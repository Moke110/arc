#!/usr/bin/env python
"""Squidiff adapter (ticket 15).

Author recipe (`.scratch/stage2-03-bench-models/research/03-squidiff-cmonge-decant.md`
§1, cited to the vendored source):

- one control-conditioned diffusion MLP; encoder input is
  ``[control expression ; rFCFP4(1024) x log10(dose_uM + 1)]`` and the denoiser
  is ``x_t + t-embedding + 60-d semantic latent -> 2048 x 3 -> 5,000``
  (``Squidiff/MLPModel.py:137-207``).
- epsilon-prediction MSE with uniform timestep sampling, 1000 diffusion steps,
  linear beta schedule (``Squidiff/diffusion.py:735-810``,
  ``train_squidiff.py:112,129``).
- AdamW lr 1e-4, linear anneal to 0 over ``lr_anneal_steps``, wd 0.0, batch 64,
  **fixed 100,000 steps, no early stopping**, EMA 0.9999
  (``train_squidiff.py:109-135``, ``Squidiff/train_util.py:107-119,188-286``).
- inference precomputes the semantic latent once and runs DDIM from pure noise
  (``sample_squidiff.py:148-158``).

Benchmark-specific choices (all in ``models/squidiff/FIDELITY.md``): 5,000
genes instead of the authors' 200 (a 90x extrapolation), the fold's training
wells instead of the full SciPlex atlas, random init (no public weights), a
10% checkpoint grid evaluated with the author loss on the cold valid drugs, and
one matched DMSO cell per well as the paired control (the corpus has no
``paired_control_index``).
"""

from __future__ import annotations

import os
import signal
import sys
import time
from pathlib import Path

# Model code is vendored read-only under `pert_model/`; never leave __pycache__.
sys.dont_write_bytecode = True
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from repo_root import EXECUTION, ROOT  # noqa: E402

from adapter_api import Adapter, run_cli  # noqa: E402
from adapters.common import (  # noqa: E402
    N_GENES,
    Controls,
    case_summary,
    ensure_float32,
    load_json,
    load_npz,
    log,
    morgan_fingerprints,
    write_json,
)

SQUIDIFF_SRC = ROOT / "pert_model/Squidiff/source"
sys.path.insert(0, str(SQUIDIFF_SRC))

from Squidiff.nn import update_ema  # noqa: E402
from Squidiff.script_util import create_gaussian_diffusion, create_model  # noqa: E402


def _smoke() -> bool:
    return os.environ.get("BENCH_SMOKE", "").strip() not in ("", "0", "false", "False")


def _grid(total_steps: int, slots: int) -> list[int]:
    """The ticket-09 fixed checkpoint grid: ``slots`` evenly spaced steps."""

    steps = sorted({max(1, int(round(total_steps * i / slots))) for i in range(1, slots + 1)})
    return steps


def _save_atomic(payload, path: Path) -> None:
    """Write a checkpoint through a temporary file, then rename it."""

    path = Path(path)
    tmp = Path(str(path) + ".tmp")
    torch.save(payload, tmp)
    os.replace(tmp, path)


class SquidiffAdapter(Adapter):
    name = "squidiff"

    def __init__(self) -> None:
        self.cfg = load_json(EXECUTION / "models/squidiff/config.json")

    # -- model / diffusion ---------------------------------------------

    def _build(self, cfg: dict, device: torch.device) -> torch.nn.Module:
        """The authors' ``create_model`` with the 5,000-gene, drug-structure flags.

        ``num_layers``, the 2048-d hidden width, the 60-d latent and the 2048-d
        time embedding are author-fixed defaults inside ``create_model`` /
        ``MLPModel``; only ``gene_size``/``output_dim`` and the drug flags are
        passed, exactly as in ``train_squidiff.py``.
        """

        model = create_model(
            int(cfg["gene_size"]),
            num_layers=int(cfg["num_layers"]),
            output_dim=int(cfg["output_dim"]),
            class_cond=bool(cfg["class_cond"]),
            use_encoder=bool(cfg["use_encoder"]),
            use_drug_structure=bool(cfg["use_drug_structure"]),
            drug_dimension=int(cfg["drug_dimension"]),
            comb_num=int(cfg["comb_num"]),
        )
        return model.to(device)

    def _diffusion(self, cfg: dict, timestep_respacing: str = ""):
        return create_gaussian_diffusion(
            steps=int(cfg["diffusion_steps"]),
            noise_schedule=str(cfg["noise_schedule"]),
            timestep_respacing=timestep_respacing,
            use_encoder=bool(cfg["use_encoder"]),
        )

    # -- inputs ---------------------------------------------------------

    def _condition(
        self, case: dict[str, np.ndarray], controls: Controls, cfg: dict
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Paired DMSO cell + the authors' drug condition.

        Pairing rule: every well takes the DMSO cell on its own plate selected by
        the well's 0-based ``replicate`` index (``Controls.cell_for``), i.e. one
        matched control cell per treated well.  The corpus has no
        ``paired_control_index``, so this is the closest reproducible rule.

        Drug condition: ``rFCFP4(1024) * log10(dose_uM + 1)`` built with the
        shared PRnet fingerprint helper -- byte-for-byte the authors'
        ``Drug_dose_encoder``.
        """

        control = controls.cell_for(case["control_group"], case["replicate"])
        fingerprints = morgan_fingerprints(
            case["smiles"],
            n_bits=int(cfg["drug_dimension"]),
            radius=int(cfg["fingerprint"]["radius"]),
            use_features=True,
        )
        dose = np.log10(case["dose_um"].astype(np.float64) + 1.0).astype(np.float32)
        return (
            torch.from_numpy(control.astype(np.float32)),
            torch.from_numpy((fingerprints * dose[:, None]).astype(np.float32)),
        )

    # -- validation loss (author epsilon-MSE, deterministic) -------------

    def _valid_loss(
        self,
        model: torch.nn.Module,
        diffusion,
        control: torch.Tensor,
        drug: torch.Tensor,
        target: torch.Tensor,
        timesteps: list[int],
        noise: torch.Tensor,
        batch_size: int,
        device: torch.device,
    ) -> float:
        model.eval()
        total, seen = 0.0, 0
        with torch.no_grad():
            for start in range(0, len(target), batch_size):
                stop = min(start + batch_size, len(target))
                x = target[start:stop].to(device)
                cond = {
                    "group": None,
                    "drug_dose": drug[start:stop].to(device),
                    "control_feature": control[start:stop].to(device),
                }
                noise_batch = noise[start:stop].to(device)
                for timestep in timesteps:
                    t = torch.full((x.shape[0],), int(timestep), device=device, dtype=torch.long)
                    terms = diffusion.training_losses(
                        model, x, t, model_kwargs=dict(cond), noise=noise_batch
                    )
                    total += float(terms["loss"].mean().item()) * x.shape[0]
                    seen += x.shape[0]
        model.train()
        return total / max(seen, 1)

    # -- parameter swapping for EMA evaluation --------------------------

    def _swap_in(self, model: torch.nn.Module, values: list[torch.Tensor]) -> list[torch.Tensor]:
        backup = [p.detach().clone() for p in model.parameters()]
        with torch.no_grad():
            for param, value in zip(model.parameters(), values):
                param.copy_(value)
        return backup

    def _swap_out(self, model: torch.nn.Module, backup: list[torch.Tensor]) -> None:
        with torch.no_grad():
            for param, value in zip(model.parameters(), backup):
                param.copy_(value)
        del backup

    def _state_dict(self, model: torch.nn.Module, values: list[torch.Tensor]) -> dict:
        state = {name: tensor.detach().clone() for name, tensor in model.state_dict().items()}
        for (name, _), value in zip(model.named_parameters(), values):
            state[name] = value.detach().clone()
        return state

    # -- Adapter API ----------------------------------------------------

    def fit(self, train, valid, config, run_dir, controls=None):
        if controls is None:
            raise RuntimeError("Squidiff needs --controls (the shared DMSO cells)")
        smoke = _smoke()
        cfg = dict(self.cfg)
        if smoke:
            cfg.update(cfg.get("smoke", {}))
        seed = int(config.get("seed", cfg["seed"]))
        torch.manual_seed(seed)
        np.random.seed(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        torch.use_deterministic_algorithms(True, warn_only=True)
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        train_case = load_npz(train)
        valid_case = load_npz(valid)
        control_sets = Controls(controls)
        log(
            f"device={device} smoke={smoke} train={case_summary(train_case)} "
            f"valid={case_summary(valid_case)}"
        )

        control_train, drug_train = self._condition(train_case, control_sets, cfg)
        target_train = torch.from_numpy(np.asarray(train_case["log1p_cp10k"], dtype=np.float32))
        control_valid, drug_valid = self._condition(valid_case, control_sets, cfg)
        target_valid = torch.from_numpy(np.asarray(valid_case["log1p_cp10k"], dtype=np.float32))
        del train_case, valid_case

        model = self._build(cfg, device)
        model.train()
        n_params = sum(p.numel() for p in model.parameters())
        log(f"model parameters={n_params:,}")

        optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=float(cfg["optimizer"]["lr"]),
            weight_decay=float(cfg["optimizer"]["weight_decay"]),
        )
        ema_rate = float(cfg["ema_rate"])
        ema_values = [p.detach().clone() for p in model.parameters()]

        diffusion_train = self._diffusion(cfg)
        batch_size = int(cfg["batch_size"])
        total_steps = int(cfg["total_steps"])
        slots = int(cfg["checkpoint_grid"])
        grid = _grid(total_steps, slots)
        log_every = max(1, total_steps // 20)

        valid_timesteps = [
            int(round(i * (int(cfg["diffusion_steps"]) - 1) / (int(cfg["valid_timesteps"]) - 1)))
            for i in range(int(cfg["valid_timesteps"]))
        ]
        noise_generator = torch.Generator(device="cpu").manual_seed(seed + 1)
        valid_noise = torch.randn(target_valid.shape, generator=noise_generator)

        checkpoint = Path(run_dir) / "checkpoint"
        checkpoint.mkdir(parents=True, exist_ok=True)
        best_path = checkpoint / "model.pt"
        best_valid = float("inf")
        best_record: dict = {}
        history: list[dict] = []

        stop = {"requested": False}
        previous_handler = signal.getsignal(signal.SIGTERM)

        def _request_stop(signum, frame):  # pragma: no cover - signal path
            stop["requested"] = True

        signal.signal(signal.SIGTERM, _request_stop)

        # The authors' shuffled-stream data order: a fresh permutation whenever
        # the cursor is exhausted (`train_util.py:195` rebuilds the loader each
        # step, which is the same distribution without the per-step reshuffle).
        order_generator = torch.Generator(device="cpu").manual_seed(seed + 2)
        order = torch.randperm(len(target_train), generator=order_generator).tolist()
        cursor = 0
        started = time.time()
        lr_now = float(cfg["optimizer"]["lr"])
        steps_run = 0

        try:
            for step in range(total_steps):
                steps_run = step + 1
                if stop["requested"]:
                    steps_run = step
                    log(f"SIGTERM at step {step}; stopping and keeping the best checkpoint")
                    break
                if cursor + batch_size > len(order):
                    order = torch.randperm(len(target_train), generator=order_generator).tolist()
                    cursor = 0
                index = order[cursor : cursor + batch_size]
                cursor += batch_size

                x = target_train[index].to(device)
                cond = {
                    "group": None,
                    "drug_dose": drug_train[index].to(device),
                    "control_feature": control_train[index].to(device),
                }
                t = torch.randint(
                    0, int(cfg["diffusion_steps"]), (x.shape[0],), device=device
                )
                terms = diffusion_train.training_losses(model, x, t, model_kwargs=cond)
                loss = terms["loss"].mean()

                for group in optimizer.param_groups:
                    group["lr"] = lr_now
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                optimizer.step()
                update_ema(ema_values, list(model.parameters()), rate=ema_rate)
                # The authors anneal after the step with the pre-increment index
                # (`train_util.py:280-286`).
                lr_now = float(cfg["optimizer"]["lr"]) * (
                    1.0 - step / max(total_steps, 1)
                )

                if step % log_every == 0 or (step + 1) in grid:
                    log(
                        f"step {step + 1}/{total_steps} loss={float(loss.item()):.6f} "
                        f"lr={lr_now:.2e}"
                    )

                if (step + 1) in grid:
                    candidates = {"raw": list(model.parameters()), "ema": ema_values}
                    scores = {}
                    for name, values in candidates.items():
                        backup = self._swap_in(model, values)
                        scores[name] = self._valid_loss(
                            model,
                            diffusion_train,
                            control_valid,
                            drug_valid,
                            target_valid,
                            valid_timesteps,
                            valid_noise,
                            batch_size,
                            device,
                        )
                        self._swap_out(model, backup)
                    entry = {"step": step + 1, "valid_loss": scores}
                    history.append(entry)
                    log(
                        f"grid step {step + 1}: valid raw={scores['raw']:.6f} "
                        f"ema={scores['ema']:.6f} (best {best_valid:.6f})"
                    )
                    chosen = min(scores, key=scores.get)
                    if scores[chosen] < best_valid:
                        best_valid = scores[chosen]
                        values = candidates[chosen]
                        state = self._state_dict(model, values)
                        _save_atomic(
                            {
                                "model": "squidiff",
                                "state_dict": state,
                                "config": cfg,
                                "seed": seed,
                                "n_genes": N_GENES,
                                "step": step + 1,
                                "weights": chosen,
                                "valid_loss": best_valid,
                                "valid_timesteps": valid_timesteps,
                            },
                            best_path,
                        )
                        best_record = {
                            "step": step + 1,
                            "weights": chosen,
                            "valid_loss": best_valid,
                        }
                        del state
        finally:
            signal.signal(signal.SIGTERM, previous_handler)

        if not best_path.exists():
            # Only reachable if SIGTERM arrived before the first grid point.
            best_valid = self._valid_loss(
                model,
                diffusion_train,
                control_valid,
                drug_valid,
                target_valid,
                valid_timesteps,
                valid_noise,
                batch_size,
                device,
            )
            _save_atomic(
                {
                    "model": "squidiff",
                    "state_dict": self._state_dict(model, ema_values),
                    "config": cfg,
                    "seed": seed,
                    "n_genes": N_GENES,
                    "step": 0,
                    "weights": "ema",
                    "valid_loss": best_valid,
                    "valid_timesteps": valid_timesteps,
                },
                best_path,
            )
            best_record = {"step": 0, "weights": "ema", "valid_loss": best_valid}

        wall = time.time() - started
        peak_gpu_mb = (
            torch.cuda.max_memory_allocated(device) / (1024.0 ** 2) if device.type == "cuda" else 0.0
        )
        write_json(
            Path(run_dir) / "fit_history.json",
            {
                "history": history,
                "best": best_record,
                "grid": grid,
                "total_steps": total_steps,
                "steps_run": steps_run,
                "smoke": smoke,
                "wall_s": wall,
                "n_parameters": n_params,
                "checkpoint_bytes": best_path.stat().st_size,
                "peak_gpu_mb": peak_gpu_mb,
            },
        )
        log(
            f"checkpoint {best_path} ({best_path.stat().st_size / 1e6:.0f} MB, "
            f"{best_record['weights']} @ step {best_record['step']}, "
            f"valid {best_valid:.6f}) in {wall:.1f}s, peak GPU {peak_gpu_mb:.0f} MB"
        )
        return checkpoint

    def predict(self, cases, checkpoint, run_dir, controls=None):
        if controls is None:
            raise RuntimeError("Squidiff needs --controls (the shared DMSO cells)")
        stored = torch.load(Path(checkpoint) / "model.pt", map_location="cpu", weights_only=True)
        cfg = dict(stored["config"])
        seed = int(stored.get("seed", cfg["seed"]))
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model = self._build(cfg, device)
        model.load_state_dict(stored["state_dict"])
        model.eval()

        sample_steps = int(cfg.get("sample_steps", cfg["diffusion_steps"]))
        respacing = "" if sample_steps >= int(cfg["diffusion_steps"]) else f"ddim{sample_steps}"
        diffusion = self._diffusion(cfg, timestep_respacing=respacing)
        log(
            f"sampling with {diffusion.num_timesteps} DDIM steps "
            f"(weights={stored.get('weights')} @ step {stored.get('step')})"
        )

        case = load_npz(cases)
        control_sets = Controls(controls)
        control, drug = self._condition(case, control_sets, cfg)
        out = np.empty((len(control), N_GENES), dtype=np.float32)
        batch_size = int(cfg["batch_size"])
        noise_generator = torch.Generator(device=device).manual_seed(seed + 3)
        with torch.no_grad():
            for start in range(0, len(control), batch_size):
                stop = min(start + batch_size, len(control))
                control_batch = control[start:stop].to(device)
                drug_batch = drug[start:stop].to(device)
                latent = model.encoder(
                    control_batch,
                    label=None,
                    drug_dose=drug_batch,
                    control_feature=control_batch,
                )
                noise = torch.randn(
                    (control_batch.shape[0], N_GENES),
                    generator=noise_generator,
                    device=device,
                )
                generated = diffusion.ddim_sample_loop(
                    model,
                    shape=(control_batch.shape[0], N_GENES),
                    noise=noise,
                    model_kwargs={"z_mod": latent},
                )
                out[start:stop] = generated.detach().cpu().numpy()
        # Δ is the generated absolute log1p_cp10k minus the matched plate DMSO mean.
        delta = out - control_sets.mean_for(case["control_group"])
        if device.type == "cuda":
            log(f"peak GPU {torch.cuda.max_memory_allocated(device) / (1024.0 ** 2):.0f} MB")
        log(f"predict {delta.shape} finite={np.isfinite(delta).all()}")
        return ensure_float32(delta)


if __name__ == "__main__":
    raise SystemExit(run_cli(SquidiffAdapter()))
