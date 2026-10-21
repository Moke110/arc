#!/usr/bin/env python
"""PRnet adapter (ticket 13).

Author recipe (from `.scratch/stage2-03-bench-models/research/01-...`):

- architecture: encoder ``x + comb -> 128``, ``z=64``; drug adaptor
  ``1024 rFCFP x log10(dose+1) -> 128 -> 64``; decoder ``64+64+10 -> 128 -> 2*18164``;
  dropout 0.05; Gaussian NLL loss (author `loss=['GUSS']`).
- Adam lr 1e-3, wd 1e-8; ReduceLROnPlateau x0.5, patience 10; batch 512;
  100 epochs with early stopping after 20 epochs without improvement.

The benchmark-specific changes (all recorded in `models/prnet/FIDELITY.md`):
the model is trained on the frozen cp10k log1p space instead of PRnet's
per-cell median-normalised counts (the test bundle has no library sizes, so an
exact inverse is impossible), and the input control is the matched plate DMSO
mean (the authors pool every dose-0 cell into every split, which would leak
across folds).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import torch

BENCH = Path(__file__).resolve().parents[2]
ROOT = BENCH.parents[1]
sys.path.insert(0, str(BENCH))

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

PRNET_SRC = ROOT / "pert_model/prnet"
sys.path.insert(0, str(PRNET_SRC))

from models.PRnet import PRnet  # noqa: E402


def _smoke() -> bool:
    return os.environ.get("BENCH_SMOKE", "").strip() not in ("", "0", "false", "False")


class PRnetAdapter(Adapter):
    name = "prnet"

    def __init__(self) -> None:
        self.cfg = load_json(BENCH / "models/prnet/config.json")

    # -- model ---------------------------------------------------------

    def _build(self, device: torch.device) -> torch.nn.Module:
        model = PRnet(
            None,
            x_dimension=self.cfg["x_dimension"],
            hidden_layer_sizes=list(self.cfg["hidden_layer_sizes"]),
            z_dimension=self.cfg["z_dimension"],
            adaptor_layer_sizes=list(self.cfg["adaptor_layer_sizes"]),
            comb_dimension=self.cfg["comb_dimension"],
            comb_num=self.cfg["comb_num"],
            drug_dimension=self.cfg["drug_dimension"],
            dr_rate=self.cfg["dr_rate"],
        ).get_PGM()
        model = model.to(device)
        # The author applies this weight init after moving to the device.
        for module in model.modules():
            if isinstance(module, torch.nn.Linear):
                module.weight.data.normal_(0, 0.02)
                if module.bias is not None:
                    module.bias.data.zero_()
            elif isinstance(module, torch.nn.BatchNorm1d):
                module.weight.data.normal_(1, 0.02)
                if module.bias is not None:
                    module.bias.data.zero_()
        return model

    # -- tensors -------------------------------------------------------

    def _features(
        self,
        case: dict[str, np.ndarray],
        controls: Controls,
        with_target: bool = True,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor | None]:
        control = controls.mean_for(case["control_group"])
        fingerprints = morgan_fingerprints(
            case["smiles"],
            n_bits=self.cfg["drug_dimension"],
            radius=2,
            use_features=True,
        )
        dose = np.log10(case["dose_um"].astype(np.float64) + 1.0).astype(np.float32)
        label = fingerprints * dose[:, None]
        target = None
        if with_target:
            target = np.asarray(case["log1p_cp10k"], dtype=np.float32)
        return (
            torch.from_numpy(control),
            torch.from_numpy(label),
            torch.from_numpy(target) if target is not None else None,
        )

    def _loss(
        self,
        model: torch.nn.Module,
        control: torch.Tensor,
        label: torch.Tensor,
        target: torch.Tensor,
        device: torch.device,
        noise_dim: int,
        generator: torch.Generator | None = None,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        if generator is None:
            noise = torch.zeros(control.shape[0], noise_dim)
        else:
            noise = torch.randn(control.shape[0], noise_dim, generator=generator)
        noise = noise.to(device)
        reconstruction = model(control, label, noise)
        dim = reconstruction.shape[1] // 2
        mean = reconstruction[:, :dim]
        var = torch.nn.functional.softplus(reconstruction[:, dim:])
        loss = torch.nn.functional.gaussian_nll_loss(mean, target, var, reduction="mean")
        return loss, mean, var

    # -- Adapter API ---------------------------------------------------

    def fit(self, train, valid, config, run_dir, controls=None):
        if controls is None:
            raise RuntimeError("PRnet needs --controls (the shared DMSO means)")
        cfg = dict(self.cfg)
        if _smoke():
            cfg.update(cfg.get("smoke", {}))
        seed = int(config.get("seed", cfg["seed"]))
        torch.manual_seed(seed)
        np.random.seed(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        train_case = load_npz(train)
        valid_case = load_npz(valid)
        control_sets = Controls(controls)
        log(f"device={device} train={case_summary(train_case)} valid={case_summary(valid_case)}")

        x_train, c_train, y_train = self._features(train_case, control_sets)
        x_valid, c_valid, y_valid = self._features(valid_case, control_sets)
        del train_case, valid_case

        model = self._build(device)
        optimizer = torch.optim.Adam(
            model.parameters(),
            lr=float(cfg["optimizer"]["lr"]),
            weight_decay=float(cfg["optimizer"]["weight_decay"]),
        )
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            optimizer,
            mode="min",
            factor=float(cfg["scheduler"]["factor"]),
            patience=int(cfg["scheduler"]["patience"]),
            min_lr=float(cfg["scheduler"]["min_lr"]),
        )
        batch_size = int(cfg["batch_size"])
        max_epochs = int(cfg["max_epochs"])
        patience = int(cfg["early_stop_patience"])
        noise_dim = int(cfg["noise_dim"])
        generator = torch.Generator(device="cpu").manual_seed(seed)
        max_batches = int(cfg["smoke"].get("max_batches", 10**9)) if _smoke() else 10**9

        checkpoint = Path(run_dir) / "checkpoint"
        checkpoint.mkdir(parents=True, exist_ok=True)
        best_path = checkpoint / "model.pt"
        best_valid = float("inf")
        stale = 0
        history: list[dict] = []

        for epoch in range(max_epochs):
            model.train()
            order = torch.randperm(len(x_train), generator=generator).numpy()
            total, seen = 0.0, 0
            for start in range(0, len(order), batch_size):
                if start // batch_size >= max_batches:
                    break
                idx = order[start : start + batch_size]
                control = x_train[idx].to(device)
                label = c_train[idx].to(device)
                target = y_train[idx].to(device)
                loss, _, _ = self._loss(
                    model, control, label, target, device, noise_dim, generator
                )
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                optimizer.step()
                total += float(loss.item()) * len(idx)
                seen += len(idx)

            model.eval()
            valid_total, valid_seen = 0.0, 0
            with torch.no_grad():
                for start in range(0, len(x_valid), batch_size):
                    idx = slice(start, min(start + batch_size, len(x_valid)))
                    control = x_valid[idx].to(device)
                    label = c_valid[idx].to(device)
                    target = y_valid[idx].to(device)
                    loss, _, _ = self._loss(
                        model, control, label, target, device, noise_dim, None
                    )
                    n = control.shape[0]
                    valid_total += float(loss.item()) * n
                    valid_seen += n

            train_loss = total / max(seen, 1)
            valid_loss = valid_total / max(valid_seen, 1)
            scheduler.step(valid_loss)
            history.append(
                {
                    "epoch": epoch,
                    "train_loss": train_loss,
                    "valid_loss": valid_loss,
                    "lr": float(optimizer.param_groups[0]["lr"]),
                }
            )
            log(
                f"epoch {epoch + 1}/{max_epochs} train={train_loss:.6f} "
                f"valid={valid_loss:.6f} lr={optimizer.param_groups[0]['lr']:.2e}"
            )

            if valid_loss < best_valid:
                best_valid = valid_loss
                stale = 0
                torch.save(
                    {
                        "model": "prnet",
                        "state_dict": model.state_dict(),
                        "config": cfg,
                        "seed": seed,
                        "n_genes": N_GENES,
                        "best_valid_loss": best_valid,
                        "epoch": epoch,
                    },
                    best_path,
                )
            else:
                stale += 1
                if not _smoke() and stale > patience:
                    log(f"early stop after {stale} epochs without improvement")
                    break

        if not best_path.exists():
            torch.save(
                {
                    "model": "prnet",
                    "state_dict": model.state_dict(),
                    "config": cfg,
                    "seed": seed,
                    "n_genes": N_GENES,
                    "best_valid_loss": best_valid,
                    "epoch": max_epochs - 1,
                },
                best_path,
            )

        write_json(
            Path(run_dir) / "fit_history.json",
            {"history": history, "best_valid_loss": best_valid, "smoke": _smoke()},
        )
        log(f"checkpoint {best_path} (best valid {best_valid:.6f})")
        return checkpoint

    def predict(self, cases, checkpoint, run_dir, controls=None):
        if controls is None:
            raise RuntimeError("PRnet needs --controls (the shared DMSO means)")
        stored = torch.load(Path(checkpoint) / "model.pt", map_location="cpu")
        cfg = stored["config"]
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model = self._build(device)
        model.load_state_dict(stored["state_dict"])
        model.eval()

        case = load_npz(cases)
        control_sets = Controls(controls)
        x, c, _ = self._features(case, control_sets, with_target=False)
        out = np.empty((len(x), N_GENES), dtype=np.float32)
        batch_size = int(cfg["batch_size"])
        with torch.no_grad():
            for start in range(0, len(x), batch_size):
                stop = min(start + batch_size, len(x))
                control = x[start:stop].to(device)
                label = c[start:stop].to(device)
                noise = torch.zeros(control.shape[0], int(cfg["noise_dim"]), device=device)
                reconstruction = model(control, label, noise)
                mean = reconstruction[:, :N_GENES]
                out[start:stop] = (mean - control).cpu().numpy()
        log(f"predict {out.shape} finite={np.isfinite(out).all()}")
        return ensure_float32(out)


if __name__ == "__main__":
    raise SystemExit(run_cli(PRnetAdapter()))
