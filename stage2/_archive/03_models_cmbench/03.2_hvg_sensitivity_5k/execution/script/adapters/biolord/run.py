#!/usr/bin/env python
"""biolord adapter (ticket 14).

Author recipe (research/02-biolord-perturbnet-scvidr.md §1, ticket 09):

- single joint objective; no pretraining and no released weights
  (`pert_model/biolord/author_log.json`);
- ``ordered_attributes_keys=["rdkit2d_dose"]`` and
  ``categorical_attributes_keys=["cell_type"]`` (SI Table 4);
- tuned sci-Plex config (SI Table 5): tau=1e4, sigma=20, decoder 4096x4,
  n_latent 256, ordered/categorical attribute latents 256/3, attribute NN
  2048x2, one Adam per module (latent 1e-4, decoder 1e-4, attribute 1e-2,
  wd 1e-4/1e-4/4e-8), cosine schedule to 1e-5, batch 512, 500 epochs.

biolord is **transductive**: a cell's baseline is a per-row latent code learned
at registration time (`_module.py:322-329`), not an amortised encoder.  The
adapter therefore registers one AnnData at fit time holding every training
treated well plus all 192 shared DMSO control cells, keeps the row order, and
at prediction reuses the **control rows** as the basal state: their
``rdkit2d_dose`` attribute is overwritten with the held-out drug's vector and
``module.get_expression`` is called on the registered AnnData
(``_model.py:414-416``, ``_module.py:574-599``).  This is the same construction
``_compute_prediction`` uses internally; it avoids the categorical-only
counterfactual helper (`_model.py:825-829`).

Benchmark deviations (all in ``models/biolord/FIDELITY.md``): 5,000 rank-sum HVG
genes instead of the authors' 7,561/4,999, 2,523 train wells instead of a
sci-Plex-sized corpus, random init, a 50-epoch checkpoint grid on the cold
valid drugs (ticket 09) instead of the authors' per-epoch Lightning early
stopping, and RDKit2D descriptors from descriptastorus standardised on the
training drugs instead of chemprop's pre-normalised ``rdkit_2d_normalized``
block.  The predicted rows are registered DMSO cells, which is the
``ideal_fit: partial`` transductive caveat.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from repo_root import EXECUTION, ROOT  # noqa: E402
sys.dont_write_bytecode = True

from adapter_api import Adapter, run_cli  # noqa: E402
from adapters.common import (  # noqa: E402
    N_GENES,
    Controls,
    apply_standardizer,
    case_summary,
    ensure_float32,
    fit_standardizer,
    load_json,
    load_npz,
    log,
    rdkit2d_descriptors,
    write_json,
)

BIOLORD_SRC = ROOT / "pert_model/biolord/source/src"
sys.path.insert(0, str(BIOLORD_SRC))

import anndata as ad  # noqa: E402
import torch  # noqa: E402
import torch.nn.functional as F  # noqa: E402

import biolord  # noqa: E402

ATTRIBUTE = "rdkit2d_dose"
CELL_TYPE_KEY = "cell_type"
CELL_TYPE = "ipsc_cm"
DMSO_SMILES = "CS(C)=O"
#: basename of the checkpoint state written by the grid callback
STATE_NAME = "best.pt"
LAST_NAME = "last.pt"


def _smoke() -> bool:
    return os.environ.get("BENCH_SMOKE", "").strip() not in ("", "0", "false", "False")


class BiolordAdapter(Adapter):
    name = "biolord"

    def __init__(self) -> None:
        self.cfg = load_json(EXECUTION / "models/biolord/config.json")

    # ------------------------------------------------------------------
    # Attributes and AnnData

    @staticmethod
    def _dose(case: dict[str, np.ndarray]) -> np.ndarray:
        """Raw dose in uM: the authors concatenate raw dosage (SI Table 4)."""

        return np.asarray(case["dose_um"], dtype=np.float32)

    def _fit_standardizer(self, train: dict[str, np.ndarray]) -> dict:
        """Standardise RDKit2D columns on the unique training drugs only."""

        smiles, first = np.unique(np.asarray(train["smiles"]).astype(str), return_index=True)
        descriptors = rdkit2d_descriptors(smiles)
        return fit_standardizer(descriptors, min_std=float(self.cfg["descriptor"]["min_std"]))

    def _attributes(
        self,
        smiles: np.ndarray,
        dose_um: np.ndarray,
        stats: dict,
    ) -> np.ndarray:
        descriptors = apply_standardizer(rdkit2d_descriptors(smiles), stats)
        dose = np.asarray(dose_um, dtype=np.float32).reshape(-1, 1)
        return np.concatenate([descriptors, dose], axis=1).astype(np.float32)

    def _control_attribute(self, stats: dict) -> np.ndarray:
        """DMSO row attribute: the vehicle descriptor at dose 0.

        Only the control rows ever carry it while training; at prediction the
        attribute is overwritten with the held-out drug's vector, so this
        choice only affects how well the DMSO rows reconstruct themselves.
        """

        return self._attributes(
            np.asarray([DMSO_SMILES]), np.zeros(1, dtype=np.float32), stats
        )[0]

    def _build_adata(
        self,
        train: dict[str, np.ndarray],
        controls: Controls,
        stats: dict,
    ) -> tuple[ad.AnnData, dict]:
        n_treated = len(train["well_id"])
        treated_x = np.asarray(train["log1p_cp10k"], dtype=np.float32)
        control_cells = np.asarray(controls.cells, dtype=np.float32)
        n_groups, n_reps = control_cells.shape[:2]
        control_x = control_cells.reshape(n_groups * n_reps, N_GENES)

        X = np.concatenate([treated_x, control_x], axis=0)
        treated_attr = self._attributes(train["smiles"], self._dose(train), stats)
        control_attr = np.repeat(
            self._control_attribute(stats)[None, :], n_groups * n_reps, axis=0
        )
        obsm = {ATTRIBUTE: np.concatenate([treated_attr, control_attr], axis=0)}

        obs = {
            CELL_TYPE_KEY: np.full(X.shape[0], CELL_TYPE, dtype=object),
            "well_id": np.concatenate(
                [
                    np.asarray(train["well_id"]).astype(str),
                    np.asarray(
                        [
                            f"control::{controls.groups[g]}::{r}"
                            for g in range(n_groups)
                            for r in range(n_reps)
                        ]
                    ),
                ]
            ),
            "control_group": np.concatenate(
                [
                    np.asarray(train["control_group"]).astype(str),
                    np.asarray(
                        [controls.groups[g] for g in range(n_groups) for _ in range(n_reps)]
                    ),
                ]
            ),
            "dose_um": np.concatenate(
                [self._dose(train), np.zeros(n_groups * n_reps, dtype=np.float32)]
            ),
        }
        adata = ad.AnnData(X=X, obs=obs, obsm=obsm)
        adata.obs[CELL_TYPE_KEY] = adata.obs[CELL_TYPE_KEY].astype("category")

        control_rows = {}
        for g in range(n_groups):
            for r in range(n_reps):
                row = n_treated + g * n_reps + r
                control_rows[f"{controls.groups[g]}|{r}"] = int(row)
        if adata.n_vars != N_GENES:
            raise RuntimeError(f"AnnData has {adata.n_vars} genes, expected {N_GENES}")
        return adata, control_rows

    # ------------------------------------------------------------------
    # Model

    def _make_model(self, adata: ad.AnnData, cfg: dict) -> "biolord.Biolord":
        biolord.Biolord.setup_anndata(
            adata,
            ordered_attributes_keys=[ATTRIBUTE],
            categorical_attributes_keys=[CELL_TYPE_KEY],
        )
        module_params = dict(cfg["module_params"])
        module_params["seed"] = int(cfg["seed"])
        return biolord.Biolord(
            adata=adata,
            n_latent=int(cfg["n_latent"]),
            module_params=module_params,
            train_classifiers=False,
            split_key=None,
        )

    @staticmethod
    def _load_state(model: "biolord.Biolord", path: Path, device: torch.device) -> None:
        stored = torch.load(Path(path), map_location="cpu", weights_only=False)
        model.module.load_state_dict(stored["state_dict"])
        model.module.to(device)
        model.module.eval()
        model.is_trained_ = True

    @staticmethod
    def _control_rows_for(
        case: dict[str, np.ndarray], control_rows: dict[str, int]
    ) -> np.ndarray:
        groups = [str(g) for g in case["control_group"]]
        reps = np.asarray(case["replicate"], dtype=np.int32)
        n_reps = max(int(k.rsplit("|", 1)[1]) for k in control_rows) + 1
        reps = np.mod(reps, n_reps)
        try:
            return np.asarray([control_rows[f"{g}|{int(r)}"] for g, r in zip(groups, reps)])
        except KeyError as exc:  # fail closed on a plate that was not registered
            raise KeyError(f"case control cell is not in the training AnnData: {exc}") from exc

    # ------------------------------------------------------------------
    # Adapter API

    def fit(self, train, valid, config, run_dir, controls=None):
        if controls is None:
            raise RuntimeError("biolord needs --controls (the 192 shared DMSO cells)")
        cfg = dict(self.cfg)
        if _smoke():
            cfg.update(cfg.get("smoke", {}))
        seed = int(config.get("seed", cfg["seed"]))
        cfg["seed"] = seed
        torch.manual_seed(seed)
        np.random.seed(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        train_case = load_npz(train)
        valid_case = load_npz(valid)
        control_sets = Controls(controls)
        log(
            f"device={device} train={case_summary(train_case)} "
            f"valid={case_summary(valid_case)} smoke={_smoke()}"
        )

        stats = self._fit_standardizer(train_case)
        adata, control_rows = self._build_adata(train_case, control_sets, stats)
        attribute_width = int(adata.obsm[ATTRIBUTE].shape[1])
        log(
            f"registered {adata.n_obs} rows "
            f"({len(train_case['well_id'])} treated + {len(control_rows)} controls), "
            f"attribute width {attribute_width}"
        )

        checkpoint = Path(run_dir) / "checkpoint"
        checkpoint.mkdir(parents=True, exist_ok=True)
        np.savez(
            checkpoint / "standardizer.npz",
            mean=np.asarray(stats["mean"]),
            std=np.asarray(stats["std"]),
            keep=np.asarray(stats["keep"]),
        )
        adata.write_h5ad(checkpoint / "train_adata.h5ad", compression=None)
        meta = {
            "model": "biolord",
            "n_genes": N_GENES,
            "attribute": ATTRIBUTE,
            "cell_type": CELL_TYPE,
            "attribute_width": attribute_width,
            "n_treated": int(len(train_case["well_id"])),
            "n_controls": int(len(control_rows)),
            "control_rows": control_rows,
            "control_order": sorted(control_rows, key=control_rows.get),
            "treated_well_ids": [str(w) for w in train_case["well_id"]],
            "stats_min_std": float(self.cfg["descriptor"]["min_std"]),
            "dose_transform": "raw_uM",
            "descriptor": self.cfg["descriptor"]["kind"],
        }
        write_json(checkpoint / "meta.json", meta)

        model = self._make_model(adata, cfg)
        model.module.to(device)

        def _state() -> dict[str, torch.Tensor]:
            return {
                key: value.detach().cpu().clone()
                for key, value in model.module.state_dict().items()
            }

        def _save(path: Path) -> None:
            tmp = Path(str(path) + ".tmp")
            torch.save(
                {
                    "model": "biolord",
                    "state_dict": _state(),
                    "config": cfg,
                    "seed": seed,
                    "n_genes": N_GENES,
                    "attribute_width": attribute_width,
                },
                tmp,
            )
            os.replace(tmp, path)

        # Always leave a usable checkpoint for a timed-out fit.
        _save(checkpoint / LAST_NAME)

        # Cold-valid author loss: Gaussian NLL + tau*MSE + lambda*||z_unknown||^2
        # evaluated on counterfactual predictions for the valid drugs from each
        # well's matched DMSO control row.
        valid_rows = self._control_rows_for(valid_case, control_rows)
        valid_attr = self._attributes(
            valid_case["smiles"], self._dose(valid_case), stats
        )
        valid_target = torch.from_numpy(
            np.asarray(valid_case["log1p_cp10k"], dtype=np.float32)
        )
        tau = float(cfg["module_params"]["reconstruction_penalty"])
        lam = float(cfg["module_params"]["unknown_attribute_penalty"])
        batch_size = int(cfg["batch_size"])
        grid = int(cfg["checkpoint_grid"])
        history: list[dict] = []
        best = {"loss": float("inf"), "epoch": -1}

        @torch.no_grad()
        def valid_loss() -> float:
            model.module.eval()
            nll_sum = 0.0
            mse_sum = 0.0
            pen_sum = 0.0
            elements = 0
            for start in range(0, len(valid_rows), batch_size):
                stop = min(start + batch_size, len(valid_rows))
                rows = valid_rows[start:stop]
                dataset = model.get_dataset(adata, indices=rows)
                dataset[ATTRIBUTE] = torch.from_numpy(valid_attr[start:stop])
                means, variances = model.module.get_expression(dataset)
                target = valid_target[start:stop].to(means.device)
                nll = F.gaussian_nll_loss(
                    means, target, variances, reduction="none"
                )
                mse = F.mse_loss(means, target, reduction="none")
                nll_sum += float(nll.sum().item())
                mse_sum += float(mse.sum().item())
                elements += int(nll.numel())
                codes = model.module.latent_codes(
                    torch.as_tensor(rows, dtype=torch.long, device=means.device)
                )
                pen_sum += float((codes**2).sum(dim=1).mean().item()) * len(rows)
            model.module.train()
            return nll_sum / max(elements, 1) + tau * mse_sum / max(elements, 1) + lam * pen_sum / max(
                len(valid_rows), 1
            )

        from lightning.pytorch.callbacks import Callback

        class _Grid(Callback):
            def on_train_epoch_end(self, trainer, pl_module):  # noqa: ANN001
                epoch = int(trainer.current_epoch) + 1
                if epoch % grid:
                    return
                loss = float(valid_loss())
                finite = np.isfinite(loss)
                history.append(
                    {
                        "epoch": epoch,
                        "valid_loss": loss,
                        "finite": bool(finite),
                    }
                )
                _save(checkpoint / LAST_NAME)
                if finite and loss < best["loss"]:
                    best["loss"] = loss
                    best["epoch"] = epoch
                    _save(checkpoint / STATE_NAME)
                log(
                    f"grid epoch {epoch} valid_author_loss={loss:.6f} "
                    f"best={best['loss']:.6f}@{best['epoch']}"
                )

        trainer_kwargs = {
            "callbacks": [_Grid()],
            "default_root_dir": str(Path(run_dir) / "lightning"),
            # The registered AnnData is entirely training data; validation is
            # the cold-drug counterfactual loss computed by the callback.
            "limit_val_batches": 0,
            "enable_progress_bar": True,
        }
        if cfg.get("limit_train_batches"):
            trainer_kwargs["limit_train_batches"] = int(cfg["limit_train_batches"])
        model.train(
            max_epochs=int(cfg["max_epochs"]),
            batch_size=batch_size,
            accelerator="gpu" if device.type == "cuda" else "cpu",
            device=1,  # "1 device"; parse_device_args maps it to cuda:0 / cpu
            train_size=1.0,
            validation_size=0.0,
            early_stopping=False,
            enable_checkpointing=False,
            plan_kwargs=dict(cfg["plan_kwargs"]),
            num_workers=0,
            **trainer_kwargs,
        )

        if not (checkpoint / STATE_NAME).exists():
            log("no finite grid checkpoint; keeping the last state")
            _save(checkpoint / STATE_NAME)
            best = {"loss": float("nan"), "epoch": int(cfg["max_epochs"])}
        write_json(
            Path(run_dir) / "fit_history.json",
            {
                "grid": grid,
                "history": history,
                "best": best,
                "smoke": _smoke(),
                "n_registered": int(adata.n_obs),
                "attribute_width": attribute_width,
                "epoch_history": getattr(model.training_plan, "epoch_history", {}),
            },
        )
        size_mb = (checkpoint / STATE_NAME).stat().st_size / 1e6
        log(
            f"checkpoint {checkpoint} best valid {best['loss']:.6f} "
            f"at epoch {best['epoch']} ({size_mb:.1f} MB state)"
        )
        return checkpoint

    def predict(self, cases, checkpoint, run_dir, controls=None):
        if controls is None:
            raise RuntimeError("biolord needs --controls (the 192 shared DMSO cells)")
        checkpoint = Path(checkpoint)
        state_path = checkpoint / STATE_NAME
        if not state_path.exists():
            state_path = checkpoint / LAST_NAME
        stored = torch.load(state_path, map_location="cpu", weights_only=False)
        cfg = dict(stored["config"])
        meta = load_json(checkpoint / "meta.json")
        stats_npz = np.load(checkpoint / "standardizer.npz", allow_pickle=False)
        try:
            stats = {
                "mean": np.asarray(stats_npz["mean"]),
                "std": np.asarray(stats_npz["std"]),
                "keep": np.asarray(stats_npz["keep"]),
            }
        finally:
            stats_npz.close()
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        adata = ad.read_h5ad(checkpoint / "train_adata.h5ad")
        model = self._make_model(adata, cfg)
        self._load_state(model, state_path, device)

        case = load_npz(cases)
        control_sets = Controls(controls)
        rows = self._control_rows_for(case, meta["control_rows"])
        attribute = self._attributes(case["smiles"], self._dose(case), stats)
        if attribute.shape[1] != int(meta["attribute_width"]):
            raise RuntimeError(
                f"attribute width {attribute.shape[1]} != checkpoint {meta['attribute_width']}"
            )

        batch_size = int(cfg["batch_size"])
        out = np.empty((len(rows), N_GENES), dtype=np.float32)
        with torch.no_grad():
            for start in range(0, len(rows), batch_size):
                stop = min(start + batch_size, len(rows))
                dataset = model.get_dataset(adata, indices=rows[start:stop])
                dataset[ATTRIBUTE] = torch.from_numpy(attribute[start:stop])
                means, _ = model.module.get_expression(dataset)
                out[start:stop] = means.detach().cpu().numpy()
        basal = control_sets.mean_for(case["control_group"])
        delta = out - basal
        log(
            f"predict {delta.shape} finite={np.isfinite(delta).all()} "
            f"source=registered DMSO cells basal=plate DMSO mean"
        )
        return ensure_float32(delta)


if __name__ == "__main__":
    raise SystemExit(run_cli(BiolordAdapter()))
