#!/usr/bin/env python
"""CPA adapter (ticket 13).

Author recipe (research/01-prrnnet-cpa-chemcpa §2, ticket 09):

- CPA (scvi-tools 0.20.3): Gaussian reconstruction on normalised/log data,
  ``n_latent=128``, encoder/decoder 512x3, dropout 0.1, ``logsigm`` doser,
  3 Adam optimizers (AE 5e-4 / adversary 3e-4 / doser 3e-4, wd 1e-6/4e-7/4e-7),
  StepLR(45, 0.9), batch 128, 2,000 epochs, early stop patience 10.
- Cold drugs: ``use_rdkit_embeddings=True`` (frozen 2,048-bit Morgan + shared
  Linear).  ``setup_anndata`` runs on the full fold AnnData (train + valid +
  test metadata rows) so every drug has an embedding row; only training rows
  are optimised.

Benchmark deviations (recorded in ``models/cpa/FIDELITY.md``): no 978/5k HVG
selection (5,000 genes), the cold drug's per-drug dose scaler row is left at
its β=1/bias=0 initialisation because it was never trained, and the author's
whole-object HVG step is skipped as the gene space is fixed by the corpus.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from repo_root import EXECUTION, ROOT  # noqa: E402
# The fake `cpa-tools` dist-info stub that used to make the import work was
# deleted in ticket 06; the adapter supplies the version fallback itself.
sys.path.insert(0, str(ROOT / "pert_model/CPA/source"))

from adapter_api import Adapter, run_cli  # noqa: E402
from adapters.common import (  # noqa: E402
    N_GENES,
    Controls,
    case_summary,
    ensure_float32,
    load_json,
    load_npz,
    log,
    write_json,
)

import cpa  # noqa: E402


def _smoke() -> bool:
    return os.environ.get("BENCH_SMOKE", "").strip() not in ("", "0", "false", "False")


def _var_names() -> pd.Index:
    return pd.Index([f"gene_{i:05d}" for i in range(N_GENES)], dtype=object)


class CPAAdapter(Adapter):
    name = "cpa"

    def __init__(self) -> None:
        self.cfg = load_json(EXECUTION / "models/cpa/config.json")

    # -- fold AnnData --------------------------------------------------

    def _fold_adata(
        self,
        train: Path,
        valid: Path,
        test: Path,
        controls: Controls,
    ):
        import anndata as ad

        cases = {
            "train": load_npz(train),
            "valid": load_npz(valid),
            "test": load_npz(test),
        }
        log(
            "fold cases: "
            + " ".join(f"{role}={case_summary(case)['n_cases']}" for role, case in cases.items())
        )

        blocks = []
        obs = []
        for role in ("train", "valid", "test"):
            case = cases[role]
            n = len(case["well_id"])
            if role == "test":
                expression = np.zeros((n, N_GENES), dtype=np.float32)
            else:
                expression = np.asarray(case["log1p_cp10k"], dtype=np.float32)
            blocks.append(expression)
            obs.append(
                pd.DataFrame(
                    {
                        "condition": np.asarray(case["drug_key"]).astype(str),
                        "dose": np.asarray(case["dose_um"], dtype=np.float32),
                        "smiles": np.asarray(case["smiles"]).astype(str),
                        "cell_type": np.asarray(case["cell_type"]).astype(str),
                        "split": role,
                        "well_id": np.asarray(case["well_id"]).astype(str),
                    }
                )
            )

        # Controls are never held out: they live in the training split.
        cells = controls.cells.reshape(-1, N_GENES)
        blocks.append(cells.astype(np.float32))
        n_controls = cells.shape[0]
        obs.append(
            pd.DataFrame(
                {
                    "condition": np.full(n_controls, "DMSO", dtype=object),
                    "dose": np.zeros(n_controls, dtype=np.float32),
                    "smiles": np.full(n_controls, "", dtype=object),
                    "cell_type": np.full(n_controls, "ipsc_cm", dtype=object),
                    "split": np.full(n_controls, "train", dtype=object),
                    "well_id": np.asarray(
                        [f"control_{i:04d}" for i in range(n_controls)], dtype=object
                    ),
                }
            )
        )

        adata = ad.AnnData(
            X=np.vstack(blocks),
            obs=pd.concat(obs, ignore_index=True),
            var=pd.DataFrame(index=_var_names()),
        )
        adata.obs_names = adata.obs["well_id"].astype(str).to_numpy()
        return adata

    def _pred_adata(self, case: dict[str, np.ndarray], controls: Controls):
        import anndata as ad

        basal = controls.mean_for(case["control_group"])
        obs = pd.DataFrame(
            {
                "condition": np.asarray(case["drug_key"]).astype(str),
                "dose": np.asarray(case["dose_um"], dtype=np.float32),
                "smiles": np.asarray(case["smiles"]).astype(str),
                "cell_type": np.asarray(case["cell_type"]).astype(str),
                "split": np.full(len(case["well_id"]), "test", dtype=object),
                "well_id": np.asarray(case["well_id"]).astype(str),
            }
        )
        adata = ad.AnnData(
            X=np.asarray(basal, dtype=np.float32),
            obs=obs,
            var=pd.DataFrame(index=_var_names()),
        )
        adata.obs_names = adata.obs["well_id"].astype(str).to_numpy()
        return adata

    # -- Adapter API ---------------------------------------------------

    def fit(self, train, valid, config, run_dir, controls=None):
        if controls is None:
            raise RuntimeError("CPA needs --controls (the shared DMSO means)")
        cfg = dict(self.cfg)
        if _smoke():
            cfg.update(cfg.get("smoke", {}))
        seed = int(config.get("seed", cfg["seed"]))
        torch.manual_seed(seed)
        np.random.seed(seed)
        control_sets = Controls(controls)

        # CPA keeps the drug dictionary as class attributes; clear any state
        # from a previous fold in the same process (ticket 01 pitfall).
        cpa.CPA.pert_encoder = None
        cpa.CPA.covars_encoder = None
        cpa.CPA.pert_smiles_map = None

        adata = self._fold_adata(train, valid, Path(config["cases"]["test"]), control_sets)
        log(f"fold AnnData {adata.shape} (train/valid/test/control rows)")

        cpa.CPA.setup_anndata(
            adata,
            perturbation_key="condition",
            control_group="DMSO",
            dosage_key="dose",
            categorical_covariate_keys=["cell_type"],
            smiles_key="smiles",
            is_count_data=False,
            max_comb_len=1,
        )

        model = cpa.CPA(
            adata,
            split_key="split",
            train_split="train",
            valid_split="valid",
            test_split="test",
            use_rdkit_embeddings=True,
            n_latent=int(cfg["n_latent"]),
            recon_loss=str(cfg["recon_loss"]),
            doser_type=str(cfg["doser_type"]),
            n_hidden_encoder=int(cfg["n_hidden_encoder"]),
            n_layers_encoder=int(cfg["n_layers_encoder"]),
            n_hidden_decoder=int(cfg["n_hidden_decoder"]),
            n_layers_decoder=int(cfg["n_layers_decoder"]),
            dropout_rate_encoder=float(cfg["dropout_rate_encoder"]),
            dropout_rate_decoder=float(cfg["dropout_rate_decoder"]),
            seed=seed,
        )
        log(
            f"CPA n_genes={adata.n_vars} n_perts={len(model.pert_encoder)} "
            f"params={sum(p.numel() for p in model.module.parameters()):,}"
        )

        checkpoint = Path(run_dir) / "checkpoint"
        checkpoint.mkdir(parents=True, exist_ok=True)
        model.train(
            max_epochs=int(cfg["max_epochs"]),
            batch_size=int(cfg["batch_size"]),
            use_gpu=torch.cuda.is_available(),
            check_val_every_n_epoch=int(cfg["check_val_every_n_epoch"]),
            early_stopping_patience=int(cfg["early_stop_patience"]),
            plan_kwargs=dict(cfg["plan"]),
            save_path=str(checkpoint),
        )
        write_json(
            Path(run_dir) / "fit_history.json",
            {
                "smoke": _smoke(),
                "epoch_history_csv": model.epoch_history.to_csv(index=False),
            },
        )
        log(f"checkpoint {checkpoint}")
        return checkpoint

    def predict(self, cases, checkpoint, run_dir, controls=None):
        if controls is None:
            raise RuntimeError("CPA needs --controls (the shared DMSO means)")
        case = load_npz(cases)
        control_sets = Controls(controls)
        pred_adata = self._pred_adata(case, control_sets)
        model = cpa.CPA.load(
            dir_path=str(checkpoint),
            adata=pred_adata,
            use_gpu=torch.cuda.is_available(),
        )
        model.predict(pred_adata, batch_size=128, n_samples=1, return_mean=True)
        prediction = np.asarray(pred_adata.obsm["CPA_pred"], dtype=np.float32)
        basal = pred_adata.X
        if hasattr(basal, "toarray"):
            basal = basal.toarray()
        delta = prediction - np.asarray(basal, dtype=np.float32)
        log(f"predict {delta.shape} finite={np.isfinite(delta).all()}")
        return ensure_float32(delta)


if __name__ == "__main__":
    raise SystemExit(run_cli(CPAAdapter()))
