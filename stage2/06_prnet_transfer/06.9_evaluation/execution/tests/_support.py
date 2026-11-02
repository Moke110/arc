"""Shared synthetic support for the 06.9 offline tests (stdlib unittest only).

Nothing here touches the real frozen input, the real checkpoints or any GPU: the tests
build a tiny *contract-shaped* frozen directory (``manifest.json`` + ``locked_config.json``,
which is all ``build_units.py`` reads), fake candidate run directories with a valid
config-digest chain, and synthetic evaluator output directories whose
``aggregate.json`` is consistent with their ``per_condition.csv`` — exactly the invariants
``load_unit_eval`` checks. The synthetic fixture is a test aid, never evidence.
"""

from __future__ import annotations

import atexit
import hashlib
import json
import shutil
import sys
import tempfile
from pathlib import Path

EXECUTION = Path(__file__).resolve().parents[1]
if str(EXECUTION) not in sys.path:
    sys.path.insert(0, str(EXECUTION))

import common  # noqa: E402

SEEDS = [20261021, 20261022, 20261023]
LR = 0.001
N_FOLDS = 5
K = 200

_TMP_ROOT: Path | None = None


def tmp_root() -> Path:
    global _TMP_ROOT
    if _TMP_ROOT is None:
        _TMP_ROOT = Path(tempfile.mkdtemp(prefix="t09-evaluation-"))
        atexit.register(shutil.rmtree, _TMP_ROOT, True)
    return _TMP_ROOT


def clean_dir(path: Path) -> Path:
    if Path(path).exists():
        shutil.rmtree(path)
    Path(path).mkdir(parents=True, exist_ok=True)
    return Path(path)


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")


def write_bytes(path: Path, payload: bytes = b"checkpoint-bytes") -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return hashlib.sha256(payload).hexdigest()


# --------------------------------------------------------------------------- #
# A contract-shaped frozen input (only the two identity files are read)
# --------------------------------------------------------------------------- #
ARTIFACTS = ("cm_transfer.h5ad", "tahoe_transfer.h5ad", "folds.parquet")


def make_data_dir(root: Path, artifacts=ARTIFACTS) -> dict:
    data_dir = clean_dir(Path(root) / "data")
    hashes: dict[str, str] = {}
    for name in artifacts:
        hashes[name] = write_bytes(data_dir / name, f"frozen-{name}".encode())
    write_json(data_dir / "manifest.json", {"artifacts": {name: digest for name, digest in hashes.items()}})
    write_json(
        data_dir / "locked_config.json",
        {"metric": {"name": "direction_auprc", "k": K}},
    )
    frozen = common.frozen_identity(data_dir)
    return {
        "dir": data_dir,
        "frozen": frozen,
        "data_identity": common.data_identity_sha256(frozen),
        "structure": common.pipeline_modules()["config"].read_locked_config(
            data_dir / "locked_config.json"
        )[0]["structure"],
    }


# --------------------------------------------------------------------------- #
# Fake candidate run directories with a valid config-digest chain
# --------------------------------------------------------------------------- #
def make_candidate_run(
    root: Path,
    *,
    fold: int,
    seed: int,
    init: str,
    data: dict,
    lr: float = LR,
    init_source: Path | None = None,
    init_source_sha: str | None = None,
    structure: dict | None = None,
    status: str = "completed",
    stop_reason: str = "early_stop_no_auprc_improvement_20",
) -> dict:
    run_dir = Path(root) / f"run_f{fold}_s{seed}_{lr}_{init}"
    clean_dir(run_dir)
    config_mod = common.pipeline_modules()["config"]
    structure = structure or dict(data["structure"])
    config = {
        "stage": "finetune",
        "fold": int(fold),
        "seed": int(seed),
        "structure": structure,
        "recipe": {"lr": float(lr), "batch_size": 512, "max_epochs": 100},
        "metric": {"name": "direction_auprc", "k": K, "plateau_mode": "max"},
        "loss": "gaussian_nll",
        "optimizer": "adam",
        "prediction_noise": "zero",
        "train_drugs": "main",
        "init": {"mode": "pretrained" if init == "tahoe" else "random"},
    }
    digest = config_mod.config_fingerprint(config)
    # mirror train.py: config.json is the resolved config plus the digest and run metadata
    config["config_sha256"] = digest
    config["run"] = {
        "unit_id": f"run_f{fold}_s{seed}_{lr}_{init}",
        "wall_s": 41.8,
        "peak_rss_mb": 1330.0,
        "torch": "2.0.1+cu118",
    }
    write_json(run_dir / "config.json", config)
    checkpoint_sha = write_bytes(
        run_dir / "checkpoint" / "model.pt", f"weights-{fold}-{seed}-{init}".encode()
    )
    init_block = (
        {
            "mode": "pretrained",
            "source_path": str(init_source) if init_source else "/nowhere/model.pt",
            "checkpoint_sha256": init_source_sha,
            "checkpoint_epoch": 12,
            "checkpoint_seed": 20261021,
            "parameter_count": 7185000,
        }
        if init == "tahoe"
        else {"mode": "random", "parameter_count": 7185000}
    )
    manifest = {
        "schema": "stage2/06.4_pipeline/run/v1",
        "kind": "run",
        "stage": "finetune",
        "status": status,
        "fold": int(fold),
        "seed": int(seed),
        "role": None,
        "config": config,
        "config_sha256": digest,
        "init": init_block,
        "structure": structure,
        "data": {
            "dir": str(data["dir"]),
            "artifacts": dict(data["frozen"]["artifacts"]),
            "gene_axis": {"n_genes": 18164, "measured": 18164},
        },
        "checkpoint": {"path": "checkpoint/model.pt", "sha256": checkpoint_sha},
        "best_epoch": 7,
        "best_valid_auprc": 0.09,
        "epochs_run": 28,
        "stop_reason": stop_reason,
        "timings": {"wall_s": 41.8},
        "peak_rss_mb": 1330.0,
    }
    write_json(run_dir / "manifest.json", manifest)
    write_json(
        run_dir / "fit_history.json",
        {"history": [{"epoch": 1, "train_loss": 1.0, "valid_loss": 1.1, "valid_auprc": 0.05}]},
    )
    return {
        "run_dir": run_dir,
        "checkpoint_sha256": checkpoint_sha,
        "config_sha256": digest,
        "checkpoint_bytes": (run_dir / "checkpoint" / "model.pt").stat().st_size,
        "best_valid_auprc": 0.09,
        "best_epoch": 7,
        "epochs_run": 28,
        "stop_reason": stop_reason,
    }


def make_selection_csv(path: Path, runs: list[dict], data: dict, init: str) -> Path:
    import pandas as pd

    rows = []
    for run in runs:
        rows.append(
            {
                "fold": int(run["fold"]),
                "seed": int(run["seed"]),
                "seed_index": SEEDS.index(int(run["seed"])),
                "lr_selected": float(run["lr"]),
                "unit_id": f"cm_f{run['fold']}_s{run['seed']}_lr{common.format_lr_tag(run['lr'])}_{init}",
                "config_sha256": run["result"]["config_sha256"],
                "checkpoint_worktree_relpath": str(
                    run["result"]["run_dir"] / "checkpoint" / "model.pt"
                ),
                "checkpoint_sha256": run["result"]["checkpoint_sha256"],
                "checkpoint_bytes": run["result"]["checkpoint_bytes"],
                "best_valid_auprc": run["result"]["best_valid_auprc"],
                "best_epoch": run["result"]["best_epoch"],
                "epochs_run": run["result"]["epochs_run"],
                "stop_reason": run["result"]["stop_reason"],
                "data_identity": data["data_identity"],
            }
        )
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def make_selection_inputs(root: Path, *, inits=("random", "tahoe"), lr_by_fold=None) -> dict:
    """Build a frozen dir + candidate runs + both selection CSVs."""
    root = Path(root)
    data = make_data_dir(root)
    pretrain_ckpt = root / "pretrain" / "checkpoint" / "model.pt"
    pretrain_sha = write_bytes(pretrain_ckpt, b"pretrain-weights")
    lr_by_fold = lr_by_fold or {fold: LR for fold in range(N_FOLDS)}
    outputs = {"data": data}
    for init in inits:
        runs = []
        for fold in range(N_FOLDS):
            for seed in SEEDS:
                lr = float(lr_by_fold[fold])
                result = make_candidate_run(
                    root / "runs" / init,
                    fold=fold, seed=seed, init=init, data=data, lr=lr,
                    init_source=pretrain_ckpt if init == "tahoe" else None,
                    init_source_sha=pretrain_sha if init == "tahoe" else None,
                )
                runs.append({"fold": fold, "seed": seed, "lr": lr, "result": result})
        outputs[f"{init}_record"] = make_selection_csv(
            root / "records" / f"{init}_selected_checkpoints.csv", runs, data, init
        )
        outputs[f"{init}_runs"] = runs
    outputs["pretrain_checkpoint"] = pretrain_ckpt
    return outputs


# --------------------------------------------------------------------------- #
# Synthetic evaluator output directories
# --------------------------------------------------------------------------- #
def make_unit(
    *, fold: int, seed: int, init: str, lr: float = LR, data_identity: str = "d" * 64
) -> dict:
    return {
        "unit_id": f"e9_f{fold}_s{seed}_lr{common.format_lr_tag(lr)}_{init}",
        "init_name": init,
        "fold": fold,
        "seed": seed,
        "seed_index": SEEDS.index(seed),
        "lr": lr,
        "lr_tag": common.format_lr_tag(lr),
        "checkpoint": {"sha256": f"{fold}{seed}{init}"[:8].ljust(64, "0"), "bytes": 10},
        "config_sha256": "c" * 64,
        "data_identity": data_identity,
        "best_valid_auprc": 0.09,
        "run_dir": f"/tmp/run_{fold}_{seed}_{init}",
        "init": {"mode": "pretrained" if init == "tahoe" else "random"},
        "training": {"best_epoch": 7, "epochs_run": 28, "stop_reason": "early_stop", "train_wall_s": 40.0, "train_peak_rss_mb": 1300.0},
    }


N_GENES = 400


def condition_metrics(
    *, fold: int, seed: int, init: str, drugs: list[str], doses=(0.2, 1.0, 10.0)
) -> dict:
    """Deterministic AUPRC values with a small Tahoe advantage (fold/seed dependent)."""
    advantage = 0.01 if init == "tahoe" else 0.0
    values: dict[str, float] = {}
    for drug_index, drug in enumerate(drugs):
        for dose_index, dose in enumerate(doses):
            condition = f"{drug}|{dose:g}"
            base = 0.08 + 0.01 * ((drug_index + dose_index + fold + seed_index(seed)) % 4)
            values[condition] = base + advantage
    return values


def seed_index(seed: int) -> int:
    return SEEDS.index(int(seed)) if int(seed) in SEEDS else 0


def write_unit_eval(
    out: Path,
    unit: dict,
    metric_values: dict[str, float],
    *,
    drugs: list[str],
    doses=(0.2, 1.0, 10.0),
    rmse: float = 0.14,
    true_positive_genes: int = 250,
    data_artifacts: dict | None = None,
) -> Path:
    """Write an evaluator-shaped ``eval/`` directory consistent with its aggregate."""
    import numpy as np
    import pandas as pd

    unit_dir = Path(out) / "units" / unit["unit_id"]
    eval_dir = unit_dir / "eval"
    clean_dir(eval_dir)
    conditions = sorted(metric_values)
    n_conditions = len(conditions)
    frame_rows = []
    for index, condition in enumerate(conditions):
        frame_rows.append(
            {
                "condition": condition,
                "n_replicates": 1,
                f"auprc_{K}": metric_values[condition],
                f"auprc_up_{K}": metric_values[condition],
                f"auprc_down_{K}": metric_values[condition],
                f"n_up_{K}": K,
                f"n_down_{K}": 50,
                "rmse": rmse + 0.001 * index,
                "pcc": 0.25,
                "spearman": 0.2,
                "pred_constant": False,
                "true_constant": False,
                "regression_nan_reason": "",
            }
        )
    frame = pd.DataFrame(frame_rows)
    frame.to_csv(eval_dir / "per_condition.csv", index=False, float_format="%.10g")

    true_delta = np.zeros((n_conditions, N_GENES), dtype=np.float32)
    true_delta[:, :true_positive_genes] = 1.0
    true_delta[:, true_positive_genes : true_positive_genes + 100] = -1.0
    pred_delta = np.zeros_like(true_delta)
    np.savez_compressed(
        eval_dir / "deltas.npz",
        conditions=np.asarray(conditions, dtype=str),
        true_delta=true_delta,
        pred_delta=pred_delta,
        condition_ids=np.arange(n_conditions, dtype=np.int32),
        rows=np.arange(n_conditions, dtype=np.int64),
        k=np.asarray([100, 200, 500], dtype=np.int32),
    )
    macro_values = np.asarray(list(metric_values.values()), dtype=np.float64)
    aggregate = {
        "stage": "finetune",
        "role": "test",
        "fold": int(unit["fold"]),
        "n_conditions": n_conditions,
        "n_rows": n_conditions,
        "k_list": [100, 200, 500],
        "headline_k": K,
        "by_k": {
            str(K): {
                "macro": {"mean": float(macro_values.mean()), "sd": float(macro_values.std(ddof=1)), "nan": 0, "n": int(macro_values.size)},
                "up": {"mean": float(macro_values.mean()), "sd": 0.0, "nan": 0, "n": int(macro_values.size)},
                "down": {"mean": float(macro_values.mean()), "sd": 0.0, "nan": 0, "n": int(macro_values.size)},
            }
        },
        "regression": {
            "rmse": {"mean": rmse, "sd": 0.01, "nan": 0, "n": n_conditions},
            "pcc": {"mean": 0.25, "sd": 0.01, "nan": 0, "n": n_conditions},
            "spearman": {"mean": 0.2, "sd": 0.01, "nan": 0, "n": n_conditions},
        },
    }
    write_json(eval_dir / "aggregate.json", aggregate)
    write_json(
        eval_dir / "manifest.json",
        {
            "kind": "evaluate",
            "stage": "finetune",
            "fold": int(unit["fold"]),
            "role": "test",
            "checkpoint": {"sha256": unit["checkpoint"]["sha256"]},
            "data": {"artifacts": dict(data_artifacts or {})},
            "predictions": {"pred_expr_sha256": "p" * 64, "rows": n_conditions},
            "truth": {"source": "synthetic"},
        },
    )
    return eval_dir


def write_synthetic_units_json(
    path: Path,
    out: Path,
    *,
    data: dict,
    units: list[dict],
    phase: str = "test",
    complete: bool = True,
    fold_reference: dict | None = None,
    created_utc: str = "2026-01-01T00:00:00Z",
) -> Path:
    payload = {
        "schema": "stage2/06.9_evaluation/units/v1",
        "kind": "evaluation_units",
        "created_utc": created_utc,
        "ticket": "#9",
        "frozen": data["frozen"],
        "data_identity": data["data_identity"],
        "locked_config": {
            "structure": data["structure"], "n_folds": N_FOLDS, "seeds": SEEDS,
            "cm_lrs": [1e-4, 3e-4, 1e-3], "metric_k": K,
        },
        "requirements": {
            "n_folds": N_FOLDS, "seeds": SEEDS, "lrs": [1e-4, 3e-4, 1e-3],
            "inits": ["random", "tahoe"], "n_units": 30, "k": K, "complete": complete,
        },
        "complete": complete,
        "missing_inits": [] if complete else ["tahoe"],
        "sources": {
            "random_record": {"label": "random-init (#6)", "sha256": "a" * 64},
            "tahoe_record": {"label": "Tahoe-init (#7)", "sha256": "b" * 64} if complete else None,
            "baselines": None,
        },
        "fold_reference": fold_reference or {},
        "gates": {
            "unit_count_ok": True,
            "checkpoint_sha256_verified": True,
            "checkpoint_bytes_match": True,
            "config_digest_chain_ok": True,
            "data_identity_uniform": True,
            "init_provenance_ok": True,
            "no_test_artifacts_in_candidate_runs": True,
            "three_seeds_per_fold_init": True,
            "one_lr_per_fold_init": True,
        },
        "units": units,
    }
    write_json(path, payload)
    return path


def write_synthetic_population(
    root: Path,
    *,
    drugs_per_fold: int = 4,
    doses=(0.2, 1.0, 10.0),
    out: Path | None = None,
    advantage: float = 0.0,
    folds=(0, 1, 2, 3, 4),
    seeds=tuple(SEEDS),
    inits=("random", "tahoe"),
    phase: str = "test",
) -> dict:
    """Write 30 evaluator-shaped unit dirs + a units.json + score records."""
    out = Path(out) if out is not None else Path(root) / "out"
    data = make_data_dir(Path(root))
    units: list[dict] = []
    drugs_by_fold: dict[int, list[str]] = {}
    for fold in folds:
        drugs = [f"DRUG{fold:02d}{index:02d}" for index in range(drugs_per_fold)]
        drugs_by_fold[fold] = drugs
        for seed in seeds:
            for init in inits:
                unit = make_unit(fold=fold, seed=seed, init=init, data_identity=data["data_identity"])
                values = {}
                for drug_index, drug in enumerate(drugs):
                    for dose_index, dose in enumerate(doses):
                        values[f"{drug}|{dose:g}"] = (
                            0.10
                            + 0.01 * ((drug_index + dose_index + fold + seed_index(seed)) % 4)
                            + (advantage if init == "tahoe" else 0.0)
                        )
                write_unit_eval(
                    out, unit, values, drugs=drugs, doses=doses,
                    data_artifacts=data["frozen"]["artifacts"],
                )
                # per-unit score record for budgets + status
                write_json(
                    out / "units" / unit["unit_id"] / "record.json",
                    {
                        "unit_id": unit["unit_id"], "phase": "test", "status": "completed",
                        "fold": fold, "seed": seed, "init": init, "lr": LR,
                        "checkpoint_sha256": unit["checkpoint"]["sha256"],
                        "steps": {
                            "predict": {"exit_code": 0, "wall_s": 2.0, "peak_rss_mb": 1200.0},
                            "evaluate": {"exit_code": 0, "wall_s": 3.0, "peak_rss_mb": 1100.0},
                        },
                        "gpu": {"peak_mb": 2500},
                    },
                )
                units.append(unit)
    fold_reference = {
        str(fold): {
            "n_valid_rows": drugs_per_fold * len(doses),
            "n_test_rows": drugs_per_fold * len(doses),
            "n_test_conditions": drugs_per_fold * len(doses),
            "n_test_drugs": drugs_per_fold,
            "test_row_sha256": "t" * 64,
            "gene_axis_sha256": "g" * 64,
        }
        for fold in folds
    }
    units_path = write_synthetic_units_json(
        out / "record" / "units.json", out, data=data, units=units, fold_reference=fold_reference
    )
    # score summary + runs (gates use the runs file to prove the ordering)
    write_json(
        out / "record" / "score_summary.json",
        {
            "schema": "stage2/06.9_evaluation/score/v1/summary", "phase": phase,
            "units_planned": len(units), "units_completed": len(units), "units_failed": 0,
            "failures": [], "wall_s": 150.0, "host_peak_rss_mb": 1300.0,
            "gpu_peak_mb": 2600, "started_utc": "2026-01-02T00:00:00Z",
        },
    )
    runs = "unit_id,phase,status,fold,seed,init,lr,predict_exit,predict_wall_s,predict_peak_rss_mb,evaluate_exit,evaluate_wall_s,evaluate_peak_rss_mb,gpu_peak_mb,eval_manifest_sha256,error\n"
    runs += "\n".join(
        f"{unit['unit_id']},test,completed,{unit['fold']},{unit['seed']},{unit['init_name']},"
        f"{LR},0,2.0,1200.0,0,3.0,1100.0,2500,"
        for unit in units
    )
    (out / "record" / "score_runs.csv").write_text(runs, encoding="utf-8")
    return {
        "out": out, "units": units, "units_path": units_path, "data": data,
        "drugs_by_fold": drugs_by_fold,
        "fold_reference": fold_reference,
    }


def write_synthetic_baselines(root: Path, *, drugs_by_fold: dict, doses=(0.2, 1.0, 10.0)) -> Path:
    """#8-shaped per_drug_dose.csv + summary_table.csv + fold_coverage.csv."""
    import pandas as pd

    baselines = Path(root) / "baselines"
    baselines.mkdir(parents=True, exist_ok=True)
    rows = []
    for fold, drugs in drugs_by_fold.items():
        for drug in drugs:
            for dose in doses:
                for name, value in (("no_change", 0.011), ("delta_mean", 0.15), ("ridge_fp", 0.13), ("delta_mean_dose", 0.14)):
                    rows.append(
                        {
                            "fold": fold, "baseline": name, "drug_key": drug,
                            "dose_um": dose, "condition": f"{drug}|{dose:g}",
                            "n_replicates": 1, f"auprc_{K}": value,
                            f"auprc_up_{K}": value, f"auprc_down_{K}": value,
                            f"n_up_{K}": K, f"n_down_{K}": 50,
                            "rmse": 0.14, "pcc": 0.24, "spearman": 0.19,
                            "pred_constant": name == "no_change",
                            "true_constant": False, "regression_nan_reason": "",
                        }
                    )
    pd.DataFrame(rows).to_csv(baselines / "per_drug_dose.csv", index=False)
    summary_rows = []
    for name, value, role in (
        ("no_change", 0.011, "primary"),
        ("delta_mean", 0.15, "primary"),
        ("delta_mean_dose", 0.14, "auxiliary"),
        ("ridge_fp", 0.13, "auxiliary"),
    ):
        summary_rows.append(
            {
                "scope": "macro_over_folds", "fold": "ALL", "baseline": name, "role": role,
                "auprc": value, "auprc_up": value, "auprc_down": value, "rmse": 0.14,
                "pcc": 0.24, "spearman": 0.19,
            }
        )
    pd.DataFrame(summary_rows).to_csv(baselines / "summary_table.csv", index=False)
    coverage = []
    for fold, drugs in drugs_by_fold.items():
        coverage.append(
            {
                "fold": fold, "n_train_rows": 100, "n_train_conditions": 30,
                "n_train_drugs": 10, "n_valid_rows": len(drugs) * len(doses),
                "n_test_rows": len(drugs) * len(doses),
                "n_test_conditions": len(drugs) * len(doses), "n_test_drugs": len(drugs),
                "n_test_conditions_seen_in_train": 0, "n_plates_train": 3,
                "n_scored_conditions": len(drugs) * len(doses),
                "n_scored_replicates": len(drugs) * len(doses),
                "scored_conditions_unique": True, "test_rows_disjoint_from_train": True,
                "test_row_sha256": "t" * 64, "gene_axis_sha256": "g" * 64,
            }
        )
    pd.DataFrame(coverage).to_csv(baselines / "fold_coverage.csv", index=False)
    return baselines
