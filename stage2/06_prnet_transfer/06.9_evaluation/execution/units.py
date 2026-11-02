"""The 30-unit test manifest for ticket #9.

Ticket #9 evaluates exactly the checkpoints #6 and #7 **selected** — 5 folds × 3 seeds ×
2 initialisations — after every selection is frozen. This module turns the two
``selected_checkpoints.csv`` records into one resolved unit manifest and verifies, before
any test row is read:

* the checkpoint file exists and its sha256/bytes match the record;
* the run's resolved config reproduces the recorded ``config_sha256`` (three-way chain:
  ``config.json`` -> embedded digest -> candidate manifest -> selected record);
* every unit's candidate manifest used the frozen bytes (recorded ``data_identity`` equals
  the identity recomputed from the frozen input);
* the split of the 15 selected checkpoints is exactly 5 folds × 3 seeds per init, with one
  common LR per (fold, init) — the rule #6/#7 locked;
* no candidate run directory holds a prediction/evaluation artifact (selection never read
  test expression);
* Tahoe-init units carry a strict warm-start provenance (``init.mode == "pretrained"`` with
  a recorded source checkpoint sha256), random-init units carry ``init.mode == "random"``.

The manifest is *evidence*: it is written once, then read by ``score.py``/``write_record.py``
and never recomputed from test results.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd

import common

__all__ = [
    "SCHEMA",
    "UnitError",
    "SelectionRecord",
    "read_selection_record",
    "read_candidate_run",
    "build_unit_manifest",
    "expected_fold_counts",
    "load_unit_manifest",
]

SCHEMA = "stage2/06.9_evaluation/units/v1"

_LR_COLUMNS = ("lr_selected", "lr")
_PATH_COLUMNS = ("checkpoint_worktree_relpath", "checkpoint_path")
_FORBIDDEN_RUN_DIRS = ("prediction", "predictions", "pred", "evaluation", "eval", "test")


class UnitError(common.EvaluationError):
    """The unit manifest cannot be built or verified."""


@dataclass
class SelectionRecord:
    """One init branch's frozen selection record (``selected_checkpoints.csv``)."""

    label: str
    init: str
    path: Path
    sha256: str
    rows: list[dict[str, Any]]
    data_identity: str
    lr_by_fold: dict[int, float]
    seeds: list[int]
    source: str = ""

    @property
    def n_rows(self) -> int:
        return len(self.rows)

    def record(self) -> dict[str, Any]:
        return {
            "label": self.label,
            "init": self.init,
            "source": self.source,
            "path": common.relpath(self.path),
            "sha256": self.sha256,
            "n_rows": self.n_rows,
            "data_identity": self.data_identity,
            "lr_by_fold": {str(k): v for k, v in sorted(self.lr_by_fold.items())},
            "seeds": list(self.seeds),
        }


def _require_columns(frame: pd.DataFrame, path: Path, required: Sequence[str]) -> None:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise UnitError(
            f"selection record {path} is missing required columns {missing}; "
            f"present: {list(frame.columns)}"
        )


def read_selection_record(
    path: Path | str,
    *,
    label: str,
    init: str,
    expected_folds: int = 5,
    expected_seeds: Sequence[int] | None = None,
) -> SelectionRecord:
    """Read + cheaply validate one ``selected_checkpoints.csv`` (no checkpoint IO yet)."""
    path = Path(path)
    if not path.exists():
        raise UnitError(
            f"{label} selection record is missing: {path}\n"
            f"(ticket #9 evaluates only frozen selections; build it in the {label} substep "
            "first, or pass --allow-incomplete for a harness dry run)"
        )
    frame = common.read_table(
        path,
        numeric_columns=(
            "fold", "seed", "seed_index", "lr_selected", "lr",
            "checkpoint_bytes", "best_valid_auprc", "best_epoch", "epochs_run",
        ),
        string_columns=("unit_id", "config_sha256", "checkpoint_sha256",
                        "checkpoint_worktree_relpath", "checkpoint_path",
                        "stop_reason", "data_identity"),
    )
    if frame.empty:
        raise UnitError(f"{label} selection record is empty: {path}")
    _require_columns(
        frame,
        path,
        ("fold", "seed", "config_sha256", "checkpoint_sha256", "best_valid_auprc"),
    )
    lr_column = next((c for c in _LR_COLUMNS if c in frame.columns), None)
    path_column = next((c for c in _PATH_COLUMNS if c in frame.columns), None)
    if lr_column is None:
        raise UnitError(f"{path}: no learning-rate column in {_LR_COLUMNS}")
    if path_column is None:
        raise UnitError(f"{path}: no checkpoint-path column in {_PATH_COLUMNS}")

    folds = sorted({int(value) for value in frame["fold"]})
    if folds != list(range(int(expected_folds))):
        raise UnitError(
            f"{label} selection record covers folds {folds}, expected "
            f"{list(range(int(expected_folds)))}"
        )
    rows: list[dict[str, Any]] = []
    lr_by_fold: dict[int, float] = {}
    seen: set[tuple[int, int]] = set()
    identities: set[str] = set()
    for _, raw in frame.iterrows():
        row = raw.to_dict()
        fold = int(row["fold"])
        seed = int(row["seed"])
        if (fold, seed) in seen:
            raise UnitError(f"{label}: duplicate (fold={fold}, seed={seed}) in {path}")
        seen.add((fold, seed))
        lr = float(row[lr_column])
        if not math.isfinite(lr) or lr <= 0:
            raise UnitError(f"{label}: fold {fold} seed {seed} has invalid lr {lr!r}")
        if fold in lr_by_fold and abs(lr_by_fold[fold] - lr) > 0:
            raise UnitError(
                f"{label}: fold {fold} has more than one selected LR "
                f"({lr_by_fold[fold]} vs {lr}); the common-LR rule is violated"
            )
        lr_by_fold[fold] = lr
        identity = str(row.get("data_identity", "") or "")
        if identity:
            identities.add(identity)
        row["_lr"] = lr
        row["_checkpoint_relpath"] = str(row[path_column])
        row["_record_line"] = int(len(rows))
        rows.append(row)
    if len(identities) > 1:
        raise UnitError(
            f"{label}: the selected checkpoints record {len(identities)} different data "
            f"identities; all units must share one frozen input"
        )
    seeds = sorted({int(row["seed"]) for row in rows})
    if expected_seeds is not None and seeds != sorted(int(s) for s in expected_seeds):
        raise UnitError(
            f"{label}: seeds {seeds} do not match the frozen locked config "
            f"{sorted(int(s) for s in expected_seeds)}"
        )
    per_fold = {
        fold: sum(1 for row in rows if int(row["fold"]) == fold) for fold in folds
    }
    bad = {fold: count for fold, count in per_fold.items() if count != len(seeds)}
    if bad:
        raise UnitError(
            f"{label}: every fold must keep all {len(seeds)} seeds; counts {per_fold}"
        )
    rows.sort(key=lambda item: (int(item["fold"]), int(item["seed"])))
    return SelectionRecord(
        label=label,
        init=init,
        path=path,
        sha256=common.sha256_file(path),
        rows=rows,
        data_identity=next(iter(identities), ""),
        lr_by_fold=lr_by_fold,
        seeds=seeds,
    )


def _read_json(path: Path) -> Any:
    return common.load_json(path)


def read_candidate_run(
    run_dir: Path | str,
    *,
    fold: int,
    seed: int,
    lr: float,
    config_sha256: str,
    checkpoint_path: Path,
    checkpoint_sha256: str,
    checkpoint_bytes: int | None,
    init: str,
    expected_structure: Mapping[str, Any] | None = None,
    expected_k: int | None = None,
    frozen_artifacts: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Verify one candidate run directory against its selection row.

    Returns the provenance block stored in the unit manifest; raises on any mismatch.
    """
    run_dir = Path(run_dir)
    manifest_file = run_dir / "manifest.json"
    config_file = run_dir / "config.json"
    if not manifest_file.exists():
        raise UnitError(f"candidate run has no manifest.json: {run_dir}")
    if not config_file.exists():
        raise UnitError(f"candidate run has no config.json: {run_dir}")
    manifest = _read_json(manifest_file)
    config = _read_json(config_file)
    modules = common.pipeline_modules()
    fingerprint = modules["config"].config_fingerprint

    if manifest.get("status") != "completed":
        raise UnitError(
            f"{run_dir}: candidate manifest status is {manifest.get('status')!r}, "
            "not 'completed'"
        )
    if manifest.get("stage") != "finetune":
        raise UnitError(f"{run_dir}: candidate stage is {manifest.get('stage')!r}")
    if int(manifest.get("fold", -1)) != int(fold):
        raise UnitError(f"{run_dir}: candidate fold {manifest.get('fold')} != {fold}")
    if int(manifest.get("seed", -1)) != int(seed):
        raise UnitError(f"{run_dir}: candidate seed {manifest.get('seed')} != {seed}")
    if str(manifest.get("role") or ""):
        raise UnitError(f"{run_dir}: a finetune candidate manifest records role {manifest.get('role')!r}")

    # config digest chain: config.json -> its embedded digest -> candidate manifest -> record
    digest = fingerprint(config)
    embedded = config.get("config_sha256")
    if embedded and embedded != digest:
        raise UnitError(f"{run_dir}: config.json embedded digest disagrees with its content")
    if manifest.get("config_sha256") != digest:
        raise UnitError(
            f"{run_dir}: candidate manifest config_sha256 {manifest.get('config_sha256')} "
            f"!= recomputed {digest}"
        )
    # train.py writes config.json as the resolved config *plus* `config_sha256` and a
    # `run` metadata block; config_fingerprint() hashes only the canonical config keys,
    # so the digest is stable, but the manifest's own config block must agree too.
    manifest_config = manifest.get("config") or {}
    if manifest_config and fingerprint(manifest_config) != digest:
        raise UnitError(
            f"{run_dir}: manifest.config does not reproduce the recorded config digest "
            f"{digest}"
        )
    if str(config_sha256) != digest:
        raise UnitError(
            f"selected record config_sha256 {config_sha256} != recomputed {digest} "
            f"for {run_dir}"
        )
    if int(config.get("fold", -1)) != int(fold) or int(config.get("seed", -1)) != int(seed):
        raise UnitError(f"{run_dir}: config fold/seed disagree with the record")
    run_recipe = (config.get("recipe") or {})
    if "lr" in run_recipe and abs(float(run_recipe["lr"]) - float(lr)) > 0:
        raise UnitError(
            f"{run_dir}: config recipe lr {run_recipe['lr']} != selected lr {lr}"
        )
    if str(config.get("train_drugs", "main")) != "main":
        raise UnitError(f"{run_dir}: candidate trained with train_drugs={config.get('train_drugs')!r}")
    if expected_structure is not None:
        if dict(config.get("structure") or {}) != dict(expected_structure):
            raise UnitError(f"{run_dir}: candidate structure disagrees with the locked config")
    metric = config.get("metric") or {}
    if expected_k is not None and int(metric.get("k", -1)) != int(expected_k):
        raise UnitError(f"{run_dir}: candidate metric.k {metric.get('k')} != {expected_k}")

    # checkpoint bytes
    if not checkpoint_path.exists():
        raise UnitError(f"selected checkpoint file is missing: {checkpoint_path}")
    actual_sha = common.sha256_file(checkpoint_path)
    if actual_sha != str(checkpoint_sha256):
        raise UnitError(
            f"checkpoint sha256 mismatch for {checkpoint_path}: record "
            f"{checkpoint_sha256}, file {actual_sha}"
        )
    actual_bytes = checkpoint_path.stat().st_size
    if checkpoint_bytes not in (None, "") and int(checkpoint_bytes) != actual_bytes:
        raise UnitError(
            f"checkpoint size mismatch for {checkpoint_path}: record {checkpoint_bytes}, "
            f"file {actual_bytes}"
        )

    # frozen input identity
    if frozen_artifacts is not None:
        run_artifacts = dict((manifest.get("data") or {}).get("artifacts") or {})
        if run_artifacts and run_artifacts != dict(frozen_artifacts):
            raise UnitError(
                f"{run_dir}: candidate manifest data artifacts disagree with the frozen input"
            )
    gene_axis = dict((manifest.get("data") or {}).get("gene_axis") or {})
    if gene_axis and int(gene_axis.get("n_genes", -1)) <= 0:
        raise UnitError(f"{run_dir}: candidate manifest carries no usable gene axis")

    # initialisation provenance
    init_block = dict(manifest.get("init") or {})
    mode = str(init_block.get("mode", ""))
    if init == "random" and mode != "random":
        raise UnitError(f"{run_dir}: expected a random-init run, manifest says {mode!r}")
    if init == "tahoe":
        if mode != "pretrained":
            raise UnitError(f"{run_dir}: expected a Tahoe warm-start, manifest says {mode!r}")
        source = init_block.get("source_path")
        source_sha = init_block.get("checkpoint_sha256")
        if not source or not source_sha:
            raise UnitError(f"{run_dir}: Tahoe-init provenance is incomplete: {init_block}")
        source_path = Path(str(source))
        init_verified = False
        if source_path.exists():
            actual = common.sha256_file(source_path)
            if actual != str(source_sha):
                raise UnitError(
                    f"{run_dir}: warm-start checkpoint sha256 {actual} != recorded {source_sha}"
                )
            init_verified = True
        init_provenance = {
            "mode": mode,
            "source_path": str(source),
            "checkpoint_sha256": str(source_sha),
            "checkpoint_epoch": init_block.get("checkpoint_epoch"),
            "checkpoint_seed": init_block.get("checkpoint_seed"),
            "source_file_present": source_path.exists(),
            "source_sha256_verified": init_verified,
            "parameter_count": init_block.get("parameter_count"),
        }
    else:
        init_provenance = {
            "mode": mode,
            "parameter_count": init_block.get("parameter_count"),
        }

    # selection never read test expression: no prediction/evaluation tree in the run dir
    forbidden = [
        str(entry.relative_to(run_dir))
        for entry in run_dir.iterdir()
        if entry.name.lower() in _FORBIDDEN_RUN_DIRS
    ]
    if forbidden:
        raise UnitError(
            f"{run_dir}: candidate run holds prediction/evaluation artifacts {forbidden}; "
            "a selected checkpoint must never have read the test split"
        )
    if not (run_dir / "fit_history.json").exists():
        raise UnitError(f"{run_dir}: candidate run has no fit_history.json (no validation curve)")
    fit_history_sha = common.sha256_file(run_dir / "fit_history.json")
    history = _read_json(run_dir / "fit_history.json")
    entries = history.get("history") if isinstance(history, dict) else None
    if not entries:
        raise UnitError(f"{run_dir}: fit_history.json holds no epochs")

    return {
        "run_dir": str(run_dir),
        "run_dir_worktree_relpath": _worktree_relpath(run_dir),
        "candidate_manifest_sha256": common.sha256_file(manifest_file),
        "candidate_config_sha256": common.sha256_file(config_file),
        "fit_history_sha256": fit_history_sha,
        "candidate_status": str(manifest.get("status")),
        "init": init_provenance,
        "best_epoch": manifest.get("best_epoch"),
        "epochs_run": manifest.get("epochs_run"),
        "stop_reason": manifest.get("stop_reason"),
        "train_wall_s": (manifest.get("timings") or {}).get("wall_s"),
        "train_peak_rss_mb": manifest.get("peak_rss_mb"),
        "structure": dict(config.get("structure") or {}),
        "n_genes": int((config.get("structure") or {}).get("n_genes", -1)),
        "parameter_count": init_block.get("parameter_count"),
    }


def _worktree_relpath(path: Path) -> str:
    """Path relative to whichever repository root contains it (record-portable)."""
    resolved = path.resolve()
    for root in common.repo_roots():
        try:
            return str(resolved.relative_to(root))
        except ValueError:
            continue
    return str(resolved)


def expected_fold_counts(baselines_dir: Path | str | None) -> dict[str, Any]:
    """Per-fold reference sizes/hashes from #8's record (the identical test population)."""
    if baselines_dir is None:
        return {}
    coverage = Path(baselines_dir) / "fold_coverage.csv"
    if not coverage.exists():
        return {}
    frame = common.read_table(
        coverage,
        numeric_columns=(
            "fold", "n_valid_rows", "n_test_rows", "n_test_conditions",
            "n_test_drugs",
        ),
        string_columns=("test_row_sha256", "gene_axis_sha256", "c_star_sha256"),
    )
    reference: dict[str, Any] = {}
    for _, row in frame.iterrows():
        fold = int(row["fold"])
        reference[str(fold)] = {
            "n_valid_rows": int(row["n_valid_rows"]),
            "n_test_rows": int(row["n_test_rows"]),
            "n_test_conditions": int(row["n_test_conditions"]),
            "n_test_drugs": int(row["n_test_drugs"]),
            "test_row_sha256": str(row.get("test_row_sha256", "")),
            "gene_axis_sha256": str(row.get("gene_axis_sha256", "")),
            "c_star_sha256": str(row.get("c_star_sha256", "")),
        }
    reference["source"] = {
        "path": common.relpath(coverage),
        "sha256": common.sha256_file(coverage),
    }
    return reference


def _gates(
    units: Sequence[Mapping[str, Any]], requirements: Mapping[str, Any], records: Sequence[SelectionRecord]
) -> dict[str, Any]:
    by_key = {(int(u["fold"]), int(u["seed"]), str(u["init_name"])): u for u in units}
    expected = {
        (fold, seed, init)
        for fold in range(int(requirements["n_folds"]))
        for seed in requirements["seeds"]
        for init in requirements["inits"]
    }
    missing = sorted(expected - set(by_key))
    extra = sorted(set(by_key) - expected)
    seeds_per_fold_init = {
        f"{fold}:{init}": sorted(
            int(u["seed"]) for u in units
            if int(u["fold"]) == fold and u["init_name"] == init
        )
        for fold in range(int(requirements["n_folds"]))
        for init in requirements["inits"]
    }
    lr_per_fold_init = {
        f"{fold}:{init}": sorted(
            {float(u["lr"]) for u in units
             if int(u["fold"]) == fold and u["init_name"] == init}
        )
        for fold in range(int(requirements["n_folds"]))
        for init in requirements["inits"]
    }
    return {
        "unit_count": len(units),
        "unit_count_ok": len(units) == int(requirements["n_units"]),
        "no_missing_unit": not missing,
        "no_extra_unit": not extra,
        "missing_units": [list(item) for item in missing],
        "extra_units": [list(item) for item in extra],
        "seeds_per_fold_init": seeds_per_fold_init,
        "three_seeds_per_fold_init": all(
            value == sorted(int(s) for s in requirements["seeds"])
            for value in seeds_per_fold_init.values()
        ),
        "one_lr_per_fold_init": all(len(value) == 1 for value in lr_per_fold_init.values()),
        "lr_per_fold_init": lr_per_fold_init,
        "checkpoint_sha256_verified": all(u["checkpoint"]["sha256_verified"] for u in units),
        "checkpoint_bytes_match": all(u["checkpoint"]["bytes_match"] for u in units),
        "config_digest_chain_ok": all(u["verification"]["config_digest_chain"] for u in units),
        "data_identity_uniform": len({u["data_identity"] for u in units}) == 1,
        "init_provenance_ok": all(
            (u["init"]["mode"] == "random") if u["init_name"] == "random"
            else (u["init"]["mode"] == "pretrained" and bool(u["init"]["checkpoint_sha256"]))
            for u in units
        ),
        "warm_start_sources_verified": all(
            u["init"]["source_sha256_verified"]
            for u in units
            if u["init_name"] == "tahoe"
        ),
        "warm_start_sources_present": all(
            u["init"]["source_file_present"]
            for u in units
            if u["init_name"] == "tahoe"
        ),
        "no_test_artifacts_in_candidate_runs": all(
            u["verification"]["no_test_artifacts"] for u in units
        ),
        "records_complete": all(record.n_rows > 0 for record in records),
    }


def build_unit_manifest(
    *,
    data_dir: Path | str,
    random_record: Path | str,
    tahoe_record: Path | str | None = None,
    baselines_dir: Path | str | None = None,
    allow_incomplete: bool = False,
) -> dict[str, Any]:
    """Build + verify the unit manifest; fails closed unless it is complete (by default)."""
    data_dir = Path(data_dir)
    frozen = common.frozen_identity(data_dir)
    identity = common.data_identity_sha256(frozen)
    modules = common.pipeline_modules()
    locked, notes = modules["config"].read_locked_config(data_dir / "locked_config.json")
    n_folds = int(locked.get("n_folds", 5))
    seeds = [int(s) for s in locked.get("seeds", [])]
    lrs = [float(v) for v in (locked.get("lrs") or {}).get("cm", [])]
    structure = dict(locked.get("structure") or {})
    k = int((locked.get("metric") or {}).get("k", 200))
    if not seeds:
        raise UnitError("the frozen locked config records no seeds")

    random_rec = read_selection_record(
        random_record, label="random-init (#6)", init="random",
        expected_folds=n_folds, expected_seeds=seeds,
    )
    tahoe_rec: SelectionRecord | None = None
    if tahoe_record is not None and Path(tahoe_record).exists():
        tahoe_rec = read_selection_record(
            tahoe_record, label="Tahoe-init (#7)", init="tahoe",
            expected_folds=n_folds, expected_seeds=seeds,
        )
    elif not allow_incomplete:
        raise UnitError(
            f"Tahoe-init selection record is missing: {tahoe_record}\n"
            "ticket #9 cannot evaluate 30 units until #7 freezes its 15 checkpoints; "
            "the harness can be exercised with --allow-incomplete (random branch only)"
        )
    records = [rec for rec in (random_rec, tahoe_rec) if rec is not None]
    for record in records:
        if record.data_identity and record.data_identity != identity:
            raise UnitError(
                f"{record.label}: recorded data_identity {record.data_identity} != "
                f"recomputed frozen identity {identity}"
            )
        recorded_lrs = set(record.lr_by_fold.values())
        if not recorded_lrs.issubset(set(lrs)):
            raise UnitError(
                f"{record.label}: selected LRs {sorted(recorded_lrs)} are not a subset of "
                f"the frozen candidate LRs {lrs}"
            )

    units: list[dict[str, Any]] = []
    for record in records:
        for row in record.rows:
            fold = int(row["fold"])
            seed = int(row["seed"])
            lr = float(row["_lr"])
            relpath = row["_checkpoint_relpath"]
            checkpoint_path, root = common.resolve_checkpoint_path(relpath)
            run_dir = checkpoint_path.parent.parent
            provenance = read_candidate_run(
                run_dir,
                fold=fold,
                seed=seed,
                lr=lr,
                config_sha256=str(row["config_sha256"]),
                checkpoint_path=checkpoint_path,
                checkpoint_sha256=str(row["checkpoint_sha256"]),
                checkpoint_bytes=(
                    int(row["checkpoint_bytes"]) if "checkpoint_bytes" in row and not _blank(row["checkpoint_bytes"]) else None
                ),
                init=record.init,
                expected_structure=structure,
                expected_k=k,
                frozen_artifacts=frozen["artifacts"],
            )
            unit_id = (
                f"e9_f{fold}_s{seed}_lr{common.format_lr_tag(lr)}_{record.init}"
            )
            units.append(
                {
                    "unit_id": unit_id,
                    "init_name": record.init,
                    "fold": fold,
                    "seed": seed,
                    "seed_index": int(row.get("seed_index", record.seeds.index(seed))),
                    "lr": lr,
                    "lr_tag": common.format_lr_tag(lr),
                    "checkpoint": {
                        "path": str(checkpoint_path),
                        "worktree_relpath": _worktree_relpath(checkpoint_path),
                        "resolved_root": str(root),
                        "sha256": str(row["checkpoint_sha256"]),
                        "sha256_verified": True,
                        "bytes": checkpoint_path.stat().st_size,
                        "bytes_match": (
                            True
                            if "checkpoint_bytes" not in row or _blank(row["checkpoint_bytes"])
                            else int(row["checkpoint_bytes"]) == checkpoint_path.stat().st_size
                        ),
                        "record_sha256": str(row["checkpoint_sha256"]),
                    },
                    "config_sha256": str(row["config_sha256"]),
                    "data_identity": identity,
                    "best_valid_auprc": float(row["best_valid_auprc"]),
                    "run_dir": provenance["run_dir"],
                    "run_dir_worktree_relpath": provenance["run_dir_worktree_relpath"],
                    "init": provenance["init"],
                    "training": {
                        "best_epoch": provenance["best_epoch"],
                        "epochs_run": provenance["epochs_run"],
                        "stop_reason": provenance["stop_reason"],
                        "train_wall_s": provenance["train_wall_s"],
                        "train_peak_rss_mb": provenance["train_peak_rss_mb"],
                    },
                    "structure": provenance["structure"],
                    "n_genes": provenance["n_genes"],
                    "parameter_count": provenance["parameter_count"],
                    "record_label": record.label,
                    "record_sha256": record.sha256,
                    "verification": {
                        "config_digest_chain": True,
                        "no_test_artifacts": True,
                        "candidate_manifest_sha256": provenance["candidate_manifest_sha256"],
                        "fit_history_sha256": provenance["fit_history_sha256"],
                    },
                    "manifest_sha256": provenance["candidate_manifest_sha256"],
                    "fit_history_sha256": provenance["fit_history_sha256"],
                }
            )

    units.sort(key=lambda item: (item["init_name"] != "random", item["fold"], item["seed_index"]))
    requirements = {
        "n_folds": n_folds,
        "seeds": seeds,
        "lrs": lrs,
        "inits": ["random", "tahoe"],
        "n_units": n_folds * len(seeds) * 2,
        "k": k,
        "complete": tahoe_rec is not None,
    }
    gates = _gates(units, requirements, records)
    missing_inits = [
        init for init in ("random", "tahoe") if init not in {u["init_name"] for u in units}
    ]
    # `warm_start_sources_*` are recorded as warnings: the warm-start weights are
    # gitignored local files and their absence is reported, never silently accepted as
    # "verified". The unit is still evaluable because its own checkpoint is verified and
    # the candidate manifest records the strict-load provenance.
    soft_gates = {"warm_start_sources_verified", "warm_start_sources_present"}
    complete = tahoe_rec is not None and all(
        value
        for key, value in gates.items()
        if isinstance(value, bool) and key not in soft_gates
    )
    payload = {
        "schema": SCHEMA,
        "kind": "evaluation_units",
        "created_utc": common.now_utc(),
        "ticket": "#9 (Moke110/arc)",
        "frozen": frozen,
        "data_identity": identity,
        "locked_config": {
            "path": common.relpath(data_dir / "locked_config.json"),
            "sha256": frozen["locked_config_sha256"],
            "structure": structure,
            "n_folds": n_folds,
            "seeds": seeds,
            "cm_lrs": lrs,
            "metric_k": k,
            "mapping": notes.get("mapping", []),
            "defaulted": notes.get("defaulted", []),
        },
        "requirements": requirements,
        "complete": bool(complete),
        "missing_inits": missing_inits,
        "sources": {
            "random_record": random_rec.record(),
            "tahoe_record": tahoe_rec.record() if tahoe_rec is not None else None,
            "baselines": common.relpath(baselines_dir) if baselines_dir else None,
        },
        "fold_reference": expected_fold_counts(baselines_dir),
        "gates": gates,
        "units": units,
    }
    return payload


def _blank(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return str(value).strip() == ""


def load_unit_manifest(path: Path | str, *, require_complete: bool = False) -> dict[str, Any]:
    """Read a written unit manifest and validate its shape."""
    path = Path(path)
    if not path.exists():
        raise UnitError(
            f"unit manifest is missing: {path}\nrun `make units` (build_units.py) first"
        )
    payload = common.load_json(path)
    if payload.get("schema") != SCHEMA:
        raise UnitError(f"{path}: expected schema {SCHEMA}, got {payload.get('schema')!r}")
    units = payload.get("units")
    if not isinstance(units, list) or not units:
        raise UnitError(f"{path}: manifest holds no units")
    if require_complete and not payload.get("complete"):
        raise UnitError(
            f"{path}: the unit manifest is incomplete ({payload.get('missing_inits')}); "
            "the 30-unit test evaluation needs both init branches"
        )
    return payload


def fold_reference_for(payload: Mapping[str, Any], fold: int) -> dict[str, Any] | None:
    reference = payload.get("fold_reference") or {}
    entry = reference.get(str(fold))
    return dict(entry) if isinstance(entry, Mapping) else None
