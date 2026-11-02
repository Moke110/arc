"""Shared test support for ``06.8_baselines`` (stdlib unittest only).

The tests build the **synthetic frozen-input fixture of ticket #4**
(``06.4_pipeline/execution/tools/fixture.py``) because ticket #3's real frozen
input may legitimately be absent while Phase A is under construction.  The
fixture is a test aid, never evidence: it writes the exact file names / column
semantics the pipeline loader (and therefore the baselines) expects, on the real
18,164-gene axis.

The evaluator itself is never re-implemented here — tests call
``ptlib.evaluate.run_evaluation`` and the ``baselines.py`` writer, i.e. the same
code #9 will use.
"""

from __future__ import annotations

import atexit
import hashlib
import json
import shutil
import sys
import tempfile
from pathlib import Path

EXECUTION = Path(__file__).resolve().parents[1]  # .../06.8_baselines/execution
for path in (str(EXECUTION),):
    if path not in sys.path:
        sys.path.insert(0, path)


def find_repo_root(start: Path | None = None) -> Path:
    path = Path(start).resolve() if start is not None else Path(__file__).resolve()
    if path.is_file():
        path = path.parent
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    raise RuntimeError(f"no .git at or above {path}")


ROOT = find_repo_root()
PIPELINE_EXECUTION = ROOT / "stage2" / "06_prnet_transfer" / "06.4_pipeline" / "execution"
if str(PIPELINE_EXECUTION) not in sys.path:
    sys.path.insert(0, str(PIPELINE_EXECUTION))

_TMP_ROOT: Path | None = None
_FIXTURE_DIR: Path | None = None


def tmp_root() -> Path:
    global _TMP_ROOT
    if _TMP_ROOT is None:
        _TMP_ROOT = Path(tempfile.mkdtemp(prefix="t08-baselines-"))
        atexit.register(shutil.rmtree, _TMP_ROOT, True)
    return _TMP_ROOT


def _purge_bytecode(root: Path) -> None:
    for cache in list(root.rglob("__pycache__")):
        shutil.rmtree(cache, ignore_errors=True)


def fixture_dir() -> Path:
    """Build (once per process) and return the synthetic frozen-input directory."""
    global _FIXTURE_DIR
    if _FIXTURE_DIR is not None:
        return _FIXTURE_DIR
    from tools.fixture import build_fixture

    target = tmp_root() / "fixture"
    build_fixture(target)
    _purge_bytecode(PIPELINE_EXECUTION)
    _FIXTURE_DIR = target
    return target


def copy_fixture_to(src: Path, dst: Path) -> Path:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    return dst


# --------------------------------------------------------------------------- #
# frozen-input helpers
# --------------------------------------------------------------------------- #
def open_frozen(data_dir: Path | None = None):
    """Open a fixture frozen input (hash-verified) with a per-process cache."""
    from ptlib import data as data_mod

    target = Path(data_dir) if data_dir is not None else fixture_dir()
    key = str(target)
    cached = _FROZEN_CACHE.get(key)
    if cached is not None:
        return cached
    frozen = data_mod.open_frozen_input(target)
    _FROZEN_CACHE[key] = frozen
    return frozen


_FROZEN_CACHE: dict[str, object] = {}


def folds_table(frozen):
    from ptlib import data as data_mod

    return data_mod.load_folds(frozen)


def source_for(frozen, fold: int = 0):
    from ptlib import data as data_mod

    return data_mod.CMPredictSource(frozen, folds_table(frozen), int(fold), frozen.structures())


def truth_delta_for_rows(frozen, rows):
    """The evaluator's own truth Δ path (never a second implementation)."""
    import numpy as np

    from ptlib import evaluate as evaluate_mod

    treated, baseline, info = evaluate_mod._truth_for_rows(
        frozen, "finetune", np.asarray(rows, dtype=np.int64), None
    )
    return treated.astype(np.float64) - baseline.astype(np.float64), info


def condition_level_truth(frozen, source, rows) -> tuple:
    """Condition-mean truth Δ broadcast back to rows (condition order preserved)."""
    import numpy as np

    from ptlib import metrics as metrics_mod

    rows = np.asarray(rows, dtype=np.int64)
    row_truth, info = truth_delta_for_rows(frozen, rows)
    keys = np.asarray(source.condition_key(rows), dtype=object)
    unique, inverse = np.unique(keys, return_inverse=True)
    inverse = np.asarray(inverse, dtype=np.int64)
    condition_truth = metrics_mod.aggregate_rows(row_truth, inverse, int(unique.size))
    return condition_truth[inverse], condition_truth, unique, inverse, info


# --------------------------------------------------------------------------- #
# fixture mutation (leakage tests)
# --------------------------------------------------------------------------- #
def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rewrite_manifest(data_dir: Path) -> dict[str, str]:
    data_dir = Path(data_dir)
    manifest_file = data_dir / "manifest.json"
    payload = json.loads(manifest_file.read_text(encoding="utf-8"))
    names = sorted(payload.get("artifacts", {}).keys())
    hashes: dict[str, str] = {}
    for name in names:
        digest = _sha256(data_dir / name)
        hashes[name] = digest
        record = payload["artifacts"][name]
        if isinstance(record, dict):
            record["sha256"] = digest
        else:
            payload["artifacts"][name] = digest
    manifest_file.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return hashes


def poison_cm_targets(src: Path, dst: Path, rows, delta: float = 1.0) -> Path:
    """Copy *src* and shift the log1p targets (and raw counts) of *rows*."""
    import h5py
    import numpy as np

    dst = copy_fixture_to(Path(src), Path(dst))
    rows = np.sort(np.asarray(rows, dtype=np.int64))
    with h5py.File(dst / "cm_transfer.h5ad", "r+") as handle:
        for key in ("X", "layers/log1p_cp10k"):
            dataset = handle[key]
            for row in rows:
                dataset[row, :] = dataset[row, :].astype(np.float64) + delta
    rewrite_manifest(dst)
    return dst


def test_rows_for_fold(data_dir: Path, fold: int):
    """``(train, valid, test)`` treated row positions from metadata only."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    obs = ad.read_h5ad(Path(data_dir) / "cm_transfer.h5ad", backed="r").obs.copy()
    folds = pd.read_parquet(Path(data_dir) / "folds.parquet")
    role_map = dict(
        zip(
            folds[folds["fold"] == fold]["drug_key"].astype(str),
            folds[folds["fold"] == fold]["role"].astype(str),
        )
    )
    is_control = obs["is_control"].to_numpy(dtype=bool)
    qc = obs["qc_pass"].to_numpy(dtype=bool)
    drug = obs["drug_key"].astype(str).to_numpy()
    role_of_row = np.array([role_map.get(key, "") for key in drug], dtype=object)
    rows = {}
    for role in ("train", "valid", "test"):
        mask = (~is_control) & qc & (role_of_row == role)
        rows[role] = np.flatnonzero(mask)
    return rows["train"], rows["valid"], rows["test"]
