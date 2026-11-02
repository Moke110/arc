"""Shared test support: path setup, one synthetic fixture per process, helpers.

The fixture is the 18,164-gene synthetic frozen input from
``tools/fixture.py``.  It is built once per test process under a temp directory;
tests that mutate input artifacts copy it first.  The fixture is a **test aid**,
never evidence.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

EXECUTION = Path(__file__).resolve().parents[1]
if str(EXECUTION) not in sys.path:
    sys.path.insert(0, str(EXECUTION))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

_FIXTURE_DIR: Path | None = None
_TMP_ROOT: Path | None = None


def tmp_root() -> Path:
    global _TMP_ROOT
    if _TMP_ROOT is None:
        _TMP_ROOT = Path(tempfile.mkdtemp(prefix="ptlib-tests-"))
    return _TMP_ROOT


def fixture_dir() -> Path:
    """Build (once per process) and return the synthetic frozen-input directory."""
    global _FIXTURE_DIR
    if _FIXTURE_DIR is not None:
        return _FIXTURE_DIR
    override = os.environ.get("PTLIB_FIXTURE_DIR")
    if override:
        _FIXTURE_DIR = Path(override)
        return _FIXTURE_DIR
    from tools.fixture import build_fixture

    target = tmp_root() / "fixture"
    build_fixture(target)
    _FIXTURE_DIR = target
    return target


def copy_fixture(name: str) -> Path:
    """Copy the fixture into a fresh directory (for mutation/poisoning tests)."""
    source = fixture_dir()
    target = tmp_root() / name
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(source, target)
    return target


def artifact_hashes(data_dir: Path) -> dict[str, str]:
    path = Path(data_dir) / "manifest.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    artifacts = payload.get("artifacts", payload)
    out: dict[str, str] = {}
    for name, record in artifacts.items():
        out[name] = record["sha256"] if isinstance(record, dict) else record
    return out


def rewrite_manifest(data_dir: Path) -> dict[str, str]:
    """Recompute the fixture manifest after a deliberate mutation."""
    data_dir = Path(data_dir)
    manifest_file = data_dir / "manifest.json"
    payload = json.loads(manifest_file.read_text(encoding="utf-8"))
    names = sorted(payload.get("artifacts", {}).keys())
    hashes: dict[str, str] = {}
    for name in names:
        digest = hashlib.sha256()
        with open(data_dir / name, "rb") as handle:
            for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
                digest.update(block)
        hashes[name] = digest.hexdigest()
    payload["artifacts"] = {name: {"sha256": digest} for name, digest in hashes.items()}
    manifest_file.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return hashes


def poison_cm_targets(src: Path, dst: Path, rows: np.ndarray, delta: float = 1.0) -> Path:
    """Copy *src* and add *delta* to the log1p targets (and raw counts) of *rows*."""
    import h5py

    dst = copy_fixture_to(src, dst)
    rows = np.sort(np.asarray(rows, dtype=np.int64))
    with h5py.File(dst / "cm_transfer.h5ad", "r+") as handle:
        for key in ("X", "layers/log1p_cp10k"):
            dataset = handle[key]
            for row in rows:
                dataset[row, :] = dataset[row, :] + delta
    rewrite_manifest(dst)
    return dst


def copy_fixture_to(src: Path, dst: Path) -> Path:
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst)
    return dst


def cm_obs(data_dir: Path) -> pd.DataFrame:
    import anndata as ad

    return ad.read_h5ad(Path(data_dir) / "cm_transfer.h5ad").obs


def tahoe_obs(data_dir: Path) -> pd.DataFrame:
    import anndata as ad

    return ad.read_h5ad(Path(data_dir) / "tahoe_transfer.h5ad").obs


def test_rows_for_fold(data_dir: Path, fold: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(train_rows, valid_rows, test_rows)`` treated positions for a fold."""
    obs = cm_obs(data_dir).reset_index(drop=True)
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
    rows = {"train": [], "valid": [], "test": []}
    for role in rows:
        mask = (~is_control) & qc & (role_of_row == role)
        rows[role] = np.flatnonzero(mask)
    return rows["train"], rows["valid"], rows["test"]


def membership_rows(data_dir: Path, fold: int) -> np.ndarray:
    frame = pd.read_parquet(Path(data_dir) / "consensus_membership.parquet")
    return np.asarray(frame[frame["fold"] == fold]["row"], dtype=np.int64)


def small_train_options(**overrides):
    """CPU training options at fixture scale (fast, no GPU)."""
    from ptlib.train import TrainOptions

    base = dict(
        stage="finetune",
        out_dir=tmp_root() / "unused-run",
        data_dir=fixture_dir(),
        seed=7,
        device="cpu",
        fold=0,
        epochs=1,
        min_epochs=1,
        batch_size=16,
        model_card=False,
        smoke=True,
        command=["test"],
    )
    base.update(overrides)
    return TrainOptions(**base)


_SMALL_FIXTURE_DIR: Path | None = None


def small_fixture_dir() -> Path:
    """A 512-gene fixture for tests that only need the mechanics (queue, cache).

    The full 18,164-gene fixture is used where the real matrix shapes / frozen
    metric matter; the small one keeps subprocess-heavy tests affordable.
    """
    global _SMALL_FIXTURE_DIR
    if _SMALL_FIXTURE_DIR is not None:
        return _SMALL_FIXTURE_DIR
    from tools.fixture import build_fixture

    target = tmp_root() / "fixture-small"
    build_fixture(target, n_genes=512, n_main_drugs=10, n_excluded=2)
    _SMALL_FIXTURE_DIR = target
    return target
