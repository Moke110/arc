"""Per-well case bundles for the Stage 2 benchmark (ticket 05/12).

A case bundle is one `.npz` per phase (train/valid/test).  Every row is one
**well**.  Train and valid carry expression and observed Δ; test carries
metadata only and is checked for forbidden expression keys at serialization
time.

The bundles are model-independent, so the runner builds them once per fold and
all adapters read the same files.  They live under
`output/_cases/fold_<k>/`.
"""

from __future__ import annotations

import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

# Keys that must never appear in a test bundle.
FORBIDDEN_TEST_KEYS = {
    "X",
    "log1p_cp10k",
    "delta",
    "well_row",
    "delta_row",
    "counts_row",
    "lib_size",
    "n_genes_detected",
}

META_KEYS = (
    "well_id",
    "group_id",
    "drug_key",
    "compound_name",
    "smiles",
    "dose_um",
    "dose_log",
    "time_h",
    "cell_type",
    "plate",
    "control_group",
    "replicate",
    "pharm_atc_moa",
)

EXPRESSION_KEYS = ("X", "log1p_cp10k", "delta", "well_row", "delta_row")


def save_npz(path: Path, arrays: dict[str, np.ndarray]) -> None:
    """Deterministic `.npz` writer (fixed 1980 timestamp, sorted keys)."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_STORED) as archive:
        for name in sorted(arrays):
            buffer = io.BytesIO()
            array = np.asarray(arrays[name])
            if array.ndim > 0:
                array = np.ascontiguousarray(array)
            np.lib.format.write_array(buffer, array, allow_pickle=False)
            info = zipfile.ZipInfo(name + ".npy", date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_STORED
            info.external_attr = 0o600 << 16
            archive.writestr(info, buffer.getvalue())


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


class BenchData:
    """Lazy access to the frozen `data/` artifacts."""

    def __init__(self, data_dir: Path):
        self.data_dir = Path(data_dir)
        self.wells = pd.read_parquet(self.data_dir / "wells.parquet")
        self.groups = pd.read_parquet(self.data_dir / "groups.parquet")
        self.splits = pd.read_parquet(self.data_dir / "splits.parquet")
        self.counts = np.load(self.data_dir / "counts.npy", mmap_mode="r")
        self.log1p = np.load(self.data_dir / "log1p_cp10k.npy", mmap_mode="r")
        self.delta = np.load(self.data_dir / "delta_log1p_cp10k.npy", mmap_mode="r")
        self.group_delta = np.load(self.data_dir / "group_delta.npy", mmap_mode="r")
        self.controls = np.load(self.data_dir / "controls.npz", allow_pickle=False)

        self.well_row = {str(w): i for i, w in enumerate(self.wells["well"])}
        self.group_row = {str(g): i for i, g in enumerate(self.groups["group_id"])}
        self.group_meta = self.groups.set_index("group_id", drop=False)
        self.manifest = json.loads((self.data_dir / "manifest.json").read_text())

        if len(self.wells) != 4478 or len(self.groups) != 1380:
            raise RuntimeError("data artifacts are not the ticket-11 frozen build")
        if len(self.splits) != 21030:
            raise RuntimeError("splits.parquet is not the 5-fold long table")

    # -- case construction -------------------------------------------------

    def case_frame(self, fold: int, role: str) -> pd.DataFrame:
        if role not in {"train", "valid", "test"}:
            raise ValueError(role)
        sub = self.splits[(self.splits["fold"] == fold) & (self.splits["role"] == role)]
        if not len(sub):
            raise RuntimeError(f"no {role} rows for fold {fold}")
        sub = sub.sort_values("well", kind="mergesort").reset_index(drop=True)

        meta = self.groups[
            ["group_id", "drug_key", "dose_um", "compound_name"]
        ].rename(columns={"compound_name": "group_compound_name"})
        sub = sub.merge(
            meta,
            on=["drug_key", "dose_um"],
            how="left",
            validate="many_to_one",
        )
        if sub["group_id"].isna().any():
            raise RuntimeError(f"fold {fold} {role}: missing group metadata")

        well_meta = self.wells[["well", "time_h", "cell_type", "delta_row"]]
        sub = sub.merge(well_meta, on="well", how="left", validate="one_to_one")
        if sub["delta_row"].isna().any():
            raise RuntimeError(f"fold {fold} {role}: missing well metadata")

        sub["compound_name"] = sub["group_compound_name"].fillna("").astype(str)
        sub["smiles"] = sub["smiles"].fillna("").astype(str)
        sub["well_row"] = sub["well"].map(self.well_row).astype(np.int32)
        sub["dose_log"] = np.log10(sub["dose_um"].astype(np.float64) + 1.0)
        return sub

    def build_arrays(self, fold: int, role: str) -> dict[str, np.ndarray]:
        frame = self.case_frame(fold, role)
        arrays: dict[str, np.ndarray] = {
            "well_id": frame["well"].astype(str).to_numpy(dtype=str),
            "group_id": frame["group_id"].astype(str).to_numpy(dtype=str),
            "drug_key": frame["drug_key"].astype(str).to_numpy(dtype=str),
            "compound_name": frame["compound_name"].astype(str).to_numpy(dtype=str),
            "smiles": frame["smiles"].astype(str).to_numpy(dtype=str),
            "dose_um": frame["dose_um"].astype(np.float64).to_numpy(),
            "dose_log": frame["dose_log"].astype(np.float64).to_numpy(),
            "time_h": frame["time_h"].astype(np.float64).to_numpy(),
            "cell_type": frame["cell_type"].astype(str).to_numpy(dtype=str),
            "plate": frame["plate"].astype(str).to_numpy(dtype=str),
            "control_group": frame["control_group"].astype(str).to_numpy(dtype=str),
            "replicate": frame["replicate"].astype(np.int32).to_numpy(),
            "pharm_atc_moa": frame["pharm_atc_moa"].astype(bool).to_numpy(),
        }
        if role != "test":
            well_row = frame["well_row"].to_numpy(dtype=np.int32)
            delta_row = frame["delta_row"].to_numpy(dtype=np.int32)
            arrays["well_row"] = well_row
            arrays["delta_row"] = delta_row
            arrays["X"] = np.asarray(self.counts[well_row], dtype=np.int32)
            arrays["log1p_cp10k"] = np.asarray(self.log1p[well_row], dtype=np.float32)
            arrays["delta"] = np.asarray(self.delta[delta_row], dtype=np.float32)
        return arrays

    # -- serialization -----------------------------------------------------

    def write_cases(self, fold: int, out_dir: Path, force: bool = False) -> dict[str, Path]:
        out_dir = Path(out_dir)
        paths = {
            role: out_dir / f"{role}.npz" for role in ("train", "valid", "test")
        }
        if not force and all(p.exists() for p in paths.values()):
            return paths
        for role, path in paths.items():
            save_npz(path, self.build_arrays(fold, role))
        return paths

    def verify_cases(self, fold: int, case_paths: dict[str, Path]) -> dict:
        """Fail-closed leakage and integrity checks for one fold's bundles."""

        loaded = {
            role: np.load(path, allow_pickle=False)
            for role, path in case_paths.items()
        }
        try:
            test_keys = set(loaded["test"].files)
            bad = test_keys & FORBIDDEN_TEST_KEYS
            if bad:
                raise RuntimeError(f"test bundle contains forbidden keys: {sorted(bad)}")
            missing = set(META_KEYS) - test_keys
            if missing:
                raise RuntimeError(f"test bundle missing metadata: {sorted(missing)}")

            for role in ("train", "valid"):
                missing = (set(META_KEYS) | set(EXPRESSION_KEYS)) - set(loaded[role].files)
                if missing:
                    raise RuntimeError(f"{role} bundle missing keys: {sorted(missing)}")
                x = loaded[role]["X"]
                log1p = loaded[role]["log1p_cp10k"]
                delta = loaded[role]["delta"]
                if x.dtype != np.int32 or log1p.dtype != np.float32 or delta.dtype != np.float32:
                    raise RuntimeError(f"{role} bundle dtype mismatch")
                if x.shape != log1p.shape or x.shape != delta.shape:
                    raise RuntimeError(f"{role} bundle shape mismatch")
                if int(x.min()) < 0 or not np.isfinite(log1p).all() or not np.isfinite(delta).all():
                    raise RuntimeError(f"{role} bundle invalid values")
                if not np.array_equal(
                    delta, np.asarray(self.delta[loaded[role]["delta_row"]], dtype=np.float32)
                ):
                    raise RuntimeError(f"{role} delta does not match the frozen Δ")
                if not np.array_equal(
                    x, np.asarray(self.counts[loaded[role]["well_row"]], dtype=np.int32)
                ):
                    raise RuntimeError(f"{role} X does not match counts.npy")

            train_ids = set(map(str, loaded["train"]["well_id"]))
            valid_ids = set(map(str, loaded["valid"]["well_id"]))
            test_ids = set(map(str, loaded["test"]["well_id"]))
            if not (train_ids | valid_ids | test_ids) == set(
                self.splits[self.splits["fold"] == fold]["well"].astype(str)
            ):
                raise RuntimeError("case wells do not partition the fold")
            if (train_ids & valid_ids) or (train_ids & test_ids) or (valid_ids & test_ids):
                raise RuntimeError("case well ids overlap")

            train_groups = set(map(str, loaded["train"]["group_id"]))
            valid_groups = set(map(str, loaded["valid"]["group_id"]))
            test_groups = set(map(str, loaded["test"]["group_id"]))
            if train_groups & test_groups or valid_groups & test_groups:
                raise RuntimeError("test groups leak into train/valid")
            if set(map(str, loaded["train"]["drug_key"])) & set(map(str, loaded["test"]["drug_key"])):
                raise RuntimeError("test drug keys leak into train")
            if set(map(str, loaded["valid"]["drug_key"])) & set(map(str, loaded["test"]["drug_key"])):
                raise RuntimeError("test drug keys leak into valid")

            return {
                "fold": fold,
                "n_train": int(len(loaded["train"]["well_id"])),
                "n_valid": int(len(loaded["valid"]["well_id"])),
                "n_test": int(len(loaded["test"]["well_id"])),
                "n_train_groups": len(train_groups),
                "n_valid_groups": len(valid_groups),
                "n_test_groups": len(test_groups),
                "test_keys": sorted(test_keys),
            }
        finally:
            for npz in loaded.values():
                npz.close()

    def cases_manifest(self, fold: int, case_paths: dict[str, Path]) -> dict:
        manifest = {"fold": fold, "files": {}}
        for role, path in case_paths.items():
            with np.load(path, allow_pickle=False) as npz:
                manifest["files"][role] = {
                    "path": str(path),
                    "sha256": sha256_file(path),
                    "bytes": path.stat().st_size,
                    "keys": sorted(npz.files),
                    "n_cases": int(len(npz["well_id"])),
                }
        return manifest
