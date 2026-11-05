"""Offline contract tests for `stage2/03_models_cmbench/03.1_cold_drug_18k/output/input/`.

These tests never touch a GPU and never rebuild the data.  They load the
frozen artifacts, re-run the twelve split/integrity gates from
`build_data.run_gates`, and spot-check the Δ/group/true-set tables against
independent recomputations.

Run:
    python -m unittest discover -s stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/tests -v
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
CODE = HERE.parent.parent / "script"  # execution/script
sys.path.insert(0, str(CODE))

from repo_root import BENCH, ROOT  # noqa: E402

DATA = BENCH / "output" / "input"

import build_data  # noqa: E402
import metrics  # noqa: E402


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


class DataContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = json.loads((DATA / "manifest.json").read_text())
        cls.wells = pd.read_parquet(DATA / "wells.parquet")
        cls.groups = pd.read_parquet(DATA / "groups.parquet")
        cls.splits = pd.read_parquet(DATA / "splits.parquet")
        cls.counts = np.load(DATA / "counts.npy", mmap_mode="r")
        cls.log1p = np.load(DATA / "log1p_cp10k.npy", mmap_mode="r")
        cls.delta = np.load(DATA / "delta_log1p_cp10k.npy", mmap_mode="r")
        cls.group_delta = np.load(DATA / "group_delta.npy", mmap_mode="r")
        cls.controls = np.load(DATA / "controls.npz", allow_pickle=False)
        cls.true_sets = {
            k: np.load(DATA / f"true_sets_k{k}.npz", allow_pickle=False)
            for k in (100, 200, 500)
        }
        cls.genes = pd.read_csv(
            ROOT / "stage2/_frozen/cm-data-v1/corpus/v1/genes.csv", usecols=["symbol"]
        )["symbol"].astype(str).to_numpy()

    # -- manifest ---------------------------------------------------------

    def test_artifacts_match_manifest(self) -> None:
        for name, info in self.manifest["artifacts"].items():
            path = DATA / name
            self.assertTrue(path.exists(), name)
            self.assertEqual(sha256_file(path), info["sha256"], name)
            self.assertEqual(path.stat().st_size, info["bytes"], name)

    def test_manifest_source_and_gene_order(self) -> None:
        source = self.manifest["source"]
        self.assertEqual(source["shape"], [4478, 18164])
        self.assertEqual(self.manifest["n_genes"], 18164)
        self.assertEqual(self.manifest["n_groups"], 1380)
        self.assertEqual(self.manifest["n_treated"], 4206)
        order = hashlib.sha256("\n".join(self.genes).encode()).hexdigest()
        self.assertEqual(self.manifest["gene_order_sha256"], order)

    # -- wells / controls / delta ----------------------------------------

    def test_wells_and_delta_row(self) -> None:
        self.assertEqual(len(self.wells), 4478)
        self.assertTrue(self.wells["well"].is_unique)
        treated = self.wells["delta_row"].to_numpy() >= 0
        self.assertEqual(int(treated.sum()), 4206)
        self.assertEqual(int((~treated).sum()), 272)
        self.assertTrue(self.wells.loc[treated, "delta_row"].is_unique)
        self.assertEqual(
            self.wells.loc[treated, "delta_row"].to_numpy().min(), 0
        )
        self.assertEqual(
            self.wells.loc[treated, "delta_row"].to_numpy().max(), 4205
        )
        # Controls are exactly the 272 non-treated rows, and no treated row is
        # a Media/vehicle control.
        self.assertTrue(self.wells.loc[~treated, "is_control"].all())
        self.assertFalse(self.wells.loc[treated, "is_control"].any())

    def test_controls_bundle(self) -> None:
        cg = self.controls["control_groups"].astype(str)
        cells = self.controls["cells"]
        mean = self.controls["mean"]
        wells = self.controls["wells"]
        self.assertEqual(cg.shape, (16,))
        self.assertEqual(cells.shape, (16, 12, 18164))
        self.assertEqual(mean.shape, (16, 18164))
        self.assertEqual(wells.shape, (16, 12))
        self.assertEqual(cells.dtype, np.float32)
        self.assertEqual(mean.dtype, np.float32)
        self.assertTrue(np.isfinite(mean).all())
        # Mean is the float32 mean of the 12 DMSO cells.
        recomputed = cells.mean(axis=1, dtype=np.float64).astype(np.float32)
        self.assertTrue(np.array_equal(recomputed, mean))
        self.assertEqual(set(cg) & set(self.wells["control_group"].astype(str)), set(cg))
        self.assertTrue(set(self.wells["control_group"].astype(str)) <= set(cg))

    def test_delta_is_control_subtracted_log1p(self) -> None:
        cg = list(self.controls["control_groups"].astype(str))
        mean = self.controls["mean"]
        treated_rows = np.flatnonzero(self.wells["delta_row"].to_numpy() >= 0)
        # Sample across plates and both single/multi-plate groups.
        sample = treated_rows[:: max(1, len(treated_rows) // 50)][:50]
        for row in sample:
            group = str(self.wells.iloc[row]["control_group"])
            expected = (
                self.log1p[row].astype(np.float64)
                - mean[cg.index(group)].astype(np.float64)
            ).astype(np.float32)
            got = self.delta[int(self.wells.iloc[row]["delta_row"])]
            self.assertTrue(np.array_equal(expected, got), msg=str(row))

    # -- groups -----------------------------------------------------------

    def test_groups(self) -> None:
        g = self.groups
        self.assertEqual(len(g), 1380)
        self.assertTrue(g["group_id"].is_unique)
        self.assertEqual(g["n_wells"].value_counts().to_dict(), {3: 1326, 4: 48, 6: 6})
        pfoa = g[g["drug_key"] == "SNGREZUHAYWORS-UHFFFAOYSA-N"]
        self.assertEqual(len(pfoa), 3)
        self.assertTrue((pfoa["n_wells"] == 6).all())
        self.assertTrue((pfoa["n_plates"] == 2).all())
        self.assertTrue(pfoa["plates"].str.contains("Plate5").all())
        self.assertTrue(pfoa["plates"].str.contains("Plate9").all())
        tam = g[g["drug_key"] == "NKANXQFJJICGDU-QPLCGJKRSA-N"]
        self.assertEqual(len(tam), 3)
        self.assertTrue(tam["smiles_is_parent"].all())

    def test_group_delta_is_mean_of_well_delta(self) -> None:
        for i in (0, 123, 969, 1379):
            g = self.groups.iloc[i]
            rows = self.wells[
                (self.wells["drug_key"] == g["drug_key"])
                & (self.wells["dose_um"] == g["dose_um"])
            ]["delta_row"].to_numpy()
            expected = self.delta[rows].mean(axis=0, dtype=np.float64).astype(np.float32)
            self.assertTrue(np.array_equal(expected, self.group_delta[i]), msg=str(i))

    # -- true sets --------------------------------------------------------

    def test_true_sets_shapes_and_sizes(self) -> None:
        ids = self.groups["group_id"].astype(str).to_numpy()
        for k, npz in self.true_sets.items():
            self.assertEqual(npz["labels"].shape, (1380, 18164))
            self.assertEqual(npz["labels"].dtype, np.int8)
            self.assertTrue(np.array_equal(npz["group_ids"].astype(str), ids))
            self.assertTrue(np.array_equal(npz["n_up"], np.full(1380, k, np.int32)))
            self.assertTrue(np.array_equal(npz["n_down"], np.full(1380, k, np.int32)))
            self.assertEqual(int(np.count_nonzero(npz["labels"] == 1)), 1380 * k)
            self.assertEqual(int(np.count_nonzero(npz["labels"] == -1)), 1380 * k)
            # Up and down are disjoint in every row.
            self.assertEqual(int(np.count_nonzero(npz["labels"] == 1, axis=1).min()), k)
            self.assertEqual(int(np.count_nonzero(npz["labels"] == -1, axis=1).min()), k)

    def test_true_sets_match_metric_code(self) -> None:
        ids = self.groups["group_id"].astype(str).to_numpy()
        for k, npz in self.true_sets.items():
            for i in (0, 123, 500, 969, 1379):
                labels = npz["labels"][i]
                up = np.flatnonzero(labels == 1)
                down = np.flatnonzero(labels == -1)
                ref_up, ref_down = metrics.true_set_indices(
                    self.group_delta[i].astype(np.float64), k
                )
                self.assertTrue(np.array_equal(np.sort(up), np.sort(ref_up)))
                self.assertTrue(np.array_equal(np.sort(down), np.sort(ref_down)))
                self.assertEqual(ids[i], str(npz["group_ids"][i]))

    # -- splits -----------------------------------------------------------

    def test_splits_rerun_all_gates(self) -> None:
        gates = build_data.run_gates(
            obs=self.wells,
            X=self.counts,
            log1p=self.log1p,
            var_names=self.genes,
            groups=self.groups,
            splits=self.splits,
            control_groups=self.controls["control_groups"],
            cells=self.controls["cells"],
            control_mean=self.controls["mean"],
            delta=self.delta,
            delta_row=self.wells["delta_row"].to_numpy(),
            group_delta=self.group_delta,
            true_sets={k: {name: npz[name] for name in npz.files} for k, npz in self.true_sets.items()},
        )
        self.assertTrue(all(gates.values()), gates)

    def test_split_shape_and_roles(self) -> None:
        s = self.splits
        self.assertEqual(len(s), 21030)
        self.assertFalse(s.duplicated(["well", "fold"]).any())
        self.assertEqual(set(s["role"]), {"train", "valid", "test"})
        for k in range(5):
            f = s[s["fold"] == k]
            self.assertEqual(f.loc[f["role"] == "test", "drug_key"].nunique(), 92)
            self.assertEqual(f.loc[f["role"] == "valid", "drug_key"].nunique(), 92)
            self.assertEqual(f.loc[f["role"] == "train", "drug_key"].nunique(), 276)
            self.assertEqual(f["well"].nunique(), 4206)

    def test_no_quarantined_paths(self) -> None:
        self.assertEqual(self.counts.shape[1], 18164)
        self.assertEqual(self.log1p.shape[1], 18164)
        self.assertNotIn("978", " ".join(self.manifest["artifacts"]))
        self.assertNotIn("978", self.manifest["source"]["path"])


if __name__ == "__main__":
    unittest.main()
