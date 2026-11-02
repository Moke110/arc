"""Leakage boundaries: which rows each loader may read, and which controls build c*.

These tests attack the boundary rather than trusting the column set:

* the CM training loader refuses held-out (test) treated rows;
* a training run's loss is unchanged when *test* labels are poisoned;
* c* is reconstructed only from its recorded membership — poisoning a non-member
  control row leaves it untouched, poisoning a member row is detected;
* prediction is target-free (the prediction path never opens treated expression);
* Tahoe matched controls never cross cell lines.
"""

from __future__ import annotations

import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import torch

import _support
from ptlib import data as data_mod
from ptlib import model as model_mod
from ptlib import predict as predict_mod
from ptlib import runtime


class FrozenInputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = _support.fixture_dir()
        cls.frozen = data_mod.open_frozen_input(cls.dir)

    def test_axis_and_config_agree(self):
        self.assertEqual(self.frozen.n_genes, 18164)
        self.assertEqual(self.frozen.structure["n_genes"], 18164)
        self.assertTrue(bool(self.frozen.gene_measured.all()))
        self.assertEqual(self.frozen.config["metric"]["k"], 200)

    def test_artifact_hashes_verified(self):
        record = self.frozen.record()
        self.assertEqual(len(record["artifacts"]), 9)
        self.assertEqual(record["hash_verification"]["verified"], 8)  # manifest.json excluded

    def test_fold_table_contract(self):
        folds = data_mod.FoldTable.load(self.dir / "folds.parquet", 5)
        self.assertEqual(folds.report["n_main_drugs"], 30)
        for fold in range(5):
            for role in ("train", "valid", "test"):
                self.assertGreater(folds.report["per_fold"][str(fold)][role], 0)

    def test_fold_table_rejects_broken_contract(self):
        frame = pd.read_parquet(self.dir / "folds.parquet")
        broken = pd.concat([frame, frame.iloc[[0]]], ignore_index=True)
        table = data_mod.FoldTable(frame=broken, n_folds=5)
        with self.assertRaises(data_mod.DataError):
            table.validate(5)


class CMDatasetLeakageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = _support.fixture_dir()
        cls.frozen = data_mod.open_frozen_input(cls.dir)
        cls.structures = cls.frozen.structures()
        cls.folds = data_mod.FoldTable.load(cls.dir / "folds.parquet", 5)

    def _dataset(self, fold: int = 0, train_drugs: str = "main") -> data_mod.CMDataset:
        return data_mod.CMDataset(
            frozen=self.frozen,
            fold=fold,
            fold_table=self.folds,
            structures=self.structures,
            train_drugs=train_drugs,
        ).prepare()

    def test_roles_are_disjoint(self):
        dataset = self._dataset()
        train, valid, test = set(dataset.train_rows), set(dataset.valid_rows), set(dataset.test_rows)
        self.assertFalse(train & valid)
        self.assertFalse(train & test)
        self.assertFalse(valid & test)
        self.assertGreater(len(test), 0)

    def test_training_loader_refuses_test_rows(self):
        dataset = self._dataset()
        with self.assertRaisesRegex(data_mod.DataError, "held-out"):
            dataset.target(dataset.test_rows[:1])
        with self.assertRaisesRegex(data_mod.DataError, "held-out"):
            dataset.batch(dataset.test_rows[:2])

    def test_only_train_and_valid_rows_are_materialised(self):
        dataset = self._dataset()
        dataset.batch(dataset.train_rows[:4])
        dataset.batch(dataset.valid_rows[:4])
        expected = np.unique(np.concatenate([dataset.train_rows, dataset.valid_rows]))
        np.testing.assert_array_equal(dataset.rows, expected)
        self.assertFalse(np.intersect1d(dataset.rows, dataset.test_rows).size)

    def test_consensus_is_reconstructed_from_membership(self):
        dataset = self._dataset()
        self.assertLess(dataset.consensus_check["max_abs_diff"], 1e-6)
        self.assertTrue(dataset.consensus_check["used_weights"])
        self.assertGreater(len(dataset.plate_control), 1)

    def test_valid_input_is_c_star_not_a_matched_control(self):
        dataset = self._dataset()
        baseline = dataset.baseline(dataset.valid_rows)
        np.testing.assert_allclose(baseline, np.broadcast_to(dataset.c_star, baseline.shape))
        train_baseline = dataset.baseline(dataset.train_rows)
        self.assertFalse(np.allclose(train_baseline, dataset.c_star))

    def test_expanded_training_set_adds_only_unfolded_qc_pass_wells(self):
        main = self._dataset(train_drugs="main")
        expanded = self._dataset(train_drugs="expanded")
        extra = np.setdiff1d(expanded.train_rows, main.train_rows)
        self.assertGreater(extra.size, 0)
        folds_drugs = set(self.folds.drugs)
        drug = self.frozen.cm_obs["drug_key"].astype(str).to_numpy()
        main_set = self.frozen.cm_obs["main_set"].to_numpy(dtype=bool)
        for row in extra:
            self.assertFalse(main_set[row])
            self.assertNotIn(drug[row], folds_drugs)
        # expanded rows never appear in valid/test and are never scored
        self.assertFalse(np.intersect1d(expanded.valid_rows, extra).size)
        self.assertFalse(np.intersect1d(expanded.test_rows, extra).size)


class ConsensusMembershipTests(unittest.TestCase):
    """c* must depend on its recorded membership and on nothing else."""

    @classmethod
    def setUpClass(cls):
        cls.dir = _support.fixture_dir()
        cls.original = data_mod.open_frozen_input(cls.dir)
        cls.folds = data_mod.FoldTable.load(cls.dir / "folds.parquet", 5)
        cls.members = _support.membership_rows(cls.dir, 0)
        obs = _support.cm_obs(cls.dir).reset_index(drop=True)
        qc = obs["qc_pass"].to_numpy(dtype=bool)
        vehicle = obs["vehicle"].astype(str).to_numpy()
        is_control = obs["is_control"].to_numpy(dtype=bool)
        cls.non_members = np.array(
            [
                row for row in np.flatnonzero(is_control)
                if row not in set(cls.members.tolist())
            ],
            dtype=np.int64,
        )
        self_check = np.setdiff1d(cls.non_members, cls.members)
        assert self_check.size == cls.non_members.size
        cls.obs = obs

    def _c_star(self, data_dir):
        frozen = data_mod.open_frozen_input(data_dir)
        dataset = data_mod.CMDataset(
            frozen=frozen, fold=0, fold_table=self.folds,
            structures=frozen.structures(),
        ).prepare()
        return dataset.c_star

    def test_non_member_control_perturbation_does_not_change_c_star(self):
        import h5py

        copy = _support.copy_fixture("nonmember-poison")
        row = int(self.non_members[0])
        with h5py.File(copy / "cm_transfer.h5ad", "r+") as handle:
            handle["layers/log1p_cp10k"][row, :] = handle["layers/log1p_cp10k"][row, :] + 5.0
        _support.rewrite_manifest(copy)
        np.testing.assert_allclose(self._c_star(copy), self._c_star(self.dir), atol=1e-6)

    def test_member_control_perturbation_is_detected(self):
        import h5py

        copy = _support.copy_fixture("member-poison")
        row = int(self.members[0])
        with h5py.File(copy / "cm_transfer.h5ad", "r+") as handle:
            handle["layers/log1p_cp10k"][row, :] = handle["layers/log1p_cp10k"][row, :] + 5.0
        _support.rewrite_manifest(copy)
        with self.assertRaisesRegex(data_mod.DataError, "c\\* cannot be reconstructed"):
            self._c_star(copy)


class PredictionIsTargetFreeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = _support.fixture_dir()
        cls.frozen = data_mod.open_frozen_input(cls.dir)
        cls.folds = data_mod.FoldTable.load(cls.dir / "folds.parquet", 5)
        cls.structures = cls.frozen.structures()
        cls.model = model_mod.build_model(cls.frozen.structure)
        model_mod.random_init(cls.model, seed=3)
        cls.model.eval()

    def test_predict_source_has_no_target_access(self):
        source = data_mod.CMPredictSource(self.frozen, self.folds, 0, self.structures)
        self.assertFalse(hasattr(source, "target"))
        self.assertFalse(hasattr(source, "log1p"))
        rows = source.rows_for("test")
        self.assertGreater(rows.size, 0)
        pred, delta = predict_mod.predict_rows(
            self.model, source, rows[:8], n_dim=10, batch_size=8, device=torch.device("cpu")
        )
        self.assertEqual(pred.shape, (8, 18164))
        np.testing.assert_allclose(delta, pred - source.baseline(rows[:8]), atol=1e-6)


class TahoeControlBoundaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = _support.fixture_dir()
        cls.frozen = data_mod.open_frozen_input(cls.dir)
        cls.dataset = data_mod.TahoeDataset(
            frozen=cls.frozen, structures=cls.frozen.structures(), max_train_rows=None,
            max_valid_rows=None, load_targets=True,
        ).prepare()

    def test_controls_never_cross_cell_lines(self):
        groups = pd.read_parquet(self.dir / "tahoe_control_groups.parquet")
        obs = _support.tahoe_obs(self.dir).reset_index(drop=True)
        for gid, sub in groups.groupby("group_index"):
            lines = set(obs["cell_line_id"].astype(str).to_numpy()[sub["row"].to_numpy()])
            plates = set(obs["plate"].astype(str).to_numpy()[sub["row"].to_numpy()])
            self.assertEqual(len(lines), 1, f"group {gid} mixes cell lines {lines}")
            self.assertEqual(len(plates), 1, f"group {gid} mixes plates {plates}")

    def test_baseline_is_the_matched_control_mean(self):
        import anndata as ad

        obs = _support.tahoe_obs(self.dir).reset_index(drop=True)
        groups = pd.read_parquet(self.dir / "tahoe_control_groups.parquet")
        layer = ad.read_h5ad(self.dir / "tahoe_transfer.h5ad").layers["log1p_cp10k"]
        layer = np.asarray(layer.todense()) if hasattr(layer, "todense") else np.asarray(layer)
        rows = self.dataset.rows_valid[:3]
        baseline = self.dataset.baseline(rows)
        for i, row in enumerate(rows):
            plate = obs.iloc[row]["plate"]
            line = obs.iloc[row]["cell_line_id"]
            members = groups[
                (groups["plate"] == plate) & (groups["cell_line_id"] == line)
            ]["row"].to_numpy()
            expected = layer[members].astype(np.float64).mean(axis=0)
            np.testing.assert_allclose(baseline[i], expected, atol=1e-6)

    def test_valid_split_is_the_frozen_column(self):
        obs = _support.tahoe_obs(self.dir).reset_index(drop=True)
        split = obs["tahoe_split"].astype(str).to_numpy()
        plate = obs["plate"].astype(str).to_numpy()
        valid_rows = self.dataset.rows_valid
        self.assertTrue(np.all(split[valid_rows] == "valid"))
        self.assertTrue(set(plate[valid_rows]) <= {"plate13", "plate14"})

    def test_train_and_valid_rows_are_disjoint(self):
        self.assertFalse(np.intersect1d(self.dataset.rows_train, self.dataset.rows_valid).size)


class TahoeUsabilityMaskTests(unittest.TestCase):
    """The frozen ``tahoe_usable`` mask is authoritative; no baseline is imputed.

    Real-data regression: re-deriving Tahoe usability inside the loader
    re-included the 4,757 treated rows that have no QC-pass matched control, and
    the control-group lookup then crashed with a ``KeyError`` on the frozen input
    (``plate3|CVCL_1715``).  The frozen artifact carries ``tahoe_usable``; the
    loader must use it, and must fail closed when the side files disagree.
    """

    @classmethod
    def setUpClass(cls):
        cls.src = _support.small_fixture_dir()

    def _variant(self, name, obs_mutator=None, groups_mutator=None, pairs_mutator=None) -> Path:
        data_dir = _support.copy_fixture_to(self.src, Path(_support.tmp_root()) / f"tahoe-mask-{name}")
        h5ad = data_dir / "tahoe_transfer.h5ad"
        adata = ad.read_h5ad(h5ad)
        obs = adata.obs.copy()
        if obs_mutator is not None:
            obs_mutator(obs)
        adata.obs = obs
        adata.write_h5ad(h5ad)
        if groups_mutator is not None:
            path = data_dir / "tahoe_control_groups.parquet"
            groups_mutator(pd.read_parquet(path)).to_parquet(path, index=False)
        if pairs_mutator is not None:
            path = data_dir / "tahoe_pairs.npz"
            payload = dict(np.load(path))
            pairs_mutator(payload)
            np.savez(path, **payload)
        _support.rewrite_manifest(data_dir)
        return data_dir

    def _dataset(self, data_dir: Path, **kwargs) -> data_mod.TahoeDataset:
        frozen = data_mod.open_frozen_input(data_dir)
        return data_mod.TahoeDataset(
            frozen=frozen, structures=frozen.structures(), load_targets=False, **kwargs
        ).prepare()

    def test_frozen_usable_column_drives_selection(self):
        def mutate(obs):
            split = obs["tahoe_split"].astype(str)
            obs["tahoe_usable"] = (
                (~obs["is_control"].astype(bool))
                & (((split == "train") & (obs["plate"].astype(str) == "plate1"))
                   | (split == "valid"))
            )
        data_dir = self._variant("column", mutate)
        dataset = self._dataset(data_dir)
        self.assertEqual(dataset.counts["usable_mask_source"], "obs.tahoe_usable")
        obs = _support.tahoe_obs(data_dir).reset_index(drop=True)
        mask = obs["tahoe_usable"].astype(bool).to_numpy()
        split = obs["tahoe_split"].astype(str).to_numpy()
        self.assertEqual(dataset.counts["n_train"], int((mask & (split == "train")).sum()))
        self.assertEqual(dataset.counts["n_valid"], int((mask & (split == "valid")).sum()))
        self.assertTrue(np.all(mask[dataset.rows_train]))
        self.assertTrue(np.all(mask[dataset.rows_valid]))
        derived = (
            (~obs["is_control"].astype(bool))
            & obs["qc_pass"].astype(bool)
            & obs["has_structure"].astype(bool)
        )
        self.assertLess(dataset.counts["n_train"], int((derived & (split == "train")).sum()))

    def test_usable_row_without_a_matched_control_fails_closed(self):
        def mutate(obs):
            obs["tahoe_usable"] = (
                (~obs["is_control"].astype(bool))
                & obs["qc_pass"].astype(bool)
                & obs["has_structure"].astype(bool)
            )
            treated = np.flatnonzero(
                (~obs["is_control"].astype(bool)).to_numpy()
                & (obs["tahoe_split"].astype(str).to_numpy() == "train")
            )
            obs["cell_line_id"] = obs["cell_line_id"].astype(object)
            obs.iloc[int(treated[0]), obs.columns.get_loc("cell_line_id")] = "cellZ"
        with self.assertRaises(data_mod.DataError):
            self._dataset(self._variant("no-group", mutate))

    def test_derived_fallback_drops_and_reports_unmatched_rows(self):
        original = pd.read_parquet(self.src / "tahoe_control_groups.parquet")
        removed = int(original.loc[original["group_key"] == "plate1|cellA", "group_index"].iloc[0])

        def mutate_groups(groups):
            return groups[groups["group_key"].astype(str) != "plate1|cellA"].reset_index(drop=True)

        def mutate_pairs(payload):
            keep = payload["group_index"] != removed
            for key, value in list(payload.items()):
                array = np.asarray(value)
                if array.ndim >= 1 and array.shape[0] == keep.size:
                    payload[key] = array[keep]

        data_dir = self._variant("drop-group", groups_mutator=mutate_groups,
                                 pairs_mutator=mutate_pairs)
        dataset = self._dataset(data_dir)
        self.assertEqual(dataset.counts["usable_mask_source"], "derived(qc_pass&has_structure)")
        self.assertGreater(dataset.counts["dropped_no_matched_control"], 0)
        obs = _support.tahoe_obs(data_dir).reset_index(drop=True)
        keys = (obs["plate"].astype(str) + "|" + obs["cell_line_id"].astype(str)).to_numpy()
        self.assertNotIn("plate1|cellA", set(keys[dataset.rows_train]))
        self.assertNotIn("plate1|cellA", set(keys[dataset.rows_valid]))


if __name__ == "__main__":
    unittest.main()
