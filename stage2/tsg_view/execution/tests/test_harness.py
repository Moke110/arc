"""Public-boundary tests for the shared view harness (product writer + primitives)."""
import sys
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
EXECUTION = HERE.parent
for path in (str(EXECUTION), str(HERE)):
    if path not in sys.path:
        sys.path.insert(0, path)

import tsg_view_fixtures as fixtures  # noqa: E402
import views  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
GENES = np.asarray([0, 3], dtype=np.int64)


def normalized(counts):
    counts = np.asarray(counts, dtype=np.int64)
    return np.log1p(counts * 10000.0 / counts.sum())


class HarnessTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        views.bind(REPO)

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.fixture = fixtures.grid_fixture(self.root)
        self.accepted, self.accepted_manifest, _ = fixtures.build_view(self.fixture)

    # -- helpers ---------------------------------------------------------- #
    @contextmanager
    def corpus(self, source="tenv"):
        corpus = views.modules()["corpus"].open_corpus(self.fixture["corpus"])
        reader = corpus.source(source)
        try:
            yield reader
        finally:
            reader.close()

    def accepted_obs(self, source="tenv"):
        with views.modules()["view"].open_view(self.accepted) as release:
            with release.source(source) as reader:
                return reader.obs.copy()

    def context(self, view, source="tenv", size=2):
        return views.Context(view=view, source=source,
                             output_dir=self.root / "output",
                             staging_dir=self.root / "staging" / f"{view}-{source}",
                             repo_root_path=REPO, size=size, fixture=True,
                             corpus_dir=self.fixture["corpus"],
                             axis_dir=self.fixture["axis"],
                             accepted_dir=self.accepted, provenance_root=None)

    def counts(self, source="tenv"):
        return self.fixture["sources"][source]["counts"]

    # -- condition model -------------------------------------------------- #
    def test_condition_table_has_one_row_per_condition(self):
        conditions, row_to_condition = views.condition_table(self.accepted_obs(), "tenv")
        self.assertEqual(len(conditions), 9)
        self.assertEqual(int((row_to_condition >= 0).sum()), 16)
        self.assertEqual(sorted(conditions.member_count.tolist()),
                         [1, 1, 2, 2, 2, 2, 2, 2, 2])
        first = conditions.iloc[0]
        self.assertEqual(first["cell_line"], "hiPSC-CM")
        self.assertEqual(views.merge_profile_id("tenv", first["cell_line"],
                                                first["drug_parent_inchikey"],
                                                first["dose_um"]),
                         "merge::tenv::hiPSC-CM::P1::1")
        # A cross-plate condition records both matched control groups.
        cross = conditions[conditions.plates.apply(len) == 2]
        self.assertEqual(len(cross), 1)
        self.assertEqual(sorted(cross.iloc[0]["control_groups"]),
                         ["tenv::Plate1::DMSO", "tenv::Plate2::DMSO"])

    # -- pooling ---------------------------------------------------------- #
    def test_pooled_profile_is_the_normalized_sum_of_member_counts(self):
        conditions, row_to_condition = views.condition_table(self.accepted_obs(), "tenv")
        with self.corpus() as reader:
            pooled = views.pooled_condition_counts(reader, GENES, row_to_condition, conditions)
            group_sums = views.vehicle_group_sums(reader, GENES)
            control, retained = views.pooled_control(conditions, group_sums)
            blocks = list(views.pooled_condition_blocks(pooled, control))
        x = np.concatenate([block[1]["X"] for block in blocks])
        retained = np.asarray(retained)
        counts = self.counts()
        # Condition 0 is P1 at 1 uM with members rows 0 and 1 on Plate1.
        expected_counts = np.asarray(counts[0])[GENES] + np.asarray(counts[1])[GENES]
        expected = normalized(expected_counts).astype(np.float32)
        np.testing.assert_allclose(x[0], expected, rtol=0, atol=1e-6)
        # The pooled profile is not the mean of the members' normalized values.
        member_mean = (normalized(np.asarray(counts[0])[GENES])
                       + normalized(np.asarray(counts[1])[GENES])) / 2
        self.assertFalse(np.allclose(x[0], member_mean, rtol=0, atol=1e-3))
        self.assertEqual(int(pooled[0].sum()), int(expected_counts.sum()))
        self.assertEqual(int(retained[0]), int(np.asarray(counts[16])[GENES].sum()
                                              + np.asarray(counts[17])[GENES].sum()))
        np.testing.assert_array_equal(pooled[0], expected_counts)

    def test_pooled_control_sums_the_vehicle_union(self):
        conditions, _ = views.condition_table(self.accepted_obs(), "tenv")
        with self.corpus() as reader:
            group_sums = views.vehicle_group_sums(reader, GENES)
            control, retained = views.pooled_control(conditions, group_sums)
        counts = self.counts()
        # The cross-plate condition (P1 at 100 uM) unions Plate1 and Plate2 vehicles.
        cross = int(np.flatnonzero(conditions.plates.apply(len).to_numpy() == 2)[0])
        self.assertEqual(conditions.iloc[cross]["drug_parent_inchikey"], "P1")
        self.assertEqual(float(conditions.iloc[cross]["dose_um"]), 100.0)
        plate1 = np.asarray(counts[16])[GENES] + np.asarray(counts[17])[GENES]
        plate2 = np.asarray(counts[18])[GENES] + np.asarray(counts[19])[GENES]
        expected = normalized(plate1 + plate2).astype(np.float32)
        np.testing.assert_allclose(control[cross], expected, rtol=0, atol=1e-6)
        self.assertEqual(int(retained[cross]), int((plate1 + plate2).sum()))
        # It is not the mean of the two plates' independently normalized controls.
        plates_mean = (normalized(plate1) + normalized(plate2)) / 2
        self.assertFalse(np.allclose(control[cross], plates_mean, rtol=0, atol=1e-3))
        # A single-plate condition still differs from any single well's normalization.
        plate1_only = normalized(np.asarray(counts[16])[GENES]).astype(np.float32)
        self.assertNotEqual(control[0][0], plate1_only[0])

    # -- solo draw -------------------------------------------------------- #
    def test_solo_draw_is_stable_and_independent_of_row_order(self):
        obs = self.accepted_obs()
        conditions, _ = views.condition_table(obs, "tenv")
        first = views.solo_choice(conditions, rule="draw", source="tenv")
        second = views.solo_choice(conditions, rule="draw", source="tenv")
        np.testing.assert_array_equal(first, second)
        chosen_ids = [conditions.iloc[i]["profile_ids"][list(
            conditions.iloc[i]["members"]).index(int(row))] for i, row in enumerate(first)]
        # Reordering the view rows (as a different fold layout would) cannot change the draw.
        shuffled = obs.sample(frac=1.0, random_state=19)
        shuffled_conditions, _ = views.condition_table(shuffled, "tenv")
        shuffled_choice = views.solo_choice(shuffled_conditions, rule="draw", source="tenv")
        shuffled_ids = [shuffled_conditions.iloc[i]["profile_ids"][list(
            shuffled_conditions.iloc[i]["members"]).index(int(row))]
            for i, row in enumerate(shuffled_choice)]
        self.assertEqual(chosen_ids, shuffled_ids)
        # Singletons always keep their only member.
        singletons = conditions.member_count == 1
        for position in np.flatnonzero(singletons.to_numpy()):
            self.assertEqual(int(first[position]), int(conditions.iloc[position]["members"][0]))
        # A different source or seed changes the draw for multi-member conditions.
        other_seed = views.solo_choice(conditions, rule="draw", seed=1, source="tenv")
        self.assertFalse(np.array_equal(first, other_seed))
        # The member table can carry a per-condition identity (as solo/pooled-label need).
        table = views.member_table(conditions, chosen_ids)
        self.assertEqual(sorted(table.profile_id.unique()), sorted(chosen_ids))
        self.assertEqual(len(table), int(conditions.member_count.sum()))

    def test_solo_deepest_picks_the_largest_retained_total(self):
        conditions, _ = views.condition_table(self.accepted_obs(), "tenv")
        chosen = views.solo_choice(conditions, rule="deepest", source="tenv")
        for position, row in enumerate(chosen):
            condition = conditions.iloc[position]
            totals = np.asarray(condition["library_totals"], dtype=np.int64)
            members = list(condition["members"])
            profile_ids = list(condition["profile_ids"])
            winner = sorted(range(len(members)),
                            key=lambda index: (-int(totals[index]), profile_ids[index]))[0]
            self.assertEqual(int(row), int(members[winner]))

    # -- product writer --------------------------------------------------- #
    def test_delta_product_is_an_alias_without_an_expression_payload(self):
        context = self.context("full-pooled-label")
        with views.modules()["view"].open_view(self.accepted) as release:
            with release.source("tenv") as reader:
                obs = reader.obs.copy()
                blocks = []
                for start in range(0, len(obs), views.BATCH):
                    stop = min(start + views.BATCH, len(obs))
                    blocks.append((start, {"delta": reader.X(start, stop)
                                           - reader.control(start, stop)}))
        rows = obs.copy()
        rows["profile_id"] = rows.index.astype(str)
        alias = {
            "release": "fixture-view-2",
            "manifest_sha256": views.sha256(self.accepted / "manifest.json"),
            "payload_sha256": {"tenv": self.accepted_manifest["outputs"]["tenv.h5ad"]["sha256"],
                               "tahoe_bulk": self.accepted_manifest["outputs"][
                                   "tahoe_bulk.h5ad"]["sha256"]},
            "control": {"file": "tenv.h5ad", "sha256": self.accepted_manifest["outputs"][
                "tenv.h5ad"]["sha256"]},
        }
        manifest = context.write_product(
            rows=rows, blocks=iter(blocks), mode="delta", alias=alias,
            recipe=views.RECIPE["full-pooled-label"],
            control_rule=views.CONTROL_RULE["full-pooled-label"], delta_kind="pooled",
            members=None, control_sources=None)
        staged = sorted(path.name for path in context.staging_dir.iterdir())
        self.assertEqual(staged, ["tenv.delta.npy", "tenv.json", "tenv.profile_ids.csv"])
        self.assertEqual(manifest["payload"]["file"], "tenv.delta.npy")
        self.assertEqual(manifest["payload"]["kind"], "delta")
        self.assertEqual(manifest["delta_kind"], "pooled")
        self.assertEqual(manifest["payload"]["layers"], ["delta"])
        self.assertIsNone(manifest["members"])
        self.assertEqual(manifest["alias"]["control"]["sha256"],
                         self.accepted_manifest["outputs"]["tenv.h5ad"]["sha256"])
        delta = np.load(context.staging_dir / "tenv.delta.npy")
        self.assertEqual(delta.shape, (len(rows), 2))
        self.assertEqual(delta.dtype, np.float32)
        views.modules()["release"].check_relative_paths(manifest)

    def test_payload_writer_requires_delta_to_equal_x_minus_control(self):
        context = self.context("merge")
        rows = pd.DataFrame({
            "profile_id": ["merge::tenv::hiPSC-CM::P1::1"],
            "role": ["treated"], "cell_line": ["hiPSC-CM"], "cell_line_id": [""],
            "dose_um": [1.0], "drug_parent_inchikey": ["P1"],
            "drug_parent_smiles": ["CCO"], "library_total": [10],
            "plate": ["Plate1"], "compound_name": ["compound"], "vehicle": [""],
            "source_item": ["tsg_view:merge"], "source_item_manifest_sha256": ["0" * 64],
            "source_row": [-1], "dose_original_value": [1.0],
            "dose_original_unit": ["uM"], "source_inchikey": [""], "n_cells": [0],
            "matched_control_group": ["tenv::Plate1::DMSO"]})
        x = np.full((1, 2), 0.5, dtype=np.float32)
        control = np.full((1, 2), 0.25, dtype=np.float32)
        good = [(0, {"X": x, "control": control})]
        context.write_product(rows=rows, blocks=iter(good), mode="payload",
                              recipe=views.RECIPE["merge"],
                              control_rule=views.CONTROL_RULE["merge"],
                              members=pd.DataFrame({"profile_id": rows["profile_id"],
                                                    "member_profile_id": ["tenv:p00"]}),
                              control_sources=pd.DataFrame({
                                  "profile_id": rows["profile_id"],
                                  "control_rule": ["pooled-sum"],
                                  "control_groups": ["tenv::Plate1::DMSO"],
                                  "pooled_retained_total": [12]}))
        with self.assertRaisesRegex(views.HarnessError, "delta must equal X - control"):
            fresh = views.Context(view="merge", source="tenv",
                                  output_dir=self.root / "out", staging_dir=self.root / "bad",
                                  repo_root_path=REPO, size=2, fixture=True,
                                  corpus_dir=self.fixture["corpus"],
                                  axis_dir=self.fixture["axis"],
                                  accepted_dir=self.accepted, provenance_root=None)
            fresh.write_product(rows=rows, mode="payload",
                                blocks=iter([(0, {"X": x, "control": control,
                                                  "delta": np.full((1, 2), 0.75,
                                                                   dtype=np.float32)})]),
                                recipe=views.RECIPE["merge"],
                                control_rule=views.CONTROL_RULE["merge"],
                                members=pd.DataFrame({"profile_id": rows["profile_id"],
                                                      "member_profile_id": ["tenv:p00"]}),
                                control_sources=pd.DataFrame({
                                    "profile_id": rows["profile_id"],
                                    "control_rule": ["pooled-sum"],
                                    "control_groups": ["tenv::Plate1::DMSO"],
                                    "pooled_retained_total": [12]}))

    def test_writer_rejects_a_mode_that_does_not_match_the_view(self):
        context = self.context("full")
        with self.assertRaisesRegex(views.HarnessError, "must publish in alias mode"):
            context.write_product(rows=pd.DataFrame({
                "profile_id": ["x"], "role": ["treated"], "cell_line": ["c"],
                "cell_line_id": [""], "dose_um": [1.0], "drug_parent_inchikey": ["P"],
                "drug_parent_smiles": ["C"], "library_total": [1]}),
                blocks=None, mode="payload", recipe="r", control_rule="c")


if __name__ == "__main__":
    unittest.main()
