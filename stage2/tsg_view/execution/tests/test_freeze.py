"""Freeze seam: a view product becomes a verified Step 04-shaped input freeze."""
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
EXECUTION = HERE.parent
for path in (str(EXECUTION), str(HERE)):
    if path not in sys.path:
        sys.path.insert(0, path)

import build  # noqa: E402
import freeze  # noqa: E402
import inputs  # noqa: E402
import tsg_view_fixtures as fixtures  # noqa: E402
import views  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
SIZE = 5000


class FreezeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        views.bind(REPO)

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.fixture = fixtures.dense_fixture(self.root)
        self.accepted, self.accepted_manifest, _ = fixtures.build_view(self.fixture, size=SIZE)
        self.population = fixtures.treated_population(self.accepted, "tenv")
        self.sealed = fixtures.write_sealed(self.root, "tenv", self.population)
        self.views_root = self.root / "views"
        self.genes = np.asarray(self.fixture["axes"][SIZE], dtype=np.int64)

    # -- helpers ---------------------------------------------------------- #
    def accepted_obs(self, source="tenv"):
        with views.modules()["view"].open_view(self.accepted) as release:
            with release.source(source) as reader:
                return reader.obs.copy()

    def product_context(self, view):
        return views.Context(view=view, source="tenv", output_dir=self.views_root,
                             staging_dir=self.root / f"staging-{view}",
                             repo_root_path=REPO, size=SIZE, fixture=True,
                             corpus_dir=self.fixture["corpus"],
                             axis_dir=self.fixture["axis"],
                             accepted_dir=self.accepted, provenance_root=None)

    def build_full(self):
        return build.build_product(
            view="full", source="tenv", output=self.views_root,
            staging=self.root / "staging-build", repo_root_path=REPO, size=SIZE,
            fixture=True, corpus=self.fixture["corpus"], axis=self.fixture["axis"],
            accepted=self.accepted, provenance_root=None)

    def build_merge(self):
        """Assemble the replicate-pooled product with the shared primitives only."""
        context = self.product_context("merge")
        obs = self.accepted_obs()
        conditions, row_to_condition = views.condition_table(obs, "tenv")
        corpus = views.modules()["corpus"].open_corpus(self.fixture["corpus"])
        reader = corpus.source("tenv")
        try:
            pooled = views.pooled_condition_counts(reader, self.genes, row_to_condition,
                                                   conditions)
            group_sums = views.vehicle_group_sums(reader, self.genes)
        finally:
            reader.close()
        control, retained = views.pooled_control(conditions, group_sums)
        rows = views.merge_rows(conditions, "tenv", retained,
                                context.source_manifests()["accepted_view_sha256"])
        members = views.member_table(
            conditions, lambda condition: views.merge_profile_id(
                "tenv", condition.cell_line, condition.drug_parent_inchikey,
                condition.dose_um))
        controls = views.control_source_table(conditions, rows["profile_id"].tolist(),
                                             views.CONTROL_RULE["merge"], retained)
        context.write_product(rows=rows, blocks=views.pooled_condition_blocks(pooled, control),
                              mode="payload", recipe=views.RECIPE["merge"],
                              control_rule=views.CONTROL_RULE["merge"], members=members,
                              control_sources=controls)
        context.install()
        return context, rows, pooled, control

    def freeze(self, destination=None, **kwargs):
        options = dict(view="full", source="tenv", views_root=self.views_root,
                       destination=destination or self.root / "freezes/full/tenv",
                       repo_root_path=REPO, sealed=self.sealed,
                       axis_dir=self.fixture["axis"], accepted_dir=self.accepted)
        options.update(kwargs)
        return freeze.freeze_product(**options)

    def verified(self, directory):
        return inputs.Inputs(directory, sealed=self.sealed, repo_root=REPO)

    def tamper(self, directory, name, edit):
        """Re-write one freeze file and its manifest digest, as a hostile edit would."""
        directory = Path(directory)
        copied = directory.parent / f"{directory.name}-tampered"
        if copied.exists():
            shutil.rmtree(copied)
        shutil.copytree(directory, copied)
        path = copied / name
        edit(path)
        manifest = json.loads((copied / "manifest.json").read_text())
        manifest["outputs"][name] = views.sha256(path)
        (copied / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        return copied

    # -- full view freeze ------------------------------------------------- #
    def test_full_freeze_is_byte_equal_to_the_sealed_splits(self):
        self.build_full()
        receipt = self.freeze()
        directory = self.root / "freezes/full/tenv"
        self.assertTrue(receipt["verified"])
        self.assertTrue(receipt["splits_byte_equal_to_sealed"])
        self.assertEqual(receipt["splits_sha256"], views.sha256(self.sealed / "splits.csv"))
        self.assertEqual(sorted(path.name for path in directory.iterdir()),
                         sorted(inputs.OUTPUTS + ("manifest.json",)))
        data = self.verified(directory)
        self.assertEqual(len(data.population), len(self.population))
        self.assertEqual(data.config["delta_kind"], "expression_minus_control")
        self.assertEqual(data.config["n_conditions"], 9)
        self.assertEqual(data.expression.shape, (len(self.population), SIZE))
        self.assertEqual(data.control.shape, (len(self.population), SIZE))
        self.assertEqual(len(data.all_obs), len(self.population) + 4)
        self.assertTrue(data.controls_sources["pooled_retained_total"].isna().all())
        self.assertTrue(data.controls_sources["control_rule"].eq(
            "mean-of-normalized-vehicles").all())
        self.assertEqual(views.sha256(directory / "morgan.npy"),
                         views.sha256(self.sealed / "morgan.npy"))
        # Every member of a condition shares one partition in every scenario x fold.
        keys = ["cell_line", "drug_parent_inchikey", "dose_um"]
        grouped = data.splits.assign(
            condition=data.population.reset_index(drop=True)[keys].astype(str).agg(
                "|".join, axis=1).to_numpy()[data.splits["row"]]).groupby(
            ["scenario", "fold", "condition"])["partition"].nunique()
        self.assertTrue((grouped == 1).all())

    def test_full_freeze_rejects_a_changed_partition(self):
        self.build_full()
        self.freeze()
        directory = self.root / "freezes/full/tenv"

        sealed = pd.read_csv(self.sealed / "splits.csv")
        row = sealed[sealed.partition == "test"].iloc[0]
        target = (row.scenario, row.fold, row.drug_parent, row.dose_um)

        def edit(path):
            splits = pd.read_csv(path)
            mask = ((splits.scenario == target[0]) & (splits.fold == target[1])
                    & (splits.drug_parent == target[2]) & (splits.dose_um == target[3]))
            self.assertTrue(mask.any())
            splits.loc[mask, "partition"] = "train"
            splits.to_csv(path, index=False)

        tampered = self.tamper(directory, "splits.csv", edit)
        with self.assertRaisesRegex(views.HarnessError, "changed partition"):
            self.verified(tampered)

    def test_full_freeze_rejects_changed_parents_and_genes(self):
        self.build_full()
        self.freeze()
        directory = self.root / "freezes/full/tenv"
        def edit_parents(path):
            parents = pd.read_csv(path)
            parents.loc[0, "drug_parent_smiles"] = "CCC"
            parents.to_csv(path, index=False)

        with self.assertRaisesRegex(views.HarnessError, "parents.csv differs"):
            self.verified(self.tamper(directory, "parents.csv", edit_parents))
        with self.assertRaisesRegex(views.HarnessError, "genes.csv differs"):
            self.verified(self.tamper(directory, "genes.csv",
                                      lambda path: path.write_text(
                                          path.read_text().replace("S", "T", 1))))

    def test_verifier_rejects_a_sealed_reference_that_changed(self):
        self.build_full()
        self.freeze()
        directory = self.root / "freezes/full/tenv"
        sealed = self.root / "sealed/tenv"
        splits = sealed / "splits.csv"
        original = splits.read_bytes()
        splits.write_bytes(original + b"\n")
        try:
            with self.assertRaisesRegex(views.HarnessError,
                                        "sealed Step 04 freeze changed"):
                self.verified(directory)
        finally:
            splits.write_bytes(original)

    def test_freeze_refuses_a_non_empty_destination(self):
        self.build_full()
        destination = self.root / "freezes/full/tenv"
        destination.mkdir(parents=True)
        (destination / "stale.txt").write_text("stale\n")
        with self.assertRaisesRegex(views.HarnessError, "new empty candidate"):
            self.freeze(destination=destination)

    # -- derived view freeze ---------------------------------------------- #
    def test_merge_freeze_has_one_row_per_condition_with_pooled_control(self):
        self.build_merge()
        receipt = self.freeze(view="merge", destination=self.root / "freezes/merge/tenv")
        self.assertTrue(receipt["verified"])
        directory = self.root / "freezes/merge/tenv"
        data = self.verified(directory)
        self.assertEqual(len(data.population), 9)
        self.assertEqual(len(data.all_obs), 9 + 4)
        self.assertEqual(data.expression.shape, (9, SIZE))
        self.assertEqual(float(data.population["dose_um"].nunique()), 3.0)
        self.assertTrue(data.controls_sources["control_rule"].eq("pooled-sum").all())
        self.assertTrue(data.population["profile_id"].str.startswith("merge::tenv::").all())
        # Splits stay at the sealed condition level: 9 conditions x 6 blocks.
        self.assertEqual(len(data.splits), 9 * 6)
        sealed = inputs.partition_map(pd.read_csv(self.sealed / "splits.csv"))
        proven = inputs.partition_map(data.splits)
        joined = sealed.merge(proven, on=["scenario", "fold", "drug_parent", "dose_um"],
                              suffixes=("_sealed", "_view"))
        self.assertTrue((joined.partition_sealed == joined.partition_view).all())
        # Pooled control is a normalize-of-sum composition, not a mean of normalized rows.
        recovered = np.expm1(np.asarray(data.control, dtype=np.float64)).sum(axis=1)
        np.testing.assert_allclose(recovered, 10000.0, rtol=1e-5, atol=5.0)


if __name__ == "__main__":
    unittest.main()
