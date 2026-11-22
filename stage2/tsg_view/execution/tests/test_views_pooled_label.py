"""Build/freeze seam: `full-pooled-label` keeps full inputs and replaces only the label.

The fixture grid carries nine conditions with singleton, same-plate and cross-plate repeats,
so a pooled label can be told apart from each member's own ``X - control``: the view must
keep every full-view row (treated and vehicle), keep the full expression and control bytes,
and give every treated row of one condition the same Δ — that condition's `merge` Δ.
"""
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

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
CONDITION_KEYS = ["cell_line", "drug_parent_inchikey", "dose_um"]


def read_product(path):
    """Read a materialized product's obs and layers straight from its declared layout."""
    with h5py.File(path, "r") as handle:
        obs = read_elem(handle["obs"]).reset_index(drop=True)
        return (obs,
                np.asarray(handle["X"], dtype=np.float32),
                np.asarray(handle["layers/control"], dtype=np.float32),
                np.asarray(handle["layers/delta"], dtype=np.float32))


class PooledLabelBuildTest(unittest.TestCase):
    """The build seam: alias the full rows, transport only the pooled Δ."""

    @classmethod
    def setUpClass(cls):
        views.bind(REPO)

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.fixture = fixtures.grid_fixture(self.root)
        self.accepted, self.accepted_manifest, _ = fixtures.build_view(self.fixture, size=2)
        self.output = self.root / "views"

    # -- helpers ---------------------------------------------------------- #
    def build(self, view, source="tenv", output=None):
        return build.build_product(
            view=view, source=source, output=output or self.output,
            staging=self.root / "staging", repo_root_path=REPO, size=2, fixture=True,
            corpus=self.fixture["corpus"], axis=self.fixture["axis"],
            accepted=self.accepted, provenance_root=None)

    def prepare(self, source="tenv"):
        """The two products `full-pooled-label` consumes: `full` and `merge`."""
        self.build("full", source)
        self.build("merge", source)

    def manifest(self, view, source="tenv"):
        return json.loads((self.output / view / f"{source}.json").read_text())

    def delta(self, source="tenv"):
        manifest = self.manifest("full-pooled-label", source)
        return np.load(self.output / "full-pooled-label" / manifest["payload"]["file"])

    def merge_product(self, source="tenv"):
        manifest = self.manifest("merge", source)
        return read_product(self.output / "merge" / manifest["payload"]["file"])

    def accepted_vector(self, source="tenv"):
        """The accepted fixture view's full-row obs, X and control."""
        view = views.modules()["view"]
        with view.open_view(self.accepted) as release:
            reader = release.source(source)
            try:
                obs = reader.obs.copy()
                x = np.asarray(reader.X(0, len(obs)), dtype=np.float32)
                control = np.asarray(reader.control(0, len(obs)), dtype=np.float32)
            finally:
                reader.close()
        obs["profile_id"] = obs.index.astype(str)
        return obs.reset_index(drop=True), x, control

    def merge_row_per_row(self, obs, source, merge_obs):
        index = {str(pid): position
                 for position, pid in enumerate(merge_obs["profile_id"])}
        targets = np.full(len(obs), -1, dtype=np.int64)
        for position, row in enumerate(obs.itertuples()):
            if row.role != "treated":
                continue
            profile = views.merge_profile_id(source, row.cell_line,
                                             row.drug_parent_inchikey, row.dose_um)
            targets[position] = index[profile]
        return targets

    # -- expression / control identity to `full` --------------------------- #
    def test_rows_expression_and_control_alias_the_full_product(self):
        self.prepare("tenv")
        receipt = self.build("full-pooled-label", "tenv")
        full = self.manifest("full", "tenv")
        manifest = self.manifest("full-pooled-label", "tenv")
        obs, x, control = self.accepted_vector("tenv")

        self.assertEqual(receipt["view"], "full-pooled-label")
        self.assertEqual(receipt["mode"], "delta")
        # Same row space as `full`: every treated and vehicle row, same identity and order.
        self.assertEqual(receipt["rows"], full["rows"])
        self.assertEqual(receipt["treated"], full["treated"])
        self.assertEqual(receipt["rows"], len(obs))
        self.assertEqual(manifest["rows"], full["rows"])
        self.assertEqual(manifest["profile_ids"]["sha256"],
                         full["profile_ids"]["sha256"])
        self.assertEqual(manifest["profile_ids"]["bytes"], full["profile_ids"]["bytes"])
        # X and control are never copied: they are the accepted payload's digests.
        self.assertEqual(manifest["alias"]["payload_sha256"],
                         full["alias"]["payload_sha256"])
        self.assertEqual(manifest["alias"]["manifest_sha256"],
                         full["alias"]["manifest_sha256"])
        self.assertEqual(manifest["alias"]["payload_sha256"]["tenv"],
                         views.sha256(self.accepted / "tenv.h5ad"))
        self.assertEqual(manifest["alias"]["manifest_sha256"],
                         views.sha256(self.accepted / "manifest.json"))
        # The delta payload is row-aligned to the full rows and the frozen axis.
        self.assertEqual(manifest["payload"]["kind"], "delta")
        self.assertEqual(manifest["payload"]["layers"], ["delta"])
        self.assertEqual(manifest["payload"]["rows"], len(obs))
        self.assertEqual(manifest["payload"]["cols"], 2)
        self.assertEqual(manifest["payload"]["dtype"], "float32")
        # Vehicle rows have no condition: their Δ stays the full view's X - control.
        vehicle = np.flatnonzero(obs["role"].to_numpy() == "vehicle")
        self.assertGreater(len(vehicle), 0)
        delta = self.delta()
        np.testing.assert_array_equal(delta[vehicle], (x - control)[vehicle])

    # -- the label is the merge Δ of the row's own condition ---------------- #
    def test_every_treated_row_equals_its_conditions_merge_delta(self):
        self.prepare("tenv")
        self.build("full-pooled-label", "tenv")
        obs, _, _ = self.accepted_vector("tenv")
        merge_obs, _, _, merge_delta = self.merge_product("tenv")
        delta = self.delta()
        targets = self.merge_row_per_row(obs, "tenv", merge_obs)
        treated = np.flatnonzero(obs["role"].to_numpy() == "treated")
        self.assertEqual(int((targets[treated] >= 0).sum()), len(treated))
        for row in treated:
            np.testing.assert_array_equal(delta[row], merge_delta[targets[row]],
                                          err_msg=f"row {row}")

    def test_rows_sharing_a_condition_share_one_pooled_delta(self):
        self.prepare("tenv")
        self.build("full-pooled-label", "tenv")
        obs, x, control = self.accepted_vector("tenv")
        delta = self.delta()
        treated = obs[obs["role"] == "treated"]
        ids = obs["profile_id"].to_numpy()
        grouped = list(treated.groupby(CONDITION_KEYS, sort=True))
        self.assertEqual(len(grouped), 9)
        repeated = 0
        for _, members in grouped:
            positions = [int(np.flatnonzero(ids == pid)[0])
                         for pid in members["profile_id"]]
            self.assertGreaterEqual(len(positions), 1)
            for position in positions[1:]:
                np.testing.assert_array_equal(delta[position], delta[positions[0]])
            if len(positions) > 1:
                repeated += 1
        self.assertEqual(repeated, 7)
        # The pooled label is not each row's own X - control on at least one repeated row.
        own = x - control
        self.assertTrue(any(not np.array_equal(delta[row], own[row])
                            for row in treated.index))
        # ... and it equals the merge Δ, not the members' mean Δ.
        merge_obs, _, _, merge_delta = self.merge_product("tenv")
        targets = self.merge_row_per_row(obs, "tenv", merge_obs)
        rows = np.flatnonzero(ids == "tenv:p00")[0].item()
        member_rows = np.flatnonzero(ids == "tenv:p01")
        self.assertTrue(len(member_rows))
        mean = ((own[rows] + own[member_rows[0]]) / 2.0).astype(np.float32)
        self.assertTrue(np.array_equal(delta[rows], merge_delta[targets[rows]]))
        self.assertFalse(np.array_equal(delta[rows], mean))

    # -- manifest provenance ------------------------------------------------ #
    def test_manifest_records_full_merge_and_per_condition_counts(self):
        self.prepare("tenv")
        self.build("full-pooled-label", "tenv")
        manifest = self.manifest("full-pooled-label", "tenv")
        full = self.manifest("full", "tenv")
        merge = self.manifest("merge", "tenv")
        merge_path = self.output / "merge" / "tenv.json"
        merge_payload = self.output / "merge" / merge["payload"]["file"]
        obs, _, _ = self.accepted_vector("tenv")

        self.assertEqual(manifest["schema"], "tsg-view-product-v1")
        self.assertEqual(manifest["product"], "tenv-full-pooled-label-5k-view")
        self.assertEqual(manifest["delta_kind"], "pooled")
        self.assertEqual(manifest["control_rule"],
                         views.CONTROL_RULE["full-pooled-label"])
        self.assertEqual(manifest["recipe"], views.RECIPE["full-pooled-label"])

        # The full provenance: the `full` product manifest and accepted payload digests.
        self.assertEqual(manifest["full"]["product"], full["product"])
        self.assertEqual(manifest["full"]["manifest"]["sha256"],
                         views.sha256(self.output / "full" / "tenv.json"))
        self.assertEqual(manifest["full"]["payload"]["sha256"],
                         full["alias"]["payload_sha256"]["tenv"])
        self.assertEqual(manifest["full"]["rows"], full["rows"])

        # The merge Δ provenance: the exact payload/manifest bytes and the pooling rule.
        label = manifest["pooled_label"]["label_source"]
        self.assertEqual(label["view"], "merge")
        self.assertEqual(label["product"], merge["product"])
        self.assertEqual(label["manifest"]["sha256"], views.sha256(merge_path))
        self.assertEqual(label["payload"]["sha256"], views.sha256(merge_payload))
        self.assertEqual(label["payload"]["sha256"], merge["payload"]["sha256"])
        self.assertEqual(label["payload"]["layers"], ["control", "delta"])
        self.assertEqual(label["layer"], "delta")
        self.assertEqual(label["recipe"], merge["recipe"])
        self.assertEqual(label["control_rule"], merge["control_rule"])
        self.assertEqual(manifest["pooled_label"]["delta_kind"], "pooled")

        # Per-condition member counts reconcile to the treated rows they label.
        histogram = {key: value for key, value in
                     manifest["pooled_label"]["member_count_histogram"].items()
                     if value}
        expected_histogram = {key: value for key, value in
                              merge["special_conditions"]["member_count_histogram"].items()
                              if value}
        self.assertEqual(histogram, expected_histogram)
        self.assertEqual(sum(int(key) * value for key, value in histogram.items()),
                         manifest["treated"])
        table = manifest["pooled_label"]["conditions_table"]
        self.assertEqual(table["rows"], 9)
        self.assertEqual(table["sha256"],
                         views.sha256(self.output / "full-pooled-label" / table["file"]))
        frame = pd.read_csv(self.output / "full-pooled-label" / table["file"])
        self.assertEqual(list(frame.columns),
                         ["condition_key", "merge_profile_id", "cell_line",
                          "drug_parent_inchikey", "dose_um", "member_count",
                          "member_profile_ids", "profile_ids"])
        self.assertEqual(int(frame["member_count"].sum()), manifest["treated"])
        self.assertEqual(set(frame["merge_profile_id"]),
                         set(self.merge_product("tenv")[0]["profile_id"]))
        # The members table maps each treated product row to its condition's members.
        members_table = pd.read_csv(self.output / "full-pooled-label"
                                    / manifest["members"]["file"])
        self.assertEqual(list(members_table.columns), ["profile_id", "member_profile_id"])
        treated_ids = set(obs[obs["role"] == "treated"]["profile_id"])
        self.assertEqual(set(members_table["profile_id"]), treated_ids)
        counts = members_table.groupby("profile_id").size()
        # Each treated row names every member of its condition: sum = sum(member_count**2).
        self.assertEqual(int(counts.sum()),
                         int((frame["member_count"] ** 2).sum()))
        self.assertEqual(int(counts.max()), 2)

    def test_delta_payload_is_deterministic(self):
        self.prepare("tenv")
        first = self.build("full-pooled-label", "tenv")
        before = first["files"]["tenv.delta.npy"]["sha256"]
        second = self.build("full-pooled-label", "tenv")
        self.assertEqual(second["files"]["tenv.delta.npy"]["sha256"], before)
        self.assertEqual(second["manifest"]["sha256"], first["manifest"]["sha256"])

    # -- fail closed -------------------------------------------------------- #
    def test_build_requires_the_merge_product(self):
        self.build("full", "tenv")
        with self.assertRaisesRegex(views.HarnessError, "merge product"):
            self.build("full-pooled-label", "tenv")

    def test_build_rejects_a_tampered_merge_payload(self):
        self.prepare("tenv")
        tampered = self.root / "views-tampered"
        shutil.copytree(self.output, tampered)
        with (tampered / "merge" / "tenv.h5ad").open("ab") as handle:
            handle.write(b"tamper")
        with self.assertRaisesRegex(views.HarnessError, "merge payload digest"):
            self.build("full-pooled-label", "tenv", output=tampered)


class PooledLabelFreezeTest(unittest.TestCase):
    """The freeze seam: `expression = full control + pooled Δ`, sealed splits unchanged."""

    @classmethod
    def setUpClass(cls):
        views.bind(REPO)

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.fixture = fixtures.dense_fixture(self.root)
        self.accepted, self.accepted_manifest, _ = fixtures.build_view(
            self.fixture, size=SIZE)
        self.output = self.root / "views"
        self.population = fixtures.treated_population(self.accepted, "tenv")
        self.sealed = fixtures.write_sealed(self.root, "tenv", self.population)
        for view in ("full", "merge", "full-pooled-label"):
            self.build(view)

    def build(self, view, source="tenv"):
        return build.build_product(
            view=view, source=source, output=self.output,
            staging=self.root / f"staging-{view}", repo_root_path=REPO, size=SIZE,
            fixture=True, corpus=self.fixture["corpus"], axis=self.fixture["axis"],
            accepted=self.accepted, provenance_root=None)

    def freeze(self, view, destination):
        return freeze.freeze_product(
            view=view, source="tenv", views_root=self.output, destination=destination,
            repo_root_path=REPO, sealed=self.sealed, axis_dir=self.fixture["axis"],
            accepted_dir=self.accepted)

    def accepted_vector(self):
        view = views.modules()["view"]
        with view.open_view(self.accepted) as release:
            with release.source("tenv") as reader:
                obs = reader.obs.copy()
                x = np.asarray(reader.X(0, len(obs)), dtype=np.float32)
                control = np.asarray(reader.control(0, len(obs)), dtype=np.float32)
        obs["profile_id"] = obs.index.astype(str)
        return obs.reset_index(drop=True), x, control

    def test_freeze_pools_the_training_target_and_keeps_the_sealed_splits(self):
        full = self.freeze("full", self.root / "freezes/full/tenv")
        receipt = self.freeze("full-pooled-label",
                              self.root / "freezes/full-pooled-label/tenv")
        directory = self.root / "freezes/full-pooled-label/tenv"
        self.assertTrue(receipt["verified"])
        self.assertEqual(receipt["delta_kind"], "pooled")
        self.assertTrue(receipt["splits_byte_equal_to_sealed"])
        self.assertEqual(receipt["splits_sha256"], views.sha256(self.sealed / "splits.csv"))
        # Row count and splits are the `full` view's, byte for byte.
        self.assertEqual(receipt["rows"], full["rows"])
        self.assertEqual(receipt["rows"], len(self.population))
        self.assertEqual(receipt["splits_sha256"], full["splits_sha256"])

        data = inputs.Inputs(directory, sealed=self.sealed, repo_root=REPO)
        self.assertEqual(data.config["delta_kind"], "pooled")
        self.assertEqual(data.config["n_treated"], len(self.population))
        self.assertEqual(len(data.population), len(self.population))
        self.assertEqual(data.expression.shape, (len(self.population), SIZE))
        self.assertEqual(data.control.shape, (len(self.population), SIZE))

        obs, x, control = self.accepted_vector()
        treated = np.flatnonzero(obs["role"].to_numpy() == "treated")
        delta_manifest = json.loads(
            (self.output / "full-pooled-label" / "tenv.json").read_text())
        pooled = np.load(self.output / "full-pooled-label"
                         / delta_manifest["payload"]["file"])[treated]
        # control is the full view's control; expression is control + the pooled label.
        np.testing.assert_array_equal(np.asarray(data.control), control[treated])
        expected = (control[treated].astype(np.float64)
                    + pooled.astype(np.float64)).astype(np.float32)
        np.testing.assert_array_equal(np.asarray(data.expression), expected)
        # The training target's Δ is the pooled label, not each row's own X - control.
        recovered = (np.asarray(data.expression, dtype=np.float64)
                     - np.asarray(data.control, dtype=np.float64)).astype(np.float32)
        # The float32 round-trip of `expression` perturbs Δ by ~1e-7 (plan C7.10), so the
        # frozen target is asserted at that tolerance, never by array equality.
        np.testing.assert_allclose(recovered, pooled, rtol=0, atol=1e-6)
        own = (x - control)[treated]
        self.assertTrue(any(not np.array_equal(a, b) for a, b in zip(recovered, own)))
        # Every condition's rows still share one partition in every scenario x fold.
        keys = ["cell_line", "drug_parent_inchikey", "dose_um"]
        frame = data.population.reset_index(drop=True)[keys].astype(str).agg(
            "|".join, axis=1).to_numpy()
        grouped = data.splits.assign(condition=frame[data.splits["row"]]).groupby(
            ["scenario", "fold", "condition"])["partition"].nunique()
        self.assertTrue((grouped == 1).all())


if __name__ == "__main__":
    unittest.main()
