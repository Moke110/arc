"""Build/verify seam: `solo` and `solo-deepest` keep one full-view profile per condition.

The fixture is the shared multi-condition grid: nine tenv conditions with two singleton and
seven multi-member conditions, so a draw, a deepest pick and the singleton-equality rule are
all observable without touching real payloads.
"""
import json
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
import tsg_view_fixtures as fixtures  # noqa: E402
import views  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
CONDITION_KEYS = ["cell_line", "drug_parent_inchikey", "dose_um"]


class SoloFixture(unittest.TestCase):
    """Shared fixture view and build helpers for the two member-choice views."""

    @classmethod
    def setUpClass(cls):
        views.bind(REPO)

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.fixture = fixtures.grid_fixture(self.root)
        self.accepted, self.accepted_manifest, _ = fixtures.build_view(
            self.fixture, size=2)

    # -- helpers ---------------------------------------------------------- #
    def build(self, view="solo", source="tenv", output=None):
        return build.build_product(
            view=view, source=source, output=output or self.root / "views",
            staging=self.root / "staging", repo_root_path=REPO, size=2, fixture=True,
            corpus=self.fixture["corpus"], axis=self.fixture["axis"],
            accepted=self.accepted, provenance_root=None)

    def product(self, view, source="tenv"):
        directory = self.root / "views" / view
        manifest = json.loads((directory / f"{source}.json").read_text())
        payload = directory / manifest["payload"]["file"]
        with h5py.File(payload, "r") as handle:
            obs = read_elem(handle["obs"])
            x = np.asarray(handle["X"], dtype=np.float32)
            control = np.asarray(handle["layers/control"], dtype=np.float32)
            delta = np.asarray(handle["layers/delta"], dtype=np.float32)
        obs = pd.DataFrame(obs)
        obs["profile_id"] = obs["profile_id"].astype(str)
        return manifest, obs.reset_index(drop=True), x, control, delta

    def accepted_view(self, source="tenv"):
        """The accepted fixture view's obs, X and control, indexed by profile id."""
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
        return obs, x, control

    def members(self, manifest):
        directory = self.root / "views" / manifest["view"]
        table = pd.read_csv(directory / manifest["members"]["file"],
                            dtype={"profile_id": str, "member_profile_id": str})
        return table

class SoloViewTest(SoloFixture):
    def test_solo_keeps_exactly_one_full_view_row_per_condition(self):
        receipt = self.build("solo")
        self.assertEqual(receipt["view"], "solo")
        self.assertEqual(receipt["rows"], 9)
        self.assertEqual(receipt["treated"], 9)
        self.assertEqual(receipt["conditions"], 9)
        self.assertEqual(receipt["mode"], "payload")
        manifest, obs, x, control, delta = self.product("solo")
        self.assertEqual(manifest["schema"], "tsg-view-product-v1")
        self.assertEqual(manifest["product"], "tenv-solo-5k-view")
        self.assertEqual(manifest["rows"], 9)
        self.assertEqual(manifest["conditions"], 9)
        self.assertEqual(manifest["seed"], views.SOLO_SEED)
        self.assertIn("drawn", manifest["selection"])
        self.assertEqual(manifest["payload"]["layers"], ["control", "delta"])
        self.assertEqual(sorted(manifest["payload"].keys() - {"file", "sha256", "bytes",
                                                             "rows", "cols", "dtype",
                                                             "layers"}),
                         [])
        # One row per condition, so nine unique condition keys.
        keys = np.asarray([views.condition_key(row.cell_line, row.drug_parent_inchikey,
                                               row.dose_um)
                           for row in obs.itertuples()])
        self.assertEqual(len(set(keys.tolist())), 9)
        self.assertEqual(len(obs), 9)
        # The chosen row is an accepted full-view row with its own identity and values.
        accepted_obs, accepted_x, accepted_control = self.accepted_view()
        position = {profile: row for row, profile in enumerate(accepted_obs["profile_id"])}
        for row in range(len(obs)):
            profile = obs.loc[row, "profile_id"]
            self.assertIn(profile, position)
            source_row = position[profile]
            np.testing.assert_array_equal(x[row], accepted_x[source_row])
            np.testing.assert_array_equal(control[row], accepted_control[source_row])
            self.assertEqual(int(obs.loc[row, "library_total"]),
                             int(accepted_obs.iloc[source_row]["library_total"]))
            self.assertEqual(obs.loc[row, "matched_control_group"],
                             accepted_obs.iloc[source_row]["matched_control_group"])
            np.testing.assert_array_equal(delta[row], x[row] - control[row])
        # Every product row is a member of its own condition, and the member table lists
        # all 16 treated fixture profiles under their chosen identity.
        conditions, _ = views.condition_table(accepted_obs, "tenv")
        member_of = {}
        for condition in conditions.itertuples():
            for member in condition.members:
                member_of[int(member)] = str(condition.condition_key)
        table = self.members(manifest)
        self.assertEqual(len(table), 16)
        self.assertEqual(sorted(table["profile_id"].unique()), sorted(obs["profile_id"]))
        for profile, condition_key in zip(obs["profile_id"], keys):
            self.assertEqual(member_of[position[profile]], condition_key)
        # The manifest reports the draw and the discarded total.
        self.assertEqual(receipt["evidence"]["discarded_members"], 16 - 9)
        self.assertEqual(receipt["evidence"]["source_conditions"], 9)
        views.modules()["release"].check_relative_paths(manifest)

    def test_solo_draw_is_stable_across_runs(self):
        self.build("solo", output=self.root / "views-a")
        self.build("solo", output=self.root / "views-b")
        first = json.loads((self.root / "views-a" / "solo" / "tenv.json").read_text())
        second = json.loads((self.root / "views-b" / "solo" / "tenv.json").read_text())
        self.assertEqual(first, second)
        self.assertEqual(first["payload"]["sha256"], second["payload"]["sha256"])
        self.assertEqual(
            views.sha256(self.root / "views-a" / "solo" / "tenv.h5ad"),
            views.sha256(self.root / "views-b" / "solo" / "tenv.h5ad"))


if __name__ == "__main__":
    unittest.main()
