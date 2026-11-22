"""Build/verify seam: `solo` and `solo-deepest` keep one full-view profile per condition.

The fixture is the shared multi-condition grid: nine tenv conditions with two singleton and
seven multi-member conditions, so a draw, a deepest pick and the singleton-equality rule are
all observable without touching real payloads.
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


class SoloDeepestViewTest(SoloFixture):
    def test_solo_deepest_takes_the_largest_retained_total(self):
        receipt = self.build("solo-deepest")
        self.assertEqual(receipt["rows"], 9)
        self.assertEqual(receipt["conditions"], 9)
        manifest, obs, x, control, delta = self.product("solo-deepest")
        self.assertIsNone(manifest["seed"])
        self.assertIn("largest retained selected-axis", manifest["selection"])
        accepted_obs, accepted_x, accepted_control = self.accepted_view()
        conditions, _ = views.condition_table(accepted_obs, "tenv")
        position = {profile: row for row, profile in enumerate(accepted_obs["profile_id"])}
        for row in range(len(obs)):
            key = views.condition_key(obs.loc[row, "cell_line"],
                                      obs.loc[row, "drug_parent_inchikey"],
                                      obs.loc[row, "dose_um"])
            condition = conditions.loc[conditions["condition_key"] == key].iloc[0]
            totals = np.asarray(condition["library_totals"], dtype=np.int64)
            members = [int(member) for member in condition["members"]]
            profile_ids = [str(profile) for profile in condition["profile_ids"]]
            winner = sorted(range(len(members)),
                            key=lambda index: (-int(totals[index]), profile_ids[index]))[0]
            self.assertEqual(obs.loc[row, "profile_id"], profile_ids[winner])
            self.assertEqual(int(obs.loc[row, "library_total"]), int(totals[winner]))
            np.testing.assert_array_equal(
                x[row], accepted_x[position[profile_ids[winner]]])
            np.testing.assert_array_equal(
                control[row], accepted_control[position[profile_ids[winner]]])
            np.testing.assert_array_equal(delta[row], x[row] - control[row])
        self.assertEqual(receipt["evidence"]["discarded_members"], 16 - 9)

    def test_solo_deepest_equals_solo_on_single_member_conditions(self):
        self.build("solo")
        self.build("solo-deepest")
        solo = pd.read_csv(self.root / "views/solo/tenv.profile_ids.csv",
                           dtype={"profile_id": str})
        deepest = pd.read_csv(self.root / "views/solo-deepest/tenv.profile_ids.csv",
                              dtype={"profile_id": str})
        _, solo_obs, _, _, _ = self.product("solo")
        _, deepest_obs, _, _, _ = self.product("solo-deepest")
        accepted_obs, _, _ = self.accepted_view()
        conditions, _ = views.condition_table(accepted_obs, "tenv")
        singles = {views.condition_key(row.cell_line, row.drug_parent_inchikey, row.dose_um)
                   for row in conditions[conditions["member_count"] == 1].itertuples()}
        self.assertEqual(len(singles), 2)
        chosen = {}
        for name, obs in (("solo", solo_obs), ("solo-deepest", deepest_obs)):
            chosen[name] = {views.condition_key(row.cell_line, row.drug_parent_inchikey,
                                                row.dose_um): row.profile_id
                            for row in obs.itertuples()}
        for key in singles:
            self.assertEqual(chosen["solo"][key], chosen["solo-deepest"][key])
        self.assertEqual(len(solo), len(deepest))
        self.assertEqual(sorted(solo["row"]), list(range(9)))

    def test_deepest_tie_breaks_to_the_smallest_profile_id(self):
        # Two members with the same retained total: the smaller profile_id must win,
        # independent of the member order the condition table happens to carry.
        frame = pd.DataFrame([{
            "condition_key": "L|P|1", "members": np.asarray([3, 7], dtype=np.int64),
            "profile_ids": np.asarray(["z-last", "a-first"]), "member_count": 2,
            "library_totals": np.asarray([100, 100], dtype=np.int64)}])
        chosen = views.solo_choice(frame, rule="deepest", source="tenv")
        self.assertEqual([int(row) for row in chosen], [7])
        reversed_frame = pd.DataFrame([{
            "condition_key": "L|P|1", "members": np.asarray([7, 3], dtype=np.int64),
            "profile_ids": np.asarray(["a-first", "z-last"]), "member_count": 2,
            "library_totals": np.asarray([100, 100], dtype=np.int64)}])
        self.assertEqual(
            [int(row) for row in views.solo_choice(reversed_frame, rule="deepest",
                                                   source="tenv")], [7])


class SoloChoiceIndependenceTest(SoloFixture):
    """The draw is a function of the condition table and the frozen seed alone."""

    def chosen_ids(self, conditions):
        chosen = views.solo_choice(conditions, rule="draw", source="tenv")
        ids = []
        for position, row in enumerate(chosen):
            condition = conditions.iloc[position]
            members = [int(member) for member in condition["members"]]
            ids.append(str(condition["profile_ids"][members.index(int(row))]))
        return ids

    def test_draw_is_independent_of_fold_layout_and_the_split_assignment(self):
        accepted_obs, _, _ = self.accepted_view()
        conditions, _ = views.condition_table(accepted_obs, "tenv")
        expected = self.chosen_ids(conditions)
        rng = np.random.default_rng(5)
        # A different fold layout and a different split assignment are columns the choice
        # never reads; the same held-parent/dose grid with arbitrary partitions is enough
        # to prove the chosen member cannot move with them.
        for scramble in (0, 1):
            layout = accepted_obs.copy()
            layout["scenario"] = "cold-drug-dose"
            layout["fold"] = rng.integers(0, 3, len(layout))
            layout["partition"] = rng.choice(
                ["train", "validation", "test", "excluded-held-parent",
                 "excluded-held-dose"], len(layout))
            if scramble:
                layout = layout.sample(frac=1.0, random_state=11)
            relaid, _ = views.condition_table(layout, "tenv")
            self.assertEqual(self.chosen_ids(relaid), expected)
        # The draw is a real draw, not "first member in profile order".
        first_member = []
        for condition in conditions.itertuples():
            first_member.append(str(condition.profile_ids[0]))
        self.assertNotEqual(expected, first_member)

    def test_training_seed_is_not_an_input_to_the_choice(self):
        accepted_obs, _, _ = self.accepted_view()
        conditions, _ = views.condition_table(accepted_obs, "tenv")
        chosen = views.solo_choice(conditions, rule="draw", source="tenv")
        # The only seed the interface accepts is the frozen draw seed; a different draw
        # seed moves the choice, which proves the seed (not the fit) owns it.
        other = views.solo_choice(conditions, rule="draw", seed=1, source="tenv")
        self.assertFalse(np.array_equal(chosen, other))
        self.assertEqual([int(row) for row in views.solo_choice(conditions, rule="draw",
                                                                seed=views.SOLO_SEED,
                                                                source="tenv")],
                         [int(row) for row in chosen])


class SoloManifestReconciliationTest(SoloFixture):
    def test_manifest_members_and_retained_totals_reconcile_to_the_corpus(self):
        receipt = self.build("solo")
        manifest, obs, _, _, _ = self.product("solo")
        accepted_obs, accepted_x, accepted_control = self.accepted_view()
        conditions, _ = views.condition_table(accepted_obs, "tenv")
        # profile_ids.csv is the product row space, in payload order.
        profile_ids = pd.read_csv(self.root / "views/solo/tenv.profile_ids.csv",
                                  dtype={"profile_id": str})
        self.assertEqual(list(profile_ids["row"]), list(range(manifest["rows"])))
        self.assertEqual(list(profile_ids["profile_id"]), list(obs["profile_id"]))
        self.assertEqual(manifest["profile_ids"]["bytes"],
                         (self.root / "views/solo/tenv.profile_ids.csv").stat().st_size)
        self.assertEqual(manifest["profile_ids"]["sha256"],
                         views.sha256(self.root / "views/solo/tenv.profile_ids.csv"))
        # The member table lists every treated fixture profile under its chosen identity.
        table = self.members(manifest)
        self.assertEqual(manifest["members"]["rows"], len(table))
        self.assertEqual(len(table), int(conditions["member_count"].sum()))
        self.assertEqual(manifest["treated"], manifest["rows"])
        # Per condition: exactly one chosen id, the discarded count is member_count - 1,
        # and the product row states both numbers next to the declared member table.
        self.assertIn("member_count", obs.columns)
        self.assertIn("discarded_members", obs.columns)
        discarded = []
        for condition in conditions.itertuples():
            chosen = [profile for profile in table["profile_id"].unique()
                      if profile in set(condition.profile_ids)]
            self.assertEqual(len(chosen), 1)
            self.assertIn(chosen[0], list(condition.profile_ids))
            listed = set(table.loc[table["profile_id"] == chosen[0], "member_profile_id"])
            self.assertEqual(listed, set(condition.profile_ids))
            row = obs.loc[obs["profile_id"] == chosen[0]].iloc[0]
            self.assertEqual(int(row["member_count"]), len(listed))
            self.assertEqual(int(row["discarded_members"]), len(listed) - 1)
            discarded.append(len(listed) - 1)
        self.assertEqual(sum(discarded), manifest["members"]["rows"] - manifest["rows"])
        self.assertEqual(sum(discarded), receipt["evidence"]["discarded_members"])
        # Retained totals are the chosen members' own full-view totals.
        retained = manifest["retained_totals"]
        totals = np.asarray([accepted_obs.loc[profile, "library_total"]
                             for profile in obs["profile_id"]], dtype=np.int64)
        self.assertEqual(retained["sum"], int(totals.sum()))
        self.assertEqual(retained["min"], int(totals.min()))
        self.assertEqual(retained["max"], int(totals.max()))
        # The control provenance names the chosen member's own vehicle group only.
        controls = pd.read_csv(self.root / "views/solo/tenv.control_sources.csv",
                               dtype={"profile_id": str})
        self.assertEqual(list(controls["profile_id"]), list(obs["profile_id"]))
        self.assertTrue(controls["control_rule"].eq(
            views.CONTROL_RULE["solo"]).all())
        groups = {profile: group for profile, group in
                  zip(accepted_obs["profile_id"], accepted_obs["matched_control_group"])}
        for row in controls.itertuples():
            self.assertEqual(row.control_groups, groups[row.profile_id])
        views.modules()["release"].check_relative_paths(manifest)


class SoloFreezeTest(unittest.TestCase):
    """Freeze seam: both member-choice products become verified Step 04-shaped inputs."""

    SIZE = 5000

    @classmethod
    def setUpClass(cls):
        views.bind(REPO)

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.root = Path(self._tmp.name)
        self.fixture = fixtures.dense_fixture(self.root)
        self.accepted, _, _ = fixtures.build_view(self.fixture, size=self.SIZE)
        self.population = fixtures.treated_population(self.accepted, "tenv")
        self.sealed = fixtures.write_sealed(self.root, "tenv", self.population)
        self.views_root = self.root / "views"

    def build(self, view):
        return build.build_product(
            view=view, source="tenv", output=self.views_root,
            staging=self.root / "staging", repo_root_path=REPO, size=self.SIZE,
            fixture=True, corpus=self.fixture["corpus"], axis=self.fixture["axis"],
            accepted=self.accepted, provenance_root=None)

    def freeze(self, view):
        return freeze.freeze_product(
            view=view, source="tenv", views_root=self.views_root,
            destination=self.root / "freezes" / view / "tenv", repo_root_path=REPO,
            sealed=self.sealed, axis_dir=self.fixture["axis"], accepted_dir=self.accepted)

    def accepted_rows(self, source="tenv"):
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

    def test_both_views_freeze_to_the_sealed_condition_fold_mapping(self):
        sealed_map = inputs.partition_map(pd.read_csv(self.sealed / "splits.csv"))
        accepted_obs, accepted_x, accepted_control = self.accepted_rows()
        position = {profile: row for row, profile in enumerate(accepted_obs["profile_id"])}
        for view in ("solo", "solo-deepest"):
            with self.subTest(view=view):
                self.build(view)
                receipt = self.freeze(view)
                self.assertTrue(receipt["verified"])
                self.assertFalse(receipt["splits_byte_equal_to_sealed"])
                data = inputs.Inputs(self.root / "freezes" / view / "tenv",
                                     sealed=self.sealed, repo_root=REPO)
                self.assertEqual(data.view, view)
                self.assertEqual(len(data.population), 9)
                self.assertEqual(len(data.splits), 9 * 6)
                self.assertEqual(data.config["delta_kind"], "expression_minus_control")
                self.assertTrue(data.controls_sources["control_rule"].eq(
                    "mean-of-normalized-vehicles").all())
                # The sealed condition -> partition mapping is preserved exactly.
                proven = inputs.partition_map(data.splits)
                joined = sealed_map.merge(
                    proven, on=["scenario", "fold", "drug_parent", "dose_um"],
                    suffixes=("_sealed", "_view"), indicator=True)
                self.assertTrue((joined["_merge"] == "both").all())
                self.assertTrue(
                    (joined["partition_sealed"] == joined["partition_view"]).all())
                self.assertEqual(set(zip(joined["drug_parent"], joined["dose_um"])),
                                 set(zip(sealed_map["drug_parent"],
                                         sealed_map["dose_um"])))
                # Each frozen row is the chosen full-view row verbatim.
                for row, profile in enumerate(data.population["profile_id"]):
                    source_row = position[str(profile)]
                    np.testing.assert_array_equal(
                        np.asarray(data.expression[row], dtype=np.float32),
                        accepted_x[source_row])
                    np.testing.assert_array_equal(
                        np.asarray(data.control[row], dtype=np.float32),
                        accepted_control[source_row])
                self.assertEqual(data.config["n_conditions"], 9)
                self.assertEqual(views.sha256(self.root / "freezes" / view / "tenv"
                                              / "parents.csv"),
                                 views.sha256(self.sealed / "parents.csv"))

    def test_solo_freeze_rejects_a_condition_that_moved_partition(self):
        self.build("solo")
        self.freeze("solo")
        directory = self.root / "freezes/solo/tenv"
        sealed = pd.read_csv(self.sealed / "splits.csv")
        row = sealed[sealed.partition == "test"].iloc[0]
        target = (row.scenario, row.fold, row.drug_parent, row.dose_um)

        copied = directory.parent / "tenv-tampered"
        shutil.copytree(directory, copied)
        splits = pd.read_csv(copied / "splits.csv")
        mask = ((splits.scenario == target[0]) & (splits.fold == target[1])
                & (splits.drug_parent == target[2]) & (splits.dose_um == target[3]))
        self.assertTrue(mask.any())
        splits.loc[mask, "partition"] = "train"
        splits.to_csv(copied / "splits.csv", index=False)
        manifest = json.loads((copied / "manifest.json").read_text())
        manifest["outputs"]["splits.csv"] = views.sha256(copied / "splits.csv")
        (copied / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        with self.assertRaisesRegex(views.HarnessError, "changed partition"):
            inputs.Inputs(copied, sealed=self.sealed, repo_root=REPO)


if __name__ == "__main__":
    unittest.main()
