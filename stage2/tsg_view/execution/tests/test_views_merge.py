"""Build seam: the `merge` view pools replicates and controls by normalize-of-sum.

Worked examples come from the literal fixture counts in ``tsg_view_fixtures``: a condition's
expression is the log-normalized **sum** of its members' retained integer counts, never the
mean of its members' normalized values; its control is the log-normalized sum of the union
of its matched vehicle wells, never a mean of the plate-level controls.
"""
import json
import math
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
from views import merge  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
AXIS = fixtures.grid_sources()["axes"][2]


def read_product(path):
    """Read a materialized product's obs and layers straight from its declared layout."""
    import h5py
    from anndata.io import read_elem

    with h5py.File(path, "r") as handle:
        obs = read_elem(handle["obs"])
        return (obs.reset_index(drop=True),
                np.asarray(handle["X"], dtype=np.float32),
                np.asarray(handle["layers/control"], dtype=np.float32),
                np.asarray(handle["layers/delta"], dtype=np.float32))


def normalized(counts):
    counts = np.asarray(counts, dtype=np.float64)
    return np.log1p(counts * (10000.0 / counts.sum()))


class MergeViewBuildTest(unittest.TestCase):
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

    def build(self, source="tenv"):
        return build.build_product(
            view="merge", source=source, output=self.output,
            staging=self.root / f"staging-{source}", repo_root_path=REPO, size=2,
            fixture=True, corpus=self.fixture["corpus"], axis=self.fixture["axis"],
            accepted=self.accepted, provenance_root=None)

    def product(self, source="tenv"):
        return read_product(self.output / "merge" / f"{source}.h5ad")

    def manifest(self, source="tenv"):
        return json.loads((self.output / "merge" / f"{source}.json").read_text())

    def counts(self, source="tenv"):
        return self.fixture["sources"][source]["counts"]

    def accepted_obs(self, source="tenv"):
        with views.modules()["view"].open_view(self.accepted) as release:
            with release.source(source) as reader:
                obs = reader.obs.copy()
        obs["profile_id"] = obs.index.astype(str)
        return obs

    # -- one row per condition -------------------------------------------- #
    def test_one_materialized_row_per_condition_with_exact_delta(self):
        receipt = self.build("tenv")
        obs, x, control, delta = self.product("tenv")
        self.assertEqual(receipt["rows"], 9)
        self.assertEqual(receipt["conditions"], 9)
        self.assertEqual(len(obs), 9)
        self.assertTrue(obs["role"].eq("treated").all())
        self.assertEqual(obs["profile_id"].nunique(), 9)
        self.assertTrue(obs["profile_id"].str.startswith("merge::tenv::").all())
        self.assertEqual(x.shape, (9, len(AXIS)))
        self.assertFalse(obs["profile_id"].duplicated().any())
        np.testing.assert_array_equal(delta, x - control)

    # -- pooling identity -------------------------------------------------- #
    def test_pooled_expression_is_the_normalized_sum_not_the_member_mean(self):
        self.build("tenv")
        obs, x, control, _ = self.product("tenv")
        profile = "merge::tenv::hiPSC-CM::P1::1"
        row = int(np.flatnonzero(obs["profile_id"].to_numpy() == profile)[0])
        # Members are fixture rows 0 and 1: selected counts [4, 6] total 10 and [7, 11] 18.
        pooled = np.asarray([4 + 7, 6 + 11], dtype=np.float64)          # 11 and 17, total 28
        expected = np.log1p(pooled * (10000.0 / pooled.sum()))
        np.testing.assert_allclose(x[row], expected.astype(np.float32), rtol=0, atol=1e-6)
        member_mean = (np.log1p(np.asarray([4, 6]) * 1000.0)
                       + np.log1p(np.asarray([7, 11]) * (10000.0 / 18.0))) / 2
        self.assertFalse(np.allclose(x[row], member_mean, rtol=0, atol=1e-3))
        # The retained total is the sum of the members' selected-axis totals.
        self.assertEqual(int(obs["library_total"].iloc[row]), 28)

    def test_pooled_control_sums_the_vehicle_union_not_the_plate_mean(self):
        self.build("tenv")
        obs, _, control, _ = self.product("tenv")
        profile = "merge::tenv::hiPSC-CM::P1::100"
        row = int(np.flatnonzero(obs["profile_id"].to_numpy() == profile)[0])
        # Plate1 vehicles [6, 9] + [18, 3]; Plate2 vehicles [3, 12] + [9, 7] -> [36, 31] / 67.
        union = np.asarray([36, 31], dtype=np.float64)
        expected = np.log1p(union * (10000.0 / union.sum()))
        np.testing.assert_allclose(control[row], expected.astype(np.float32), rtol=0,
                                   atol=1e-6)
        plate1 = np.log1p(np.asarray([24, 12]) * (10000.0 / 36.0))
        plate2 = np.log1p(np.asarray([12, 19]) * (10000.0 / 31.0))
        self.assertFalse(np.allclose(control[row], (plate1 + plate2) / 2, rtol=0,
                                     atol=1e-3))
        # A singleton condition still deepens its control away from any one vehicle well.
        single = int(np.flatnonzero(
            obs["profile_id"].to_numpy() == "merge::tenv::hiPSC-CM::P3::1")[0])
        self.assertEqual(int(obs["library_total"].iloc[single]), 16)
        self.assertNotEqual(float(control[single][0]),
                            float(normalized(np.asarray([6, 9]))[0]))

    def test_cross_plate_repeat_pools_into_one_profile(self):
        self.build("tenv")
        obs, _, _, _ = self.product("tenv")
        repeats = obs[obs["member_count"] > 1] if "member_count" in obs else obs
        self.assertTrue("member_count" in obs)                # materialized obs records it
        self.assertEqual(int(obs["member_count"].max()), 2)
        self.assertEqual(obs["profile_id"].nunique(), len(obs))
        cross = obs[obs["matched_control_group"].str.count(";") == 1]
        self.assertEqual(len(cross), 1)
        self.assertEqual(cross["profile_id"].iloc[0], "merge::tenv::hiPSC-CM::P1::100")
        self.assertEqual(sorted(cross["plate"].iloc[0].split(";")), ["Plate1", "Plate2"])
        manifest = self.manifest("tenv")
        members = pd.read_csv(self.output / "merge" / manifest["members"]["file"])
        self.assertEqual(len(members[members.profile_id == "merge::tenv::hiPSC-CM::P1::100"]),
                         2)

    # -- manifest reconciliation ------------------------------------------- #
    def test_manifest_reconciles_members_totals_and_counts(self):
        self.build("tenv")
        manifest = self.manifest("tenv")
        obs, x, _, control = self.product("tenv")
        accepted = self.accepted_obs("tenv")
        treated = accepted[accepted.role == "treated"]
        grouped = list(treated.groupby(["cell_line", "drug_parent_inchikey", "dose_um"]))
        self.assertEqual(manifest["rows"], len(obs))
        self.assertEqual(manifest["conditions"], len(grouped))
        self.assertEqual(manifest["treated"], len(obs))
        self.assertEqual(manifest["view"], "merge")
        self.assertEqual(manifest["payload"]["layers"], ["control", "delta"])
        self.assertEqual(manifest["control_rule"], views.CONTROL_RULE["merge"])
        members = pd.read_csv(self.output / "merge" / manifest["members"]["file"])
        self.assertEqual(list(members.columns), ["profile_id", "member_profile_id"])
        self.assertEqual(len(members), len(treated))
        self.assertEqual(set(members.profile_id), set(obs["profile_id"]))
        self.assertEqual(set(members.member_profile_id), set(treated["profile_id"]))
        controls = pd.read_csv(self.output / "merge" / manifest["control_sources"]["file"])
        self.assertEqual(list(controls.profile_id), list(obs["profile_id"]))
        self.assertTrue(controls.control_rule.eq("pooled-sum").all())
        # Every pooled retained total reconciles to the sum of the members' totals ...
        expected_totals = treated.groupby(["cell_line", "drug_parent_inchikey", "dose_um"])[
            "library_total"].sum()
        expected_total = np.asarray([
            int(expected_totals.loc[(line, parent, float(dose))])
            for line, parent, dose in zip(obs["cell_line"], obs["drug_parent_inchikey"],
                                          obs["dose_um"])])
        np.testing.assert_array_equal(obs["library_total"].to_numpy(dtype=np.int64),
                                      expected_total)
        # ... and every control's retained total reconciles to the union of its vehicle wells.
        identity = self.fixture["sources"]["tenv"]["identity"]
        counts = self.counts("tenv")
        vehicle_wells = {}
        for offset, entry in enumerate(identity[len(treated):]):
            vehicle_wells.setdefault(entry["group"], []).append(len(treated) + offset)
        for position, profile in enumerate(obs["profile_id"]):
            member_positions = [int(pid.rsplit("p", 1)[1]) for pid in
                                members[members.profile_id == profile].member_profile_id]
            used = {identity[row]["group"] for row in member_positions}
            wells = [well for group in used for well in vehicle_wells[group]]
            expected = sum(int(np.asarray(counts[well], dtype=np.int64)[AXIS].sum())
                           for well in wells)
            self.assertEqual(int(controls.pooled_retained_total.iloc[position]), expected)
        # The manifest pins the companion tables it declares.
        for record in (manifest["members"], manifest["control_sources"],
                       manifest["profile_ids"]):
            path = self.output / "merge" / record["file"]
            self.assertEqual(views.sha256(path), record["sha256"])

    def test_manifest_names_the_tenv_traps_and_member_histogram(self):
        self.build("tenv")
        special = self.manifest("tenv")["special_conditions"]
        self.assertEqual(len(special["four_well_reference_compounds"]), 0)   # no traps in
        self.assertEqual(len(special["six_well_salt_pairs"]), 0)             # the fixture
        histogram = {int(key): int(value) for key, value in
                     special["member_count_histogram"].items()}
        self.assertEqual(histogram, {1: 2, 2: 7})
        self.assertEqual(special["cross_plate_conditions"], 1)
        self.assertEqual(special["source_has_named_traps"], False)
        # The pooled-identity probe travelled with the manifest.
        probe = self.manifest("tenv")["pooling_identity"]
        self.assertGreater(probe["checked_conditions"], 0)
        self.assertGreater(probe["differing_rows"], 0)
        self.assertGreater(probe["max_abs_difference"], 1e-6)

    def test_special_conditions_names_the_reference_compounds_and_salt_pairs(self):
        rows = []
        compounds = [f"Reference{i:02d}" for i in range(16)]
        for index, name in enumerate(compounds):
            for dose in (0.2, 1.0, 10.0):
                rows.append({"cell_line": "hiPSC-CM", "drug_parent_inchikey": f"PARENT{index}",
                             "dose_um": dose, "member_count": 4,
                             "drug_parent_smiles": "CCO",
                             "compound_names": [name], "plates": [f"Plate{index + 1}"],
                             "control_groups": [f"tenv::Plate{index + 1}::DMSO"],
                             "profile_ids": [f"tenv:p{index}-{dose}-{k}" for k in range(4)]})
        pairs = [("Tamoxifen", "Tamoxifen citrate"), ("Ammonium perfluorooctanoate",
                                                      "Perfluorooctanoic acid"),
                 ("Perfluorooctanesulfonic acid", "Potassium perfluorooctanesulfonate")]
        for index, names in enumerate(pairs):
            for dose in (0.2, 1.0, 10.0):
                rows.append({"cell_line": "hiPSC-CM", "drug_parent_inchikey": f"SALT{index}",
                             "dose_um": dose, "member_count": 6,
                             "drug_parent_smiles": "CCO",
                             "compound_names": list(names), "plates": ["Plate8"],
                             "control_groups": ["tenv::Plate8::DMSO"],
                             "profile_ids": [f"tenv:s{index}-{dose}-{k}" for k in range(6)]})
        rows.append({"cell_line": "hiPSC-CM", "drug_parent_inchikey": "PLAIN",
                     "dose_um": 1.0, "member_count": 3, "drug_parent_smiles": "CCO",
                     "compound_names": ["Plain"],
                     "plates": ["Plate1"], "control_groups": ["tenv::Plate1::DMSO"],
                     "profile_ids": ["a", "b", "c"]})
        special = merge.special_conditions(pd.DataFrame(rows), "tenv")
        four = special["four_well_reference_compounds"]
        six = special["six_well_salt_pairs"]
        self.assertEqual(len(four), 48)
        self.assertEqual(len(six), 9)
        self.assertEqual({entry["compound_names"][0] for entry in four}, set(compounds))
        self.assertEqual({name for entry in six for name in entry["compound_names"]},
                         {name for pair in pairs for name in pair})
        self.assertTrue(all(len(entry["compound_names"]) == 2 for entry in six))
        self.assertEqual({entry["drug_parent_inchikey"] for entry in six},
                         {"SALT0", "SALT1", "SALT2"})
        self.assertTrue(all(entry["member_count"] == 6 for entry in six))
        self.assertEqual(special["member_count_histogram"],
                         {"1": 0, "2": 0, "3": 1, "4": 48, "5": 0, "6": 9})
        # Only tenv carries the named traps; a Tahoe condition table names none.
        tahoe = merge.special_conditions(pd.DataFrame(rows), "tahoe_bulk")
        self.assertEqual(tahoe["four_well_reference_compounds"], [])
        self.assertEqual(tahoe["six_well_salt_pairs"], [])
        self.assertFalse(tahoe["source_has_named_traps"])

    # -- freeze seam ------------------------------------------------------- #
    def test_merge_freeze_keeps_the_sealed_condition_partitions(self):
        fixture = fixtures.dense_fixture(self.root / "dense")
        accepted, _, _ = fixtures.build_view(fixture, size=5000)
        output = self.root / "views-frozen"
        build.build_product(view="merge", source="tenv", output=output,
                            staging=self.root / "staging-frozen", repo_root_path=REPO,
                            size=5000, fixture=True, corpus=fixture["corpus"],
                            axis=fixture["axis"], accepted=accepted, provenance_root=None)
        population = fixtures.treated_population(accepted, "tenv")
        sealed = fixtures.write_sealed(self.root / "dense", "tenv", population)
        destination = self.root / "freezes/merge/tenv"
        receipt = freeze.freeze_product(view="merge", source="tenv", views_root=output,
                                        destination=destination, repo_root_path=REPO,
                                        sealed=sealed, axis_dir=fixture["axis"],
                                        accepted_dir=accepted)
        self.assertTrue(receipt["verified"])
        data = inputs.Inputs(destination, sealed=sealed, repo_root=REPO)
        self.assertEqual(len(data.population), 9)
        self.assertEqual(len(data.splits), 9 * 6)
        self.assertTrue(data.controls_sources["control_rule"].eq("pooled-sum").all())
        sealed_map = inputs.partition_map(pd.read_csv(sealed / "splits.csv"))
        view_map = inputs.partition_map(data.splits)
        joined = sealed_map.merge(view_map, on=["scenario", "fold", "drug_parent", "dose_um"],
                                  suffixes=("_sealed", "_view"))
        self.assertEqual(len(joined), len(sealed_map))
        self.assertTrue((joined.partition_sealed == joined.partition_view).all())
        # Every condition's members land in exactly one partition per scenario x fold.
        members = pd.read_csv(output / "merge" / "tenv.members.csv")
        full_population = fixtures.treated_population(accepted, "tenv")
        parent_of = dict(zip(full_population.profile_id,
                             full_population.drug_parent_inchikey))
        for profile, group in members.groupby("profile_id"):
            parents = {parent_of[member] for member in group.member_profile_id}
            self.assertEqual(len(parents), 1, profile)
        self.assertEqual(set(data.population["drug_parent_inchikey"]),
                         set(pd.read_csv(sealed / "parents.csv")["drug_parent_inchikey"]))
        # expression/control are the materialized pooled rows, delta exact float32.
        _, product_x, product_control, product_delta = read_product(
            output / "merge" / "tenv.h5ad")
        np.testing.assert_array_equal(np.asarray(data.expression), product_x)
        np.testing.assert_array_equal(np.asarray(data.control), product_control)
        np.testing.assert_array_equal(product_delta, product_x - product_control)
        np.testing.assert_array_equal(
            np.asarray(data.delta(range(len(data.population)))),
            np.asarray(data.expression, dtype=np.float64)
            - np.asarray(data.control, dtype=np.float64))


if __name__ == "__main__":
    unittest.main()
