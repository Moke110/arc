"""Truth seam: frozen primary/condition-mean truth and the noise diagnostics (plan C3).

The fixtures are literal-count corpora, so the pooling recipe, the per-member delta and the
diagnostic quantities can be recomputed in the test from those literals rather than from
the implementation.
"""
import json
import subprocess
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
import truth  # noqa: E402
import tsg_view_fixtures as fixtures  # noqa: E402
import views  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
#: The materialised full freeze lives only in the checkout that built the payload (a
#: worktree shares git objects, not ignored payloads), so it is named explicitly there.
FULL_FREEZE = fixtures.payload_root(REPO) / "stage2/tsg_view/execution/freezes/full/tenv"
SIZE = 3
GENES = [0, 3, 4]
CELL_LINE = {"tenv": "hiPSC-CM", "tahoe_bulk": "line1"}
PLATES = {"tenv": ("Plate1", "Plate2", "Plate3"), "tahoe_bulk": ("pl1", "pl2", "pl3")}
VEHICLE_COUNTS = [[5, 0, 0, 9, 0, 0], [7, 0, 0, 3, 0, 0]]


def _counts(index, plate):
    """Six-gene literal counts; genes 0 and 3 carry the varying signal."""
    return [1 + (index * 3 + plate) % 5, 40, 60, 1 + (index * 5 + 2 * plate) % 7, 50, 30]


def _source_spec(source):
    """Nine 3-dose x 2-member parents plus one singleton parent, one cell line."""
    identity, counts = [], []
    plates = PLATES[source]
    index = 0
    for parent_index in range(9):
        parent = f"P{parent_index + 1}"
        for dose_index, dose in enumerate((1.0, 10.0, 100.0)):
            for member in range(2):
                plate = plates[(parent_index + dose_index + member) % len(plates)]
                if parent == "P1" and dose == 1.0:
                    # a cross-plate repeat: the two members sit on different plates
                    plate = plates[member]
                identity.append({
                    "role": "treated", "plate": plate, "cell_line": CELL_LINE[source],
                    "cell_line_id": "", "dose": dose, "parent": parent,
                    "smiles": "CCO", "group": f"{source}::{plate}::DMSO"})
                counts.append(_counts(index, (parent_index + dose_index) % 5))
                index += 1
    for dose in (1.0, 10.0, 100.0):
        identity.append({"role": "treated", "plate": plates[0],
                         "cell_line": CELL_LINE[source], "cell_line_id": "",
                         "dose": dose, "parent": "SINGLE", "smiles": "CCO",
                         "group": f"{source}::{plates[0]}::DMSO"})
        counts.append(_counts(index, 1))
        index += 1
    for plate in plates:
        for well in VEHICLE_COUNTS:
            identity.append({"role": "vehicle", "plate": plate,
                             "cell_line": CELL_LINE[source], "cell_line_id": "",
                             "group": f"{source}::{plate}::DMSO", "vehicle": "DMSO"})
            counts.append(list(well))
    return {"counts": counts, "identity": identity}


def _fixture(root):
    spec = {"sources": {source: _source_spec(source) for source in views.SOURCES},
            "axes": {SIZE: GENES}}
    written = fixtures.write_release(root, spec["sources"], spec["axes"])
    accepted, manifest, _ = fixtures.build_view(written, size=SIZE)
    populations = {source: fixtures.treated_population(accepted, source)
                   for source in views.SOURCES}
    sealed = {source: fixtures.write_sealed(root, source, populations[source])
              for source in views.SOURCES}
    return {"written": written, "accepted": accepted, "populations": populations,
            "sealed": sealed, **spec}


class TruthTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        views.bind(REPO)
        cls._tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls._tmp.name)
        cls.fixture = _fixture(cls.root)
        cls.truth_root = cls.root / "truth"
        cls.receipt = cls.build(cls.truth_root)

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    @classmethod
    def build(cls, output, source="tenv", views_root=None):
        fixture = cls.fixture
        return truth.build_truth(
            source=source, output_dir=output, repo_root_path=REPO, size=SIZE,
            fixture=True, corpus_dir=fixture["written"]["corpus"],
            axis_dir=fixture["written"]["axis"], accepted_dir=fixture["accepted"],
            sealed_dir=fixture["sealed"][source], views_root=views_root)

    # -- helpers ---------------------------------------------------------- #
    def output(self, source="tenv"):
        return self.truth_root / source

    def sealed_population(self, source="tenv"):
        return pd.read_csv(self.fixture["sealed"][source] / "population.csv",
                           dtype={"profile_id": str})

    def sealed_splits(self, source="tenv"):
        return pd.read_csv(self.fixture["sealed"][source] / "splits.csv")

    def population(self, kind="primary", fold=0, source="tenv"):
        return pd.read_csv(self.output(source) / kind / f"fold-{fold}" / "population.csv",
                           dtype={"profile_id": str, "members": str})

    def payload(self, name, fold=0, kind="primary", source="tenv"):
        return np.load(self.output(source) / kind / f"fold-{fold}" / name)

    def noise_floor(self, source="tenv"):
        return pd.read_csv(self.output(source) / "noise_floor.csv",
                           dtype={"condition_key": str})

    def accepted_delta(self, source="tenv"):
        view = views.modules()["view"]
        with view.open_view(self.fixture["accepted"]) as release:
            with release.source(source) as reader:
                obs = reader.obs.copy()
                delta = np.concatenate(
                    [reader.delta(start, min(start + 256, len(obs)))
                     for start in range(0, len(obs), 256)], axis=0)
                control = np.concatenate(
                    [reader.control(start, min(start + 256, len(obs)))
                     for start in range(0, len(obs), 256)], axis=0)
        obs["profile_id"] = obs.index.astype(str)
        obs["view_row"] = np.arange(len(obs), dtype=np.int64)
        return obs, delta, control

    def corpus(self, source="tenv"):
        return views.modules()["corpus"].open_corpus(self.fixture["written"]["corpus"]
                                                     ).source(source)

    def normalise(self, counts):
        from bases import view as view_api

        counts = np.asarray(counts, dtype=np.int64).reshape(1, -1)
        return view_api.log1p_cp10k(counts, view_api.library_totals(counts))[0]

    # -- primary truth ---------------------------------------------------- #
    def test_primary_population_is_the_sealed_test_population(self):
        splits = self.sealed_splits()
        folds = sorted(set(splits[splits.scenario == truth.SCENARIO]["fold"]))
        for fold in folds:
            population = self.population("primary", fold)
            self.assertEqual(list(population.columns), list(truth.PRIMARY_COLUMNS))
            block = splits[(splits.scenario == truth.SCENARIO) & (splits.fold == fold)
                           & (splits.partition == "test")]
            self.assertEqual(population["profile_id"].tolist(),
                             block["profile_id"].astype(str).tolist())
            self.assertTrue(population["scenario"].eq(truth.SCENARIO).all())
            self.assertTrue(population["partition"].eq("test").all())
            self.assertTrue(population["fold"].eq(fold).all())
            self.assertEqual(self.payload("delta.npy", fold).shape,
                             (len(population), SIZE))
            self.assertEqual(self.payload("control.npy", fold).shape,
                             (len(population), SIZE))
            self.assertEqual(self.payload("delta.npy", fold).dtype, np.float32)
        total = sum(len(self.population("primary", fold)) for fold in folds)
        self.assertEqual(total, 19)  # nine two-member parents + one singleton parent

    def test_primary_delta_is_full_view_minus_control(self):
        obs, delta, control = self.accepted_delta()
        views_by_profile = dict(zip(obs["profile_id"], obs["view_row"]))
        for fold in range(3):
            population = self.population("primary", fold)
            rows = np.asarray([views_by_profile[profile]
                               for profile in population["profile_id"]], dtype=np.int64)
            published_delta = self.payload("delta.npy", fold)
            published_control = self.payload("control.npy", fold)
            self.assertTrue(np.array_equal(published_control, control[rows]))
            self.assertTrue(np.array_equal(published_delta, delta[rows]))
            self.assertTrue(np.array_equal(published_delta,
                                           (self.accepted_x(rows) - control[rows])))
        self.assertTrue(self.receipt["delta_equals_view_delta_layer"])

    def accepted_x(self, rows):
        view = views.modules()["view"]
        with view.open_view(self.fixture["accepted"]) as release:
            with release.source("tenv") as reader:
                return np.asarray(reader.X(int(rows.min()), int(rows.max()) + 1),
                                  dtype=np.float32)[rows - int(rows.min())]

    def test_primary_fingerprints_are_copied_byte_identically(self):
        sealed = self.fixture["sealed"]["tenv"]
        for fold in range(3):
            directory = self.output() / "primary" / f"fold-{fold}"
            self.assertEqual(views.sha256(directory / "parents.csv"),
                             views.sha256(sealed / "parents.csv"))
            self.assertEqual(views.sha256(directory / "morgan.npy"),
                             views.sha256(sealed / "morgan.npy"))

    def test_primary_manifest_records_the_fold_provenance(self):
        manifest = json.loads((self.output() / "primary" / "fold-0" / "manifest.json")
                              .read_text())
        self.assertEqual(manifest["schema"], truth.FOLD_SCHEMA)
        self.assertEqual(manifest["kind"], "primary")
        self.assertEqual(manifest["split_provenance"]["scenario"], truth.SCENARIO)
        for name in ("population.csv", "delta.npy", "control.npy", "parents.csv",
                     "morgan.npy"):
            self.assertEqual(manifest["outputs"][name]["sha256"],
                             views.sha256(self.output() / "primary" / "fold-0" / name))

    # -- condition-mean truth --------------------------------------------- #
    def test_condition_mean_only_replicated_test_conditions(self):
        splits = self.sealed_splits()
        for fold in range(3):
            population = self.population("condition_mean", fold)
            self.assertEqual(list(population.columns), list(truth.CONDITION_COLUMNS))
            self.assertTrue((population["count"] > 1).all())
            self.assertTrue(population["partition"].eq("test").all())
            members = [member for value in population["members"]
                       for member in str(value).split(";")]
            self.assertEqual(len(members), int(population["count"].sum()))
            self.assertEqual(len(set(members)), len(members))
            block = splits[(splits.scenario == truth.SCENARIO) & (splits.fold == fold)
                           & (splits.partition == "test")]
            expected = block.groupby(["drug_parent", "dose_um"]).size()
            expected = expected[expected > 1]
            self.assertEqual(len(population), len(expected))
            expected_profiles = {f"merge::tenv::{CELL_LINE['tenv']}::{parent}::{dose:g}"
                                 for parent, dose in expected.index}
            self.assertEqual(set(population["profile_id"]), expected_profiles)
            primary = set(self.population("primary", fold)["profile_id"])
            self.assertTrue(set(members) <= primary,
                            "a condition-mean member is absent from the fold's primary truth")
            self.assertEqual(self.payload("delta.npy", fold, kind="condition_mean").shape,
                             (len(population), SIZE))
            self.assertFalse((self.output() / "condition_mean" / f"fold-{fold}"
                              / "control.npy").exists())

    def test_condition_mean_partition_purity(self):
        splits = self.sealed_splits()
        for fold in range(3):
            population = self.population("condition_mean", fold)
            block = splits[(splits.scenario == truth.SCENARIO) & (splits.fold == fold)]
            partitions = dict(zip(block["profile_id"].astype(str), block["partition"]))
            for condition, members in zip(population["profile_id"], population["members"]):
                for member in str(members).split(";"):
                    self.assertEqual(partitions[member], "test",
                                     f"{condition} mixes a non-test member {member}")

    def hand_pooled_delta(self, condition):
        """Pool a condition from the fixture literals (never through the truth code)."""
        source = "tenv"
        reader = self.corpus(source)
        try:
            obs = reader.obs
            obs["profile_id"] = obs.index.astype(str)
            members = obs[(obs.role == "treated")
                          & (obs.drug_parent_inchikey == condition["drug_parent"])
                          & (obs.dose_um == float(condition["dose_um"]))].index.tolist()
            rows = [int(np.flatnonzero(obs.index.to_numpy() == member)[0])
                    for member in members]
            raw = np.asarray(reader.counts(0, len(obs)), dtype=np.int64)[:, GENES]
            groups = sorted({str(group) for group in
                             obs["matched_control_group"].to_numpy()[rows]})
            pooled = raw[rows].sum(axis=0)
            pooled_x = np.float32(self.normalise(pooled))
            group_rows = [int(position) for group in groups
                          for position in reader.control_rows(group)]
            vehicle = np.asarray(reader.counts(0, len(obs)), dtype=np.int64)
            pooled_control = np.float32(self.normalise(
                vehicle[group_rows][:, GENES].sum(axis=0)))
            member_controls = []
            for row in rows:
                local = [int(position) for position in
                         reader.control_rows(str(obs["matched_control_group"].iloc[row]))]
                member_controls.append(
                    np.asarray([self.normalise(vehicle[position][GENES])
                                for position in local], dtype=np.float64).mean(axis=0))
            member_x = np.asarray([self.normalise(raw[row]) for row in rows],
                                  dtype=np.float64)
            member_delta = (member_x.astype(np.float32)
                            - np.asarray(member_controls, dtype=np.float32))
        finally:
            reader.close()
        return (pooled_x - pooled_control).astype(np.float32), members, member_delta

    def test_condition_mean_delta_matches_the_pooling_recipe(self):
        population = pd.concat([self.population("condition_mean", fold) for fold in range(3)],
                               ignore_index=True)
        delta = np.concatenate([self.payload("delta.npy", fold, kind="condition_mean")
                                for fold in range(3)], axis=0)
        condition = population.iloc[0].to_dict()
        expected, _members, _member_delta = self.hand_pooled_delta(condition)
        position = int(np.flatnonzero(
            population["profile_id"].to_numpy() == condition["profile_id"])[0])
        self.assertTrue(np.array_equal(delta[position], expected))
        # the normalize-of-sum rule must differ from a mean of normalized members
        reader = self.corpus("tenv")
        try:
            obs = reader.obs
            obs["profile_id"] = obs.index.astype(str)
            rows = [int(np.flatnonzero(obs.index.to_numpy() == member)[0])
                    for member in str(condition["members"]).split(";")]
            raw = np.asarray(reader.counts(0, len(obs)), dtype=np.int64)[:, GENES]
            mean_of_normalised = np.asarray([self.normalise(raw[row]) for row in rows]
                                            ).mean(axis=0)
            pooled = np.float32(self.normalise(raw[rows].sum(axis=0)))
        finally:
            reader.close()
        self.assertGreater(float(np.max(np.abs(pooled - mean_of_normalised))), 1e-6)

    def test_condition_mean_delta_equals_merge_product(self):
        views_root = self.root / "views-merge"
        build.build_product(
            view="merge", source="tenv", output=views_root,
            staging=self.root / "staging-merge", repo_root_path=REPO, size=SIZE,
            fixture=True, corpus=self.fixture["written"]["corpus"],
            axis=self.fixture["written"]["axis"], accepted=self.fixture["accepted"],
            provenance_root=None)
        receipt = self.build(self.root / "truth-merge", views_root=views_root)
        self.assertTrue(receipt["merge_cross_check"]["available"])
        self.assertTrue(receipt["merge_cross_check"]["equal"])
        merge_manifest = json.loads((views_root / "merge" / "tenv.json").read_text())
        table = pd.read_csv(views_root / "merge" / merge_manifest["profile_ids"]["file"],
                            dtype={"profile_id": str})
        lookup = dict(zip(table["profile_id"], table["row"]))
        import h5py

        with h5py.File(views_root / "merge" / merge_manifest["payload"]["file"], "r") as handle:
            merged = np.asarray(handle["layers/delta"])
        for fold in range(3):
            population = self.population("condition_mean", fold, source="tenv")
            published = self.payload("delta.npy", fold, kind="condition_mean")
            for row, identity in enumerate(population["profile_id"]):
                self.assertTrue(np.array_equal(published[row], merged[lookup[identity]]))

    # -- diagnostics ------------------------------------------------------ #
    def test_noise_floor_and_repeat_correlation_match_hand_computed_deltas(self):
        population = pd.concat([self.population("condition_mean", fold) for fold in range(3)],
                               ignore_index=True)
        condition = population.iloc[0].to_dict()
        _delta, _members, member_delta = self.hand_pooled_delta(condition)
        frame = self.noise_floor()
        row = frame[frame["condition_key"] ==
                    f"{condition['cell_line']}|{condition['drug_parent']}|"
                    f"{float(condition['dose_um']):g}"]
        self.assertEqual(len(row), 1)
        deltas = np.asarray(member_delta, dtype=np.float64)
        expected_headline = float(np.sqrt(np.mean(np.var(deltas, axis=0, ddof=0))))
        expected_sample_sd = float(np.sqrt(np.mean(np.var(deltas, axis=0, ddof=1))))
        expected_rho = float(np.mean([
            np.corrcoef(deltas[left], deltas[right])[0, 1]
            for left in range(len(deltas)) for right in range(left + 1, len(deltas))]))
        self.assertAlmostEqual(float(row["noise_floor"].iloc[0]), expected_headline,
                               places=6)
        self.assertAlmostEqual(float(row["noise_floor_sample_sd"].iloc[0]),
                               expected_sample_sd, places=6)
        self.assertGreater(expected_sample_sd, expected_headline)
        self.assertAlmostEqual(float(row["rho_mean"].iloc[0]), expected_rho, places=6)
        self.assertEqual(int(row["member_count"].iloc[0]), len(deltas))
        self.assertEqual(int(row["weight_profile"].iloc[0]), len(deltas))

    def test_diagnostics_cover_every_replicated_condition_of_the_source(self):
        frame = self.noise_floor()
        conditions = self.fixture["populations"]["tenv"].copy()
        conditions["profile_id"] = conditions.index.astype(str)
        counts = conditions.groupby(["cell_line", "drug_parent_inchikey",
                                     "dose_um"]).size()
        replicated = counts[counts > 1]
        self.assertEqual(len(frame), len(replicated))
        self.assertTrue((frame["member_count"] > 1).all())
        self.assertEqual(int(frame["weight_profile"].sum()), int(replicated.sum()))

    def test_diagnostics_contract_and_reference_blocks(self):
        diagnostics = json.loads((self.output() / "diagnostics.json").read_text())
        self.assertEqual(diagnostics["schema"], truth.DIAGNOSTICS_SCHEMA)
        self.assertEqual(diagnostics["label"], "heuristic diagnostic")
        self.assertEqual(diagnostics["pcc_ceiling"]["label"], "heuristic diagnostic")
        self.assertEqual(diagnostics["pcc_ceiling"]["weighting"],
                         "member-profile count per replicated condition")
        self.assertEqual(diagnostics["pcc_ceiling"]["input_conditions"],
                         int(diagnostics["condition_population"][
                             "replicated_conditions"]))
        self.assertEqual(diagnostics["noise_floor"]["adr_reference_median"],
                         truth.ADR_REFERENCE["tenv"]["noise_floor_median"])
        self.assertIn("ddof=0", diagnostics["noise_floor_definition"])
        self.assertIn("ddof=1", diagnostics["noise_floor_sample_sd_definition"])
        self.assertIn("sensitivity only",
                      diagnostics["noise_floor_sample_sd_definition"])
        self.assertAlmostEqual(
            diagnostics["noise_floor"]["median"],
            float(self.noise_floor()["noise_floor"].median()), places=8)
        self.assertTrue(np.isfinite(diagnostics["pcc_ceiling"]["value"]))
        manifest = json.loads((self.output() / "manifest.json").read_text())
        self.assertEqual(manifest["schema"], truth.SCHEMA)
        self.assertTrue(manifest["view_independent"])
        self.assertEqual(manifest["diagnostics"]["label"], "heuristic diagnostic")
        for fold in range(3):
            self.assertEqual(manifest["counts"]["primary"][f"fold-{fold}"],
                             len(self.population("primary", fold)))
            outputs = manifest["primary"][f"fold-{fold}"]["outputs"]
            for name, record in outputs.items():
                self.assertEqual(record["sha256"],
                                 views.sha256(self.output() / "primary" / f"fold-{fold}"
                                              / name))

    def test_diagnostic_helpers_handle_the_identity_reductions(self):
        self.assertAlmostEqual(truth.noise_floor(np.zeros((2, SIZE))), 0.0)
        self.assertAlmostEqual(truth.noise_floor_sample_sd(np.ones((3, SIZE))), 0.0)
        deltas = np.vstack([np.arange(SIZE, dtype=np.float64),
                            np.arange(SIZE, dtype=np.float64) * 2.0])
        self.assertGreater(truth.noise_floor_sample_sd(deltas),
                           truth.noise_floor(deltas))
        self.assertTrue(np.isnan(truth.repeat_correlation(np.ones((3, SIZE)))))
        replicate = np.vstack([np.arange(SIZE) + 1.0] * 3)
        self.assertAlmostEqual(truth.repeat_correlation(replicate), 1.0)
        self.assertTrue(np.isnan(truth.repeat_correlation(np.zeros((2, SIZE)))))
        base = np.arange(SIZE, dtype=np.float64) + 1.0
        block = np.vstack([base, base * 2.0, base * 0.5])
        permuted = block[[2, 0, 1]]
        self.assertAlmostEqual(truth.noise_floor(block), truth.noise_floor(permuted),
                               places=12)
        self.assertAlmostEqual(truth.repeat_correlation(block),
                               truth.repeat_correlation(permuted), places=12)
        self.assertAlmostEqual(truth.pcc_ceiling([1.0], [7]), 1.0)
        self.assertAlmostEqual(truth.pcc_ceiling([0.25, 0.81], [1.0, 3.0]),
                               (0.5 + 2.7) / 4.0)
        with self.assertRaises(views.HarnessError):
            truth.condition_diagnostics(np.zeros((1, SIZE)))
        with self.assertRaises(views.HarnessError):
            truth.noise_floor(np.zeros(SIZE))

    def test_condition_mean_is_deterministic_byte_for_byte(self):
        first = self.truth_root / "tenv"
        second = self.root / "truth-again"
        self.build(second)
        names = sorted(path.relative_to(first) for path in first.rglob("*")
                       if path.is_file())
        again = sorted(path.relative_to(second / "tenv") for path in
                       (second / "tenv").rglob("*") if path.is_file())
        self.assertEqual(names, again)
        for name in names:
            self.assertEqual(views.sha256(first / name),
                             views.sha256(second / "tenv" / name), str(name))

    def test_truth_module_has_no_model_or_fit_dependency(self):
        source = (EXECUTION / "truth.py").read_text()
        for forbidden in ("torch", "learner", "matrix", "np.random", "random."):
            self.assertNotIn(forbidden, source, forbidden)
        # Check the import in a fresh interpreter: the suite's discovery may already have
        # imported the fitting stack through other test modules.
        command = [sys.executable, "-c",
                   ("import sys; sys.path.insert(0, %r); import truth; "
                    "assert not ({'torch', 'tsg_round1.learner'} & set(sys.modules)), "
                    "'truth import pulled in the fitting stack'" % str(EXECUTION))]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    @unittest.skipUnless(
        (FULL_FREEZE / "parents.csv").is_file(),
        fixtures.dependency_message("a materialised full freeze", FULL_FREEZE))
    def test_fixture_build_does_not_probe_the_materialised_full_freeze(self):
        """Fixture isolation: main-checkout geometry must not leak real freezes in.

        The regression scenario is the main checkout, where
        ``stage2/tsg_view/execution/freezes/full/<source>`` exists next to the real
        payloads. A fixture build must still use its synthetic sealed directory, never
        compare against the materialised full freeze.
        """
        main = fixtures.payload_root(REPO)
        fixture = self.fixture
        output = self.root / "truth-main-geometry"
        receipt = truth.build_truth(
            source="tenv", output_dir=output, repo_root_path=main, size=SIZE,
            fixture=True, corpus_dir=fixture["written"]["corpus"],
            axis_dir=fixture["written"]["axis"], accepted_dir=fixture["accepted"],
            sealed_dir=fixture["sealed"]["tenv"])
        manifest = json.loads((output / "tenv" / "manifest.json").read_text())
        self.assertFalse(receipt["merge_cross_check"]["available"],
                         "a fixture build must not cross-check materialised views")
        fingerprints = manifest["provenance"]["fingerprints"]
        self.assertTrue(fingerprints["byte_identical_to_sealed_freeze"])
        self.assertEqual(fingerprints["parents_sha256"],
                         views.sha256(fixture["sealed"]["tenv"] / "parents.csv"))
        self.assertEqual(receipt["counts"]["primary"]["fold-0"],
                         len(self.population("primary", 0)))

    def test_published_truth_reproduces_the_adr_reference_medians(self):
        """The committed real-data products must land on ADR 0082's reference medians."""
        published = HERE.parents[1] / "output/truth"
        if not (published / "tenv/diagnostics.json").is_file():
            self.skipTest("the materialised truth product is not in this checkout")
        for source in ("tenv", "tahoe_bulk"):
            diagnostics = json.loads((published / source / "diagnostics.json").read_text())
            reference = truth.ADR_REFERENCE[source]
            headline = diagnostics["noise_floor"]
            self.assertLessEqual(
                abs(headline["median"] - reference["noise_floor_median"]), 0.001,
                f"{source}: noise floor {headline['median']} misses the ADR reference "
                f"{reference['noise_floor_median']}")
            correlation = diagnostics["repeat_correlation"]
            self.assertLessEqual(
                abs(correlation["rho_median_pair"]
                    - reference["repeat_correlation_median"]), 0.001,
                f"{source}: repeat-correlation median {correlation['rho_median_pair']} "
                f"misses the ADR reference {reference['repeat_correlation_median']}")
            self.assertGreater(diagnostics["noise_floor_sample_sd"]["median"],
                               headline["median"])
            self.assertTrue(np.isfinite(diagnostics["pcc_ceiling"]["value"]))
            axis = {"tenv": (0.35, 0.42), "tahoe_bulk": (0.40, 0.55)}[source]
            self.assertTrue(axis[0] <= diagnostics["pcc_ceiling"]["value"] <= axis[1],
                            f"{source}: ceiling {diagnostics['pcc_ceiling']['value']} "
                            f"outside {axis}")

    def test_frozen_expectations_match_the_ticket(self):
        self.assertEqual(truth.EXPECTED_TEST_COUNTS["tenv"], (469, 468, 465))
        self.assertEqual(truth.EXPECTED_TEST_COUNTS["tahoe_bulk"],
                         (6643, 6650, 7950))
        self.assertEqual(truth.EXPECTED_CONDITIONS, {"tenv": 1377,
                                                     "tahoe_bulk": 52777})
        self.assertEqual(truth.EXPECTED_REPLICATED, {"tenv": 1377,
                                                     "tahoe_bulk": 9741})
        self.assertEqual(truth.ADR_REFERENCE["tenv"]["noise_floor_median"], 0.343)
        self.assertEqual(truth.ADR_REFERENCE["tahoe_bulk"]["noise_floor_median"], 0.070)


if __name__ == "__main__":
    unittest.main()
