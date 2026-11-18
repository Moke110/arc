"""Model-facing views: selected-axis normalization, control/delta, replay, index space.

Every expected number is hand-worked from the fixture's integer counts, never
recomputed through the implementation. Tests target the public build/reader
boundary named by ADR 0077 and issue #244.
"""
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import view

import view_fixtures


def hand_normalized(counts, totals):
    """log1p(CP10K) from the given per-row library totals, hand-worked."""
    counts = np.asarray(counts, dtype=np.float64)
    totals = np.asarray(totals, dtype=np.float64)
    return np.log1p(counts * (10000.0 / totals[:, None]))


class SelectedAxisNormalizationTest(unittest.TestCase):
    """A fixture where the full-axis and view-local denominators differ."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        fixture = view_fixtures.normalization_fixture(self.root)
        self.corpus = fixture["corpus"]
        self.axis = fixture["axis"]
        self.sources = fixture["sources"]
        self.release = self.root / "view"
        view.build(self.release, corpus_dir=self.corpus, axis_dir=self.axis,
                   size=2, fixture=True)
        self.opened = view.open_view(self.release)
        self.addCleanup(self.opened.close)

    def counted(self, name):
        return np.asarray(self.sources[name]["counts"], dtype=np.int64)

    def test_x_uses_selected_integer_totals(self):
        tenv = self.opened.source("tenv")
        # Fixture genes 0,3: [12,24], [6,6], [18,18], [3,3].
        expected = hand_normalized(self.counted("tenv")[:, [0, 3]],
                                   [36, 12, 36, 6])
        np.testing.assert_allclose(tenv.X(), expected, rtol=0, atol=1e-6)
        self.assertEqual(tenv.library_totals().tolist(), [36, 12, 36, 6])

    def test_control_is_the_mean_of_the_named_matched_controls(self):
        genes = [0, 3]
        tenv = self.opened.source("tenv")
        counts = self.counted("tenv")
        normalized = hand_normalized(counts[:, genes], counts[:, genes].sum(axis=1))
        expected = np.broadcast_to(normalized[[1, 2]].mean(axis=0), (3, 2))
        np.testing.assert_allclose(tenv.control()[[0, 1, 2]], expected,
                                   rtol=0.0, atol=1e-6)
        # A singleton vehicle group is its own control (the Media well).
        np.testing.assert_allclose(tenv.control()[3], normalized[3],
                                   rtol=0.0, atol=1e-6)

    def test_delta_is_x_minus_control_exactly(self):
        for name in self.opened.sources:
            source = self.opened.source(name)
            self.addCleanup(source.close)
            np.testing.assert_array_equal(source.delta(), source.X() - source.control())

    def test_obs_carries_the_corpus_row_identity_and_control_correspondence(self):
        tenv = self.opened.source("tenv")
        counts = self.counted("tenv")
        for column in ("role", "source_item", "plate", "cell_line", "compound_name",
                       "vehicle", "dose_um", "dose_original_value",
                       "drug_parent_inchikey", "drug_parent_smiles",
                       "matched_control_group", "library_total"):
            self.assertIn(column, tenv.obs.columns)
        self.assertEqual(tenv.obs.role.tolist(), self.sources["tenv"]["role"])
        self.assertEqual(tenv.obs.matched_control_group.tolist(),
                         self.sources["tenv"]["matched_control_group"])
        self.assertEqual(tenv.obs.library_total.tolist(), [36, 12, 36, 6])
        self.assertEqual(tenv.matched_control_rows(0).tolist(), [1, 2])

    def test_var_is_the_axis_ascending_and_never_renumbered(self):
        tenv = self.opened.source("tenv")
        self.assertEqual(tenv.gene_index.tolist(), [0, 3])
        self.assertEqual(tenv.var.gene_index.tolist(), [0, 3])
        self.assertEqual(tenv.shape, (4, 2))


class ControlScopeTest(unittest.TestCase):
    """A row's control is the mean of its own named group, never another plate's."""

    def test_control_is_scoped_to_the_rows_matched_group(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sources = {
                "tenv": {
                    "counts": [[4, 4, 0], [2, 2, 0], [8, 0, 8], [4, 0, 4]],
                    "role": ["treated", "vehicle", "treated", "vehicle"],
                    "matched_control_group": ["tenv::P1::DMSO", "tenv::P1::DMSO",
                                              "tenv::P2::DMSO", "tenv::P2::DMSO"],
                },
                "tahoe_bulk": {
                    "counts": [[1, 1, 1], [1, 0, 0]],
                    "role": ["treated", "vehicle"],
                    "matched_control_group": ["tahoe::pl1::L1", "tahoe::pl1::L1"],
                },
            }
            fixture = view_fixtures.write_release(root, sources, {2: [0, 2]})
            release = root / "view"
            view.build(release, corpus_dir=fixture["corpus"],
                       axis_dir=fixture["axis"], size=2, fixture=True)
            with view.open_view(release) as opened:
                tenv = opened.source("tenv")
                counts = np.asarray(sources["tenv"]["counts"], dtype=np.float64)
                totals = np.array([4, 2, 16, 8])
                normalized = hand_normalized(counts[:, [0, 2]], totals)
                np.testing.assert_allclose(tenv.control()[0], normalized[1],
                                           rtol=0.0, atol=1e-6)
                np.testing.assert_allclose(tenv.control()[2], normalized[3],
                                           rtol=0.0, atol=1e-6)
                self.assertFalse(np.allclose(tenv.control()[0], tenv.control()[2]))


class ViewReleaseTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        fixture = view_fixtures.normalization_fixture(self.root)
        self.corpus = fixture["corpus"]
        self.axis = fixture["axis"]
        self.release = self.root / "view"
        self.kwargs = dict(corpus_dir=self.corpus, axis_dir=self.axis, size=2,
                           fixture=True)
        view.build(self.release, **self.kwargs)

    def manifest(self):
        return json.loads((self.release / "manifest.json").read_text())

    def test_verify_accepts_the_built_release(self):
        self.assertEqual(view.verify(self.release, fixture=True)["result"], "PASS")

    def test_build_twice_produces_identical_digests(self):
        second = self.root / "second"
        view.build(second, **self.kwargs)
        view.verify(second, fixture=True)
        self.assertEqual(json.loads((second / "manifest.json").read_text()),
                         self.manifest())

    def test_replay_reproduces_the_same_delta(self):
        second = self.root / "replayed"
        view.build(second, **self.kwargs)
        with view.open_view(self.release) as accepted, view.open_view(second) as replayed:
            for name in accepted.sources:
                np.testing.assert_array_equal(accepted.source(name).X(),
                                              replayed.source(name).X())
                np.testing.assert_array_equal(accepted.source(name).control(),
                                              replayed.source(name).control())
                np.testing.assert_array_equal(accepted.source(name).delta(),
                                              replayed.source(name).delta())
        result = view.replay_release(self.release, root=self.root / "staging",
                                     **self.kwargs)
        self.assertEqual(result["result"], "PASS")

    def test_manifest_has_no_absolute_paths(self):
        text = (self.release / "manifest.json").read_text()
        self.assertNotIn(str(self.root), text)
        self.assertNotIn("/tmp/", text)
        from bases import release as release_module
        release_module.check_relative_paths(self.manifest())

    def test_no_view_truth_directory_is_created(self):
        for path in self.release.rglob("*"):
            self.assertNotIn("eval_truth", path.name)
        self.assertFalse((self.release / "2k_eval_truth").exists())
        self.assertFalse((self.release / "5k_eval_truth").exists())

    def test_tampered_payload_fails_verification(self):
        import h5py
        from bases import release as release_module
        with h5py.File(self.release / "tenv.h5ad", "a") as handle:
            handle["X"][0, 0] = np.float32(handle["X"][0, 0] + 1.0)
        with self.assertRaises(release_module.ReleaseError):
            view.verify(self.release, fixture=True)

    def test_delta_that_is_not_x_minus_control_fails_verification(self):
        import h5py
        from bases import release as release_module
        payload = self.release / "tenv.h5ad"
        with h5py.File(payload, "a") as handle:
            handle["layers/delta"][0, 0] = np.float32(
                handle["layers/delta"][0, 0] + 0.5)
        manifest = self.manifest()
        manifest["outputs"]["tenv.h5ad"]["sha256"] = release_module.sha256(payload)
        (self.release / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        with self.assertRaises(view.ViewError):
            view.verify(self.release, fixture=True)

    def test_verification_rejects_wrong_control_even_with_consistent_delta(self):
        import h5py
        from bases.release import sha256
        payload = self.release / "tenv.h5ad"
        with h5py.File(payload, "a") as handle:
            handle["layers/control"][0] = 0
            handle["layers/delta"][0] = handle["X"][0]
        manifest = self.manifest()
        manifest["outputs"]["tenv.h5ad"]["sha256"] = sha256(payload)
        (self.release / "manifest.json").write_text(json.dumps(manifest))
        with self.assertRaisesRegex(view.ViewError, "matched.*vehicle"):
            view.verify(self.release, fixture=True)

    def test_reader_windows_are_bounded_slices(self):
        with view.open_view(self.release) as opened:
            tenv = opened.source("tenv")
            np.testing.assert_array_equal(tenv.X(1, 3), tenv.X()[1:3])
            np.testing.assert_array_equal(tenv.control(1, 3), tenv.control()[1:3])
            np.testing.assert_array_equal(tenv.delta(1, 3), tenv.delta()[1:3])
            for start, x, control, delta in tenv.batches(2):
                np.testing.assert_array_equal(x - control, delta)


class ReleasedAxisSpaceTest(unittest.TestCase):
    """Synthetic axes preserve exact sizes and original indices; production is 5k only."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.tmp.name)
        fixture, axes = view_fixtures.full_size_fixture(cls.root)
        cls.corpus = fixture["corpus"]
        cls.axis = fixture["axis"]
        cls.axes = axes
        cls.releases = {}
        for size in (2000, 5000):
            output = cls.root / f"view-{size}"
            view.build(output, corpus_dir=cls.corpus, axis_dir=cls.axis, size=size,
                       fixture=True)
            cls.releases[size] = output

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_each_view_holds_exactly_its_axis_genes_ascending(self):
        for size, output in self.releases.items():
            with view.open_view(output) as opened:
                self.assertEqual(opened.size, size)
                self.assertEqual(opened.axis_size, 5010)
                for name in opened.sources:
                    source = opened.source(name)
                    self.assertEqual(source.shape[1], size)
                    expected = self.axes[size].tolist()
                    self.assertEqual(source.gene_index.tolist(), expected)
                    self.assertEqual(source.var.gene_index.tolist(), expected)
                    self.assertEqual(len(set(expected)), size)
                    self.assertEqual(expected, sorted(expected))
                    self.assertNotEqual(expected, list(range(size)))

    def test_released_views_verify(self):
        for output in self.releases.values():
            self.assertEqual(view.verify(output, fixture=True)["result"], "PASS")

    def test_manifest_records_selected_denominator_and_axis(self):
        for size, output in self.releases.items():
            manifest = json.loads((output / "manifest.json").read_text())
            for name in ("tenv", "tahoe_bulk"):
                summary = manifest["sources"][name]
                self.assertGreater(summary["selected_axis_library_total"]["sum"], 0)
                self.assertNotIn("full_axis_library_total", summary)
            self.assertIn("selected", manifest["normalization"]["denominator"])
            self.assertEqual(manifest["axis"]["size"], size)


class LocalScientificContractTest(unittest.TestCase):
    def build(self, root, sources):
        fixture = view_fixtures.write_release(root, sources, {2: [0, 2]})
        output = root / "view"
        view.build(output, corpus_dir=fixture["corpus"], axis_dir=fixture["axis"],
                   size=2, fixture=True)
        return output, fixture

    def sources(self):
        return {name: {
            "counts": [[2, 900, 6], [1, 100, 1], [9, 700, 3]],
            "role": ["treated", "vehicle", "vehicle"],
            "matched_control_group": [group] * 3,
        } for name, group in (("tenv", "tenv::P1::DMSO"),
                              ("tahoe_bulk", "tahoe::P1::L1"))}

    def test_each_unequal_library_vehicle_is_normalized_before_mean(self):
        with tempfile.TemporaryDirectory() as tmp:
            output, _ = self.build(Path(tmp), self.sources())
            # Vehicles have selected totals 2 and 12, ratios (1/2,1/2), (3/4,1/4).
            expected = np.array([(math.log1p(5000) + math.log1p(7500)) / 2,
                                 (math.log1p(5000) + math.log1p(2500)) / 2])
            with view.open_view(output) as opened:
                for name in opened.sources:
                    source = opened.source(name)
                    np.testing.assert_allclose(source.X(0, 1)[0],
                                               [math.log1p(2500), math.log1p(7500)], atol=1e-6)
                    np.testing.assert_allclose(source.control(0, 1)[0], expected, atol=1e-6)
                    self.assertEqual(source.library_totals().tolist(), [8, 2, 12])
                    self.assertEqual(source.matched_control_rows(0).tolist(), [1, 2])

    def test_unselected_counts_cannot_change_x_control_or_delta(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sources = self.sources()
            original, _ = self.build(root / "first", sources)
            for item in sources.values():
                for row in item["counts"]:
                    row[1] = 2000000000
            altered, _ = self.build(root / "altered", sources)
            with view.open_view(original) as first, view.open_view(altered) as other:
                for name in first.sources:
                    for field in ("X", "control", "delta", "library_totals"):
                        np.testing.assert_array_equal(getattr(first.source(name), field)(),
                                                      getattr(other.source(name), field)())

    def test_zero_selected_total_reports_profile_and_preserves_accepted(self):
        from bases.release import sha256
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output, _ = self.build(root / "original", self.sources())
            digest = sha256(output / "manifest.json")
            sources = self.sources()
            sources["tenv"]["counts"][0] = [0, 100, 0]
            inputs = view_fixtures.write_release(root / "invalid", sources, {2: [0, 2]})
            with self.assertRaisesRegex(view.ViewError, "tenv.*tenv:p0.*selected-axis"):
                view.publish_release(output, corpus_dir=inputs["corpus"],
                                     axis_dir=inputs["axis"], size=2, fixture=True,
                                     root=root / "staging")
            self.assertEqual(sha256(output / "manifest.json"), digest)
            self.assertEqual(view.verify(output, fixture=True)["result"], "PASS")

    def test_public_totals_identify_the_normalization_and_axis(self):
        with tempfile.TemporaryDirectory() as tmp:
            output, _ = self.build(Path(tmp), self.sources())
            with view.open_view(output) as opened:
                totals = opened.library_totals()
                self.assertEqual(set(totals.normalization_schema), {"stage2-view-v3"})
                self.assertEqual(set(totals.denominator_axis_sha256),
                                 {opened.manifest()["axis"]["sha256"]})

    def test_only_5k_is_a_production_view(self):
        with tempfile.TemporaryDirectory() as tmp:
            for size in (2000, 3000, 4000):
                with self.assertRaisesRegex(view.ViewError, "5000"):
                    view.build(Path(tmp) / str(size), corpus_dir="missing", axis_dir="missing", size=size)

    def test_missing_input_payload_cannot_replace_accepted(self):
        from bases.release import sha256
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output, inputs = self.build(root, self.sources())
            before = sha256(output / "manifest.json")
            (inputs["corpus"] / "tenv.h5ad").unlink()
            with self.assertRaises(FileNotFoundError):
                view.publish_release(output, corpus_dir=inputs["corpus"],
                                     axis_dir=inputs["axis"], size=2, fixture=True,
                                     root=root / "staging")
            self.assertEqual(sha256(output / "manifest.json"), before)

    def test_manifest_full_axis_generation_and_axis_mismatch_are_rejected(self):
        for field in ("schema", "axis"):
            with tempfile.TemporaryDirectory() as tmp:
                output, _ = self.build(Path(tmp), self.sources())
                manifest = json.loads((output / "manifest.json").read_text())
                if field == "schema":
                    manifest["schema"] = "stage2-view-v2"
                else:
                    manifest["axis"]["sha256"] = "0" * 64
                (output / "manifest.json").write_text(json.dumps(manifest))
                with self.assertRaises(view.ViewError):
                    view.open_view(output).manifest()

    def test_reader_rejects_stale_payload_digest(self):
        from bases.release import ReleaseError
        with tempfile.TemporaryDirectory() as tmp:
            output, _ = self.build(Path(tmp), self.sources())
            with (output / "tenv.h5ad").open("ab") as handle:
                handle.write(b"stale")
            with self.assertRaises(ReleaseError):
                view.open_view(output).source("tenv")

    def test_reader_rejects_old_and_mixed_metadata_and_axis(self):
        import h5py
        from anndata.io import read_elem, write_elem
        from bases.release import sha256
        for alteration in ("schema", "normalization", "axis"):
            with self.subTest(alteration=alteration), tempfile.TemporaryDirectory() as tmp:
                output, _ = self.build(Path(tmp), self.sources())
                payload = output / "tenv.h5ad"
                with h5py.File(payload, "a") as handle:
                    uns = read_elem(handle["uns"])
                    if alteration == "schema":
                        uns["schema"] = "stage2-view-v2"
                    elif alteration == "axis":
                        uns["axis"]["sha256"] = "0" * 64
                    else:
                        uns["normalization"]["denominator"] = "full-axis"
                    write_elem(handle, "uns", uns)
                manifest = json.loads((output / "manifest.json").read_text())
                manifest["outputs"]["tenv.h5ad"]["sha256"] = sha256(payload)
                (output / "manifest.json").write_text(json.dumps(manifest))
                with self.assertRaises(view.ViewError):
                    view.open_view(output).source("tenv")
                with self.assertRaises(view.ViewError):
                    view.verify(output, fixture=True)


if __name__ == "__main__":
    unittest.main()
