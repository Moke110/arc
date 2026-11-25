"""Public contracts used by the Stage 2 consumer migration."""
import math
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from stage2_util import metrics


class MetricTests(unittest.TestCase):
    def test_correlations_and_errors(self):
        row = metrics.vector_metrics([1, 2, 3], [3, 2, 1])
        self.assertAlmostEqual(row["pcc"], -1.0)
        self.assertAlmostEqual(row["spearman"], -1.0)
        self.assertAlmostEqual(row["rmse"], math.sqrt(8 / 3))
        self.assertAlmostEqual(row["r2"], -3.0)

    def test_constant_masks_match_accepted_scoring(self):
        row = metrics.vector_metrics([1, 2, 3], [2, 2, 2])
        self.assertEqual(row["pcc"], 0.0)
        self.assertEqual(row["spearman"], 0.0)
        self.assertTrue(row["constant_prediction"])
        row = metrics.vector_metrics([2, 2, 2], [1, 2, 3])
        self.assertTrue(math.isnan(row["pcc"]))
        self.assertTrue(math.isnan(row["r2"]))

    def test_non_vector_metrics_reject_invalid_input(self):
        for truth, prediction in [([], []), ([[1, 2]], [[3, 4]]), ([1], [float("nan")])]:
            with self.subTest(truth=truth):
                with self.assertRaises(ValueError):
                    metrics.vector_metrics(truth, prediction)


class SeedTests(unittest.TestCase):
    def test_prefix_matches_accepted_seeds(self):
        from stage2_util import seeds

        self.assertEqual(seeds.read_seeds(), (17, 42, 2026))
        self.assertEqual(seeds.read_seeds(2), (17, 42))
        with self.assertRaises(ValueError):
            seeds.read_seeds(4)


class PlotTests(unittest.TestCase):
    def test_render_from_csv_and_palette(self):
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib import pyplot as plt
        from stage2_util import plotting

        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            values = directory / "points.csv"
            values.write_text("epoch,loss\n0,2.5\n1,1.5\n")
            fig, axis = plt.subplots()
            rows = plotting.read_rows(values)
            line, = plotting.line(axis, [int(row["epoch"]) for row in rows],
                                  [float(row["loss"]) for row in rows], hue=220, depth=0)
            self.assertEqual(list(line.get_ydata()), [2.5, 1.5])
            self.assertEqual(tuple(round(channel, 3) for channel in line.get_color()),
                             (.048, .176, .432, .90))
            png = directory / "plot.png"
            plotting.save_figure(fig, png)
            self.assertTrue(png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))
            plt.close(fig)


class IdentityTests(unittest.TestCase):
    def test_frozen_gene_axis_preserves_index(self):
        from stage2_util import identity

        genes = identity.load_gene_table()
        self.assertEqual(len(genes), 18164)
        self.assertEqual(genes[0]["hgnc_id"], "HGNC:5")
        self.assertEqual(genes[5091]["symbol"], "FAM153B")
        self.assertEqual(genes[5091]["gene_index"], 5091)


class RepresentationTests(unittest.TestCase):
    def test_ordered_parent_identities_collapse_repeated_profiles(self):
        import pandas as pd
        from stage2_util.representations import PARENT_COLUMNS, ordered_parents

        observations = pd.DataFrame(
            [("B", "CCC"), ("A", "CCO"), ("B", "CCC")], columns=PARENT_COLUMNS
        )
        self.assertEqual(ordered_parents(observations).values.tolist(),
                         [["A", "CCO"], ["B", "CCC"]])
        for invalid in ((None, "CCO"), ("", "CCO"), ("A", "CCC")):
            with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError, "parent"):
                ordered_parents(pd.DataFrame([("A", "CCO"), invalid], columns=PARENT_COLUMNS))

    def test_digest_bound_parent_order_and_vectors(self):
        import hashlib
        import json
        import numpy as np
        from stage2_util.representations import RepresentationTable, WIDTHS

        def digest(path):
            return hashlib.sha256(path.read_bytes()).hexdigest()

        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            parents = directory / "parents.csv"
            parents.write_text("drug_parent_inchikey,drug_parent_smiles\n"
                               "LFQSCWFLJHTTHZ-UHFFFAOYSA-N,CCO\n"
                               "OKKJLVBELUTLKV-UHFFFAOYSA-N,CO\n")
            conformers = directory / "e3fp_conformers.csv"
            conformers.write_text("drug_parent_inchikey,conformers\n"
                                  "LFQSCWFLJHTTHZ-UHFFFAOYSA-N,1\n"
                                  "OKKJLVBELUTLKV-UHFFFAOYSA-N,2\n")
            for arm, width in WIDTHS.items():
                values = np.zeros((2, width), dtype=np.float32)
                values[1, 1] = 1
                np.save(directory / f"{arm}.npy", values)
            manifest = directory / "manifest.json"
            manifest.write_text(json.dumps({
                "schema": "tsg-drug-emb-representations-v1",
                "source": "tenv",
                "inputs_sha256": "e" * 64,
                "widths": WIDTHS,
                "files": {path.name: digest(path) for path in directory.iterdir()},
            }))
            pin = digest(manifest)
            table = RepresentationTable(directory, source="tenv", manifest_sha256=pin)
            self.assertEqual(table.parents["drug_parent_smiles"].tolist(), ["CCO", "CO"])
            keys = ["OKKJLVBELUTLKV-UHFFFAOYSA-N", "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"]
            reordered = table.parents.iloc[[1, 0]]
            self.assertEqual(table.verified_parent_keys(reordered), keys)
            altered = reordered.assign(drug_parent_smiles=["CC", "CCO"])
            with self.assertRaisesRegex(ValueError, "structures differ"):
                table.verified_parent_keys(altered)
            self.assertEqual(table.vectors("ECFP4C", keys)[:, 1].tolist(), [1, 0])
            self.assertEqual(table.width("KPGT"), 2304)
            parents.write_text(parents.read_text().replace("CCO", "CC"))
            with self.assertRaisesRegex(ValueError, "parents.csv"):
                RepresentationTable(directory, source="tenv", manifest_sha256=pin)
