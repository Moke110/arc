import json
import sys
import tempfile
import unittest
from pathlib import Path

STEP = Path(__file__).resolve().parents[2]
STAGE2 = STEP.parents[1]
sys.path.insert(0, str(STEP / "output/api"))
sys.path.insert(0, str(STAGE2 / "01_data_raw/output/api"))
from prediction_axis import open_axes


class PredictionAxisTest(unittest.TestCase):
    def test_canonical_axis_id_selects_the_published_gene_order(self):
        from prediction_axis import open_axis

        selected = open_axis("denoise-3k")
        self.assertEqual(selected.sha256, "fccdb1e8290cc93e7e7af97f3b2076626fc678f59e38c4ff56aacf48efa9a94b")
        self.assertEqual(selected.gene_indices, tuple(sorted(set(selected.gene_indices))))
        self.assertEqual(len(selected.gene_indices), 3000)
        with self.assertRaisesRegex(ValueError, "prediction axis"):
            open_axis("denoise-deg-3k")

    def test_denoise_deg_release_has_distinct_frozen_prediction_axes(self):
        release = open_axes(family="denoise")
        self.assertEqual(release.sizes, (1000, 2000, 3000, 4000, 5000))
        axis = release.axis(3000)
        self.assertEqual(axis.sha256, "fccdb1e8290cc93e7e7af97f3b2076626fc678f59e38c4ff56aacf48efa9a94b")
        self.assertEqual(axis.gene_indices, tuple(sorted(set(axis.gene_indices))))
        self.assertNotEqual(axis.sha256, open_axes().axis(3000).sha256)

    def test_release_has_original_ordered_axes(self):
        release = open_axes()
        self.assertEqual(release.sizes, (2000, 3000, 4000, 5000))
        evidence = json.loads((STEP / "output/release/evidence/gene-axis/accepted_manifest.json").read_text())
        for size in release.sizes:
            with self.subTest(size=size):
                axis = release.axis(size)
                original = f"stage2/01_bases/output/gene-axis/gene_axis_{size}.csv"
                self.assertEqual(axis.sha256, evidence["outputs"][original]["sha256"])
                self.assertEqual(len(axis.gene_indices), size)
                self.assertEqual(tuple(sorted(set(axis.gene_indices))), axis.gene_indices)
                self.assertLess(axis.gene_indices[-1], 18164)
        self.assertEqual(release.axis(2000).gene_indices[:3], (4, 35, 55))
        self.assertEqual(release.axis(2000).symbols[:3], ("A3GALT2", "ABCA4", "ABCC3"))

    def test_native_projection_preserves_integer_values_and_prediction_axis_order(self):
        import numpy as np
        import pandas as pd
        from anndata import AnnData
        from native_corpus import open_corpus, seal_corpus
        from prediction_axis import COMMON_SHA256, COMMON_SIZE

        axis = open_axes().axis(2000)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obs = pd.DataFrame({"role": ["treated", "vehicle"],
                                "matched_control_group": ["plate", "plate"]}, index=["drug", "vehicle"])
            var = pd.DataFrame(index=[f"native-{column}" for column in range(COMMON_SIZE)])
            matrix = np.zeros((2, COMMON_SIZE), dtype=np.int32)
            matrix[0, [18159, 18128, 18108]] = [149, 73, 211]
            mapping = {gene: COMMON_SIZE - 1 - gene for gene in range(COMMON_SIZE)}
            AnnData(matrix, obs=obs, var=var).write_h5ad(root / "reversed.h5ad")
            seal_corpus(root, {"reversed": mapping}, COMMON_SIZE, COMMON_SHA256, {}, {})

            with open_corpus(root) as corpus:
                projected = axis.project(corpus, "reversed")
                block = projected.counts(0, 1)
                self.assertEqual(block.shape, (1, 2000))
                self.assertEqual(block.dtype, np.int32)
                self.assertEqual(block[0, :3].tolist(), [149, 73, 211])
                self.assertFalse(block[0, 3:].any())
                self.assertEqual([start for start, _ in projected.batches(1)], [0, 1])

    def test_source_missing_a_common_gene_cannot_use_prediction_axes(self):
        import numpy as np
        import pandas as pd
        from anndata import AnnData
        from native_corpus import open_corpus, seal_corpus
        from prediction_axis import COMMON_SHA256, COMMON_SIZE

        axis = open_axes().axis(2000)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obs = pd.DataFrame({"role": ["vehicle"], "matched_control_group": ["plate"]},
                               index=["vehicle"])
            var = pd.DataFrame(index=[f"native-{gene}" for gene in range(COMMON_SIZE - 1)])
            matrix = np.zeros((1, COMMON_SIZE - 1), dtype=np.int32)
            # Index 0 is outside G_2000, so the source covers the selected genes.
            mapping = {gene: gene - 1 for gene in range(1, COMMON_SIZE)}
            AnnData(matrix, obs=obs, var=var).write_h5ad(root / "partial.h5ad")
            seal_corpus(root, {"partial": mapping}, COMMON_SIZE, COMMON_SHA256, {}, {})
            with open_corpus(root) as corpus, self.assertRaisesRegex(ValueError, "fixed common gene axis"):
                axis.project(corpus, "partial")

    def test_verify_recomputes_axis_selection_from_owner_evidence(self):
        sys.path.insert(0, str(STEP / "execution"))
        import build_axes

        result = build_axes.verify_release(STEP / "output/release")
        self.assertEqual(result.sizes, (2000, 3000, 4000, 5000))

    def test_publish_without_native_corpus_keeps_the_previous_release(self):
        from unittest.mock import patch

        sys.path.insert(0, str(STEP / "execution"))
        import build_axes

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "accepted"
            output.mkdir()
            marker = output / "manifest.json"
            marker.write_text("previous release")
            with patch.object(build_axes, "NATIVE_CORPUS", root / "missing"):
                with self.assertRaisesRegex(ValueError, "requires a native corpus manifest"):
                    build_axes.run("publish", root, output)
            self.assertEqual(marker.read_text(), "previous release")
