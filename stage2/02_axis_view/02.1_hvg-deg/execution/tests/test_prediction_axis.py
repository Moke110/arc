import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

STEP = Path(__file__).resolve().parents[2]
STAGE2 = STEP.parents[1]
sys.path.insert(0, str(STEP / "output/api"))
sys.path.insert(0, str(STAGE2 / "01_corpus/output/api"))
from prediction_axis import open_axes


class PredictionAxisTest(unittest.TestCase):
    def test_release_has_original_ordered_axes(self):
        release = open_axes()
        self.assertEqual(release.sizes, (2000, 3000, 4000, 5000))
        for size in release.sizes:
            with self.subTest(size=size):
                axis = release.axis(size)
                original = STAGE2 / "01_bases/output/gene-axis" / f"gene_axis_{size}.csv"
                self.assertEqual(axis.sha256, hashlib.sha256(original.read_bytes()).hexdigest())
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
