import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from anndata import AnnData

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "output/api"))
from native_corpus import CorpusError, open_corpus, seal_corpus


class NativeCorpusTest(unittest.TestCase):
    def test_preserves_counts_coverage_and_unmatched_observations(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            obs = pd.DataFrame({
                "role": ["treated", "vehicle", "treated"],
                "matched_control_group": ["plate-1", "plate-1", "plate-2"],
            }, index=["drug-1", "vehicle-1", "drug-2"])
            var = pd.DataFrame(index=["native-a", "native-b", "native-c"])
            AnnData(np.array([[0, 7, 2], [4, 0, 8], [1, 2, 0]], dtype=np.int32),
                    obs=obs, var=var).write_h5ad(root / "mini.h5ad")
            seal_corpus(root, {"mini": {0: 1, 2: 0}}, 3, "1" * 64, {}, {})
            with open_corpus(root) as corpus:
                data = corpus.source("mini")
                self.assertEqual(data.obs.index.tolist(), ["drug-1", "vehicle-1", "drug-2"])
                self.assertEqual(data.matched_control_rows(0), (1,))
                self.assertEqual(data.matched_control_rows(2), ())
                self.assertEqual(data.obs.match_status.tolist(), ["matched", "vehicle", "unmatched"])
                self.assertEqual(data.native_columns([2, 0]), (0, 1))
                self.assertEqual(data.counts(0, 1).tolist(), [[0, 7, 2]])
                with self.assertRaises(CorpusError):
                    data.native_columns([1])
