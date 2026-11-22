import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from build_corpus import open_corpus, seal_corpus
from cardiac_sources import build_count_table


class CardiacSourceTableTest(unittest.TestCase):
    def test_distinct_versioned_features_are_retained_without_claiming_ambiguous_coverage(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "counts.csv"
            source.write_text(
                "feature,well1,well2\n"
                "ENSG0001.1,11,0\n"
                "ENSG0001.2,2,13\n"
                "ENSG0002.4,0,1\n"
            )

            def metadata(wells):
                return pd.DataFrame({
                    "role": ["treated", "vehicle"],
                    "matched_control_group": ["plate", "plate"],
                }, index=wells)

            common_axis = pd.DataFrame({"ensembl_id": ["ENSG0001", "ENSG0002", "ENSG0003"]})
            coverage = build_count_table(root, "example", source, metadata, common_axis, "ensembl_id", ",")
            self.assertEqual(coverage, {1: 2})
            seal_corpus(root, {"example": coverage}, 3, "1" * 64, {}, {})
            with open_corpus(root) as corpus:
                data = corpus.source("example")
                self.assertEqual(data.var.source_gene_id.tolist(), [
                    "ENSG0001.1", "ENSG0001.2", "ENSG0002.4",
                ])
                self.assertEqual(data.counts(0, 2).tolist(), [[11, 2, 0], [0, 13, 1]])
                self.assertFalse(data.shared_axis_complete)
                self.assertEqual(data.matched_control_rows(0), (1,))

    def test_rejects_duplicate_publisher_well_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "duplicate.csv"
            source.write_text("feature,well1,well1\nENSG0001,1,2\n")
            with self.assertRaisesRegex(ValueError, "duplicate source observation"):
                build_count_table(root, "example", source, lambda _: None,
                                  pd.DataFrame(), "ensembl_id", ",")
