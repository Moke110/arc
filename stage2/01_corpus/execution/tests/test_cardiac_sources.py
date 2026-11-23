import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
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

    def test_rejects_duplicate_first_well_when_feature_header_is_omitted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "duplicate.tsv"
            source.write_text("well1\twell1\nENSG0001\t1\t2\n")
            with self.assertRaisesRegex(ValueError, "duplicate source observation"):
                build_count_table(root, "example", source, lambda _: None,
                                  pd.DataFrame(), "ensembl_id")

    def test_pooled_primary_counts_follow_gene_ids_and_keep_original_contexts(self):
        from cardiac_sources import merge_primary_tables

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            axis = pd.DataFrame({"symbol": ["g1", "g2"]})
            for name, rows, well_experiment in (
                ("first", "gene,well,control\ng1,3,8\ng2,7,4\n", "exp1"),
                ("second", "gene,well,existing\ng2,12,7\ng1,9,3\n", "exp2"),
            ):
                counts = root / f"{name}.csv"
                counts.write_text(rows)

                def metadata(wells):
                    return pd.DataFrame({
                        "Plate": "plate", "Cell": "cell", "Dish": "dish", "Time": "24",
                        "Well": ["well", "well" if name == "second" else "control"],
                        "Experiment": [well_experiment, "exp1"],
                        "State": ["drug", "drug" if name == "second" else "CTRL-DMSO"],
                        "role": ["treated", "treated" if name == "second" else "vehicle"],
                        "matched_control_group": "plate::cell::24::DMSO",
                    }, index=wells)

                build_count_table(root, name, counts, metadata, axis, "symbol", ",")
            tables = {"first": root / "first.h5ad", "second": root / "second.h5ad"}
            mapping = merge_primary_tables(root, tables, axis)
            seal_corpus(root, {"dtoxs_primary": mapping}, 2, "1" * 64, {}, {},
                        supplements={"dtoxs_primary": ["source_observations.csv"]})
            with open_corpus(root) as corpus:
                source = corpus.source("dtoxs_primary")
                self.assertEqual(source.counts(0, 3).tolist(), [[3, 7], [8, 4], [9, 12]])
                self.assertEqual(source.matched_control_rows(0), (1,))
                self.assertEqual(source.matched_control_rows(2), ())
                self.assertEqual(json.loads(source.obs.control_contexts.iloc[0]), [
                    ["first", "plate::cell::24::DMSO"], ["second", "plate::cell::24::DMSO"],
                ])
                ledger = pd.read_csv(root / "dtoxs_primary/source_observations.csv")
                self.assertEqual(ledger.source_profile_id.tolist(), ["well", "control", "well", "existing"])
                self.assertEqual(ledger.profile_id.nunique(), 3)

            with h5py.File(tables["second"], "r+") as handle:
                handle["X"][1, 0] = 44
            with self.assertRaisesRegex(ValueError, "conflicting counts"):
                merge_primary_tables(root, tables, axis)
