"""Per-(source, class) DEG table and pooled-class arithmetic tests (#189).

The fixtures are hand-worked: an eight-gene frozen axis, a two-stratum
multi-contrast class DEG, and a tiny DToxS Calc-table archive. The collapse and
rank rules themselves are covered by ``test_deg_rank.py``; here the tests assert
that the cardiotoxicity products consume that contract correctly and that the
pooled-class arithmetic is exactly "member wells + the same-plate CTRL wells".
"""
import io
import json
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import cardiotox, deg_rank, release


def axis_fixture():
    return pd.DataFrame({
        "gene_index": np.arange(8, dtype=np.int64),
        "symbol": [f"G{index}" for index in range(8)],
        "hgnc_id": [f"HGNC:{index}" for index in range(8)],
        "ensembl_id": [f"ENSG{index:011d}" for index in range(8)],
    })


def resolver_fixture():
    return deg_rank.SymbolResolver(axis_fixture())


def contrast_rows():
    return pd.DataFrame({
        "contrast": ["class|a", "class|a", "class|b", "class|b", "class|a"],
        "gene_symbol": ["G1", "G2", "G1", "G2", "NOPE"],
        "log2fc": [2.0, 0.5, -3.0, 1.0, 9.0],
        "padj": [1e-4, 0.5, 1e-3, np.nan, 1e-6],
        "stat": [5.0, 0.2, -8.0, 0.1, 99.0],
        "p_value": [1e-5, 0.5, 1e-8, 0.9, 1e-7],
        "base_mean": [10.0, 2.0, 11.0, 2.1, 5.0],
    })


class ClassDegTest(unittest.TestCase):
    def setUp(self):
        self.axis = axis_fixture()
        self.resolver = resolver_fixture()

    def test_multi_contrast_collapse_and_rank_contract(self):
        table, rank, stats = cardiotox.class_deg(contrast_rows(), self.resolver,
                                                 axis_size=8)
        self.assertEqual(stats["n_measured"], 2)
        self.assertEqual(stats["n_unmeasured"], 6)
        self.assertEqual(stats["rank_norm_denominator"], "n_measured")
        gene1 = table[table.gene_index.eq(1)].iloc[0]
        # |2| * -log10(1e-4) = 8 for class|a; |-3| * -log10(1e-3) = 9 for class|b.
        self.assertEqual(gene1.contrast, "class|b")
        np.testing.assert_allclose(gene1.effect, 9.0)
        self.assertEqual(int(gene1.direction), -1)
        self.assertEqual(int(gene1["rank"]), 1)
        np.testing.assert_allclose(gene1.rank_norm, 0.5)
        # The winning contrast's statistic/p-value/expression level are carried.
        np.testing.assert_allclose(gene1.stat, -8.0)
        np.testing.assert_allclose(gene1.p_value, 1e-8)
        np.testing.assert_allclose(gene1.base_mean, 11.0)
        gene2 = table[table.gene_index.eq(2)].iloc[0]
        self.assertEqual(gene2.contrast, "class|a")  # NA padj is non-significant
        self.assertFalse(bool(gene2.significant))
        self.assertEqual(int(gene2.direction), 0)

    def test_unmeasured_genes_get_n_measured_plus_one(self):
        table, rank, stats = cardiotox.class_deg(contrast_rows(), self.resolver,
                                                 axis_size=8)
        self.assertEqual(len(rank), 8)
        self.assertEqual(rank["rank"].tolist(), [3, 1, 2, 3, 3, 3, 3, 3])
        # n_measured = 2, unmeasured rank 3: rank_norm = 3/2 = 1.5 (ADR 0066).
        np.testing.assert_allclose(rank.rank_norm.tolist(),
                                   [1.5, 0.5, 1.0, 1.5, 1.5, 1.5, 1.5, 1.5])
        self.assertEqual(int(stats["unmeasured_rank"]), 3)

    def test_alignment_losses_are_recorded(self):
        table, _, stats = cardiotox.class_deg(contrast_rows(), self.resolver,
                                              axis_size=8)
        self.assertEqual(stats["alignment"]["n_source_symbols"], 5)
        self.assertEqual(stats["alignment"]["n_resolved"], 4)
        self.assertEqual(stats["alignment"]["unmapped"], 1)

    def test_deg_is_deterministic(self):
        first = cardiotox.class_deg(contrast_rows(), self.resolver, axis_size=8)
        second = cardiotox.class_deg(contrast_rows(), self.resolver, axis_size=8)
        pd.testing.assert_frame_equal(first[0], second[0])
        pd.testing.assert_frame_equal(first[1], second[1])
        self.assertEqual(first[2], second[2])


class PublishedTableTest(unittest.TestCase):
    def setUp(self):
        self.axis = axis_fixture()
        table, rank, _ = cardiotox.class_deg(contrast_rows(), resolver_fixture(),
                                             axis_size=8)
        self.table, self.rank = table, rank
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_deg_table_has_the_contract_schema(self):
        path = self.root / "deg.csv"
        cardiotox.write_deg_table(path, self.table, self.axis)
        frame = pd.read_csv(path)
        self.assertEqual(list(frame.columns), cardiotox.DEG_COLUMNS)
        self.assertEqual(len(frame), self.table.shape[0])
        self.assertEqual(frame.gene_symbol.tolist()[:2], ["G1", "G2"])
        self.assertTrue(frame.direction.dtype.kind in "iu")
        self.assertTrue(frame["rank"].dtype.kind in "iu")

    def test_rank_table_covers_the_full_axis(self):
        path = self.root / "rank.csv"
        cardiotox.write_rank_table(path, self.rank)
        frame = pd.read_csv(path)
        self.assertEqual(list(frame.columns), cardiotox.RANK_COLUMNS)
        self.assertEqual(frame.gene_index.tolist(), list(range(8)))
        self.assertEqual(len(frame), 8)


def _calc_header(controls, treated):
    normalised = [f"{well}.Norm" for well in controls + treated]
    return (controls + treated + normalised
            + ["logFC", "logCPM", "PValue", "FDR", "Regulation"])


def write_calc_archive(path, members):
    """members: {code: (controls, treated)} -> GSE217421 Calc table archive."""
    with tarfile.open(path, "w:gz") as archive:
        for code, (controls, treated) in members.items():
            header = _calc_header(list(controls), list(treated))
            row = ["GENE1"] + ["1"] * len(header)
            payload = ("\t".join(header) + "\n" + "\t".join(row) + "\n").encode()
            info = tarfile.TarInfo(
                f"Human-MSN01-03R-CM-Hour-48-Plate-0-Calc.CTRL-{code}.tsv")
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))


class PooledStrataTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_member_wells_are_pooled_against_one_shared_control_set(self):
        archive = self.root / "deg.tar.gz"
        controls = ["CTRL.A", "CTRL.B", "CTRL.C"]
        write_calc_archive(archive, {
            "DAU": (controls, ["DAU.A", "DAU.B"]),
            "DOX": (controls, ["DOX.A", "DOX.B", "DOX.C"]),
        })
        strata = cardiotox.dtox_pooled_strata(archive, {"DAU", "DOX"})
        self.assertEqual(len(strata), 1)
        stratum = next(iter(strata.values()))
        self.assertEqual(stratum["controls"], controls)
        self.assertEqual(stratum["treated"]["DAU"], ["DAU.A", "DAU.B"])
        self.assertEqual(sum(len(wells) for wells in stratum["treated"].values()), 5)

    def test_member_tables_must_share_the_control_wells(self):
        archive = self.root / "deg.tar.gz"
        write_calc_archive(archive, {
            "DAU": (["CTRL.A", "CTRL.B", "CTRL.C"], ["DAU.A"]),
            "DOX": (["CTRL.A", "CTRL.B", "CTRL.D"], ["DOX.A"]),
        })
        with self.assertRaises(cardiotox.CardiotoxError):
            cardiotox.dtox_pooled_strata(archive, {"DAU", "DOX"})

    def test_missing_members_are_an_error(self):
        archive = self.root / "deg.tar.gz"
        write_calc_archive(archive, {"DAU": (["CTRL.A"], ["DAU.A"])})
        with self.assertRaises(cardiotox.CardiotoxError):
            cardiotox.dtox_pooled_strata(archive, {"GEF"})


class CatalogShapeTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.axis = axis_fixture()
        table, rank, _ = cardiotox.class_deg(contrast_rows(), resolver_fixture(),
                                             axis_size=8)
        entry = cardiotox._dataset_entry(
            "dtoxs_ipsc", "atc_l4:L01EA", "BCR-ABL tyrosine kinase inhibitors",
            ["Bosutinib", "Dasatinib", "Imatinib"], ["BOS", "DAS", "IMA"],
            [{"contrast": "x|a", "stratum": "a", "n_case": 7, "n_control": 6}],
            {"n_measured": 2, "n_significant": 0, "alignment": {"n_resolved": 2},
             "informative": True, "informativeness_rule": "fixture",
             "n_distinct_effect": 2, "informativeness_reason": None},
            "x.deg.csv", "x.rank.csv", {}, {"name": "test", "implementation": "t"})
        entry["category"] = "cardiotox_other"
        self.record = {"deg": "x.deg.csv", "rank_name": "x.rank.csv",
                       "table": table, "rank": rank, "axis": self.axis,
                       "entry": entry}

    def tearDown(self):
        self.temporary.cleanup()

    def test_canonical_catalog_shape_and_digest(self):
        cardiotox.write_category_product(
            self.root / "cat", "cardiotox_other", [self.record], {}, {},
            "0" * 64)
        catalog = json.loads((self.root / "cat/catalog.json").read_text())
        self.assertEqual(catalog["schema"], cardiotox.CATALOG_SCHEMA)
        self.assertEqual(catalog["category"], "cardiotox_other")
        self.assertEqual(catalog["tier"], 2)
        self.assertEqual(catalog["evidence_type"], "bulk")
        self.assertIn("axis", catalog)
        self.assertIsInstance(catalog["datasets"], list)
        record = catalog["datasets"][0]
        self.assertEqual(record["dataset_id"], "dtoxs_ipsc")
        for dataset in catalog["datasets"]:
            counted = (int(dataset["n_case"]) > 0
                       and int(dataset["n_control"]) > 0)
            explained = bool(dataset.get("n_case_basis")
                             and dataset.get("n_control_basis"))
            self.assertTrue(counted or explained,
                            f"{dataset['dataset_id']} lacks n_case/n_control")
        self.assertIn("method", record)
        self.assertNotIn("test_method", record)
        self.assertIn("class_id", record)
        self.assertIn("n_measured", record)
        # The digest is verifiable from the catalog payload alone.
        verified = cardiotox.catalog_digest(self.root / "cat/catalog.json")
        self.assertEqual(verified, catalog["catalog_sha256"])
        manifest = json.loads((self.root / "cat/manifest.json").read_text())
        self.assertEqual(manifest["catalog_sha256"], verified)
        self.assertIn("catalog.json", manifest["outputs"])
        self.assertEqual(cardiotox.CategoryDeg(self.root / "cat",
                                               "cardiotox_other").catalog_digest,
                         verified)


class DatasetEntryTest(unittest.TestCase):
    def test_constant_arm_size_keeps_the_per_contrast_count(self):
        entry = cardiotox._dataset_entry(
            "gse244740", "atc_l4:L01DB", "Anthracyclines and related substances",
            ["Doxorubicin", "Idarubicin", "Mitoxantrone"],
            ["Doxorubicin", "Idarubicin", "Mitoxantrone"],
            [{"contrast": "x|plate", "stratum": "plate", "n_case": 72,
              "n_control": 22}],
            {"n_measured": 10, "n_significant": 1, "alignment": {},
             "informative": True, "informativeness_rule": "fixture",
             "n_distinct_effect": 2, "informativeness_reason": None},
            "x.deg.csv", "x.rank.csv", {}, {"name": "t"})
        self.assertEqual(entry["n_case"], 72)
        self.assertEqual(entry["n_control"], 22)
        self.assertIn("per-contrast", entry["n_case_basis"])


    def test_catalog_entry_records_the_required_fusion_fields(self):
        entry = cardiotox._dataset_entry(
            "dtoxs_ipsc", "atc_l4:L01EA", "BCR-ABL tyrosine kinase inhibitors",
            ["Bosutinib", "Dasatinib", "Imatinib"], ["BOS", "DAS", "IMA"],
            [{"contrast": "x|a", "stratum": "a", "n_case": 7, "n_control": 6},
             {"contrast": "x|b", "stratum": "b", "n_case": 9, "n_control": 8}],
            {"n_measured": 15000, "n_significant": 400,
             "alignment": {"n_resolved": 15000}, "informative": True,
             "informativeness_rule": "fixture", "n_distinct_effect": 2,
             "informativeness_reason": None}, "d.deg.csv", "d.rank.csv", {},
            {"test": "edgeR"})
        self.assertEqual(entry["dataset_id"], "dtoxs_ipsc")
        self.assertEqual(entry["class_id"], "atc_l4:L01EA")
        self.assertEqual(entry["class_slug"], "atc_l4-L01EA")
        self.assertEqual(entry["n_small_molecule_drugs"], 3)
        self.assertEqual(entry["n_contrasts"], 2)
        # Varying per-stratum arm sizes are summed over the disjoint strata, not null.
        self.assertEqual(entry["n_case"], 16)
        self.assertEqual(entry["n_control"], 14)
        self.assertTrue(entry["n_case_basis"])
        self.assertTrue(entry["n_control_basis"])
        self.assertEqual(entry["contrasts"][0]["n_case"], 7)
        self.assertEqual(entry["method"], {"test": "edgeR"})
        release.check_relative_paths(entry)


if __name__ == "__main__":
    unittest.main()
