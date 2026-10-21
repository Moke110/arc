"""Unit and integration tests for the Stage 2 CM corpus build.

Run from the repository root:

    python -m unittest discover -s stage2/02_cm_data/tests -v
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from cmdata import paths  # noqa: E402
from cmdata.annotation import CompoundAnnotation, HgncIndex  # noqa: E402
from cmdata.drugs import CompoundError, _validate_key, canonical_smiles  # noqa: E402
from cmdata.genes import (compute_universe, map_probe_space,  # noqa: E402
                          map_symbol_space, split_probe)
from cmdata.sources import _parse_conc_um  # noqa: E402
from cmdata.text import normalize_name  # noqa: E402
from cmdata.views import _replicate_index  # noqa: E402


def fake_hgnc() -> HgncIndex:
    table = pd.DataFrame([
        {"hgnc_id": "HGNC:1", "symbol": "A1BG", "alias_symbol": "",
         "prev_symbol": "", "entrez_id": "1", "ensembl_gene_id": "ENSG1"},
        {"hgnc_id": "HGNC:2", "symbol": "A2M", "alias_symbol": "A2M-ALIAS",
         "prev_symbol": "OLD2M", "entrez_id": "2", "ensembl_gene_id": "ENSG2"},
        {"hgnc_id": "HGNC:3", "symbol": "JPT2", "alias_symbol": "HN1L",
         "prev_symbol": "", "entrez_id": "90861", "ensembl_gene_id": "ENSG3"},
        {"hgnc_id": "HGNC:4", "symbol": "SHARED", "alias_symbol": "",
         "prev_symbol": "", "entrez_id": "904", "ensembl_gene_id": "ENSG4"},
        {"hgnc_id": "HGNC:5", "symbol": "OTHER", "alias_symbol": "DUP",
         "prev_symbol": "", "entrez_id": "905", "ensembl_gene_id": "ENSG5"},
        {"hgnc_id": "HGNC:6", "symbol": "THIRD", "alias_symbol": "DUP",
         "prev_symbol": "", "entrez_id": "906", "ensembl_gene_id": "ENSG6"},
    ])
    return HgncIndex(table)


class TextTests(unittest.TestCase):
    def test_normalize_name(self):
        self.assertEqual(normalize_name("Acetylsalicylic_acid"), "Acetylsalicylic acid")
        self.assertEqual(normalize_name("  Acetaminophen  "), "Acetaminophen")
        self.assertEqual(normalize_name("Mitomycin_C"), "Mitomycin C")
        self.assertEqual(normalize_name("5_Fluorouracil"), "5 Fluorouracil")
        self.assertEqual(normalize_name(None), "")

    def test_split_probe(self):
        self.assertEqual(split_probe("A1BG_25586"), ("A1BG", "25586"))
        self.assertEqual(split_probe("HLA_A_123"), ("HLA_A", "123"))
        self.assertEqual(split_probe("NOSUFFIX"), ("NOSUFFIX", ""))


class GeneTests(unittest.TestCase):
    def setUp(self):
        self.index = fake_hgnc()

    def test_symbol_current_alias_previous(self):
        self.assertEqual(self.index.resolve_symbol("A1BG"), ("HGNC:1", "direct"))
        self.assertEqual(self.index.resolve_symbol("A2M-ALIAS"), ("HGNC:2", "alias"))
        self.assertEqual(self.index.resolve_symbol("OLD2M"), ("HGNC:2", "previous"))
        self.assertEqual(self.index.resolve_symbol("JPT2"), ("HGNC:3", "direct"))
        self.assertEqual(self.index.resolve_symbol("HN1L"), ("HGNC:3", "alias"))

    def test_ambiguous_alias_is_dropped(self):
        self.assertEqual(self.index.resolve_symbol("DUP"), (None, "ambiguous"))

    def test_entrez_and_ensembl(self):
        self.assertEqual(self.index.resolve_entrez("90861")[0], "HGNC:3")
        self.assertEqual(self.index.resolve_ensembl("ENSG2.5")[0], "HGNC:2")

    def test_entrez_verifies_reassigned_symbol(self):
        mapping = map_symbol_space("dtoxs_primary", ["ALIAS_FOR_A2M"],
                                   self.index, {"ALIAS_FOR_A2M": "2"})
        self.assertEqual(mapping.hgnc_ids[0], "HGNC:2")
        self.assertEqual(mapping.match_types[0], "entrez")
        self.assertTrue(mapping.verified[0])

    def test_probe_space_collapses_by_symbol(self):
        mapping = map_probe_space("temp_env",
                                  ["A1BG_25586", "A1BG_1", "A2M_2", "NOPE_9"],
                                  self.index)
        self.assertEqual(list(mapping.hgnc_ids[:3]), ["HGNC:1", "HGNC:1", "HGNC:2"])
        self.assertIsNone(mapping.hgnc_ids[3])
        self.assertEqual(mapping.counts["unmapped"], 1)

    def test_intersection_and_order(self):
        m1 = map_symbol_space("a", ["A1BG", "A2M", "JPT2"], self.index)
        m2 = map_symbol_space("b", ["A1BG", "A2M", "HN1L"], self.index)
        m3 = map_probe_space("c", ["JPT2_1", "A2M_2", "A1BG_3"], self.index)
        universe = compute_universe(self.index, [m1, m2, m3])
        self.assertEqual(universe.symbols, ["A1BG", "A2M", "JPT2"])
        self.assertEqual(universe.var.loc["JPT2", "match_type"], "alias")
        self.assertEqual(universe.var.loc["A1BG", "match_type"], "direct")


class DrugIdentityTests(unittest.TestCase):
    def test_connectivity_only_key_rejected(self):
        with self.assertRaises(CompoundError):
            _validate_key("AOJJSUZBOXZQNB", "Doxorubicin")

    def test_full_key_accepted(self):
        key = "AOJJSUZBOXZQNB-TZSSRYMLSA-N"
        self.assertEqual(_validate_key(key, "Doxorubicin"), key)
        self.assertNotEqual(key, "AOJJSUZBOXZQNB-VTZDEGQISA-N")

    def test_canonical_smiles_retains_stereo(self):
        smile = "C[C@@H](O)c1ccccc1"
        canonical = canonical_smiles(smile)
        self.assertIn("@@", canonical)

    def test_structure_first_lookup(self):
        table = pd.DataFrame([
            {"lookup_name": "Amiodarone", "dataset_inchikey": "",
             "inchikey": "IYIKLHRQXLHMJQ-UHFFFAOYSA-N",
             "parent_inchikey": "IYIKLHRQXLHMJQ-UHFFFAOYSA-N",
             "smiles": "CCCCc1oc2ccccc2c1C(=O)c1cc(I)c(OCCN(CC)CC)c(I)c1",
             "pubchem_cid": "2157", "chembl_id": "CHEMBL633", "atc": "C01BD01",
             "moa": "", "mesh_classes": "Anti-Arrhythmia Agents",
             "verification": "verified", "notes": ""},
            {"lookup_name": "Amiodarone", "dataset_inchikey": "QMKSCTINHJSWPD-UHFFFAOYSA-N",
             "inchikey": "QMKSCTINHJSWPD-UHFFFAOYSA-N",
             "parent_inchikey": "QMKSCTINHJSWPD-UHFFFAOYSA-N",
             "smiles": "CCCCCc1oc2ccccc2c1C(=O)c1cc(I)c(OCCN(CC)CC)c(I)c1",
             "pubchem_cid": "16218975", "chembl_id": "", "atc": "",
             "moa": "", "mesh_classes": "",
             "verification": "verified", "notes": "homolog"},
        ])
        annotation = CompoundAnnotation(table)
        homolog = annotation.get("Amiodarone", "QMKSCTINHJSWPD-UHFFFAOYSA-N")
        name_only = annotation.get("Amiodarone", "")
        self.assertEqual(homolog["parent_inchikey"], "QMKSCTINHJSWPD-UHFFFAOYSA-N")
        self.assertEqual(name_only["parent_inchikey"], "IYIKLHRQXLHMJQ-UHFFFAOYSA-N")


class SourceParsingTests(unittest.TestCase):
    def test_concentration_units(self):
        self.assertAlmostEqual(_parse_conc_um("50nM"), 0.05)
        self.assertAlmostEqual(_parse_conc_um("2uM"), 2.0)
        self.assertAlmostEqual(_parse_conc_um("2 uM"), 2.0)
        self.assertAlmostEqual(_parse_conc_um("0.5uM /10uM"), 0.5)
        self.assertEqual(_parse_conc_um("0"), 0.0)

    def test_replicate_is_within_group_and_deterministic(self):
        meta = pd.DataFrame({
            "well": ["c", "a", "b", "d"],
            "control_group": ["g", "g", "g", "g"],
            "is_control": [False, False, False, False],
            "vehicle": ["DMSO"] * 4,
            "drug_key": ["K", "K", "K", "J"],
            "dose_um": [1.0, 1.0, 1.0, 2.0],
        })
        first = _replicate_index(meta)
        second = _replicate_index(meta.sample(frac=1, random_state=0))
        # Reordering the frame reorders the result, but each well keeps its index
        # within the dose/drug group.
        by_well = {well: int(value) for well, value in zip(meta["well"], first)}
        shuffled = meta.sample(frac=1, random_state=0).reset_index(drop=True)
        by_well_shuffled = {well: int(value) for well, value in zip(shuffled["well"], second)}
        self.assertEqual(by_well, by_well_shuffled)
        self.assertEqual(by_well["a"], 0)
        self.assertEqual(by_well["b"], 1)
        self.assertEqual(by_well["c"], 2)
        self.assertEqual(by_well["d"], 0)


class NamingTests(unittest.TestCase):
    def test_dataset_ids_map_to_accessions(self):
        self.assertEqual(set(paths.SOURCE_DATASETS),
                         {"dtoxs_primary", "dtoxs_ipsc", "screen_ctox",
                          "temp_env"})
        self.assertEqual(paths.SOURCE_ACCESSIONS,
                         {"screen_ctox": "GSE244740",
                          "temp_env": "GSE262419"})
        self.assertEqual(paths.source_accession("dtoxs_primary"), "")
        self.assertEqual(paths.stage1_dataset_name("screen_ctox"), "gse244740")
        self.assertEqual(paths.stage1_dataset_name("temp_env"), "gse262419")
        self.assertEqual(paths.stage1_dataset_name("dtoxs_primary"),
                         "dtoxs_primary")

    def test_raw_folders_keep_accession(self):
        self.assertEqual(paths.SCREEN_CTOX_DIR.name, "gse244740")
        self.assertEqual(paths.TEMP_ENV_DIR.name, "gse262419")

    def test_view_names_use_dataset_ids(self):
        for name, datasets, _ in paths.VIEWS:
            self.assertNotIn("gse", name)
            stem = name.removesuffix(".h5ad")
            if len(datasets) == 1:
                self.assertTrue(stem == datasets[0]
                                or stem.startswith(datasets[0] + "_"))
            else:
                self.assertIn("_cm", stem)

    def test_expected_views_match_views(self):
        self.assertEqual(set(paths.EXPECTED_VIEWS),
                         {name for name, _, _ in paths.VIEWS})


@unittest.skipUnless(paths.GENES_CSV.exists() and paths.BUNDLE_DIR.exists(),
                     "built corpus/bundle not present")
class BuiltBundleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import anndata as ad

        cls.genes = pd.read_csv(paths.GENES_CSV)
        cls.manifest = pd.read_csv(paths.BUNDLE_DIR / "manifest.csv")
        cls.views = {
            name: ad.read_h5ad(paths.BUNDLE_DIR / name, backed="r")
            for name in cls.manifest["view"]
        }

    def test_gene_axis(self):
        self.assertEqual(self.genes["symbol"].is_unique, True)
        self.assertEqual(len(self.genes), 18164)
        for adata in self.views.values():
            self.assertEqual(list(adata.var_names), self.genes["symbol"].tolist())

    def test_contract_columns(self):
        required = {"well", "dataset_id", "compound_name", "drug_key", "inchikey",
                    "smiles", "dose_um", "dose_source", "time_h", "is_control",
                    "vehicle", "plate", "control_group", "cell_type", "donor",
                    "replicate", "atc", "moa", "mesh_classes", "pharm_atc_moa",
                    "pharm_mesh", "annotation_source", "lib_size",
                    "n_genes_detected"}
        for adata in self.views.values():
            self.assertTrue(required.issubset(adata.obs.columns))
            self.assertEqual(list(adata.layers), ["log1p_cp10k"])
            self.assertEqual(adata.X.dtype, np.int32)
            for forbidden in ("pert_idx", "cell_idx", "dose_idx"):
                self.assertNotIn(forbidden, adata.obs.columns)

    def test_manifest_counts(self):
        for row in self.manifest.itertuples(index=False):
            adata = self.views[row.view]
            self.assertEqual(adata.n_obs, row.n_obs)
            self.assertEqual(adata.n_vars, 18164)
            controls = int(adata.obs["is_control"].sum())
            self.assertEqual(controls, row.n_controls)
            self.assertEqual(adata.n_obs - controls, row.n_perturbed)
            self.assertEqual(adata.obs.loc[~adata.obs["is_control"], "drug_key"].nunique(),
                             row.n_drugs)
            self.assertEqual(adata.uns["control_mean"].shape[0], row.n_control_groups)

    def test_control_rows_are_clean(self):
        for adata in self.views.values():
            controls = adata.obs[adata.obs["is_control"]]
            self.assertTrue((controls["dose_um"] == 0).all())
            self.assertTrue((controls["drug_key"] == "").all())
            self.assertTrue(controls["control_group"].isin(
                adata.obs["control_group"]).all())


if __name__ == "__main__":
    unittest.main()
