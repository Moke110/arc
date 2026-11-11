"""The release seam is exercised using independent, worked count examples."""

import sys
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
import anndata as ad
import h5py
from anndata.io import read_elem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bases.corpora import publish, verify_training_source


class CorpusReleaseTest(unittest.TestCase):
    def fixture(self, *, missing_control=False, invalid_dose=False, incompatible=False):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        native = root / "native"
        (native / "tables").mkdir(parents=True)
        (native / "tables/genes.csv").write_text(
            "symbol,hgnc_id,ensembl_id\nA,HGNC:1,E1\nB,HGNC:2,E2\n"
        )
        (native / "tables/gene_crosswalk.csv").write_text(
            "dataset_id,source_gene_id,hgnc_id\ntahoe,E1,HGNC:1\ntahoe,1,HGNC:1\ntahoe,E2,HGNC:2\ntahoe,2,HGNC:2\ntahoe,E3,\ntahoe,3,\n"
        )
        var = pd.DataFrame(
            {
                "hgnc_id": ["HGNC:1", "HGNC:2", "HGNC:3"],
                "ensembl_id": ["E1", "E2", "E3"],
                "token_id": [1, 2, 3],
            },
            index=["A", "B", "off"],
        )
        counts = np.array(
            [[1, 3, 900], [3, 1, 0], [2, 2, 100], [4, 0, 20], [9, 1, 0]], dtype="int32"
        )
        for source in ["tenv", "tahoe"]:
            obs = pd.DataFrame(
                {
                    "plate": ["P", "P", "P", "Q" if missing_control else "P", "P"],
                    "well": ["c1", "c2", "t1", "t2", "media"],
                    "is_control": [True, True, False, False, True],
                    "is_blank": [False] * 5,
                    "dose_um": [0, 0, 1, -1 if invalid_dose else 1, 0],
                    "compound_name": ["DMSO", "DMSO", "ethanol", "ethanol", "Media"],
                    "drug": ["DMSO_TF", "DMSO_TF", "ethanol", "ethanol", "OTHER"],
                    "smiles": ["", "", "CCO", "CCO", ""],
                    "inchikey": [
                        "",
                        "",
                        "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                        "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                        "",
                    ],
                    "drug_key": [
                        "",
                        "",
                        "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                        "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
                        "",
                    ],
                    "vehicle": ["DMSO"] * 5,
                    "time_h": [24, 24, 24, 48 if incompatible else 24, 24],
                    "cell_line_id": ["L"] * 5,
                    "cell_name": ["line"] * 5,
                },
                index=pd.Index(
                    ["c1", "c2", "t1", "t2", "media"],
                    name="well" if source == "tenv" else "profile_id",
                ),
            )
            ad.AnnData(counts, obs=obs, var=var).write_h5ad(native / f"{source}.h5ad")
        drugs = root / "drugs.parquet"
        pd.DataFrame({"drug": ["ethanol"], "canonical_smiles": ["CCO"]}).to_parquet(
            drugs
        )
        return native, root / "shared", drugs, root / "release"

    def test_separate_raw_locations_and_table_binding(self):
        from bases import raw
        native, shared, drugs, release = self.fixture()
        tahoe = native.parent / "tahoe-counts-v1"
        tahoe.mkdir()
        (native / "tahoe.h5ad").rename(tahoe / "tahoe.h5ad")
        tables = native.parent / "identity-tables"
        (native / "tables").rename(tables)
        bindings = {"tenv.h5ad": native / "tenv.h5ad",
                    "tahoe.h5ad": tahoe / "tahoe.h5ad",
                    "genes.csv": tables / "genes.csv",
                    "gene_crosswalk.csv": tables / "gene_crosswalk.csv"}
        raw.build(bindings, shared, fixture=True)
        publish(bindings, shared, drugs, release, fixture=True)
        self.assertEqual(publish(bindings, shared, drugs, release,
                                 fixture=True, verify=True)["result"], "PASS")
        for source in ("tenv", "tahoe"):
            check = verify_training_source(bindings, release, source, tables,
                                           drugs, fixture=True)
            self.assertEqual(check["result"], "PASS")

    def test_per_well_projection_control_mean_and_label_isolation(self):
        native, shared, drugs, release = self.fixture()
        result = publish(native, shared, drugs, release, fixture=True)
        self.assertEqual(result["sources"]["tenv"]["samples"], 2)
        self.assertEqual(
            ad.read_h5ad(release / "tenv.h5ad").shape, (2, 2)
        )  # tiny fixture
        with h5py.File(release / "tenv.h5ad") as f:
            obs = read_elem(f["obs"])
            self.assertEqual(list(obs.cell_identity), ["hiPSC-CM"] * 2)
            self.assertEqual(list(obs.source_profile_id), ["t1", "t2"])
            self.assertEqual(set(f["layers"]), {"control", "delta"})
            np.testing.assert_allclose(f["X"][0], [8.517393, 8.517393], rtol=1e-6)
            np.testing.assert_allclose(
                f["layers/control"][0], [8.373619, 8.373619], rtol=1e-6
            )
            self.assertEqual(obs.matched_control_ids.iloc[0], '["c1","c2"]')
        self.assertEqual(
            publish(native, shared, drugs, release, fixture=True, verify=True)[
                "result"
            ],
            "PASS",
        )

    def test_missing_control_excludes_well_with_explicit_accounting(self):
        native, shared, drugs, release = self.fixture(missing_control=True)
        result = publish(native, shared, drugs, release, fixture=True)
        self.assertEqual(
            result["sources"]["tenv"]["excluded_by_reason"],
            {"missing_matched_control": 1},
        )
        self.assertEqual(result["sources"]["tenv"]["samples"], 1)

    def test_invalid_dose_excludes_well_instead_of_grouping_it(self):
        native, shared, drugs, release = self.fixture(invalid_dose=True)
        result = publish(native, shared, drugs, release, fixture=True)
        self.assertEqual(
            result["sources"]["tenv"]["excluded_by_reason"], {"invalid_treated_dose": 1}
        )

    def test_parent_dose_conflicting_conditions_fail_before_publication(self):
        native, shared, drugs, release = self.fixture(incompatible=True)
        with self.assertRaisesRegex(
            ValueError, "incompatible or unresolved tenv parent-dose"
        ):
            publish(native, shared, drugs, release, fixture=True)
        self.assertFalse(release.exists())

    def test_salt_and_free_acid_keep_wells_and_auditable_parent_remapping(self):
        native, shared, drugs, release = self.fixture()
        data = ad.read_h5ad(native / "tenv.h5ad")  # tiny fixture only
        data.obs["smiles"] = data.obs["smiles"].astype(object)
        data.obs["drug_key"] = data.obs["drug_key"].astype(object)
        data.obs.loc["t1", "smiles"] = "CC(=O)O"
        data.obs.loc["t2", "smiles"] = "CC(=O)[O-].[Na+]"
        data.obs.loc["t1", "drug_key"] = "QTBSBXVTEAMEQO-UHFFFAOYSA-N"
        data.obs.loc["t2", "drug_key"] = "QTBSBXVTEAMEQO-UHFFFAOYSA-M"
        data.write_h5ad(native / "tenv.h5ad")
        result = publish(native, shared, drugs, release, fixture=True)
        self.assertEqual(result["sources"]["tenv"]["samples"], 2)
        self.assertEqual(result["sources"]["tenv"]["parent_dose_groups"], 1)
        with h5py.File(release / "tenv.h5ad") as file:
            obs = read_elem(file["obs"])
            self.assertEqual(obs.drug_parent.nunique(), 1)
            self.assertEqual(obs.source_parent_key.nunique(), 2)
            self.assertIn("[Na+]", obs.source_smiles.iloc[1])

    def test_source_replay_joins_stable_profiles_after_raw_row_reordering(self):
        native, shared, drugs, release = self.fixture()
        publish(native, shared, drugs, release, fixture=True)
        raw = ad.read_h5ad(native / "tenv.h5ad")  # tiny fixture only
        raw = raw[[3, 1, 4, 0, 2]].copy()
        raw.write_h5ad(native / "tenv.h5ad")
        result = verify_training_source(native, release, "tenv", native / "tables", fixture=True)
        self.assertEqual(result["result"], "PASS")
        self.assertEqual(result["rows_replayed"], 2)
