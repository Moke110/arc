"""TOX source measurement/admission/control boundary fixtures (ADR 0072)."""

import sys
import io
import gzip
import tarfile
import tempfile
import unittest
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api
from bases import tox_sources


class SourcePreparationTest(unittest.TestCase):
    def test_technical_repeat_is_summed_without_inventing_donor(self):
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp)
            root = repo / "disease_data/GSE313746/source"
            root.mkdir(parents=True)
            columns = [
                "DOX_tx0_1",
                "VEH_tx0_1",
                "DOX_tx0_6",
                "VEH_tx0_6",
                "DOX_tx0_6R",
                "VEH_tx0_6R",
            ]
            pd.DataFrame(
                [[101, 10, 20, 30, 40, 7, 9]], columns=["Entrez_ID"] + columns
            ).to_csv(
                root / "GSE313746_Counts_RNA_EMPfortmiller.txt.gz", sep=" ", index=False
            )
            ann = repo / "pert_data/Annotation/source"
            ann.mkdir(parents=True)
            pd.DataFrame(
                {
                    "entrez_id": ["101"],
                    "symbol": ["G1"],
                    "status": ["Approved"],
                    "locus_group": ["protein-coding gene"],
                    "hgnc_id": ["HGNC:1"],
                }
            ).to_csv(ann / "hgnc_complete_set.txt", sep="\t", index=False)
            drugs = pd.DataFrame(
                [
                    {
                        "source_id": "DOX",
                        "drug_name": "Doxorubicin",
                        "parent_inchikey": "KEY",
                    }
                ]
            )
            value = tox_sources.prepare_source(repo, "gse313746", drugs)
            part = value["strata"][0]
            self.assertEqual(len(part["samples"]), 4)
            self.assertEqual(sorted(set(part["samples"].donor)), ["donor_1", "donor_6"])
            self.assertEqual(part["matrix"]["DOX_tx0_6"].iloc[0], 37)
            self.assertEqual(part["matrix"]["VEH_tx0_6"].iloc[0], 49)
            self.assertEqual(part["design"]["formula"], "~0 + group + donor")
            self.assertEqual(len(value["audit"]), 6)
            repeat = next(
                row for row in value["audit"] if row["sample"] == "DOX_tx0_6R"
            )
            self.assertEqual(repeat["analysis_sample"], "DOX_tx0_6")
            self.assertEqual(repeat["matched_controls"], ["VEH_tx0_6", "VEH_tx0_6R"])

    def test_fpkm_keeps_unmodified_drug_and_excludes_coexposure(self):
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp)
            root = repo / "disease_data/GSE157282/source"
            root.mkdir(parents=True)
            with tarfile.open(root / "GSE157282_RAW.tar", "w") as archive:
                for name, value in [
                    ("Control1", 3),
                    ("Control2", 7),
                    ("Dox1", 15),
                    ("Dox2", 31),
                    ("DHGDox1", 255),
                ]:
                    data = gzip.compress(f"ENSG1\tG1\tchr1:1-2\t{value}\n".encode())
                    member = tarfile.TarInfo(f"GSM_{name}_S1.txt.gz")
                    member.size = len(data)
                    archive.addfile(member, io.BytesIO(data))
            drugs = pd.DataFrame(
                [
                    {
                        "source_id": "DOX",
                        "drug_name": "Doxorubicin",
                        "parent_inchikey": "KEY",
                    }
                ]
            )
            value = tox_sources.prepare_source(repo, "gse157282", drugs)
            part = value["strata"][0]
            self.assertEqual(part["design"]["kind"], "expression")
            self.assertEqual(part["matrix"]["GSM_Dox1_S1"].iloc[0], 4)
            self.assertEqual(len(part["samples"]), 4)
            excluded = next(row for row in value["audit"] if "DHGDox" in row["sample"])
            self.assertEqual(excluded["role"], "excluded")
            self.assertIn("co-exposure", excluded["reason"])

    def test_plate_vehicle_never_uses_blank_or_other_plate(self):
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp)
            root = repo / "pert_data/GSE262419/source"
            root.mkdir(parents=True)
            pd.DataFrame(
                {
                    "Plate_ID": [1, 1, 1, 2, 2],
                    "Well_ID": ["B2", "B3", "B4", "B2", "B3"],
                    "Chemical_name": ["VEH", "Known", "Blank", "VEH", "Known"],
                    "Chemical_Concentration_uM": [0, 1, 0, 0, 2],
                }
            ).to_csv(root / "GSE262419_hash.csv.gz", index=False)
            for plate in [1, 2]:
                values = {
                    "Genes": ["G1_001"],
                    f"Plate{plate}-B02": [10],
                    f"Plate{plate}-B03": [20],
                }
                if plate == 1:
                    values["Plate1-B04"] = [999]
                pd.DataFrame(values).to_csv(
                    root / f"GSE262419_Plate{plate:02}.csv.gz", index=False
                )
            value = tox_sources.prepare_source(
                repo,
                "gse262419",
                pd.DataFrame(
                    [
                        {
                            "source_id": "Known",
                            "drug_name": "Known",
                            "parent_inchikey": "KEY",
                        }
                    ]
                ),
            )
            self.assertEqual(len(value["strata"]), 2)
            treated = [r for r in value["audit"] if r["role"] == "treated"]
            self.assertEqual(treated[0]["matched_controls"], ["Plate1-B02"])
            self.assertEqual(treated[1]["matched_controls"], ["Plate2-B02"])
            self.assertNotIn("Plate1-B04", value["strata"][0]["matrix"])

    def test_sequencing_flowcell_is_adjusted_without_excluding_crossflow_dose(self):
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp)
            root = repo / "disease_data/GSE242692/source"
            root.mkdir(parents=True)
            names = [
                "c1_DOX_DMSO",
                "c1_DOX_1uM",
                "c1_DOX_03uM",
                "c2_DOX_DMSO",
                "c2_DOX_1uM",
            ]
            pd.DataFrame([[10] * 5], columns=names, index=["G1"]).to_csv(
                root / "GSE242692_countsMatrix_raw.txt.gz", sep="\t", index_label="gene"
            )
            fields = {
                "title": ["c1 0 µM", "c1 1 µM", "c1 0.3 µM", "c2 0 µM", "c2 1 µM"],
                "geo_accession": ["GSM1", "GSM2", "GSM3", "GSM4", "GSM5"],
                "description": ["c1", "c1", "c1", "c2", "c2"],
                "characteristics_ch1": [
                    "flowcell: Run1",
                    "flowcell: Run1",
                    "flowcell: Run2",
                    "flowcell: Run2",
                    "flowcell: Run2",
                ],
            }
            text = "".join(
                "!Sample_"
                + key
                + "\t"
                + "\t".join('"' + v + '"' for v in values)
                + "\n"
                for key, values in fields.items()
            )
            (root / "GSE242692_series_matrix.txt.gz").write_bytes(
                gzip.compress(text.encode())
            )
            value = tox_sources.prepare_source(
                repo,
                "gse242692",
                pd.DataFrame(
                    [
                        {
                            "source_id": "DOX",
                            "drug_name": "Doxorubicin",
                            "parent_inchikey": "KEY",
                        }
                    ]
                ),
            )
            self.assertEqual(len(value["strata"]), 1)
            self.assertEqual(len(value["strata"][0]["samples"]), 5)
            row = next(r for r in value["audit"] if r["sample"] == "c1_DOX_03uM")
            self.assertEqual(row["role"], "treated")
            self.assertEqual(row["matched_controls"], ["c1_DOX_DMSO"])
            self.assertEqual(row["flowcell"], "Run2")
            self.assertEqual(
                value["method"]["design_rank"], value["method"]["design_columns"]
            )

    def test_eht_uses_admitted_identity_without_a_doxorubicin_requirement(self):
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp)
            root = repo / "disease_data/GSE232331/source"
            root.mkdir(parents=True)
            names = ["VS_5uMACLA_080321_A1", "VS_DOXO_080321_A2", "VS_DMSO_080321_A3"]
            pd.DataFrame([[100, 200, 50]], columns=names, index=["G1"]).to_csv(
                root / "GSE232331_hg38-counts.tsv.gz", sep="\t", index_label="gene"
            )
            value = tox_sources.prepare_source(
                repo,
                "gse232331",
                pd.DataFrame(
                    [
                        {
                            "source_id": "ACLA",
                            "drug_name": "Aclarubicin",
                            "parent_inchikey": "KEY",
                        }
                    ]
                ),
            )
            treated = [r for r in value["audit"] if r["role"] == "treated"]
            self.assertEqual([r["source_id"] for r in treated], ["ACLA"])
            self.assertEqual(
                value["strata"][0]["samples"].batch.tolist(),
                ["date_080321", "date_080321"],
            )


if __name__ == "__main__":
    unittest.main()
