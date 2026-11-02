"""Tests for evidence extraction, identity mapping and the 467->460 collapse."""

import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import evidence  # noqa: E402

COMPOUNDS = evidence.ROOT / "stage2/02_cm_data/output/datasets-1/compounds.csv"


def pugview_record(sections, title="Example"):
    return json.dumps({"Record": {"RecordTitle": title, "Section": sections}})


def leaf(heading, strings):
    return {
        "TOCHeading": heading,
        "Information": [{"Value": {"StringWithMarkup": [{"String": s}]}} for s in strings],
    }


class TestNormalizeName(unittest.TestCase):
    def test_case_space_and_punctuation_insensitive(self):
        self.assertEqual(evidence.normalize_name("PK 11195"), "pk11195")
        self.assertEqual(evidence.normalize_name("3'-Azido-3'-deoxythymidine"),
                         "3azido3deoxythymidine")

    def test_salt_and_parent_are_not_merged(self):
        self.assertNotEqual(
            evidence.normalize_name("Diltiazem"),
            evidence.normalize_name("Diltiazem hydrochloride"),
        )


class TestExtractPubchem(unittest.TestCase):
    def test_extracts_documented_fields(self):
        record = pugview_record(
            [
                {
                    "TOCHeading": "Drug and Medication Information",
                    "Section": [
                        {
                            "TOCHeading": "Drug Development Summary",
                            "Section": [leaf("Max Phase", ["Phase 2"])],
                        },
                        leaf("Drug Classes", ["Anti-Inflammatory Agents"]),
                        leaf("Drug Indication", ["Pain and fever"]),
                    ],
                },
                {
                    "TOCHeading": "Pharmacology and Biochemistry",
                    "Section": [
                        leaf("ATC Code", ["N02BA01", "N - Nervous system"]),
                        leaf("Mechanism of Action", ["Blocks prostaglandin synthesis"]),
                    ],
                },
                {
                    "TOCHeading": "Use and Manufacturing",
                    "Section": [
                        {
                            "TOCHeading": "Uses",
                            "Section": [leaf("Use Classification", ["Pharmaceuticals"])],
                        }
                    ],
                },
            ]
        )
        out = evidence.extract_pubchem(record)
        self.assertEqual(out["pubchem_max_phase"], "Phase 2")
        self.assertEqual(out["atc_pubchem"], ["N02BA01", "N - Nervous system"])
        self.assertEqual(out["drug_classes"], ["Anti-Inflammatory Agents"])
        self.assertEqual(out["drug_indications"], ["Pain and fever"])
        self.assertEqual(out["mechanism_texts"], ["Blocks prostaglandin synthesis"])
        self.assertEqual(out["use_classification"], ["Pharmaceuticals"])

    def test_boilerplate_is_not_treated_as_use_evidence(self):
        record = pugview_record(
            [
                {
                    "TOCHeading": "Drug and Medication Information",
                    "Section": [
                        leaf(
                            "Therapeutic Uses",
                            ["ClinicalTrials.gov is a registry and results database of "
                             "publicly and privately supported clinical studies."],
                        )
                    ],
                }
            ]
        )
        out = evidence.extract_pubchem(record)
        self.assertEqual(out["therapeutic_uses"], [])

    def test_bare_uses_statements_are_extracted_separately(self):
        # The HSDB/EPA bare "Uses" statements document what the compound is used
        # for; v1.3 reads them as use evidence (v1.2 ignored them).
        record = pugview_record(
            [
                {
                    "TOCHeading": "Use and Manufacturing",
                    "Section": [
                        {
                            "TOCHeading": "Uses",
                            "Information": [
                                {
                                    "Value": {
                                        "StringWithMarkup": [
                                            {"String": "Used as a broad-spectrum "
                                                       "fungicide for cereals."}
                                        ]
                                    }
                                }
                            ],
                            "Section": [
                                leaf("Use Classification", ["Pesticides"]),
                                leaf("Industry Uses", ["Farming (Pesticides)"]),
                            ],
                        }
                    ],
                }
            ]
        )
        out = evidence.extract_pubchem(record)
        self.assertEqual(
            out["uses_texts"], ["Used as a broad-spectrum fungicide for cereals."]
        )
        self.assertEqual(out["use_classification"], ["Pesticides"])
        self.assertEqual(out["industry_uses"], ["Farming (Pesticides)"])

    def test_uses_texts_does_not_duplicate_subsection_leaves(self):
        record = pugview_record(
            [
                {
                    "TOCHeading": "Use and Manufacturing",
                    "Section": [
                        {
                            "TOCHeading": "Uses",
                            "Section": [leaf("Use Classification", ["Insecticides"])],
                        }
                    ],
                }
            ]
        )
        out = evidence.extract_pubchem(record)
        self.assertEqual(out["uses_texts"], [])
        self.assertEqual(out["use_classification"], ["Insecticides"])

    def test_empty_record_yields_no_claims(self):
        out = evidence.extract_pubchem(pugview_record([]))
        self.assertEqual(out["pubchem_max_phase"], "")
        self.assertEqual(out["mechanism_texts"], [])
        self.assertEqual(out["uses_texts"], [])


class TestLocalTimestamps(unittest.TestCase):
    def test_file_timestamp_is_the_file_mtime_not_now(self):
        # A local evidence source's retrieval timestamp must be a property of the
        # evidence file, so an offline rebuild yields a byte-identical record.
        from datetime import datetime, timezone

        path = evidence.BROAD_PHASE
        first = evidence.file_timestamp(path)
        second = evidence.file_timestamp(path)
        self.assertEqual(first, second)
        self.assertEqual(
            first,
            datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc)
            .astimezone()
            .isoformat(timespec="seconds"),
        )


class TestAtcCodes(unittest.TestCase):
    def test_pubchem_atc_token_is_split_from_its_name(self):
        codes = evidence.pubchem_atc_codes(
            ["QP53 - Ectoparasiticides, insecticides and repellents", "N02BA01"]
        )
        self.assertEqual(codes, ["QP53", "N02BA01"])

    def test_q_codes_are_veterinary_not_human_use(self):
        codes = evidence.pubchem_atc_codes(["QP53", "A01AD05"])
        human = [c for c in codes if not c.startswith("Q") and len(c) >= 5]
        self.assertEqual(human, ["A01AD05"])


class TestFrozenKeys(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.keys = evidence.frozen_keys(COMPOUNDS)
        cls.by_key = {row["drug_key"]: row for row in cls.keys}

    def test_460_structures_from_467_names(self):
        self.assertEqual(len(self.keys), 460)
        self.assertEqual(sum(row["n_name_records"] for row in self.keys), 467)
        self.assertEqual(len(self.by_key), 460)

    def test_multi_name_structures_keep_every_name(self):
        aspirin = self.by_key["BSYNRYMUTXBXSQ-UHFFFAOYSA-N"]
        self.assertEqual(aspirin["n_name_records"], 2)
        self.assertIn("Aspirin", aspirin["all_names"])
        self.assertIn("Acetylsalicylic acid", aspirin["all_names"])

    def test_no_structure_is_merged_or_dropped(self):
        self.assertTrue(all(row["inchikey"] for row in self.keys))
        self.assertTrue(all(row["old_tier_all_smiles"] for row in self.keys))
        self.assertEqual(
            sum(1 for row in self.keys if row["old_tier_atc_moa"]), 113
        )
        self.assertEqual(
            sum(1 for row in self.keys if row["old_tier_atc_moa_mesh"]), 252
        )


class TestIdentityRecord(unittest.TestCase):
    # The registered salt keeps its own InChIKey; drug_key is the desalted parent.
    row = {
        "drug_key": "HSUGRBWQSSZJOP-RTWAWAEBSA-N",
        "inchikey": "HSUGRBWQSSZJOP-RTWAWAEBSA-X",
        "compound_name": "Diltiazem hydrochloride",
        "smiles": "Cl.COC(=O)[C@H](c1ccccc1)N1CCN(CC1)c1ccc(OC)cc1",
    }

    def test_unresolved_when_no_source_matches(self):
        out = evidence.identity_record(
            self.row, {"cids": [], "http_status": 0}, {"broad_hub_member": False}, {}
        )
        self.assertFalse(out["identity_resolved"])
        self.assertIn("no pubchem cid", out["identity_notes"])

    def test_resolved_by_broad_membership_alone(self):
        out = evidence.identity_record(
            self.row, {"cids": []}, {"broad_hub_member": True}, {}
        )
        self.assertTrue(out["identity_resolved"])

    def test_salt_parent_and_name_mismatch_are_recorded(self):
        out = evidence.identity_record(
            self.row,
            {"cids": [39186], "synonyms": ["Diltiazem"]},
            {"broad_hub_member": True},
            {},
        )
        self.assertIn("salt/parent", out["identity_notes"])
        self.assertFalse(out["name_in_pubchem_synonyms"])


if __name__ == "__main__":
    unittest.main()
