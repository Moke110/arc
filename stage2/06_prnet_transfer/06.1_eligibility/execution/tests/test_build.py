"""Completeness and consistency gates on the written eligibility record."""

import csv
import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import build_eligibility  # noqa: E402
import evidence  # noqa: E402
import rules  # noqa: E402

OUT = evidence.STEP / "output"


class TestWrittenRecord(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table_path = OUT / "eligibility_table.csv"
        if not cls.table_path.exists():
            raise unittest.SkipTest(f"{cls.table_path} not built yet")
        with open(cls.table_path, newline="", encoding="utf-8") as handle:
            cls.rows = list(csv.DictReader(handle))
        cls.summary = json.loads((OUT / "summary.json").read_text(encoding="utf-8"))
        cls.handoff = json.loads((OUT / "handoff.json").read_text(encoding="utf-8"))
        cls.manifest = json.loads((OUT / "manifest.json").read_text(encoding="utf-8"))

    def test_all_gates_pass(self):
        ok, lines = build_eligibility.check_artifacts(OUT)
        self.assertTrue(ok, "\n".join(lines))

    def test_460_rows_one_per_structure(self):
        self.assertEqual(len(self.rows), 460)
        self.assertEqual(len({row["drug_key"] for row in self.rows}), 460)
        frozen = evidence.frozen_keys(
            evidence.ROOT / build_eligibility.COMPOUNDS_CSV
        )
        self.assertEqual({row["drug_key"] for row in self.rows},
                         {row["drug_key"] for row in frozen})

    def test_required_columns_are_populated(self):
        required = [
            "drug_key",
            "all_names",
            "inchikey",
            "smiles",
            "decision",
            "decision_class",
            "clinical_phase",
            "moa",
            "targets",
            "atc_codes",
            "mesh_frozen",
            "evidence_urls",
            "retrieved_at",
            "rationale",
            "old_tier_atc_moa",
            "old_tier_atc_moa_mesh",
            "old_tier_all_smiles",
        ]
        for row in self.rows:
            for column in required:
                self.assertIn(column, row)
            for column in ("drug_key", "all_names", "inchikey", "decision",
                           "decision_class", "rationale", "retrieved_at"):
                self.assertTrue(row[column].strip(), f"{row['drug_key']} {column} empty")

    def test_every_include_and_exclude_is_traceable(self):
        for row in self.rows:
            if row["decision"] in {"include", "exclude"}:
                self.assertGreaterEqual(int(row["evidence_ok_urls"]), 1, row["drug_key"])
                self.assertTrue(row["evidence_urls"].strip(), row["drug_key"])
                self.assertTrue(row["rationale"].strip(), row["drug_key"])

    def test_uncertain_is_not_silently_included(self):
        uncertain = {row["drug_key"] for row in self.rows if row["decision"] == "uncertain"}
        self.assertEqual(uncertain, set(self.summary["uncertain_keys"]))
        self.assertEqual(len(uncertain), len(self.handoff["uncertain"]))
        include = {row["drug_key"] for row in self.rows if row["decision"] == "include"}
        self.assertFalse(uncertain & include)

    def test_rule_hash_is_the_locked_one(self):
        self.assertEqual(self.summary["rule"]["sha256"], rules.rule_sha256())
        self.assertEqual(self.handoff["rule"]["sha256"], rules.rule_sha256())
        self.assertEqual(self.manifest["rule"]["sha256"], rules.rule_sha256())
        self.assertEqual(self.summary["rule"]["version"], rules.RULE_VERSION)

    def test_counts_partition_the_table(self):
        counts = self.summary["counts"]
        self.assertEqual(counts["total_keys"], 460)
        self.assertEqual(counts["total_name_records"], 467)
        self.assertEqual(
            counts["include"] + counts["exclude"] + counts["uncertain"], 460
        )
        self.assertEqual(counts["include"], len(self.handoff["include"]))
        self.assertEqual(counts["exclude"], len(self.handoff["exclude"]))
        self.assertEqual(counts["uncertain"], len(self.handoff["uncertain"]))

    def test_class_and_decision_are_consistent(self):
        allowed = {
            "approved": "include",
            "clinical": "include",
            "mechanism_tool": "include",
            "other_use": "exclude",
            "unknown": {"uncertain"},
        }
        for row in self.rows:
            expected = allowed[row["decision_class"]]
            if isinstance(expected, str):
                self.assertEqual(row["decision"], expected, row["drug_key"])
            else:
                self.assertIn(row["decision"], expected, row["drug_key"])

    def test_hash_manifest_matches_files(self):
        hashes = {item["path"]: item["sha256"] for item in self.manifest["files"]}
        for name in ("eligibility_table.csv", "handoff.json", "summary.json"):
            self.assertEqual(hashes[name], evidence.sha256_file(OUT / name), name)


if __name__ == "__main__":
    unittest.main()
