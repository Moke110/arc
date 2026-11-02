"""Completeness and consistency gates on the written eligibility record.

v1.3 writes a versioned artifact set: ``output/v1.2/`` keeps the locked v1.2
record byte-identical, ``output/v1.3/`` is the revised record, and
``output/current.json`` points ticket #3 at the current handoff.
"""

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
V13 = OUT / "v1.3"
V12 = OUT / "v1.2"


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


class TestWrittenRecord(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table_path = V13 / "eligibility_table.csv"
        if not cls.table_path.exists():
            raise unittest.SkipTest(f"{cls.table_path} not built yet")
        cls.rows = read_csv(cls.table_path)
        cls.summary = json.loads((V13 / "summary.json").read_text(encoding="utf-8"))
        cls.handoff = json.loads((V13 / "handoff.json").read_text(encoding="utf-8"))
        cls.manifest = json.loads((V13 / "manifest.json").read_text(encoding="utf-8"))
        cls.current = json.loads((OUT / "current.json").read_text(encoding="utf-8"))
        cls.change = read_csv(V13 / "change_report.csv")

    def test_all_gates_pass(self):
        ok, lines = build_eligibility.check_artifacts(OUT)
        self.assertTrue(ok, "\n".join(lines))

    def test_460_rows_one_per_structure(self):
        self.assertEqual(len(self.rows), 460)
        self.assertEqual(len({row["drug_key"] for row in self.rows}), 460)
        frozen = evidence.frozen_keys(evidence.ROOT / build_eligibility.COMPOUNDS_CSV)
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
            "decision_precedence",
            "clinical_phase",
            "moa",
            "targets",
            "atc_codes",
            "mesh_frozen",
            "evidence_urls",
            "retrieved_at",
            "rationale",
            "human_use_basis",
            "other_use_basis",
            "uses_texts",
            "v12_decision",
            "v12_decision_class",
            "old_tier_atc_moa",
            "old_tier_atc_moa_mesh",
            "old_tier_all_smiles",
        ]
        for row in self.rows:
            for column in required:
                self.assertIn(column, row)
            for column in ("drug_key", "all_names", "inchikey", "decision",
                           "decision_class", "decision_precedence", "rationale",
                           "retrieved_at"):
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

    def test_rule_hash_is_the_locked_v13_one(self):
        self.assertEqual(self.summary["rule"]["sha256"], rules.rule_sha256())
        self.assertEqual(self.handoff["rule"]["sha256"], rules.rule_sha256())
        self.assertEqual(self.manifest["rule"]["sha256"], rules.rule_sha256())
        self.assertEqual(self.current["rule"]["sha256"], rules.rule_sha256())
        self.assertEqual(self.summary["rule"]["version"], rules.RULE_VERSION)
        self.assertEqual(self.summary["rule"]["version"], "eligibility-rule-v1.3")

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
            "unknown": "uncertain",
        }
        for row in self.rows:
            self.assertEqual(row["decision"], allowed[row["decision_class"]], row["drug_key"])

    def test_hash_manifest_matches_files(self):
        hashes = {item["path"]: item["sha256"] for item in self.manifest["files"]}
        for name in ("eligibility_table.csv", "handoff.json", "summary.json",
                     "change_report.csv", "old_tier_comparison.csv"):
            self.assertEqual(hashes[name], evidence.sha256_file(V13 / name), name)

    def test_handoff_is_the_v2_interface(self):
        self.assertEqual(self.handoff["schema"], "prnet-eligibility-handoff-v2")
        self.assertEqual(self.current["handoff"]["schema"], "prnet-eligibility-handoff-v2")
        self.assertEqual(self.current["handoff"]["path"], "v1.3/handoff.json")
        self.assertEqual(
            self.current["handoff"]["sha256"], evidence.sha256_file(V13 / "handoff.json")
        )

    def test_v12_record_is_preserved_byte_identical(self):
        v12_manifest = json.loads((V12 / "manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(v12_manifest["rule"]["version"], "eligibility-rule-v1.2")
        self.assertEqual(
            v12_manifest["counts"],
            {"total_keys": 460, "total_name_records": 467,
             "include": 126, "exclude": 193, "uncertain": 141},
        )
        for item in v12_manifest["files"]:
            self.assertEqual(
                item["sha256"], evidence.sha256_file(V12 / item["path"]), item["path"]
            )
        # v1.3 manifest carries the predecessor identity it was built from.
        self.assertEqual(self.manifest["predecessor"]["artifact_version"], "v1.2")
        self.assertEqual(
            self.manifest["predecessor"]["rule"]["sha256"],
            v12_manifest["rule"]["sha256"],
        )
        self.assertEqual(
            {item["path"]: item["sha256"] for item in self.manifest["predecessor"]["files"]},
            {item["path"]: item["sha256"] for item in v12_manifest["files"]},
        )


class TestPrecedenceAndChangeReport(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = V13 / "eligibility_table.csv"
        if not path.exists():
            raise unittest.SkipTest(f"{path} not built yet")
        cls.rows = read_csv(path)
        cls.change = read_csv(V13 / "change_report.csv")
        cls.summary = json.loads((V13 / "summary.json").read_text(encoding="utf-8"))
        cls.by_key = {row["drug_key"]: row for row in cls.rows}

    def test_precedence_label_matches_the_decision(self):
        expected = {
            "include": {"human_drug_use", "mechanism_tool"},
            "exclude": {"environmental_or_industrial_use"},
            "uncertain": {"identity_unresolved", "insufficient_evidence"},
        }
        for row in self.rows:
            self.assertIn(row["decision_precedence"], expected[row["decision"]],
                          row["drug_key"])

    def test_every_include_records_its_human_use_basis(self):
        for row in self.rows:
            if row["decision"] == "include":
                self.assertTrue(row["human_use_basis"].strip(), row["drug_key"])

    def test_every_exclude_records_its_other_use_basis(self):
        for row in self.rows:
            if row["decision"] == "exclude":
                self.assertTrue(row["other_use_basis"].strip(), row["drug_key"])

    def test_mechanism_tool_includes_carry_a_use_marker(self):
        # A mechanism-tool include must name the use/tool marker next to the
        # mechanism; a bare mechanism text is never enough.
        for row in self.rows:
            if row["decision_class"] == "mechanism_tool":
                self.assertIn("mechanism:", row["human_use_basis"], row["drug_key"])
                self.assertIn("use:", row["human_use_basis"], row["drug_key"])

    def test_change_report_lists_every_moved_row(self):
        moved = {
            row["drug_key"]: row
            for row in self.rows
            if (row["v12_decision"], row["v12_decision_class"])
            != (row["decision"], row["decision_class"])
        }
        reported = {row["drug_key"]: row for row in self.change}
        self.assertEqual(set(moved), set(reported))
        self.assertEqual(self.summary["change_report"]["moved"], len(moved))
        for key, row in moved.items():
            change = reported[key]
            self.assertEqual(change["v12_decision"], row["v12_decision"])
            self.assertEqual(change["v12_decision_class"], row["v12_decision_class"])
            self.assertEqual(change["v13_decision"], row["decision"])
            self.assertEqual(change["v13_decision_class"], row["decision_class"])
            self.assertTrue(change["deciding_evidence"].strip(), key)
            self.assertTrue(change["rationale_v13"].strip(), key)

    def test_change_report_counts_are_in_the_summary(self):
        blocks = self.summary["change_report"]
        self.assertEqual(blocks["moved"], len(self.change))
        self.assertEqual(
            blocks["by_movement"],
            {
                "include->exclude": sum(
                    1 for r in self.change
                    if r["v12_decision"] == "include" and r["v13_decision"] == "exclude"
                ),
                "exclude->include": sum(
                    1 for r in self.change
                    if r["v12_decision"] == "exclude" and r["v13_decision"] == "include"
                ),
                "uncertain->exclude": sum(
                    1 for r in self.change
                    if r["v12_decision"] == "uncertain" and r["v13_decision"] == "exclude"
                ),
                "uncertain->include": sum(
                    1 for r in self.change
                    if r["v12_decision"] == "uncertain" and r["v13_decision"] == "include"
                ),
                "include->uncertain": sum(
                    1 for r in self.change
                    if r["v12_decision"] == "include" and r["v13_decision"] == "uncertain"
                ),
                "exclude->uncertain": sum(
                    1 for r in self.change
                    if r["v12_decision"] == "exclude" and r["v13_decision"] == "uncertain"
                ),
            },
        )

    def test_rule_text_states_the_precedence(self):
        self.assertIn(
            "human drug use > documented environmental/industrial use > mechanism text",
            rules.RULE_TEXT,
        )
        self.assertIn("Mechanism text alone never satisfies eligibility", rules.RULE_TEXT)


if __name__ == "__main__":
    unittest.main()
