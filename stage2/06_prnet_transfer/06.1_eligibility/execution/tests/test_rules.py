"""Tests for the locked eligibility rule engine (rules.py).

Expectations come from the frozen annotation format (``compounds.csv`` ``atc``
field: ``CODE Name; CODE Name``) and from the decision rule text, not from the
implementation.
"""

import unittest

from rules import (
    RULE_VERSION,
    atc_level_counts,
    decide,
    parse_atc_codes,
    rule_sha256,
)


class TestParseAtc(unittest.TestCase):
    def test_complete_and_name(self):
        out = parse_atc_codes("N05BA04 Benzodiazepine derivatives")
        self.assertEqual(out["complete"], ["N05BA04"])
        self.assertEqual(out["partial"], [])

    def test_multiple_codes_with_names(self):
        out = parse_atc_codes(
            "L02AA03 Estrogens; G03CA01 Natural and semisynthetic estrogens, plain"
        )
        self.assertEqual(set(out["complete"]), {"L02AA03", "G03CA01"})

    def test_partial_five_character_code(self):
        out = parse_atc_codes("A03FA Propulsives")
        self.assertEqual(out["complete"], [])
        self.assertEqual(out["partial"], ["A03FA"])

    def test_partial_is_not_the_prefix_of_a_complete_code(self):
        # N05BA is the L4 prefix of N05BA04 and must not be reported as an
        # independent partial annotation.
        out = parse_atc_codes("N05BA04 Benzodiazepine derivatives")
        self.assertEqual(out["partial"], [])

    def test_empty_and_nan(self):
        self.assertEqual(parse_atc_codes(""), {"complete": [], "partial": []})
        self.assertEqual(parse_atc_codes(None), {"complete": [], "partial": []})


class TestAtcLevels(unittest.TestCase):
    def test_complete_only_projects_upwards(self):
        # L1..L5 = prefix lengths 1,3,4,5,7 as in the audit.
        counts = atc_level_counts(["A01AA01", "A01AA02"])
        self.assertEqual(counts["complete"]["L1"], 1)
        self.assertEqual(counts["complete"]["L2"], 1)
        self.assertEqual(counts["complete"]["L3"], 1)
        self.assertEqual(counts["complete"]["L4"], 1)
        self.assertEqual(counts["complete"]["L5"], 2)
        # partial reading still has only the complete codes here
        self.assertEqual(counts["complete_plus_partial"]["L4"], 1)

    def test_partial_adds_its_own_levels(self):
        counts = atc_level_counts(["A01AA01"], ["B05CA"])
        self.assertEqual(counts["complete"]["L1"], 1)
        self.assertEqual(counts["complete"]["L5"], 1)
        self.assertEqual(counts["complete_plus_partial"]["L1"], 2)
        self.assertEqual(counts["complete_plus_partial"]["L4"], 2)
        # a level-4 annotation cannot invent a level-5 substance code
        self.assertEqual(counts["complete_plus_partial"]["L5"], 1)

    def test_two_codes_sharing_a_group_count_once(self):
        counts = atc_level_counts(["A01AA01", "A01AA02"])
        self.assertEqual(counts["complete"]["L5"], 2)  # substances differ
        self.assertEqual(counts["complete"]["L3"], 1)  # A01A shared


class TestRuleIdentity(unittest.TestCase):
    def test_identity_unresolved_is_uncertain_even_with_a_phase(self):
        out = decide({"identity_resolved": False, "pubchem_max_phase": "Approved"})
        self.assertEqual(out["decision"], "uncertain")
        self.assertEqual(out["decision_class"], "unknown")


class TestRuleClasses(unittest.TestCase):
    def test_pubchem_approved_is_include_approved(self):
        out = decide({"identity_resolved": True, "pubchem_max_phase": "Approved"})
        self.assertEqual(out["decision"], "include")
        self.assertEqual(out["decision_class"], "approved")

    def test_broad_launched_is_include_approved(self):
        out = decide({"identity_resolved": True, "broad_phase": "Launched"})
        self.assertEqual(out["decision_class"], "approved")

    def test_phase_two_is_clinical(self):
        out = decide({"identity_resolved": True, "pubchem_max_phase": "Phase 2"})
        self.assertEqual(out["decision"], "include")
        self.assertEqual(out["decision_class"], "clinical")

    def test_withdrawn_is_clinical_include_with_flag(self):
        out = decide({"identity_resolved": True, "broad_phase": "Withdrawn"})
        self.assertEqual(out["decision"], "include")
        self.assertEqual(out["decision_class"], "clinical")
        self.assertTrue(out["flags"]["withdrawn"])

    def test_research_tool_needs_mechanism_and_tool_use(self):
        out = decide(
            {
                "identity_resolved": True,
                "broad_moa": "Hedgehog signalling pathway inhibitor",
                "broad_target": "SMO",
                "broad_hub_member": True,
            }
        )
        self.assertEqual(out["decision"], "include")
        self.assertEqual(out["decision_class"], "mechanism_tool")

    def test_mechanism_without_use_evidence_is_unknown(self):
        out = decide(
            {
                "identity_resolved": True,
                "mechanism_texts": ["Inhibits sodium channels"],
            }
        )
        self.assertEqual(out["decision"], "uncertain")
        self.assertEqual(out["decision_class"], "unknown")

    def test_pesticide_mechanism_alone_does_not_qualify(self):
        out = decide(
            {
                "identity_resolved": True,
                "mechanism_texts": ["Organophosphates poison insects by inhibiting "
                                    "acetylcholinesterase"],
                "drug_classes": ["Cholinesterase Inhibitors"],
                "use_classification": ["INSECTICIDES"],
            }
        )
        self.assertEqual(out["decision"], "exclude")
        self.assertEqual(out["decision_class"], "other_use")
        self.assertEqual(
            out["rationale_category"], "pesticide_or_agrochemical"
        )

    def test_human_atc_beats_industrial_use(self):
        out = decide(
            {
                "identity_resolved": True,
                "atc_human_codes": ["A01AD05"],
                "mechanism_texts": ["Counter-irritation"],
                "use_classification": ["Plasticisers"],
            }
        )
        self.assertEqual(out["decision"], "include")
        self.assertEqual(out["decision_class"], "mechanism_tool")

    def test_pharmaceutical_use_string_keeps_a_drug_with_a_fragrance_use(self):
        # methyl salicylate: approved and used as fragrance/flavouring, but the
        # documented pharmaceutical use keeps D1 from excluding it.
        out = decide(
            {
                "identity_resolved": True,
                "pubchem_max_phase": "Approved",
                "use_classification": [
                    "Fragrance Ingredients",
                    "Pharmaceuticals -> Animal Drugs -> Approved in Taiwan",
                ],
            }
        )
        self.assertEqual(out["decision"], "include")
        self.assertEqual(out["decision_class"], "approved")

    def test_generic_pharmaceuticals_string_does_not_promote_a_pesticide_tool(self):
        # Dicloran: "Pharmaceuticals" use + fungicide mechanism, no ATC, no Broad.
        out = decide(
            {
                "identity_resolved": True,
                "mechanism_texts": ["As a fungicide, it delays germination"],
                "use_classification": ["Pharmaceuticals"],
            }
        )
        self.assertEqual(out["decision"], "uncertain")

    def test_unknown_is_uncertain_never_excluded(self):
        out = decide({"identity_resolved": True})
        self.assertEqual(out["decision"], "uncertain")
        self.assertEqual(out["decision_class"], "unknown")

    def test_qp_atc_is_not_human_pharmacological_use(self):
        out = decide(
            {
                "identity_resolved": True,
                "atc_human_codes": [],
                "use_classification": ["Insecticides"],
            }
        )
        self.assertEqual(out["decision"], "exclude")


class TestRuleProvenance(unittest.TestCase):
    def test_rule_version_and_hash_are_stable(self):
        self.assertTrue(RULE_VERSION.startswith("eligibility-rule-"))
        self.assertRegex(rule_sha256(), r"^[0-9a-f]{64}$")
        self.assertEqual(rule_sha256(), rule_sha256())


if __name__ == "__main__":
    unittest.main()
