"""Tests for the locked eligibility rule engine (rules.py).

Expectations come from the frozen annotation format (``compounds.csv`` ``atc``
field: ``CODE Name; CODE Name``) and from the v1.3 decision rule text, not from
the implementation.

The v1.3 precedence under test (ticket #1 review, 2026-11-02):

    documented human drug use  >  documented environmental/industrial use
                               >  mechanism text

A mechanism text never admits a compound on its own, and a bare approval /
clinical-phase string is not "human drug use" when the same record documents an
environmental or industrial use (the ChEMBL-derived Max Phase marks
agrochemicals "Approved").  Corroboration is a human ATC code, a human
therapeutic indication, an FDA-approved-drug listing or the WHO essential
medicines list.
"""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rules import (
    PRECEDENCE_DOMAIN,
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
        self.assertEqual(out["precedence"], "identity_unresolved")


class TestHumanDrugUsePrecedence(unittest.TestCase):
    """P1: documented human drug use outranks environmental/industrial use."""

    def test_approval_with_human_indication_beats_insecticide_use(self):
        # Lindane: FDA-approved scabicide whose record also documents insecticide
        # use.  The human indication corroborates the approval, so human use wins.
        out = decide(
            {
                "identity_resolved": True,
                "pubchem_max_phase": "Approved",
                "broad_phase": "Launched",
                "broad_indication": "scabies",
                "use_classification": ["INSECTICIDES"],
            }
        )
        self.assertEqual(out["decision"], "include")
        self.assertEqual(out["decision_class"], "approved")
        self.assertEqual(out["precedence"], "human_drug_use")
        self.assertIn("scabies", out["human_use_basis"])

    def test_bare_status_on_a_documented_pesticide_is_not_human_use(self):
        # Carbaryl: Max Phase "Approved" is ChEMBL-derived and the recorded use is
        # an insecticide with no human ATC/indication.  A bare status string is not
        # "documented human drug use", so the pesticide use excludes it.
        out = decide(
            {
                "identity_resolved": True,
                "pubchem_max_phase": "Approved",
                "broad_phase": "Launched",
                "use_classification": ["INSECTICIDES"],
            }
        )
        self.assertEqual(out["decision"], "exclude")
        self.assertEqual(out["decision_class"], "other_use")
        self.assertEqual(out["precedence"], "environmental_or_industrial_use")
        self.assertTrue(out["other_use_basis"])

    def test_clinical_phase_with_human_indication_beats_food_additive_use(self):
        # Quercetin: documented clinical phase + human indications + a food-additive
        # use.  v1.2 ordered the other-use test first and excluded it; v1.3 decides
        # by the human evidence.
        out = decide(
            {
                "identity_resolved": True,
                "pubchem_max_phase": "Phase 3",
                "broad_indication": "fatigue|allergic rhinitis|drowsiness",
                "use_classification": ["Food Additives -> COLOUR"],
            }
        )
        self.assertEqual(out["decision"], "include")
        self.assertEqual(out["decision_class"], "clinical")
        self.assertEqual(out["precedence"], "human_drug_use")

    def test_approved_phase_conflict_is_approved(self):
        # Either source documenting approval is enough (unchanged from v1.2).
        out = decide({"identity_resolved": True, "broad_phase": "Launched"})
        self.assertEqual(out["decision_class"], "approved")

    def test_clinical_phase_alone_is_include(self):
        out = decide({"identity_resolved": True, "pubchem_max_phase": "Phase 2"})
        self.assertEqual(out["decision"], "include")
        self.assertEqual(out["decision_class"], "clinical")
        self.assertEqual(out["precedence"], "human_drug_use")

    def test_withdrawn_is_clinical_include_with_flag(self):
        out = decide({"identity_resolved": True, "broad_phase": "Withdrawn"})
        self.assertEqual(out["decision"], "include")
        self.assertEqual(out["decision_class"], "clinical")
        self.assertTrue(out["flags"]["withdrawn"])

    def test_approved_salt_without_marker_and_without_other_use_is_include(self):
        # A sparse salt-form record (moxifloxacin hydrochloride) has no other-use
        # evidence; the status claim is not contradicted, so it stays approved.
        out = decide({"identity_resolved": True, "pubchem_max_phase": "Approved"})
        self.assertEqual(out["decision"], "include")
        self.assertEqual(out["decision_class"], "approved")


class TestEnvironmentalUseOutranksMechanism(unittest.TestCase):
    """P2: documented environmental/industrial use excludes before mechanism."""

    def test_mechanism_tool_with_documented_herbicide_use_is_excluded(self):
        # Daminozide: Broad Hub Preclinical entry with a KDM2A MoA, but its record
        # documents plant-growth-regulator (pesticide) use.
        out = decide(
            {
                "identity_resolved": True,
                "broad_phase": "Preclinical",
                "broad_moa": "KDM2A inhibitor",
                "broad_target": "KDM2A|KDM7A|PHF8",
                "broad_hub_member": True,
                "use_classification": ["Environmental transformation -> Pesticides"],
            }
        )
        self.assertEqual(out["decision"], "exclude")
        self.assertEqual(out["decision_class"], "other_use")
        self.assertEqual(out["precedence"], "environmental_or_industrial_use")

    def test_bare_uses_text_counts_as_use_evidence(self):
        # The HSDB/EPA bare "Uses" statements now count as use evidence.
        out = decide(
            {
                "identity_resolved": True,
                "broad_phase": "Preclinical",
                "broad_moa": "mitochondrial electron transport inhibitor",
                "broad_hub_member": True,
                "uses_texts": [
                    "Used as broad-spectrum fungicide for cereals, fruits, vegetables"
                ],
            }
        )
        self.assertEqual(out["decision"], "exclude")
        self.assertEqual(out["decision_class"], "other_use")

    def test_human_atc_alone_does_not_rescue_industrial_use(self):
        # Dibutyl phthalate: Broad Preclinical tool entry + plasticizer use + a
        # (repellent) ATC code.  In v1.2 the ATC marker escaped the other-use test;
        # v1.3 requires documented human drug use (status or indication), not an
        # ATC annotation alone.
        out = decide(
            {
                "identity_resolved": True,
                "broad_phase": "Preclinical",
                "broad_target": "TRPA1",
                "broad_hub_member": True,
                "atc_human_codes": ["P03BX03"],
                "use_classification": ["PLASTICIZERS"],
            }
        )
        self.assertEqual(out["decision"], "exclude")
        self.assertEqual(out["decision_class"], "other_use")

    def test_q_atc_is_not_human_pharmacological_use(self):
        out = decide(
            {
                "identity_resolved": True,
                "atc_human_codes": [],
                "use_classification": ["Insecticides"],
            }
        )
        self.assertEqual(out["decision"], "exclude")


class TestMechanismTool(unittest.TestCase):
    """P3: a mechanism tool needs a mechanism *and* documented tool/human use."""

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
        self.assertEqual(out["precedence"], "mechanism_tool")

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
        self.assertEqual(out["rationale_category"], "pesticide_or_agrochemical")

    def test_mechanism_plus_human_atc_and_no_other_use_is_mechanism_tool(self):
        out = decide(
            {
                "identity_resolved": True,
                "atc_human_codes": ["S01EB10"],
                "mechanism_texts": ["Cholinesterase inhibitor"],
            }
        )
        self.assertEqual(out["decision"], "include")
        self.assertEqual(out["decision_class"], "mechanism_tool")

    def test_unknown_is_uncertain_never_excluded(self):
        out = decide({"identity_resolved": True})
        self.assertEqual(out["decision"], "uncertain")
        self.assertEqual(out["decision_class"], "unknown")
        self.assertEqual(out["precedence"], "insufficient_evidence")


class TestHumanUseMarkers(unittest.TestCase):
    """A corroborating marker must document *human drug* use."""

    def test_pharmaceutical_use_string_is_not_a_human_use_marker(self):
        # The generic PubChem "Pharmaceuticals" use-classification string is not
        # evidence of human drug use (v1.2 already refused it as tool evidence).
        out = decide(
            {
                "identity_resolved": True,
                "pubchem_max_phase": "Approved",
                "use_classification": ["Pharmaceuticals"],
                "uses_texts": ["Used as a pesticide."],
            }
        )
        self.assertEqual(out["decision"], "exclude")
        self.assertEqual(out["decision_class"], "other_use")

    def test_veterinary_indication_is_not_a_human_use_marker(self):
        # Fipronil: "flea control" is a veterinary indication, not human use.
        out = decide(
            {
                "identity_resolved": True,
                "broad_phase": "Launched",
                "broad_indication": "flea control",
                "use_classification": ["Insecticides"],
            }
        )
        self.assertEqual(out["decision"], "exclude")
        self.assertEqual(out["decision_class"], "other_use")

    def test_cosmetic_indication_is_not_a_human_use_marker(self):
        # Sulisobenzone/methylparaben: "sunscreen lotion"/"cosmetic" indications are
        # cosmetic, not therapeutic, evidence.
        out = decide(
            {
                "identity_resolved": True,
                "broad_phase": "Launched",
                "broad_indication": "sunscreen lotion",
                "use_classification": ["Cosmetics -> Uv absorber"],
            }
        )
        self.assertEqual(out["decision"], "exclude")
        self.assertEqual(out["decision_class"], "other_use")

    def test_pubchem_drug_indication_text_is_not_a_human_use_marker(self):
        # Fenoxaprop-ethyl: PubChem carries a mis-mapped human indication for this
        # herbicide.  The rule records the PubChem indication but does not use it as
        # human-use evidence (deviation D5).
        out = decide(
            {
                "identity_resolved": True,
                "pubchem_max_phase": "Phase 2",
                "drug_indications": ["Chronic angina and Coronary Artery Disease"],
                "uses_texts": ["Used as selective post-emergence herbicide"],
            }
        )
        self.assertEqual(out["decision"], "exclude")
        self.assertEqual(out["decision_class"], "other_use")

    def test_allergenic_patch_test_indication_is_not_a_human_use_marker(self):
        # 1,3-Diphenylguanidine: "approved for use within allergenic patch tests" is
        # a diagnostic reagent, not a documented human drug use.
        out = decide(
            {
                "identity_resolved": True,
                "drug_indications": [
                    "Diphenylguanidine is approved for use within allergenic "
                    "epicutaneous patch tests"
                ],
                "use_classification": ["PLASTICIZERS"],
            }
        )
        self.assertEqual(out["decision"], "exclude")


class TestRuleProvenance(unittest.TestCase):
    def test_rule_version_and_hash_are_stable(self):
        self.assertTrue(RULE_VERSION.startswith("eligibility-rule-"))
        self.assertEqual(RULE_VERSION, "eligibility-rule-v1.3")
        self.assertRegex(rule_sha256(), r"^[0-9a-f]{64}$")
        self.assertEqual(rule_sha256(), rule_sha256())

    def test_precedence_domain_is_the_documented_set(self):
        self.assertEqual(
            PRECEDENCE_DOMAIN,
            {
                "identity_unresolved",
                "human_drug_use",
                "environmental_or_industrial_use",
                "mechanism_tool",
                "insufficient_evidence",
            },
        )


if __name__ == "__main__":
    unittest.main()
