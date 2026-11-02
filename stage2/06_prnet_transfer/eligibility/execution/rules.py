"""Locked pharmacological eligibility rule for the temp_env 460 structures.

The rule is a pure function of documented evidence; it was written and hashed
before the final table was built and before any model was trained (no model is
trained in this ticket).  ``RULE_TEXT`` is the canonical, human-readable rule and
``rule_sha256()`` its hash, both recorded in the table, the summary and the
handoff.  Editing any clause of the text changes the hash and therefore the
version identity of every artifact that cites it.

Evidence used here is only what ``evidence.py`` extracted from documented
sources (PubChem, the Broad Repurposing Hub, BindingDB) plus the identity
resolution result.  Frozen ``compounds.csv`` ATC/MoA/MeSH strings are annotation
under audit, never eligibility evidence.
"""

import hashlib
import re

RULE_VERSION = "eligibility-rule-v1.2"

RULE_TEXT = """\
eligibility-rule-v1.2 - locked 2026-11-02, before the final table and before any
model was trained (this ticket trains nothing). v1.0 was a pre-final draft: it
took a documented clinical phase as sufficient, and PubChem's ChEMBL-derived Max
Phase marks several agrochemicals "Approved"/"Phase 2" (e.g. Carbaryl, Dieldrin),
so v1.1 runs the documented-other-use test first and adds an explicit human-use
escape. v1.2 additionally drops a PubChem "Drug Classes" string as sufficient
research-tool evidence, because "Cholinesterase Inhibitors"/"Insecticides" class
strings let agrochemicals in through D4. No model result existed at any revision.

Scope: every drug_key of the frozen stage2/02_cm_data/output/datasets-1/compounds.csv
with temp_env_all-smiles=True: 460 structure identities carried by 467 name records.

Evidence sources (per-row URLs and retrieval timestamps are kept in the table):
S1 PubChem PUG-REST + PUG-View record resolved from the frozen InChIKey: structure
   identity, synonyms, ATC, Max Phase, Drug Classes, Drug Indication, Therapeutic
   Uses, Mechanism of Action, Drug Targets, Use Classification.
S2 Broad Institute Repurposing Hub (clinical_phase, moa, target), matched by
   InChIKey or PubChem CID (never by name alone for a status claim).
S3 BindingDB 2024-09 ligand-target pairs, matched by Ligand InChI Key; curated
   binding evidence used only as supplementary target evidence.
The frozen compounds.csv ATC/MoA/MeSH strings are annotation under audit, not
eligibility evidence.  ChEMBL was attempted on 2026-11-02 and was unreachable
(HTTP 000 and HTTP 500); no ChEMBL field is used as evidence.

Decision procedure, in order (first match wins):
D0 identity gate: a key is resolved when the frozen InChIKey maps to a PubChem CID
   or to a Broad/BindingDB InChIKey record.  Unresolved identity => uncertain, no
   matter what other evidence exists.
D1 other_use: a documented non-pharmacological use (pesticide/agrochemical,
   industrial, cosmetic/fragrance/flavour, food additive) AND no human-use marker
   => exclude.  The documented use string is the reason.
D2 approved: PubChem Max Phase == "Approved" OR Broad clinical_phase == "Launched"
   => include.
D3 clinical: a documented clinical phase (Early Phase 1, Phase 1/2/3/4,
   Phase 1/Phase 2, Phase 2/Phase 3, Withdrawn) in either source => include.
D4 mechanism_tool: mechanism/target evidence AND (research-tool use evidence OR a
   human-use marker) => include.
D5 unknown: nothing above is documented => uncertain.

Two different markers are used because they answer different questions:
- documented pharmaceutical use (escape from D1) = a 7-character human ATC code
  (not the Q veterinary group) OR a PubChem Use Classification string naming
  pharmaceutical use (e.g. "Pharmaceuticals -> Animal Drugs").  This keeps a real
  drug that also carries a fragrance/food use (methyl salicylate, coumarin) out of
  the other-use exclude.
- human-use marker (D4 qualifier) = a 7-character human ATC code only.  The generic
  "Pharmaceuticals" use string is not enough to promote an agrochemical with a
  mechanism text into a mechanism tool (Dicloran, Prometon, Cyanazine).
A PubChem Drug Indication string is neither marker: when PubChem has one for a
pesticide it describes the pesticidal use (e.g. Allethrin).  Neither marker
substitutes for D2/D3: an approved drug whose PubChem record is sparse (e.g. a salt
form) is still included by D2, and the table flags such rows with
status_without_human_use_marker.  A compound whose only use strings are
agrochemical/industrial is excluded by D1 even if an aggregator records an approval
for it (Carbaryl, Lindane).
Mechanism/target evidence = Broad moa or target, PubChem Mechanism of Action text,
PubChem Drug Targets, or a BindingDB ligand-target pair.
Research-tool use evidence = Broad Repurposing Hub membership (a curated research
repurposing library). A PubChem Drug Classes / MeSH pharmacological-class string
alone is not research-tool evidence: it names a pharmacological action, not a use,
and pesticide class strings are exactly the leakage v1.1 admitted.

Decisions: include for D2/D3/D4, exclude for D1, uncertain for D0/D5.
"unknown" is never an exclude reason; a pesticide/herbicide mechanism alone does
not qualify a compound (D1 runs before D4).
"""

# Documented non-pharmacological use categories.  A term only counts when the
# string came from PubChem's Use Classification / Industry Uses evidence.
OTHER_USE_CATEGORIES = (
    (
        "pesticide_or_agrochemical",
        (
            "insecticide",
            "herbicide",
            "fungicide",
            "pesticide",
            "rodenticide",
            "nematicide",
            "acaricide",
            "molluscicide",
            "agrochemical",
            "plant growth regulator",
            "defoliant",
            "fumigant",
        ),
    ),
    (
        "industrial_chemical",
        (
            "plasticiser",
            "plasticizer",
            "flame retardant",
            "fluorosurfactant",
            "surfactant",
            "solvent",
            "lubricant",
            "dye",
            "pigment",
            "photographic",
            "explosive",
            "propellant",
            "refrigerant",
            "fuel",
            "adhesive",
            "coating",
            "textile",
            "plating",
            "corrosion inhibitor",
            "antifreeze",
            "hydraulic fluid",
            "ink",
            "polymer",
            "rubber",
            "preservative",
        ),
    ),
    (
        "cosmetic_or_fragrance",
        (
            "fragrance",
            "flavour",
            "flavor",
            "cosmetic",
            "perfume",
            "skin conditioning",
            "humectant",
        ),
    ),
    (
        "food_additive",
        (
            "food additive",
            "flavouring agent",
            "flavoring agent",
            "sweetener",
        ),
    ),
)

ATC_COMPLETE = re.compile(r"\b[A-Z][0-9]{2}[A-Z]{2}[0-9]{2}\b")
ATC_PARTIAL = re.compile(r"\b[A-Z][0-9]{2}[A-Z]{2}\b(?![0-9])")
ATC_GROUP = re.compile(r"\b[A-Z][0-9]{2}\b(?![A-Z0-9])")

APPROVED_PHASES = {"approved", "launched"}
CLINICAL_PHASES = {
    "early phase 1",
    "phase 1",
    "phase 2",
    "phase 3",
    "phase 4",
    "phase 1/phase 2",
    "phase 2/phase 3",
    "phase 1/ phase 2",
    "phase 2/ phase 3",
    "withdrawn",
}
NON_CLINICAL_PHASES = {"preclinical", "unknown", ""}

# ATC levels as reported in the audit: L1 = 1-character anatomical group,
# L2 = 3-character therapeutic subgroup, L3 = 4-character pharmacological
# subgroup, L4 = 5-character chemical subgroup, L5 = 7-character substance.
ATC_LEVELS = (("L1", 1), ("L2", 3), ("L3", 4), ("L4", 5), ("L5", 7))


def rule_sha256():
    """sha256 of the canonical rule text."""
    return hashlib.sha256(RULE_TEXT.encode("utf-8")).hexdigest()


def parse_atc_codes(text):
    """Split an ATC annotation string into complete and partial codes.

    The frozen field is ``CODE Name; CODE Name``; a code may be complete
    (7 characters, e.g. ``N02BA01``) or partial (5 characters, e.g. ``A03FA``).
    A 5-character code that is the prefix of a complete code in the same string
    is not reported again.
    """
    if text is None:
        return {"complete": [], "partial": []}
    text = str(text)
    if not text.strip():
        return {"complete": [], "partial": []}
    complete = ATC_COMPLETE.findall(text)
    partial = [code for code in ATC_PARTIAL.findall(text) if code not in complete]
    return {"complete": complete, "partial": partial}


def atc_level_counts(complete_codes, partial_codes=()):
    """Distinct prefix counts at levels L1/L3/L4/L5/L7 for both ATC readings.

    ``complete`` projects the 7-character substance codes upwards (the historical
    reading); ``complete_plus_partial`` also adds the 5-character annotations at
    the levels they support.  A 5-character annotation cannot create a level-7
    substance code, so ``L7`` is identical in both readings.
    """
    complete = {code for code in complete_codes if len(code) == 7}
    partial = {code for code in partial_codes if len(code) == 5}
    out = {"complete": {}, "complete_plus_partial": {}}
    for label, length in ATC_LEVELS:
        out["complete"][label] = len({code[:length] for code in complete})
    for label, length in ATC_LEVELS:
        codes = complete if length == 7 else complete | partial
        out["complete_plus_partial"][label] = len({code[:length] for code in codes})
    return out


def _norm_phase(value):
    return str(value or "").strip().lower().replace("  ", " ")


def _other_use_term(strings):
    """First (category, term, string) documented non-pharmacological use hit."""
    for raw in strings:
        text = str(raw or "").lower()
        if not text:
            continue
        for category, terms in OTHER_USE_CATEGORIES:
            for term in terms:
                if term in text:
                    return category, term, str(raw)
    return None


def decide(evidence):
    """Map one evidence record to (decision, class, rationale, flags).

    ``evidence`` is the plain dict produced by ``evidence.build_evidence``; only
    the documented keys are read.  Pure function: no I/O, no globals mutated.
    """
    def strings(key):
        value = evidence.get(key) or []
        if isinstance(value, str):
            return [value]
        return [str(item) for item in value if item]

    identity_resolved = bool(evidence.get("identity_resolved"))
    max_phase = _norm_phase(evidence.get("pubchem_max_phase"))
    broad_phase = _norm_phase(evidence.get("broad_phase"))
    atc_human = strings("atc_human_codes")
    indications = strings("drug_indications")
    drug_classes = strings("drug_classes")
    mechanism = (
        strings("mechanism_texts")
        + strings("drug_targets")
        + strings("bindingdb_targets")
    )
    broad_mechanism = [
        value
        for value in (evidence.get("broad_moa"), evidence.get("broad_target"))
        if value
    ]
    mechanism = mechanism + broad_mechanism
    tool_use = bool(evidence.get("broad_hub_member"))
    use_strings = strings("use_classification") + strings("industry_uses")
    human_pharm_use = bool(atc_human)
    pharmaceutical_use = human_pharm_use or any(
        "pharmaceutical" in item.lower() for item in use_strings
    )
    other_use = _other_use_term(use_strings)

    approved = max_phase in APPROVED_PHASES or broad_phase in APPROVED_PHASES
    clinical = (
        (max_phase in CLINICAL_PHASES and max_phase not in APPROVED_PHASES)
        or (broad_phase in CLINICAL_PHASES and broad_phase not in APPROVED_PHASES)
    )
    withdrawn = "withdrawn" in {max_phase, broad_phase}
    mechanism_evidence = bool(mechanism)

    flags = {
        "identity_resolved": identity_resolved,
        "approved": approved,
        "clinical": clinical,
        "withdrawn": withdrawn,
        "mechanism_evidence": mechanism_evidence,
        "tool_use": tool_use,
        "human_pharm_use": human_pharm_use,
        "pharmaceutical_use": pharmaceutical_use,
        "other_use": bool(other_use),
    }

    if not identity_resolved:
        return {
            "decision": "uncertain",
            "decision_class": "unknown",
            "rationale_category": "identity_unresolved",
            "rationale": "structure identity unresolved: frozen InChIKey maps to no "
            "PubChem/Broad/BindingDB record",
            "flags": flags,
        }

    if other_use and not pharmaceutical_use:
        category, term, raw = other_use
        return {
            "decision": "exclude",
            "decision_class": "other_use",
            "rationale_category": category,
            "rationale": f"other use: documented use {raw!r} (term {term!r}) and no "
            "human-use marker",
            "flags": flags,
        }

    if approved:
        return {
            "decision": "include",
            "decision_class": "approved",
            "rationale_category": "documented_approval",
            "rationale": "approved: documented status "
            + (f"PubChem Max Phase={evidence.get('pubchem_max_phase')}"
               if max_phase in APPROVED_PHASES
               else f"Broad clinical_phase={evidence.get('broad_phase')}"),
            "flags": flags,
        }

    if clinical:
        phase = (
            evidence.get("pubchem_max_phase")
            if max_phase in CLINICAL_PHASES
            else evidence.get("broad_phase")
        )
        return {
            "decision": "include",
            "decision_class": "clinical",
            "rationale_category": "documented_clinical_phase",
            "rationale": f"clinical: documented phase {phase}"
            + (" (withdrawn)" if withdrawn else ""),
            "flags": flags,
        }

    if mechanism_evidence and (tool_use or human_pharm_use):
        if broad_mechanism:
            basis = "Broad moa/target " + "/".join(str(v) for v in broad_mechanism)
        elif strings("mechanism_texts"):
            basis = "PubChem Mechanism of Action text"
        elif strings("drug_targets"):
            basis = "PubChem Drug Targets"
        else:
            basis = "BindingDB ligand-target pair"
        use = (
            "Broad Hub research-library membership"
            if evidence.get("broad_hub_member")
            else ("PubChem Drug Classes" if drug_classes else "human drug indication")
        )
        return {
            "decision": "include",
            "decision_class": "mechanism_tool",
            "rationale_category": "mechanism_plus_use",
            "rationale": f"mechanism_tool: mechanism from {basis}; use evidence from {use}",
            "flags": flags,
        }

    gap = (
        "mechanism evidence present but no research-tool or human pharmacological "
        "use evidence"
        if mechanism_evidence
        else "no documented status, mechanism or use"
    )
    return {
        "decision": "uncertain",
        "decision_class": "unknown",
        "rationale_category": "insufficient_evidence",
        "rationale": f"unknown: {gap}",
        "flags": flags,
    }
