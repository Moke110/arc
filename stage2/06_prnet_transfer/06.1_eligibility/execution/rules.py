"""Locked pharmacological eligibility rule for the temp_env 460 structures.

The rule is a pure function of documented evidence; ``RULE_TEXT`` is the
canonical, human-readable rule and ``rule_sha256()`` its hash, both recorded in
every artifact.  Editing any clause of the text changes the hash and therefore
the version identity of every artifact that cites it.

v1.3 is a revision of the locked v1.2 rule (ticket #1 review, 2026-11-02).  It
makes the precedence explicit — documented human drug use outranks documented
environmental/industrial use, which in turn outranks mechanism text — and it
re-decides every affected row.  No model has been trained: the revision is still
pre-registration, and the v1.2 record is preserved unchanged beside v1.3.

Evidence used here is only what ``evidence.py`` extracted from documented
sources (PubChem, the Broad Repurposing Hub, BindingDB) plus the identity
resolution result.  Frozen ``compounds.csv`` ATC/MoA/MeSH strings are annotation
under audit, never eligibility evidence.
"""

import hashlib
import re

RULE_VERSION = "eligibility-rule-v1.3"

RULE_TEXT = """\
eligibility-rule-v1.3 - locked 2026-11-02, before the frozen input is built and
before any model is trained (this ticket trains nothing). Revision of the locked
eligibility-rule-v1.2 (sha256 e1af8cdb...), reopened by the ticket #1 review,
which found two rule-consistency defects in v1.2:
  (a) the mechanism-tool clause admitted agrochemical/industrial compounds
      (Ethoxyquin, Dibutyl phthalate, 2-Phenylphenol, Daminozide,
      Trifloxystrobin, Cyfluthrin, Paraoxon) on a mechanism text alone, which
      this ticket forbids ("a pesticide mechanism does not automatically satisfy
      eligibility"); v1.2 had read Broad Repurposing Hub membership as
      research-tool use and a bare ATC code as a human-use marker;
  (b) the approved/clinical-versus-other-use precedence was not stated and was
      applied inconsistently to natural products with registered human evidence
      (Quercetin was excluded as food additive while Genistein was clinical).
v1.3 states the precedence, encodes it, and re-decides every one of the 460
drug_keys; the v1.2 record is kept byte-identical beside it. No model result
existed at any revision, so this is pre-registration development.

Scope: every drug_key of the frozen stage2/02_cm_data/output/datasets-1/compounds.csv
with temp_env_all-smiles=True: 460 structure identities carried by 467 name records.

Evidence sources (per-row URLs and retrieval timestamps are kept in the table):
S1 PubChem PUG-REST + PUG-View record resolved from the frozen InChIKey: structure
   identity, synonyms, ATC, Max Phase, Drug Classes, Drug Indication, Therapeutic
   Uses, Mechanism of Action, Drug Targets, Use Classification, Industry Uses, and
   the bare Use and Manufacturing/Uses statements (HSDB/EPA "Uses" text).
S2 Broad Institute Repurposing Hub (clinical_phase, moa, target, indication,
   disease_area), matched by InChIKey or PubChem CID (never by name alone for a
   status claim).
S3 BindingDB 2024-09 ligand-target pairs, matched by Ligand InChI Key; curated
   binding evidence used only as supplementary target evidence.  In the v1.3
   build the join produced no usable per-key target pairs (recorded in
   summary.json); mechanism evidence therefore comes from S1 and S2.
The frozen compounds.csv ATC/MoA/MeSH strings are annotation under audit, not
eligibility evidence.  ChEMBL was attempted on 2026-11-02 and was unreachable
(HTTP 000 and HTTP 500); no ChEMBL field is used as evidence.

Documented human drug use (the highest-precedence evidence) is any of:
  H1 approval: PubChem Max Phase == "Approved" OR Broad clinical_phase == "Launched";
  H2 clinical phase: Early Phase 1, Phase 1/2/3/4, Phase 1/Phase 2,
     Phase 2/Phase 3, or Withdrawn in either source;
  H3 human indication: a Broad `indication` entry that is a human therapeutic
     indication (not itself an environmental/industrial use, not a cosmetic or
     fragrance use, not a veterinary/animal indication);
plus corroboration by a human-use marker: a 7-character human ATC code (not the
Q veterinary group), a human therapeutic indication (H3), a PubChem
"FDA Approved Drugs" listing, or a PubChem "WHO Essential Medicines" listing.
H1/H2 outrank a documented environmental/industrial use ONLY when corroborated by
a human-use marker: a bare status string is not human drug use when the same
record documents an environmental/industrial use, because PubChem Max Phase is
ChEMBL-derived and marks agrochemicals "Approved"/"Phase 2" (Carbaryl, Dieldrin).
H1/H2 are sufficient on their own when no environmental/industrial use is
documented (this keeps approved salt forms with sparse records, e.g. moxifloxacin
hydrochloride, included by status as in v1.2).

Documented environmental/industrial use (the middle precedence) is a
non-pharmacological use in PubChem Use Classification, PubChem Industry Uses, or
the bare PubChem Uses statements, naming a pesticide/agrochemical (insecticide,
herbicide, fungicide, pesticide, rodenticide, nematicide, acaricide,
molluscicide, agrochemical, plant growth regulator, defoliant, fumigant), an
industrial chemical (plasticiser, flame retardant, surfactant, solvent,
lubricant, dye, pigment, photographic, explosive, propellant, refrigerant, fuel,
adhesive, coating, textile, plating, corrosion inhibitor, antifreeze, hydraulic
fluid, ink, polymer, rubber, preservative), a cosmetic/fragrance/flavour use, or
a food additive.  When a record documents several such uses, the reported category
follows the declaration order pesticide/agrochemical, industrial,
cosmetic/fragrance, food additive, and the matched string is kept in the row.

Decision procedure, in order (first match wins). The precedence is
human drug use > documented environmental/industrial use > mechanism text:
D0 identity gate: a key is resolved when the frozen InChIKey maps to a PubChem CID
   or to a Broad/BindingDB InChIKey record.  Unresolved identity => uncertain, no
   matter what other evidence exists.
D1 human drug use: H1 or H2 is documented AND (no environmental/industrial use is
   documented OR a human-use marker corroborates it) => include, class approved
   for H1 else clinical.  This is where documented human evidence outranks a
   documented food-additive/fragrance/industrial use (Quercetin, Lindane, Eugenol).
D2 environmental/industrial use: a documented use above AND no D1 match =>
   exclude other_use, whatever mechanism text, Broad Hub membership or ATC code
   the record also carries.  Mechanism text alone never satisfies eligibility.
D3 mechanism_tool: mechanism/target evidence AND (Broad Hub research-library
   membership OR a human-use marker) => include.  A mechanism tool therefore
   needs a documented use/tool marker in addition to the mechanism; a bare
   mechanism text is not enough.
D4 unknown: nothing above is documented => uncertain.

Decisions: include for D1/D3, exclude for D2, uncertain for D0/D4.
"unknown" is never an exclude reason. A PubChem Drug Indication / Therapeutic
Uses string and a historical or experimental treatment note are recorded in the
table but are NOT human-use markers in v1.3: PubChem carries diagnostic
patch-test reagents ("approved for use within allergenic patch tests"), mis-mapped
entries (the herbicide fenoxaprop-ethyl carries a human angina indication) and
historical notes (paraoxon "formerly used for glaucoma"), so a status claim
corroborated only by such a string stays uncorroborated. The generic PubChem
"Pharmaceuticals" use-classification string is likewise not a human-use marker.
An ATC code alone is annotation, not a use: it does not escape D2 and is not a
substitute for H1/H2/H3. Neither marker substitutes for H1/H2 when an
environmental/industrial use is documented.
"""

# The precedence labels recorded per row and enforced by the build gates.
PRECEDENCE_DOMAIN = {
    "identity_unresolved",
    "human_drug_use",
    "environmental_or_industrial_use",
    "mechanism_tool",
    "insufficient_evidence",
}

# Documented non-pharmacological use categories.  A term only counts when the
# string came from PubChem's Use Classification / Industry Uses / Uses evidence.
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

# A Broad `indication` string carrying one of these is not human drug-use
# evidence: it is a cosmetic/fragrance use or a veterinary/animal indication.
# ("lice", "pediculosis" and "scabies" are human indications and stay markers.)
NON_HUMAN_INDICATION_TERMS = (
    "sunscreen",
    "lotion",
    "cosmetic",
    "perfume",
    "fragrance",
    "flavour",
    "flavor",
    "flea",
    "tick",
    "veterinary",
    "veterinarian",
    "animal",
    "livestock",
    "cattle",
    "poultry",
    "swine",
    "canine",
    "feline",
    "equine",
    "fescue",
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


def _strings(evidence, key):
    value = evidence.get(key) or []
    if isinstance(value, str):
        return [value]
    return [str(item) for item in value if item]


def _other_use_term(strings):
    """Documented non-pharmacological use hit, by category priority.

    Categories are scanned in declaration order (pesticide/agrochemical first,
    then industrial, cosmetic/fragrance, food additive) across all of the
    record's use strings, so a compound that documents both a pesticide and a
    food use is attributed to the pesticide (Ethoxyquin, 2-Phenylphenol) rather
    than to whichever string happened to come first in the record.
    """
    texts = [(str(raw or ""), str(raw or "").lower()) for raw in strings]
    for category, terms in OTHER_USE_CATEGORIES:
        for raw, text in texts:
            if not text:
                continue
            for term in terms:
                if term in text:
                    return category, term, raw
    return None


def _human_indication(strings):
    """First Broad indication entry that documents a human therapeutic use."""
    for raw in strings:
        for part in re.split(r"[|;]", str(raw or "")):
            part = part.strip()
            if not part:
                continue
            if _other_use_term([part]):
                continue
            lowered = part.lower()
            if any(term in lowered for term in NON_HUMAN_INDICATION_TERMS):
                continue
            return part
    return None


def human_use_marker(evidence):
    """Return the corroborating human-use marker string ('' when none).

    A marker says "this record documents a human drug use": a human ATC code, a
    human therapeutic indication (Broad), an FDA-approved-drug listing, or the
    WHO essential-medicines listing.  ATC is a marker here but is not a status:
    on its own it cannot turn a documented environmental/industrial use into an
    include (rule D2).
    """
    parts = []
    codes = _strings(evidence, "atc_human_codes")
    if codes:
        parts.append("atc:" + ",".join(codes[:3]))
    indication = _human_indication(_strings(evidence, "broad_indication"))
    if indication:
        parts.append(f"indication:{indication}")
    if _strings(evidence, "fda_approved_drugs"):
        parts.append("fda_approved_drugs")
    if _strings(evidence, "who_essential_medicines"):
        parts.append("who_essential_medicines")
    return "; ".join(parts)


def _status_basis(evidence):
    """Human-readable status evidence; returns (approved, clinical, withdrawn, basis)."""
    approved_parts, clinical_parts = [], []
    for label, value in (
        ("PubChem Max Phase", evidence.get("pubchem_max_phase")),
        ("Broad clinical_phase", evidence.get("broad_phase")),
    ):
        normalized = _norm_phase(value)
        if not normalized:
            continue
        if normalized in APPROVED_PHASES:
            approved_parts.append(f"{label}={value}")
        elif normalized in CLINICAL_PHASES:
            clinical_parts.append(f"{label}={value}")
    approved = bool(approved_parts)
    clinical = bool(clinical_parts) and not approved
    withdrawn = bool(
        {"withdrawn"} & {_norm_phase(evidence.get("pubchem_max_phase")),
                         _norm_phase(evidence.get("broad_phase"))}
    )
    if approved:
        basis = "approval: " + "; ".join(approved_parts)
    elif clinical:
        basis = "clinical phase: " + "; ".join(clinical_parts)
    else:
        basis = ""
    return approved, clinical, withdrawn, basis


def _mechanism_basis(evidence):
    """Human-readable mechanism evidence string ('' when none)."""
    if evidence.get("broad_moa"):
        basis = "Broad moa " + str(evidence["broad_moa"])
        if evidence.get("broad_target"):
            basis += f" (target {evidence['broad_target']})"
        return basis
    if _strings(evidence, "mechanism_texts"):
        return "PubChem Mechanism of Action text"
    if _strings(evidence, "drug_targets"):
        return "PubChem Drug Targets"
    if _strings(evidence, "bindingdb_targets"):
        return "BindingDB ligand-target pair"
    return ""


def _mechanism_evidence(evidence):
    """Mechanism/target evidence from the documented sources."""
    return bool(
        evidence.get("broad_moa")
        or evidence.get("broad_target")
        or _strings(evidence, "mechanism_texts")
        or _strings(evidence, "drug_targets")
        or _strings(evidence, "bindingdb_targets")
    )


def decide(evidence):
    """Map one evidence record to a v1.3 decision dict.

    ``evidence`` is the plain dict produced by ``evidence.build_evidence``; only
    the documented keys are read.  Pure function: no I/O, no globals mutated.

    The returned dict carries the decision, its class, the precedence label that
    produced it, the per-row evidence string for the deciding clause, and flags.
    """
    identity_resolved = bool(evidence.get("identity_resolved"))
    approved, clinical, withdrawn, status_basis = _status_basis(evidence)
    marker = human_use_marker(evidence)
    other_use = _other_use_term(
        _strings(evidence, "use_classification")
        + _strings(evidence, "industry_uses")
        + _strings(evidence, "uses_texts")
    )
    mechanism = _mechanism_evidence(evidence)
    mechanism_basis = _mechanism_basis(evidence)
    tool_use = bool(evidence.get("broad_hub_member"))
    atc_human = bool(_strings(evidence, "atc_human_codes"))
    human_indication = _human_indication(_strings(evidence, "broad_indication"))

    other_use_basis = ""
    if other_use:
        category, term, raw = other_use
        other_use_basis = f"{category} ({term}): {raw}"

    flags = {
        "identity_resolved": identity_resolved,
        "approved": approved,
        "clinical": clinical,
        "withdrawn": withdrawn,
        "mechanism_evidence": mechanism,
        "tool_use": tool_use,
        "human_use_marker": bool(marker),
        "atc_human": atc_human,
        "human_indication": human_indication is not None,
        "other_use": bool(other_use),
    }
    base = {
        "flags": flags,
        "human_use_basis": "",
        "other_use_basis": other_use_basis,
    }

    if not identity_resolved:
        return {
            **base,
            "decision": "uncertain",
            "decision_class": "unknown",
            "precedence": "identity_unresolved",
            "rationale_category": "identity_unresolved",
            "rationale": "identity unresolved: frozen InChIKey maps to no "
            "PubChem/Broad/BindingDB record",
        }

    # D1 human drug use outranks environmental/industrial use.
    human_status = approved or clinical
    if human_status and (not other_use or marker):
        basis = status_basis
        if marker:
            basis += f"; corroborated by {marker}"
        if other_use:
            basis += (
                f"; documented {other_use[0]} use ({other_use[1]!r}) is outranked "
                "by documented human drug use"
            )
        return {
            **base,
            "decision": "include",
            "decision_class": "approved" if approved else "clinical",
            "precedence": "human_drug_use",
            "rationale_category": (
                "documented_approval" if approved else "documented_clinical_phase"
            ),
            "rationale": ("approved: " if approved else "clinical: ") + basis,
            "human_use_basis": basis,
        }

    # D2 environmental/industrial use, with no documented human drug use.
    if other_use:
        category, term, raw = other_use
        if human_status:
            reason = (
                "status claim is not corroborated as human drug use (no human ATC/"
                "indication/FDA/WHO marker)"
            )
        elif mechanism:
            reason = "a mechanism text does not satisfy eligibility"
        elif atc_human:
            reason = "an ATC code alone is not a human drug use"
        else:
            reason = "no documented human drug use"
        return {
            **base,
            "decision": "exclude",
            "decision_class": "other_use",
            "precedence": "environmental_or_industrial_use",
            "rationale_category": category,
            "rationale": f"other use: documented use {raw!r} (term {term!r}); "
            + reason,
        }

    # D3 mechanism tool: mechanism AND a documented tool/human-use marker.
    if mechanism and (tool_use or marker):
        use = (
            "Broad Hub research-library membership"
            if tool_use
            else f"human-use marker {marker}"
        )
        basis = f"mechanism: {mechanism_basis}; use: {use}"
        return {
            **base,
            "decision": "include",
            "decision_class": "mechanism_tool",
            "precedence": "mechanism_tool",
            "rationale_category": "mechanism_plus_use",
            "rationale": f"mechanism_tool: {basis}",
            "human_use_basis": basis,
        }

    gap = (
        "mechanism evidence present but no research-tool or human-use marker"
        if mechanism
        else "no documented status, mechanism or use"
    )
    return {
        **base,
        "decision": "uncertain",
        "decision_class": "unknown",
        "precedence": "insufficient_evidence",
        "rationale_category": "insufficient_evidence",
        "rationale": f"unknown: {gap}",
    }
