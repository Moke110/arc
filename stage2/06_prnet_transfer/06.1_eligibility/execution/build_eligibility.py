"""Build the 460-structure pharmacological eligibility audit for temp_env.

Reads the frozen ``compounds.csv`` (identity + annotation under audit), collects
documented evidence (PubChem, Broad Repurposing Hub, BindingDB), applies the
locked rule in ``rules.py``, and writes the evidence table, coverage summary, the
machine-readable handoff for ticket #3, the v1.2 -> v1.3 change report and the
hashed artifact manifest.

Artifacts are versioned by rule version under ``output/<version>/``; the frozen
v1.2 record is never overwritten (ADR 0030) and ``output/current.json`` is the
machine-readable pointer to the current record.

    python execution/build_eligibility.py --out-root stage2/06_prnet_transfer
    python execution/build_eligibility.py --offline          # cached evidence only
    python execution/build_eligibility.py --check            # gates on the record
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, OrderedDict
from pathlib import Path

import evidence
import rules

ARTIFACT_VERSION = "v1.3"
PREVIOUS_VERSION = "v1.2"
SCHEMA_TABLE = "prnet-eligibility-table-v2"
SCHEMA_SUMMARY = "prnet-eligibility-summary-v2"
SCHEMA_HANDOFF = "prnet-eligibility-handoff-v2"
SCHEMA_MANIFEST = "prnet-eligibility-manifest-v2"
SCHEMA_CHANGE_REPORT = "prnet-eligibility-change-report-v1"
SCHEMA_CURRENT = "prnet-eligibility-current-v1"
COMPOUNDS_CSV = Path("stage2/02_cm_data/output/datasets-1/compounds.csv")
DEFAULT_BINDINGDB = evidence.DEFAULT_BINDINGDB

# Files that belong to one artifact version.
VERSION_FILES = (
    "eligibility_table.csv",
    "old_tier_comparison.csv",
    "evidence_manifest.csv",
    "summary.json",
    "handoff.json",
    "change_report.csv",
    "manifest.json",
)

TABLE_COLUMNS = [
    "drug_key",
    "all_names",
    "n_name_records",
    "compound_name",
    "inchikey",
    "smiles",
    "pubchem_cid",
    "decision",
    "decision_class",
    "decision_precedence",
    "rationale_category",
    "rationale",
    "human_use_basis",
    "other_use_basis",
    "v12_decision",
    "v12_decision_class",
    # documented status
    "clinical_phase",
    "clinical_phase_source",
    "phase_conflict",
    "pubchem_max_phase",
    "broad_phase",
    "broad_match",
    "withdrawn",
    # mechanism / targets
    "moa",
    "moa_sources",
    "moa_frozen",
    "targets",
    "target_sources",
    "drug_classes",
    # ATC / MeSH annotation under audit
    "atc_codes",
    "atc_frozen_complete",
    "atc_frozen_partial",
    "atc_pubchem",
    "atc_human_codes",
    "mesh_frozen",
    "mesh_corroborated",
    # use evidence
    "use_classification",
    "industry_uses",
    "uses_texts",
    # PubChem indication text, recorded but not human-use evidence (deviation D5)
    "drug_indications",
    "therapeutic_uses",
    # flags and provenance
    "flag_mechanism_evidence",
    "flag_tool_use",
    "flag_human_use_marker",
    "flag_atc_human",
    "flag_human_indication",
    "flag_other_use",
    "status_without_human_use_marker",
    "identity_resolved",
    "identity_notes",
    "name_in_pubchem_synonyms",
    "smiles_key",
    "bindingdb_n_pairs",
    "bindingdb_targets",
    "old_tier_atc_moa",
    "old_tier_atc_moa_mesh",
    "old_tier_all_smiles",
    "evidence_urls",
    "retrieved_at",
    "evidence_ok_urls",
]


def _join(items, sep="; "):
    return sep.join(str(item) for item in items if str(item).strip())


def _truncate(text, limit=400):
    text = str(text or "").strip()
    return text if len(text) <= limit else text[: limit - 1].rstrip() + "…"


def _mesh_corroborated(frozen_mesh, haystack):
    tokens = [token.strip() for token in str(frozen_mesh or "").split(";") if token.strip()]
    if not tokens:
        return ""
    found = 0
    for token in tokens:
        normalized = evidence.normalize_name(token)
        if len(normalized) >= 4 and normalized in haystack:
            found += 1
    return f"{found}/{len(tokens)}"


def assemble(keys, evidence_map, previous):
    """Apply the locked rule and flatten each key into one table row.

    ``previous`` is the frozen v1.2 table keyed by drug_key; its decision and
    class are carried as ``v12_*`` columns so every later change is traceable.
    """
    rows = []
    for key in keys:
        ev = evidence_map[key["drug_key"]]
        verdict = rules.decide(ev)
        flags = verdict["flags"]
        prior = previous.get(key["drug_key"], {})

        frozen_complete = rules.parse_atc_codes(key["atc_frozen"])["complete"]
        frozen_partial = rules.parse_atc_codes(key["atc_frozen"])["partial"]
        pubchem_codes = list(ev.get("atc_pubchem_codes", []))
        haystack = " ".join(
            str(item).lower()
            for field in ("drug_classes", "therapeutic_uses", "mechanism_texts")
            for item in (ev.get(field) or [])
        )
        for item in ev.get("drug_indications") or []:
            haystack += " " + str(item).lower()

        moa, moa_sources = "", []
        if ev.get("broad_moa"):
            moa, moa_sources = ev["broad_moa"], moa_sources + ["broad_repurposing_hub"]
        elif ev.get("mechanism_texts"):
            moa = _truncate(ev["mechanism_texts"][0])
            moa_sources = ["pubchem_mechanism_of_action"]

        targets, target_sources = [], []
        if ev.get("broad_target"):
            targets.extend(part for part in ev["broad_target"].split("|") if part.strip())
            target_sources.append("broad_repurposing_hub")
        if ev.get("drug_targets"):
            targets.extend(ev["drug_targets"])
            target_sources.append("pubchem_drug_targets")
        if ev.get("bindingdb_targets"):
            targets.extend(ev["bindingdb_targets"])
            target_sources.append("bindingdb")
        targets = list(OrderedDict.fromkeys(t for t in targets if str(t).strip()))

        def phase_category(value):
            normalized = str(value or "").strip().lower()
            if normalized in rules.APPROVED_PHASES:
                return "approved"
            if normalized in rules.CLINICAL_PHASES:
                return "clinical"
            return "other"

        phase_conflict = bool(
            ev.get("pubchem_max_phase")
            and ev.get("broad_phase")
            and phase_category(ev.get("pubchem_max_phase"))
            != phase_category(ev.get("broad_phase"))
        )

        phase = ev.get("pubchem_max_phase") or ev.get("broad_phase") or ""
        phase_source = (
            "pubchem_max_phase"
            if ev.get("pubchem_max_phase")
            else ("broad_repurposing_hub" if ev.get("broad_phase") else "")
        )

        rows.append(
            {
                "drug_key": key["drug_key"],
                "all_names": key["all_names"],
                "n_name_records": key["n_name_records"],
                "compound_name": key["compound_name"],
                "inchikey": key["inchikey"],
                "smiles": key["smiles"],
                "pubchem_cid": ev.get("pubchem_cid") or "",
                "decision": verdict["decision"],
                "decision_class": verdict["decision_class"],
                "decision_precedence": verdict["precedence"],
                "rationale_category": verdict["rationale_category"],
                "rationale": verdict["rationale"],
                "human_use_basis": verdict["human_use_basis"],
                "other_use_basis": _truncate(verdict["other_use_basis"], 400),
                "v12_decision": prior.get("decision", ""),
                "v12_decision_class": prior.get("decision_class", ""),
                # documented status
                "clinical_phase": phase,
                "clinical_phase_source": phase_source,
                "phase_conflict": phase_conflict,
                "pubchem_max_phase": ev.get("pubchem_max_phase") or "",
                "broad_phase": ev.get("broad_phase") or "",
                "broad_match": ev.get("broad_match") or "",
                "withdrawn": flags["withdrawn"],
                # mechanism / targets
                "moa": moa,
                "moa_sources": _join(moa_sources),
                "moa_frozen": key["moa_frozen"],
                "targets": _join(targets[:25]),
                "target_sources": _join(target_sources),
                "drug_classes": _join(ev.get("drug_classes") or []),
                # ATC / MeSH annotation under audit
                "atc_codes": _join(frozen_complete + frozen_partial),
                "atc_frozen_complete": _join(frozen_complete),
                "atc_frozen_partial": _join(frozen_partial),
                "atc_pubchem": _join(pubchem_codes),
                "atc_human_codes": _join(ev.get("atc_human_codes") or []),
                "mesh_frozen": key["mesh_frozen"],
                "mesh_corroborated": _mesh_corroborated(key["mesh_frozen"], haystack),
                # use evidence
                "use_classification": _join(ev.get("use_classification") or []),
                "industry_uses": _join(ev.get("industry_uses") or []),
                "uses_texts": _truncate(_join(ev.get("uses_texts") or [], " | "), 800),
                # indication text: recorded, not human-use evidence
                "drug_indications": _truncate(_join(ev.get("drug_indications") or [], " | "), 500),
                "therapeutic_uses": _truncate(
                    _join(ev.get("therapeutic_uses") or [], " | "), 500
                ),
                # flags and provenance
                "flag_mechanism_evidence": flags["mechanism_evidence"],
                "flag_tool_use": flags["tool_use"],
                "flag_human_use_marker": flags["human_use_marker"],
                "flag_atc_human": flags["atc_human"],
                "flag_human_indication": flags["human_indication"],
                "flag_other_use": flags["other_use"],
                "status_without_human_use_marker": (
                    verdict["decision_class"] in {"approved", "clinical"}
                    and not flags["human_use_marker"]
                ),
                "identity_resolved": flags["identity_resolved"],
                "identity_notes": ev.get("identity_notes") or "",
                "name_in_pubchem_synonyms": ev.get("name_in_pubchem_synonyms"),
                "smiles_key": ev.get("smiles_key") or "",
                "bindingdb_n_pairs": ev.get("bindingdb_n_pairs", 0),
                "bindingdb_targets": _join((ev.get("bindingdb_targets") or [])[:10]),
                "old_tier_atc_moa": key["old_tier_atc_moa"],
                "old_tier_atc_moa_mesh": key["old_tier_atc_moa_mesh"],
                "old_tier_all_smiles": key["old_tier_all_smiles"],
                "evidence_urls": _join([u["url"] for u in ev.get("evidence_urls", [])], " | "),
                "retrieved_at": ev.get("retrieved_at") or "",
                "evidence_ok_urls": ev.get("evidence_ok_urls", 0),
            }
        )
    return rows


def _movement(previous_decision, decision):
    return f"{previous_decision}->{decision}"


def summarize(rows, keys, sources, predecessor):
    total = len(rows)
    include = [r for r in rows if r["decision"] == "include"]
    exclude = [r for r in rows if r["decision"] == "exclude"]
    uncertain = [r for r in rows if r["decision"] == "uncertain"]

    classes = Counter(r["decision_class"] for r in rows)
    by_decision = Counter(r["decision"] for r in rows)
    class_by_decision = Counter((r["decision"], r["decision_class"]) for r in rows)
    precedence_counts = Counter(r["decision_precedence"] for r in rows)

    frozen_complete, frozen_partial = set(), set()
    for row in rows:
        frozen_complete.update(filter(None, row["atc_frozen_complete"].split("; ")))
        frozen_partial.update(filter(None, row["atc_frozen_partial"].split("; ")))
    frozen_levels = rules.atc_level_counts(frozen_complete, frozen_partial)

    pubchem_complete, pubchem_partial = set(), set()
    for row in rows:
        for code in filter(None, row["atc_pubchem"].split("; ")):
            (pubchem_complete if len(code) == 7 else pubchem_partial).add(code)
    pubchem_levels = rules.atc_level_counts(pubchem_complete, pubchem_partial)

    atc_keys_frozen = sum(1 for r in rows if r["atc_frozen_complete"] or r["atc_frozen_partial"])
    atc_keys_pubchem = sum(1 for r in rows if r["atc_pubchem"])
    atc_agree = sum(
        1
        for r in rows
        if r["atc_frozen_complete"]
        and set(r["atc_frozen_complete"].split("; ")) & set(r["atc_pubchem"].split("; "))
    )

    frozen_moa = [key["moa_frozen"] for key in keys if key["moa_frozen"]]
    frozen_fragments = set()
    for value in frozen_moa:
        for fragment in value.split(";"):
            fragment = fragment.strip()
            if fragment:
                frozen_fragments.add(fragment)

    old = {}
    for tier, column in (
        ("atc-moa", "old_tier_atc_moa"),
        ("atc-moa-mesh", "old_tier_atc_moa_mesh"),
        ("all-smiles", "old_tier_all_smiles"),
    ):
        old_keys = {r["drug_key"] for r in rows if r[column]}
        old.update(
            {
                tier: {
                    "old_n": len(old_keys),
                    "old_include_new_include": sum(
                        1 for r in rows if r[column] and r["decision"] == "include"
                    ),
                    "old_include_new_exclude": sum(
                        1 for r in rows if r[column] and r["decision"] == "exclude"
                    ),
                    "old_include_new_uncertain": sum(
                        1 for r in rows if r[column] and r["decision"] == "uncertain"
                    ),
                    "false_negative": sum(
                        1 for r in rows if not r[column] and r["decision"] == "include"
                    ),
                    "false_positive": sum(
                        1 for r in rows if r[column] and r["decision"] == "exclude"
                    ),
                    "false_negative_keys": [
                        r["drug_key"] for r in rows if not r[column] and r["decision"] == "include"
                    ],
                    "false_positive_keys": [
                        r["drug_key"] for r in rows if r[column] and r["decision"] == "exclude"
                    ],
                }
            }
        )

    moved = [
        r
        for r in rows
        if (r["v12_decision"], r["v12_decision_class"])
        != (r["decision"], r["decision_class"])
    ]
    movement_counts = Counter(
        _movement(r["v12_decision"], r["decision"]) for r in moved
    )
    by_movement = {
        key: movement_counts.get(key, 0)
        for key in (
            "include->exclude",
            "exclude->include",
            "uncertain->exclude",
            "uncertain->include",
            "include->uncertain",
            "exclude->uncertain",
        )
    }
    retrieved = sorted({r["retrieved_at"] for r in rows if r["retrieved_at"]})
    generated_at = retrieved[-1] if retrieved else evidence.now_iso()

    return {
        "schema": SCHEMA_SUMMARY,
        "artifact_version": ARTIFACT_VERSION,
        "rule": {"version": rules.RULE_VERSION, "sha256": rules.rule_sha256()},
        "predecessor": predecessor,
        # Deterministic evidence cutoff (latest retrieval in the record), so an
        # offline rebuild reproduces every digest.
        "generated_at": generated_at,
        "source": sources,
        "counts": {
            "total_keys": total,
            "total_name_records": sum(int(r["n_name_records"]) for r in rows),
            "include": len(include),
            "exclude": len(exclude),
            "uncertain": len(uncertain),
        },
        "decision_class_counts": dict(sorted(classes.items())),
        "decision_counts": dict(sorted(by_decision.items())),
        "decision_by_class": {
            f"{decision}::{klass}": count
            for (decision, klass), count in sorted(class_by_decision.items())
        },
        "precedence_counts": dict(sorted(precedence_counts.items())),
        "human_use": {
            "keys_with_human_use_marker": sum(
                1 for r in rows if r["flag_human_use_marker"]
            ),
            "keys_with_human_atc": sum(1 for r in rows if r["flag_atc_human"]),
            "keys_with_human_indication": sum(
                1 for r in rows if r["flag_human_indication"]
            ),
            "status_without_human_use_marker": sum(
                1 for r in rows if r["status_without_human_use_marker"]
            ),
        },
        "status": {
            "phase_conflicts": sum(1 for r in rows if r["phase_conflict"]),
            "approved_by_pubchem": sum(
                1 for r in rows
                if str(r["pubchem_max_phase"]).strip().lower() in rules.APPROVED_PHASES
            ),
            "approved_by_broad": sum(
                1 for r in rows
                if str(r["broad_phase"]).strip().lower() in rules.APPROVED_PHASES
            ),
            "clinical": sum(1 for r in rows if r["decision_class"] == "clinical"),
            "withdrawn": sum(1 for r in rows if r["withdrawn"]),
            "status_without_human_use_marker": sum(
                1 for r in rows if r["status_without_human_use_marker"]
            ),
            "excluded_with_documented_status": [
                r["drug_key"]
                for r in rows
                if r["decision"] == "exclude"
                and str(r["clinical_phase"]).strip().lower()
                in rules.APPROVED_PHASES | rules.CLINICAL_PHASES
            ],
        },
        "identity": {
            "resolved": sum(1 for r in rows if r["identity_resolved"]),
            "unresolved": [
                r["drug_key"] for r in rows if not r["identity_resolved"]
            ],
            "pubchem_cid_resolved": sum(1 for r in rows if r["pubchem_cid"] != ""),
            "name_in_synonyms_true": sum(
                1 for r in rows if r["name_in_pubchem_synonyms"] is True
            ),
            "name_in_synonyms_false": sum(
                1 for r in rows if r["name_in_pubchem_synonyms"] is False
            ),
            "name_unchecked": sum(
                1 for r in rows if r["name_in_pubchem_synonyms"] is None
            ),
            "salt_or_parent_of_drug_key": sum(
                1 for r in rows if "salt/parent" in r["identity_notes"]
            ),
            "smiles_key_mismatch": sum(
                1 for r in rows if "frozen smiles InChIKey differs" in r["identity_notes"]
            ),
        },
        "atc": {
            "frozen": {
                "keys_with_atc": atc_keys_frozen,
                "complete_codes": len(frozen_complete),
                "partial_codes": len(frozen_partial),
                "levels_complete_only": frozen_levels["complete"],
                "levels_complete_plus_partial": frozen_levels["complete_plus_partial"],
            },
            "pubchem_verified": {
                "keys_with_atc": atc_keys_pubchem,
                "complete_codes": len(pubchem_complete),
                "partial_codes": len(pubchem_partial),
                "levels_complete_only": pubchem_levels["complete"],
                "levels_complete_plus_partial": pubchem_levels["complete_plus_partial"],
            },
            "keys_where_frozen_and_pubchem_share_a_complete_code": atc_agree,
        },
        "moa": {
            "keys_with_frozen_moa": len(frozen_moa),
            "distinct_frozen_moa_strings": len(set(frozen_moa)),
            "distinct_frozen_moa_fragments": len(frozen_fragments),
            "keys_with_documented_mechanism": sum(
                1 for r in rows if r["flag_mechanism_evidence"]
            ),
            "keys_with_broad_moa": sum(
                1 for r in rows if "broad_repurposing_hub" in r["moa_sources"]
            ),
            "keys_with_pubchem_mechanism": sum(
                1 for r in rows if "pubchem_mechanism_of_action" in r["moa_sources"]
            ),
        },
        "mesh": {
            "keys_with_frozen_mesh": sum(1 for r in rows if r["mesh_frozen"]),
            "keys_corroborated_by_pubchem": sum(
                1 for r in rows if r["mesh_corroborated"].startswith(("1/", "2/", "3/"))
            ),
        },
        "old_tier": old,
        "change_report": {
            "schema": SCHEMA_CHANGE_REPORT,
            "moved": len(moved),
            "class_only_changes": sum(
                1 for r in moved if r["v12_decision"] == r["decision"]
            ),
            "by_movement": by_movement,
            "by_previous_decision": dict(
                sorted(Counter(r["v12_decision"] for r in moved).items())
            ),
            "moved_keys": [r["drug_key"] for r in moved],
        },
        "sources": sources.get("evidence_sources", {}),
        "uncertain_keys": [r["drug_key"] for r in uncertain],
        "excluded_keys": [r["drug_key"] for r in exclude],
    }


def write_table(rows, path):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TABLE_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in TABLE_COLUMNS})


def write_old_tier_comparison(rows, path):
    fieldnames = [
        "drug_key",
        "compound_name",
        "old_tier_atc_moa",
        "old_tier_atc_moa_mesh",
        "new_decision",
        "new_decision_class",
        "comparison_atc_moa",
        "rationale",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            if row["old_tier_atc_moa"] and row["decision"] == "include":
                comparison = "agrees_include"
            elif row["old_tier_atc_moa"] and row["decision"] == "exclude":
                comparison = "false_positive"
            elif row["old_tier_atc_moa"] and row["decision"] == "uncertain":
                comparison = "old_include_now_uncertain"
            elif not row["old_tier_atc_moa"] and row["decision"] == "include":
                comparison = "false_negative"
            elif not row["old_tier_atc_moa"] and row["decision"] == "exclude":
                comparison = "agrees_exclude"
            else:
                comparison = "old_exclude_now_uncertain"
            writer.writerow(
                {
                    "drug_key": row["drug_key"],
                    "compound_name": row["compound_name"],
                    "old_tier_atc_moa": row["old_tier_atc_moa"],
                    "old_tier_atc_moa_mesh": row["old_tier_atc_moa_mesh"],
                    "new_decision": row["decision"],
                    "new_decision_class": row["decision_class"],
                    "comparison_atc_moa": comparison,
                    "rationale": row["rationale"],
                }
            )


CHANGE_COLUMNS = [
    "drug_key",
    "compound_name",
    "v12_decision",
    "v12_decision_class",
    "v13_decision",
    "v13_decision_class",
    "decision_precedence",
    "deciding_evidence",
    "human_use_basis",
    "other_use_basis",
    "rationale_v13",
    "evidence_urls",
]


def write_change_report(rows, path):
    """One row per decision that moved from v1.2, with the evidence that decided it."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CHANGE_COLUMNS)
        writer.writeheader()
        for row in rows:
            if (row["v12_decision"], row["v12_decision_class"]) == (
                row["decision"],
                row["decision_class"],
            ):
                continue
            deciding = (
                row["human_use_basis"]
                or row["other_use_basis"]
                or row["rationale"]
            )
            writer.writerow(
                {
                    "drug_key": row["drug_key"],
                    "compound_name": row["compound_name"],
                    "v12_decision": row["v12_decision"],
                    "v12_decision_class": row["v12_decision_class"],
                    "v13_decision": row["decision"],
                    "v13_decision_class": row["decision_class"],
                    "decision_precedence": row["decision_precedence"],
                    "deciding_evidence": deciding,
                    "human_use_basis": row["human_use_basis"],
                    "other_use_basis": row["other_use_basis"],
                    "rationale_v13": row["rationale"],
                    "evidence_urls": row["evidence_urls"],
                }
            )


def _handoff_entry(row):
    return {
        "drug_key": row["drug_key"],
        "decision": row["decision"],
        "decision_class": row["decision_class"],
        "decision_precedence": row["decision_precedence"],
        "compound_name": row["compound_name"],
        "all_names": row["all_names"],
        "inchikey": row["inchikey"],
        "pubchem_cid": row["pubchem_cid"],
        "clinical_phase": row["clinical_phase"],
        "moa": row["moa"],
        "targets": row["targets"],
        "atc_codes": row["atc_codes"],
        "human_use_basis": row["human_use_basis"],
        "other_use_basis": row["other_use_basis"],
        "evidence_urls": row["evidence_urls"],
        "retrieved_at": row["retrieved_at"],
        "rationale": row["rationale"],
    }


def write_handoff(rows, summary, path):
    payload = {
        "schema": SCHEMA_HANDOFF,
        "artifact_version": ARTIFACT_VERSION,
        "rule": {"version": rules.RULE_VERSION, "sha256": rules.rule_sha256()},
        "predecessor": summary["predecessor"],
        "generated_at": summary["generated_at"],
        "counts": summary["counts"],
        "decision_class_counts": summary["decision_class_counts"],
        "precedence_counts": summary["precedence_counts"],
        "include": [_handoff_entry(r) for r in rows if r["decision"] == "include"],
        "uncertain": [_handoff_entry(r) for r in rows if r["decision"] == "uncertain"],
        "exclude": [_handoff_entry(r) for r in rows if r["decision"] == "exclude"],
    }
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def write_evidence_manifest(rows, path):
    fieldnames = [
        "drug_key",
        "decision",
        "decision_class",
        "decision_precedence",
        "n_urls",
        "n_ok_urls",
        "sources",
        "retrieved_at",
        "identity_resolved",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            urls = [url for url in row["evidence_urls"].split(" | ") if url]
            sources = []
            for url in urls:
                if "pubchem" in url:
                    sources.append("pubchem")
                elif "broadinstitute" in url:
                    sources.append("broad")
                elif "bindingdb" in url:
                    sources.append("bindingdb")
            writer.writerow(
                {
                    "drug_key": row["drug_key"],
                    "decision": row["decision"],
                    "decision_class": row["decision_class"],
                    "decision_precedence": row["decision_precedence"],
                    "n_urls": len(urls),
                    "n_ok_urls": row["evidence_ok_urls"],
                    "sources": "; ".join(sorted(set(sources))),
                    "retrieved_at": row["retrieved_at"],
                    "identity_resolved": row["identity_resolved"],
                }
            )


def _read_table(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return {row["drug_key"]: row for row in csv.DictReader(handle)}


def _read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _write_json(payload, path):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def build(out_dir, compounds_csv, allow_network=True, bindingdb=DEFAULT_BINDINGDB,
          jobs=3):
    out_dir = Path(out_dir)
    version_dir = out_dir / ARTIFACT_VERSION
    previous_dir = out_dir / PREVIOUS_VERSION
    version_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = out_dir / "cache"

    previous_table = previous_dir / "eligibility_table.csv"
    previous_manifest = previous_dir / "manifest.json"
    if not previous_table.exists() or not previous_manifest.exists():
        raise SystemExit(
            f"frozen v1.2 record missing under {previous_dir}: the v1.3 build "
            "carries the v1.2 decisions as its change-report baseline"
        )
    previous = _read_table(previous_table)
    previous_manifest_payload = _read_json(previous_manifest)

    compounds = evidence.ROOT / compounds_csv
    if not compounds.exists():
        raise SystemExit(f"frozen compounds table not found: {compounds}")
    keys = evidence.frozen_keys(compounds)
    if allow_network:
        evidence.prefetch(keys, cache_dir, allow_network=True, jobs=jobs)
    key_rows = evidence.build_evidence(
        keys, cache_dir, allow_network=allow_network, bindingdb_path=bindingdb
    )
    evidence_map, binding_info = key_rows
    rows = assemble(keys, evidence_map, previous)

    chembl = evidence.probe_chembl(cache_dir, allow_network=allow_network)
    binding_available = bool(binding_info.get("bindingdb", {}).get("available"))
    usable_pairs = sum(
        1
        for row in rows
        if row["target_sources"] and "bindingdb" in row["target_sources"]
    )
    sources = {
        "compounds_csv": str(compounds_csv),
        "compounds_csv_sha256": evidence.sha256_file(compounds),
        "evidence_sources": {
            "pubchem": "https://pubchem.ncbi.nlm.nih.gov/ (PUG-REST + PUG-View)",
            "broad_repurposing_hub": "repo-hub.broadinstitute.org (local Biomni parquet)",
            "bindingdb": {
                "path": binding_info.get("bindingdb", {}).get("path"),
                "available": binding_available,
                "keys_with_target_evidence": usable_pairs,
                "join_note": (
                    "attempted; the InChIKey join yielded no per-key ligand-target "
                    "pairs in this build, so BindingDB is recorded but not used as "
                    "mechanism evidence (deviation D5)"
                ),
            },
            "chembl_probe": chembl,
            "frozen_annotation": "compounds.csv atc/moa/mesh_classes (verified, not evidence)",
        },
        "cache_dir": str(cache_dir),
    }
    predecessor = {
        "artifact_version": PREVIOUS_VERSION,
        "rule": previous_manifest_payload["rule"],
        "files": previous_manifest_payload["files"],
    }
    summary = summarize(rows, keys, sources, predecessor)

    table_path = version_dir / "eligibility_table.csv"
    old_path = version_dir / "old_tier_comparison.csv"
    evidence_path = version_dir / "evidence_manifest.csv"
    summary_path = version_dir / "summary.json"
    handoff_path = version_dir / "handoff.json"
    change_path = version_dir / "change_report.csv"
    write_table(rows, table_path)
    write_old_tier_comparison(rows, old_path)
    write_evidence_manifest(rows, evidence_path)
    write_handoff(rows, summary, handoff_path)
    write_change_report(rows, change_path)
    _write_json(summary, summary_path)

    manifest = {
        "schema": SCHEMA_MANIFEST,
        "artifact_version": ARTIFACT_VERSION,
        "rule": summary["rule"],
        "predecessor": predecessor,
        "counts": summary["counts"],
        "generated_at": summary["generated_at"],
        "files": [],
    }
    for path in (table_path, old_path, evidence_path, summary_path, handoff_path,
                 change_path):
        manifest["files"].append(
            {
                "path": path.name,
                "sha256": evidence.sha256_file(path),
                "bytes": path.stat().st_size,
            }
        )
    _write_json(manifest, version_dir / "manifest.json")

    current = {
        "schema": SCHEMA_CURRENT,
        "artifact_version": ARTIFACT_VERSION,
        "rule": summary["rule"],
        "handoff": {
            "path": f"{ARTIFACT_VERSION}/handoff.json",
            "schema": SCHEMA_HANDOFF,
            "sha256": evidence.sha256_file(handoff_path),
        },
        "table": {
            "path": f"{ARTIFACT_VERSION}/eligibility_table.csv",
            "sha256": evidence.sha256_file(table_path),
        },
        "manifest": {
            "path": f"{ARTIFACT_VERSION}/manifest.json",
            "sha256": evidence.sha256_file(version_dir / "manifest.json"),
        },
        "previous": {
            "artifact_version": PREVIOUS_VERSION,
            "manifest": f"{PREVIOUS_VERSION}/manifest.json",
            "rule": predecessor["rule"],
        },
    }
    _write_json(current, out_dir / "current.json")
    return rows, summary, manifest


def check_artifacts(out_dir, compounds_csv=COMPOUNDS_CSV):
    """Re-run the audit gates on the written record; returns (ok, lines)."""
    out_dir = Path(out_dir)
    version_dir = out_dir / ARTIFACT_VERSION
    previous_dir = out_dir / PREVIOUS_VERSION
    lines, ok = [], True

    def gate(name, passed, detail):
        nonlocal ok
        ok = ok and passed
        lines.append(f"{'PASS' if passed else 'FAIL'}  {name}: {detail}")

    table_path = version_dir / "eligibility_table.csv"
    if not table_path.exists():
        return False, [f"FAIL  table: {table_path} missing"]
    rows = list(csv.DictReader(open(table_path, newline="", encoding="utf-8")))
    keys = evidence.frozen_keys(evidence.ROOT / compounds_csv)
    expected = {key["drug_key"] for key in keys}

    # -- row completeness -------------------------------------------------
    gate("row_count", len(rows) == 460, f"{len(rows)} rows (expected 460)")
    gate("key_uniqueness", len({r["drug_key"] for r in rows}) == len(rows), "drug_key unique")
    gate(
        "key_set_matches_frozen",
        {r["drug_key"] for r in rows} == expected,
        f"{len(expected)} frozen keys, table matches",
    )
    gate(
        "name_records",
        sum(int(r["n_name_records"]) for r in rows) == 467,
        f"{sum(int(r['n_name_records']) for r in rows)} name records (expected 467)",
    )
    valid_decisions = {"include", "exclude", "uncertain"}
    valid_classes = {"approved", "clinical", "mechanism_tool", "other_use", "unknown"}
    gate(
        "decision_domain",
        all(r["decision"] in valid_decisions for r in rows)
        and all(r["decision_class"] in valid_classes for r in rows),
        "decision and decision_class in domain",
    )

    # -- rule v1.3 identity and precedence --------------------------------
    rule_ok = True
    for artifact in ("summary.json", "handoff.json", "manifest.json",
                     "current.json"):
        path = out_dir / artifact if artifact == "current.json" else version_dir / artifact
        if not path.exists():
            rule_ok = False
            continue
        payload = _read_json(path)
        rule_ok = rule_ok and payload.get("rule") == {
            "version": rules.RULE_VERSION,
            "sha256": rules.rule_sha256(),
        }
    gate(
        "rule_identity_v13",
        rule_ok,
        f"every artifact carries {rules.RULE_VERSION} sha256 {rules.rule_sha256()[:12]}",
    )
    expected_precedence = {
        "include": {"human_drug_use", "mechanism_tool"},
        "exclude": {"environmental_or_industrial_use"},
        "uncertain": {"identity_unresolved", "insufficient_evidence"},
    }
    bad_precedence = [
        r["drug_key"]
        for r in rows
        if r["decision_precedence"] not in expected_precedence.get(r["decision"], set())
    ]
    gate(
        "precedence_recorded",
        not bad_precedence,
        f"{len(bad_precedence)} rows with a precedence label inconsistent with the decision",
    )
    no_basis = [
        r["drug_key"] for r in rows if r["decision"] == "include" and not r["human_use_basis"]
    ]
    gate(
        "include_records_human_use",
        not no_basis,
        f"{len(no_basis)} includes without a recorded human-use basis",
    )
    no_use_basis = [
        r["drug_key"] for r in rows if r["decision"] == "exclude" and not r["other_use_basis"]
    ]
    gate(
        "exclude_records_other_use",
        not no_use_basis,
        f"{len(no_use_basis)} excludes without a recorded other-use basis",
    )
    mechanism_only = [
        r["drug_key"]
        for r in rows
        if r["decision"] == "include"
        and r["decision_class"] == "mechanism_tool"
        and "use:" not in r["human_use_basis"]
    ]
    gate(
        "mechanism_alone_never_includes",
        not mechanism_only,
        f"{len(mechanism_only)} mechanism tools without a use/tool marker",
    )

    # -- existing traceability gates -------------------------------------
    missing = [
        r["drug_key"]
        for r in rows
        if r["decision"] in {"include", "exclude"}
        and (not r["rationale"].strip() or int(r["evidence_ok_urls"]) < 1)
    ]
    gate("traceable_decisions", not missing, f"{len(missing)} include/exclude rows untraceable")
    unknown_excluded = [
        r["drug_key"] for r in rows if r["decision"] == "exclude" and r["decision_class"] == "unknown"
    ]
    gate("unknown_never_excluded", not unknown_excluded, f"{len(unknown_excluded)} unknown excludes")
    other_use_excluded = [
        r["drug_key"]
        for r in rows
        if r["decision"] == "exclude" and r["flag_other_use"] != "True"
    ]
    gate("exclude_requires_other_use", not other_use_excluded, f"{len(other_use_excluded)} bad excludes")
    identity_gate = [
        r["drug_key"]
        for r in rows
        if r["identity_resolved"] != "True" and r["decision"] != "uncertain"
    ]
    gate("unresolved_identity_uncertain", not identity_gate, f"{len(identity_gate)} violations")
    gate(
        "old_tier_flags",
        all(r["old_tier_atc_moa"] in {"True", "False"} for r in rows)
        and all(r["old_tier_atc_moa_mesh"] in {"True", "False"} for r in rows)
        and all(r["old_tier_all_smiles"] in {"True", "False"} for r in rows),
        "boolean old-tier flags on every row",
    )
    gate(
        "retrieval_timestamps",
        all(r["retrieved_at"].strip() for r in rows),
        "every row has a retrieval timestamp",
    )

    # -- v1.2 record preserved (ADR 0030: no silent overwrite) -----------
    previous_ok = True
    previous_detail = "v1.2 manifest hashes recomputed on the frozen record"
    previous_manifest_path = previous_dir / "manifest.json"
    if not previous_manifest_path.exists():
        previous_ok = False
        previous_detail = f"{previous_manifest_path} missing"
    else:
        previous_manifest = _read_json(previous_manifest_path)
        previous_ok = previous_manifest["rule"]["version"] == "eligibility-rule-v1.2"
        for item in previous_manifest["files"]:
            path = previous_dir / item["path"]
            if not path.exists() or evidence.sha256_file(path) != item["sha256"]:
                previous_ok = False
                previous_detail = f"v1.2 file changed or missing: {item['path']}"
                break
    gate("v12_preserved", previous_ok, previous_detail)

    # -- change report completeness --------------------------------------
    change_path = version_dir / "change_report.csv"
    change_ok, change_detail = False, f"{change_path} missing"
    previous_by_key = _read_table(previous_dir / "eligibility_table.csv") \
        if (previous_dir / "eligibility_table.csv").exists() else {}
    if change_path.exists():
        change = list(csv.DictReader(open(change_path, newline="", encoding="utf-8")))
        reported = {r["drug_key"] for r in change}
        moved = {
            r["drug_key"]
            for r in rows
            if (r["v12_decision"], r["v12_decision_class"])
            != (r["decision"], r["decision_class"])
        }
        change_ok = reported == moved
        change_detail = f"{len(reported)} reported / {len(moved)} moved rows"
        by_key = {r["drug_key"]: r for r in rows}
        mismatched = [
            r["drug_key"]
            for r in change
            if not r["deciding_evidence"].strip()
            or not r["rationale_v13"].strip()
            or r["drug_key"] not in by_key
            or r["v13_decision"] != by_key[r["drug_key"]]["decision"]
            or r["v13_decision_class"] != by_key[r["drug_key"]]["decision_class"]
        ]
        if mismatched:
            change_ok = False
            change_detail += f"; {len(mismatched)} rows without deciding evidence"
        stale = [
            r["drug_key"]
            for r in change
            if r["drug_key"] in previous_by_key
            and r["v12_decision"] != previous_by_key[r["drug_key"]]["decision"]
        ]
        if stale:
            change_ok = False
            change_detail += f"; {len(stale)} rows disagree with the frozen v1.2 table"
    gate("change_report_complete", change_ok, change_detail)

    # -- manifest hashes and the current pointer -------------------------
    manifest_path = version_dir / "manifest.json"
    handoff_path = version_dir / "handoff.json"
    manifest_ok, manifest_detail = False, f"{manifest_path} or {handoff_path} missing"
    if manifest_path.exists() and handoff_path.exists():
        manifest = _read_json(manifest_path)
        hashes = {item["path"]: item["sha256"] for item in manifest["files"]}
        manifest_ok = True
        for name in ("eligibility_table.csv", "handoff.json", "summary.json",
                     "old_tier_comparison.csv", "evidence_manifest.csv",
                     "change_report.csv"):
            if hashes.get(name) != evidence.sha256_file(version_dir / name):
                manifest_ok = False
                manifest_detail = f"{name} sha256 does not match manifest"
                break
        else:
            manifest_detail = "all tracked artifact sha256 match the manifest"
        predecessor = manifest.get("predecessor", {})
        previous_manifest = _read_json(previous_manifest_path) \
            if previous_manifest_path.exists() else {}
        if not (
            predecessor.get("artifact_version") == PREVIOUS_VERSION
            and predecessor.get("rule") == previous_manifest.get("rule")
            and predecessor.get("files") == previous_manifest.get("files")
        ):
            manifest_ok = False
            manifest_detail = "v1.3 manifest predecessor block does not match v1.2"
    gate("manifest_hashes", manifest_ok, manifest_detail)

    handoff = _read_json(handoff_path) if handoff_path.exists() else {}
    summary = _read_json(version_dir / "summary.json") \
        if (version_dir / "summary.json").exists() else {}
    gate(
        "handoff_schema_v2",
        handoff.get("schema") == SCHEMA_HANDOFF
        and handoff.get("counts") == summary.get("counts"),
        f"schema {SCHEMA_HANDOFF}, counts equal summary",
    )
    gate(
        "handoff_partition",
        len(handoff.get("include", []))
        + len(handoff.get("uncertain", []))
        + len(handoff.get("exclude", []))
        == len(rows),
        "include + uncertain + exclude = every row",
    )
    current_path = out_dir / "current.json"
    current_ok, current_detail = False, f"{current_path} missing"
    if current_path.exists():
        current = _read_json(current_path)
        current_ok = (
            current.get("artifact_version") == ARTIFACT_VERSION
            and current.get("rule") == summary.get("rule")
            and current.get("handoff", {}).get("path") == f"{ARTIFACT_VERSION}/handoff.json"
            and current.get("handoff", {}).get("sha256") == evidence.sha256_file(handoff_path)
            and current.get("table", {}).get("sha256") == evidence.sha256_file(table_path)
            and current.get("previous", {}).get("artifact_version") == PREVIOUS_VERSION
        )
        current_detail = "current.json points at the v1.3 record and matches its digests"
    gate("current_pointer", current_ok, current_detail)
    return ok, lines


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-root",
        default=str(evidence.STEP),
        help="step directory that owns output/ (default: this step)",
    )
    parser.add_argument("--compounds", default=str(COMPOUNDS_CSV))
    parser.add_argument("--bindingdb", default=str(DEFAULT_BINDINGDB))
    parser.add_argument("--offline", action="store_true", help="use cached evidence only")
    parser.add_argument("--jobs", type=int, default=3, help="parallel PubChem fetches")
    parser.add_argument("--check", action="store_true", help="verify the written record")
    args = parser.parse_args(argv)
    out_dir = Path(args.out_root) / "output"

    if args.check:
        ok, lines = check_artifacts(out_dir, Path(args.compounds))
        print("\n".join(lines))
        return 0 if ok else 1

    rows, summary, manifest = build(
        out_dir,
        Path(args.compounds),
        allow_network=not args.offline,
        bindingdb=Path(args.bindingdb),
        jobs=args.jobs,
    )
    print(
        f"wrote {len(rows)} rows: "
        f"include={summary['counts']['include']} "
        f"exclude={summary['counts']['exclude']} "
        f"uncertain={summary['counts']['uncertain']} "
        f"(rule {summary['rule']['version']} sha256 {summary['rule']['sha256'][:12]})"
    )
    print(
        f"v1.2 -> v1.3: {summary['change_report']['moved']} rows moved "
        f"{summary['change_report']['by_movement']}"
    )
    print(f"artifacts: {out_dir / ARTIFACT_VERSION}")
    ok, lines = check_artifacts(out_dir, Path(args.compounds))
    print("\n".join(lines))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
