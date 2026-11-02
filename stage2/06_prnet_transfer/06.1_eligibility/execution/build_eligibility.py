"""Build the 460-structure pharmacological eligibility audit for temp_env.

Reads the frozen ``compounds.csv`` (identity + annotation under audit), collects
documented evidence (PubChem, Broad Repurposing Hub, BindingDB), applies the
locked rule in ``rules.py``, and writes the evidence table, coverage summary, the
machine-readable handoff for ticket #3, and the hashed artifact manifest.

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
from datetime import datetime, timezone
from pathlib import Path

import evidence
import rules

SCHEMA_TABLE = "prnet-eligibility-table-v1"
SCHEMA_SUMMARY = "prnet-eligibility-summary-v1"
SCHEMA_HANDOFF = "prnet-eligibility-handoff-v1"
COMPOUNDS_CSV = Path("stage2/02_cm_data/output/datasets-1/compounds.csv")
DEFAULT_BINDINGDB = evidence.DEFAULT_BINDINGDB

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
    "rationale_category",
    "rationale",
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
    # flags and provenance
    "flag_mechanism_evidence",
    "flag_tool_use",
    "flag_human_pharm_use",
    "flag_other_use",
    "flag_pharmaceutical_use",
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


def assemble(keys, evidence_map):
    """Apply the locked rule and flatten each key into one table row."""
    rows = []
    for key in keys:
        ev = evidence_map[key["drug_key"]]
        verdict = rules.decide(ev)
        flags = verdict["flags"]

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
                "rationale_category": verdict["rationale_category"],
                "rationale": verdict["rationale"],
                "clinical_phase": phase,
                "clinical_phase_source": phase_source,
                "phase_conflict": phase_conflict,
                "pubchem_max_phase": ev.get("pubchem_max_phase") or "",
                "broad_phase": ev.get("broad_phase") or "",
                "broad_match": ev.get("broad_match") or "",
                "withdrawn": flags["withdrawn"],
                "moa": moa,
                "moa_sources": _join(moa_sources),
                "moa_frozen": key["moa_frozen"],
                "targets": _join(targets[:25]),
                "target_sources": _join(target_sources),
                "drug_classes": _join(ev.get("drug_classes") or []),
                "atc_codes": _join(frozen_complete + frozen_partial),
                "atc_frozen_complete": _join(frozen_complete),
                "atc_frozen_partial": _join(frozen_partial),
                "atc_pubchem": _join(pubchem_codes),
                "atc_human_codes": _join(ev.get("atc_human_codes") or []),
                "mesh_frozen": key["mesh_frozen"],
                "mesh_corroborated": _mesh_corroborated(key["mesh_frozen"], haystack),
                "use_classification": _join(ev.get("use_classification") or []),
                "industry_uses": _join(ev.get("industry_uses") or []),
                "flag_mechanism_evidence": flags["mechanism_evidence"],
                "flag_tool_use": flags["tool_use"],
                "flag_human_pharm_use": flags["human_pharm_use"],
                "flag_other_use": flags["other_use"],
                "flag_pharmaceutical_use": flags["pharmaceutical_use"],
                "status_without_human_use_marker": (
                    verdict["decision_class"] in {"approved", "clinical"}
                    and not flags["human_pharm_use"]
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


def summarize(rows, keys, sources):
    total = len(rows)
    include = [r for r in rows if r["decision"] == "include"]
    exclude = [r for r in rows if r["decision"] == "exclude"]
    uncertain = [r for r in rows if r["decision"] == "uncertain"]

    classes = Counter(r["decision_class"] for r in rows)
    by_decision = Counter(r["decision"] for r in rows)
    class_by_decision = Counter((r["decision"], r["decision_class"]) for r in rows)

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

    return {
        "schema": SCHEMA_SUMMARY,
        "rule": {"version": rules.RULE_VERSION, "sha256": rules.rule_sha256()},
        "generated_at": evidence.now_iso(),
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
            "withdrawn": sum(1 for r in rows if r["withdrawn"] == "True"),
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
            "unresolved": [r["drug_key"] for r in rows if not r["identity_resolved"]],
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
            "keys_with_broad_moa": sum(1 for r in rows if "broad_repurposing_hub" in r["moa_sources"]),
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


def write_handoff(rows, summary, path):
    def entry(row):
        return {
            "drug_key": row["drug_key"],
            "decision": row["decision"],
            "decision_class": row["decision_class"],
            "compound_name": row["compound_name"],
            "all_names": row["all_names"],
            "inchikey": row["inchikey"],
            "pubchem_cid": row["pubchem_cid"],
            "clinical_phase": row["clinical_phase"],
            "moa": row["moa"],
            "targets": row["targets"],
            "atc_codes": row["atc_codes"],
            "evidence_urls": row["evidence_urls"],
            "retrieved_at": row["retrieved_at"],
            "rationale": row["rationale"],
        }

    payload = {
        "schema": SCHEMA_HANDOFF,
        "rule": {"version": rules.RULE_VERSION, "sha256": rules.rule_sha256()},
        "generated_at": summary["generated_at"],
        "counts": summary["counts"],
        "decision_class_counts": summary["decision_class_counts"],
        "include": [entry(r) for r in rows if r["decision"] == "include"],
        "uncertain": [entry(r) for r in rows if r["decision"] == "uncertain"],
        "exclude": [entry(r) for r in rows if r["decision"] == "exclude"],
    }
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def write_evidence_manifest(rows, path):
    fieldnames = [
        "drug_key",
        "decision",
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
                    "n_urls": len(urls),
                    "n_ok_urls": row["evidence_ok_urls"],
                    "sources": "; ".join(sorted(set(sources))),
                    "retrieved_at": row["retrieved_at"],
                    "identity_resolved": row["identity_resolved"],
                }
            )


def build(out_dir, compounds_csv, allow_network=True, bindingdb=DEFAULT_BINDINGDB,
          jobs=3):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = out_dir / "cache"

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
    rows = assemble(keys, evidence_map)

    chembl = evidence.probe_chembl(cache_dir, allow_network=allow_network)
    sources = {
        "compounds_csv": str(compounds_csv),
        "compounds_csv_sha256": evidence.sha256_file(compounds),
        "evidence_sources": {
            "pubchem": "https://pubchem.ncbi.nlm.nih.gov/ (PUG-REST + PUG-View)",
            "broad_repurposing_hub": "repo-hub.broadinstitute.org (local Biomni parquet)",
            "bindingdb": {
                "path": binding_info.get("path"),
                "available": bool(binding_info.get("available")),
            },
            "chembl_probe": chembl,
            "frozen_annotation": "compounds.csv atc/moa/mesh_classes (verified, not evidence)",
        },
        "cache_dir": str(cache_dir),
    }
    summary = summarize(rows, keys, sources)

    table_path = out_dir / "eligibility_table.csv"
    old_path = out_dir / "old_tier_comparison.csv"
    evidence_path = out_dir / "evidence_manifest.csv"
    summary_path = out_dir / "summary.json"
    handoff_path = out_dir / "handoff.json"
    write_table(rows, table_path)
    write_old_tier_comparison(rows, old_path)
    write_evidence_manifest(rows, evidence_path)
    write_handoff(rows, summary, handoff_path)
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, ensure_ascii=False)
        handle.write("\n")

    manifest = {
        "schema": "prnet-eligibility-manifest-v1",
        "rule": summary["rule"],
        "counts": summary["counts"],
        "generated_at": summary["generated_at"],
        "files": [],
    }
    for path in (table_path, old_path, evidence_path, summary_path, handoff_path):
        manifest["files"].append(
            {
                "path": path.name,
                "sha256": evidence.sha256_file(path),
                "bytes": path.stat().st_size,
            }
        )
    with open(out_dir / "manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    return rows, summary, manifest


def check_artifacts(out_dir, compounds_csv=COMPOUNDS_CSV):
    """Re-run the audit gates on the written record; returns (ok, lines)."""
    out_dir = Path(out_dir)
    lines, ok = [], True

    def gate(name, passed, detail):
        nonlocal ok
        ok = ok and passed
        lines.append(f"{'PASS' if passed else 'FAIL'}  {name}: {detail}")

    table_path = out_dir / "eligibility_table.csv"
    if not table_path.exists():
        return False, [f"FAIL  table: {table_path} missing"]
    with open(table_path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    keys = evidence.frozen_keys(evidence.ROOT / compounds_csv)
    expected = {key["drug_key"] for key in keys}

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

    manifest_path = out_dir / "manifest.json"
    handoff_path = out_dir / "handoff.json"
    if manifest_path.exists() and handoff_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        hashes = {item["path"]: item["sha256"] for item in manifest["files"]}
        gate(
            "handoff_hash",
            hashes.get("handoff.json") == evidence.sha256_file(handoff_path),
            "handoff sha256 matches manifest",
        )
        gate(
            "table_hash",
            hashes.get("eligibility_table.csv") == evidence.sha256_file(table_path),
            "table sha256 matches manifest",
        )
    else:
        gate("manifest", False, "manifest.json or handoff.json missing")

    handoff = json.loads(handoff_path.read_text(encoding="utf-8")) if handoff_path.exists() else {}
    gate(
        "handoff_partition",
        len(handoff.get("include", []))
        + len(handoff.get("uncertain", []))
        + len(handoff.get("exclude", []))
        == len(rows),
        "include + uncertain + exclude = every row",
    )
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
    print(f"artifacts: {out_dir}")
    ok, lines = check_artifacts(out_dir, Path(args.compounds))
    print("\n".join(lines))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
