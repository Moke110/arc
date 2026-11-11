"""Reproducible metadata/structure audit of the migrated bases layout (#156).

Successor of ``log/155-reproduce.py`` (the #155 record keeps its own bytes and old
paths). This audit reads the stable products, the moved acceptance records, the
dissolved audit evidence, the preserved historical seals and the published
figure. It never opens ``*.h5ad``, ``*.npy``, ``*.npz`` or ``*.bin`` payloads, so
it is cheap enough to re-run on every change; the one-off payload re-hash is
``execution/verify_payload_digests.py`` (gate E1).

    python3 stage2/01_bases/execution/audit_metadata_bindings.py \
        --step-root /home/chang/arc/stage2/01_bases \
        --write /path/to/log/156-metadata-binding-audit.json

Exit status: 0 when every check is PASS/EXPECTED_GAP, 1 otherwise.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

SCHEMA = "bases-156-metadata-binding-audit-v1"
EXPRESSION = {".h5ad", ".npy", ".npz", ".bin"}
STABLE_DIRS = {"gene-identity", "tenv-raw-counts", "tahoe-raw-counts",
               "shared-expression", "training-corpora", "tenv-truth",
               "tenv-truth-successor", "tenv-hvg", "historical-training-release"}
OLD_NAMES = ("identity-v1", "native-v1", "shared-v1", "shared-v2", "training-v1",
             "training-v2", "truth-tenv-v1", "truth-tenv-v2", "tahoe-counts-v1",
             "hvg-tenv-v1", "hvg-metrics-diagnostic-v1", "build_record")
DIAGNOSTIC_EVIDENCE = {"cohorts.csv": "tenv-hvg/cohorts.csv"}
FIGURE_PNG = "35d14266283b296c1585af3d58f0bacb15a3d7f854cde4972c433375bbdd911f"

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--step-root", type=Path,
                    default=Path(__file__).resolve().parents[1])
parser.add_argument("--figure", type=Path,
                    help="figure root; defaults to step-root/figure (a split checkout can "
                         "carry the tracked code/figure and the payload data in different trees)")
parser.add_argument("--code-root", type=Path,
                    help="checkout that carries the published output/api implementation")
parser.add_argument("--write", type=Path)
args = parser.parse_args()
STEP = args.step_root.resolve()
OUT = STEP / "output"
REC = STEP / "log/records"
HEALTH = STEP / "log/health"
DIAG = STEP / "log/diagnostic"
HIST = STEP / "log/historical"
SEALS = HIST / "seals"
FIG = (args.figure or (STEP / "figure")).resolve()
CODE = (args.code_root or STEP).resolve()

checks, structure = [], []
manifests = {}


def label(path):
    """Step-relative when possible; a split checkout may keep code/figure elsewhere."""
    try:
        return str(Path(path).relative_to(STEP))
    except ValueError:
        return str(path)


def sha(path):
    value = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            value.update(chunk)
    return value.hexdigest()


def add(target, name, status, detail):
    target.append({"name": name, "status": status, "detail": detail})


def digest_map(target, base, declared, *, prefix=""):
    matched, mismatch, absent, skipped, problems = 0, 0, 0, 0, []
    for name, value in declared.items():
        if Path(name).suffix in EXPRESSION:
            skipped += 1
            continue
        if value is None:
            continue
        path = base / name
        if not path.is_file():
            absent += 1
            problems.append({"file": name, "problem": "declared but absent"})
            continue
        actual = sha(path)
        if actual == value:
            matched += 1
        else:
            mismatch += 1
            problems.append({"file": name, "declared": value, "actual": actual})
    add(target, f"{prefix}{base.name}: declared non-expression digests",
        "PASS" if mismatch == 0 and absent == 0 else "FAIL",
        {"matched": matched, "mismatched": mismatch, "absent": absent,
         "expression_skipped": skipped, "problems": problems[:6]})
    return mismatch + absent


def load(relative):
    path = OUT / relative
    manifests[relative] = sha(path)
    return json.loads(path.read_text())


def rows(path):
    return list(csv.DictReader(open(path, newline="")))


# --- identity -----------------------------------------------------------------
idn = load("gene-identity/manifest.json")
add(checks, "gene-identity: genes_sha256 -> genes.csv",
    "PASS" if sha(OUT / "gene-identity/genes.csv") == idn["genes_sha256"] else "FAIL",
    {"genes_sha256": idn["genes_sha256"]})
digest_map(checks, OUT / "gene-identity/tables",
           {name: value["sha256"] for name, value in idn["inputs"].items()
            if name in ("genes.csv", "gene_crosswalk.csv")})
digest_map(checks, HIST / "shared-v1-mappings",
           {name: value["sha256"] for name, value in idn["inputs"].items()
            if name in ("tahoe_mapping.csv", "tenv_mapping.csv")})
annotation = STEP.parents[1] / "pert_data/Annotation/source/hgnc_complete_set.txt"
add(checks, "gene-identity: HGNC annotation input seal",
    "PASS" if sha(annotation) == idn["inputs"]["hgnc_complete_set.txt"]["sha256"] else "FAIL",
    {"path": str(annotation)})

# --- raw counts ---------------------------------------------------------------
nat = load("tenv-raw-counts/tenv_manifest.json")
digest_map(checks, OUT / "tenv-raw-counts", nat["files"])
tk = load("tahoe-raw-counts/manifest.json")
add(checks, "tahoe-raw-counts: mapping_sha256 -> token_mapping.csv",
    "PASS" if sha(OUT / "tahoe-raw-counts/token_mapping.csv") == tk["mapping_sha256"] else "FAIL",
    {"mapping_sha256": tk["mapping_sha256"],
     "source_shards_declared": len(tk["inputs"]["source_files"]),
     "source_shard_rehash": "not run (immutable external source)"})

# --- expression releases ------------------------------------------------------
sh = load("shared-expression/manifest.json")
digest_map(checks, OUT / "shared-expression", dict(sh["files"]))
digest_map(checks, OUT / "gene-identity/tables",
           {name: value for name, value in sh["inputs"].items()
            if name in ("genes.csv", "gene_crosswalk.csv")})
tr = load("training-corpora/manifest.json")
digest_map(checks, OUT / "training-corpora", dict(tr["files"]))
h1 = load("historical-training-release/manifest.json")
digest_map(checks, OUT / "historical-training-release",
           {name: value for name, value in h1["files"].items() if name != "tahoe_bulk.h5ad"})
add(checks, "historical-training-release: rejected Tahoe payload still present",
    "PASS" if (OUT / "historical-training-release/tahoe_bulk.h5ad").is_file()
    and h1["files"]["tahoe_bulk.h5ad"] == "6e3ac2ec880974fe1da69d522eb1e8007a2851b2886510e040719002e4fe548e"
    else "FAIL",
    {"digest": h1["files"]["tahoe_bulk.h5ad"],
     "verdict": "REJECTED_SUPERSEDED; never an input and never opened by this audit"})

# --- truth --------------------------------------------------------------------
truth = load("tenv-truth/manifest.json")
digest_map(checks, OUT / "tenv-truth", dict(truth["files"]))
successor = json.loads((OUT / "tenv-truth-successor/manifest.json").read_text())
add(checks, "tenv-truth-successor: worker_sha256 -> deseq2_contrasts.R",
    "PASS" if sha(OUT / "tenv-truth-successor/deseq2_contrasts.R") == successor["worker_sha256"] else "FAIL",
    {"worker_sha256": successor["worker_sha256"]})
digest_map(checks, OUT / "tenv-truth", dict(successor["source_truth_files"]))
binding_rows = []
for key, item in successor["bindings"].items():
    path = Path(item["path"])
    row = {"binding": key, "recorded_path": item["path"],
           "path_is_resolved_absolute": str(path.resolve()) == str(path),
           "exists": path.exists(), "suffix": path.suffix}
    if path.is_file() and path.suffix not in EXPRESSION:
        row["digest_matches"] = sha(path) == item["sha256"]
    elif path.is_file():
        row["digest_matches"] = "not rehashed (expression payload)"
    binding_rows.append(row)
add(checks, "tenv-truth-successor: recorded binding paths are live absolute dependencies",
    "PASS" if all(row["path_is_resolved_absolute"] and row["exists"] for row in binding_rows) else "FAIL",
    binding_rows)
add(checks, "tenv-truth-successor: bound training release is training-corpora",
    "PASS" if successor["bindings"]["training_manifest"]["path"]
    == str((OUT / "training-corpora/manifest.json").resolve()) else "FAIL",
    {"path": successor["bindings"]["training_manifest"]["path"]})

# --- HVG ----------------------------------------------------------------------
hv = load("tenv-hvg/manifest.json")
declared = {"cohorts.csv": hv["cohorts"]["sha256"]}
for population in hv["populations"].values():
    declared[population["ranking"]] = population["sha256"]
digest_map(checks, OUT / "tenv-hvg", declared)

# --- audit evidence and records ----------------------------------------------
acceptance = json.loads((DIAG / "acceptance.json").read_text())
for name, value in acceptance["files"].items():
    target = (OUT / DIAGNOSTIC_EVIDENCE[name]) if name in DIAGNOSTIC_EVIDENCE else (DIAG / name)
    add(checks, f"log/diagnostic: acceptance.json declares {name}",
        "PASS" if target.is_file() and sha(target) == value else "FAIL",
        {"resolved": str(target.relative_to(STEP)), "declared": value,
         "actual": sha(target) if target.is_file() else None})
visual = json.loads((DIAG / "visual_qa.json").read_text())
add(checks, "log/diagnostic: visual QA seal still matches the accepted PNG",
    "PASS" if visual["sha256"] == FIGURE_PNG and sha(DIAG / "diagnostic_pr.png") == FIGURE_PNG else "FAIL",
    {"sha256": visual["sha256"], "result": visual["result"]})

# --- figure -------------------------------------------------------------------
values = FIG / "data/diagnostic_tenv_pr.csv"
figure_rows = rows(values) if values.is_file() else []
curves = [row for row in figure_rows if row["record"] == "curve"]
summaries = [row for row in figure_rows if row["record"] == "summary"]
report = json.loads((DIAG / "diagnostic_pr.json").read_text())
expected = {(item["label"], metric): item[metric] for item in report["summaries"]
            for metric in ("AUPRC", "AUPRC_positive", "AUPRC_negative")}
observed = {(row["series"], row["panel"]): float(row["macro_ap"]) for row in summaries}
add(structure, "figure: one complete plotted-values table for the diagnostic panel",
    "PASS" if (len(curves) == 7222 and len(summaries) == 9
               and {row["record"] for row in figure_rows} == {"figure", "panel", "summary", "curve"})
    else "FAIL",
    {"rows": len(figure_rows), "curve_rows": len(curves), "summary_rows": len(summaries),
     "columns": list(figure_rows[0]) if figure_rows else []})
add(structure, "figure: displayed macro AP values are in the CSV and equal the accepted panel",
    "PASS" if observed == expected else "FAIL",
    {"mismatches": [{"key": list(key), "csv": observed.get(key), "accepted": expected[key]}
                    for key in expected if observed.get(key) != expected[key]][:5]})
add(structure, "figure: fixed-group definition, uncertainty and reference statement are present",
    "PASS" if all(row.get("group_definition") for row in figure_rows if row["record"] == "figure")
    and all(row.get("uncertainty") for row in figure_rows if row["record"] == "figure")
    and all(row.get("reference") for row in figure_rows if row["record"] == "figure") else "FAIL",
    {"group_definition": next((row["group_definition"] for row in figure_rows
                               if row["record"] == "figure"), None)})
add(structure, "figure: published PNG is the accepted byte-identical panel",
    "PASS" if (FIG / "diagnostic_tenv_pr.png").is_file()
    and sha(FIG / "diagnostic_tenv_pr.png") == FIGURE_PNG else "FAIL",
    {"sha256": sha(FIG / "diagnostic_tenv_pr.png") if (FIG / "diagnostic_tenv_pr.png").is_file() else None})
add(structure, "figure: a renderer exists and the step root no longer publishes duplicates",
    "PASS" if (FIG / "script/render_diagnostic_tenv_pr.py").is_file()
    and not (STEP / "diagnostic_tenv_pr.png").exists()
    and not (STEP / "diagnostic_tenv_pr.csv").exists() else "FAIL",
    {"renderer": label(FIG / "script/render_diagnostic_tenv_pr.py")})

# --- layout -------------------------------------------------------------------
stable = {path.name for path in OUT.iterdir() if path.is_dir()}
published_api = (CODE / "output/api/bases").is_dir() or (OUT / "api/bases").is_dir()
add(structure, "layout: the migrated tree holds exactly the stable products and publishes output/api",
    "PASS" if stable == STABLE_DIRS and published_api else "FAIL",
    {"missing": sorted(STABLE_DIRS - stable), "extra": sorted(stable - STABLE_DIRS),
     "public_api": published_api, "checked_in": str(CODE)})
preserved = SEALS / "truth-successor-envelope-v1/manifest.json"
add(structure, "layout: old product/record names survive only under log/historical",
    "PASS" if preserved.is_file() and not (STEP / "execution/build_record").exists() else "FAIL",
    {"preserved_envelope": str(preserved.relative_to(STEP)) if preserved.is_file() else None})
payloads = [str(path.relative_to(STEP)) for path in (STEP / "log").rglob("*")
            if path.suffix in EXPRESSION]
add(structure, "layout: no expression payload under log/", "PASS" if not payloads else "FAIL",
    {"payloads": payloads[:5]})
add(structure, "layout: repo-health records are under log/health",
    "PASS" if (HEALTH / "final_health/summary.json").is_file() else "FAIL",
    {"summary": str((HEALTH / "final_health/summary.json").relative_to(STEP))
     if (HEALTH / "final_health/summary.json").is_file() else None})

# --- cross-bindings -----------------------------------------------------------
def same(label, left, right):
    add(checks, f"cross: {label}", "PASS" if left == right else "FAIL",
        {"expected": left, "actual": right})


add(checks, "cross: historical-training-release.inputs.shared_manifest -> deleted shared-v1 seal",
    "EXPECTED_GAP",
    {"declared": h1["inputs"]["shared_manifest"],
     "detail": "shared-v1 has no manifest.json: its expression payload and seal were deleted with the "
               "historical release; only the mapping CSVs survive under log/historical/shared-v1-mappings, "
               "so the historical release stays a seal, not a reconstruction input"})
same("training-corpora.inputs.shared_manifest -> shared-expression-manifest",
     tr["inputs"]["shared_manifest"], manifests["shared-expression/manifest.json"])
same("tenv-truth.provenance.training_manifest_sha256 -> historical-training-release manifest",
     truth["provenance"]["training_manifest_sha256"],
     manifests["historical-training-release/manifest.json"])
same("tenv-truth.provenance.native_counts_sha256 -> tenv-raw-counts tenv.h5ad",
     truth["provenance"]["native_counts_sha256"], nat["files"]["tenv.h5ad"])
same("tenv-truth.provenance.tenv_training_sha256 -> historical tenv.h5ad",
     truth["provenance"]["tenv_training_sha256"], h1["files"]["tenv.h5ad"])
same("tenv-truth.provenance.worker_sha256 -> successor worker seal",
     truth["provenance"]["worker_sha256"], successor["worker_sha256"])
same("tenv-truth manifest group_count -> groups.csv CSV records",
     truth["group_count"], len(rows(OUT / "tenv-truth/groups.csv")))
same("truth-complete-verify.manifest_sha256 -> tenv-truth/manifest.json",
     json.loads((REC / "truth_complete_verify.json").read_text())["manifest_sha256"],
     manifests["tenv-truth/manifest.json"])
same("tenv-hvg.training_sha256 -> training-corpora tenv.h5ad",
     hv["training_sha256"], tr["files"]["tenv.h5ad"])
same("tenv-hvg.cohorts.sha256 -> tenv-hvg/cohorts.csv",
     hv["cohorts"]["sha256"], sha(OUT / "tenv-hvg/cohorts.csv"))
historical_successor = json.loads((REC / "truth_successor_acceptance.json").read_text())
superseded = sha(preserved) if preserved.is_file() else None
same("historical truth-successor acceptance binds the preserved superseded envelope",
     historical_successor["manifest_sha256"], superseded)
hvg_acceptance = json.loads((REC / "hvg_metrics_actual_acceptance.json").read_text())["acceptance"]
same("hvg acceptance.hvg_manifest_sha256 -> tenv-hvg/manifest.json",
     hvg_acceptance["hvg_manifest_sha256"], manifests["tenv-hvg/manifest.json"])
same("hvg acceptance.truth_manifest_sha256 -> tenv-truth/manifest.json",
     hvg_acceptance["truth_manifest_sha256"], manifests["tenv-truth/manifest.json"])
same("hvg acceptance.training_sha256 -> training-corpora tenv.h5ad",
     hvg_acceptance["training_sha256"], tr["files"]["tenv.h5ad"])
same("hvg acceptance.successor_manifest_sha256 -> preserved superseded envelope",
     hvg_acceptance["successor_manifest_sha256"], superseded)
acceptance_files = {}
for name, value in hvg_acceptance["files"].items():
    target = (OUT / DIAGNOSTIC_EVIDENCE[name]) if name in DIAGNOSTIC_EVIDENCE else (DIAG / name)
    acceptance_files[name] = {"path": str(target.relative_to(STEP)),
                              "digest_matches": target.is_file() and sha(target) == value}
add(checks, "hvg acceptance: its recorded diagnostic files are unchanged evidence",
    "PASS" if all(item["digest_matches"] for item in acceptance_files.values()) else "FAIL",
    acceptance_files)

# --- structural / scientific invariants ---------------------------------------
axis = rows(OUT / "gene-identity/genes.csv")
add(structure, "identity axis: 18164 rows, gene_index 0..18163, unique symbol/HGNC",
    "PASS" if (len(axis) == 18164 and [int(row["gene_index"]) for row in axis] == list(range(18164))
               and len({row["symbol"] for row in axis}) == 18164
               and len({row["hgnc_id"] for row in axis}) == 18164) else "FAIL",
    {"rows": len(axis), "index_5091": axis[5091]["symbol"]})
tables = rows(OUT / "gene-identity/tables/genes.csv")
add(structure, "frozen axis order identical in gene-identity tables and the identity product",
    "PASS" if [row["symbol"] for row in tables] == [row["symbol"] for row in axis] else "FAIL",
    {"rows": len(tables),
     "shared_expression_genes_equals_tables_bytes":
     sha(OUT / "gene-identity/tables/genes.csv") == sha(OUT / "shared-expression/genes.csv")})

groups = rows(OUT / "tenv-truth/groups.csv")
ids = [row["group_id"] for row in groups]
add(structure, "tenv-truth groups: 1377 unique group ids",
    "PASS" if len(groups) == 1377 and len(set(ids)) == 1377 else "FAIL",
    {"rows": len(groups), "statuses": sorted({row["status"] for row in groups}),
     "sum_treated": sum(int(row["n_treated"]) for row in groups),
     "sum_controls": sum(int(row["n_controls"]) for row in groups),
     "tested_empty": sum(1 for row in groups
                         if row["status"] == "eligible" and int(row["n_significant"]) == 0),
     "multiplate": sum(1 for row in groups if int(row["n_plates"]) > 1)})
files = list((OUT / "tenv-truth/groups").glob("*.json"))
add(structure, "tenv-truth per-group payloads cover exactly the group table",
    "PASS" if {path.stem for path in files} == set(ids) else "FAIL", {"payloads": len(files)})
contrast = rows(OUT / "tenv-truth/contrast_samples.csv")
add(structure, "tenv-truth contrast membership: 4206 treated wells plus controls",
    "PASS" if sum(1 for row in contrast if row["condition"] == "treated") == 4206 else "FAIL",
    {"rows": len(contrast),
     "treated": sum(1 for row in contrast if row["condition"] == "treated"),
     "control": sum(1 for row in contrast if row["condition"] == "control")})

cohorts = rows(OUT / "tenv-hvg/cohorts.csv")
per_fold = {}
for row in cohorts:
    per_fold.setdefault(row["fold"], {}).setdefault(row["role"], 0)
    per_fold[row["fold"]][row["role"]] += 1
roles = {}
for row in cohorts:
    roles.setdefault(row["sample_id"], set()).add((row["fold"], row["role"]))
complete = all(len(value) == 3 and sum(1 for fold, role in value if role == "test") == 1
               for value in roles.values())
add(structure, "HVG cohorts: every well in exactly three folds with exactly one test role",
    "PASS" if complete and len(roles) == 4206 else "FAIL",
    {"samples": len(roles), "rows": len(cohorts), "per_fold_role": per_fold})
for population in ("descriptive", "fold_0", "fold_1", "fold_2"):
    ranks = sorted(int(row["rank"]) for row in rows(OUT / f"tenv-hvg/{population}.csv"))
    add(structure, f"tenv-hvg/{population}: rank permutation over 18164 genes",
        "PASS" if ranks == list(range(1, 18165)) else "FAIL",
        {"rows": len(ranks), "min": ranks[0], "max": ranks[-1]})

v1 = {row["symbol"]: row for row in rows(HIST / "shared-v1-mappings/tahoe_mapping.csv")}
v2 = {row["symbol"]: row for row in rows(OUT / "shared-expression/tahoe_mapping.csv")}
symbol_delta = sum(1 for key in v1 if v1[key]["native_symbol"] != v2[key]["native_symbol"])
token_delta = sum(1 for key in v1 if v1[key]["native_token_id"] != v2[key]["native_token_id"])
add(structure, "historical vs current Tahoe mapping: same tokens, annotated symbol vintage differs by 229",
    "PASS" if symbol_delta == 229 and token_delta == 0 else "FAIL",
    {"rows": len(v1), "native_symbol_differences": symbol_delta,
     "native_token_id_differences": token_delta,
     "example": {"symbol": "ACKR5", "v1": v1["ACKR5"]["native_symbol"],
                 "v2": v2["ACKR5"]["native_symbol"]}})
add(structure, "historical vs current tenv mapping: byte-identical evidence",
    "PASS" if sha(HIST / "shared-v1-mappings/tenv_mapping.csv")
    == sha(OUT / "shared-expression/tenv_mapping.csv") else "FAIL",
    {"rows": len(rows(HIST / "shared-v1-mappings/tenv_mapping.csv"))})

payload = {"schema": SCHEMA, "ticket": 156, "parent": 154, "step_root": str(STEP),
           "scope": "stable products, moved records, dissolved audit evidence, preserved seals "
                    "and the published figure; expression payloads never opened",
           "expression_payloads_opened": 0,
           "manifest_sha256": manifests, "checks": checks, "structure": structure}
failed = [item for item in checks + structure if item["status"] not in ("PASS", "EXPECTED_GAP")]
payload["summary"] = {"checks": len(checks), "structure": len(structure),
                      "non_pass": [item["name"] for item in failed],
                      "expected_gaps": [item["name"] for item in checks + structure
                                        if item["status"] == "EXPECTED_GAP"]}
if args.write:
    args.write.parent.mkdir(parents=True, exist_ok=True)
    args.write.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(json.dumps(payload["summary"], indent=2))
for item in failed:
    print("NON-PASS:", item["name"], json.dumps(item["detail"])[:600])
sys.exit(0 if not failed else 1)
