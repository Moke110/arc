#!/usr/bin/env python
"""Run one biomni A1 agent per drug from the report final-drug set.

Final-drug set (archived report step): the `10_filter` survivors whose
integrated-score percentile (computed over all screened drugs in
`11_integrate`) exceeds 0.5.
Each drug runs in an isolated subprocess with a wall-clock timeout. Results are
written per-drug to cache/ as JSON, so the batch is resumable (already-done
drugs are skipped). Assembled at the end into cache/results_all.json.

The agent query is the checked-in prompt.md template with {DRUG}/{SMILES}
substituted. The biomni agent backend (/home/chang/cc-research) is an external
dependency and is not vendored.

Usage:
  python report/run_drug_agents.py --per-drug-timeout 1500
"""
import argparse
import json
import os
import re
import sys
import time
from multiprocessing import Process
from pathlib import Path

from dotenv import load_dotenv

import pandas as pd
from scipy.stats import percentileofscore

WORKFLOW_ROOT = Path(__file__).resolve().parent.parent
REPORT_DIR = WORKFLOW_ROOT / "report"
CACHE_DIR = REPORT_DIR / "cache"
PROMPT_PATH = REPORT_DIR / "prompt.md"

# External biomni agent backend (not vendored); override with $BIOMNI_CC_DIR.
CC = os.environ.get("BIOMNI_CC_DIR", "/home/chang/cc-research")
load_dotenv(os.path.join(CC, ".env"))
sys.path.insert(0, os.path.join(CC, "biomni"))

INTEGRATED = WORKFLOW_ROOT / "11_integrate" / "output" / "integrated_scores.csv"
PATH_FILTERED = WORKFLOW_ROOT / "10_filter" / "output" / "path_filtered.csv"

FIELDS = [
    ("ENGLISH_NAME", "English name (with common brand names)"),
    ("CHINESE_NAME", "Chinese name (中文名, standard transliteration if any)"),
    ("US_FDA_APPROVAL", "US FDA approval status (approved / withdrawn / not found; year and indication if known)"),
    ("CHINA_NMPA_APPROVAL", "China NMPA approval status (approved / not found)"),
    ("CLASSIFICATION", "Pharmacological classification (drug class, mechanism, ATC if known)"),
    ("CLINICAL_SCOPE", "Clinical application scope (approved indications and main clinical uses)"),
    ("MIRI_AMI_CLINICAL", "Connection to MIRI/AMI: clinical evidence (trials, post-MI use, epidemiology)"),
    ("MIRI_AMI_PRECLINICAL", "Connection to MIRI/AMI: preclinical evidence (animal models, mechanism, cardioprotection)"),
]


def load_final_drugs():
    """Final-report drugs: 10_filter survivors with integrated pct > 0.5.

    Percentile of each survivor's integrated score vs ALL screened drugs
    (11_integrate integrated_scores.csv), kind='rank' — the same statistic the
    step-11 'integrated pct >= 0.5' top set uses.
    """
    integrated = pd.read_csv(INTEGRATED)
    path_filtered = pd.read_csv(PATH_FILTERED)
    merged = (path_filtered[["drug_name"]]
              .merge(integrated, on="drug_name", how="left")
              .sort_values("integrated", ascending=False)
              .reset_index(drop=True))
    full = integrated["integrated"].values
    pcts = [percentileofscore(full, v, kind="rank") / 100.0
            for v in merged["integrated"].values]
    return merged[pd.Series(pcts).values >= 0.5].copy()


def make_query(drug, smi):
    template = PROMPT_PATH.read_text()
    return template.replace("{DRUG}", drug).replace("{SMILES}", smi)


def run_one(drug, smi, out_path):
    t0 = time.time()
    record = {"drug": drug, "smiles": smi, "ok": False, "time_s": None,
              "solution": None, "error": None}
    try:
        from biomni.agent import A1
        data_path = os.environ.get("BIOMNI_PATH", os.path.join(CC, "biomni", "data"))
        agent = A1(path=data_path, expected_data_lake_files=[])
        log, content = agent.go(make_query(drug, smi))
        content = content or ""
        m = re.search(r"<solution>(.*?)</solution>", content, re.DOTALL | re.IGNORECASE)
        solution = m.group(1).strip() if m else content.strip()
        record["ok"] = True
        record["solution"] = solution
        record["steps"] = len(log)
    except Exception as e:
        record["error"] = f"{type(e).__name__}: {e}"
    record["time_s"] = round(time.time() - t0, 1)
    with open(out_path, "w") as f:
        json.dump(record, f, ensure_ascii=False, indent=1)
    print(f"[{record['time_s']}s] {drug}: ok={record['ok']} "
          f"err={'N/A' if record['ok'] else record['error'][:80]}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-drug-timeout", type=int, default=1500)
    ap.add_argument("--drugs", default=None,
                    help="Comma-separated drug names to restrict to (skip others)")
    args = ap.parse_args()

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    df = load_final_drugs()
    drugs = df[["drug_name", "smiles"]].values.tolist()
    restrict = set(args.drugs.split(",")) if args.drugs else None
    if restrict:
        drugs = [d for d in drugs if d[0] in restrict]

    results = {}
    for drug, smi in drugs:
        out_path = CACHE_DIR / f"{drug}.json"
        if out_path.exists():
            with open(out_path) as f:
                rec = json.load(f)
            print(f"skip (cached) {drug} ok={rec.get('ok')}", flush=True)
            results[drug] = rec
            continue
        p = Process(target=run_one, args=(drug, smi, str(out_path)))
        p.start()
        p.join(timeout=args.per_drug_timeout)
        if p.is_alive():
            p.terminate()
            p.join()
            rec = {"drug": drug, "smiles": smi, "ok": False,
                   "time_s": args.per_drug_timeout, "solution": None,
                   "error": f"wall-clock timeout > {args.per_drug_timeout}s"}
            with open(out_path, "w") as f:
                json.dump(rec, f, ensure_ascii=False, indent=1)
            print(f"[TIMEOUT {args.per_drug_timeout}s] {drug}", flush=True)
            results[drug] = rec
            continue
        with open(out_path) as f:
            rec = json.load(f)
        results[drug] = rec

    with open(CACHE_DIR / "results_all.json", "w") as f:
        json.dump(results, f, ensure_ascii=False, indent=1)
    ok = sum(1 for r in results.values() if r.get("ok"))
    print(f"\nDONE {ok}/{len(drugs)} ok. results_all.json written.")


if __name__ == "__main__":
    main()
