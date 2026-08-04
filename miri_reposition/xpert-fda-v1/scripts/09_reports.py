#!/usr/bin/env python3
"""
Stage 09 — gpt-researcher reports (steps 4-6 of the screening pipeline).

This stage drives the literature/regulatory evidence reports via the
gpt-researcher skill. It also produces the deterministic downstream artifacts
that gate the reports:

  data/reports/4-regulation.csv       — FDA+NMPA regulatory status (agent-driven)
  data/reports/4-approved.csv         — FDA+NMPA dual-approved subset
  data/reports/5-ir_report.csv        — IR + MI evidence scores (agent-driven)
  data/reports/5-report_filtered.csv  — drugs with no MIRI evidence
  data/reports/6-deep_candidates.csv  — no-MIRI drugs with integrated pct >= 0.8

The evidence-gathering reports (4-regulation, 5_efficacy_reports,
6-reposition_report) require the gpt-researcher skill + API keys in .env and
are executed by the main agent (see screen.ipynb). This script performs the
deterministic CSV joins/filters and prints the prompts to run; pass
--run-gptr with API keys configured to invoke gpt-researcher directly.

Usage:
    python scripts/09_reports.py                       # deterministic artifacts only
    python scripts/09_reports.py --run-gptr            # also run gpt-researcher reports
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import percentileofscore

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import WORKFLOW_ROOT, WORKSPACE_ROOT  # noqa: E402

SCORES = WORKFLOW_ROOT / "data" / "scores" / "1-scores.csv"
PATH_FILTERED = WORKFLOW_ROOT / "data" / "ranked" / "3-path_filtered.csv"
REPORTS_DIR = WORKFLOW_ROOT / "data" / "reports"
PROMPTS_DIR = WORKFLOW_ROOT / "prompts"
GTR = WORKSPACE_ROOT / ".agents/skills/gpt-researcher/run_research.py"


def prepare_reports(ir_report_csv):
    """Build the deterministic report-gating artifacts from the IR report CSV."""
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    ir = pd.read_csv(ir_report_csv)
    no_miri = ir[ir["miri_score"] == 0].sort_values("integrated", ascending=False)
    out_cols = [c for c in ["drug_name", "drug_name_cn", "integrated"] if c in no_miri.columns]
    no_miri[out_cols].to_csv(REPORTS_DIR / "5-report_filtered.csv", index=False)

    scores = pd.read_csv(SCORES)
    integrated_full = scores["integrated"].values
    no_miri = no_miri.copy()
    no_miri["integrated_pct"] = [
        percentileofscore(integrated_full, v, kind="rank") / 100.0
        for v in no_miri["integrated"].values
    ]
    deep = no_miri[no_miri["integrated_pct"] >= 0.8].sort_values(
        "integrated", ascending=False)
    deep.to_csv(REPORTS_DIR / "6-deep_candidates.csv", index=False)

    print(f"  no-MIRI drugs: {len(no_miri)}, deep candidates (pct>=0.8): {len(deep)}")
    return no_miri, deep


def run_gptr(prompt_path, output_path, report_type="research_report"):
    """Invoke gpt-researcher for one report stage."""
    prompt = Path(prompt_path).read_text()
    cmd = [sys.executable, str(GTR), "--output", str(output_path),
           "--report-type", report_type, prompt]
    print(f"Running gpt-researcher -> {output_path}")
    subprocess.run(cmd, cwd=WORKSPACE_ROOT, check=True)


def main():
    parser = argparse.ArgumentParser(description="Generate MIRI screening reports (steps 4-6)")
    parser.add_argument("--run-gptr", action="store_true",
                        help="Run gpt-researcher reports (requires API keys in .env)")
    parser.add_argument("--ir-report-csv",
                        default=str(REPORTS_DIR / "5-ir_report.csv"))
    parser.add_argument("--path-filtered", default=str(PATH_FILTERED))
    parser.add_argument("--approved-csv", default=str(REPORTS_DIR / "4-approved.csv"))
    args = parser.parse_args()

    REPORTS_DIR.mkdir(parents=True, exist_ok=True)

    # Step 4 artifact: approval filtering happens agent-side; if a regulation
    # CSV already exists, derive the dual-approved subset.
    reg_csv = REPORTS_DIR / "4-regulation.csv"
    if reg_csv.exists():
        reg = pd.read_csv(reg_csv)
        approved = reg[reg["fda_approval"].str.contains("Approved", na=False)]
        approved.to_csv(REPORTS_DIR / "4-approved.csv", index=False)
        print(f"  Dual-approved drugs: {len(approved)} -> 4-approved.csv")

    # Step 5b / 6a: deterministic filters once 5-ir_report.csv exists.
    if Path(args.ir_report_csv).exists():
        prepare_reports(args.ir_report_csv)
    else:
        print(f"  (no {args.ir_report_csv} yet — run the IR evidence step in screen.ipynb)")

    print("\nPrompts for gpt-researcher steps (run via screen.ipynb or --run-gptr):")
    for p in sorted(PROMPTS_DIR.glob("*.txt")):
        print(f"  {p}")

    if args.run_gptr:
        for prompt, out in [
            (PROMPTS_DIR / "4-regulation.txt", REPORTS_DIR / "4-regulation.md"),
            (PROMPTS_DIR / "5-ir-evidence.txt", REPORTS_DIR / "5_efficacy_reports.md"),
            (PROMPTS_DIR / "6-reposition-report.txt", REPORTS_DIR / "6-reposition_report.md"),
        ]:
            if prompt.exists():
                run_gptr(prompt, out)

    print(f"\nStage 09 complete — artifacts in {REPORTS_DIR}")


if __name__ == "__main__":
    main()
