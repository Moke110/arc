#!/usr/bin/env python3
"""
Regenerate the v2-vs-v1 FDA divergence report (ticket 07, gate F).

Compares the v2 FDA screening tail (steps 06-11) against v1 on the aligned
FDA-extension HG rows and writes ``11_integrate/divergence_vs_v1.md``.
Divergence is the expected outcome of HG-alignment, not a pass/fail.

Usage:
    python _lib/compare_v1.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import WORKFLOW_ROOT, WORKSPACE_ROOT  # noqa: E402

V1 = WORKSPACE_ROOT / "miri_reposition" / "xpert-fda-v1"
OUT = WORKFLOW_ROOT / "11_integrate" / "divergence_vs_v1.md"


def main():
    deg1 = pd.read_csv(V1 / "06_perturbed_transcriptome/output/xpert_fda-cm_deg.csv")
    deg2 = pd.read_csv(WORKFLOW_ROOT / "06_perturbed_transcriptome/output/xpert_fda-cm_deg.csv")
    gc = deg1.columns[3:]

    d = np.abs(deg1[gc].values.astype(float) - deg2[gc].values.astype(float))
    pearson = float(np.mean([
        pearsonr(deg1[gc].values[i].astype(float),
                 deg2[gc].values[i].astype(float)).statistic
        for i in range(len(deg1))]))

    r1 = pd.read_csv(V1 / "11_integrate/output/ranked_drugs.csv")
    r2 = pd.read_csv(WORKFLOW_ROOT / "11_integrate/output/ranked_drugs.csv")
    i1 = pd.read_csv(V1 / "11_integrate/output/integrated_scores.csv")
    i2 = pd.read_csv(WORKFLOW_ROOT / "11_integrate/output/integrated_scores.csv")
    mi = i1[["drug_name", "integrated"]].merge(
        i2[["drug_name", "integrated"]], on="drug_name", suffixes=("_v1", "_v2"))

    s1, s2 = set(r1["drug_name"]), set(r2["drug_name"])
    n_wcs1 = len(pd.read_csv(V1 / "10_filter/output/wcs_filtered.csv"))
    n_wcs2 = len(pd.read_csv(WORKFLOW_ROOT / "10_filter/output/wcs_filtered.csv"))

    def table(df):
        return "\n".join(
            f"| {i} | {row['drug_name']} ({row['integrated']:+.3f}) |"
            for i, row in enumerate(df.head(10).iterrows(), start=1))

    lines = [
        "# FDA divergence report — v2 (HG-aligned) vs v1",
        "",
        "Captured after re-running the v2 FDA screening tail (steps 06-11) on the "
        "author-pinned FDA-extension HG rows. Divergence is the point of the "
        "HG-aligned feature, not a failure.",
        "",
        "## Predicted DEG divergence (step 06)",
        "",
        f"- Per-drug Pearson r between v1 and v2 FDA DEG: **{pearson:.3f}** (mean)",
        f"- Max |diff|: {d.max():.3f}; mean |diff|: {d.mean():.4f}",
        "",
        "## Screening funnel",
        "",
        "| Stage | v1 | v2 |",
        "|---|---|---:|",
        f"| FDA drugs screened | {len(deg1)} | {len(deg2)} |",
        f"| WCS < 0 | {n_wcs1} | {n_wcs2} |",
        f"| WCS < 0 + >= 6/8 critical pathways | {len(r1)} | {len(r2)} |",
        "",
        "## Integrated ranking divergence (step 11)",
        "",
        f"- Ranked candidates: v1 **{len(r1)}**, v2 **{len(r2)}**; overlap **{len(s1 & s2)} / {len(s1 | s2)}**.",
        f"- Integrated-score correlation across all drugs: **{mi['integrated_v1'].corr(mi['integrated_v2']):.2f}**.",
        f"- Mean |delta integrated|: {(mi['integrated_v1'] - mi['integrated_v2']).abs().mean():.2f}; "
        f"max: {(mi['integrated_v1'] - mi['integrated_v2']).abs().max():.2f}.",
        "",
        "Top 10 candidates:",
        "",
        "| # | v1 | v2 |",
        "|---:|---|---|",
    ]
    for i, (a, b) in enumerate(zip(
            r1.head(10).itertuples(), r2.head(10).itertuples()), start=1):
        lines.append(f"| {i} | {a.drug_name} ({a.integrated:+.3f}) | {b.drug_name} ({b.integrated:+.3f}) |")

    OUT.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
