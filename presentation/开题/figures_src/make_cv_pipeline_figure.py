#!/usr/bin/env python3
"""Slide-4 figure: cardiovascular trials lag other therapeutic areas.

Counts of *registered interventional studies* by start year for four therapeutic
areas, from the ClinicalTrials.gov API v2 (real registry data; not fabricated).

Output -> presentation/figures/cv_pipeline_by_area.png
"""

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/home/chang/arc")
OUT = ROOT / "presentation/figures"
OUT.mkdir(parents=True, exist_ok=True)

YEARS = list(range(2008, 2025))
AREAS = [
    ("Oncology", "cancer", "#C0392B"),
    ("Cardiovascular", "cardiovascular", "#1F4E79"),
    ("Diabetes", "diabetes", "#E8A13C"),
    ("Infectious disease", "infectious disease", "#4A8C5A"),
]
BASE = ("https://clinicaltrials.gov/api/v2/studies?"
        "query.cond={cond}&filter.advanced=AREA%5BStartDate%5DRANGE%5B{y}-01-01,{y}-12-31%5D"
        "&countTotal=true&pageSize=1&fields=NCTId")


def count(cond, year):
    url = BASE.format(cond=urllib.parse.quote(cond), y=year)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.load(r).get("totalCount", 0)
        except Exception:
            time.sleep(2 * (attempt + 1))
    return None


def main():
    plt.rcParams.update({"figure.dpi": 200, "savefig.dpi": 200,
                         "font.size": 16, "axes.labelsize": 18,
                         "xtick.labelsize": 14, "ytick.labelsize": 14,
                         "legend.fontsize": 15, "axes.linewidth": 1.3})
    fig, ax = plt.subplots(figsize=(10, 6))
    for label, cond, color in AREAS:
        vals = []
        for y in YEARS:
            c = count(cond, y)
            vals.append(c)
            time.sleep(0.3)
        ax.plot(YEARS, vals, "-o", color=color, lw=2.4, ms=6, label=label)
    ax.set_xlabel("Trial start year")
    ax.set_ylabel("Registered interventional studies")
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    out = OUT / "cv_pipeline_by_area.png"
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Saved:", out.relative_to(ROOT))


if __name__ == "__main__":
    main()
