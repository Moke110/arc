#!/usr/bin/env python
"""Aggregate the per-model `FIDELITY.md` tables into `results/fidelity.csv`.

The per-model ledgers are the source of truth (ticket 09); this script only
widens them into one machine-readable table with a `model` column.  It fails
closed on a table row whose cell count does not match its header.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODELS = ["prnet", "cpa", "chemcpa", "biolord", "perturbnet", "squidiff", "cmonge"]
COLUMNS = [
    "model", "component", "author_setting", "run_setting", "class",
    "reason", "expected_effect", "ticket",
]


def parse(path: Path, model: str) -> tuple[list[dict], list[str]]:
    rows: list[dict] = []
    problems: list[str] = []
    header: list[str] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith("|"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) < 6:
            continue
        if cells[0].lower().startswith("component"):
            header = [cell.lower().replace(" ", "_") for cell in cells]
            continue
        if set("".join(cells)) <= set("-: "):
            continue
        if header is None:
            continue
        if len(cells) != len(header):
            problems.append(f"{model}: {len(cells)} cells vs {len(header)} header cells: {line[:80]}")
            continue
        row = dict(zip(header, cells))
        row["model"] = model
        rows.append(row)
    if header is None:
        problems.append(f"{model}: no table header found")
    return rows, problems


def main() -> int:
    rows: list[dict] = []
    problems: list[str] = []
    for model in MODELS:
        parsed, bad = parse(HERE.parent / "models" / model / "FIDELITY.md", model)
        rows.extend(parsed)
        problems.extend(bad)
    if problems:
        print("FAILED:", file=sys.stderr)
        for problem in problems:
            print(" -", problem, file=sys.stderr)
        return 1
    out = HERE.parent / "results" / "fidelity.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} deviation rows -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
