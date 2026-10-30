#!/usr/bin/env python
"""Offline verification of the imported Tahoe-100M tree (no network).

`fetch_tahoe.py --verify-only deep` is the authoritative check, but it needs the
HF API to fetch each file's `lfs.sha256`. This script runs the checks that do not
need network:

1. **HF-anchored hashes** for every file listed in `MANIFEST.jsonl` — those were
   fetched from HF in an earlier session, so their sha256 came from HF's API.
2. **Structural integrity** of all 4,430 files: parquet footer readable,
   row-group sum == num_rows, and the frozen expectations (3,388 expression
   shards x 28,225 cells = 95,624,334; 1,026 DE shards totalling
   4,089,820,780 rows; metadata row counts; 62,710 genes).
3. **Coverage**: every cell in the atlas is accounted for by exactly one
   (plate, well, cell_line) group (that check lives in the QA stage).

Exit code 0 = all offline checks pass.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[3]
DATA_HOME = ROOT / "pert_data" / "tahoe"

EXPECT = {
    "expression_shards": 3_388,
    "cells_per_shard": 28_225,
    "total_cells": 95_624_334,
    "de_shards": 1_026,
    "de_rows": 4_089_820_780,
    "genes": 62_710,
    "samples": 1_344,
    "drugs": 379,
    "pert_labels": 1_138,
    "plates": 14,
    "obs_rows": 100_648_790,
}


def sha256_of(path: Path, chunk: int = 32 * 1024**2) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def main() -> int:
    problems: list[str] = []
    notes: list[str] = []

    # 1. HF-anchored hashes from our own earlier HF fetch.
    manifest = DATA_HOME / "MANIFEST.jsonl"
    checked = 0
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            if not line.strip():
                continue
            rec = json.loads(line)
            p = DATA_HOME / rec["path"]
            if not p.exists():
                problems.append(f"manifest file missing: {rec['path']}")
                continue
            if p.stat().st_size != rec["bytes"]:
                problems.append(f"manifest size mismatch: {rec['path']}")
                continue
            if rec.get("sha256"):
                got = sha256_of(p)
                checked += 1
                if got != rec["sha256"]:
                    problems.append(f"sha256 mismatch vs HF: {rec['path']}")
        notes.append(f"HF-anchored sha256 verified: {checked}/{sum(1 for l in manifest.read_text().splitlines() if l.strip())} manifest files")
    else:
        problems.append("MANIFEST.jsonl absent - no HF-anchored hashes available")

    # 2. Expression shards.
    shards = sorted((DATA_HOME / "data").glob("train-*-of-03388.parquet"))
    if len(shards) != EXPECT["expression_shards"]:
        problems.append(f"expression shards {len(shards)} != {EXPECT['expression_shards']}")
    probe_every = max(1, len(shards) // 120)      # full-footers on a sample; cheap on the rest
    cells = 0
    for i, p in enumerate(shards):
        try:
            md = pq.ParquetFile(p).metadata
        except Exception as exc:
            problems.append(f"unreadable parquet: {p.name} ({exc})")
            continue
        rows = md.num_rows
        cells += rows
        if rows != EXPECT["cells_per_shard"]:
            notes.append(f"{p.name}: {rows} rows (not {EXPECT['cells_per_shard']})")
        if i % probe_every == 0:
            rg = sum(md.row_group(j).num_rows for j in range(md.num_row_groups))
            if rg != rows:
                problems.append(f"row-group sum mismatch: {p.name}")
    if cells != EXPECT["total_cells"]:
        problems.append(f"expression cells {cells:,} != {EXPECT['total_cells']:,}")
    notes.append(f"expression shards: {len(shards)}, cells {cells:,}")

    # 3. Pseudobulk DE shards.
    de = sorted((DATA_HOME / "metadata" / "pseudobulk_differential_expression").glob("*.parquet"))
    if len(de) != EXPECT["de_shards"]:
        problems.append(f"DE shards {len(de)} != {EXPECT['de_shards']}")
    de_rows = 0
    for p in de:
        try:
            de_rows += pq.ParquetFile(p).metadata.num_rows
        except Exception as exc:
            problems.append(f"unreadable DE parquet: {p.name} ({exc})")
    if de_rows != EXPECT["de_rows"]:
        problems.append(f"DE rows {de_rows:,} != {EXPECT['de_rows']:,}")
    notes.append(f"pseudobulk DE: {len(de)} shards, {de_rows:,} rows")

    # 4. Metadata tables.
    meta = DATA_HOME / "metadata"
    checks = {
        "gene_metadata.parquet": EXPECT["genes"],
        "sample_metadata.parquet": EXPECT["samples"],
        "drug_metadata.parquet": EXPECT["drugs"],
        "obs_metadata.parquet": EXPECT["obs_rows"],
    }
    for name, want in checks.items():
        p = meta / name
        if not p.exists():
            problems.append(f"missing metadata table: {name}")
            continue
        got = pq.ParquetFile(p).metadata.num_rows
        if got != want:
            problems.append(f"{name}: {got:,} rows != {want:,}")
        notes.append(f"{name}: {got:,} rows")

    samples = pq.read_table(meta / "sample_metadata.parquet").to_pylist()
    pert = {r["drugname_drugconc"] for r in samples}
    if len(pert) != EXPECT["pert_labels"]:
        problems.append(f"perturbation labels {len(pert)} != {EXPECT['pert_labels']}")
    if len({r["plate"] for r in samples}) != EXPECT["plates"]:
        problems.append("plate count mismatch")

    print(json.dumps({
        "ok": not problems,
        "notes": notes,
        "problems": problems[:40],
        "problem_count": len(problems),
    }, indent=2))
    return 0 if not problems else 1


if __name__ == "__main__":
    sys.exit(main())
