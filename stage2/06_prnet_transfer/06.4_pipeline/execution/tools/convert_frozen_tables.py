#!/usr/bin/env python
"""Build the derived CSV table cache for the pinned GPU env (no parquet engine).

The pinned training env has no parquet engine, so the frozen ``folds.parquet`` /
``consensus_membership.parquet`` / ``tahoe_control_groups.parquet`` tables are
converted once to CSV under the pipeline's gitignored ``output/input_cache/``.
Each entry records the source parquet's frozen sha256; the loader refuses a
cache entry built from a different source.  Run this with the ``arc`` env (which
has a parquet engine):

    /home/chang/miniconda3/envs/arc/bin/python tools/convert_frozen_tables.py \
        --data-dir <06.3_input/output>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ptlib import manifest as manifest_mod  # noqa: E402
from ptlib import paths  # noqa: E402
from ptlib import tables  # noqa: E402

TABLES = ("folds.parquet", "consensus_membership.parquet", "tahoe_control_groups.parquet")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Derived CSV cache for the frozen parquet tables.")
    parser.add_argument("--data-dir", type=Path, default=paths.DEFAULT_FROZEN_INPUT)
    parser.add_argument("--cache-dir", type=Path, default=tables.DEFAULT_TABLE_CACHE)
    args = parser.parse_args(argv)
    available, engine = tables.parquet_engine_available()
    if not available:
        raise SystemExit(
            "no parquet engine in this interpreter; run this with the `arc` env "
            "(/home/chang/miniconda3/envs/arc/bin/python)"
        )
    raw = manifest_mod.read_json(args.data_dir / "manifest.json")
    hashes = manifest_mod.normalize_artifact_manifest(raw)
    entries = {}
    for name in TABLES:
        if name not in hashes:
            raise SystemExit(f"frozen manifest has no hash for {name}")
        entry = tables.convert_parquet(args.data_dir / name, args.cache_dir, hashes[name])
        entries[name] = entry
        print(f"{name}: {entry['rows']} rows -> {entry['csv']}")
    print(json.dumps({"engine": engine, "cache_dir": str(args.cache_dir),
                      "tables": sorted(entries)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
