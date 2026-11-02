#!/usr/bin/env python
"""Write the synthetic frozen-input fixture (test/smoke aid; not a deliverable)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tools.fixture import build_fixture  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build the 06.4 synthetic frozen-input fixture.")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--n-genes", type=int, default=18164)
    parser.add_argument("--n-main-drugs", type=int, default=30)
    parser.add_argument("--n-excluded", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20261021)
    ns = parser.parse_args(argv)
    meta = build_fixture(
        ns.out,
        n_genes=ns.n_genes,
        n_main_drugs=ns.n_main_drugs,
        n_excluded=ns.n_excluded,
        seed=ns.seed,
    )
    counts = {
        "dir": meta["dir"],
        "n_genes": meta["n_genes"],
        "n_folds": meta["n_folds"],
        "n_main_drugs": len(meta["main_drugs"]),
        "n_excluded_drugs": len(meta["excluded_drugs"]),
        "n_artifacts": len(meta["artifact_hashes"]),
    }
    print(json.dumps(counts, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
