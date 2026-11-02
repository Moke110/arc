#!/usr/bin/env python
"""Build + verify the derived Tahoe log1p-CP10k memmap cache (ticket #5 memory rule).

    tools/build_cache.py --data-dir <06.3_input/output> --cache <out>/cache/tahoe_log1p_cp10k_v1.npy

Hashes the frozen h5ad first (fail closed), streams the CSR layer into a dense
float32 ``.npy``, records the full-array checksum + a random-sample equality
check against the frozen file, and writes the sidecar the trainer verifies.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ptlib import cache as cache_mod  # noqa: E402
from ptlib import manifest as manifest_mod  # noqa: E402
from ptlib import paths  # noqa: E402
from ptlib import runtime  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build the derived Tahoe log1p memmap cache.")
    parser.add_argument("--data-dir", type=Path, default=paths.DEFAULT_FROZEN_INPUT)
    parser.add_argument("--cache", type=Path,
                        default=paths.PIPELINE / "output" / "cache" / cache_mod.DEFAULT_CACHE_NAME)
    parser.add_argument("--verify", choices=("metadata", "deep"), default="deep",
                        help="deep (default) also verifies the existing cache checksum")
    ns = parser.parse_args(argv)
    runtime.setup_logging(None)
    raw = manifest_mod.read_json(ns.data_dir / "manifest.json")
    hashes = manifest_mod.normalize_artifact_manifest(raw)
    source_sha = hashes.get("tahoe_transfer.h5ad")
    if not source_sha:
        raise SystemExit("frozen manifest has no tahoe_transfer.h5ad hash")
    if ns.cache.exists():
        cache = cache_mod.open_cache(
            ns.cache, source_sha,
            source_h5ad=ns.data_dir / "tahoe_transfer.h5ad",
            verify=ns.verify,
        )
        print(json.dumps({"status": "existing-ok", "path": str(ns.cache),
                          "sha256": cache.record.get("file_sha256"),
                          "sample_max_abs_diff": cache.sample_max_abs_diff}, indent=2))
        return 0
    record = cache_mod.build_cache(
        ns.data_dir / "tahoe_transfer.h5ad", ns.cache, source_sha, log=None
    )
    print(json.dumps(record, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
