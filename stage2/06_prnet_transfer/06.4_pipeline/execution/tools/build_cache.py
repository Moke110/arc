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
    parser.add_argument("--source-h5ad", type=Path, default=None,
                        help="default <data-dir>/tahoe_transfer.h5ad; e.g. the 04 pretrain view "
                             "for a pre-build before ticket #3 freezes the input")
    parser.add_argument("--source-sha256", default=None,
                        help="recorded source hash; computed from the file when omitted")
    parser.add_argument("--cache", type=Path,
                        default=paths.PIPELINE / "output" / "cache" / cache_mod.DEFAULT_CACHE_NAME)
    parser.add_argument("--verify", choices=("metadata", "deep"), default="deep",
                        help="deep (default) also verifies the existing cache checksum")
    parser.add_argument("--rebind-frozen", action="store_true",
                        help="existing cache: prove it is byte-identical to the frozen layer "
                             "and record the frozen source sha256 (pre-built from the 04 view)")
    ns = parser.parse_args(argv)
    runtime.setup_logging(None)
    source_h5ad = ns.source_h5ad or (ns.data_dir / "tahoe_transfer.h5ad")
    if not source_h5ad.exists():
        raise SystemExit(f"Tahoe source h5ad not found: {source_h5ad}")
    source_sha = ns.source_sha256
    if not source_sha and (ns.data_dir / "manifest.json").exists() and ns.source_h5ad is None:
        raw = manifest_mod.read_json(ns.data_dir / "manifest.json")
        source_sha = manifest_mod.normalize_artifact_manifest(raw).get("tahoe_transfer.h5ad")
    if not source_sha:
        print(f"hashing source {source_h5ad} …")
        source_sha = manifest_mod.sha256_file(source_h5ad)
    print(f"source sha256={source_sha}")
    if ns.cache.exists():
        if ns.rebind_frozen:
            record = cache_mod.rebind_cache(ns.cache, source_h5ad, source_sha)
            print(json.dumps({"status": "rebound", "path": str(ns.cache),
                              "layer_sha256": record.get("layer_sha256"),
                              "source_sha256": record.get("source_sha256"),
                              "rebound_from": record.get("rebound_from"),
                              "sample_max_abs_diff": record.get("sample_max_abs_diff")},
                             indent=2, default=str))
            return 0
        cache = cache_mod.open_cache(
            ns.cache, source_sha, source_h5ad=source_h5ad, verify=ns.verify,
        )
        print(json.dumps({"status": "existing-ok", "path": str(ns.cache),
                          "sha256": cache.record.get("file_sha256"),
                          "sample_max_abs_diff": cache.sample_max_abs_diff}, indent=2))
        return 0
    record = cache_mod.build_cache(source_h5ad, ns.cache, source_sha, log=None)
    print(json.dumps(record, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
