#!/usr/bin/env python3
"""Publish or replay the Stage 2 raw-count release on the frozen gene axis.

ADR 0066: ``output/raw_counts/`` is the replayable ``source/`` -> axis bridge,
not a model-facing product (integer counts only; no normalization, no delta, no
axis subsetting). This driver only binds paths: the release schema, readers and
manifest live once in ``output/api/bases/raw_counts.py``.

    python stage2/01_bases/execution/build_raw_counts.py publish
    python stage2/01_bases/execution/build_raw_counts.py verify

``publish`` writes a candidate under ``execution/checkpoints/`` (never inside
``output/``), verifies it, and only then replaces ``output/raw_counts/``.
``verify`` rebuilds from the same pinned inputs into a scratch candidate and
compares every digest; it never writes the accepted release.

Default mode rehomes the source-accepted payloads (ADR 0066 *Removals*): their
bytes are re-checked against the recorded source-replay acceptance records and
the Tahoe token mapping is re-derived from the immutable source metadata.
``--full-rebuild`` instead rebuilds both payloads from ``source/`` through the
tenv/Tahoe builders and replays them from source; that is the stronger, much
longer mode (recorded: Tahoe build 1 h 35 m, independent source replay 2 h 35 m
under the 10 GiB cap).
"""
import argparse
import json
import sys
from pathlib import Path

import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import raw_counts

STEP = Path(__file__).resolve().parents[1]
REPO = next(path for path in Path(__file__).resolve().parents if (path / ".git").exists())


def build_kwargs(args, mode):
    return {
        "mode": mode,
        "data_root": args.data_root,
        "axis_path": args.axis,
        "crosswalk_path": args.crosswalk,
        "tenv_payload": args.tenv_payload,
        "tahoe_payload": args.tahoe_payload,
        "tenv_acceptance": args.tenv_acceptance,
        "tahoe_acceptance": args.tahoe_acceptance,
        "gene_metadata": args.gene_metadata,
        "checkpoints_root": args.staging / "raw-counts-tahoe",
        "repo_root": REPO,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("publish", "verify"))
    parser.add_argument("--full-rebuild", action="store_true",
                        help="rebuild both payloads from immutable source/ instead of rehoming")
    parser.add_argument("--data-root", type=Path, default=REPO / "pert_data")
    parser.add_argument("--axis", type=Path,
                        default=STEP / "output/gene-identity/tables/genes.csv")
    parser.add_argument("--crosswalk", type=Path,
                        default=STEP / "output/gene-identity/tables/gene_crosswalk.csv")
    parser.add_argument("--tenv-payload", type=Path,
                        default=STEP / "output/tenv-raw-counts")
    parser.add_argument("--tahoe-payload", type=Path,
                        default=STEP / "output/tahoe-raw-counts")
    parser.add_argument("--tenv-acceptance", type=Path,
                        default=REPO / raw_counts.TENV_ACCEPTANCE)
    parser.add_argument("--tahoe-acceptance", type=Path,
                        default=REPO / raw_counts.TAHOE_ACCEPTANCE)
    parser.add_argument("--gene-metadata", type=Path,
                        default=REPO / "pert_data/Tahoe-100M/source/metadata/gene_metadata.parquet")
    parser.add_argument("--staging", type=Path, default=STEP / "execution/checkpoints")
    parser.add_argument("--output", type=Path, default=STEP / "output/raw_counts")
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    mode = "source-rebuild" if args.full_rebuild else "source-accepted-replay"
    kwargs = build_kwargs(args, mode)
    if args.action == "publish":
        raw_counts.publish_release(args.output, root=args.staging, **kwargs)
    else:
        raw_counts.replay_release(args.output, root=args.staging, **kwargs)
    manifest = json.loads((args.output / "manifest.json").read_text())
    result = {
        "action": args.action,
        "mode": mode,
        "manifest_sha256": raw_counts.sha256(args.output / "manifest.json"),
        "counts": manifest["counts"],
        "outputs": sorted(manifest["outputs"]),
        "provenance": {source: manifest["sources"][source]["provenance"]["mode"]
                       for source in ("tenv", "tahoe")},
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(text)
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
