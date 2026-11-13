#!/usr/bin/env python3
"""Publish or replay the Stage 2 raw corpus (ADR 0066).

``output/corpus/`` holds integer counts on the frozen 18,164 ``gene_index`` axis
with controls as ordinary rows: ``tenv`` 4,206 treated + 272 vehicle and
``tahoe_bulk`` 65,318 treated + 1,400 ``DMSO_TF``, each row carrying its role,
source identity, plate, cell line, desalted drug-parent structure, dose, source
item digest and ``matched_control_group``. No normalization, no delta, no HVG
subsetting.

This driver only binds paths; row metadata, matched-control correspondence and
the reader live once in ``output/api/bases/corpus.py``.

    python stage2/01_bases/execution/build_corpus.py publish
    python stage2/01_bases/execution/build_corpus.py verify
"""
import argparse
import json
import sys
from pathlib import Path

import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import corpus

STEP = Path(__file__).resolve().parents[1]
REPO = next(path for path in Path(__file__).resolve().parents if (path / ".git").exists())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("publish", "verify"))
    parser.add_argument("--raw-counts", type=Path,
                        default=STEP / "output/raw_counts")
    parser.add_argument("--axis", type=Path,
                        default=STEP / "output/gene-identity/tables/genes.csv")
    parser.add_argument("--drug-metadata", type=Path,
                        default=REPO / "pert_data/Tahoe-100M/source/metadata/drug_metadata.parquet")
    parser.add_argument("--staging", type=Path, default=STEP / "execution/checkpoints")
    parser.add_argument("--output", type=Path, default=STEP / "output/corpus")
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    kwargs = {"raw_counts_dir": args.raw_counts, "axis_path": args.axis,
              "drug_metadata": args.drug_metadata, "repo_root": REPO}
    if args.action == "publish":
        corpus.publish_release(args.output, root=args.staging, **kwargs)
    else:
        corpus.replay_release(args.output, root=args.staging, **kwargs)
    manifest = json.loads((args.output / "manifest.json").read_text())
    result = {
        "action": args.action,
        "manifest_sha256": corpus.sha256(args.output / "manifest.json"),
        "counts": manifest["counts"],
        "sources": {source: manifest["sources"][source]["counts"]
                    for source in manifest["sources"]},
        "controls": manifest["controls"],
        "outputs": sorted(manifest["outputs"]),
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(text)
    sys.stdout.write(text)


if __name__ == "__main__":
    main()
