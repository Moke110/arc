"""Build and replay one tenv training release on a published ordered axis."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from build_training import (FOLD_RECEIPT, FOLD_ROOT, STAGE, TENV_RECIPE, _build_dataset,
                            _verify_dataset, load_folds, open_axes, open_corpus,
                            publication_lock, publish_candidate, release)
from training_data import TENV_DEFAULT, open_tenv_training, sha256

STAGING = STAGE.parent / "tmp/stage2-tenv-training-staging"


def verify_tenv(root, corpus, axis, fold_axis, fold_root, fold_receipt):
    """Check the pinned inputs, complete raw ledgers and bounded expression reads."""
    folds = load_folds(fold_root, fold_receipt, "tenv", fold_axis.sha256)
    source_record = corpus.manifest["sources"]["tenv"]
    expected_inputs = {
        "raw_source": "tenv", "raw_sha256": corpus.manifest_sha256,
        "raw_dataset_sha256": source_record["manifest_sha256"],
        "raw_counts_sha256": source_record["sha256"], "axis_sha256": axis.sha256,
        "folds": {**folds.evidence, "axis_sha256": fold_axis.sha256},
    }
    with open_tenv_training(root) as data:
        if data.manifest["inputs"] != expected_inputs or data.manifest["recipe"] != {
            **TENV_RECIPE, "vehicle": "DMSO",
        }:
            raise ValueError("tenv training inputs or recipe differ")
        _verify_dataset(data, corpus.source("tenv").obs, axis)
        return {"result": "PASS", "manifest_sha256": data.manifest_sha256,
                "counts": data.manifest["counts"]}


def publish_tenv(destination, corpus, axis, fold_axis, fold_root=FOLD_ROOT,
                 fold_receipt=FOLD_RECEIPT, *, staging=STAGING, expected_sha256=None):
    """Build from raw, verify and atomically install, optionally requiring an exact digest."""
    destination = Path(destination)
    if destination.resolve() == TENV_DEFAULT.resolve() and expected_sha256 is None:
        raise ValueError("stable publication requires a verified candidate manifest digest")
    if destination.resolve() == TENV_DEFAULT.resolve() and sha256(Path(fold_receipt)) != sha256(FOLD_RECEIPT):
        raise ValueError("stable publication requires the original Step 07 fold receipt")

    def build(candidate):
        _build_dataset(candidate, "tenv", corpus, axis, fold_root, fold_receipt, tenv_fold_axis=fold_axis)

    def verify(candidate):
        result = verify_tenv(candidate, corpus, axis, fold_axis, fold_root, fold_receipt)
        if expected_sha256 is not None and result["manifest_sha256"] != expected_sha256:
            raise ValueError("tenv training manifest differs from the verified candidate")
        return result

    with publication_lock(destination) as output:
        stage = Path(staging).resolve()
        if stage.is_relative_to(output):
            raise ValueError("training staging must be outside the accepted output")
        if output.exists() and not (output / "manifest.json").is_file():
            raise ValueError("refusing to replace unrecognized tenv training output")
        return publish_candidate(release, output, build, verify, staging=stage)


def replay_tenv(destination, corpus, axis, fold_axis, fold_root=FOLD_ROOT,
                fold_receipt=FOLD_RECEIPT, *, staging=STAGING):
    """Rebuild outside the accepted path and compare manifest and all output digests."""
    with publication_lock(destination) as output:
        stage = Path(staging).resolve()
        if stage.is_relative_to(output):
            raise ValueError("training staging must be outside the accepted output")
        args = (corpus, axis, fold_axis, fold_root, fold_receipt)
        accepted = verify_tenv(output, *args)
        candidate = release.candidate(output, stage)
        try:
            _build_dataset(candidate, "tenv", corpus, axis, fold_root, fold_receipt,
                           tenv_fold_axis=fold_axis)
            if verify_tenv(candidate, *args) != accepted:
                raise ValueError("replayed tenv training manifest differs from the accepted release")
        finally:
            shutil.rmtree(candidate)
        return accepted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("publish", "verify"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--staging", type=Path, default=STAGING)
    parser.add_argument("--raw", type=Path, default=STAGE / "01_data_raw/output")
    parser.add_argument("--axis", type=Path, help="published HVG-DEG axis directory")
    parser.add_argument("--size", type=int, default=5000)
    parser.add_argument("--fold-root", type=Path, default=FOLD_ROOT)
    parser.add_argument("--fold-receipt", type=Path, default=FOLD_RECEIPT)
    parser.add_argument("--expected-sha256", help="manifest digest of the verified candidate")
    args = parser.parse_args()
    if args.command == "verify" and args.expected_sha256 is not None:
        parser.error("--expected-sha256 applies to publication only")
    axis = open_axes(args.axis).axis(args.size)
    fold_axis = open_axes().axis(5000)
    if args.output.resolve() == TENV_DEFAULT.resolve() and axis.sha256 != fold_axis.sha256:
        parser.error("stable 03.1 publication requires the accepted HVG-DEG 5k axis")
    with open_corpus(args.raw) as corpus:
        arguments = (args.output, corpus, axis, fold_axis, args.fold_root, args.fold_receipt)
        if args.command == "publish":
            output = publish_tenv(*arguments, staging=args.staging, expected_sha256=args.expected_sha256)
            result = {"result": "PASS", "manifest_sha256": sha256(output / "manifest.json")}
        else:
            result = replay_tenv(*arguments, staging=args.staging)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
