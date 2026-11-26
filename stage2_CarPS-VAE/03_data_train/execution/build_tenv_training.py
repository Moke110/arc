"""Build and replay one tenv training release on a published ordered axis."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from build_training import (FOLD_RECEIPT, FOLD_ROOT, STAGE, TENV_RECIPE, _build_dataset,
                            _verify_dataset, load_folds, open_corpus,
                            publication_lock, publish_candidate, release)
from prediction_axis import open_axis
from training_data import (TENV_AXIS_STEPS, TENV_BASELINE_AXIS, TENV_DEFAULT,
                           open_tenv_training, sha256, tenv_axis_output)

STAGING = STAGE.parent / "tmp/stage2-tenv-training-staging"


def _verify_reference_population(data, reference):
    with open_tenv_training(reference) as base:
        if {k: v for k, v in base.manifest["inputs"].items() if k != "axis_sha256"} != {
                k: v for k, v in data.manifest["inputs"].items() if k != "axis_sha256"}:
            raise ValueError("tenv reference raw inputs or frozen fold provenance differ")
        for label, actual, expected in (
            ("control population", data.controls.drop(columns="library_total"),
             base.controls.drop(columns="library_total")),
            ("treatment population or fold assignment", data.treatments.drop(columns="library_total"),
             base.treatments.drop(columns="library_total")),
            ("raw membership ledger", data.membership(), base.membership()),
            ("raw exclusion ledger", data.exclusions(), base.exclusions()),
        ):
            if not actual.equals(expected):
                raise ValueError(f"tenv {label} differs from {reference} "
                                 f"({len(actual)} rows versus {len(expected)})")
        return base.manifest_sha256


def verify_tenv(root, corpus, axis, fold_axis, fold_root, fold_receipt, *, reference=None):
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
        result = {"result": "PASS", "manifest_sha256": data.manifest_sha256,
                  "counts": data.manifest["counts"]}
        if reference is not None:
            result["reference_manifest_sha256"] = _verify_reference_population(data, reference)
        return result


def publish_tenv(destination, corpus, axis, fold_axis, fold_root=FOLD_ROOT,
                 fold_receipt=FOLD_RECEIPT, *, staging=STAGING, expected_sha256=None,
                 reference=None):
    """Build from raw, verify and atomically install, optionally requiring an exact digest."""
    destination = Path(destination)
    for axis_name in TENV_AXIS_STEPS:
        if destination.resolve() != tenv_axis_output(axis_name).resolve():
            continue
        if expected_sha256 is None:
            raise ValueError("stable publication requires a verified candidate manifest digest")
        if sha256(Path(fold_receipt)) != sha256(FOLD_RECEIPT):
            raise ValueError("stable publication requires the original Step 07 fold receipt")
        if axis != open_axis(axis_name):
            raise ValueError(f"stable tenv publication requires the accepted {axis_name} axis")
        if axis_name != TENV_BASELINE_AXIS:
            reference = TENV_DEFAULT if reference is None else Path(reference)
            if sha256(reference / "manifest.json") != sha256(TENV_DEFAULT / "manifest.json"):
                raise ValueError("stable tenv publication requires the accepted HVG-DEG 5k reference")
        break

    def build(candidate):
        _build_dataset(candidate, corpus, axis, fold_axis, fold_root, fold_receipt)

    def verify(candidate):
        result = verify_tenv(candidate, corpus, axis, fold_axis, fold_root, fold_receipt,
                             reference=reference)
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
                fold_receipt=FOLD_RECEIPT, *, staging=STAGING, reference=None):
    """Rebuild outside the accepted path and compare manifest and all output digests."""
    with publication_lock(destination) as output:
        stage = Path(staging).resolve()
        if stage.is_relative_to(output):
            raise ValueError("training staging must be outside the accepted output")
        args = (corpus, axis, fold_axis, fold_root, fold_receipt)
        accepted = verify_tenv(output, *args, reference=reference)
        candidate = release.candidate(output, stage)
        try:
            _build_dataset(candidate, corpus, axis, fold_axis, fold_root, fold_receipt)
            if verify_tenv(candidate, *args, reference=reference) != accepted:
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
    parser.add_argument("--axis", choices=TENV_AXIS_STEPS, default=TENV_BASELINE_AXIS)
    parser.add_argument("--axis-root", type=Path, help="override the published directory for the selected axis")
    parser.add_argument("--fold-root", type=Path, default=FOLD_ROOT)
    parser.add_argument("--fold-receipt", type=Path, default=FOLD_RECEIPT)
    parser.add_argument("--reference", type=Path, help="HVG-DEG 5k release used for population equality")
    parser.add_argument("--expected-sha256", help="manifest digest of the verified candidate")
    args = parser.parse_args()
    if args.command == "verify" and args.expected_sha256 is not None:
        parser.error("--expected-sha256 applies to publication only")
    axis = open_axis(args.axis, args.axis_root)
    fold_axis = open_axis(TENV_BASELINE_AXIS)
    reference = args.reference
    if reference is None and args.axis != TENV_BASELINE_AXIS:
        reference = TENV_DEFAULT
    with open_corpus(args.raw) as corpus:
        arguments = (args.output, corpus, axis, fold_axis, args.fold_root, args.fold_receipt)
        if args.command == "publish":
            output = publish_tenv(*arguments, staging=args.staging, expected_sha256=args.expected_sha256,
                                  reference=reference)
            result = {"result": "PASS", "manifest_sha256": sha256(output / "manifest.json")}
        else:
            result = replay_tenv(*arguments, staging=args.staging, reference=reference)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
