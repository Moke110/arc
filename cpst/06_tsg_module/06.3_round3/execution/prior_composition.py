"""Verify all 12 reused fits and preserve their accepted Round 2 lineage."""

import argparse
from itertools import product
import json
from pathlib import Path

import pandas as pd

from composition_model import ROUND_TWO, axis_run
from composition_roster import REUSED_SEED, ROUND_TWO_SHAPES, SHAPES
from architecture_roster import coordinates as round_two_coordinates
from matrix_architecture import SCHEMA as ROUND_TWO_MATRIX_SCHEMA
from prior_architecture import PRIOR_MATRIX as ROUND_ONE_MATRIX, load_originals
from prior_shape import verify_matrix
from tsg_axis.metrics import METRICS
from verify_architecture import read_run as read_round_two

ROUND_TWO_MATRIX = ROUND_TWO / "output/formal/matrix"
RELEASE_RECEIPT = ROUND_TWO / "log/353-formal.json"
RELEASE_RECEIPT_SHA256 = "aa66defedbbb127f41d3db88a859a469fd53792e035941732d02d8b38ae0b872"
REUSED_SHAPES = {ROUND_TWO_SHAPES[shape]: shape for shape in SHAPES[:4]}
PUBLISHED_SCORES = "figure/data/fold-scores.csv"
PUBLISHED_BANDWIDTHS = "log/bandwidths.csv"


def load_reused(round_two_matrix=None, round_one_matrix=None):
    round_two_matrix = ROUND_TWO_MATRIX if round_two_matrix is None else round_two_matrix
    round_one_matrix = ROUND_ONE_MATRIX if round_one_matrix is None else round_one_matrix
    if axis_run.digest(RELEASE_RECEIPT) != RELEASE_RECEIPT_SHA256:
        raise ValueError("Step 06.2 receipt differs from the accepted #353 release")
    release = axis_run.read_json(RELEASE_RECEIPT)
    report = round_two_matrix.parent / "report"
    for file, field in ((round_two_matrix / "manifest.json", "matrix_manifest_sha256"),
                        (round_two_matrix / "preflight.json", "preflight_sha256"),
                        (report / "manifest.json", "report_manifest_sha256")):
        if axis_run.digest(file) != release[field]:
            raise ValueError(f"Step 06.2 source identity differs from #353: {file}")
    report_receipt = axis_run.read_json(report / "manifest.json")
    if report_receipt["outputs"] != release["report_outputs_sha256"]:
        raise ValueError("Step 06.2 report output identities differ from #353")
    axis_run.check_outputs(report, report_receipt)

    paths, fits = verify_matrix(round_two_matrix, round_two_coordinates(False, include_original=False),
                                ROUND_TWO_MATRIX_SCHEMA)
    originals, original_lineage = load_originals(round_one_matrix)
    if original_lineage != release["prior_release"]:
        raise ValueError("sealed Step 06.1 originals differ from the accepted #353 lineage")

    published = {name: release["outputs"][name] for name in (PUBLISHED_SCORES, PUBLISHED_BANDWIDTHS)}
    axis_run.check_outputs(ROUND_TWO, {"outputs": published})
    if any(report_receipt["outputs"][Path(name).name] != digest for name, digest in published.items()):
        raise ValueError("Step 06.2 report and published reused-score evidence differ")

    scores = pd.read_csv(ROUND_TWO / PUBLISHED_SCORES, float_precision="round_trip")
    scores = scores.loc[scores["shape"].isin(REUSED_SHAPES)].copy()
    if list(scores[["shape", "fold", "seed", "metric"]].itertuples(index=False, name=None)) != list(
        product(REUSED_SHAPES, range(3), (REUSED_SEED,), METRICS)
    ):
        raise ValueError("Step 06.2 reused-score evidence has the wrong coordinate roster")
    scores["shape"] = scores["shape"].map(REUSED_SHAPES)
    scores = scores.reset_index(drop=True)
    bandwidths = pd.read_csv(ROUND_TWO / PUBLISHED_BANDWIDTHS, float_precision="round_trip")

    verified, reused_fits = {}, []
    for prior, shape in REUSED_SHAPES.items():
        for fold in range(3):
            key = prior, fold
            if key in originals:
                run = originals[key]
                fit = original_lineage["coordinates"][fold]
                ticket = original_lineage["sealed_source_ticket"]
            else:
                found, run = read_round_two(paths[key], False)
                if found != key:
                    raise ValueError(f"sealed Step 06.2 run does not match its coordinate: {key}")
                fit, ticket = fits[key], release["ticket"]
            _, config, replay, _, _ = run
            accepted = scores.loc[scores["shape"].eq(shape) & scores["fold"].eq(fold)]
            if (config["seed"] != REUSED_SEED
                    or not accepted["run_sha256"].eq(replay["registration_sha256"]).all()
                    or fit["registration_sha256"] != replay["registration_sha256"]):
                raise ValueError(f"sealed Step 06.2 run has mismatched scoring evidence: {key}")
            verified[shape, fold, REUSED_SEED] = run
            reused_fits.append({**fit, "shape": shape, "source_shape": fit.get("source_shape", prior),
                                "seed": REUSED_SEED, "source_ticket": ticket})

    return verified, {
        "ticket": release["ticket"],
        "release_receipt_sha256": RELEASE_RECEIPT_SHA256,
        "release_revision": "2173902d59dc0daec24116cda2a6bb272b525318",
        "matrix_code_revision": release["matrix_code_revision"],
        "matrix_manifest_sha256": release["matrix_manifest_sha256"],
        "report_manifest_sha256": release["report_manifest_sha256"],
        "published_files_sha256": published,
        "original_lineage": original_lineage,
        "coordinates": reused_fits,
    }, scores, bandwidths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--round-two-matrix", type=Path, default=ROUND_TWO_MATRIX)
    parser.add_argument("--round-one-matrix", type=Path, default=ROUND_ONE_MATRIX)
    args = parser.parse_args()
    _, provenance, _, _ = load_reused(args.round_two_matrix, args.round_one_matrix)
    print(json.dumps(provenance, sort_keys=True))


if __name__ == "__main__":
    main()
