"""Verify Round 3's two retained references against the accepted publication."""

import argparse
from itertools import product
import json
from pathlib import Path
import sys

import pandas as pd

from measured_roster import SOURCE_SHAPES, coordinates
from measured_scores import SCORE_FIELDS

STEP = Path(__file__).resolve().parents[1]
ROUND_THREE = STEP.with_name("06.3_round3")
sys.path.insert(0, str(ROUND_THREE / "execution"))

from composition_model import axis_run
from prior_composition import load_reused
from prior_unnormalized import load_unnormalized
from tsg_axis.metrics import METRICS

RELEASE_REVISION = "28f5eea4593170466ff6ca4110def1ef3b125cce"
RECEIPT = "log/356-formal.json"
RECEIPT_SHA256 = "cc53a143c5c5a121ef27050e5bc6122469f2ea73f0ce83b332e00ffd26a5ebb2"
PUBLISHED = {
    "figure/data/fold-scores.csv": "fold-scores.csv",
    "log/bandwidths.csv": "bandwidths.csv",
    "log/runs.csv": "runs.csv",
}
SOURCE_SCORE_FIELDS = [*SCORE_FIELDS, "run_sha256", "value", "n_conditions", "n_eligible", "n_undefined"]


def verify_release(release):
    """Bind the Round 3 matrix and scoring artifacts to the accepted #356 receipt."""
    receipt_path = release / RECEIPT
    if axis_run.digest(receipt_path) != RECEIPT_SHA256:
        raise ValueError("Step 06.3 source receipt differs from the accepted #356 release")
    receipt = axis_run.read_json(receipt_path)
    formal = release / "output/formal"
    artifacts = {
        formal / "matrix/manifest.json": "matrix_manifest_sha256",
        formal / "matrix/preflight.json": "preflight_sha256",
        formal / "report/manifest.json": "report_manifest_sha256",
    }
    for path, field in artifacts.items():
        if axis_run.digest(path) != receipt[field]:
            raise ValueError(f"Step 06.3 source artifact differs from the accepted release: {path}")
    report = axis_run.read_json(formal / "report/manifest.json")
    if report["outputs"] != receipt["report_outputs_sha256"] or report["prior_release"] != receipt["prior_release"]:
        raise ValueError("Step 06.3 source scores or lineage differ from the accepted report")
    hashes = {path: receipt["outputs"][path] for path in PUBLISHED}
    axis_run.check_outputs(release, {"outputs": hashes})
    axis_run.check_outputs(formal / "report", {"outputs": {PUBLISHED[name]: digest for name, digest in hashes.items()}})
    return receipt, hashes


def load_sources(release=ROUND_THREE, *, round_one_matrix=None, round_two_matrix=None, prior_unnormalized=None):
    release = Path(release)
    receipt, hashes = verify_release(release)
    if round_one_matrix is None:
        round_one_matrix = release.with_name("06.1_round1") / "output/formal/prior-347/matrix"
    if round_two_matrix is None:
        round_two_matrix = release.with_name("06.2_round2") / "output/formal/matrix"
    if prior_unnormalized is None:
        prior_unnormalized = release / "output/formal/prior-unnormalized"
    round_two, round_two_lineage, _, _ = load_reused(round_two_matrix, round_one_matrix)
    first_attempt, first_attempt_lineage = load_unnormalized(prior_unnormalized)
    round_two_lineage["coordinates"] = [entry for entry in round_two_lineage["coordinates"]
                                         if entry["shape"] in ("0_original", "1_output-delta", "2_output-delta-no-x2")]
    lineages = {"round_two": round_two_lineage, "unnormalized_round_three": first_attempt_lineage}
    if lineages != receipt["prior_release"] or set(round_two) & set(first_attempt):
        raise ValueError("source fit lineage differs from the accepted Round 3 lineage")

    sources = round_two | first_attempt
    all_lineage = {
        (entry["shape"], entry["fold"], entry["seed"]): entry
        for branch in lineages.values() for entry in branch["coordinates"]
    }
    selected = tuple(key for key in coordinates(False) if key[0] in SOURCE_SHAPES)
    names = {prior: shape for shape, prior in SOURCE_SHAPES.items()}
    tables = [pd.read_csv(release / name, float_precision="round_trip") for name in PUBLISHED]
    scores, bandwidths, runs = tables
    retained = pd.MultiIndex.from_tuples([
        (SOURCE_SHAPES[shape], fold, seed) for shape, fold, seed in selected])
    for table in (scores, runs):
        table.drop(index=table.index[~pd.MultiIndex.from_frame(table[SCORE_FIELDS[:3]]).isin(retained)], inplace=True)
        table["shape"] = table["shape"].map(names)
    if (list(runs[SCORE_FIELDS[:3]].itertuples(index=False, name=None)) != list(selected)
            or list(scores[SCORE_FIELDS].itertuples(index=False, name=None)) != [
                (*key, metric) for key, metric in product(selected, METRICS)]):
        raise ValueError("source evidence does not contain one accepted row per reference and metric")

    verified, lineage = {}, []
    for accepted in runs.itertuples(index=False):
        key = accepted.shape, int(accepted.fold), int(accepted.seed)
        shape, fold, seed = key
        source_key = SOURCE_SHAPES[shape], fold, seed
        path, config, replay, train, trace = sources[source_key]
        fit = all_lineage[source_key]
        score_rows = scores.loc[scores["shape"].eq(shape) & scores.fold.eq(fold) & scores.seed.eq(seed)]
        if (not score_rows.run_sha256.eq(accepted.run_sha256).all()
                or replay["registration_sha256"] != accepted.run_sha256
                or fit["registration_sha256"] != accepted.run_sha256
                or fit["run_manifest_sha256"] != accepted.manifest_sha256
                or axis_run.digest(path / "manifest.json") != accepted.manifest_sha256
                or replay["outputs"]["checkpoint.pt"] != accepted.checkpoint_sha256
                or replay["outputs"]["learning_curve.csv"] != accepted.trace_sha256
                or config["shape"] != accepted.source_shape
                or config["protocol"] != accepted.protocol
                or (train["best_epoch"], train["epochs"]) != (accepted.best_epoch, accepted.stopping_epoch)):
            raise ValueError(f"source reference differs from its accepted fit or score evidence: {key}")
        verified[key] = path, config, replay, train, trace
        lineage.append({**fit, "shape": shape, "round_three_shape": source_key[0],
                        "round_three_fit_origin": accepted.source})

    ancestry = {
        name: {**branch, "coordinates": [
            fit for fit in branch["coordinates"]
            if (names.get(fit["shape"]), fit["fold"], fit["seed"]) in verified]}
        for name, branch in lineages.items()
    }
    provenance = {
        "ticket": receipt["ticket"],
        "release_revision": RELEASE_REVISION,
        "release_receipt_sha256": RECEIPT_SHA256,
        "matrix_code_revision": receipt["matrix_code_revision"],
        "matrix_manifest_sha256": receipt["matrix_manifest_sha256"],
        "report_manifest_sha256": receipt["report_manifest_sha256"],
        "published_files_sha256": hashes,
        "ancestry": ancestry,
        "coordinates": lineage,
    }
    return verified, provenance, scores.reset_index(drop=True), bandwidths


def verify_source_scores(scores, bandwidths, accepted_scores, accepted_bandwidths):
    coordinates = pd.MultiIndex.from_frame(accepted_scores[SCORE_FIELDS[:3]])
    selection = pd.MultiIndex.from_frame(scores[SCORE_FIELDS[:3]]).isin(coordinates)
    accepted = accepted_scores[SOURCE_SCORE_FIELDS].set_index(SCORE_FIELDS).sort_index()
    reproduced = scores.loc[selection, SOURCE_SCORE_FIELDS].set_index(SCORE_FIELDS).sort_index()
    if not reproduced.equals(accepted):
        raise ValueError("rescored source predictions or eligibility differ from the accepted Round 3 scores")
    if not bandwidths.sort_values("fold").reset_index(drop=True).equals(
            accepted_bandwidths[bandwidths.columns].sort_values("fold").reset_index(drop=True)):
        raise ValueError("rescored fold bandwidths differ from the accepted Round 3 values")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-release", type=Path, default=ROUND_THREE)
    parser.add_argument("--round-one-matrix", type=Path)
    parser.add_argument("--round-two-matrix", type=Path)
    parser.add_argument("--prior-unnormalized", type=Path)
    args = vars(parser.parse_args())
    release = args.pop("source_release")
    _, provenance, scores, _ = load_sources(release, **args)
    print(json.dumps({"verified_score_rows": len(scores), **provenance}, sort_keys=True))


if __name__ == "__main__":
    main()
