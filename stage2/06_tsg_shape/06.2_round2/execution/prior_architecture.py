"""Verify the three sealed originals against the published Step 06.1 comparison."""

from itertools import product

import pandas as pd

from architecture_model import STEP, axis_run
from architecture_roster import CONTROL
from evaluate_shape import read_run as read_round_one
from prior_shape import load_prior
from tsg_axis.metrics import METRICS

PRIOR_STEP = STEP.with_name("06.1_round1")
PRIOR_MATRIX = PRIOR_STEP / "output/formal/prior-347/matrix"
RELEASE_RECEIPT = "log/350-formal.json"
RELEASE_RECEIPT_SHA256 = "a3d4703b86ecfab316b0becbea9960ee2e523b50b7d155943103c2d4af039885"
PUBLISHED_SCORES = "figure/data/fold-scores.csv"
PUBLISHED_BANDWIDTHS = "log/bandwidths.csv"
ORIGINALS = tuple(product((CONTROL,), range(3)))
SCORE_FIELDS = ("fold", "seed", "run_sha256", "metric", "value", "n_conditions",
                "n_eligible", "n_undefined")


def load_originals(matrix):
    release_path = PRIOR_STEP / RELEASE_RECEIPT
    if axis_run.digest(release_path) != RELEASE_RECEIPT_SHA256:
        raise ValueError("Step 06.1 receipt differs from the accepted #350 release")
    release = axis_run.read_json(release_path)
    paths, fits, provenance = load_prior(matrix)
    if provenance != release["prior_release"]:
        raise ValueError("original coordinates differ from the #350 source lineage")

    originals = {}
    for coordinate in ORIGINALS:
        key, evidence = read_round_one(paths[coordinate], False, prior=True)
        path, config, manifest, training, _ = evidence
        if key != coordinate:
            raise ValueError(f"sealed original does not match its fold: {path}")
        trace = pd.read_csv(path / "learning_curve.csv", float_precision="round_trip")
        originals[key] = (path, config, manifest, training, trace)

    published = {name: release["outputs"][name] for name in (PUBLISHED_SCORES, PUBLISHED_BANDWIDTHS)}
    axis_run.check_outputs(PRIOR_STEP, {"outputs": published})
    return originals, {
        "ticket": 350,
        "publication_receipt_sha256": RELEASE_RECEIPT_SHA256,
        "publication_report_sha256": release["report_manifest_sha256"],
        "sealed_source_ticket": provenance["ticket"],
        "sealed_source_revision": provenance["release_revision"],
        "sealed_source_receipt_sha256": provenance["release_receipt_sha256"],
        "sealed_source_matrix_sha256": provenance["matrix_manifest_sha256"],
        "sealed_source_report_sha256": provenance["report_manifest_sha256"],
        "published_files_sha256": published,
        "coordinates": [fits[coordinate] for coordinate in ORIGINALS],
    }


def verify_original_scores(scores, published):
    def original_rows(frame):
        return frame.loc[frame["shape"].eq(CONTROL), list(SCORE_FIELDS)].reset_index(drop=True)

    found, accepted = original_rows(scores), original_rows(published)
    if (list(zip(found["fold"], found["metric"])) != list(product(range(3), METRICS))
            or not found.equals(accepted)):
        raise ValueError("Step 06.2 original scores or eligibility differ from the accepted Step 06.1 publication")
