"""Exercise the real Stage 2 releases with the local archive absent."""

import json
from pathlib import Path
import sys

import numpy as np

STAGE = Path(__file__).resolve().parents[2]
RAW = STAGE / "01_data_raw/output"
TRAINING = STAGE / "03_data_train"
TAHOE = TRAINING / "03.2_tahoe_hvg-5k"
sys.path[:0] = [str(path) for path in (
    RAW / "api", STAGE / "02_gene_axis/02.1_hvg-deg/output/api",
    TRAINING / "execution", TRAINING / "output/api",
    TAHOE / "execution", TAHOE / "output/api", STAGE / "util/output/api",
)]

from build_tenv_training import FOLD_RECEIPT, FOLD_ROOT, verify_tenv
from native_corpus import open_corpus
from prediction_axis import open_axis, open_axes
from publish_tahoe import verify_release
from tahoe_expression import open_tahoe
from tenv_axes import TENV_AXIS_STEPS, open_tenv_axis_training, tenv_axis_output
from training_data import TENV_BASELINE_AXIS


def verify_active_inputs():
    if (STAGE / "history").exists():
        raise ValueError("archive-isolated smoke requires an unavailable Stage 2 history")
    axes = {}
    with open_corpus(RAW) as corpus:
        source = corpus.source("tenv")
        treated = np.flatnonzero(source.obs.role.eq("treated") & source.obs.match_status.eq("matched"))
        row = int(treated[0])
        controls = source.matched_control_rows(row)
        matches = source.obs.iloc[list(controls)]
        if (not len(matches) or not matches.plate.eq(source.obs.plate.iloc[row]).all()
                or not matches.vehicle.eq("DMSO").all()):
            raise ValueError("native tenv treatment lacks a same-plate DMSO match")
        for index in (row, *controls):
            counts = source.counts(int(index), int(index) + 1)
            if not np.issubdtype(counts.dtype, np.integer) or (counts < 0).any():
                raise ValueError("native tenv count window must be nonnegative integer")
        native = source.counts(row, row + 1)
        for family in ("hvg", "denoise"):
            release = open_axes(family=family)
            for size in release.sizes:
                axis = release.axis(size)
                projected = axis.project(corpus, "tenv").counts(row, row + 1)
                np.testing.assert_array_equal(projected, native[:, source.native_columns(axis.gene_indices)])
                axes[f"{family}-{size // 1000}k"] = axis.sha256
        fold_axis = open_axis(TENV_BASELINE_AXIS)
        baseline = tenv_axis_output(TENV_BASELINE_AXIS)
        tenv = {}
        for axis_id in TENV_AXIS_STEPS:
            axis = open_axis(axis_id)
            root = tenv_axis_output(axis_id)
            tenv[axis_id] = verify_tenv(root, corpus, axis, fold_axis, FOLD_ROOT, FOLD_RECEIPT,
                                       reference=None if root == baseline else baseline)
            with open_tenv_axis_training(axis_id) as training:
                positions = training.treatments.control_row.iloc[:4].to_numpy(dtype=int)
                matched = np.stack([training.read_controls(index, index + 1)[0] for index in positions])
                np.testing.assert_array_equal(training.read_treatments(0, 4, "delta"),
                                              training.read_treatments(0, 4) - matched)
                tenv[axis_id]["folds"] = training.manifest["inputs"]["folds"]
        raw = {"manifest_sha256": corpus.manifest_sha256, "dataset_count": len(corpus.manifest["sources"]),
               "tenv_treatment_row": row, "tenv_matched_vehicle_rows": list(map(int, controls)),
               "native_count_dtype": str(native.dtype), "native_width": native.shape[1]}
    tahoe_root = TAHOE / "output/release"
    tahoe = verify_release(tahoe_root)
    with open_tahoe(tahoe_root) as view:
        for start in (0, view.shape[0] - 4):
            stop = start + 4
            np.testing.assert_array_equal(view.read(start, stop, "delta"),
                                          view.read(start, stop) - view.read(start, stop, "control"))
    return {"result": "PASS", "history_available": False, "raw": raw,
            "axes": axes, "training": {"tenv": tenv, "tahoe": tahoe}}


if __name__ == "__main__":
    print(json.dumps(verify_active_inputs(), sort_keys=True, indent=2))
