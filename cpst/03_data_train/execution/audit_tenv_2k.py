"""Audit the new tenv training releases directly against bounded native count windows."""

import argparse
import json
from pathlib import Path

import numpy as np

from build_training import open_corpus
from prediction_axis import open_axis
from tenv_axes import NEW_TENV_AXIS_STEPS, TENV_AXIS_STEPS, open_tenv_axis_training

STEP = Path(__file__).resolve().parents[1]
RAW = STEP.parent / "01_data_raw/output"
BATCH_ROWS = 64


def audit_counts(training, source, axis):
    np.testing.assert_array_equal(training.genes.gene_index, axis.gene_indices)
    np.testing.assert_array_equal(training.genes.symbol, axis.symbols)
    members = training.membership().set_index("profile_id").reindex(source.obs.index)
    offset = len(training.controls)
    target = members["row"].fillna(-1).to_numpy(dtype=np.int64)
    target[members.role.eq("treatment").to_numpy()] += offset
    columns = source.native_columns(axis.gene_indices)
    pooled = np.zeros((offset + len(training.treatments), len(columns)), dtype=np.int64)
    for start in range(0, len(target), BATCH_ROWS):
        stop = min(start + BATCH_ROWS, len(target))
        rows = target[start:stop]
        retained = rows >= 0
        counts = source.counts(start, stop)
        np.add.at(pooled, rows[retained], counts[np.ix_(retained, columns)])

    totals = pooled.sum(axis=1, dtype=np.int64)
    np.testing.assert_array_equal(totals[:offset], training.controls.library_total)
    np.testing.assert_array_equal(totals[offset:], training.treatments.library_total)
    controls = training.read_controls(0, offset)
    expected = np.log1p(10000 * pooled[:offset] / totals[:offset, None]).astype(np.float32)
    np.testing.assert_array_equal(expected, controls)
    np.testing.assert_array_equal(
        expected.mean(axis=0, dtype=np.float64).astype(np.float32), training.control(),
    )
    for start in range(0, len(training.treatments), BATCH_ROWS):
        stop = min(start + BATCH_ROWS, len(training.treatments))
        block = slice(offset + start, offset + stop)
        expected = np.log1p(10000 * pooled[block] / totals[block, None]).astype(np.float32)
        matched = controls[training.treatments.control_row.iloc[start:stop].to_numpy()]
        np.testing.assert_array_equal(expected, training.read_treatments(start, stop))
        np.testing.assert_array_equal(expected - matched, training.read_treatments(start, stop, "delta"))
    return {"controls": offset, "treatments": len(training.treatments)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=RAW)
    parser.add_argument("--training", type=Path, default=STEP, help="root of the two 2k releases")
    parser.add_argument("--established", type=Path, default=STEP, help="root of the six established releases")
    args = parser.parse_args()
    report = {}
    with open_corpus(args.raw) as corpus:
        source = corpus.source("tenv")
        for axis_id in TENV_AXIS_STEPS:
            new = axis_id in NEW_TENV_AXIS_STEPS
            with open_tenv_axis_training(axis_id, args.training if new else args.established) as data:
                axis = open_axis(axis_id)
                if data.manifest["inputs"]["axis_sha256"] != axis.sha256:
                    raise ValueError(f"{axis_id} training view uses a different prediction axis")
                evidence = {"manifest_sha256": data.manifest_sha256, "axis_sha256": axis.sha256}
                if new:
                    evidence["independently_audited_rows"] = audit_counts(data, source, axis)
                report[axis_id] = evidence
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
