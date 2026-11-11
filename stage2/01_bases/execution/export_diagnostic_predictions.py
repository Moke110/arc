"""Export one diagnostic fold's measured delta rows for the public evaluator (#156 E5).

Gate E5 proves the standalone evaluator works end to end through the stable
paths on the accepted successor envelope. There are no surviving model
predictions, so the predictions are the measured delta layer of the current
training release for the requested fold/role. Rows are read in bounded blocks
through h5py; the whole matrix is never materialized.

    python3 stage2/01_bases/execution/export_diagnostic_predictions.py \
        --training stage2/01_bases/output/training-corpora \
        --cohorts stage2/01_bases/output/tenv-hvg/cohorts.csv \
        --fold 0 --role test --output <empty candidate dir>
"""
import argparse
import json
from pathlib import Path
import sys

from anndata.io import read_elem
import h5py
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--training", type=Path, required=True,
                        help="training release directory (its tenv.h5ad is read)")
    parser.add_argument("--cohorts", type=Path, required=True)
    parser.add_argument("--fold", type=int, required=True)
    parser.add_argument("--role", default="test", choices=("train", "validation", "test"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--block", type=int, default=64)
    args = parser.parse_args()
    output = args.output
    if output.exists():
        parser.error(f"refusing existing output directory: {output}")
    assignments = pd.read_csv(args.cohorts, keep_default_na=False)
    assignments = assignments[assignments.fold.eq(args.fold) & assignments.role.eq(args.role)]
    if assignments.sample_id.duplicated().any() or assignments.empty:
        raise ValueError("requested fold/role is empty or repeats a sample")
    payload = args.training / "tenv.h5ad"
    output.mkdir(parents=True)
    with h5py.File(payload, "r") as handle:
        obs = read_elem(handle["obs"]).reset_index(drop=True)
        positions = pd.Series(np.arange(len(obs)), index=obs.sample_id.astype(str))
        wanted = assignments[["sample_id"]].copy()
        wanted["_position"] = wanted.sample_id.map(positions)
        if wanted._position.isna().any():
            raise ValueError("requested samples are absent from the training release")
        wanted = wanted.sort_values("_position").reset_index(drop=True)
        rows = wanted._position.to_numpy(dtype=int)
        matrix = handle["layers/delta"]
        if matrix.shape != handle["X"].shape:
            raise ValueError("training release lacks the complete delta layer")
        predictions = np.empty((len(rows), matrix.shape[1]), dtype=np.float32)
        for start in range(0, len(rows), args.block):
            part = rows[start:start + args.block]
            predictions[start:start + len(part)] = matrix[part]
    np.save(output / "predictions.npy", predictions)
    samples = obs.iloc[rows][["sample_id", "drug_parent", "dose_um"]].reset_index(drop=True)
    samples.to_csv(output / "samples.csv", index=False)
    record = {"schema": "bases-156-diagnostic-predictions-v1",
              "source": "measured delta layer of the current training release",
              "training": str(payload.resolve()), "cohorts": str(args.cohorts.resolve()),
              "fold": args.fold, "role": args.role,
              "rows": len(rows), "columns": int(matrix.shape[1]),
              "predictions": str((output / "predictions.npy").resolve()),
              "samples": str((output / "samples.csv").resolve())}
    (output / "export.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
