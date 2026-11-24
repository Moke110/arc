"""Verify the archived accepted payloads and replay a sealed prediction."""
import argparse
from hashlib import file_digest
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

STEP = Path(__file__).resolve().parents[1]
ARCHIVE = STEP / "historical/04_model_view"
sys.path.insert(0, str(ARCHIVE / "output/api"))

from tenv5k.inputs import Inputs
from tenv5k.protocol import score_prediction


def digest(path):
    with path.open("rb") as stream:
        return file_digest(stream, "sha256").hexdigest()


def verify(archive):
    original_receipt = archive / "log/231-canonical-preservation.json"
    receipt = json.loads(original_receipt.read_bytes())
    if receipt["schema"] != "tenv5k-canonical-preservation-v1":
        raise ValueError("unknown sealed single-dose preservation record")
    trees = {}
    for tree in receipt["trees"]:
        parts = Path(tree["destination"]).parts
        relative = Path(*parts[parts.index("02_models_tenv-5k") + 1:])
        root = archive / relative
        expected = tree["relative_file_sha256"]
        assets = {path.relative_to(root).as_posix(): digest(path)
                  for path in root.rglob("*") if path.is_file()}
        if assets != expected or len(assets) != tree["files"]:
            raise ValueError(f"sealed single-dose assets differ: {relative}")
        trees[relative.as_posix()] = len(assets)
    axis_local = Inputs(archive / "execution/checkpoints/axis-local-245-input-smoke")
    data = Inputs(archive / "output/inputs")
    run = archive / "output/runs/cold-drug/transigen/0/42"
    rows = data.rows("cold-drug", 0, "test")
    ids = pd.read_csv(run / "identity.csv").profile_id
    prediction = np.load(run / "prediction.npy", mmap_mode="r", allow_pickle=False)
    wells, parents = score_prediction(data.obs.iloc[rows], data.expression[rows], prediction,
                                     data.control[rows], prediction_ids=ids)
    for name, scores in (("per_well", wells), ("per_drug", parents)):
        if scores.to_csv(index=False).encode() != (run / f"{name}.csv").read_bytes():
            raise ValueError(f"archived fixed-prediction playback differs: {name}")
    return {"result": "PASS", "trees": trees, "tested_wells": len(wells), "tested_parents": len(parents),
            "fixed_prediction_sha256": digest(run / "prediction.npy"),
            "original_receipt_sha256": digest(original_receipt),
            "axis_local_manifest_sha256": digest(axis_local.directory / "manifest.json")}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=ARCHIVE)
    args = parser.parse_args()
    print(json.dumps(verify(args.archive), sort_keys=True))
