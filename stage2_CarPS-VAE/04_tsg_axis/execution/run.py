"""Register, fit and exactly replay one cold-drug TranSiGen axis coordinate."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import resource
import secrets
import sys
from tempfile import TemporaryDirectory
import time

import numpy as np
import pandas as pd
import torch

STEP = Path(__file__).resolve().parents[1]
REPO = next(parent for parent in STEP.parents if (parent / ".git").exists())
STAGE = REPO / "stage2_CarPS-VAE"
AXIS_API = STAGE / "02_gene_axis/02.1_hvg-deg/output/api"
sys.path[:0] = [str(STEP / "output/api"), str(AXIS_API), str(STAGE)]

from prediction_axis import axis_coordinates, open_axes
from tsg_axis.inputs import open_inputs, score_fixed_prediction
from tsg_axis.protocol import MAX_EPOCHS, PATIENCE, PROTOCOL, RUN_SCHEMA
from tsg_axis import model as kpgt
from tsg_axis.model import DoseView
from tsg_axis.learning import EarlyStopping, fixed_rng
from stage2_util.derivation.release import check_outputs, read_json, sha256 as digest, write_json
from stage2_util.publication import publication_lock
from training_data import TENV_AXIS_STEPS
from util.training import transigen
from util.training.transigen import (
    PREDICTION_FILES, PROCESS, TRAINING_FILES, file_hashes, fit, predict_test,
    prediction, record_predictions, validate, verified_receipt,
)

SEED = 42
VALIDATION_SEED = 22617
TEST_SEED = 100042
PARTITIONS = ("train", "valid", "test")
SMOKE_LIMITS = dict(zip(PARTITIONS, (64, 24, 24), strict=True))


def code_hashes():
    api = STEP / "output/api/tsg_axis"
    paths = [
        Path(__file__),
        *sorted((STAGE / "util/training").glob("*.py")),
        AXIS_API / "prediction_axis.py",
        *(api / name for name in ("inputs.py", "protocol.py", "model.py", "native_import.py", "learning.py")),
        STAGE / "03_data_train/output/api/training_data.py",
        STAGE / "util/output/api/stage2_util/representations.py",
    ]
    return {path.relative_to(REPO).as_posix(): digest(path) for path in paths}


def bind(args, data, vectors):
    family, width = axis_coordinates(args["axis"])
    axes = open_axes(family=family)
    axis_sha = axes.axis(width).sha256
    if (axis_sha != data.release.manifest["inputs"]["axis_sha256"]
            or axis_sha != data.release.manifest["outputs"]["genes.csv"] or len(data.genes) != width):
        raise ValueError("training release lacks the registered ordered prediction axis")
    rows = {name: data.rows(args["fold"], name) for name in PARTITIONS}
    parents = {name: set(data.parents[part]) for name, part in rows.items()}
    if (any(not part.size for part in rows.values())
            or sum(map(len, parents.values())) != len(set().union(*parents.values()))):
        raise ValueError("cold-drug partitions require disjoint nonempty parent populations")
    if args["smoke"]:
        rows = {name: part[:SMOKE_LIMITS[name]] for name, part in rows.items()}
    transform = kpgt.fit_transform(data.doses[rows["train"]])
    model = kpgt.configuration(transform)
    model.update(n_genes=width, method_identity="tenv-axis-transigen-kpgt-latent-numeric-v2",
                 control_input="matched plate DMSO for every partition",
                 validation_objective="native TranSiGen MSE + beta-weighted KLD, batch sums per condition",
                 adaptations=[f"{width} gene interface" if item == "5000 gene interface" else item
                              for item in model["adaptations"]
                              if item != "shared splits and fixed validation treated-expression MSE"]
                 + ["shared cold-drug splits and fixed validation native MSE + KLD"])
    model.update(batch_size=64, predictive_head_dim=532, dose_vector_dim=32)
    native = file_hashes(Path(args["source"]), ("src/model.py", "src/train_TranSiGen.py", "src/utils.py"))
    if (native["src/model.py"] != kpgt.SOURCE_DIGEST
            or native["src/train_TranSiGen.py"] != model["config_script_sha256"]):
        raise ValueError("registered native TranSiGen source differs from the pinned author release")
    hardware = (torch.cuda.get_device_name() if args["device"] == "cuda" else platform.machine())
    config = {"schema": RUN_SCHEMA, "protocol": PROTOCOL, "axis_sha256": axis_sha,
            "axis_family_manifest_sha256": axes.manifest_sha256,
            "training_sha256": data.release.manifest_sha256,
            "raw_sha256": data.release.manifest["inputs"]["raw_sha256"],
            "folds": data.release.manifest["inputs"]["folds"],
            "kpgt_manifest_sha256": vectors.manifest_sha256, "kpgt_vectors_sha256": vectors.manifest["files"]["KPGT.npy"],
            "native_sha256": native,
            "code_sha256": code_hashes(), "model": model, "partitions": {key: part.tolist() for key, part in rows.items()},
            "conditions": {key: list(data.condition_ids(part)) for key, part in rows.items()},
            "environment": {"python": platform.python_version(), "torch": str(torch.__version__),
                            "numpy": np.__version__, "cuda": torch.version.cuda,
                            "torch_threads": torch.get_num_threads(), "hardware": hardware},
            "scenario": "cold-drug", "seed": SEED,
            "validation_seed": VALIDATION_SEED, "test_seed": TEST_SEED,
            "max_epochs": 2 if args["smoke"] else MAX_EPOCHS, "patience": PATIENCE}
    return transigen.bind_experiment(args, config)


def register(output, args):
    return transigen.register_run(output, args, bind)


def build_model(args, config, data):
    native = kpgt.adapter(args["source"], config["model"]["transform"])
    return native, native.build(config["seed"], args["device"], n_genes=len(data.genes))


def execute(output, phase, paths):
    return transigen.execute_phase(output, phase, paths, bind=bind, build_model=build_model, fit=fit)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="phase", required=True)
    registration = commands.add_parser("register")
    registration.add_argument("--axis", choices=TENV_AXIS_STEPS, required=True)
    registration.add_argument("--fold", choices=range(3), type=int, required=True)
    registration.add_argument("--training", type=Path, required=True)
    registration.add_argument("--representations", type=Path, default=STAGE / "util/output/representations")
    registration.add_argument("--source", type=Path, default=REPO / "pert_model/TranSiGen/source")
    registration.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    registration.add_argument("--smoke", action="store_true")
    for phase in ("train", "predict", "replay"):
        command = commands.add_parser(phase)
        for field in ("training", "representations", "source"):
            command.add_argument(f"--{field}", type=Path, help="same pinned input at a new location")
    for child in commands.choices.values():
        child.add_argument("--output", type=Path, required=True)
    args = vars(parser.parse_args())
    phase, output = args.pop("phase"), args.pop("output").resolve()
    if phase == "register":
        options = {key: str(value.resolve()) if isinstance(value, Path) else value for key, value in args.items()}
        print(json.dumps({"registration_sha256": register(output, options)}))
    else:
        paths = {key: str(value.resolve()) for key, value in args.items() if value is not None}
        print(json.dumps(execute(output, phase, paths), sort_keys=True))


if __name__ == "__main__":
    main()
