"""The deterministic ``Mean`` reference for the tenv comparison (#165).

ADR 0063 keeps the retired ``delta_mean`` reference under the name ``Mean`` with
its definition unchanged.  The retired implementation
(``git show 1543781^:stage2/03_modes_tenv/execution/script/baselines.py``)
computed **one global prediction vector**: the equal-weight mean over the
training-partition drug-dose groups of each group's mean measured Δ, and copied
that same vector to every test well.  It did not look up the test well's own
group.  A cold-drug test group is therefore absent from the training partition by
construction and still receives the same global reference — that is the recorded
semantics, not a fallback.  ADR 0063 words the reference as "the mean Δ of the
same group in the training partition"; the retired code is the authority for the
exact arithmetic and this module records the interpretation it implements.

This module is presentation-free and deterministic: same corpus, fold and inputs
give the same bytes.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
STEP_ROOT = SCRIPT_DIR.parent.parent
EXPECTED_GENES = 18164
SCHEMA = "stage2-tenv-mean-reference-v1"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def group_mean_prediction(train_rows, train_delta):
    """Equal-weight mean over training drug-dose groups of each group's mean Δ.

    Byte-compatible with the retired ``delta_mean``: groups sorted by ``group_id``,
    float64 accumulation, then a float32 prediction vector.
    """
    train_delta = np.asarray(train_delta, dtype=np.float64)
    if train_delta.ndim != 2 or len(train_delta) != len(train_rows) or not len(train_rows):
        raise ValueError("training responses must be one nonempty row per training well")
    groups = sorted({str(row["group_id"]) for row in train_rows})
    if not groups:
        raise ValueError("Mean requires at least one training-partition group")
    delta = np.zeros(train_delta.shape[1], dtype=np.float64)
    for group in groups:
        members = [index for index, row in enumerate(train_rows)
                   if str(row["group_id"]) == group]
        delta += np.mean(train_delta[members], axis=0, dtype=np.float64) / len(groups)
    return delta.astype(np.float32)


def publish(directory, *, fold, seed, case_ids, prediction, source):
    """Write a complete ``Mean`` run directory the shared driver validates."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    case_ids = [str(value) for value in case_ids]
    prediction = np.ascontiguousarray(prediction, dtype=np.float32)
    if prediction.shape != (len(case_ids), EXPECTED_GENES):
        raise ValueError(f"Mean prediction shape {prediction.shape} is not "
                         f"{(len(case_ids), EXPECTED_GENES)}")
    if not np.isfinite(prediction).all():
        raise ValueError("Mean prediction must be finite")
    np.save(directory / "prediction.npy", prediction)
    write_json(directory / "case_ids.json", case_ids)
    identity = {
        "schema": SCHEMA,
        "identity": {
            "model": "Mean",
            "deterministic_reference": True,
            "fold": int(fold),
            "seed": int(seed),
            "formula": ("equal-weight mean over training-partition drug-dose groups of "
                        "each group's mean measured delta; one global vector for every "
                        "test well"),
            "cold_drug_semantics": ("a test group absent from the training partition "
                                    "receives the same global reference; no per-group "
                                    "lookup and no fallback exist"),
            "retired_implementation": "stage2/03_modes_tenv/execution/script/baselines.py@1543781^",
        },
        "source": source,
        "prediction_sha256": sha256(directory / "prediction.npy"),
        "case_ids_sha256": sha256(directory / "case_ids.json"),
        "trainable_parameters": 0,
        "group_weighting": "equal training drug-dose groups",
    }
    write_json(directory / "identity.json", identity)
    write_json(directory / "checkpoint.json", {
        "history": [],
        "deterministic_reference": True,
        "selection": "no fitted checkpoint",
        "training_seconds": 0.0,
        "identity": identity["identity"],
    })
    write_json(directory / "terminal.json", {
        "schema": SCHEMA,
        "status": "complete",
        "phase": "complete",
        "prediction": "prediction.npy",
        "training_seconds": 0.0,
        "identity": identity["identity"],
    })
    return identity


def build(*, corpus, fold_product, fold_row, seed_row, output, declared_folds=2,
          declared_seeds=2):
    """Build one fold's Mean run from the accepted corpus and fold product."""
    sys.path.insert(0, str(SCRIPT_DIR))
    import pandas as pd
    from evaluation_driver import EXPECTED_GENES as GENES, fold_cases, read_layer_rows
    readers = public_readers()
    index = readers["fold_seed"].FoldSeedIndex.load(fold_product)
    folds = index.select_folds(declared_folds)
    seeds = index.select_seeds(declared_seeds)
    fold = int(folds.iloc[fold_row]["fold"])
    seed = int(seeds.iloc[seed_row]["seed"])
    cohorts = pd.read_csv(Path(fold_product) / "cohorts.csv", keep_default_na=False,
                          dtype={"sample_id": str, "group_id": str, "role": str})
    part = cohorts[cohorts.fold.astype(int).eq(fold)]
    train = part[part.role.eq("train")].sort_values("sample_id")
    test = fold_cases(cohorts, [fold])[fold]
    train_rows = [{"sample_id": str(row.sample_id), "group_id": str(row.group_id)}
                  for row in train.itertuples()]
    obs, var = readers["corpora"].training_header(corpus)
    if len(var) != GENES:
        raise ValueError(f"corpus gene axis is {len(var)}, not {GENES}")
    position = {str(value): row for row, value in enumerate(obs["sample_id"])}
    missing = [row["sample_id"] for row in train_rows if row["sample_id"] not in position]
    if missing:
        raise ValueError(f"{len(missing)} training well(s) missing from the corpus")
    rows = np.asarray([position[row["sample_id"]] for row in train_rows], dtype=np.int64)
    train_delta = read_layer_rows(corpus, "layers/delta", rows, GENES)
    delta = group_mean_prediction(train_rows, train_delta)
    prediction = np.broadcast_to(delta, (len(test), GENES))
    return publish(output, fold=fold, seed=seed,
                   case_ids=[row["sample_id"] for row in test],
                   prediction=prediction,
                   source={"corpus_sha256": sha256(corpus),
                           "cohorts_sha256": sha256(Path(fold_product) / "cohorts.csv"),
                           "folds_sha256": sha256(Path(fold_product) / "folds.csv"),
                           "seeds_sha256": sha256(Path(fold_product) / "seeds.csv"),
                           "fold_row": int(fold_row), "seed_row": int(seed_row),
                           "declared": {"folds": declared_folds, "seeds": declared_seeds}})


def public_readers():
    """Bind the accepted public implementation under ``01_bases/output/api``."""
    root = STEP_ROOT.parents[1]  # repo root: stage2/02_models_tenv -> repo
    api = root / "stage2/01_bases/output/api"
    if not api.is_dir():
        raise RuntimeError(f"public implementation not found: {api}")
    if str(api) not in sys.path:
        sys.path.insert(0, str(api))
    from bases import corpora, fold_seed
    return {"corpora": corpora, "fold_seed": fold_seed}


def select_identities(fold_product, fold_row, seed_row, declared_folds=2, declared_seeds=2):
    """Declared identity for one matrix cell, read through the fold/seed index."""
    readers = public_readers()
    index = readers["fold_seed"].FoldSeedIndex.load(fold_product)
    folds = index.select_folds(declared_folds)
    seeds = index.select_seeds(declared_seeds)
    return {"fold_row": int(fold_row), "fold": int(folds.iloc[fold_row]["fold"]),
            "seed_row": int(seed_row), "seed": int(seeds.iloc[seed_row]["seed"])}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    root = STEP_ROOT.parents[1]
    parser.add_argument("--corpus", default=str(
        root / "stage2/01_bases/output/training-corpora/tenv.h5ad"))
    parser.add_argument("--fold-product", default=str(
        root / "stage2/01_bases/output/fold_n_seed"))
    parser.add_argument("--output", required=True)
    parser.add_argument("--fold-row", type=int, required=True)
    parser.add_argument("--seed-row", type=int, required=True)
    parser.add_argument("--folds", type=int, default=2)
    parser.add_argument("--seeds", type=int, default=2)
    args = parser.parse_args(argv)
    identity = build(corpus=args.corpus, fold_product=args.fold_product,
                     fold_row=args.fold_row, seed_row=args.seed_row,
                     output=args.output, declared_folds=args.folds,
                     declared_seeds=args.seeds)
    print(json.dumps({"identity": identity["identity"],
                      "prediction_sha256": identity["prediction_sha256"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
