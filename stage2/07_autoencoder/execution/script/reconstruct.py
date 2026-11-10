"""Freeze plate memberships, fit a reconstruction model and replay its scores."""

import argparse
import json
from pathlib import Path

import numpy as np

from benchmark import (
    AAE_MODEL_IDS, AAE_TRANSFORMER_MODEL_IDS, Corpus, DEFAULT_MEAN, DEFAULT_RELEASE, DEFAULT_SPLIT,
    MATRIX_AE_MODEL_ID,
    REPLAY_TOLERANCE, ROLES, SOURCES, Split, checked_test_role, fit_mean, freeze_split, json_bytes,
    load_mean, publish_artifact, publish_json, score, validate_batch_size, validation_mse_matches,
    verify_artifact, verify_replay,
)
from mlp import MODEL_IDS as MLP_MODEL_IDS, fit_mlp, mlp_predictor
from pca import MODEL_IDS as PCA_MODEL_IDS, fit_pca, pca_predictor


def fit(args, corpus, split):
    if args.model in AAE_MODEL_IDS:
        from aae import fit_aae

        if args.output is None:
            raise ValueError("AAE fit requires --output for the new artifact")
        return fit_aae(corpus, split, args.output, args.batch_size, args.device, args.model)
    if args.model == MATRIX_AE_MODEL_ID:
        from matrix_ae import fit_matrix_ae

        if args.output is None:
            raise ValueError(f"{MATRIX_AE_MODEL_ID} fit requires --output for the new artifact")
        return fit_matrix_ae(corpus, split, args.output, args.batch_size, args.device)
    if args.model in MLP_MODEL_IDS:
        if args.output is None:
            raise ValueError("MLP fit requires --output for the new artifact")
        return fit_mlp(corpus, split, args.output, args.model, args.batch_size, args.device)
    batch_size = 64 if args.batch_size is None else args.batch_size
    if args.model == "PCA":
        if args.output is None:
            raise ValueError("PCA fit requires --output for the new model set")
        return fit_pca(corpus, split, args.output, batch_size, args.device)
    output = args.output or DEFAULT_MEAN
    mean = fit_mean(corpus, split, batch_size)
    if not np.isfinite(mean).all():
        raise ValueError("training-gene mean is nonfinite")
    metadata = {
        "fit_role": "train",
        "fit_profiles": {source: len(split.rows[source]["train"]) for source in SOURCES},
        "batch_size": batch_size,
        "reduction": "float64 accumulation, float32 inference",
        "latent_tokens": None,
        "latent_scalars": None,
    }
    with publish_artifact(output, corpus, split, "mean", metadata) as directory:
        np.save(directory / "mean.npy", mean, allow_pickle=False)
    return verify_artifact(output, corpus, split)


def score_artifact(args, corpus, split):
    checked_test_role(corpus, args.role, args.allow_test)
    record = verify_artifact(args.artifact, corpus, split)
    batch_size = 64 if args.batch_size is None else args.batch_size
    if record["model_id"] == "mean":
        reference, mean = load_mean(args.artifact, corpus, split)
        predict = lambda values: np.broadcast_to(mean, values.shape)
    else:
        reference, mean = load_mean(args.mean_artifact, corpus, split)
        if record["model_id"] in PCA_MODEL_IDS:
            predict = pca_predictor(args.artifact, record)
        elif record["model_id"] in MLP_MODEL_IDS:
            predict = mlp_predictor(args.artifact, record, args.device)
        elif record["model_id"] in AAE_MODEL_IDS:
            from aae import aae_predictor

            predict, batch_size = aae_predictor(
                args.artifact, record, corpus.fixture, args.batch_size, args.device,
            )
        elif record["model_id"] == MATRIX_AE_MODEL_ID:
            from matrix_ae import matrix_ae_predictor

            predict, batch_size = matrix_ae_predictor(
                args.artifact, record, corpus.fixture, args.batch_size, args.device,
            )
        else:
            raise ValueError(f"unsupported model: {record['model_id']}")
    result = score(corpus, split, args.role, predict, mean, batch_size, allow_test=args.allow_test)
    if record["model_id"] in AAE_TRANSFORMER_MODEL_IDS and args.role == "validation":
        selected_mse = record["metadata"]["selected_validation_mse"]
        if not validation_mse_matches(result["sources"], selected_mse):
            raise ValueError("validation MSE differs from the selected epoch")
    return {
        **result,
        "model_id": record["model_id"],
        "artifact_manifest_sha256": record["artifact_manifest_sha256"],
        "reference_manifest_sha256": reference["artifact_manifest_sha256"],
        "provenance": record["provenance"],
        "replay_tolerance": REPLAY_TOLERANCE,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("freeze-split", "fit", "score", "reload"):
        command = commands.add_parser(name)
        command.add_argument("--release", type=Path, default=DEFAULT_RELEASE)
        command.add_argument("--split", type=Path, default=DEFAULT_SPLIT)
        command.add_argument("--fixture", action="store_true")
        if name != "freeze-split":
            command.add_argument("--batch-size", type=int)
            command.add_argument("--device", choices=("cpu", "cuda"), default="cpu")
        if name == "fit":
            command.add_argument(
                "--model", choices=("mean", "PCA", *MLP_MODEL_IDS, *AAE_MODEL_IDS, MATRIX_AE_MODEL_ID),
                default="mean",
            )
            command.add_argument("--output", type=Path)
        if name in {"score", "reload"}:
            command.add_argument("--artifact", type=Path, default=DEFAULT_MEAN)
            command.add_argument("--mean-artifact", type=Path, default=DEFAULT_MEAN)
            command.add_argument("--role", choices=ROLES, default="validation")
            command.add_argument("--allow-test", action="store_true", help="fixture test scoring only")
        if name == "score":
            command.add_argument("--output", type=Path)
        if name == "reload":
            command.add_argument("--expected", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command in {"score", "reload"} and args.role == "test" and not args.fixture:
            raise ValueError("real test scoring requires a selection seal through compare.py")
        if args.command != "freeze-split" and args.batch_size is not None:
            validate_batch_size(args.batch_size)
        corpus = Corpus(args.release, args.fixture)
        if args.command == "freeze-split":
            result = freeze_split(corpus, args.split)
        else:
            split = Split(corpus, args.split)
            if args.command == "fit":
                result = fit(args, corpus, split)
            else:
                result = score_artifact(args, corpus, split)
                if args.command == "reload":
                    verify_replay(result, json.loads(args.expected.read_text()))
                elif args.output:
                    publish_json(args.output, result)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"error: {error}\n")
    print(json_bytes(result).decode(), end="")


if __name__ == "__main__":
    main()
