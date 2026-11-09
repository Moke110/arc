"""Freeze plate memberships, fit a reconstruction model and replay its scores."""

import argparse
import json
from pathlib import Path

import numpy as np

from benchmark import (
    Corpus, DEFAULT_MEAN, DEFAULT_RELEASE, DEFAULT_SPLIT, REPLAY_TOLERANCE, ROLES, SOURCES,
    Split, checked_test_role, fit_mean, freeze_split, json_bytes, load_mean, publish_artifact,
    publish_json, score, verify_artifact, verify_replay,
)


def fit(args, corpus, split):
    mean = fit_mean(corpus, split, args.batch_size)
    if not np.isfinite(mean).all():
        raise ValueError("training-gene mean is nonfinite")
    metadata = {
        "fit_role": "train",
        "fit_profiles": {source: len(split.rows[source]["train"]) for source in SOURCES},
        "batch_size": args.batch_size,
        "reduction": "float64 accumulation, float32 inference",
        "latent_tokens": None,
        "latent_scalars": None,
    }
    with publish_artifact(args.output, corpus, split, "mean", metadata) as directory:
        np.save(directory / "mean.npy", mean, allow_pickle=False)
    return verify_artifact(args.output, corpus, split)


def score_artifact(args, corpus, split):
    checked_test_role(corpus, args.role, args.allow_test)
    record, mean = load_mean(args.artifact, corpus, split)
    result = score(corpus, split, args.role, lambda values: np.broadcast_to(mean, values.shape),
                   mean, args.batch_size, allow_test=args.allow_test)
    return {
        **result,
        "model_id": record["model_id"],
        "artifact_manifest_sha256": record["artifact_manifest_sha256"],
        "reference_manifest_sha256": record["artifact_manifest_sha256"],
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
            command.add_argument("--batch-size", type=int, default=64)
        if name == "fit":
            command.add_argument("--output", type=Path, default=DEFAULT_MEAN)
        if name in {"score", "reload"}:
            command.add_argument("--artifact", type=Path, default=DEFAULT_MEAN)
            command.add_argument("--role", choices=ROLES, default="validation")
            command.add_argument("--allow-test", action="store_true")
        if name == "score":
            command.add_argument("--output", type=Path)
        if name == "reload":
            command.add_argument("--expected", type=Path, required=True)
    args = parser.parse_args()
    try:
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
