"""Seal the selected Step 07 artifacts, then score and replay the held-out comparison."""

import argparse
from hashlib import sha256
from io import BytesIO
import json
from math import isclose
from pathlib import Path

from benchmark import (
    AAE_MODEL_ID, Corpus, DEFAULT_MEAN, DEFAULT_RELEASE, DEFAULT_SPLIT, REPLAY_TOLERANCE,
    SOURCES, STEP, Split, checked_test_role, json_bytes, publish_json, verify_artifact,
    verify_replay,
)
from mlp import MODEL_IDS as MLP_MODEL_IDS
from pca import MODEL_IDS as PCA_MODEL_IDS
from reconstruct import score_artifact


MODEL_IDS = ("mean", *PCA_MODEL_IDS, *MLP_MODEL_IDS, AAE_MODEL_ID)
RECORDS = STEP / "execution"
DEFAULT_SEAL = STEP / "output" / "test_seal.json"


def latent_units(model_id):
    if model_id == "mean":
        return {"tokens": None, "scalars": None}
    if model_id == AAE_MODEL_ID:
        return {"tokens": 512, "scalars": 65_536}
    return {"tokens": None, "scalars": int(model_id.rsplit("_", 1)[1])}


def profile_counts(split, role):
    return {source: len(split.rows[source][role]) for source in SOURCES}


def check_score(report, model_id, role, artifact, reference_digest, split):
    if (report["model_id"] != model_id or report["role"] != role
            or report["provenance"] != split.provenance
            or report["artifact_manifest_sha256"] != artifact["artifact_manifest_sha256"]
            or report["reference_manifest_sha256"] != reference_digest
            or report["replay_tolerance"] != REPLAY_TOLERANCE):
        raise ValueError(f"{model_id}: {role} report does not match the artifact, reference and split")
    sources = report["sources"]
    if set(sources) != set(SOURCES) or {
        source: metrics["profiles"] for source, metrics in sources.items()
    } != profile_counts(split, role):
        raise ValueError(f"{model_id}: incomplete {role} score coverage")
    # Reject nonfinite JSON metrics before a report can be sealed or published.
    json_bytes(report)


def selected_models(models_path, corpus, split, validation_reports):
    artifacts = {model: verify_artifact(models_path / model, corpus, split) for model in MODEL_IDS}
    reference_digest = artifacts["mean"]["artifact_manifest_sha256"]
    selected = {}
    for model_id, artifact in artifacts.items():
        latent = latent_units(model_id)
        metadata = artifact["metadata"]
        if (artifact["model_id"] != model_id
                or metadata["fit_role"] != "train"
                or metadata["fit_profiles"] != profile_counts(split, "train")
                or metadata["latent_tokens"] != latent["tokens"]
                or metadata["latent_scalars"] != latent["scalars"]):
            raise ValueError(f"{model_id}: artifact fit population or latent units differ")
        report = validation_reports[model_id]
        check_score(report, model_id, "validation", artifact, reference_digest, split)
        neural = model_id in (*MLP_MODEL_IDS, AAE_MODEL_ID)
        if neural and (
            metadata["selection_role"] != "validation"
            or metadata["validation_profiles"] != profile_counts(split, "validation")
            or any(not isclose(report["sources"][source]["mse"],
                               metadata["selected_validation_mse"][source],
                               rel_tol=REPLAY_TOLERANCE["rtol"], abs_tol=REPLAY_TOLERANCE["atol"])
                   for source in SOURCES)
        ):
            raise ValueError(f"{model_id}: validation MSE differs from the selected epoch")
        selected[model_id] = {
            "artifact_manifest_sha256": artifact["artifact_manifest_sha256"],
            "files": artifact["files"],
            "latent": latent,
            "selected_epoch": metadata.get("selected_epoch"),
            "device": "cuda" if neural and not corpus.fixture else "cpu",
            "validation": report,
        }
    return {
        "schema_version": 1,
        "provenance": split.provenance,
        "gene_count": corpus.gene_count,
        "test_profiles": profile_counts(split, "test"),
        "replay_tolerance": REPLAY_TOLERANCE,
        "models": selected,
    }


def validation_path(model_id, directory):
    if directory is not None:
        return directory / f"{model_id}-validation.json"
    if model_id == "mean":
        return RECORDS / "mean_record" / "validation.json"
    family = model_id.split("_")[0].lower()
    return RECORDS / f"{family}_record" / f"{model_id}-validation.json"


def compare(args, corpus, split):
    checked_test_role(corpus, "test", args.allow_test)
    seal_bytes = args.seal.read_bytes()
    seal = json.loads(seal_bytes)
    validation_reports = {model: entry["validation"] for model, entry in seal["models"].items()}
    try:
        selected = selected_models(args.models, corpus, split, validation_reports)
    except (ValueError, KeyError) as error:
        raise ValueError(f"selection seal verification failed: {error}") from error
    if selected != seal:
        raise ValueError("selected artifacts differ from the test seal")

    results = {}
    reference_digest = seal["models"]["mean"]["artifact_manifest_sha256"]
    for model_id in MODEL_IDS:
        entry = seal["models"][model_id]
        request = argparse.Namespace(
            artifact=args.models / model_id, mean_artifact=args.models / "mean", role="test",
            allow_test=args.allow_test, batch_size=None, device=entry["device"],
        )
        result = score_artifact(request, corpus, split)
        check_score(result, model_id, "test", entry, reference_digest, split)
        results[model_id] = {**result, "latent": entry["latent"]}
    return {
        "schema_version": 1,
        "seal_sha256": sha256(seal_bytes).hexdigest(),
        "provenance": split.provenance,
        "gene_count": corpus.gene_count,
        "role": "test",
        "replay_tolerance": REPLAY_TOLERANCE,
        "models": results,
    }


def plot_comparison(report, output):
    from matplotlib.figure import Figure

    if report["role"] != "test" or set(report["models"]) != set(MODEL_IDS):
        raise ValueError("figure requires a complete test score report")
    models = report["models"]
    figure = Figure(figsize=(12, 7), layout="constrained")
    figure.suptitle("Whole-transcriptome reconstruction on held-out plates", fontsize=16)
    metrics = (("pcc", "Mean profile PCC"), ("rmse", "Mean profile RMSE"),
               ("relative_squared_error", "Squared error / training-gene mean"))
    colors = ("#57606a",) + ("#00748b",) * 3 + ("#9e4813",) * 3 + ("#6336a3",)
    labels = []
    for model in MODEL_IDS:
        latent = models[model]["latent"]
        if latent["tokens"] is not None:
            units = f"{latent['tokens']:,} tokens · {latent['scalars']:,} scalars"
        elif latent["scalars"] is not None:
            units = f"{latent['scalars']:,} scalars"
        else:
            units = "no latent representation"
        labels.append(f"{model}  ·  {units}")
    for source, name, subfigure in zip(SOURCES, ("Tahoe", "tenv"), figure.subfigures(2, 1)):
        scores = [models[model]["sources"][source] for model in MODEL_IDS]
        subfigure.suptitle(f"{name}  ·  {scores[0]['profiles']:,} profiles", fontsize=12)
        for index, (axis, (field, title)) in enumerate(zip(subfigure.subplots(1, 3), metrics)):
            axis.set_title(title, fontsize=10)
            axis.set_yticks(range(len(labels)), labels if index == 0 else [""] * len(labels))
            axis.set_ylim(len(labels) - 0.35, -0.65)
            axis.tick_params(axis="both", labelsize=9, length=0)
            axis.spines[["top", "right", "left"]].set_visible(False)
            axis.grid(axis="x", color="#dddddd", linewidth=0.7)
            axis.set_axisbelow(True)
            axis.margins(x=0.23)
            if scores[0][field] is not None:
                axis.axvline(scores[0][field], color="#a0a0a0", linestyle=":", linewidth=1)
            for row, (score, color) in enumerate(zip(scores, colors)):
                value = score[field]
                if value is None:
                    axis.text(0.05, row, "undefined", transform=axis.get_yaxis_transform(),
                              va="center", fontsize=8)
                    continue
                axis.scatter(value, row, color=color, s=25, zorder=3)
                axis.annotate(f"{value:.4f}", (value, row), xytext=(5, 0),
                              textcoords="offset points", va="center", fontsize=8)
    save_figure(figure, output)


def plot_aae_loss(history, output):
    from matplotlib.figure import Figure
    from matplotlib.ticker import MaxNLocator

    x_field = "optimizer_steps" if "optimizer_steps" in history[0] else "epoch"
    steps = [record[x_field] for record in history]
    tenv = [record["validation_mse"]["tenv"] for record in history]
    best_index = tenv.index(min(tenv))
    figure = Figure(figsize=(8, 4.5), layout="constrained")
    axis = figure.subplots()
    axis.set_title("AAE_512 training and validation loss")
    axis.set_xlabel("Optimizer updates" if x_field == "optimizer_steps" else "Epoch")
    axis.set_ylabel("Full-transcriptome MSE")
    axis.set_yscale("log")
    axis.xaxis.set_major_locator(MaxNLocator(integer=True))
    axis.grid(alpha=0.25)
    axis.plot(steps, [record["training_mse"] for record in history],
              color="#57606a", marker=".", label="Joint training (batch losses)")
    axis.plot(steps, [record["validation_mse"]["tahoe"] for record in history],
              color="#00748b", marker=".", label="Tahoe validation")
    axis.plot(steps, tenv, color="#6336a3", marker=".", label="tenv validation")
    axis.axvline(steps[best_index], color="#6336a3", linestyle=":", alpha=0.6,
                 label=f"Selected checkpoint ({steps[best_index]:,})")
    axis.legend()
    save_figure(figure, output)


def save_figure(figure, output):
    image = BytesIO()
    figure.savefig(image, format="png", dpi=160)
    with output.open("xb") as stream:
        stream.write(image.getbuffer())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("seal", "score", "reload"):
        command = commands.add_parser(name)
        command.add_argument("--release", type=Path, default=DEFAULT_RELEASE)
        command.add_argument("--split", type=Path, default=DEFAULT_SPLIT)
        command.add_argument("--models", type=Path, default=DEFAULT_MEAN.parent)
        command.add_argument("--seal", type=Path, default=DEFAULT_SEAL)
        command.add_argument("--fixture", action="store_true")
        if name == "seal":
            command.add_argument("--validation-dir", type=Path,
                                 help="directory of <model>-validation.json reports")
        else:
            command.add_argument("--allow-test", action="store_true")
            if name == "score":
                command.add_argument("--output", type=Path, required=True)
            else:
                command.add_argument("--expected", type=Path, required=True)
    figure = commands.add_parser("figure")
    figure.add_argument("--expected", type=Path, required=True, help="published test score report")
    figure.add_argument("--output", type=Path, required=True)
    loss = commands.add_parser("aae-loss-figure")
    loss.add_argument("--history", type=Path, required=True, help="AAE epochs.json")
    loss.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "figure":
            plot_comparison(json.loads(args.expected.read_text()), args.output)
            return
        if args.command == "aae-loss-figure":
            plot_aae_loss(json.loads(args.history.read_text()), args.output)
            return
        corpus = Corpus(args.release, args.fixture)
        split = Split(corpus, args.split)
        if args.command == "seal":
            reports = {model: json.loads(validation_path(model, args.validation_dir).read_text())
                       for model in MODEL_IDS}
            result = selected_models(args.models, corpus, split, reports)
            publish_json(args.seal, result)
        elif args.command == "score":
            if args.output.exists():
                raise ValueError(f"score report already exists: {args.output}")
            result = compare(args, corpus, split)
            publish_json(args.output, result)
        else:
            expected = json.loads(args.expected.read_text())
            result = compare(args, corpus, split)
            verify_replay(result, expected)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"error: {error}\n")
    print(json_bytes(result).decode(), end="")


if __name__ == "__main__":
    main()
