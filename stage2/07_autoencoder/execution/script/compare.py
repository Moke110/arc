"""Seal the selected Step 07 artifacts, then score and replay the held-out comparison."""

import argparse
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import subprocess

from benchmark import (
    AAE_MODEL_ID, AAE_MODEL_IDS, AAE_RECONSTRUCTION_TOKENS,
    AAE_TRANSFORMER_MODEL_ID, AAE_TRANSFORMER_MODEL_IDS, AAE_WIDTH, Corpus, DEFAULT_MEAN,
    DEFAULT_RELEASE, DEFAULT_SPLIT, MATRIX_AE_CHANNELS, MATRIX_AE_MODEL_ID, MATRIX_AE_TOKENS,
    INTERFACE_ID, MATRIX_MLP_MODEL_ID, SHARED_NEURAL_IDS, SHARED_PCA_ID,
    REPLAY_TOLERANCE, ROOT, SOURCES, STEP, Split, aae_encoder_tokens, check_checkpoint_selection,
    checked_test_role, json_bytes, publish_json, shared_extra_tokens, shared_matrix_shape,
    validation_mse_matches, verify_artifact, verify_replay,
)
from mlp import MODEL_IDS as MLP_MODEL_IDS
from pca import MODEL_IDS as PCA_MODEL_IDS
from reconstruct import score_artifact


LEGACY_MODEL_IDS = ("mean", *PCA_MODEL_IDS, *MLP_MODEL_IDS, AAE_MODEL_ID)
MODEL_IDS = (*LEGACY_MODEL_IDS, MATRIX_AE_MODEL_ID)
EXPANDED_MODEL_IDS = (*LEGACY_MODEL_IDS, *AAE_TRANSFORMER_MODEL_IDS, MATRIX_AE_MODEL_ID)
MATRIX_MLP_COMPARISON_IDS = (
    "mean", "PCA_512", "MLP_512", MATRIX_MLP_MODEL_ID, *AAE_MODEL_IDS, MATRIX_AE_MODEL_ID,
)
SHARED_COMPARISON_IDS = ("mean", SHARED_PCA_ID, *SHARED_NEURAL_IDS)
SHARED_SCHEMA = 4
COMPARISON_SCHEMAS = {
    LEGACY_MODEL_IDS: 1, MODEL_IDS: 1, EXPANDED_MODEL_IDS: 2, MATRIX_MLP_COMPARISON_IDS: 3,
    SHARED_COMPARISON_IDS: SHARED_SCHEMA,
}
NEURAL_MODEL_IDS = (*MLP_MODEL_IDS, *AAE_MODEL_IDS, MATRIX_AE_MODEL_ID, MATRIX_MLP_MODEL_ID)
EXPANDED_DISPLAY_NAMES = {AAE_MODEL_ID: "AAE_512_SA", AAE_TRANSFORMER_MODEL_ID: "AAE_512"}
MATRIX_MLP_DISPLAY_NAMES = {**EXPANDED_DISPLAY_NAMES, "MLP_512": "MLPv_512"}
RECORDS = STEP / "execution"
DEFAULT_SEAL = STEP / "output" / "test_seal.json"
EXPANDED_SEAL = STEP / "output" / "test_seal_aae_ablation.json"
MATRIX_AE_COLOR = "#287a42"


def comparison_model_ids(models):
    for model_ids in COMPARISON_SCHEMAS:
        if models.keys() == set(model_ids):
            return model_ids
    raise ValueError("unexpected comparison model set")


def display_name(model_id, schema):
    names = MATRIX_MLP_DISPLAY_NAMES if schema == 3 else EXPANDED_DISPLAY_NAMES
    return names.get(model_id, model_id)


def latent_units(model_id, show_decoder):
    if model_id == "mean":
        return {"tokens": None, "scalars": None}
    if model_id in AAE_MODEL_IDS:
        tokens = aae_encoder_tokens(model_id)
        latent = {"tokens": tokens, "scalars": tokens * AAE_WIDTH}
        if show_decoder:
            latent.update(decoder_tokens=AAE_RECONSTRUCTION_TOKENS,
                          decoder_scalars=AAE_RECONSTRUCTION_TOKENS * AAE_WIDTH)
        return latent
    if model_id == MATRIX_AE_MODEL_ID:
        return {"tokens": MATRIX_AE_TOKENS, "scalars": MATRIX_AE_TOKENS * MATRIX_AE_CHANNELS}
    if model_id == MATRIX_MLP_MODEL_ID:
        return {"tokens": AAE_RECONSTRUCTION_TOKENS, "scalars": AAE_RECONSTRUCTION_TOKENS * AAE_WIDTH}
    return {"tokens": None, "scalars": int(model_id.rsplit("_", 1)[1])}


def shared_latent_units(model_id, gene_count, fixture):
    if model_id == "mean":
        return latent_units(model_id, False)
    tokens, channels = shared_matrix_shape(gene_count, fixture)
    extra_tokens = shared_extra_tokens(model_id, fixture)
    return {"tokens": tokens, "scalars": tokens * channels, "channels": channels,
            "encoder_tokens": tokens + extra_tokens}


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


def check_matrix_checkpoint(path, artifact, device):
    metadata = artifact["metadata"]
    if (set(artifact["files"]) != {"weights.pt", "epochs.json"}
            or metadata["runtime"]["device"] != device):
        raise ValueError(f"{MATRIX_AE_MODEL_ID}: checkpoint files or fit device differ")
    history = json.loads((path / "epochs.json").read_text())
    check_checkpoint_selection(history, metadata, MATRIX_AE_MODEL_ID)


def selected_models(models_path, corpus, split, validation_reports, model_ids):
    artifacts = {model: verify_artifact(models_path / model, corpus, split) for model in model_ids}
    reference_digest = artifacts["mean"]["artifact_manifest_sha256"]
    schema = COMPARISON_SCHEMAS[model_ids]
    selected = {}
    for model_id, artifact in artifacts.items():
        latent = latent_units(model_id, schema >= 2)
        metadata = artifact["metadata"]
        if (artifact["model_id"] != model_id
                or metadata["fit_role"] != "train"
                or metadata["fit_profiles"] != profile_counts(split, "train")
                or metadata["latent_tokens"] != latent["tokens"]
                or metadata["latent_scalars"] != latent["scalars"]):
            raise ValueError(f"{model_id}: artifact fit population or latent units differ")
        report = validation_reports[model_id]
        check_score(report, model_id, "validation", artifact, reference_digest, split)
        neural = model_id in NEURAL_MODEL_IDS
        if neural and (
            metadata["selection_role"] != "validation"
            or metadata["validation_profiles"] != profile_counts(split, "validation")
            or not validation_mse_matches(report["sources"], metadata["selected_validation_mse"])
        ):
            raise ValueError(f"{model_id}: validation MSE differs from the selected epoch")
        device = "cuda" if neural and not corpus.fixture else "cpu"
        if model_id == MATRIX_AE_MODEL_ID:
            check_matrix_checkpoint(models_path / model_id, artifact, device)
        elif model_id == MATRIX_MLP_MODEL_ID:
            from matrix_mlp import check_matrix_mlp_artifact

            check_matrix_mlp_artifact(models_path / model_id, artifact, corpus.fixture, device)
        elif model_id in AAE_TRANSFORMER_MODEL_IDS:
            from aae import check_aae_artifact

            check_aae_artifact(models_path / model_id, artifact, corpus.fixture, device)
        selected[model_id] = {
            "artifact_manifest_sha256": artifact["artifact_manifest_sha256"],
            "files": artifact["files"],
            "latent": latent,
            "selected_epoch": metadata.get("selected_epoch"),
            "device": device,
            "validation": report,
        }
        if schema >= 2:
            selected[model_id]["display_name"] = display_name(model_id, schema)
    return {
        "schema_version": schema,
        "provenance": split.provenance,
        "gene_count": corpus.gene_count,
        "test_profiles": profile_counts(split, "test"),
        "replay_tolerance": REPLAY_TOLERANCE,
        "models": selected,
    }


def selected_shared_models(models_path, interface_path, corpus, split, validation_reports):
    from shared_neural import check_shared_neural_artifact

    artifacts = {model: verify_artifact(models_path / model, corpus, split)
                 for model in SHARED_COMPARISON_IDS}
    interface = verify_artifact(interface_path, corpus, split)
    _, channels = shared_matrix_shape(corpus.gene_count, corpus.fixture)
    if (set(interface["files"]) != {"weights.pt"}
            or interface["metadata"]["generator"] != [1, 64, channels]
            or interface["metadata"]["restorer"] != [channels, 64, 1]):
        raise ValueError("invalid shared expression interface architecture")
    interface_digest = interface["artifact_manifest_sha256"]
    reference_digest = artifacts["mean"]["artifact_manifest_sha256"]
    train_counts = profile_counts(split, "train")
    validation_counts = profile_counts(split, "validation")
    for model_id, artifact in {**artifacts, INTERFACE_ID: interface}.items():
        metadata = artifact["metadata"]
        if (artifact["model_id"] != model_id or metadata["fit_role"] != "train"
                or metadata["fit_profiles"] != train_counts):
            raise ValueError(f"{model_id}: artifact fit population differs from the training split")
        report = validation_reports[model_id]
        check_score(report, model_id, "validation", artifact, reference_digest, split)
        if model_id != "mean" and report["interface_manifest_sha256"] != interface_digest:
            raise ValueError(f"{model_id}: validation uses a different shared expression interface")

    device = "cpu" if corpus.fixture else "cuda"
    selected = {}
    common_training = None
    for model_id, artifact in artifacts.items():
        metadata = artifact["metadata"]
        latent = shared_latent_units(model_id, corpus.gene_count, corpus.fixture)
        entry = {
            "artifact_manifest_sha256": artifact["artifact_manifest_sha256"],
            "files": artifact["files"],
            "latent": latent,
            "trainable_parameters": 0,
            "device": "cpu" if model_id == "mean" else device,
            "validation": validation_reports[model_id],
            "display_name": model_id,
        }
        if model_id != "mean":
            if metadata["interface_manifest_sha256"] != interface_digest:
                raise ValueError(f"{model_id}: artifact uses a different shared expression interface")
            if metadata["latent_shape"] != [latent["tokens"], latent["channels"]]:
                raise ValueError(f"{model_id}: reconstruction token shape differs")
            entry["trainable_parameters"] = metadata["trainable_parameters"]
        if model_id == SHARED_PCA_ID:
            entry["fitted_parameters"] = metadata["fitted_parameters"]
        if model_id in SHARED_NEURAL_IDS:
            if (metadata["encoder_shape"] != [latent["encoder_tokens"], latent["channels"]]
                    or metadata["validation_profiles"] != validation_counts
                    or not validation_mse_matches(entry["validation"]["sources"],
                                                  metadata["selected_validation_mse"])):
                raise ValueError(f"{model_id}: encoder shape or selected validation evidence differs")
            path = models_path / model_id
            check_shared_neural_artifact(path, artifact, corpus.fixture, device)
            history = json.loads((path / "epochs.json").read_text())
            training = {
                "budget_sha256": metadata["budget_sha256"],
                "recipe": metadata["recipe"],
                "sample_order": [
                    {key: check[key] for key in ("epoch", "training_epoch", "batch_offset",
                                                 "optimizer_steps", "sample_order_sha256")}
                    for check in history
                ],
            }
            if common_training is not None and common_training != training:
                raise ValueError(f"{model_id}: neural recipe, update budget or sample order differs")
            common_training = training
            entry["selected_epoch"] = metadata["selected_epoch"]
            entry["selected_optimizer_steps"] = history[metadata["selected_epoch"] - 1]["optimizer_steps"]
        selected[model_id] = entry

    return {
        "schema_version": SHARED_SCHEMA,
        "provenance": split.provenance,
        "gene_count": corpus.gene_count,
        "fixture": corpus.fixture,
        "test_profiles": profile_counts(split, "test"),
        "reuses_published_test_plates": not corpus.fixture,
        "replay_tolerance": REPLAY_TOLERANCE,
        INTERFACE_ID: {
            "artifact_manifest_sha256": interface_digest,
            "files": interface["files"],
            "representation_shape": [corpus.gene_count, channels],
            "trainable_parameters": interface["metadata"]["trainable_parameters"],
            "validation": validation_reports[INTERFACE_ID],
        },
        "common_neural_training": common_training,
        "models": selected,
    }


def validation_path(model_id, directory):
    if directory is not None:
        return directory / f"{model_id}-validation.json"
    if model_id == "mean":
        return RECORDS / "mean_record" / "validation.json"
    if model_id in AAE_TRANSFORMER_MODEL_IDS:
        return RECORDS / "aae_ablation_record" / f"{model_id}-validation.json"
    family = model_id.split("_")[0].lower()
    return RECORDS / f"{family}_record" / f"{model_id}-validation.json"


def read_test_seal(args):
    seal_bytes = args.seal.read_bytes()
    seal = json.loads(seal_bytes)
    if not args.fixture and (seal.get("schema_version") == SHARED_SCHEMA
                             or seal["models"].keys() == set(SHARED_COMPARISON_IDS)):
        path = args.seal.resolve()
        if not path.is_relative_to(ROOT):
            raise ValueError("shared expression seal must match a file committed in HEAD")
        committed = subprocess.run(
            ["git", "-C", str(ROOT), "cat-file", "blob", f"HEAD:{path.relative_to(ROOT).as_posix()}"],
            capture_output=True, check=False,
        )
        if committed.returncode or committed.stdout != seal_bytes:
            raise ValueError("shared expression seal must match a file committed in HEAD")
    return seal_bytes, seal


def compare(args):
    seal_bytes, seal = read_test_seal(args)
    corpus = Corpus(args.release, args.fixture)
    split = Split(corpus, args.split)
    checked_test_role(corpus, "test", args.allow_test)
    try:
        model_ids = comparison_model_ids(seal["models"])
        validation_reports = {model: seal["models"][model]["validation"] for model in model_ids}
        shared = model_ids == SHARED_COMPARISON_IDS
        if shared:
            validation_reports[INTERFACE_ID] = seal[INTERFACE_ID]["validation"]
            selected = selected_shared_models(args.models, args.interface, corpus, split,
                                              validation_reports)
        else:
            selected = selected_models(args.models, corpus, split, validation_reports, model_ids)
    except (ValueError, KeyError) as error:
        raise ValueError(f"selection seal verification failed: {error}") from error
    if selected != seal:
        raise ValueError("selected artifacts differ from the test seal")

    results = {}
    reference_digest = seal["models"]["mean"]["artifact_manifest_sha256"]
    for model_id in model_ids:
        entry = seal["models"][model_id]
        request = argparse.Namespace(
            artifact=args.models / model_id, mean_artifact=args.models / "mean", role="test",
            allow_test=args.allow_test, batch_size=None, device=entry["device"],
            interface=args.interface,
        )
        result = score_artifact(request, corpus, split)
        check_score(result, model_id, "test", entry, reference_digest, split)
        if shared and model_id != "mean" and (
            result["interface_manifest_sha256"] != seal[INTERFACE_ID]["artifact_manifest_sha256"]
        ):
            raise ValueError(f"{model_id}: test score uses a different shared expression interface")
        results[model_id] = {**result, "latent": entry["latent"]}
        if seal["schema_version"] >= 2:
            results[model_id]["display_name"] = entry["display_name"]
        if shared:
            for key in ("trainable_parameters", "fitted_parameters", "selected_epoch",
                        "selected_optimizer_steps"):
                if key in entry:
                    results[model_id][key] = entry[key]
    report = {
        "schema_version": seal["schema_version"],
        "seal_sha256": sha256(seal_bytes).hexdigest(),
        "provenance": split.provenance,
        "gene_count": corpus.gene_count,
        "role": "test",
        "replay_tolerance": REPLAY_TOLERANCE,
        "models": results,
    }
    if shared:
        for key in ("fixture", "reuses_published_test_plates", INTERFACE_ID, "common_neural_training"):
            report[key] = seal[key]
    return report


def plot_comparison(report, output):
    from matplotlib.figure import Figure

    if report["role"] != "test":
        raise ValueError("figure requires a complete test score report")
    models = report["models"]
    model_ids = comparison_model_ids(models)
    schema = COMPARISON_SCHEMAS[model_ids]
    if report.get("schema_version", 1) != schema:
        raise ValueError("test score report schema version does not match its model set")
    figure = Figure(figsize=(12, 9 if schema == 2 else 7), layout="constrained")
    figure.suptitle("Whole-transcriptome reconstruction on held-out plates", fontsize=16)
    metrics = (("pcc", "Mean profile PCC"), ("rmse", "Mean profile RMSE"),
               ("relative_squared_error", "Squared error / training-gene mean"))
    colors = ("#57606a",) + ("#00748b",) * 3 + ("#9e4813",) * 3 + ("#6336a3",)
    if schema == 2:
        colors += ("#a23d86", "#4359a7")
    colors += (MATRIX_AE_COLOR,)
    if schema == 3:
        colors = ("#57606a", "#00748b", "#9e4813", "#bd6922", "#6336a3", "#a23d86",
                  "#4359a7", MATRIX_AE_COLOR)
    if schema == SHARED_SCHEMA:
        colors = ("#57606a", "#00748b", "#9e4813", "#6336a3", "#a23d86", "#4359a7")
    labels = []
    for model in model_ids:
        latent = models[model]["latent"]
        expected_latent = (shared_latent_units(model, report["gene_count"],
                                              report.get("fixture", False))
                           if schema == SHARED_SCHEMA else latent_units(model, schema >= 2))
        if latent != expected_latent:
            raise ValueError(f"{model}: unexpected latent units")
        if schema == SHARED_SCHEMA and latent["tokens"] is not None:
            units = f"{latent['tokens']:,} × {latent['channels']}"
            if latent["encoder_tokens"] > latent["tokens"]:
                units += f" ({latent['encoder_tokens']:,} encoder)"
        elif latent["tokens"] is not None:
            units = f"{latent['tokens']:,} tokens · {latent['scalars']:,} scalars"
        elif latent["scalars"] is not None:
            units = f"{latent['scalars']:,} scalars"
        else:
            units = "no latent representation"
        name = model
        if schema >= 2:
            name = display_name(model, schema)
            if models[model]["display_name"] != name:
                raise ValueError(f"{model}: unexpected display name")
        labels.append(f"{name}  ·  {units}")
    for source, name, subfigure in zip(SOURCES, ("Tahoe", "tenv"), figure.subfigures(2, 1)):
        scores = [models[model]["sources"][source] for model in model_ids]
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
            for row, (model, score, color) in enumerate(zip(model_ids, scores, colors)):
                value = score[field]
                if value is None:
                    axis.text(0.05, row, "undefined", transform=axis.get_yaxis_transform(),
                              va="center", fontsize=8)
                    continue
                marker = "D" if model == MATRIX_AE_MODEL_ID else "o"
                axis.scatter(value, row, color=color, marker=marker, s=25, zorder=3)
                axis.annotate(f"{value:.4f}", (value, row), xytext=(5, 0),
                              textcoords="offset points", va="center", fontsize=8)
    save_figure(figure, output)


def plot_attention_loss(history, output, model_id):
    from matplotlib.figure import Figure
    from matplotlib.ticker import MaxNLocator

    if not history:
        raise ValueError("loss figure requires validation history")
    x_field = "optimizer_steps"
    if model_id == AAE_MODEL_ID and x_field not in history[0]:
        x_field = "epoch"
    steps = [record[x_field] for record in history]
    tenv = [record["validation_mse"]["tenv"] for record in history]
    best_index = tenv.index(min(tenv))
    figure = Figure(figsize=(8, 4.5), layout="constrained")
    axis = figure.subplots()
    axis.set_title(f"{model_id} training and validation loss")
    axis.set_xlabel("Optimizer updates" if x_field == "optimizer_steps" else "Epoch")
    axis.set_ylabel("Full-transcriptome MSE")
    axis.set_yscale("log")
    axis.xaxis.set_major_locator(MaxNLocator(integer=True))
    axis.grid(alpha=0.25)
    axis.plot(steps, [record["training_mse"] for record in history],
              color="#57606a", marker=".", label="Joint training (batch losses)")
    axis.plot(steps, [record["validation_mse"]["tahoe"] for record in history],
              color="#00748b", marker=".", label="Tahoe validation")
    color = MATRIX_AE_COLOR if model_id == MATRIX_AE_MODEL_ID else "#6336a3"
    axis.plot(steps, tenv, color=color, marker=".", label="tenv validation")
    axis.axvline(steps[best_index], color=color, linestyle=":", alpha=0.6,
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
        command.add_argument("--interface", type=Path,
                             help=f"shared expression artifact (default: <models>/{INTERFACE_ID})")
        command.add_argument("--seal", type=Path, default=DEFAULT_SEAL)
        command.add_argument("--fixture", action="store_true")
        if name == "seal":
            roster = command.add_mutually_exclusive_group()
            roster.add_argument("--expanded-aae", action="store_true",
                                help="include both AAE ablation variants in a new comparison")
            roster.add_argument("--matrix-mlp", action="store_true",
                                help="seal the revised eight-model matrix MLP comparison")
            roster.add_argument("--shared-expression", action="store_true",
                                help="seal the six-entry shared expression comparison")
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
    loss_models = {"aae-loss-figure": AAE_MODEL_ID, "matrix-ae-loss-figure": MATRIX_AE_MODEL_ID}
    for name, model_id in loss_models.items():
        loss = commands.add_parser(name)
        loss.add_argument("--history", type=Path, required=True, help=f"{model_id} epochs.json")
        loss.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "figure":
            plot_comparison(json.loads(args.expected.read_text()), args.output)
            return
        if args.command in loss_models:
            plot_attention_loss(json.loads(args.history.read_text()), args.output, loss_models[args.command])
            return
        if (args.command == "seal" and args.expanded_aae
                and args.seal.resolve() == DEFAULT_SEAL.resolve()):
            raise ValueError("expanded comparison requires a separate --seal path")
        if (args.command == "seal" and args.matrix_mlp
                and args.seal.resolve() in {DEFAULT_SEAL.resolve(), EXPANDED_SEAL.resolve()}):
            raise ValueError("matrix MLP comparison requires a separate --seal path")
        if args.command == "seal" and args.shared_expression:
            if args.seal.resolve() in {DEFAULT_SEAL.resolve(), EXPANDED_SEAL.resolve()}:
                raise ValueError("shared expression comparison requires a separate --seal path")
            if args.validation_dir is None:
                raise ValueError("shared expression comparison requires --validation-dir")
        args.interface = args.interface or args.models / INTERFACE_ID
        if args.command == "seal":
            corpus = Corpus(args.release, args.fixture)
            split = Split(corpus, args.split)
            model_ids = MODEL_IDS
            if args.expanded_aae:
                model_ids = EXPANDED_MODEL_IDS
            elif args.matrix_mlp:
                model_ids = MATRIX_MLP_COMPARISON_IDS
            if args.shared_expression:
                reports = {
                    model: json.loads(validation_path(model, args.validation_dir).read_text())
                    for model in (*SHARED_COMPARISON_IDS, INTERFACE_ID)
                }
                result = selected_shared_models(args.models, args.interface, corpus, split, reports)
            else:
                reports = {model: json.loads(validation_path(model, args.validation_dir).read_text())
                           for model in model_ids}
                result = selected_models(args.models, corpus, split, reports, model_ids)
            publish_json(args.seal, result)
        elif args.command == "score":
            if args.output.exists():
                raise ValueError(f"score report already exists: {args.output}")
            result = compare(args)
            publish_json(args.output, result)
        else:
            expected = json.loads(args.expected.read_text())
            result = compare(args)
            verify_replay(result, expected)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"error: {error}\n")
    print(json_bytes(result).decode(), end="")


if __name__ == "__main__":
    main()
