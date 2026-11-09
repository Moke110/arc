"""Deterministic symmetric MLP reconstruction baselines."""

from benchmark import SEED


MODEL_IDS = ("MLP_32", "MLP_128", "MLP_512")
HIDDEN = 512
RECIPE = {
    "seed": SEED,
    "epochs": 20,
    "train_batch_size": 512,
    "validation_batch_size": 64,
    "max_seconds": 1800,
    "optimizer": "AdamW",
    "learning_rate": 0.001,
    "weight_decay": 0.01,
    "precision": "float32",
    "train_order": "PCG64(seed, epoch) shuffle of the joint train membership",
}


def architecture(gene_count, model_id):
    if model_id not in MODEL_IDS:
        raise ValueError(f"unsupported MLP model: {model_id}")
    latent = int(model_id.split("_")[1])
    return {
        "widths": [gene_count, HIDDEN, latent, HIDDEN, gene_count],
        "activations": ["GELU", "linear", "GELU", "linear"],
    }


def build_mlp(widths):
    from torch import nn

    inputs, hidden, latent, decoder, outputs = widths
    return nn.Sequential(
        nn.Linear(inputs, hidden), nn.GELU(),
        nn.Linear(hidden, latent),
        nn.Linear(latent, decoder), nn.GELU(),
        nn.Linear(decoder, outputs),
    )


def fit_mlp(corpus, split, output, model_id, batch_size, device):
    from neural import fit_neural

    design = architecture(corpus.gene_count, model_id)
    recipe = RECIPE.copy()
    if corpus.fixture:
        recipe.update(epochs=3, train_batch_size=2, validation_batch_size=2)
        if batch_size is not None:
            recipe.update(train_batch_size=batch_size, validation_batch_size=batch_size)
    elif batch_size is not None and batch_size != recipe["train_batch_size"]:
        raise ValueError(f"MLP fit requires the frozen batch size {recipe['train_batch_size']}")
    if recipe["train_batch_size"] < 1:
        raise ValueError("batch size must be positive")
    return fit_neural(
        corpus, split, output, model_id, lambda: build_mlp(design["widths"]), recipe, device,
        {"architecture": design, "latent_tokens": None, "latent_scalars": design["widths"][2]},
    )


def mlp_predictor(path, record, device):
    from neural import load_neural

    design = architecture(record["gene_count"], record["model_id"])
    if (record["metadata"]["architecture"] != design
            or record["metadata"]["latent_scalars"] != design["widths"][2]
            or record["metadata"]["latent_tokens"] is not None
            or set(record["files"]) != {"weights.pt", "epochs.json"}):
        raise ValueError("invalid MLP artifact")
    return load_neural(path, build_mlp(design["widths"]), device)
