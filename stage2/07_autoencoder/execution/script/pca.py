"""Joint training PCA via a streamed randomized covariance eigensolver."""

from pathlib import Path
import tempfile

import numpy as np

from benchmark import SEED, SOURCES, batches, fit_mean, publish_artifact, verify_artifact


RANKS = (32, 128, 512)
MODEL_IDS = tuple(f"PCA_{rank}" for rank in RANKS)
OVERSAMPLE = 64
POWER_PASSES = 8


def fit_pca(corpus, split, output, batch_size, device):
    import torch

    output = Path(output)
    if output.exists():
        raise ValueError(f"artifact already exists: {output}")
    fit_profiles = {source: len(split.rows[source]["train"]) for source in SOURCES}
    width = min(max(RANKS) + OVERSAMPLE, corpus.gene_count, sum(fit_profiles.values()) - 1)
    if width < max(RANKS):
        raise ValueError("PCA requires at least 512 genes and 513 training profiles")
    if device == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA is unavailable for PCA fit")
    torch.set_float32_matmul_precision("highest")
    mean = fit_mean(corpus, split, batch_size)
    center = torch.from_numpy(mean).to(device)
    generator = torch.Generator().manual_seed(SEED)
    sketch = torch.randn(corpus.gene_count, width, dtype=torch.float32, generator=generator).to(device)
    basis = torch.linalg.qr(sketch).Q

    def centered_batches():
        for source in SOURCES:
            for values in batches(corpus, split, source, "train", batch_size):
                yield torch.from_numpy(values).to(device) - center

    # Apply XᵀX to the current basis without storing X or its full covariance.
    for _ in range(POWER_PASSES):
        projected = torch.zeros_like(basis)
        for values in centered_batches():
            projected.addmm_(values.T, values @ basis)
        basis = torch.linalg.qr(projected).Q

    compressed = basis.new_zeros((width, width))
    for values in centered_batches():
        scores = values @ basis
        compressed.addmm_(scores.T, scores)
    # Diagonalize the covariance within the subspace and keep descending directions.
    _, directions = torch.linalg.eigh((compressed + compressed.T) * 0.5)
    components = (basis @ directions.flip(dims=(1,))[:, :max(RANKS)]).T.contiguous().cpu().numpy()
    if not np.isfinite(mean).all() or not np.isfinite(components).all():
        raise ValueError("nonfinite fitted PCA")

    settings = {
        "fit_role": "train",
        "fit_profiles": fit_profiles,
        "batch_size": batch_size,
        "solver": "streamed randomized covariance subspace iteration with Rayleigh–Ritz",
        "seed": SEED,
        "oversample": OVERSAMPLE,
        "sketch_width": width,
        "power_passes": POWER_PASSES,
        "center_reduction": "float64 accumulation, float32 inference",
        "solver_dtype": "float32",
        "matmul_precision": "highest",
        "fit_device": device,
        "fit_device_name": torch.cuda.get_device_name() if device == "cuda" else "cpu",
        "torch_version": torch.__version__,
        "cuda_version": torch.version.cuda if device == "cuda" else None,
        "latent_tokens": None,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent, prefix=".pca-") as staging:
        staging = Path(staging)
        records = {}
        for rank, model_id in zip(RANKS, MODEL_IDS):
            artifact = staging / model_id
            metadata = {**settings, "latent_scalars": rank}
            with publish_artifact(artifact, corpus, split, model_id, metadata) as directory:
                np.save(directory / "mean.npy", mean, allow_pickle=False)
                np.save(directory / "components.npy", components[:rank], allow_pickle=False)
            records[model_id] = verify_artifact(artifact, corpus, split)
        staging.rename(output)
    return records


def pca_predictor(artifact, record):
    rank = RANKS[MODEL_IDS.index(record["model_id"])]
    if (set(record["files"]) != {"mean.npy", "components.npy"}
            or record["metadata"]["latent_scalars"] != rank
            or record["metadata"]["latent_tokens"] is not None):
        raise ValueError("invalid PCA artifact")
    center = np.load(artifact / "mean.npy")
    components = np.load(artifact / "components.npy")
    for values, shape in ((center, (record["gene_count"],)),
                          (components, (rank, record["gene_count"]))):
        if values.dtype != np.float32 or values.shape != shape or not np.isfinite(values).all():
            raise ValueError("invalid PCA payload")

    def predict(values):
        scores = (values - center) @ components.T
        return scores @ components + center

    return predict
