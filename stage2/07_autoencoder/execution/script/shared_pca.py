"""Fit one gene-axis PCA basis, pooling centered training channels."""

from pathlib import Path
import sys

import numpy as np
import torch

from benchmark import SEED, SOURCES, batches, publish_artifact, verify_artifact
from shared_expression import dimensions, load_interface


FIT_BATCH = 32
SCORE_BATCH = 16
OVERSAMPLE = 64
POWER_PASSES = 2


class MatrixPCA:
    def __init__(self, center, components):
        self.center = center
        self.components = components

    def encode(self, representation):
        return self.components @ (representation - self.center)

    def decode(self, tokens):
        return self.components.T @ tokens + self.center


@torch.inference_mode()
def fit_shared_pca(corpus, split, output, device, interface_path):
    if Path(output).exists():
        raise ValueError(f"artifact already exists: {output}")
    record, interface = load_interface(interface_path, corpus, split, device)
    rank, channels = dimensions(corpus)
    counts = {source: len(split.rows[source]["train"]) for source in SOURCES}
    genes = corpus.gene_count
    width = min(rank + OVERSAMPLE, genes)
    torch.set_float32_matmul_precision("highest")

    def representations():
        for source in SOURCES:
            for values in batches(corpus, split, source, "train", FIT_BATCH):
                yield interface.generate(torch.from_numpy(values).to(device))

    # The center keeps the mean of each channel of each ordered gene.
    print("Shared PCA: accumulating training representation center", flush=True, file=sys.stderr)
    total = torch.zeros((genes, channels), dtype=torch.float64, device=device)
    for values in representations():
        total += values.sum(dim=0, dtype=torch.float64)
    center = (total / sum(counts.values())).float()

    def centered_channels():
        for values in representations():
            # Each (profile, channel) pair contributes one gene vector.
            pooled = (values - center).transpose(1, 2).reshape(-1, genes)
            yield pooled.half() if device == "cuda" else pooled

    def product(left, right):
        if device == "cuda":
            return torch.mm(left, right, out_dtype=torch.float32)
        return left @ right

    rng = torch.Generator().manual_seed(SEED)
    sketch = torch.randn(genes, width, generator=rng, dtype=torch.float32).to(device)
    basis = torch.linalg.qr(sketch).Q
    for index in range(POWER_PASSES):
        print(f"Shared PCA: covariance pass {index + 1}/{POWER_PASSES}", flush=True, file=sys.stderr)
        projected = torch.zeros_like(basis)
        right = basis.half() if device == "cuda" else basis
        for values in centered_channels():
            scores = product(values, right)
            left = scores.half() if device == "cuda" else scores
            projected.add_(product(values.T, left))
        basis = torch.linalg.qr(projected).Q

    print("Shared PCA: Rayleigh–Ritz projection", flush=True, file=sys.stderr)
    compressed = basis.new_zeros((width, width))
    right = basis.half() if device == "cuda" else basis
    for values in centered_channels():
        scores = product(values, right)
        compressed.addmm_(scores.T, scores)
    _, directions = torch.linalg.eigh((compressed + compressed.T) * 0.5)
    components = (basis @ directions.flip((1,))[:, :rank]).T.contiguous()
    if not (torch.isfinite(center).all() and torch.isfinite(components).all()):
        raise ValueError("nonfinite fitted matrix PCA")

    metadata = {
        "fit_role": "train",
        "fit_profiles": counts,
        "interface_manifest_sha256": record["artifact_manifest_sha256"],
        "latent_shape": [rank, channels],
        "trainable_parameters": 0,
        "fitted_parameters": center.numel() + components.numel(),
        "batch_size": FIT_BATCH,
        "score_batch_size": SCORE_BATCH,
        "solver": "streamed randomized covariance subspace iteration with Rayleigh–Ritz",
        "power_passes": POWER_PASSES,
        "sketch_width": width,
        "seed": SEED,
        "center_reduction": "float64 accumulation, float32 center",
        "covariance_products": "float16 inputs with float32 output" if device == "cuda" else "float32",
        "compressed_covariance": "float32",
        "fit_device": device,
        "torch_version": str(torch.__version__),
    }
    with publish_artifact(output, corpus, split, "PCA_1024", metadata) as directory:
        np.save(directory / "center.npy", center.cpu().numpy(), allow_pickle=False)
        np.save(directory / "components.npy", components.cpu().numpy(), allow_pickle=False)
    return verify_artifact(output, corpus, split)


def shared_pca_predictor(path, record, corpus, split, device, interface_path):
    interface_record, interface = load_interface(interface_path, corpus, split, device)
    metadata = record["metadata"]
    rank, channels = dimensions(corpus)
    if metadata["interface_manifest_sha256"] != interface_record["artifact_manifest_sha256"]:
        raise ValueError("matrix PCA is bound to a different expression interface")
    if (record["model_id"] != "PCA_1024" or metadata["latent_shape"] != [rank, channels]
            or set(record["files"]) != {"center.npy", "components.npy"}
            or metadata["score_batch_size"] != SCORE_BATCH):
        raise ValueError("invalid matrix PCA artifact")

    def load(name, shape):
        array = np.load(Path(path) / name, allow_pickle=False)
        if array.shape != shape or array.dtype != np.float32 or not np.isfinite(array).all():
            raise ValueError(f"invalid matrix PCA payload: {name}")
        return torch.from_numpy(array).to(device)

    projection = MatrixPCA(load("center.npy", (corpus.gene_count, channels)),
                           load("components.npy", (rank, corpus.gene_count)))
    torch.set_float32_matmul_precision("high" if device == "cuda" else "highest")

    @torch.inference_mode()
    def predict(values):
        representation = interface.generate(torch.from_numpy(values).to(device))
        return interface.restore(projection.decode(projection.encode(representation))).cpu().numpy()

    return predict, SCORE_BATCH
