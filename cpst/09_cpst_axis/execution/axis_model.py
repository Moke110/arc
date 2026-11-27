"""Pair the CarPS-Transcriptome axis models to the native Step 08 UniMol initialization."""

from hashlib import sha256
import json

import torch

from axis_protocol import ARM
from chem_runs import AlphaAdapter, WIDTHS

REFERENCE_WIDTH = 5000
RESIZED_PARAMETERS = frozenset({
    "encoder_x1.0.weight", "encoder_delta.0.weight", "decoder_x1.4.weight",
    "decoder_x1.4.bias", "decoder_delta.4.weight", "decoder_delta.4.bias",
})


def named_tensor_digest(state, names):
    """Hash names with their corresponding ordered initialized parameter bytes."""
    digest = sha256()
    for name in sorted(names):
        digest.update(name.encode() + b"\0" + state[name].detach().cpu().numpy().tobytes())
    return digest.hexdigest()


def torch_streams(device):
    """Record the full CPU and active accelerator generator state."""
    cuda = torch.cuda.get_rng_state(device) if device.type == "cuda" else None
    return {
        "torch_cpu_rng_sha256": sha256(torch.get_rng_state().numpy().tobytes()).hexdigest(),
        "torch_cuda_rng_sha256": sha256(cuda.numpy().tobytes()).hexdigest() if cuda is not None else None,
    }


class AxisAdapter(AlphaAdapter):
    def __init__(self, source, transform):
        super().__init__(source, transform, chemical_width=WIDTHS[ARM])
        self.epoch_streams = []

    def build(self, seed, device, *, n_genes=REFERENCE_WIDTH):
        """Retain reference generator states and copy all unchanged initialized tensors."""
        if n_genes == REFERENCE_WIDTH:
            return super().build(seed, device, n_genes=n_genes)
        reference = super().build(seed, device, n_genes=REFERENCE_WIDTH)
        cuda = [self.device] if self.device.type == "cuda" else []
        with torch.random.fork_rng(devices=cuda):
            model = super().build(seed, device, n_genes=n_genes)
        before, after = reference.state_dict(), model.state_dict()
        if before.keys() != after.keys() or {
            name for name in before if before[name].shape != after[name].shape
        } != RESIZED_PARAMETERS:
            raise ValueError("prediction-axis width must resize only the four expression projections")
        model.load_state_dict({name: after[name] if name in RESIZED_PARAMETERS else before[name]
                               for name in before})
        return model

    def train_epoch(self, model, inputs, rows, optimizer, rng):
        loss = super().train_epoch(model, inputs, rows, optimizer, rng)
        self.epoch_streams.append({
            **torch_streams(self.device),
            "numpy_rng_sha256": sha256(json.dumps(rng.bit_generator.state, sort_keys=True).encode()).hexdigest(),
        })
        return loss


def paired_initialization(model):
    """Record unchanged initialized tensors, all resized shapes and the paired stream."""
    state = model.state_dict()
    return {
        "invariant_sha256": named_tensor_digest(state, state.keys() - RESIZED_PARAMETERS),
        "axis_projection_shapes": {name: list(state[name].shape) for name in sorted(RESIZED_PARAMETERS)},
        **torch_streams(next(model.parameters()).device),
    }
