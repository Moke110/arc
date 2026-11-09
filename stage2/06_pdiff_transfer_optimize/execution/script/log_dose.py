"""Continuous log-dose covariates and the loss-selected Tahoe EMA bridge."""
import math

import torch
from torch import nn

from native_training import FORMAT
from tahoe_bridge import (CHECKPOINT_FIELDS, CHEMICAL_LOOKUP, DOSE_LOOKUP, dose_layout,
                          request_digest, state_digest)
from formal_matrix import check_binding
from training_inputs import json_hash

REFERENCE_DOSE_UM = 0.05


class LogDoseEncoder(nn.Module):
    def __init__(self, doses, embedding_dim):
        super().__init__()
        # The vocabulary is rebuilt from the local cases on fit and replay.
        # Keep it out of the tensors averaged by Lightning EMA.
        self.features = (0., *(math.log10(dose / REFERENCE_DOSE_UM) for dose in doses))
        self.projection = nn.Linear(1, embedding_dim)
        self.control = nn.Parameter(torch.empty(embedding_dim))
        nn.init.normal_(self.control)

    def forward(self, indices):
        features = self.projection.weight.new_tensor(self.features)[indices]
        positive = self.projection(features.unsqueeze(-1))
        return torch.where(indices.unsqueeze(-1) == 0, self.control, positive)


def install_encoder(model):
    """Replace the native categorical encoder before device and optimizer setup."""
    covariates = model.native.cov_encoder
    doses, indices = dose_layout(model.keys)
    if covariates.dose_indices.tolist() != indices:
        raise ValueError('log encoder dose vocabulary differs from native case order')
    # Retain the native training RNG stream outside this initialization.
    with torch.random.fork_rng(devices=[]):
        covariates.dose_encoder = LogDoseEncoder(doses, covariates.dose_encoder.embedding_dim)
    return model


def load_initialization(model, source, recipient_request):
    """Load the selected log Tahoe state and reconstruct the recipient lookups."""
    metadata = source['metadata']
    identity = metadata['identity']
    recipient_identity = recipient_request['identity']
    if ((identity['dataset'], identity['fold'], identity['arm'], identity['regime'])
            != ('tahoe', recipient_identity['fold'], 'log', 'cold-cm-drug')
            or recipient_identity['dataset'] != 'tenv'
            or identity.get('gene_sha256') != recipient_identity.get('gene_sha256')):
        raise ValueError('log Tahoe and tenv recipient identities differ')
    binding = check_binding(source['checkpoint'])
    if binding['sha256'] != metadata['checkpoint_sha256']:
        raise ValueError('log Tahoe EMA binding differs from its selection')
    checkpoint = torch.load(binding['path'], map_location='cpu', weights_only=True, mmap=True)
    if (checkpoint['format'] != FORMAT
            or any(checkpoint[key] != metadata[key] for key in CHECKPOINT_FIELDS)
            or json_hash(checkpoint['recipe']) != metadata['recipe_sha256']):
        raise ValueError('log Tahoe EMA selection metadata mismatch')
    weights = checkpoint['model_state']
    if state_digest(weights) != metadata['model_state_sha256']:
        raise ValueError('log Tahoe EMA model state checksum mismatch')

    recipient_doses, recipient_indices = dose_layout(model.keys)
    native = model.state_dict()
    lookups = {key: native[key] for key in (CHEMICAL_LOOKUP, DOSE_LOOKUP)}
    if (weights[DOSE_LOOKUP].tolist() != source['dose_indices']
            or len(weights[CHEMICAL_LOOKUP]) != len(source['dose_indices'])):
        raise ValueError('log Tahoe EMA lookups differ from the frozen source vocabulary')
    encoder = model.native.cov_encoder.dose_encoder
    if (not isinstance(encoder, LogDoseEncoder) or weights.keys() != native.keys()
            or native[DOSE_LOOKUP].tolist() != recipient_indices
            or len(native[CHEMICAL_LOOKUP]) != len(recipient_indices)):
        raise ValueError('recipient log encoder or lookups differ from its native case vocabulary')
    weights.update(lookups)
    if any(weights[key].shape != value.shape or weights[key].dtype != value.dtype
           for key, value in native.items()):
        raise ValueError('log Tahoe EMA is incompatible with the recipient model')
    model.load_state_dict(weights, strict=True)
    expected = state_digest(weights)
    if state_digest(model.state_dict()) != expected:
        raise ValueError('loaded log Tahoe state checksum mismatch')
    dose_state = {key: value for key, value in weights.items()
                  if key.startswith('native.cov_encoder.dose_encoder.')}
    return {
        'schema': 'step06-log-tahoe-bridge-v1', 'source_checkpoint': binding, 'source_selection': metadata,
        'recipient_identity': recipient_identity, 'recipient_request_sha256': request_digest(recipient_request),
        'reference_dose_um': REFERENCE_DOSE_UM,
        'recipient_doses_um': recipient_doses, 'recipient_log_features': list(encoder.features[1:]),
        'recipient_lookup_sha256': state_digest(lookups),
        'transferred_dose_sha256': state_digest(dose_state), 'mapped_state_sha256': expected,
    }
