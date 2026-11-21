"""Initialize a tenv recipient from the frozen Tahoe source model."""
import numpy as np
import torch

from tsg_round1.model import transform_doses

BRIDGE_TOLERANCE = 1e-6


def initialize_from_tahoe(source, recipient, source_transform, target_transform, doses, *, held_dose):
    """Copy the winner and preserve its dose embedding in tenv coordinates."""
    doses = tuple(sorted(set(float(dose) for dose in doses)))
    if held_dose not in doses:
        raise ValueError('dose handoff must cover the held dose')
    source_doses = transform_doses(doses, source_transform)
    target_doses = transform_doses(doses, target_transform)
    weight_scale = target_transform['std'] / source_transform['std']
    bias_shift = (target_transform['mean'] - source_transform['mean']) / source_transform['std']
    recipient.load_state_dict(source.state_dict())
    source_layer = source.dose_embedding[0]
    target_layer = recipient.dose_embedding[0]
    with torch.no_grad():
        weight = source_layer.weight.detach().to(target_layer.weight.device, dtype=torch.float64)
        bias = source_layer.bias.detach().to(target_layer.bias.device, dtype=torch.float64)
        target_layer.weight.copy_(weight * weight_scale)
        target_layer.bias.copy_(bias + weight[:, 0] * bias_shift)
        source_embedding = source.dose_embedding(
            torch.as_tensor(source_doses, device=source_layer.weight.device).reshape(-1, 1),
        ).cpu()
        recipient_embedding = recipient.dose_embedding(
            torch.as_tensor(target_doses, device=target_layer.weight.device).reshape(-1, 1),
        ).cpu()
    error = (source_embedding - recipient_embedding).abs()
    allowance = BRIDGE_TOLERANCE * (1 + source_embedding.abs())
    worst_ratio = (error / allowance).amax(dim=1).numpy()
    difference = error.amax(dim=1).numpy()
    if not np.isfinite(worst_ratio).all() or (worst_ratio > 1).any():
        raise ValueError('pre-update dose embedding differs after affine conversion')
    return dict(
        held_dose_um=held_dose, absolute_tolerance=BRIDGE_TOLERANCE,
        relative_tolerance=BRIDGE_TOLERANCE, weight_scale=weight_scale, bias_shift=bias_shift,
        source_transform=source_transform, target_transform=target_transform,
        checks=[dict(dose_um=dose, max_absolute_error=float(delta),
                     max_fraction_of_tolerance=float(fraction))
                for dose, delta, fraction in zip(doses, difference, worst_ratio, strict=True)],
    )


def verify_bridge(record, source_transform, target_transform, doses, *, held_dose):
    """Check a recorded conversion against the frozen dose coordinates and coverage."""
    doses = sorted(set(float(dose) for dose in doses))
    if (record['source_transform'] != source_transform
            or record['target_transform'] != target_transform
            or record['held_dose_um'] != held_dose
            or held_dose not in doses
            or [row['dose_um'] for row in record['checks']] != doses
            or record['absolute_tolerance'] != BRIDGE_TOLERANCE
            or record['relative_tolerance'] != BRIDGE_TOLERANCE
            or record['weight_scale'] != target_transform['std'] / source_transform['std']
            or record['bias_shift'] != (
                target_transform['mean'] - source_transform['mean']
            ) / source_transform['std']):
        raise ValueError('recorded dose handoff coordinates or coverage differ')
    errors = [row['max_absolute_error'] for row in record['checks']]
    fractions = [row['max_fraction_of_tolerance'] for row in record['checks']]
    if (not np.isfinite(errors + fractions).all() or min(errors) < 0
            or min(fractions) < 0 or max(fractions) > 1):
        raise ValueError('recorded dose handoff numerical evidence differs')
