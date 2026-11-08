"""Resolve one tenv initialization for both zero-shot Tahoe and Transfer."""
from bisect import bisect_left
import hashlib
import json
import math
from pathlib import Path

FORMAT = 'step05-tahoe-bridge-v1'
LIMITATION = ('Tahoe dose embeddings use log-dose interpolation with endpoint clamping. '
              'The 10 µM tenv dose receives the learned 5 µM Tahoe row.')
CHECKPOINT_FIELDS = ('identity', 'recipe_sha256', 'best_step', 'best_validation_loss', 'model_state_sha256')
DOSE_ROWS = 'native.cov_encoder.dose_encoder.weight'
DOSE_LOOKUP = 'native.cov_encoder.dose_indices'
CHEMICAL_LOOKUP = 'native.cov_encoder.drug_encoder.0.weight'


def file_digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def request_digest(request):
    return hashlib.sha256(json.dumps(request, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def state_digest(weights):
    """Match the digest of the native adapter's ordered state dictionary."""
    digest = hashlib.sha256()
    for name, value in weights.items():
        digest.update(name.encode())
        digest.update(memoryview(value.detach().cpu().contiguous().numpy()).cast('B'))
    return digest.hexdigest()


def dose_layout(keys):
    """Reconstruct the native positive dose vocabulary and its zero-dose control index."""
    doses = sorted({dose for _, dose in keys})
    if not doses or any(not math.isfinite(dose) or dose <= 0 for dose in doses):
        raise ValueError('positive finite dose vocabulary required')
    indices = {dose: index for index, dose in enumerate(doses, 1)}
    return doses, [indices[dose] for _, dose in keys] + [0]


def dose_mapping(source, recipient):
    records = []
    for dose in recipient:
        position = bisect_left(source, dose)
        left, right = (position, position) if position < len(source) and source[position] == dose else (
            max(position - 1, 0), min(position, len(source) - 1))
        fraction = (math.log(dose / source[left]) / math.log(source[right] / source[left])
                    if left != right else 0.)
        records.append({'recipient_um': dose, 'source_lower_um': source[left], 'source_upper_um': source[right],
                        'source_upper_weight': fraction})
    return records


def mapped_weights(model, source):
    import torch

    checkpoint_path = Path(source['path'])
    if file_digest(checkpoint_path) != source['checkpoint_sha256']:
        raise ValueError('Tahoe EMA checkpoint checksum mismatch')
    checkpoint = torch.load(checkpoint_path, map_location='cpu', weights_only=True, mmap=True)
    if (checkpoint['format'] != 'step04-ema-prediction-v1'
            or any(checkpoint[key] != source[key] for key in CHECKPOINT_FIELDS)):
        raise ValueError('Tahoe EMA selection metadata mismatch')
    weights = checkpoint['model_state']
    if state_digest(weights) != source['model_state_sha256']:
        raise ValueError('Tahoe EMA state checksum mismatch')

    native = model.state_dict()
    target_doses, target_indices = dose_layout(model.keys)
    source_doses = source['doses_um']
    source_rows = weights[DOSE_ROWS]
    if (not source_doses or source_doses != sorted(set(source_doses))
            or any(not math.isfinite(dose) or dose <= 0 for dose in source_doses)
            or len(source_rows) != len(source_doses) + 1
            or weights[DOSE_LOOKUP].tolist() != source['dose_indices']):
        raise ValueError('Tahoe dose rows differ from the source case vocabulary')
    if (native[DOSE_LOOKUP].tolist() != target_indices or len(native[DOSE_ROWS]) != len(target_doses) + 1
            or len(native[CHEMICAL_LOOKUP]) != len(target_indices)
            or weights.keys() != native.keys()):
        raise ValueError('recipient lookup differs from the native tenv vocabulary')

    rows = [source_rows[0]]
    by_dose = dict(zip(source_doses, source_rows[1:]))
    mapping = dose_mapping(source_doses, target_doses)
    for record in mapping:
        lower, upper, fraction = (record[key] for key in
                                  ('source_lower_um', 'source_upper_um', 'source_upper_weight'))
        rows.append(torch.lerp(by_dose[lower], by_dose[upper], fraction) if fraction else by_dose[lower])
    weights[DOSE_ROWS] = torch.stack(rows)
    weights[DOSE_LOOKUP] = native[DOSE_LOOKUP]
    weights[CHEMICAL_LOOKUP] = native[CHEMICAL_LOOKUP]
    if any(weights[key].shape != value.shape or weights[key].dtype != value.dtype
           for key, value in native.items()):
        raise ValueError('Tahoe state is incompatible with the recipient model')
    return weights, mapping


def bridge_record(model, source, recipient_request):
    weights, mapping = mapped_weights(model, source)
    return weights, {
        'schema': FORMAT, 'source': source,
        'recipient_identity': recipient_request['identity'],
        'recipient_request_sha256': request_digest(recipient_request),
        'dose_mapping': mapping, 'dose_limitation': LIMITATION,
        'recipient_lookup_sha256': state_digest({key: weights[key] for key in (CHEMICAL_LOOKUP, DOSE_LOOKUP)}),
        'mapped_dose_sha256': state_digest({DOSE_ROWS: weights[DOSE_ROWS]}),
        'mapped_state_sha256': state_digest(weights),
    }


def seal_initialization(model, directory, source, recipient_request):
    """Seal the mapping without duplicating the large learned source checkpoint."""
    _, record = bridge_record(model, source, recipient_request)
    path = directory / 'bridge.json'
    with path.open('x') as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return file_digest(path)


def load_initialization(model, directory, expected_sha256, recipient_request):
    """Load and verify the pinned initialization before prediction or optimizer creation."""
    path = directory / 'bridge.json'
    if file_digest(path) != expected_sha256:
        raise ValueError('Tahoe bridge record checksum mismatch')
    record = json.loads(path.read_text())
    if (record['schema'] != FORMAT or record['recipient_identity'] != recipient_request['identity']
            or record['recipient_request_sha256'] != request_digest(recipient_request)):
        raise ValueError('Tahoe bridge recipient request mismatch')
    weights, expected = bridge_record(model, record['source'], recipient_request)
    if record != expected:
        raise ValueError('Tahoe bridge mapped state differs from the sealed initialization')
    model.load_state_dict(weights, strict=True)
    if state_digest(model.state_dict()) != record['mapped_state_sha256']:
        raise ValueError('loaded Tahoe bridge state checksum mismatch')
    return record
