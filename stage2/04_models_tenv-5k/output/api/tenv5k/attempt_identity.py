"""Non-scientific execution-attempt provenance; never touches model RNG."""
import hashlib
import json
from pathlib import Path


def attempt_binding(registration):
    nonce = registration.get('attempt_nonce')
    if not isinstance(nonce, str) or len(nonce) != 64 or any(c not in '0123456789abcdef' for c in nonce):
        raise ValueError('preregistered attempt nonce missing or invalid')
    encoded = json.dumps(registration, sort_keys=True, separators=(',', ':')).encode()
    return dict(attempt_nonce=nonce, registration_sha256=hashlib.sha256(encoded).hexdigest())


def read_registration(path, *, adapter, scenario, fold, seed, device):
    if path is None:
        return None  # Standalone learner/smoke; matrix acceptance still requires binding.
    registration = json.loads(Path(path).read_text())
    if registration['run'] != dict(adapter=adapter, scenario=scenario, fold=fold, seed=seed) or registration['device'] != str(device):
        raise ValueError('attempt registration run/device mismatch')
    return attempt_binding(registration)


def verify_artifact_binding(root, registration, *, require_manifest=True):
    """Check training-produced config/checkpoint/manifest against the slot attempt."""
    root = Path(root)
    expected = attempt_binding(registration)
    config = json.loads((root/'config.json').read_text())
    if config.get('execution_binding') != expected:
        raise ValueError('training configuration attempt registration mismatch')
    manifest = root/'manifest.json'
    if require_manifest or manifest.exists():
        if json.loads(manifest.read_text()).get('execution_binding') != expected:
            raise ValueError('training manifest attempt registration mismatch')
    import torch
    saved = torch.load(root/'checkpoint.pt', map_location='cpu', weights_only=False)
    if saved.get('execution_binding') != expected or saved['config'] != config:
        raise ValueError('training checkpoint attempt registration mismatch')
    return expected
