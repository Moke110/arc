#!/usr/bin/env python3
"""Exercise the positive-dose initialization seam with tiny native models."""
import argparse
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'script'))

from dose_reset import load_dose_reset
from artifacts import sha256
from formal import use_runtime
from tahoe_bridge import (DOSE_ROWS, dose_layout, load_initialization, seal_initialization,
                          state_digest)


def replay(runtime, root):
    use_runtime(runtime)
    from native_training import FORMAT, JointPerturbDiff, np, seed_randomness, torch

    def model(doses):
        cases = [{'drug_key': 'drug', 'dose_um': str(dose)} for dose in doses]
        inputs = SimpleNamespace(n_genes=3, cases={'train': cases},
                                 morgan=np.eye(1, 2048, dtype=np.float32), drug_index={'drug': 0})
        config = {'chemical_dimension': 16, 'learning_rate': 0.0002, 'weight_decay': 0,
                  'expression_scale': 1, 'batch_groups': 2}
        return JointPerturbDiff(inputs, config, root)

    torch.set_num_threads(2)
    seed_randomness(23)
    source = model((0.05, 0.5, 5))
    with torch.no_grad():
        source.get_parameter(DOSE_ROWS).fill_(2)
    weights = source.state_dict()
    selected = {'identity': {'dataset': 'tahoe', 'seed': 23}, 'recipe_sha256': 'fixture',
                'best_step': 1, 'best_validation_loss': 0.5, 'model_state_sha256': state_digest(weights)}
    checkpoint = root / 'source.pt'
    torch.save({**selected, 'format': FORMAT, 'model_state': weights}, checkpoint)
    doses, indices = dose_layout(source.keys)
    binding = {**selected, 'path': str(checkpoint), 'checkpoint_sha256': sha256(checkpoint),
               'doses_um': doses, 'dose_indices': indices}

    request = {'identity': {'dataset': 'tenv', 'fold': 0}}
    seed_randomness(17)
    original = model((0.2, 1, 10))
    bridge_sha = seal_initialization(original, root, binding, request)
    load_initialization(original, root, bridge_sha, request)

    seed_randomness(17)
    reset = model((0.2, 1, 10))
    native_rows = reset.state_dict()[DOSE_ROWS][1:].clone()
    rng = torch.get_rng_state()
    bridge, evidence = load_dose_reset(reset, root, bridge_sha, request)
    assert torch.equal(rng, torch.get_rng_state())
    assert evidence['reset_doses_um'] == [0.2, 1.0, 10.0]
    assert torch.equal(reset.state_dict()[DOSE_ROWS][1:], native_rows)
    assert torch.equal(reset.state_dict()[DOSE_ROWS][0], original.state_dict()[DOSE_ROWS][0])
    assert all(torch.equal(value, original.state_dict()[key])
               for key, value in reset.state_dict().items() if key != DOSE_ROWS)
    assert bridge['mapped_state_sha256'] == state_digest(original.state_dict())
    assert reset.get_parameter(DOSE_ROWS).requires_grad
    print('PASS native dose reset preserves the bridge outside the three positive rows and preserves RNG')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        replay(args.runtime.resolve(), Path(temp))
