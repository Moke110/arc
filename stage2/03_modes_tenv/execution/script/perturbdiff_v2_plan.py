"""Read-only PerturbDiff v2 budget preflight; this CLI never migrates or dispatches.

Scheduler sampled training is a lower bound while an attempt is live. Such plans
are explicitly provisional and cannot authorize dispatch. An absent state is zero
only when no run was ever dispatched, which must be verified at the stop boundary.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path


def plan(queue_path, states_path):
    queue_bytes = queue_path.read_bytes()
    queue = json.loads(queue_bytes)
    jobs = queue['jobs']
    selected = [job for job in jobs if job['model'] == 'perturbdiff']
    expected = {(fold, seed) for fold in range(3) for seed in (17, 42, 2026)}
    if len(selected) != 9 or {(job['fold'], job['seed']) for job in selected} != expected:
        raise ValueError('exactly nine unique PerturbDiff fold/seed slots required')
    migrations = []
    for job in selected:
        key = f"perturbdiff_fold{job['fold']}_seed{job['seed']}"
        state_path = states_path / (key + '.json')
        state_bytes = state_path.read_bytes() if state_path.exists() else None
        state = json.loads(state_bytes) if state_bytes else {}
        if state_bytes and 'training_seconds' not in state:
            raise ValueError('historical training debit missing: ' + key)
        consumed = state.get('training_seconds', 0)
        if not isinstance(consumed, (float, int)) or not math.isfinite(consumed) or consumed < 0:
            raise ValueError('invalid historical training debit: ' + key)
        migrations.append({'run': key, 'legacy_identity': state.get('identity'),
            'legacy_status': state.get('status', 'never_dispatched_unverified'),
            'state_sha256': hashlib.sha256(state_bytes).hexdigest() if state_bytes else None,
            'historical_training_seconds_lower_bound': consumed,
            'remaining_training_seconds': max(0, 14400 - consumed),
            'initialization': 'random_init', 'legacy_reuse_as_v2': False,
            'boundary_reconciliation_required': True})
    return {'schema': 'perturbdiff-v2-preflight-v1', 'dispatch_authorized': False,
        'queue_sha256': hashlib.sha256(queue_bytes).hexdigest(),
        'protocol': {'label': 'adapted native validation stopping; pending frozen identity',
            'validation_interval': 1000, 'checkpoint_interval': 1000,
            'validation_groups': 'all', 'validation_reduction': 'group-count-weighted',
            'early_stop_min_steps': 20000, 'early_stop_patience': 10,
            'early_stop_min_delta': 0.01, 'patience_begins_at_min_steps': True,
            'earliest_early_stop_step': 30000, 'resume_supported': False,
            'cumulative_training_ceiling_seconds': 14400},
        'migrations': migrations,
        'untouched_jobs': [job for job in jobs if job['model'] != 'perturbdiff'],
        'limitations': ['Debits are sampled lower bounds until stopped boundary reconciliation.',
            'No native checkpoint or prediction hash inventory is performed by this preflight.',
            'No train completion is promised within remaining original budget.',
            'No queue, adapter, identity, checkpoint or scheduler state is modified.']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--queue', required=True, type=Path)
    parser.add_argument('--states', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = plan(args.queue, args.states)
    # Exclusive creation avoids accidentally overwriting prior planning evidence.
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write('\n')
