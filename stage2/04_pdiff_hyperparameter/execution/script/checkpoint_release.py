"""Verify the selected checkpoint release independently of retrieval scores."""
from pathlib import Path
import sys
REPO = next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())
sys.path.insert(0, str(REPO / 'stage2/03_modes_tenv/execution/script'))
from artifacts import SEEDS, sha256
FOLDS = range(3)
SLOTS = {(fold, seed) for fold in FOLDS for seed in SEEDS}
CHECKPOINT_FIELDS = ('identity', 'recipe_sha256', 'checkpoint_sha256', 'best_step',
                     'best_validation_loss', 'model_state_sha256')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_checkpoint_index(entries, folds, manifests):
    by_slot = {(row['identity']['fold'], row['identity']['seed']): row for row in entries}
    require(by_slot.keys() == SLOTS and len(entries) == len(SLOTS),
            'tenv release needs nine unique checkpoint identities')
    for (fold, seed), row in by_slot.items():
        manifest = manifests[fold, seed]
        fields = (*CHECKPOINT_FIELDS, 'run_sha256', 'replay_sha256')
        require(all(row[key] == manifest[key] for key in fields)
                and row['directory'] == f'fold-{fold}/seed-{seed}'
                and row['request_sha256'] == folds[fold]['request_sha256']
                and row['policy_sha256'] == folds[fold]['policy_sha256']
                and row['raw_recipe_sha256'] == folds[fold]['recipes'][str(seed)]['sha256'],
                f'tenv checkpoint index mismatch for fold {fold}, seed {seed}')
    return [by_slot[slot] for slot in sorted(by_slot)]


def verify_local_release(root, formal_root, index_path, checkpoints):
    require(sha256(root / 'checkpoints.json') == sha256(index_path),
            'local tenv checkpoint index differs from the recorded index')
    evidence = {}
    for row in checkpoints:
        slot = root / row['directory']
        for filename, field in [('run.json', 'run_sha256'), ('replay.json', 'replay_sha256'),
                                ('request.json', 'request_sha256'), ('recipe.json', 'raw_recipe_sha256')]:
            require(sha256(slot / filename) == row[field], f'local {slot / filename} hash mismatch')
        fold, seed = row['identity']['fold'], row['identity']['seed']
        fitted = formal_root / 'slots' / f'fold-{fold}-seed-{seed}' / 'fit/best.pt'
        checkpoint = slot / 'best.pt'
        require(checkpoint.is_file() and not checkpoint.is_symlink()
                and checkpoint.samefile(fitted) and checkpoint.stat().st_size > 0,
                f'{checkpoint} must hard-link to the verified formal fit')
        require(sha256(checkpoint) == row['checkpoint_sha256'], f'local {checkpoint} hash mismatch')
        evidence[row['directory']] = {
            'bytes': checkpoint.stat().st_size, 'same_file_as_formal_fit': True,
            'metadata_hashes_verified': True, 'checkpoint_hash_verified': True,
        }
    return evidence
