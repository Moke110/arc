"""Archive exploratory evidence and transactionally install a replayed Round 1.

This command never removes training products. cleanup-proposal.json lists old
fit directories for a separate, evidence-backed decision after publication.
"""
import argparse
import fcntl
import json
from pathlib import Path
import shutil
import uuid

from rerun_publication import ENTRIES, replay, sha

REPO = next(path for path in Path(__file__).resolve().parents if (path / '.git').exists())
STEP = REPO / 'stage2/09_xpert_optimize'


def atomic_bytes(path, payload):
    temporary = path.with_name(f'.{path.name}-{uuid.uuid4().hex}.tmp')
    try:
        temporary.write_bytes(payload)
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()


def relocated(value, old_prefix, new_prefix):
    if isinstance(value, str):
        return new_prefix + value[len(old_prefix):] if value.startswith(old_prefix + '/') else value
    if isinstance(value, list):
        return [relocated(item, old_prefix, new_prefix) for item in value]
    if isinstance(value, dict):
        return {key: relocated(item, old_prefix, new_prefix) for key, item in value.items()}
    return value


def tree_hashes(root):
    records = {}
    for path in sorted(root.rglob('*')):
        if path.is_symlink():
            raise ValueError(f'historical evidence cannot contain symlinks: {path}')
        if path.is_file():
            records[str(path.relative_to(root))] = sha(path)
    return records


def finalize(repo, staging, destination, archive, inventory, retrieval_manifest, round1_record):
    """Install only after replay; roll back both renames and metadata on error.

    Historical bindings are relocated without changing their payload hashes.
    The old inventory, retrieval manifest and committed provenance are preserved
    byte-for-byte before active binding paths are rewritten.
    """
    repo = Path(repo).absolute()
    paths = [Path(path).absolute() for path in
             (staging, destination, archive, inventory, retrieval_manifest, round1_record)]
    staging, destination, archive, inventory, retrieval_manifest, round1_record = paths
    for path in paths:
        if path.resolve() != path or not path.is_relative_to(repo):
            raise ValueError(f'expected real repository paths: {path}')
    if (staging == destination or archive.is_relative_to(destination)
            or archive.is_relative_to(staging) or staging.is_relative_to(destination)
            or destination.is_relative_to(staging)):
        raise ValueError('staging, archive and destination must be disjoint')
    if archive.exists():
        raise FileExistsError(archive)
    if {path.name for path in staging.iterdir()} != ENTRIES:
        raise ValueError('staging publication top-level entries differ')
    if not (destination / 'labels/test/genes.csv').is_file():
        raise ValueError('historical significance gene family is missing')
    replay(staging)
    checkpoints = tree_hashes(staging / 'checkpoints')
    historical = tree_hashes(destination)
    originals = {inventory: inventory.read_bytes(), retrieval_manifest: retrieval_manifest.read_bytes()}
    old_prefix = str(destination.relative_to(repo))
    new_prefix = str((archive / 'round1').relative_to(repo))
    replacements = {path: json_bytes(relocated(json.loads(payload), old_prefix, new_prefix))
                    for path, payload in originals.items()}
    archive.parent.mkdir(parents=True, exist_ok=True)
    if staging.stat().st_dev != destination.stat().st_dev or archive.parent.stat().st_dev != destination.stat().st_dev:
        raise ValueError('publication and archive require same-filesystem atomic renames')
    with (archive.parent / '.round1-finalize.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        archive.mkdir()
        (archive / 'original').mkdir()
        (archive / 'original/inventory.json').write_bytes(originals[inventory])
        (archive / 'original/retrieval-manifest.json').write_bytes(originals[retrieval_manifest])
        shutil.copytree(round1_record, archive / 'round1_record')
        evidence = {'schema': 'step09-round1-finalization-v1', 'status': 'prepared',
                    'historical_publication_sha256': historical,
                    'selected_checkpoints_sha256': checkpoints,
                    'original_inventory_sha256': sha(inventory),
                    'original_retrieval_manifest_sha256': sha(retrieval_manifest),
                    'relocation': {'old_prefix': old_prefix, 'new_prefix': new_prefix},
                    'new_publication': str(destination.relative_to(repo))}
        atomic_bytes(archive / 'finalization.json', json_bytes(evidence))
        proposal = {'status': 'proposal_only_no_deletion_performed',
                    'candidates': [{'path': str((round1_record.parent / 'tmp/round1').relative_to(repo)),
                                    'reason': 'obsolete exploratory fits; retain until archived evidence and old replay verified'}],
                    'retain': ['shared corpora', 'full-axis MLP inputs', 'new 8-run checkpoints',
                               'historical predictions bound by shared evaluator', 'unrelated work']}
        atomic_bytes(archive / 'cleanup-proposal.json', json_bytes(proposal))
        old_moved = new_moved = False
        try:
            if any(path.read_bytes() != payload for path, payload in originals.items()):
                raise ValueError('active evaluator metadata changed during preparation')
            destination.rename(archive / 'round1')
            old_moved = True
            if tree_hashes(archive / 'round1') != historical:
                raise ValueError('historical evidence changed during relocation')
            if not (archive / 'round1/labels/test/genes.csv').is_file():
                raise ValueError('archived significance gene family is missing')
            staging.rename(destination)
            new_moved = True
            for path, payload in replacements.items():
                atomic_bytes(path, payload)
            replay(destination)
            if tree_hashes(destination / 'checkpoints') != checkpoints:
                raise ValueError('selected checkpoint identities changed during installation')
            evidence['status'] = 'complete'
            evidence['relocated_inventory_sha256'] = sha(inventory)
            evidence['relocated_retrieval_manifest_sha256'] = sha(retrieval_manifest)
            atomic_bytes(archive / 'finalization.json', json_bytes(evidence))
        except BaseException:
            for path, payload in originals.items():
                atomic_bytes(path, payload)
            if new_moved:
                destination.rename(staging)
            if old_moved:
                (archive / 'round1').rename(destination)
            evidence['status'] = 'rolled_back'
            atomic_bytes(archive / 'finalization.json', json_bytes(evidence))
            raise
    return evidence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--staging', type=Path, required=True)
    args = parser.parse_args()
    finalize(REPO, args.staging, STEP / 'output/round1',
             STEP / 'execution/history/round1-exploratory',
             REPO / 'stage2/evaluation/inventory.json', STEP / 'output/retrieval/manifest.json',
             STEP / 'execution/round1_record')


if __name__ == '__main__':
    main()
