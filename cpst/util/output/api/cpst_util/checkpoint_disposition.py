"""Dispose of explicitly authorized weights after a study's full acceptance."""
import json
import os
from pathlib import Path, PurePosixPath
from tempfile import NamedTemporaryFile, mkdtemp

from .derivation.release import ReleaseError, check_outputs, output_digests, read_json, sha256
from .publication import publication_lock

SCHEMA = 'cpst-checkpoint-disposition-v1'
WEIGHT = 'checkpoint.pt'


def _inside(root, name):
    if not isinstance(name, str):
        raise ReleaseError(f'invalid disposition path: {name}')
    relative = PurePosixPath(name)
    if (not relative.parts or relative.is_absolute() or '..' in relative.parts
            or relative.as_posix() != name or any('\\' in part or ':' in part for part in relative.parts)):
        raise ReleaseError(f'invalid disposition path: {name}')
    path = root
    for part in relative.parts:
        path /= part
        if path.is_symlink():
            raise ReleaseError(f'disposition path contains a symbolic link: {name}')
    return path


def _check_pin(path, expected, label):
    if not path.is_file() or sha256(path) != expected:
        raise ReleaseError(f'{label} digest differs: {path.name}')


def _records(root, runs):
    if not runs:
        raise ReleaseError('empty checkpoint disposition roster')
    for name, digest in sorted(runs.items()):
        manifest_path = _inside(root, name)
        if manifest_path.name != 'manifest.json':
            raise ReleaseError(f'expected a sealed run manifest: {name}')
        _check_pin(manifest_path, digest, 'run manifest')
        manifest = read_json(manifest_path)
        outputs = output_digests(manifest)
        if WEIGHT not in outputs or len(outputs) < 2:
            raise ReleaseError(f'run lacks weights or retained artifacts: {name}')
        for output in outputs:
            _inside(manifest_path.parent, output)
        relative = PurePosixPath(name).with_name(WEIGHT).as_posix()
        checkpoint = _inside(root, relative)
        yield manifest, checkpoint, dict(path=relative, sha256=outputs[WEIGHT],
                                         manifest=name, manifest_sha256=digest)


def _write_receipt(path, value):
    temporary = NamedTemporaryFile(mode='w', encoding='utf-8', dir=path.parent, delete=False)
    candidate = Path(temporary.name)
    try:
        with temporary as stream:
            json.dump(value, stream, indent=2, sort_keys=True)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        candidate.replace(path)
    finally:
        candidate.unlink(missing_ok=True)


def dispose_checkpoints(root, *, receipt, acceptance, acceptance_sha256, runs, verify_full, verify_retained):
    """Dispose of an authorized subset of a newly accepted study's checkpoint weights.

    The caller supplies an authenticated acceptance receipt, a mapping of eligible
    run manifests to accepted digests, and the numbered study's full verifier.
    The full verifier must raise unless complete acceptance holds. The retained
    verifier independently rescores predictions and replays figures before completion.
    Receive the prospective weight list as its only argument. Use a new receipt
    path distinct from all study artifacts. A pending receipt never authorizes
    prediction-archive verification.
    """
    with publication_lock(root) as root:
        path = _inside(root, receipt)
        if path.exists() or not path.parent.is_dir():
            raise ReleaseError('disposition requires an unused receipt path')
        acceptance = Path(acceptance)
        _check_pin(acceptance, acceptance_sha256, 'full acceptance receipt')
        verify_full()
        _check_pin(acceptance, acceptance_sha256, 'full acceptance receipt')
        entries = []
        protected = {acceptance.resolve()}
        for manifest, checkpoint, entry in _records(root, runs):
            check_outputs(checkpoint.parent, manifest)
            size = checkpoint.stat().st_size
            if size <= 0:
                raise ReleaseError(f'checkpoint is empty: {entry["path"]}')
            entries.append({**entry, 'size': size})
            protected.add(root / entry['manifest'])
            protected.update(checkpoint.parent / name for name in output_digests(manifest))
        if path in protected:
            raise ReleaseError('disposition receipt overlaps an accepted artifact')
        quarantine = Path(mkdtemp(prefix=f'.{path.name}-weights-', dir=path.parent))
        record = dict(schema=SCHEMA, status='pending', acceptance_sha256=acceptance_sha256,
                      weights=entries, quarantine=quarantine.relative_to(root).as_posix())
        try:
            _write_receipt(path, record)
            for index, entry in enumerate(entries):
                _inside(root, entry['path']).rename(quarantine / str(index))
            _check_archived_outputs(root, runs, entries)
            verify_retained(tuple(entry.copy() for entry in entries))
            _check_pin(acceptance, acceptance_sha256, 'full acceptance receipt')
            _check_archived_outputs(root, runs, entries)
            _write_receipt(path, {**record, 'status': 'complete'})
        except BaseException:
            # Keep staged bytes until verification and receipt publication succeed.
            for index, entry in enumerate(entries):
                staged = quarantine / str(index)
                if staged.exists():
                    checkpoint = _inside(root, entry['path'])
                    if checkpoint.exists():
                        raise ReleaseError(f'cannot restore staged weights: {staged}')
                    staged.rename(checkpoint)
            quarantine.rmdir()
            raise
        for index in range(len(entries)):
            (quarantine / str(index)).unlink()
        quarantine.rmdir()
        return sha256(path)


def _check_archived_outputs(root, runs, entries):
    records = list(_records(root, runs))
    if not isinstance(entries, list) or len(records) != len(entries):
        raise ReleaseError('checkpoint disposition roster differs')
    for (manifest, checkpoint, expected), entry in zip(records, entries, strict=True):
        if (not isinstance(entry, dict) or type(entry.get('size')) is not int or entry['size'] <= 0
                or entry != {**expected, 'size': entry['size']} or checkpoint.exists()):
            raise ReleaseError(f'invalid completed checkpoint disposition: {expected["path"]}')
        retained = {name: digest for name, digest in output_digests(manifest).items() if name != WEIGHT}
        check_outputs(checkpoint.parent, {'outputs': retained})


def verify_prediction_archive(root, *, receipt, receipt_sha256, acceptance, acceptance_sha256, runs):
    """Check retained artifact bytes against a pinned completed disposition.

    The numbered study still owns independent prediction rescoring and figure
    replay. This verification does not load checkpoints or replay training.
    """
    root = Path(root).resolve(strict=True)
    path = _inside(root, receipt)
    _check_pin(path, receipt_sha256, 'disposition receipt')
    _check_pin(Path(acceptance), acceptance_sha256, 'full acceptance receipt')
    record = read_json(path)
    if (record.get('schema'), record.get('status'), record.get('acceptance_sha256')) != (
        SCHEMA, 'complete', acceptance_sha256,
    ):
        raise ReleaseError('checkpoint disposition is not complete or bound to this acceptance')
    _check_archived_outputs(root, runs, record.get('weights'))
    return record
