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


def _recover_batch(root, path, record, acceptance_pin, runs, prior_entries, predecessor, verify_retained):
    status = record.get('status')
    if ((record.get('schema'), record.get('acceptance_sha256'), record.get('previous'))
            != (SCHEMA, acceptance_pin, predecessor) or status not in ('pending', 'complete')):
        raise ReleaseError('disposition requires an unused receipt path or the same pending transaction')
    quarantine = _inside(root, record['quarantine'])
    if (quarantine.parent != path.parent or not quarantine.name.startswith(f'.{path.name}-weights-')
            or (quarantine.exists() and not quarantine.is_dir())):
        raise ReleaseError('invalid pending checkpoint quarantine')
    prior = {entry['manifest']: entry for entry in prior_entries}
    recorded = {entry['manifest']: entry for entry in record['weights']}
    if len(recorded) != len(record['weights']) or recorded.keys() != prior.keys() | runs.keys():
        raise ReleaseError('pending checkpoint disposition roster differs')
    if any(recorded[name] != entry for name, entry in prior.items()):
        raise ReleaseError('pending checkpoint disposition predecessor differs')
    moves = []
    stage_names = set()
    for index, (_, checkpoint, expected) in enumerate(_records(root, runs)):
        entry = recorded[expected['manifest']]
        if (type(entry.get('size')) is not int or entry['size'] <= 0
                or entry != {**expected, 'size': entry['size']}):
            raise ReleaseError('pending checkpoint disposition entry differs')
        staged = quarantine / str(index)
        if staged.is_symlink() or (status == 'pending' and checkpoint.exists() == staged.exists()):
            raise ReleaseError('expected one intact copy of a pending checkpoint')
        if status == 'complete':
            if checkpoint.exists():
                raise ReleaseError('completed disposition still has a checkpoint')
            if not staged.exists():
                continue
        copy = staged if staged.exists() else checkpoint
        _check_pin(copy, entry['sha256'], 'pending checkpoint')
        if copy.stat().st_size != entry['size']:
            raise ReleaseError('pending checkpoint size differs')
        if copy == staged:
            moves.append((staged, checkpoint))
            stage_names.add(staged.name)
    if quarantine.is_dir() and {item.name for item in quarantine.iterdir()} != stage_names:
        raise ReleaseError('unrecognized files in pending checkpoint quarantine')
    if status == 'complete':
        all_runs = {entry['manifest']: entry['manifest_sha256'] for entry in prior_entries} | runs
        _check_archived_outputs(root, all_runs, record['weights'])
        verify_retained(tuple(entry.copy() for entry in record['weights']))
        for staged, _ in moves:
            staged.unlink()
    else:
        for staged, checkpoint in moves:
            staged.rename(checkpoint)
    if quarantine.is_dir():
        quarantine.rmdir()
    if status == 'pending':
        path.unlink()
    return status == 'complete'


def dispose_checkpoints(root, *, receipt, acceptance, acceptance_sha256, runs, verify_full, verify_retained,
                        previous=None):
    """Dispose of an authorized subset of a newly accepted study's checkpoint weights.

    The caller supplies an authenticated acceptance receipt, a mapping of eligible
    run manifests to accepted digests, and the numbered study's full verifier.
    The full verifier must raise unless complete acceptance holds. The retained
    verifier independently rescores predictions and replays figures before completion.
    Receive the prospective weight list as its only argument. Use a new receipt
    path distinct from all study artifacts. A pending receipt never authorizes
    prediction-archive verification. Repeat the same request after an
    interruption to restore a pending batch and retry it, or finish a completed
    batch's quarantine cleanup. For additional batches, pass the previous
    receipt path, pinned digest and accepted prior run roster as `previous`.
    The before-disposal verifier must allow exactly those previously disposed
    weights to be absent and the retained verifier must check the cumulative
    roster.
    """
    with publication_lock(root) as root:
        path = _inside(root, receipt)
        if not path.parent.is_dir():
            raise ReleaseError('disposition requires an unused receipt path')
        acceptance = Path(acceptance)
        _check_pin(acceptance, acceptance_sha256, 'full acceptance receipt')
        prior_entries = []
        all_runs = runs
        predecessor = None
        if previous is not None:
            previous_path, previous_pin, previous_runs = previous
            record = verify_prediction_archive(root, receipt=previous_path, receipt_sha256=previous_pin,
                                               acceptance=acceptance, acceptance_sha256=acceptance_sha256,
                                               runs=previous_runs)
            if set(runs) & set(previous_runs):
                raise ReleaseError('new batch overlaps a previously disposed checkpoint')
            prior_entries = record['weights']
            all_runs = previous_runs | runs
            predecessor = {'receipt': previous_path, 'sha256': previous_pin}
        if path.exists() and _recover_batch(
            root, path, read_json(path), acceptance_sha256, runs, prior_entries, predecessor, verify_retained,
        ):
            return sha256(path)
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
        all_entries = sorted(prior_entries + entries, key=lambda entry: entry['manifest'])
        record = dict(schema=SCHEMA, status='pending', acceptance_sha256=acceptance_sha256,
                      weights=all_entries, previous=predecessor,
                      quarantine=quarantine.relative_to(root).as_posix())
        try:
            _write_receipt(path, record)
            for index, entry in enumerate(entries):
                _inside(root, entry['path']).rename(quarantine / str(index))
            _check_archived_outputs(root, all_runs, all_entries)
            verify_retained(tuple(entry.copy() for entry in all_entries))
            _check_pin(acceptance, acceptance_sha256, 'full acceptance receipt')
            _check_archived_outputs(root, all_runs, all_entries)
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
