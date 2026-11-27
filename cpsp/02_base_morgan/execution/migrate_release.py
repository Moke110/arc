"""Reissue the sealed Morgan release using the installed CPSP corpus identity."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from workflow import digest, encoded

ROOT = Path(__file__).resolve().parents[3]
STAGE = ROOT / 'cpsp'
EXECUTION = STAGE / '02_base_morgan/execution'
CORPUS = STAGE / '01_corpus/output'
OUTPUT = STAGE / '02_base_morgan/output'
RECORD = EXECUTION / 'release_record.json'
ORIGINAL_STAGE = 'stage3'
ORIGINAL_ROOT = Path('/home/chang/arc') / ORIGINAL_STAGE
ORIGINAL_OUTPUT = ORIGINAL_ROOT / '02_base_morgan/output'
ORIGINAL_REVISION = '2fd12768b755afe405f14a41f726e484c2f9812d'
INITIAL_REVISION = '844d849fdf7b25c39eae2c198403146dfbea0c1f'
INITIAL_RESOURCES_REVISION = '49f6839a6b7fe1f8ec128f751275814f0ee94ff1'
PREPARED = 'prepared/prepared.json'
LEDGER = 'formal/benchmark.json'
SUMMARY = 'formal/summary.json'
REPORT = 'formal/report.md'


def sha(payload):
    return hashlib.sha256(payload).hexdigest()


def check(actual, expected, label):
    if actual != expected:
        raise ValueError(f'{label} differs from the sealed release')


def git_bytes(revision, path):
    result = subprocess.run(['git', 'show', f'{revision}:{path}'], cwd=ROOT, capture_output=True)
    if result.returncode:
        raise ValueError(f'Sealed Git object unavailable: {revision}:{path}')
    return result.stdout


def formatted(data):
    return (json.dumps(data, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + '\n').encode()


def read_verified(path, expected):
    data = path.read_bytes()
    check(sha(data), expected, str(path))
    return data


def corpus_proof(original_record):
    original_path = ORIGINAL_STAGE + '/01_corpus/output/manifest.json'
    original = git_bytes(ORIGINAL_REVISION, original_path)
    original_sha = sha(original)
    check(original_sha, original_record['source_provenance']['manifest_sha256'], original_path)
    manifest = json.loads(original)
    for name, expected in manifest['code'].items():
        check(sha(git_bytes(ORIGINAL_REVISION, name)), expected, name)
    for name, expected in manifest['inputs'].items():
        read_verified(ROOT / name, expected)
    for name, expected in manifest['files'].items():
        read_verified(CORPUS / name, expected)

    name = ORIGINAL_STAGE + '/01_corpus/execution/corpus.py'
    migrated_name = 'cpsp' + name[len(ORIGINAL_STAGE):]
    manifest['code'][migrated_name] = manifest['code'].pop(name)
    for name, expected in manifest['code'].items():
        read_verified(ROOT / name, expected)
    current = CORPUS / 'manifest.json'
    check(current.read_bytes(), formatted(manifest), str(current))
    return original_sha, sha(current.read_bytes())


def code_proof(original):
    initial = original['initial_code_sha256']
    transition = original['code_transitions']
    check(initial, transition['previous_hashes'], 'Initial code transition')
    for name, expected in initial.items():
        revision = INITIAL_RESOURCES_REVISION if name.endswith('/resources.py') else INITIAL_REVISION
        check(sha(git_bytes(revision, name)), expected, name)
    for name, expected in transition['continuation_hashes'].items():
        check(sha(git_bytes(transition['commit'], name)), expected, name)
    for additional in transition['additional_transitions']:
        name = ORIGINAL_STAGE + '/02_base_morgan/execution/resources.py'
        check(sha(git_bytes(additional['commit'], name)), additional['resources_sha256'], name)

    current = {}
    for name, expected in original['current_code_sha256'].items():
        check(sha(git_bytes(ORIGINAL_REVISION, name)), expected, name)
        current_name = 'cpsp' + name[len(ORIGINAL_STAGE):]
        current[current_name] = sha((ROOT / current_name).read_bytes())
    return current


def original_metadata(output, original, migrated, corpus_shas):
    paths = [PREPARED, LEDGER, SUMMARY]
    payload = {name: (output / name).read_bytes() for name in paths}
    prepared, ledger, summary = (json.loads(payload[name]) for name in paths)
    source_sha = corpus_shas[1] if migrated else corpus_shas[0]
    check(prepared['provenance_metadata']['manifest_sha256'], source_sha, 'Prepared corpus identity')
    check(summary['source_coverage']['manifest_sha256'], source_sha, 'Summary corpus identity')
    check(ledger['prepared_fingerprint'], digest(prepared), 'Ledger prepared fingerprint')
    prepared_path = '../' + PREPARED if migrated else str(ORIGINAL_OUTPUT / PREPARED)
    check(ledger['prepared_path'], prepared_path, 'Ledger prepared path')
    for entry in ledger['runs']:
        subpath = Path('runs') / digest(entry['identity'])
        expected = str(subpath if migrated else ORIGINAL_OUTPUT / 'formal' / subpath)
        check(entry['run_path'], expected, 'Ledger run path')
        entry['run_path'] = str(ORIGINAL_OUTPUT / 'formal' / subpath)
    prepared['provenance_metadata']['manifest_sha256'] = corpus_shas[0]
    summary['source_coverage']['manifest_sha256'] = corpus_shas[0]
    ledger['prepared_path'] = str(ORIGINAL_OUTPUT / PREPARED)
    ledger['prepared_fingerprint'] = digest(prepared)

    records = [prepared, ledger, summary]
    old = {path: encoded(data) + b'\n' for path, data in zip(paths, records)}
    check(digest(prepared), original['prepared_fingerprint'], 'Prepared artifact')
    for path in (LEDGER, SUMMARY):
        check(sha(old[path]), original['aggregate_sha256'][Path(path).name], path)
    return payload, old, records


def verify_runs(output, original, ledger):
    expected = {record['directory']: record for record in original['runs']}
    directory = output / 'formal/runs'
    check(len(expected), 60, 'Formal run count')
    check({path.name for path in directory.iterdir()}, set(expected), 'Formal run directories')
    check(len(ledger['runs']), 60, 'Ledger run count')
    for entry in ledger['runs']:
        key = digest(entry['identity'])
        record = expected[key]
        check(entry['identity'], record['identity'], f'{key} run identity')
        check(entry['status'], 'completed', f'{key} run status')
        for name, expected_sha in record['artifact_sha256'].items():
            read_verified(directory / key / name, expected_sha)


def reissue(output, migrated):
    original_path = ORIGINAL_STAGE + '/02_base_morgan/execution/release_record.json'
    original_bytes = git_bytes(ORIGINAL_REVISION, original_path)
    original = json.loads(original_bytes)
    current_code = code_proof(original)
    split = ORIGINAL_STAGE + '/02_base_morgan/execution/split_manifest.json'
    check(sha(git_bytes(ORIGINAL_REVISION, split)), original['split_manifest_sha256'], split)
    read_verified(EXECUTION / 'split_manifest.json', original['split_manifest_sha256'])
    read_verified(output / REPORT, original['aggregate_sha256']['report.md'])
    corpus_shas = corpus_proof(original)
    actual, old, (prepared, ledger, summary) = original_metadata(output, original, migrated, corpus_shas)
    verify_runs(output, original, ledger)

    prepared['provenance_metadata']['manifest_sha256'] = corpus_shas[1]
    summary['source_coverage']['manifest_sha256'] = corpus_shas[1]
    ledger['prepared_fingerprint'] = digest(prepared)
    ledger['prepared_path'] = '../' + PREPARED
    for entry in ledger['runs']:
        entry['run_path'] = 'runs/' + digest(entry['identity'])
    new = {
        PREPARED: encoded(prepared) + b'\n',
        LEDGER: encoded(ledger) + b'\n',
        SUMMARY: encoded(summary) + b'\n',
    }
    check(actual, new if migrated else old, 'On-disk metadata')

    record = dict(original)
    for key in ('code_transitions', 'initial_code_sha256', 'public_replay'):
        del record[key]
    record.update(
        stage='cpsp',
        current_code_sha256=current_code,
        aggregate_sha256={Path(path).name: sha(new[path]) for path in (SUMMARY, LEDGER)}
        | {'report.md': original['aggregate_sha256']['report.md']},
        source_provenance=prepared['provenance_metadata'],
        prepared_fingerprint=ledger['prepared_fingerprint'],
        migration=dict(
            original_release=dict(revision=ORIGINAL_REVISION, path=original_path, sha256=sha(original_bytes)),
            original_code_revisions=dict(initial=INITIAL_REVISION, initial_resources=INITIAL_RESOURCES_REVISION),
            metadata={
                path: dict(original_sha256=sha(old[path]), migrated_sha256=sha(new[path])) for path in old
            } | {
                'cpsp/01_corpus/output/manifest.json': dict(
                    original_path=ORIGINAL_STAGE + '/01_corpus/output/manifest.json',
                    original_sha256=corpus_shas[0], migrated_sha256=corpus_shas[1]),
            },
        ),
    )
    record['summary'] = dict(record['summary'], source_coverage=summary['source_coverage'])
    return record, new


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    migrate = commands.add_parser('migrate', help='stage a new output directory from sealed original outputs')
    migrate.add_argument('--source', type=Path, required=True, help='directory containing original prepared/ and formal/')
    migrate.add_argument('--output', type=Path, default=OUTPUT)
    verify = commands.add_parser('verify', help='compare installed outputs and the active release record against sealed originals')
    verify.add_argument('--output', type=Path, default=OUTPUT)
    verify.add_argument('--record', type=Path, default=RECORD)
    args = parser.parse_args()
    try:
        if args.command == 'verify':
            record, _ = reissue(args.output, migrated=True)
            check(args.record.read_bytes(), formatted(record), str(args.record))
            print('CPSP release proof: 60 formal runs and original source identities preserved')
            return

        output = args.output.resolve()
        if output.exists() or output.is_relative_to(args.source.resolve()):
            raise ValueError('Choose a fresh output directory outside the source tree')
        record, metadata = reissue(args.source, migrated=False)
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
            staging = Path(temporary) / 'output'
            shutil.copytree(args.source, staging, symlinks=True)
            for name, payload in metadata.items():
                path = staging / name
                if path.is_symlink() or not path.resolve().is_relative_to(staging):
                    raise ValueError(f'Metadata must be a regular file: {name}')
                path.write_bytes(payload)
            check(reissue(staging, migrated=True)[0], record, 'Staged metadata')
            staging.rename(output)
        print(formatted(record).decode(), end='')
    except (ValueError, OSError, KeyError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
