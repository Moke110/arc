"""Public prepare/train/evaluate/verify workflow for the frozen Morgan baseline."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import random
import subprocess
import sys

SPLIT_SEED = 20261105
TRAINING_SEEDS = [17, 42, 2026]
CLASSES = [-1, 0, 1]


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def digest(value):
    return hashlib.sha256(encoded(value)).hexdigest()


def read_rows(source):
    """Accept one CSV fixture or the existing corpus release's condition CSVs."""
    from rdkit import Chem
    if source.is_dir():
        verifier = Path(__file__).resolve().parents[2] / '01_corpus/execution/corpus.py'
        result = subprocess.run([sys.executable, str(verifier), 'verify', '--output', str(source)],
                                capture_output=True, text=True)
        if result.returncode:
            raise ValueError('Frozen corpus verification failed: ' + result.stderr.strip())
    paths = sorted(source.glob('*/corpus.csv')) if source.is_dir() else [source]
    if not paths:
        raise ValueError('No condition corpus.csv inputs found')
    rows = []
    identities = set()
    for path in paths:
        with path.open(newline='') as handle:
            for record in csv.DictReader(handle):
                row = dict(record)
                dose = float(row['dose_um'])
                if not math.isfinite(dose) or dose <= 0:
                    raise ValueError('dose_um must be positive and finite')
                row['dose_um'] = dose
                direction = float(row['direction'])
                if direction not in CLASSES:
                    raise ValueError('direction must be one of -1, 0, +1')
                row['direction'] = int(direction)
                mol = Chem.MolFromSmiles(row['canonical_smiles'])
                if mol is None:
                    raise ValueError(f"Invalid parent structure: {row['canonical_smiles']}")
                row['canonical_smiles'] = Chem.MolToSmiles(mol, isomericSmiles=True)
                row['parent_id'] = Chem.MolToInchiKey(mol)
                row['measurement_id'] = digest([row[k] for k in ['condition', 'phenotype', 'source_asset', 'source_cell']])
                if row['measurement_id'] in identities:
                    raise ValueError('Duplicate source measurement identity')
                identities.add(row['measurement_id'])
                rows.append(row)
    if not rows:
        raise ValueError('Corpus contains no measurements')
    return sorted(rows, key=lambda row: row['measurement_id'])


def cold_drug(rows):
    parents = sorted({row['parent_id'] for row in rows})
    if len(parents) < 6:
        raise ValueError('cold_drug requires at least six parents for nonempty train/validation/test roles')
    random.Random(SPLIT_SEED).shuffle(parents)
    splits = []
    for index in range(3):
        test = set(parents[index::3])
        available = sorted(set(parents) - test)
        random.Random(SPLIT_SEED + index + 1).shuffle(available)
        count = max(1, math.ceil(0.2 * len(available)))
        validation = set(available[:count])
        mapping = {p: 'test' if p in test else 'validation' if p in validation else 'train' for p in parents}
        roles = [mapping[row['parent_id']] for row in rows]
        counts = {}
        for role in ['train', 'validation', 'test']:
            selected = [row for row, assigned in zip(rows, roles) if assigned == role]
            counts[role] = dict(measurements=len(selected), parents=len({r['parent_id'] for r in selected}),
                                combinations=len({(r['parent_id'], r['dose_um']) for r in selected}),
                                classes={str(c): sum(r['direction'] == c for r in selected) for c in CLASSES})
        # Each independent target must have evidence in every role.
        for target in sorted({(r['condition'], r['phenotype']) for r in rows}):
            present = {role for row, role in zip(rows, roles) if (row['condition'], row['phenotype']) == target}
            if present != {'train', 'validation', 'test'}:
                raise ValueError(f'cold_drug fold {index} cannot populate all roles for {target}; add parent coverage')
        splits.append(dict(setting='cold_drug', fold=index, parent_roles=mapping, roles=roles, counts=counts))
    return splits


SETTINGS = ['cold_drug', 'cold_dose', 'warm_start']


def reconstruct_splits(rows, settings):
    if not settings or settings != [s for s in SETTINGS if s in settings]:
        raise ValueError('Settings must be unique and in canonical cold_drug, cold_dose, warm_start order')
    splits = []
    for setting in settings:
        if setting == 'cold_drug':
            splits.extend(cold_drug(rows))
        elif setting == 'cold_dose':
            from cold_dose import cold_dose
            splits.extend(cold_dose(rows))
        else:
            from warm_start import warm_start
            splits.extend(warm_start(rows))
    return splits


def prepare(source, output, settings=None):
    import rdkit
    rows = read_rows(source)
    settings = [s for s in SETTINGS if s in (settings or ['cold_drug'])]
    splits = reconstruct_splits(rows, settings)
    expected = [dict(condition=condition, phenotype=phenotype, setting=s['setting'], fold=s['fold'], seed=seed)
                for condition, phenotype in sorted({(r['condition'], r['phenotype']) for r in rows})
                for s in splits for seed in TRAINING_SEEDS]
    payload = dict(schema_version=1, scope='+'.join(settings), settings=settings, split_seed=SPLIT_SEED, training_seeds=TRAINING_SEEDS,
                   class_order=CLASSES, rows=rows, input_fingerprint=digest(rows), splits=splits,
                   expected_runs=expected, features=dict(radius=2, bits=2048, chirality=True, rdkit_version=rdkit.__version__,
                                 effective_configuration=dict(countSimulation=False, useBondTypes=True,
                                     onlyNonzeroInvariants=False, includeRingMembership=True,
                                     includeRedundantEnvironments=False)))
    if source.is_dir():
        manifest_path = source / 'manifest.json'
        manifest = json.loads(manifest_path.read_text())
        payload['provenance_metadata'] = dict(source='frozen Heliste', authenticated=True,
            manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
            coverage=manifest.get('coverage', {}))
    else:
        payload['provenance_metadata'] = dict(source='synthetic CSV fixture', authenticated=False)
    payload['split_fingerprint'] = digest(splits)
    target = output / 'prepared.json'
    if output.exists():
        if not target.exists() or target.read_bytes() != encoded(payload) + b'\n':
            raise ValueError('Output exists with incompatible artifacts; choose a fresh output directory')
        return
    output.mkdir(parents=True)
    target.write_bytes(encoded(payload) + b'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    command = commands.add_parser('prepare')
    command.add_argument('--input', type=Path, required=True)
    command.add_argument('--output', type=Path, required=True)
    command.add_argument('--settings', nargs='+', choices=SETTINGS, default=['cold_drug'])
    command = commands.add_parser('train')
    command.add_argument('--prepared', type=Path, required=True)
    command.add_argument('--output', type=Path, required=True)
    command.add_argument('--condition', required=True)
    command.add_argument('--phenotype', required=True)
    command.add_argument('--setting', choices=SETTINGS, required=True)
    command.add_argument('--fold', type=int, required=True)
    command.add_argument('--seed', type=int, choices=TRAINING_SEEDS, required=True)
    command.add_argument('--smoke-epochs', type=int)
    command = commands.add_parser('benchmark')
    command.add_argument('--prepared', type=Path, required=True)
    command.add_argument('--output', type=Path, required=True)
    command.add_argument('--smoke-epochs', type=int)
    for name in ['evaluate', 'verify']:
        command = commands.add_parser(name)
        group = command.add_mutually_exclusive_group(required=True)
        group.add_argument('--run', type=Path)
        group.add_argument('--benchmark', type=Path)
        command.add_argument('--prepared', type=Path)
        command.add_argument('--require-complete', action='store_true')
    args = parser.parse_args()
    try:
        if args.command == 'prepare':
            prepare(args.input, args.output, args.settings)
        elif args.command == 'benchmark':
            from benchmark import benchmark
            result = benchmark(args.prepared, args.output, args.smoke_epochs)
            print(json.dumps(result, sort_keys=True, allow_nan=False))
        elif args.command in ['evaluate', 'verify'] and args.benchmark:
            from benchmark import evaluate_benchmark
            result = evaluate_benchmark(args.benchmark, args.prepared, args.require_complete)
            print(json.dumps(result, sort_keys=True, allow_nan=False))
        else:
            from training import train, evaluate, verify
            if args.command == 'train':
                result = train(args.prepared, args.output, args.condition, args.phenotype,
                               args.setting, args.fold, args.seed, smoke_epochs=args.smoke_epochs)
            elif args.command == 'evaluate':
                result = evaluate(args.run)
            else:
                result = verify(args.run, args.prepared)
            print(json.dumps(result, sort_keys=True, allow_nan=False))
    except (ValueError, KeyError, OSError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
