#!/usr/bin/env python3
"""Stage/replay a coherent ADR 0073 candidate without touching accepted products.

The workspace mirrors repository-relative input names, so every resulting pin
also resolves after acceptance. Input snapshots use hard links and are read-only
by convention; publication writes only newly built candidate product directories.
TOX is consumed from an independently measured and verified candidate.
"""
import argparse
import errno
import json
import os
from pathlib import Path
import shutil
from types import SimpleNamespace

import public_api
import build_disease_deg as deg_driver
import build_disease_fusion as fusion_driver
import build_gene_axis as axis_driver
import disease_deg_sources
from bases import disease_deg, gene_axis, release, view

NON_TOX = ('nHF', 'ICM', 'AF', 'genetic_cardiomyopathy')
STEP = Path('stage2/01_bases')


def link_or_copy(source, target):
    try:
        os.link(source, target)
    except OSError as error:
        if error.errno != errno.EXDEV:
            raise
        shutil.copyfile(source, target)


def snapshot(source, target):
    """Snapshot input bytes, using hard links when the filesystem permits."""
    source, target = Path(source), Path(target)
    if target.exists():
        raise ValueError(f'candidate input already exists: {target}')
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.is_dir():
        shutil.copytree(source, target, copy_function=link_or_copy)
    else:
        link_or_copy(source, target)


def sources():
    categories = {key: disease_deg_sources.CATEGORIES[key] for key in NON_TOX}
    categories['TOX'] = {'directory': 'tox', 'tier': 2, 'evidence_type': 'bulk'}
    return SimpleNamespace(CATEGORIES=categories,
                           datasets=lambda category: (disease_deg_sources.datasets(category)
                                                      if category != 'TOX' else []))


def prepare(repo, inputs, candidate, tox_candidate):
    if candidate.exists():
        raise ValueError('candidate exists; verify it or select a fresh candidate path')
    candidate.mkdir(parents=True)
    # Source/local data are immutable snapshots, never build destinations.
    for category in NON_TOX:
        for dataset in disease_deg_sources.datasets(category):
            item = Path(dataset['source_item'])
            for name in ('manifest.sha256', *(dataset.get(key) for key in
                          ('local_summary', 'local_deg', 'cells_per_donor'))):
                if name and not (candidate / item / name).exists():
                    snapshot(inputs / item / name, candidate / item / name)
    for name in ('gene-identity', 'tenv-hvg', 'corpus'):
        snapshot(inputs / STEP / 'output' / name, candidate / STEP / 'output' / name)
    # Pin the code actually used, rather than the accepted release's older code.
    shutil.copytree(repo / STEP / 'output/api', candidate / STEP / 'output/api',
                    copy_function=shutil.copyfile)
    snapshot(tox_candidate / 'tox', candidate / STEP / 'output/disease-deg/tox')
    # These snapshots must match all the measured TOX's original input pins.
    manifest = release.read_json(candidate / STEP / 'output/disease-deg/tox/manifest.json')
    release.check_manifest(candidate / STEP / 'output/disease-deg/tox', manifest)
    for name, record in manifest['inputs'].items():
        path = candidate / name
        if not path.exists():
            snapshot(repo / name if (repo / name).is_file() else inputs / name, path)
        if release.sha256(path) != record['sha256']:
            raise release.ReleaseError(f'TOX input generation mismatch: {name}')


def candidate_receipt(candidate):
    step = candidate / STEP
    fusion = disease_deg.DiseaseFusion(step / 'output/disease-fusion')
    catalog = fusion.catalog()
    roster = {row['directory']: {'n_deg': row['n_deg'], 'deg1000': row['deg1000'],
                                 'eligible': row['eligible'],
                                 'exclusion_reason': row['exclusion_reason']}
              for row in catalog['diseases']}
    if tuple(roster) != gene_axis.CATEGORIES:
        raise release.ReleaseError(f'candidate declared roster differs: {tuple(roster)}')
    axis = gene_axis.GeneAxis(step / 'output/gene-axis')
    manifest = axis.manifest()
    if manifest['fusion']['diseases'] != [name for name in roster if roster[name]['eligible']]:
        raise release.ReleaseError('axis participation differs from measured eligibility')
    previous = set()
    for size in gene_axis.AXIS_SIZES:
        current = set(axis.axis(size))
        if len(current) != size or not previous <= current:
            raise release.ReleaseError('candidate axes must be exact and nested')
        previous = current
    if any((step / f"output/view-{size}k").exists() for size in (2, 3, 4)):
        raise release.ReleaseError("retired small expression views in candidate")
    view.verify(step / "output/view-5k")
    manifests = {}
    for path in sorted((step / 'output').glob('**/manifest.json')):
        if any(part in ('api',) for part in path.parts):
            continue
        value = release.read_json(path)
        # Unchanged corpus inputs belong to the accepted corpus provenance.
        if path.parent.name not in ('corpus', 'gene-identity', 'tenv-hvg'):
            release.check_manifest(candidate if path.parent.name == 'gene-axis'
                                   else path.parent, value)
            for name, record in value.get('inputs', {}).items():
                expected = record['sha256'] if isinstance(record, dict) else record
                if release.sha256(candidate / name) != expected:
                    raise release.ReleaseError(f'mixed candidate input generation: {name}')
        manifests[path.relative_to(candidate).as_posix()] = release.sha256(path)
    return {'status': 'candidate-not-deployed', 'declared': roster,
            'fusion': manifest['fusion'], 'axes': manifest['axes'],
            'intersection': manifest['intersection'], 'manifests': manifests}


def run(action, candidate):
    step = candidate / STEP
    context = deg_driver.Context(step / 'output/disease-deg', candidate, disease_deg_sources)
    for category in NON_TOX:
        print(f'{action} disease-DEG {category}', flush=True)
        (deg_driver.publish_category if action == 'publish' else
         deg_driver.replay_category)(context, category)
    fcontext = fusion_driver.Context(step / 'output/disease-fusion', candidate, sources(),
                                     step / 'execution/checkpoints')
    (fusion_driver.publish if action == 'publish' else fusion_driver.replay)(fcontext)
    config = axis_driver.default_config(candidate)
    (axis_driver.publish if action == 'publish' else axis_driver.verify)(config)
    for size in view.VIEW_NAMES:
        print(f'{action} view {size}', flush=True)
        kwargs = dict(corpus_dir=step / 'output/corpus', axis_dir=step / 'output/gene-axis',
                      size=size, name=f'view-{size // 1000}k', repo_root=candidate)
        (view.publish_release if action == 'publish' else view.replay_release)(
            step / f'output/view-{size // 1000}k', root=step / 'execution/checkpoints', **kwargs)
    return candidate_receipt(candidate)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('publish', 'verify'))
    parser.add_argument('--candidate', required=True, type=Path)
    parser.add_argument('--inputs', type=Path)
    parser.add_argument('--tox-candidate', type=Path)
    parser.add_argument('--record', type=Path)
    args = parser.parse_args()
    repo = next(p for p in Path(__file__).resolve().parents if (p / '.git').exists())
    if args.action == 'publish':
        if args.inputs is None or args.tox_candidate is None:
            parser.error('publish requires --inputs and --tox-candidate')
        prepare(repo, args.inputs, args.candidate, args.tox_candidate)
    receipt = run(args.action, args.candidate)
    text = json.dumps(receipt, indent=2, sort_keys=True) + '\n'
    if args.record:
        args.record.parent.mkdir(parents=True, exist_ok=True)
        args.record.write_text(text)
    print(text, end='')


if __name__ == '__main__':
    main()
