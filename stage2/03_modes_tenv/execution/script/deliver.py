#!/usr/bin/env python
"""Replay verified evaluation artifacts and publish human-facing step results.

The frozen benchmark replay remains unchanged. This presentation-only layer adds
reference/completeness labels to its best-run figure without touching scored
values, selections or the run identity. Run this after any original replay to
restore the delivery labels; publication.json records both rendering sources.
"""
import argparse
from pathlib import Path
import shutil

import numpy as np

from artifacts import read_json, rows, sha256, write_json
from benchmark import replay
from plots import plt

TABLES = ('per_group.csv', 'per_run.csv', 'summary.csv', 'best_runs.csv', 'pr_bands.csv')
FIGURES = tuple(f'{name}.{extension}' for name in ('metrics', 'pr_panels', 'best_run')
                for extension in ('png', 'svg'))


def label_best_run(output):
    summaries = {row['model']: row for row in rows(output / 'summary.csv')}
    with np.load(output / 'pr_sources.npz', allow_pickle=False) as source:
        identifiers = source['identifiers'].tolist()
        figure, axis = plt.subplots(figsize=(10, 7), constrained_layout=True)
        for selected in rows(output / 'best_runs.csv'):
            row = summaries[selected['model']]
            completed = int(row['completed_runs'])
            expected = int(row['expected_runs'])
            incomplete = '; incomplete' if completed < expected else ''
            label = f"{selected['model']} (completed {completed}/{expected}{incomplete})"
            if row['deterministic_reference'] == 'True':
                label += ' [deterministic copies]'
            if not selected['fold']:
                axis.plot([], [], label=label + ': undefined selection')
                continue
            key = f"{selected['model']}:{selected['fold']}:{selected['seed']}"
            curves = source['curves'][identifiers.index(key)]
            axis.plot(source['recall'], (curves[0] + curves[1]) / 2,
                      label=label + f"; fold {selected['fold']}, seed {selected['seed']}")
        axis.set(xlabel='Recall', ylabel='Mean directional precision', xlim=(0, 1), ylim=(0, 1.01),
                 title='Post-hoc best test run per entry\nDifferent entries may select different folds; nine-run summaries are primary')
        axis.legend(fontsize=8)
        for extension in ('png', 'svg'):
            figure.savefig(output / f'best_run.{extension}', **({'dpi': 150} if extension == 'png' else {}))
        plt.close(figure)


def publish(output, destination):
    output = Path(output).resolve()
    destination = Path(destination).resolve()
    if output == destination:
        raise ValueError('publication destination must differ from evaluation output')
    # The original replay checks every scored source, including prediction barrier
    # metadata, before any delivery files are written. It never refits or scores.
    replay(output)
    label_best_run(output)
    destination.mkdir(parents=True, exist_ok=True)
    for name in TABLES + FIGURES:
        shutil.copyfile(output / name, destination / name)
    evaluation = read_json(output / 'evaluation.json')
    script = Path(__file__).resolve()
    manifest = {
        'evaluation_sha256': sha256(output / 'evaluation.json'),
        'source_artifacts': evaluation['artifacts'],
        'rendering_sources': {name: sha256(script.with_name(name))
                              for name in ('plots.py', 'deliver.py')},
        'presentation': 'Frozen artifact-only replay; best-run labels add completion, incomplete and deterministic-copy disclosures. No scoring or selection changes.',
        'artifacts': {name: sha256(destination / name) for name in TABLES + FIGURES},
    }
    write_json(destination / 'publication.json', manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='verified benchmark evaluation directory')
    parser.add_argument('--publish', type=Path, required=True, help='human-facing step result directory')
    args = parser.parse_args()
    publish(args.output, args.publish)


if __name__ == '__main__':
    main()
