#!/usr/bin/env python
"""Replay verified evaluation artifacts and publish human-facing step results.

The frozen benchmark replay remains unchanged. This presentation-only layer adds
reference/completeness labels to its best-run figure without touching scored
values, selections or the run identity. Run this after any original replay to
restore the delivery labels; publication.json records both rendering sources.
"""
import argparse
import csv
import os
import re
from pathlib import Path
import shutil
import tempfile
from xml.sax.saxutils import escape

import numpy as np

from artifacts import read_json, rows, sha256, write_json
from benchmark import replay
from plots import plt

TABLES = ('per_group.csv', 'per_run.csv', 'summary.csv', 'best_runs.csv', 'pr_bands.csv')
FIGURES = tuple(f'{name}.{extension}' for name in ('metrics', 'pr_panels', 'best_run')
                for extension in ('png', 'svg'))


def label_best_run(output, destination=None):
    destination = output if destination is None else destination
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
            figure.savefig(destination / f'best_run.{extension}', **({'dpi': 150} if extension == 'png' else {}))
        plt.close(figure)


def disclosure_lines(proof):
    original = proof['original_timeout']['budget']
    selected = proof['selected_attempt']
    account = selected['account']
    timing = account['timing']

    def measured(value, unit):
        return 'unavailable' if value is None else f'{value:,.2f} {unit}'

    lines = (
        f"CMonge v2 fold 0 seed 17 · original cumulative-budget timeout "
        f"({original['remaining_at_launch_seconds']:,.2f} s left of 14,400 s; 14,400 s debited).",
        f"One fresh 14,400 s training authorization from seed initialization + 1,800 s inference. "
        f"Selected attempt: {selected['status']}.",
        f"Original prior debit {original['prior_debit_seconds']:,.2f} s · "
        f"original worker wall {original['attempt_wall_seconds']:,.2f} s · "
        f"original host peak RSS {measured(original['peak_host_rss_gib'], 'GiB')}.",
        f"Selected watchdog setup/training {measured(timing['pre_inference_watchdog_seconds'], 's')} · "
        f"observed inference {measured(timing['observed_inference_seconds'], 's')}.",
        f"Selected host peak RSS {measured(account['peak_host_rss_gib'], 'GiB')} "
        "(maximum of available owned-process samples and GNU time peak).",
    )
    if 'sampled_peak_host_rss_gib' in account:
        lines += (f"Owned-process RSS samples {measured(account['sampled_peak_host_rss_gib'], 'GiB')} · "
                  f"GNU time peak RSS {measured(account['time_peak_rss_gib'], 'GiB')}.",)
    return lines


def add_figure_disclosure(directory, name, lines):
    """Append the same disclosure to both formats without changing the plots."""
    from PIL import Image, ImageDraw, ImageFont
    from matplotlib import font_manager

    png = directory / f'{name}.png'
    with Image.open(png) as source:
        original = source.convert('RGB')
    margin = 18
    font = ImageFont.truetype(font_manager.findfont('DejaVu Sans'), size=21)
    if max(font.getlength(line) for line in lines) > original.width - 2 * margin:
        raise ValueError('figure is too narrow for the budget disclosure')
    strip = 2 * margin + 27 * len(lines)
    image = Image.new('RGB', (original.width, original.height + strip), '#fff2db')
    image.paste(original, (0, 0))
    draw = ImageDraw.Draw(image)
    for i, line in enumerate(lines):
        draw.text((margin, original.height + margin + 27 * i), line, font=font, fill='#532908')
    image.save(png, dpi=(150, 150))

    svg = directory / f'{name}.svg'
    text = svg.read_text()
    root = re.search(r'<svg [^>]+>', text)
    if not root:
        raise ValueError('expected an SVG root')
    viewbox = re.search(r'viewBox="0 0 ([0-9.]+) ([0-9.]+)"', root[0])
    if not viewbox:
        raise ValueError('expected a zero-origin SVG viewBox')
    width, height = map(float, viewbox.groups())
    if max(font.getlength(line) for line in lines) * .5 > width - 20:
        raise ValueError('SVG is too narrow for the budget disclosure')
    extra = 20 + 15 * len(lines)
    header = re.sub(r'height="[0-9.]+pt"', f'height="{height + extra}pt"', root[0])
    header = header.replace(viewbox[0], f'viewBox="0 0 {width} {height + extra}"')
    annotation = f'<rect x="0" y="{height}" width="{width}" height="{extra}" fill="#fff2db"/>\n'
    for i, line in enumerate(lines):
        annotation += (f'<text x="10" y="{height + 18 + 15 * i}" fill="#532908" '
                       f'font-family="DejaVu Sans, sans-serif" font-size="10.5">{escape(line)}</text>\n')
    svg.write_text(text[:root.start()] + header + text[root.end():].replace('</svg>', annotation + '</svg>'))


def publish_supplement_tables(output, destination, proof):
    # Published CSVs retain the scorer's row order and original field strings.
    # The underlying five scorer tables remain frozen in the overlay.
    disclosure = ' '.join(disclosure_lines(proof))
    attempt = proof['selected_attempt']
    for name in TABLES:
        with (output / name).open(newline='') as source, (destination / name).open('w', newline='') as target:
            reader = csv.DictReader(source)
            fields = [*reader.fieldnames, 'cmonge_budget_exception', 'selected_attempt_id',
                      'original_attempt_status', 'supplement_attempt_status']
            writer = csv.DictWriter(target, fieldnames=fields)
            writer.writeheader()
            for row in reader:
                if row['model'] == 'cmonge':
                    row['cmonge_budget_exception'] = disclosure
                    if row.get('fold') == '0' and row.get('seed') == '17':
                        row.update(selected_attempt_id=attempt['attempt_id'],
                                   original_attempt_status='timeout', supplement_attempt_status=attempt['status'])
                writer.writerow(row)
    write_json(destination / 'slot_provenance.json', proof)

def publish(output, destination, *, supplement=None):
    output = Path(output).resolve()
    destination = Path(destination).resolve()
    if output == destination:
        raise ValueError('publication destination must differ from evaluation output')
    if supplement is None and (output / 'selection.json').exists():
        raise ValueError('use the authenticated supplement overlay publication command')
    if supplement is not None:
        return publish_supplement(output, destination, supplement)
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


def publish_supplement(output, destination, supplement):
    from learning_curves import render
    from supplement_overlay import checked_path, evaluation_proof

    record, proof, evaluation_pin = supplement
    if proof != evaluation_proof(output, record, proof['selection_sha256']):
        raise ValueError('evaluation evidence changed before replay')
    for name in FIGURES:
        if os.path.lexists(output / name):
            checked_path(output / name)
    replay(output)
    names = TABLES + FIGURES + ('slot_provenance.json', 'learning_curves.png',
                               'learning_curves.svg', 'learning_curves.json')
    destination.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.cmonge-publish-', dir=destination) as temporary:
        staged = Path(temporary)
        for name in FIGURES:
            shutil.copyfile(output / name, staged / name)
        label_best_run(output, staged)
        for name in ('metrics', 'pr_panels', 'best_run'):
            add_figure_disclosure(staged, name, disclosure_lines(proof))
        publish_supplement_tables(output, staged, proof)
        render(output, staged, final=True, supplement=supplement)
        if proof != evaluation_proof(output, record, proof['selection_sha256']):
            raise ValueError('scored artifacts changed during replay or presentation')
        script = Path(__file__).resolve()
        manifest = dict(evaluation_sha256=proof['evaluation_sha256'],
                        source_artifacts=proof['artifacts'], selection_sha256=proof['selection_sha256'],
                        evaluation_provenance_sha256=evaluation_pin,
                        presentation='Frozen artifact-only replay with supplemental budget disclosures.',
                        rendering_sources={name: sha256(script.with_name(name)) for name in
                                           ('plots.py', 'deliver.py', 'learning_curves.py', 'supplement_overlay.py')},
                        artifacts={name: sha256(staged / name) for name in names})
        write_json(staged / 'publication.json', manifest)
        for name in (*names, 'publication.json'):
            (staged / name).replace(destination / name)
    return manifest

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True, help='verified benchmark evaluation directory')
    parser.add_argument('--publish', type=Path, required=True, help='human-facing step result directory')
    args = parser.parse_args()
    publish(args.output, args.publish)


if __name__ == '__main__':
    main()
