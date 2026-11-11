"""Recompute and replay the three canonical Stage 2 retrieval metrics."""
import argparse
from collections import defaultdict
import csv
from functools import lru_cache
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, str(next(parent for parent in Path(__file__).resolve().parents if (parent / '.git').exists())))
from stage2.evaluation.inventory import REPO, ROOT, STEPS, iter_bindings, load_inventory, read, sha
from stage2.evaluation.metrics import METRICS, PooledPR, retrieval

KEYS = ('step', 'dataset', 'regime', 'role', 'model', 'fold', 'seed')


def checked(path, expected):
    path = REPO / path
    if path.resolve() != path or sha(path) != expected:
        raise ValueError(f'changed or unsafe bound artifact: {path}')
    return path


def bound(value):
    return checked(value['path'], value['sha256'])


def membership(ids):
    import hashlib
    return hashlib.sha256(json.dumps(sorted(ids), separators=(',', ':')).encode()).hexdigest()


class Labels:
    def __init__(self, dataset, label_root=None):
        if dataset != 'tenv':
            raise ValueError('Stage 2 retrieval is evaluated only on tenv')
        self.root = (Path(label_root).absolute() if label_root else ROOT / 'evaluation/output/labels') / dataset
        if not (self.root / 'manifest.json').exists():
            partial = self.root.with_name(f'{dataset}-partial')
            if (partial / 'manifest.json').exists():
                self.root = partial
        self.manifest_digest = sha(self.root / 'manifest.json')
        self.manifest = read(self.root / 'manifest.json')
        if self.manifest.get('schema') != 'stage2-significance-v1' or self.manifest.get('dataset') != dataset:
            raise ValueError('significance dataset or schema differs from the Stage 2 standard')
        if self.manifest['fdr'] != .05:
            raise ValueError('significance threshold differs from the Stage 2 standard')
        for name, digest in self.manifest['files'].items():
            checked(str((self.root / name).relative_to(REPO)), digest)
        self.groups = pd.read_csv(self.root / 'groups.csv', keep_default_na=False).set_index('group_id')
        self.index = {name: i for i, name in enumerate(self.groups.index)}
        family = pd.read_csv(self.root / 'genes.csv')
        if len(family) != 18165 or family.hgnc_id.duplicated().any() or family.iloc[-1].symbol != 'HDGFL3':
            raise ValueError('invalid significance gene family')
        self.genes = family.iloc[:18164]
        self.padj = np.load(self.root / 'padj.npy', mmap_mode='r')
        self.lfc = np.load(self.root / 'log2_fold_change.npy', mmap_mode='r')
        if self.padj.shape != (len(self.groups), 18165) or self.padj.shape != self.lfc.shape:
            raise ValueError('invalid full-CM significance release')

    def verify(self):
        checked(str((self.root / 'manifest.json').relative_to(REPO)), self.manifest_digest)
        for name, digest in self.manifest['files'].items():
            checked(str((self.root / name).relative_to(REPO)), digest)

    def truth(self, group, cases):
        metadata = self.groups.loc[group]
        if metadata.case_ids_sha256 != membership(cases):
            raise ValueError(f'significance contrast membership differs: {group}')
        index = self.index[group]
        values = self.padj[index]
        finite = np.isfinite(values)
        if (np.any(np.isinf(values)) or np.any(values[finite] < 0) or np.any(values[finite] > 1)
                or metadata.status not in ('eligible', 'untestable')
                or (metadata.status == 'eligible') != bool(finite.any())):
            raise ValueError('invalid significance eligibility or adjusted p-values')
        if metadata.status != 'eligible':
            return None
        significant = values[:18164] < .05
        direction = self.lfc[index, :18164]
        if np.any(significant & (~np.isfinite(direction) | (direction == 0))):
            raise ValueError('invalid significant-gene direction')
        return np.where(significant, np.sign(direction), 0).astype(np.int8)


def validate_prediction_parent(item):
    if 'prediction_manifest' in item:
        manifest = read(bound(item['prediction_manifest']))
        expected = manifest.get('prediction_sha256', manifest.get('files', {}).get('delta.npy'))
        if expected != item['prediction']['sha256']:
            raise ValueError('prediction manifest does not bind the selected prediction')
        if 'prediction_order' in item and manifest.get('cases_sha256') != item['prediction_order']['sha256']:
            raise ValueError('prediction manifest does not bind the case order')
        for name, key in (('cases.csv', 'cases'), ('genes.csv', 'genes')):
            if name in manifest.get('files', {}) and manifest['files'][name] != item[key]['sha256']:
                raise ValueError(f'prediction manifest does not bind {key}')
        if 'selected_run' in item:
            run = read(bound(item['selected_run']))
            if manifest.get('run_sha256') != item['selected_run']['sha256']:
                raise ValueError('prediction manifest does not bind its selected run')
            for field in ('identity', 'recipe_sha256', 'checkpoint_sha256', 'best_step',
                          'best_validation_loss', 'model_state_sha256'):
                if manifest.get(field) != run.get(field):
                    raise ValueError(f'prediction differs from the selected checkpoint: {field}')
    elif 'run' in item:
        run = read(bound(item['run']))
        if (run.get('status') != 'completed'
                or run['artifacts']['prediction.npy'] != item['prediction']['sha256']
                or run['artifacts']['case_ids.json'] != item['prediction_order']['sha256']):
            raise ValueError('selected run does not bind the completed prediction')
        if any(run['identity'][key] != item[key] for key in ('model', 'fold', 'seed')):
            raise ValueError('selected run identity differs from the frozen inventory')


def mean_rows(matrix, rows):
    total = np.zeros(matrix.shape[1], dtype=np.float64)
    for start in range(0, len(rows), 32):
        batch = matrix[rows[start:start + 32]]
        if not np.isfinite(batch).all():
            raise ValueError('nonfinite predictions')
        total += batch.sum(axis=0, dtype=np.float64)
    return total / len(rows)


@lru_cache(maxsize=8)
def training_mean(cases_path, delta_path, controls_path=None):
    cases = pd.read_csv(cases_path)
    matrix = np.load(delta_path, mmap_mode='r')
    controls = np.load(controls_path, mmap_mode='r') if controls_path else None
    total = np.zeros(matrix.shape[1], dtype=float)
    groups = cases.groupby('group_id', sort=False).indices
    for rows in groups.values():
        value = mean_rows(matrix, rows)
        if controls is not None:
            value -= controls[cases.iloc[rows].control_index].mean(axis=0, dtype=np.float64)
        total += value
    return total / len(groups)


def score(entries, destination, label_root=None, figure_destination=None):
    destination.mkdir(parents=True)
    from stage2.evaluation.continuous import ContinuousMetrics, FIELDS as CONTINUOUS
    source_hashes = {name: sha(Path(__file__).with_name(name))
                     for name in ('evaluate.py', 'metrics.py', 'figures.py', 'inventory.py', 'continuous.py')}
    label_sets = {}
    pooled, reports, dose_reports = defaultdict(PooledPR), [], []
    group_fields = [*KEYS, 'group_id', 'dose_um', 'wells', 'significance_status',
                    *METRICS, *(f'{m}_status' for m in METRICS), *(f'{m}_terminal_recall' for m in METRICS), *CONTINUOUS]
    with (destination / 'per_group.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=group_fields)
        writer.writeheader()
        for number, item in enumerate(entries):
            labels = label_sets.setdefault(item['dataset'], None)
            if labels is None:
                labels = label_sets[item['dataset']] = (Labels(item['dataset'], label_root) if label_root else Labels(item['dataset']))
            cases = pd.read_csv(bound(item['cases']))
            if cases.case_id.duplicated().any():
                raise ValueError('duplicate prediction cases')
            genes = pd.read_csv(bound(item['genes']))
            columns = labels.genes.set_index('hgnc_id').index.get_indexer(genes.hgnc_id)
            covered = columns >= 0
            if np.any(~covered & ~genes.symbol.eq('HDGFL3')):
                raise ValueError('prediction contains unexpected genes outside CM axis')
            if 'reference' in item:
                vector = (np.zeros(len(genes)) if item['reference'] == 'no_change' else
                          training_mean(str(bound(item['train_cases'])), str(bound(item['train_delta'])),
                                        str(bound(item['train_controls'])) if 'train_controls' in item else None))
                matrix = None
            else:
                validate_prediction_parent(item)
                matrix = np.load(bound(item['prediction']), mmap_mode='r')
                if matrix.shape != (len(cases), len(genes)):
                    raise ValueError('prediction shape differs from cases and genes')
                if 'prediction_order' in item:
                    order = read(bound(item['prediction_order']))
                    expected = cases[['case_id', 'group_id']].values.tolist() if order and isinstance(order[0], list) else cases.case_id.tolist()
                    if order != expected:
                        raise ValueError('prediction row order differs from cases')
            continuity = ContinuousMetrics(item, bound)
            continuous_values = {m: [] for m in CONTINUOUS}
            continuous_present = set()
            values = {m: [] for m in METRICS}
            dose_values = defaultdict(lambda: {m: [] for m in METRICS})
            dose_counts = defaultdict(int)
            terminal = {m: [] for m in METRICS}
            groups = cases.groupby('group_id', sort=False).indices
            for group, rows in groups.items():
                truth = labels.truth(group, cases.iloc[rows].case_id.tolist())
                predicted = vector if matrix is None else mean_rows(matrix, rows)
                curves = retrieval(predicted[covered], truth, columns[covered]) if truth is not None else dict.fromkeys(METRICS)
                row = {key: item[key] for key in KEYS}
                row.update(group_id=group, dose_um=cases.iloc[rows[0]].dose_um, wells=len(rows),
                           significance_status='eligible' if truth is not None else 'untestable')
                dose_counts[row['dose_um']] += 1
                dose_values[row['dose_um']]
                for metric, curve in curves.items():
                    row[metric] = curve.ap if curve is not None else None
                    row[f'{metric}_status'] = 'defined' if curve is not None else ('zero_significant' if truth is not None else 'untestable')
                    row[f'{metric}_terminal_recall'] = curve.recall[-1] if curve is not None else None
                    if curve is not None:
                        values[metric].append(curve.ap)
                        dose_values[row['dose_um']][metric].append(curve.ap)
                        terminal[metric].append(curve.recall[-1])
                    key = (*[item[k] for k in ('dataset', 'regime', 'role', 'model')], metric)
                    pooled[key].add(curve)
                continuous_row = continuity.row(group, predicted)
                row.update(continuous_row)
                continuous_present.update(continuous_row)
                for metric, value in continuous_row.items():
                    dose_values[row['dose_um']].setdefault(metric, [])
                    if value is not None:
                        continuous_values[metric].append(value)
                        dose_values[row['dose_um']][metric].append(value)
                writer.writerow(row)
            report = {key: item[key] for key in KEYS}
            report.update(groups=len(groups), covered_genes=int(covered.sum()))
            report['deterministic_reference'] = item['model'] in ('no_change', 'delta_mean')
            report.update(item.get('continuous_metrics', {}))
            for metric in (name for name in CONTINUOUS if name in continuous_present):
                numbers = continuous_values[metric]
                value = float(np.mean(numbers)) if numbers else None
                prior = item.get('continuous_metrics', {}).get(metric)
                if prior is not None and (value is None or not np.isclose(prior, value, rtol=1e-7, atol=1e-9)):
                    raise ValueError(f"{item['model']}/{item['regime']}: preserved {metric} differs from group scores")
                report[metric] = prior if prior is not None else value
                report[f'{metric}_defined_groups'] = len(numbers)
                report[f'{metric}_undefined_groups'] = len(groups) - len(numbers)
            report.update({key: item[key] for key in ('continuous_gene_count', 'continuous_provisional') if key in item})
            for metric in METRICS:
                report[metric] = float(np.mean(values[metric])) if values[metric] else None
                report[f'{metric}_defined_groups'] = len(values[metric])
                report[f'{metric}_undefined_groups'] = len(groups) - len(values[metric])
                report[f'{metric}_terminal_recall_mean'] = float(np.mean(terminal[metric])) if terminal[metric] else None
            reports.append(report)
            for dose, metrics in dose_values.items():
                dose_row = {key: item[key] for key in KEYS}
                dose_row['dose_um'] = dose
                dose_row['groups'] = dose_counts[dose]
                for metric, numbers in metrics.items():
                    dose_row[metric] = float(np.mean(numbers)) if numbers else None
                    dose_row[f'{metric}_defined_groups'] = len(numbers)
                    dose_row[f'{metric}_undefined_groups'] = dose_counts[dose] - len(numbers)
                dose_reports.append(dose_row)
            print(f"run {number + 1}/{len(entries)}: {item['dataset']}/{item['model']}/{item['fold']}/{item['seed']}", flush=True)
    pd.DataFrame(reports).to_csv(destination / 'per_run.csv', index=False)
    pd.DataFrame(dose_reports).to_csv(destination / 'per_dose.csv', index=False)
    frames = []
    for key, pool in pooled.items():
        frame = pool.frame()
        for column, value in zip(('dataset', 'regime', 'role', 'model', 'metric'), key):
            frame[column] = value
        frames.append(frame)
    pd.concat(frames, ignore_index=True).to_csv(destination / 'pr.csv', index=False)
    data = pd.DataFrame(reports)
    group = ['dataset', 'regime', 'role', 'model']
    references = data[data.deterministic_reference].drop_duplicates([*group, 'fold'])
    summary_data = pd.concat([data[~data.deterministic_reference], references], ignore_index=True)
    summary_metrics = [*METRICS, *[m for m in CONTINUOUS if m in data]]
    summary = summary_data.groupby(group, sort=False)[summary_metrics].agg(['mean', 'std', 'count'])
    summary.columns = ['_'.join(c) for c in summary.columns]
    summary.reset_index().to_csv(destination / 'summary.csv', index=False)
    from stage2.evaluation.figures import render
    figure_destination = figure_destination or destination
    render(destination, figure_destination)
    # Detect in-place input changes before publishing any replay manifest.
    for item in entries:
        for value in iter_bindings(item):
            bound(value)
    for labels in label_sets.values():
        labels.verify()
    if any(sha(Path(__file__).with_name(name)) != digest for name, digest in source_hashes.items()):
        raise ValueError('scoring source changed during evaluation')
    manifest = {'schema': 'stage2-significance-scoring-v1', 'metrics': METRICS, 'sources': source_hashes,
                'labels': {k: v.manifest_digest for k, v in label_sets.items()},
                'inventory_runs': entries,
                'files': {p.name: sha(p) for p in destination.glob('*.csv')},
                'figures': {name: sha(figure_destination / name)
                            for name in ('PR_curve.png', 'PR_curve.svg', 'AUPRC.png', 'AUPRC.svg')}}
    (destination / 'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('score', 'replay'))
    parser.add_argument('--step', type=int, choices=STEPS, required=True)
    parser.add_argument('--dataset', choices=('tenv',), default='tenv')
    parser.add_argument('--labels-root', type=Path, help='explicit evaluator label release root')
    parser.add_argument('--inventory', type=Path, default=Path(__file__).with_name('inventory.json'))
    parser.add_argument('--destination', type=Path)
    args = parser.parse_args()
    entries = [r for r in load_inventory(args.inventory)['runs'] if r['step'] == args.step and (not args.dataset or r['dataset'] == args.dataset)]
    if not entries:
        raise ValueError('empty scoring inventory')
    destination = args.destination or STEPS[args.step] / 'output/retrieval'
    figure_destination = STEPS[6] / 'output' if args.step == 6 and args.destination is None else destination
    if args.action == 'score':
        score(entries, destination, args.labels_root, figure_destination)
    else:
        from stage2.evaluation.publish import verify_existing
        verify_existing(destination, figure_destination)
        with tempfile.TemporaryDirectory(dir=destination.parent) as temp:
            replay = Path(temp) / 'replay'
            score(entries, replay, args.labels_root)
            expected = read(destination / 'manifest.json')
            actual = read(replay / 'manifest.json')
            if expected != actual or sha(destination / 'manifest.json') != sha(replay / 'manifest.json'):
                raise ValueError('retrieval replay differs')
            print('All scoring tables and figures replay byte-for-byte', flush=True)


if __name__ == '__main__':
    main()
