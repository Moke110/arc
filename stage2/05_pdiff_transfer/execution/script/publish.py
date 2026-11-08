#!/usr/bin/env python3
"""Score the frozen Step 05 matrix and replay its publication offline."""
import argparse
from collections import Counter, defaultdict
from contextlib import contextmanager
import csv
import math
from pathlib import Path
import shutil
import tempfile

import formal_matrix as matrix
import tahoe_zero_shot as tahoe
from artifacts import SEEDS, read_json, rows, sha256, write_csv, write_json
from evaluate import (GROUP_FIELDS, METRICS, RECALL_GRID, REFERENCES, RUN_FIELDS, evaluator_inputs,
                      grouped_cases, load_array, mean_delta, np, run_scores, summaries, verify_file)
from tahoe_bridge import LIMITATION
from training_inputs import json_hash
from training_policy import verify_best_checkpoint
from transfer import RECORD, REPO

SHARED = 'cold-shared'
SCORE_SCHEMA = 'step05-scores-v1'
PUBLICATION_SCHEMA = 'step05-publication-v1'
MODELS = (*matrix.ARMS, *REFERENCES)
COPY_FILES = ('per_run.csv', 'learning_points.json', 'pr_points.json', 'undefined_reasons.json')
COMPLETE = ('completed', 'reused')


@contextmanager
def fresh_directory(output):
    if output.exists() or output.is_symlink():
        raise ValueError(f'publication requires a fresh output directory: {output}')
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        pending = Path(temporary) / 'result'
        pending.mkdir()
        yield pending
        pending.rename(output)


def hashes(root):
    return {str(path.relative_to(root)): sha256(path)
            for path in sorted(root.rglob('*')) if path.is_file()}


def sealed_record(root, name, digest, schema, *, paths=None):
    files = None if paths is None else set(paths) - {name}
    if paths is not None and (name not in paths or len(set(paths.values())) != len(paths)):
        raise ValueError('invalid tracked publication paths')

    def resolve(relative):
        mapped = relative if paths is None else paths[relative]
        if Path(mapped).is_absolute():
            raise ValueError(f'artifact path must be relative: {mapped}')
        path = (root / mapped).resolve()
        if not path.is_relative_to(root.resolve()):
            raise ValueError(f'artifact is outside the record: {mapped}')
        return path

    record = read_json(verify_file(resolve(name), digest))
    if record['schema'] != schema:
        raise ValueError(f'unexpected record schema: {name}')
    if files is not None and files != set(record['files']):
        raise ValueError('tracked publication paths differ from the sealed file inventory')
    for relative, checksum in record['files'].items():
        verify_file(resolve(relative), checksum)
    return record


def checkpoint_index(plan, inventory):
    checkpoints, slots = {}, {}
    planned = {item['id']: item for item in plan['slots']}

    def register(checkpoint, history, request, recipe, stage):
        run = read_json(history['path'])
        best = verify_best_checkpoint(run, history['path'])
        digest = checkpoint['sha256']
        if (digest != run['checkpoint_sha256']
                or recipe.get('resolved_sha256', run['recipe_sha256']) != run['recipe_sha256']
                or not all(math.isfinite(item['training_loss']) for item in run['history'])):
            raise ValueError(f'checkpoint or recipe differs from its history: {history["path"]}')
        item = {
            'stage': stage, 'identity': run['identity'], 'checkpoint': checkpoint, 'history': history,
            'request': request, 'recipe': recipe, 'recipe_sha256': run['recipe_sha256'],
            'model_state_sha256': run['model_state_sha256'], 'best_step': best['step'],
            'best_validation_loss': best['validation_loss'], 'optimizer_steps': run['optimizer_steps'],
            'points': [{key: event[key] for key in ('step', 'training_loss', 'validation_loss')}
                       for event in run['history']],
        }
        if checkpoints.setdefault(digest, item) != item:
            raise ValueError(f'different training records for checkpoint {digest}')
        return digest

    for row in inventory['slots']:
        key, regime, arm, seed = (row[field] for field in ('id', 'regime', 'arm', 'seed'))
        slot = planned[key]
        link = {'regime': regime, 'arm': arm, 'seed': seed, 'fold': row['fold'], 'status': row['status']}
        if row['status'] == 'reused':
            reuse = row['reuse']
            files = {Path(item['path']).name: item for item in matrix.cm_bindings(reuse)}
            link['target'] = register(files['best.pt'], files['run.json'], files['request.json'],
                                      files['recipe.json'], 'tenv target')
        if arm in ('Tahoe', 'Transfer'):
            source = slot['source']
            if regime == 'cold-drug':
                source_fit = next(item['fit'] for item in inventory['slots']
                                  if item['id'] == f'{regime}-tahoe-{seed}')
                if source_fit['status'] == 'completed':
                    fit = source_fit['artifacts']
                    link['source_fit_seal'] = fit['fit']
                    link['source'] = register(fit['checkpoint'], fit['history'], source['request'],
                                              source['recipe'], 'Tahoe source')
            else:
                link['source'] = register({'path': source['path'], 'sha256': source['checkpoint_sha256']},
                                          source['history'], source['request'], source['recipe'],
                                          'Tahoe source')
        if arm != 'Tahoe' and row.get('fit', {}).get('status') == 'completed':
            fit = row['fit']['artifacts']
            if arm == 'Mix':
                owner = read_json(plan['worker_plans']['mix']['path'])
                request = {'artifact': plan['worker_plans']['mix'], 'embedded_request': regime,
                           'request_sha256': plan['targets'][regime]['request_sha256_by_arm']['Mix']}
                if json_hash(owner['requests'][regime]) != request['request_sha256']:
                    raise ValueError(f'{key} differs from the embedded pooled request')
            else:
                request = plan['targets'][regime]['request_files']['CM' if arm == 'CM' else 'Tahoe']
            link['target'] = register(fit['checkpoint'], fit['history'], request,
                                      slot['target_recipe'], 'tenv target')
            link['fit_seal'] = fit['fit']
            if arm == 'Transfer':
                lineage = fit['initialization']
                if lineage['source_checkpoint_sha256'] != link.get('source'):
                    raise ValueError(f'{key} source differs from the selected Tahoe checkpoint')
                link['initialization'] = lineage
        slots[key] = link
    return {'checkpoints': checkpoints, 'slots': slots}


def parent_overlap(plan, regime, target_request):
    source = next(row['source'] for row in plan['slots']
                  if row['regime'] == regime and row['arm'] == 'Tahoe')
    source_request = read_json(source['request']['path'])
    parents = tahoe.parent_blocks(source_request, tahoe.role_cases(source_request))
    source_roles = {role: {'parent_connectivity': parents[role]} for role in parents}
    overlap = tahoe.overlap({'roles': source_roles}, target_request)
    if regime == 'cold-drug' and any(overlap[role][target] for role in ('train', 'validation')
                                     for target in ('validation', 'test')):
        raise ValueError('cold-drug has a held-out tenv parent in the Tahoe fit membership')
    pool = plan['targets'][regime]['mix_source']
    pool_request = {'predictor_files': {key: item['path'] for key, item in pool.items()}}
    pool_parents = tahoe.parent_blocks(pool_request, {'train': rows(pool['train_cases']['path'])})['train']
    if pool_parents != parents['train']:
        raise ValueError(f'{regime} Mix and Tahoe source train parent memberships differ')
    return {
        'Tahoe_and_Transfer_source_roles': overlap,
        'Mix_Tahoe_train': overlap['train'],
        'counts': {role: {target: len(blocks) for target, blocks in destinations.items()}
                   for role, destinations in overlap.items()},
    }


def prediction_path(row, case_order, n_genes):
    if row['status'] not in COMPLETE:
        return None, None
    if row['status'] == 'reused':
        reuse = row['reuse']
        directory = Path(reuse['prediction_manifest']).parent
        cases, path = directory / 'cases.json', directory / 'prediction.npy'
        evidence = {'path': reuse['prediction_manifest'], 'sha256': reuse['prediction_manifest_sha256']}
    else:
        artifacts = row['prediction']['artifacts']
        cases, path = (Path(artifacts[field]['path']) for field in ('cases', 'values'))
        evidence = artifacts['prediction']
    if read_json(cases) != case_order:
        raise ValueError(f'{row["id"]} differs from the frozen test case/group order')
    return load_array(path, (len(case_order), n_genes)), evidence


def failure_reason(row):
    outcome = row['prediction']
    if outcome.get('root_failures'):
        return 'blocked by ' + ', '.join(outcome['root_failures'])
    attempts = ', '.join(
        str(item.get('reason') or item.get('error') or item.get('status') or item.get('returncode', 'unknown'))
        for item in outcome['attempts']
    )
    reason = outcome.get('reason') or outcome.get('error') or attempts or outcome['status']
    return f'{outcome["status"]}: {reason}' if reason != outcome['status'] else reason


def verify_cold_reuse(directory, runs, curves):
    """Compare the replayed shared scores with the independently qualified Step 04 replay."""
    record = RECORD / 'cold_record/replay.json'
    pin = read_json(record)['local_artifact_sha256']
    historical_runs = rows(verify_file(directory / 'per_run.csv', pin[f'{directory.name}/per_run.csv']))
    historical_curves = read_json(verify_file(directory / 'tenv/pr_curves.json',
                                            pin[f'{directory.name}/tenv/pr_curves.json']))
    shared = {(item['model'], item['seed']): item for item in runs if item['scope'] == SHARED}
    for old in historical_runs:
        arm = 'CM' if old['model'] == 'perturbdiff' else old['model']
        item = shared[arm, int(old['seed']) if old['seed'] else None]
        for key in RUN_FIELDS:
            if key == 'model':
                continue
            value = item[key]
            if ('' if value is None else str(value)) != old[key]:
                raise ValueError(f'cold replay diverges from Step 04 for {arm} {item["seed"]}: {key}')
    pr = {(item['model'], item['seed']): item['precision'] for item in curves if item['scope'] == SHARED}
    if historical_curves['recall'] != RECALL_GRID.tolist() or any(
            pr['CM' if old['model'] == 'perturbdiff' else old['model'], old['seed']] != old['precision']
            for old in historical_curves['runs']):
        raise ValueError('cold directional PR replay diverges from Step 04')


def undefined_reasons(path):
    reasons = defaultdict(Counter)
    with path.open(newline='') as stream:
        for row in csv.DictReader(stream):
            key = (row['scope'], row['model'], row['seed'])
            for metric in METRICS:
                reason = row[f'{metric}_reason']
                if reason:
                    reasons[key][metric, reason] += 1
    return [{'scope': scope, 'model': model, 'seed': int(seed) if seed else None,
             'reasons': [{'metric': metric, 'reason': reason, 'groups': count}
                         for (metric, reason), count in counts.items()]}
            for (scope, model, seed), counts in reasons.items()]


def verify_tahoe_membership(entries, request, overlap):
    """Check the recipient and exact overlap fields against verified membership."""
    for row in entries:
        if row['arm'] != 'Tahoe' or row['status'] != 'completed':
            continue
        evidence = row['prediction']['artifacts']['prediction']
        manifest = read_json(evidence['path'])
        if (manifest['recipient_identity'] != request['identity']
                or manifest['source_parent_overlap'] != overlap['Tahoe_and_Transfer_source_roles']):
            raise ValueError(f'{row["id"]} Tahoe prediction differs from the verified recipient or source overlap')


def score(runtime, formal, output, plan_digest, inventory_digest):
    verify_file(formal / 'prediction_inventory.json', inventory_digest)
    plan = matrix.read_plan(formal, plan_digest)
    if plan['runtime'] != str(runtime):
        raise ValueError('runtime differs from the frozen matrix checkout')
    matrix.verify_plan(plan)
    inventory = matrix.verify_inventory(formal, plan_digest, plan)
    learning = checkpoint_index(plan, inventory)
    runs, curves, overlap, coverage, sources = [], [], {}, {}, {}
    entries = {regime: [row for row in inventory['slots'] if row['regime'] == regime]
               for regime in matrix.REGIMES}
    recipients = {}
    for regime in matrix.REGIMES:
        recipient = plan['targets'][regime]['recipient_plan']
        directory = Path(recipient['path']).parent
        request, _, fold, warm_source = tahoe.recipient(runtime, directory, regime, recipient['sha256'])
        overlap[regime] = parent_overlap(plan, regime, request)
        verify_tahoe_membership(entries[regime], request, overlap[regime])
        recipients[regime] = request, fold, warm_source, directory
    case_orders = {
        regime: [[row['case_id'], row['group_id']] for row in rows(target['cases']['test']['path'])]
        for regime, target in plan['targets'].items()
    }
    predictions = {
        row['id']: prediction_path(row, case_orders[row['regime']], plan['targets'][row['regime']]['n_genes'])
        for row in inventory['slots']
    }
    with fresh_directory(output) as pending, tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        scratch = Path(temporary)
        with (pending / 'per_group.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=('scope', *GROUP_FIELDS))
            writer.writeheader()
            for regime in matrix.REGIMES:
                request, fold, warm_source, directory = recipients[regime]
                cases = rows(request['predictor_files']['test_cases'])
                order = [[row['case_id'], row['group_id']] for row in cases]
                target = plan['targets'][regime]
                if (order != case_orders[regime]
                        or json_hash(order) != target['cases']['test']['ordered_case_group_sha256']):
                    raise ValueError(f'{regime} differs from the frozen tenv test membership')
                indices, _ = grouped_cases(cases)
                if regime == 'warm-start':
                    with tahoe.score_truth(runtime, request, fold, warm_source, scratch) as inputs:
                        warm_inputs = inputs
                    training = rows(target['predictor_files']['train_cases']['path'])
                    train_groups, _ = grouped_cases(training)
                    delta = load_array(target['predictor_files']['train_expression']['path'],
                                       (len(training), target['n_genes']))
                    mean = np.zeros(target['n_genes'])
                    for members in train_groups.values():
                        mean += mean_delta(delta, members)
                    references = {'no_change': np.zeros(target['n_genes']), 'delta_mean': mean / len(train_groups)}
                    groups, truth, labels = warm_inputs
                elif regime == 'cold-cm-drug':
                    groups, truth, labels, references = evaluator_inputs(fold)
                # Both cold regimes use the same evaluator inputs. Score shared rows in the first.
                coverage[regime] = {
                    'expected_learned_slots': len(entries[regime]),
                    'completed_learned_slots': sum(row['status'] in COMPLETE for row in entries[regime]),
                    'reused_cm_slots': sum(row['status'] == 'reused' for row in entries[regime]),
                }
                for row in entries[regime]:
                    values, evidence = predictions.pop(row['id'])
                    sources[row['id']] = evidence
                    if row['status'] == 'reused' and regime == 'cold-drug':
                        continue
                    slot = {'scope': SHARED if row['status'] == 'reused' else regime, 'dataset': 'tenv',
                            'fold': row['fold'], 'seed': row['seed'], 'model': row['arm']}
                    status = 'completed' if evidence else row['status']
                    reason = '' if evidence else failure_reason(row)
                    runs.append(run_scores(writer, slot, status, reason, groups, truth, labels,
                                           indices, values, curves))
                if regime != 'cold-drug':
                    for arm, delta in references.items():
                        slot = {'scope': SHARED if fold is not None else regime, 'dataset': 'tenv',
                                'fold': 0, 'seed': None, 'model': arm}
                        runs.append(run_scores(writer, slot, 'completed', '', groups, truth, labels,
                                               indices, delta, curves))
                if regime == 'cold-cm-drug':
                    verify_cold_reuse(directory, runs, curves)
        write_csv(pending / 'per_run.csv', runs, ('scope', *RUN_FIELDS))
        write_json(pending / 'pr_points.json', {'recall': RECALL_GRID.tolist(), 'runs': curves,
                                                'prediction_manifests': sources})
        write_json(pending / 'learning_points.json', learning)
        write_json(pending / 'undefined_reasons.json', undefined_reasons(pending / 'per_group.csv'))
        write_json(pending / 'scores.json', {
            'schema': SCORE_SCHEMA, 'dose_limitation': LIMITATION, 'coverage': coverage, 'overlap': overlap,
            'targets': {regime: {'split': target['split'], 'cases': target['cases'],
                                 'identity': target['identity'], 'n_genes': target['n_genes'],
                                 'step04_test_reuse': regime != 'warm-start',
                                 'shared_scores': SHARED if regime != 'warm-start' else None}
                        for regime, target in plan['targets'].items()},
            'inputs': {'producer_revision': plan['producer_revision'],
                       'plan_sha256': plan_digest, 'inventory_sha256': inventory_digest,
                       'cold_cm_and_references_replay': 'exact'},
            'files': hashes(pending),
        })


def read_runs(path):
    result = []
    for row in rows(path):
        row['fold'] = int(row['fold'])
        row['seed'] = int(row['seed']) if row['seed'] else None
        row['group_count'] = int(row['group_count'])
        for metric in METRICS:
            row[metric] = float(row[metric]) if row[metric] else None
            count = f'{metric}_valid_groups'
            row[count] = int(row[count])
            if (not 0 <= row[count] <= row['group_count']
                    or (row[metric] is None) != (row[count] == 0)
                    or (row[metric] is not None and not math.isfinite(row[metric]))):
                raise ValueError(f'invalid metric or valid group count: {metric}')
        result.append(row)
    return result


def scoped_runs(runs, regime):
    scopes = (regime,) if regime == 'warm-start' else (regime, SHARED)
    selected = [row for row in runs if row['scope'] in scopes]
    expected = {(arm, seed) for arm in matrix.ARMS for seed in SEEDS}
    expected.update((arm, None) for arm in REFERENCES)
    if ({(row['model'], row['seed']) for row in selected} != expected
            or len(selected) != len(expected)
            or any((row['dataset'], row['fold']) != ('tenv', 0) for row in selected)):
        raise ValueError(f'invalid run coverage for {regime}')
    return selected


def paired_ap100(runs, regime):
    by_arm = {(row['model'], row['seed']): row for row in scoped_runs(runs, regime)}
    pairs = []
    for seed in SEEDS:
        reasons, means, counts = [], {}, {}
        for arm in ('CM', 'Transfer'):
            run = by_arm[arm, seed]
            directional = {metric: run[metric] for metric in METRICS[:2]}
            counts[arm] = {metric: run[f'{metric}_valid_groups'] for metric in directional}
            means[arm] = None
            if run['status'] != 'completed':
                reasons.append(f'{arm} {run["status"]} ({run["reason"]})')
            else:
                undefined = [metric for metric, value in directional.items() if value is None]
                if undefined:
                    reasons.append(f'{arm} undefined {", ".join(undefined)}')
                else:
                    means[arm] = sum(directional.values()) / 2
        pairs.append({'fold': 0, 'seed': seed, 'CM': means['CM'], 'Transfer': means['Transfer'],
                      'difference': None if reasons else means['Transfer'] - means['CM'],
                      'valid_group_counts': counts, 'reason': '. '.join(reasons)})
    differences = [item['difference'] for item in pairs if item['difference'] is not None]
    return {'expected_pairs': len(SEEDS), 'complete_pairs': len(differences), 'pairs': pairs,
            'mean_difference': sum(differences) / len(SEEDS) if len(differences) == len(SEEDS) else None}


COLORS = dict(zip(MODELS, ('#2166ac', '#d87a27', '#24804b', '#8255a2', '#8a8e95', '#4d5258')))
SEED_COLORS = dict(zip(SEEDS, ('#2166ac', '#d87a27', '#8255a2')))
METRIC_TITLES = ('Positive AP@100', 'Negative AP@100', 'Positive recall at P80',
                 'Negative recall at P80', 'PCC', 'Spearman', 'RMSE', 'R²')


def figure_axis(axis):
    axis.set_axisbelow(True)
    axis.grid(axis='y', color='#e4e8ed', linewidth=0.8)
    axis.spines[['top', 'right']].set_visible(False)
    axis.tick_params(labelsize=8, length=2)


def save_figure(plt, figure, path):
    figure.tight_layout(rect=(0, 0.045, 1, 0.96))
    figure.savefig(path, dpi=160)
    plt.close(figure)


def draw_learning(plt, data, path):
    figure, axes = plt.subplots(3, 4, figsize=(12.8, 8.7))
    figure.suptitle('Validation selection by regime and training stage', fontsize=15)
    titles = ('Tahoe source', 'CM target', 'Transfer target', 'Mix · pooled fit')
    for index, regime in enumerate(matrix.REGIMES):
        for column, (title, arm) in enumerate(zip(titles, ('Tahoe', 'CM', 'Transfer', 'Mix'))):
            axis = axes[index, column]
            figure_axis(axis)
            if not index:
                axis.set_title(title)
            if not column:
                axis.set_ylabel(f'{regime}\nNative loss')
            if index == len(matrix.REGIMES) - 1:
                axis.set_xlabel('Optimizer step')
            seen = False
            for seed, color in SEED_COLORS.items():
                key = data['slots'].get(f'{regime}-{arm.lower()}-{seed}', {})
                selected = key.get('source' if arm == 'Tahoe' else 'target')
                if selected is None:
                    continue
                run = data['checkpoints'][selected]
                steps = [event['step'] for event in run['points']]
                for kind, style in (('training_loss', ':'), ('validation_loss', '-')):
                    axis.plot(steps, [event[kind] for event in run['points']], style, color=color, lw=1.25)
                axis.scatter(run['best_step'], run['best_validation_loss'], marker='o', s=21, color=color,
                             edgecolor='white', linewidth=0.35, zorder=4)
                seen = True
            if not seen:
                axis.text(0.5, 0.5, 'No selected checkpoint', ha='center', transform=axis.transAxes, fontsize=9)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=color, label=f'Seed {seed}') for seed, color in SEED_COLORS.items()]
    handles += [Line2D([], [], color='#30343a', ls=style, marker=marker, label=label)
                for style, marker, label in ((':', '', 'Training'), ('-', '', 'Validation'),
                                             ('', 'o', 'Selected validation minimum'))]
    figure.legend(handles=handles, loc='lower center', ncol=6, frameon=False)
    save_figure(plt, figure, path)


def draw_metrics(plt, runs, path):
    figure, axes = plt.subplots(3, 8, figsize=(17, 7.6))
    figure.suptitle('Held-out tenv group metrics', fontsize=15)
    for index, regime in enumerate(matrix.REGIMES):
        selected = scoped_runs(runs, regime)
        for column, (metric, title) in enumerate(zip(METRICS, METRIC_TITLES)):
            axis = axes[index, column]
            figure_axis(axis)
            if not index:
                axis.set_title(title, fontsize=10)
            if not column:
                axis.set_ylabel(regime)
            axis.set_xlim(-0.5, len(MODELS) - 0.5)
            axis.set_xticks(range(len(MODELS)), [])
            for position, arm in enumerate(MODELS):
                values = [row[metric] for row in selected if row['model'] == arm and row[metric] is not None]
                if values:
                    axis.scatter(position + np.linspace(-0.12, 0.12, len(values)) if len(values) > 1 else [position],
                                 values, color=COLORS[arm], s=16, alpha=0.58)
                    axis.hlines(np.mean(values), position - 0.27, position + 0.27,
                                color=COLORS[arm], lw=1.7)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker='o', color=color, label=arm) for arm, color in COLORS.items()]
    figure.legend(handles=handles, loc='lower center', ncol=len(MODELS), frameon=False)
    save_figure(plt, figure, path)


def draw_pr(plt, record, path):
    figure, axes = plt.subplots(3, 2, figsize=(10, 10.1), sharex=True, sharey=True)
    figure.suptitle('Directional precision–recall on held-out tenv groups', fontsize=15)
    for index, regime in enumerate(matrix.REGIMES):
        selected = scoped_runs(record['runs'], regime)
        for direction, title in enumerate(('Positive changes', 'Negative changes')):
            axis = axes[index, direction]
            figure_axis(axis)
            axis.set_title(title if not index else '')
            axis.set(xlim=(0, 1), ylim=(0, 1.015))
            if direction == 0:
                axis.set_ylabel(f'{regime}\nPrecision')
            if index == len(matrix.REGIMES) - 1:
                axis.set_xlabel('Recall')
            for arm in MODELS:
                curves = [row['precision'][direction] for row in selected
                          if row['model'] == arm and row['precision'][direction] is not None]
                if len(curves) > 1:
                    for values in curves:
                        axis.plot(record['recall'], values, color=COLORS[arm], lw=0.85, alpha=0.26)
                if curves:
                    axis.plot(record['recall'], np.mean(curves, axis=0), color=COLORS[arm], lw=1.5,
                              ls='--' if arm in REFERENCES else '-')
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=color, label=arm, ls='--' if arm in REFERENCES else '-')
               for arm, color in COLORS.items()]
    figure.legend(handles=handles, loc='lower center', ncol=len(MODELS), frameon=False)
    save_figure(plt, figure, path)


def report(record, paired, runs):
    lines = [
        '# Step 05 · Tahoe transfer to tenv', '',
        '> **Dose encoding limitation.** ' + record['dose_limitation'] +
        ' Transfer starts from this mapping and can update the dose embeddings on tenv. '
        'Improving the cross-dataset encoder is a priority for follow-up.', '',
        'All arms predict held-out tenv drug–dose group means. The cold regimes reuse the Step 04 '
        'tenv fold 0 test set, which was also scored in Step 03. Warm-start uses its own held-out '
        'parent–dose combinations. Each regime is a fixed fold, so seed variation describes repeated '
        'fits on that fold.', '',
        '## Transfer minus CM', '',
        'The paired difference uses the equal mean of positive and negative AP@100. A regime mean '
        'requires all three seed pairs. Group-level scores and undefined reasons are in the full '
        'scoring output. The compact release records [undefined reason counts](undefined_reasons.json) '
        'and each metric’s valid-group counts in [per_run.csv](per_run.csv). The cold CM and '
        'deterministic references share the '
        '`cold-shared` scope.', '',
        '| Regime | Learned slots | Seed 17 | Seed 42 | Seed 2026 | Complete pairs | Mean difference |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: |',
    ]

    def delta(value):
        return f'{value:+.5f}' if value is not None else '—'

    notes = []
    for regime, comparison in paired.items():
        selected = scoped_runs(runs, regime)
        learned = [row for row in selected if row['model'] in matrix.ARMS]
        coverage = f'{sum(row["status"] == "completed" for row in learned)}/{len(learned)}'
        cells = [regime, coverage, *(delta(row['difference']) for row in comparison['pairs']),
                 f'{comparison["complete_pairs"]}/{comparison["expected_pairs"]}',
                 delta(comparison['mean_difference'])]
        lines.append('| ' + ' | '.join(cells) + ' |')
        notes.extend(f'- {regime}, seed {pair["seed"]}. {pair["reason"]}'
                     for pair in comparison['pairs'] if pair['reason'])
    if notes:
        lines += ['', *notes]
    if record['overlap']:
        lines += ['', '## Tahoe parent overlap', '',
                  'Counts reflect the realized Tahoe source memberships and the tenv roles. '
                  'Mix pools the Tahoe training membership. Tahoe and Transfer use the source training '
                  'and validation memberships.', '',
                  '| Regime | Tahoe train ∩ tenv validation | Tahoe train ∩ tenv test | '
                  'Tahoe validation ∩ tenv validation | Tahoe validation ∩ tenv test |',
                  '| --- | ---: | ---: | ---: | ---: |']
        for regime, overlap in record['overlap'].items():
            counts = overlap['counts']
            cells = [regime, *(str(counts[source][target]) for source in ('train', 'validation')
                               for target in ('validation', 'test'))]
            lines.append('| ' + ' | '.join(cells) + ' |')
    lines += ['', '## Results and checkpoint selection', '',
              '- [Eight-metric summary](summary.csv) and [held-out results](metrics.png)',
              '- [Directional precision–recall plots](pr_curves.png) and [per-run sources](pr_points.json)',
              '- [Learning curves](learning_curves.png) and [checkpoint index](checkpoints.json)',
              '- [Publication checksums](publication.json)', '']
    return '\n'.join(lines)


def render(scores, output, digest):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    record = sealed_record(scores, 'scores.json', digest, SCORE_SCHEMA)
    if not set(COPY_FILES).issubset(record['files']):
        raise ValueError('score record does not bind all publication inputs')
    runs = read_runs(scores / 'per_run.csv')
    learning = read_json(scores / 'learning_points.json')
    pr = read_json(scores / 'pr_points.json')
    paired = {regime: paired_ap100(runs, regime) for regime in matrix.REGIMES}
    with fresh_directory(output) as pending:
        for name in COPY_FILES:
            shutil.copyfile(scores / name, pending / name)
        write_json(pending / 'paired_ap100.json', paired)
        write_json(pending / 'checkpoints.json', {
            'schema': 'step05-checkpoints-v1', 'inputs': record['inputs'], 'slots': learning['slots'],
            'checkpoints': {digest: {key: value for key, value in item.items() if key != 'points'}
                            for digest, item in learning['checkpoints'].items()},
        })
        write_csv(pending / 'summary.csv', [
            {'scope': scope, **row} for scope in (*matrix.REGIMES, SHARED)
            for row in summaries(run for run in runs if run['scope'] == scope)
        ])
        (pending / 'README.md').write_text(report(record, paired, runs), encoding='utf-8')
        draw_learning(plt, learning, pending / 'learning_curves.png')
        draw_metrics(plt, runs, pending / 'metrics.png')
        draw_pr(plt, pr, pending / 'pr_curves.png')
        write_json(pending / 'publication.json', {
            'schema': PUBLICATION_SCHEMA, 'scores_sha256': digest, 'inputs': record['inputs'],
            'targets': record['targets'], 'coverage': record['coverage'], 'overlap': record['overlap'],
            'score_files': record['files'], 'files': hashes(pending),
        })


def verify_replay(scores, output, digest, *, paths=None):
    record = sealed_record(output, 'publication.json', digest, PUBLICATION_SCHEMA, paths=paths)
    if scores is None:
        return
    with tempfile.TemporaryDirectory() as temporary:
        replay = Path(temporary) / 'publication'
        render(scores, replay, record['scores_sha256'])
        verify_file(replay / 'publication.json', digest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, default=REPO)
    commands = parser.add_subparsers(dest='command', required=True)
    scorer = commands.add_parser('score')
    scorer.add_argument('formal', type=Path)
    scorer.add_argument('output', type=Path)
    scorer.add_argument('--plan-sha256', required=True)
    scorer.add_argument('--inventory-sha256', required=True)
    renderer = commands.add_parser('render')
    renderer.add_argument('scores', type=Path)
    renderer.add_argument('output', type=Path)
    renderer.add_argument('--scores-sha256', required=True)
    verifier = commands.add_parser('verify-replay')
    verifier.add_argument('scores', type=Path)
    verifier.add_argument('output', type=Path)
    verifier.add_argument('--publication-sha256', required=True)
    release = commands.add_parser('verify-release', help='verify the tracked Step 05 publication')
    release.add_argument('step', type=Path)
    release.add_argument('--scores', type=Path, help='also reproduce the seal from this scoring output')
    args = parser.parse_args()
    try:
        if args.command == 'score':
            score(args.runtime.resolve(), args.formal.resolve(), args.output,
                  args.plan_sha256, args.inventory_sha256)
        elif args.command == 'render':
            render(args.scores, args.output, args.scores_sha256)
        elif args.command == 'verify-replay':
            verify_replay(args.scores, args.output, args.publication_sha256)
        else:
            publication = read_json(args.step / 'execution/formal_record/run.json')['publication']
            verify_replay(args.scores, args.step, publication['manifest']['sha256'],
                          paths=publication['tracked_artifact_paths'])
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
