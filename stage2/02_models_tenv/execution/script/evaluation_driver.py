"""Score a fold's persisted tenv predictions against evaluator-only truth (#164).

The driver is the only evaluator in ``stage2/02_models_tenv``. It reads the
step's published run artifacts, checks that the run set is complete and that
every fold case is covered, then scores the persisted predictions through the
shared implementation in ``stage2/01_bases/output/api/bases/evaluation.py``
(ADR 0063). It writes the per-well, per-group and per-run tables plus the
long-form plot-value tables both report figures render from. No metric is
re-implemented here.

Fail-closed: a run set that is missing a declared (entry, fold, seed), a run
whose ``terminal.json`` is absent or not ``complete``, and a prediction set
that does not cover every case of the evaluated fold all raise
:class:`DriverError` instead of scoring a subset.
"""
import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
from artifacts import read_json, sha256, write_json, write_rows  # noqa: E402
import plot_tables  # noqa: E402

EXPECTED_GENES = 18164
K_MAX = 100
ROLES = ('train', 'validation', 'test')
PREDICTION_FILE = 'prediction.npy'
TERMINAL_FILE = 'terminal.json'
CASES_FILE = 'case_ids.json'
IDENTITY_FILE = 'identity.json'
CHECKPOINT_FILE = 'checkpoint.json'


class DriverError(ValueError):
    """A declared evaluation precondition was violated; never score a subset."""


def repo_root():
    for candidate in Path(__file__).resolve().parents:
        if (candidate / '.git').exists():
            return candidate
    raise RuntimeError('repository root not found')


def public_readers():
    """The accepted Stage 2 foundation readers, imported lazily."""
    root = repo_root()
    sys.path.insert(0, str(root / 'stage2/01_bases/output/api'))
    from bases import corpora, evaluation, fold_seed, significance
    return {'corpora': corpora, 'evaluation': evaluation,
            'fold_seed': fold_seed, 'significance': significance}


def _sample_ids(value, label):
    if not isinstance(value, list) or not value:
        raise DriverError(f'{label} must be a nonempty list of sample ids')
    if any(not isinstance(item, str) or not item for item in value):
        raise DriverError(f'{label} must contain nonempty sample-id strings')
    if len(set(value)) != len(value):
        raise DriverError(f'{label} contains duplicate sample ids')
    return value


@dataclass(frozen=True)
class RunArtifacts:
    entry: str
    label: str
    fold: int
    seed: int
    run_dir: Path
    case_ids: list
    prediction: np.ndarray
    terminal: dict
    identity: dict
    history: list


def _read_identity(run_dir):
    path = run_dir / IDENTITY_FILE
    if not path.is_file():
        raise DriverError(f'partially written run: missing {path.name} in {run_dir}')
    try:
        identity = read_json(path)
    except (json.JSONDecodeError, OSError) as exc:
        raise DriverError(f'invalid {path.name} in {run_dir}: {exc}') from exc
    nested = identity.get('identity')
    if not isinstance(nested, dict):
        raise DriverError(f'run identity in {run_dir} lacks an identity mapping')
    return identity


def validate_run(run_dir, expected_case_ids, gene_count=EXPECTED_GENES):
    """Check one run is complete and covers exactly the expected fold cases."""
    run_dir = Path(run_dir)
    terminal_path = run_dir / TERMINAL_FILE
    if not terminal_path.is_file():
        raise DriverError(f'partially written run: missing {TERMINAL_FILE} in {run_dir}')
    try:
        terminal = read_json(terminal_path)
    except (json.JSONDecodeError, OSError) as exc:
        raise DriverError(f'invalid {TERMINAL_FILE} in {run_dir}: {exc}') from exc
    status = terminal.get('status')
    if status != 'complete':
        raise DriverError(f'run status in {run_dir} is {status!r}, not complete')
    cases_path = run_dir / CASES_FILE
    if not cases_path.is_file():
        raise DriverError(f'partially written run: missing {CASES_FILE} in {run_dir}')
    try:
        case_ids = _sample_ids(read_json(cases_path), CASES_FILE)
    except (json.JSONDecodeError, OSError) as exc:
        raise DriverError(f'invalid {CASES_FILE} in {run_dir}: {exc}') from exc
    prediction_path = run_dir / (terminal.get('prediction') or PREDICTION_FILE)
    if not prediction_path.is_file():
        raise DriverError(f'partially written run: missing {prediction_path.name} in {run_dir}')
    prediction = np.load(prediction_path, mmap_mode='r')
    expected = _sample_ids(list(expected_case_ids), 'expected fold cases')
    missing = sorted(set(expected) - set(case_ids))
    foreign = sorted(set(case_ids) - set(expected))
    if missing or foreign:
        raise DriverError(
            f'prediction coverage mismatch in {run_dir}: '
            f'{len(missing)} fold case(s) missing, {len(foreign)} foreign case(s)')
    if prediction.shape != (len(case_ids), gene_count):
        raise DriverError(
            f'prediction shape {prediction.shape} in {run_dir} is not '
            f'{(len(case_ids), gene_count)}')
    if not np.isfinite(prediction).all():
        raise DriverError(f'nonfinite prediction in {run_dir}')
    identity = _read_identity(run_dir)
    history = []
    if (run_dir / CHECKPOINT_FILE).is_file():
        try:
            history = list(read_json(run_dir / CHECKPOINT_FILE).get('history') or [])
        except (json.JSONDecodeError, OSError) as exc:
            raise DriverError(f'invalid {CHECKPOINT_FILE} in {run_dir}: {exc}') from exc
    order = [case_ids.index(sample) for sample in expected]
    return RunArtifacts(entry=run_dir.parent.parent.name,
                        label=str(identity['identity'].get('model') or run_dir.parent.parent.name),
                        fold=int(run_dir.parent.name.removeprefix('fold_')),
                        seed=int(run_dir.name.removeprefix('seed_')),
                        run_dir=run_dir, case_ids=list(expected),
                        prediction=np.asarray(prediction, dtype=np.float64)[order],
                        terminal=terminal, identity=identity, history=history)


def _entries(value):
    entries = [item for item in str(value).split(',') if item]
    if not entries or len(set(entries)) != len(entries):
        raise DriverError('entries must be a nonempty comma-separated unique list')
    return entries


def _integers(value, label):
    try:
        values = [int(item) for item in str(value).split(',') if item != '']
    except ValueError as exc:
        raise DriverError(f'{label} must be comma-separated integers') from exc
    if not values or len(set(values)) != len(values) or min(values) < 0:
        raise DriverError(f'{label} must be a nonempty distinct list of non-negative integers')
    return sorted(values)


def load_run_set(runs_root, entries, fold_rows, seeds, expected_cases, gene_count=EXPECTED_GENES):
    """Load every declared run; a missing member is a partial run set, not a skip."""
    runs_root = Path(runs_root)
    runs = {}
    for entry in entries:
        for fold in fold_rows:
            if fold not in expected_cases:
                raise DriverError(f'no expected case axis for fold row {fold}')
            for seed in seeds:
                run_dir = runs_root / entry / f'fold_{fold}' / f'seed_{seed}'
                if not run_dir.is_dir():
                    raise DriverError(
                        f'incomplete run set: missing run {entry}/fold_{fold}/seed_{seed}')
                run = validate_run(run_dir, expected_cases[fold], gene_count)
                if run.fold != fold or run.seed != seed:
                    raise DriverError(
                        f'run identity in {run_dir} does not match fold {fold} seed {seed}')
                runs[(entry, fold, seed)] = run
    labels = {}
    for run in runs.values():
        if labels.setdefault(run.label, run.entry) != run.entry:
            raise DriverError(f'entry label {run.label!r} is shared by several entries')
    return runs


# ---------------------------------------------------------------------------
# Scoring: every metric value comes from the shared bases.evaluation module.
# ---------------------------------------------------------------------------
EXPECTED_FIELDS = ('sample_id', 'group_id', 'drug_parent', 'dose_um')
DIRECT_COLUMNS = ('PCC', 'Spearman', 'RMSE', 'R2')
CURVE_METRICS = ('Precision@k', 'Spearman@k')
SCALAR_METRICS = ('max_k@P80S40', 'Recall@max_k')


def evaluation_module():
    if 'evaluation' not in sys.modules:
        public_readers()
    from bases import evaluation
    return evaluation


@dataclass(frozen=True)
class RunScore:
    entry: str
    label: str
    fold: int
    seed: int
    expected: list
    well_results: list
    group_results: list
    group_meta: list
    direct: object
    indirect: object
    history: list


def _check_expected(expected):
    if not expected:
        raise DriverError('the evaluated fold has no test cases')
    seen = set()
    for row in expected:
        missing = [field for field in EXPECTED_FIELDS if field not in row]
        if missing:
            raise DriverError(f'expected case row lacks {missing}')
        if row['sample_id'] in seen:
            raise DriverError(f'duplicate expected sample id {row["sample_id"]}')
        seen.add(row['sample_id'])
    return expected


def _undefined_group(k_max, reason, matched=0):
    evaluation = evaluation_module()
    return evaluation.GroupIndirectMetrics(
        k=np.arange(1, k_max + 1, dtype=np.int64),
        precision=np.full(k_max, np.nan), spearman=np.full(k_max, np.nan),
        max_k=float('nan'), recall=float('nan'), matched=int(matched),
        reasons={name: reason for name in evaluation.INDIRECT_METRICS})


def score_run(expected, prediction, measured_delta, truth_rows, *, entry='run',
              label=None, fold=0, seed=0, k_max=K_MAX, history=None):
    """Score one run's persisted predictions against measured Δ and truth rows."""
    evaluation = evaluation_module()
    expected = _check_expected(list(expected))
    prediction = np.asarray(prediction, dtype=np.float64)
    measured_delta = np.asarray(measured_delta, dtype=np.float64)
    if prediction.ndim != 2 or measured_delta.shape != prediction.shape:
        raise DriverError('prediction and measured delta must share the well axis')
    if prediction.shape[0] != len(expected):
        raise DriverError('prediction rows differ from the expected fold cases')
    if prediction.shape[1] != EXPECTED_GENES:
        raise DriverError(f'prediction gene axis is not {EXPECTED_GENES}')
    if not np.isfinite(prediction).all() or not np.isfinite(measured_delta).all():
        raise DriverError('nonfinite prediction or measured delta')
    well_results = []
    for index, row in enumerate(expected):
        result = evaluation.direct_metrics(measured_delta[index], prediction[index])
        well_results.append({'sample_id': row['sample_id'], 'group_id': row['group_id'],
                             'values': result.values, 'reasons': result.reasons})
    positions = {}
    for index, row in enumerate(expected):
        positions.setdefault(row['group_id'], []).append(index)
    group_results, group_meta = [], []
    for group_id in sorted(positions):
        members = positions[group_id]
        row = expected[members[0]]
        truth = truth_rows.row(group_id)
        predicted_mean = prediction[members].mean(axis=0)
        measured_mean = measured_delta[members].mean(axis=0)
        status = str(truth.get('status') or 'eligible')
        if status == 'eligible':
            gate = evaluation.truth_gate(truth['padj'], truth['log2_fold_change'], measured_mean)
            table = evaluation.predicted_rank_table(predicted_mean)
            result = evaluation.indirect_metrics(table, gate, measured_mean, k_max=k_max)
            matched = result.matched
        else:
            gate = None
            result = _undefined_group(k_max, f'truth_{status}')
            matched = 0
        group_results.append(result)
        group_meta.append({'group_id': group_id, 'drug_parent': str(row['drug_parent']),
                           'dose_um': float(row['dose_um']), 'n_wells': len(members),
                           'gated_truth': 0 if gate is None else gate.count,
                           'direction_matched': matched, 'truth_status': status,
                           'truth_reason': str(truth.get('reason') or '')})
    return RunScore(entry=entry, label=label or entry, fold=int(fold), seed=int(seed),
                    expected=expected, well_results=well_results,
                    group_results=group_results, group_meta=group_meta,
                    direct=evaluation.aggregate_direct(
                        [type('Well', (), {'values': row['values'], 'reasons': row['reasons']})
                         for row in well_results]),
                    indirect=evaluation.aggregate_indirect(group_results),
                    history=list(history or []))


def _number(value):
    value = float(value)
    return '' if not np.isfinite(value) else value


def per_well_rows(scores):
    rows = []
    for score in scores:
        for well in score.well_results:
            for metric in DIRECT_COLUMNS:
                value = well['values'][metric]
                rows.append({'entry': score.entry, 'label': score.label, 'fold': score.fold,
                             'seed': score.seed, 'sample_id': well['sample_id'],
                             'group_id': well['group_id'], 'metric': metric,
                             'value': _number(value),
                             'status': 'defined' if np.isfinite(value) else 'undefined',
                             'reason': well['reasons'].get(metric, '')})
    return rows


def per_group_rows(scores):
    rows = []
    for score in scores:
        for result, meta in zip(score.group_results, score.group_meta):
            base = {'entry': score.entry, 'label': score.label, 'fold': score.fold,
                    'seed': score.seed, **meta}
            for metric, curve in (('Precision@k', result.precision),
                                  ('Spearman@k', result.spearman)):
                for index, value in enumerate(curve):
                    rows.append({**base, 'k': index + 1, 'metric': metric,
                                 'value': _number(value),
                                 'status': 'defined' if np.isfinite(value) else 'undefined',
                                 'reason': '' if np.isfinite(value) else result.reasons.get(metric, '')})
            for metric, value in (('max_k@P80S40', result.max_k),
                                  ('Recall@max_k', result.recall)):
                rows.append({**base, 'k': '', 'metric': metric, 'value': _number(value),
                             'status': 'defined' if np.isfinite(value) else 'undefined',
                             'reason': '' if np.isfinite(value) else result.reasons.get(metric, '')})
    return rows


def per_run_rows(scores):
    rows = []
    for score in scores:
        row = {'entry': score.entry, 'label': score.label, 'fold': score.fold,
               'seed': score.seed, 'wells': len(score.well_results),
               'groups': len(score.group_results)}
        for metric in DIRECT_COLUMNS:
            defined = int(score.direct.counts[metric])
            row[metric] = _number(score.direct.values[metric])
            row[f'{metric}_defined'] = defined
            row[f'{metric}_undefined'] = len(score.well_results) - defined
        for metric, attribute in (('max_k@P80S40', 'max_k'), ('Recall@max_k', 'recall')):
            defined = int(score.indirect.counts[metric])
            row[metric] = _number(getattr(score.indirect, attribute))
            row[f'{metric}_defined'] = defined
            row[f'{metric}_undefined'] = len(score.group_results) - defined
        precision_defined = int(score.indirect.counts['Precision@k'][0])
        row['Precision@k_groups'] = precision_defined
        row['Precision@k_groups_undefined'] = len(score.group_results) - precision_defined
        spearman_defined = (int(score.indirect.counts['Spearman@k'][2])
                            if len(score.indirect.k) >= 3 else 0)
        row['Spearman@k_groups'] = spearman_defined
        row['Spearman@k_groups_undefined'] = len(score.group_results) - spearman_defined
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Accepted inputs: fold/seed product, evaluator-only truth, measured delta.
# ---------------------------------------------------------------------------
PER_WELL_FIELDS = ['entry', 'label', 'fold', 'seed', 'sample_id', 'group_id',
                   'metric', 'value', 'status', 'reason']
PER_GROUP_FIELDS = ['entry', 'label', 'fold', 'seed', 'group_id', 'drug_parent',
                    'dose_um', 'n_wells', 'gated_truth', 'direction_matched',
                    'truth_status', 'truth_reason', 'k', 'metric', 'value',
                    'status', 'reason']
PER_RUN_FIELDS = ['entry', 'label', 'fold', 'seed', 'wells', 'groups',
                  'PCC', 'PCC_defined', 'PCC_undefined',
                  'Spearman', 'Spearman_defined', 'Spearman_undefined',
                  'RMSE', 'RMSE_defined', 'RMSE_undefined',
                  'R2', 'R2_defined', 'R2_undefined',
                  'max_k@P80S40', 'max_k@P80S40_defined', 'max_k@P80S40_undefined',
                  'Recall@max_k', 'Recall@max_k_defined', 'Recall@max_k_undefined',
                  'Precision@k_groups', 'Precision@k_groups_undefined',
                  'Spearman@k_groups', 'Spearman@k_groups_undefined']


def fold_cases(cohorts, fold_rows):
    """Test cases and their groups for each evaluated fold, sorted by sample id."""
    cases = {}
    for fold in fold_rows:
        part = cohorts[cohorts.fold.astype(int).eq(int(fold)) & cohorts.role.eq('test')]
        if part.empty:
            raise DriverError(f'fold {fold} has no test cases in the assignment')
        if part.sample_id.duplicated().any():
            raise DriverError(f'fold {fold} repeats a test sample id')
        rows = part.sort_values('sample_id')
        cases[fold] = [{'sample_id': str(row.sample_id), 'group_id': str(row.group_id),
                        'drug_parent': str(row.drug_parent),
                        'dose_um': float(row.dose_um)} for row in rows.itertuples()]
    return cases


class TruthStatistics:
    """The accepted tenv truth release: groups plus the full-axis statistics."""

    def __init__(self, root, allow_fixture=False):
        root = Path(root)
        manifest_path = root / 'manifest.json'
        if not manifest_path.is_file():
            raise DriverError(f'missing accepted truth manifest: {manifest_path}')
        self.root = root
        self.allow_fixture = allow_fixture
        self.manifest = read_json(manifest_path)
        groups_path = root / 'groups.csv'
        if not groups_path.is_file():
            raise DriverError(f'missing accepted truth groups: {groups_path}')
        self.groups = pd.read_csv(groups_path, keep_default_na=False,
                                  dtype={'group_id': str, 'status': str, 'reason': str})
        if self.groups.group_id.duplicated().any():
            raise DriverError('truth groups.csv repeats a group id')
        self.position = {group_id: index for index, group_id in enumerate(self.groups.group_id)}
        gene_count = int(self.manifest.get('gene_count') or 0)
        self.padj = np.load(root / 'padj.npy', mmap_mode='r')
        self.log2_fold_change = np.load(root / 'log2_fold_change.npy', mmap_mode='r')
        expected = (len(self.groups), gene_count)
        if self.padj.shape != expected or self.log2_fold_change.shape != expected:
            raise DriverError(f'truth statistics shapes differ from {expected}')

    def partition(self, cohorts_fold):
        significance = public_readers()['significance']
        assignments = cohorts_fold.loc[:, ['sample_id', 'role', 'group_id',
                                           'drug_parent', 'dose_um']].copy()
        return significance.PartitionTruth(self.root, assignments, 'test',
                                           allow_fixture=self.allow_fixture)

    def rows(self, partition):
        return FoldTruthRows(self, partition)


class FoldTruthRows:
    """Truth rows authorized by one fold's evaluator partition."""

    def __init__(self, statistics, partition):
        self.statistics = statistics
        self.partition = partition
        self.authorized = set(partition.memberships)

    def row(self, group_id):
        if group_id not in self.authorized:
            raise DriverError(f'truth group {group_id!r} is outside the evaluated fold')
        position = self.statistics.position.get(group_id)
        if position is None:
            raise DriverError(f'truth group {group_id!r} is missing from the accepted release')
        record = self.statistics.groups.iloc[position]
        return {'padj': np.asarray(self.statistics.padj[position]),
                'log2_fold_change': np.asarray(self.statistics.log2_fold_change[position]),
                'status': str(record.status), 'reason': str(record.reason)}


def read_layer_rows(path, dataset, positions, gene_count):
    """Read selected rows of one corpus layer in bounded sorted slabs."""
    import h5py
    positions = np.asarray(positions, dtype=np.int64)
    if positions.ndim != 1 or not len(positions):
        raise DriverError('empty corpus row selection')
    values = np.empty((len(positions), gene_count), dtype=np.float64)
    with h5py.File(path, 'r') as handle:
        layer = handle[dataset]
        if layer.shape[1] != gene_count:
            raise DriverError(f'{dataset} gene axis is not {gene_count}')
        for start in range(0, len(positions), 64):
            part = positions[start:start + 64]
            order = np.argsort(part, kind='stable')
            block = np.asarray(layer[part[order], :], dtype=np.float64)
            restored = np.empty_like(block)
            restored[order] = block
            values[start:start + len(part)] = restored
    if not np.isfinite(values).all():
        raise DriverError(f'nonfinite {dataset} value in the evaluated fold')
    return values


def load_measured_delta(corpus, case_rows, gene_count=EXPECTED_GENES):
    """The evaluated fold's measured Δ target, read from the training corpus."""
    corpora = public_readers()['corpora']
    obs, var = corpora.training_header(corpus)
    if len(var) != gene_count:
        raise DriverError(f'corpus gene axis is {len(var)}, not {gene_count}')
    sample_ids = [row['sample_id'] for row in case_rows]
    positions = obs.index.get_indexer(sample_ids)
    if (positions < 0).any():
        missing = [sample for sample, position in zip(sample_ids, positions) if position < 0]
        raise DriverError(f'{len(missing)} evaluated case(s) missing from the corpus: '
                          f'{missing[:3]}')
    return read_layer_rows(corpus, 'layers/delta', positions, gene_count)


# ---------------------------------------------------------------------------
# Figure entries: the recorded per-run values the plot-value tables carry.
# ---------------------------------------------------------------------------
def run_label(fold, seed):
    return f'fold_{int(fold)}/seed_{int(seed)}'


def figure_entries(scores, entries):
    """Plotting entries for the metrics figure; accepts any entry label."""
    plotted = []
    for index, entry in enumerate(entries):
        members = sorted((score for score in scores if score.entry == entry),
                         key=lambda score: (score.fold, score.seed))
        if not members:
            raise DriverError(f'entry {entry!r} has no scored run')
        runs = []
        for score in members:
            direct = {}
            reasons = {}
            for metric in DIRECT_COLUMNS:
                value = float(score.direct.values[metric])
                direct[metric] = None if not np.isfinite(value) else value
                if not np.isfinite(value):
                    reasons[metric] = score.direct.reasons.get(metric, 'undefined')
            indirect = {'k': [int(item) for item in score.indirect.k],
                        'precision': [float(item) for item in score.indirect.precision],
                        'spearman': [float(item) for item in score.indirect.spearman],
                        'max_k': (None if not np.isfinite(score.indirect.max_k)
                                  else float(score.indirect.max_k)),
                        'recall': (None if not np.isfinite(score.indirect.recall)
                                   else float(score.indirect.recall))}
            for metric in SCALAR_METRICS:
                if indirect['max_k' if metric == 'max_k@P80S40' else 'recall'] is None:
                    reasons[metric] = score.indirect.reasons.get(metric, 'undefined')
            for metric in CURVE_METRICS:
                values = indirect['precision' if metric == 'Precision@k' else 'spearman']
                if not np.isfinite(values).any():
                    reasons[metric] = score.indirect.reasons.get(metric, 'undefined')
            runs.append({'run': run_label(score.fold, score.seed), 'direct': direct,
                         'indirect': indirect, 'reasons': reasons})
        plotted.append({'label': members[0].label,
                        'hue': plot_tables.palette_entry_hue(index), 'runs': runs})
    return plotted


def learning_models(scores, entries):
    """Learning-curve models: only entries with a recorded training trajectory."""
    models = []
    for entry in entries:
        members = sorted((score for score in scores if score.entry == entry),
                         key=lambda score: (score.fold, score.seed))
        if not any(score.history for score in members):
            continue
        runs = []
        for score in members:
            label = run_label(score.fold, score.seed)
            if not score.history:
                runs.append({'run': label, 'epoch': [], 'train': [], 'validation': [],
                             'reason': 'no_learning_history'})
                continue
            epoch, train, validation = [], [], []
            for record in score.history:
                try:
                    epoch.append(float(record['epoch']))
                    train.append(float(record['training_loss']))
                    validation.append(float(record['validation_loss']))
                except (KeyError, TypeError, ValueError) as exc:
                    raise DriverError(f'malformed training history in {label}: {exc}') from exc
            runs.append({'run': label, 'epoch': epoch, 'train': train,
                         'validation': validation})
        models.append({'name': members[0].label, 'runs': runs})
    return models


# ---------------------------------------------------------------------------
# Driver: inputs -> tables -> plot-value CSVs + evaluation record.
# ---------------------------------------------------------------------------
def _write_csv(path, records, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_rows(path, records, fields)
    return path


def _verify_identity_hashes(run):
    recorded = run.identity.get('prediction_sha256')
    if recorded:
        actual = sha256(run.run_dir / PREDICTION_FILE)
        if actual != recorded:
            raise DriverError(f'prediction digest differs from identity.json in {run.run_dir}')
    recorded_cases = run.identity.get('case_ids_sha256')
    if recorded_cases and sha256(run.run_dir / CASES_FILE) != recorded_cases:
        raise DriverError(f'case-id digest differs from identity.json in {run.run_dir}')


def run_driver(runs_root, entries, fold_rows, seeds, fold_product, truth, corpus,
               output, figure_data=None, *, allow_fixture=False):
    """Score the declared run set and publish the tables and figure data."""
    output = Path(output)
    fold_seed = public_readers()['fold_seed']
    index = fold_seed.FoldSeedIndex.load(fold_product)
    index.select_folds(max(fold_rows) + 1)
    cohorts = pd.read_csv(Path(fold_product) / 'cohorts.csv',
                          dtype={'sample_id': str, 'drug_parent': str, 'group_id': str,
                                 'role': str})
    expected_cases = fold_cases(cohorts, fold_rows)
    statistics = TruthStatistics(truth, allow_fixture=allow_fixture)
    case_ids = {fold: [row['sample_id'] for row in rows]
                for fold, rows in expected_cases.items()}
    runs = load_run_set(runs_root, entries, fold_rows, seeds, case_ids)
    scores = []
    for (entry, fold, seed), run in sorted(runs.items()):
        _verify_identity_hashes(run)
        partition = statistics.partition(cohorts[cohorts.fold.astype(int).eq(fold)])
        measured = load_measured_delta(corpus, expected_cases[fold])
        scores.append(score_run(expected_cases[fold], run.prediction, measured,
                                statistics.rows(partition), entry=entry,
                                label=run.label, fold=fold, seed=seed,
                                history=run.history))
    well_rows = per_well_rows(scores)
    group_rows = per_group_rows(scores)
    run_rows = per_run_rows(scores)
    output.mkdir(parents=True, exist_ok=True)
    table_rows = {'per_well': (well_rows, PER_WELL_FIELDS),
                  'per_group': (group_rows, PER_GROUP_FIELDS),
                  'per_run': (run_rows, PER_RUN_FIELDS)}
    tables = {name: _write_csv(output / f'{name}.csv', rows, fields)
              for name, (rows, fields) in table_rows.items()}
    entries_values = figure_entries(scores, entries)
    metrics = plot_tables.metrics_values(entries_values)
    figure_files = {'metrics': Path(plot_tables.write_plot_values(
        metrics, output / 'metrics_values.csv'))}
    models = learning_models(scores, entries)
    if models:
        curve = plot_tables.learning_curve_values(models)
        figure_files['learning_curve'] = Path(plot_tables.write_plot_values(
            curve, output / 'learning_curve_values.csv'))
    if figure_data is not None:
        figure_data = Path(figure_data)
        for name, path in list(figure_files.items()):
            destination = figure_data / path.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(path.read_bytes())
            figure_files[name] = destination
    record = {
        'schema': 'stage2-tenv-evaluation-v1',
        'entries': list(entries),
        'seeds': list(seeds),
        'fold_rows': list(fold_rows),
        'truth': {'release': str(Path(truth).resolve()),
                  'manifest_sha256': sha256(Path(truth) / 'manifest.json'),
                  'gene_count': int(statistics.manifest['gene_count']),
                  'fixture': bool(statistics.manifest.get('fixture'))},
        'runs': [{'entry': score.entry, 'label': score.label, 'fold': score.fold,
                  'seed': score.seed, 'wells': len(score.well_results),
                  'groups': len(score.group_results),
                  'terminal_status': runs[(score.entry, score.fold, score.seed)].terminal.get('status'),
                  'prediction_sha256': runs[(score.entry, score.fold, score.seed)].identity.get('prediction_sha256'),
                  'n_gated_truth': int(sum(meta['gated_truth'] for meta in score.group_meta))}
                 for score in scores],
        'tables': {name: {'path': str(path), 'sha256': sha256(path),
                          'rows': len(table_rows[name][0])}
                   for name, path in tables.items()},
        'figures': {name: {'path': str(path), 'sha256': sha256(path)}
                    for name, path in figure_files.items()},
        'plot_value_schema': {
            'source': 'shared bases.plots' if plot_tables.shared_plot_api() else
                      'driver fallback matching the shared schema',
            'columns': list(plot_tables.PLOT_VALUE_COLUMNS),
            'metrics_schema': plot_tables.METRICS_SCHEMA,
            'learning_curve_schema': plot_tables.LEARNING_CURVE_SCHEMA},
        'learning_curve_skipped': None if models else 'no_learnable_entry_trajectory',
        'direct_metrics': list(DIRECT_COLUMNS),
        'indirect_metrics': ['Precision@k', 'Spearman@k', 'max_k@P80S40', 'Recall@max_k'],
        'k_axis': list(range(1, K_MAX + 1)),
    }
    write_json(output / 'evaluation.json', record)
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(
        description='Score persisted tenv predictions against evaluator-only truth (#164).')
    parser.add_argument('--runs-root', required=True)
    parser.add_argument('--entries', required=True)
    parser.add_argument('--fold-rows', default='0')
    parser.add_argument('--seeds', default='17')
    parser.add_argument('--fold-product', required=True)
    parser.add_argument('--truth', required=True)
    parser.add_argument('--corpus', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--figure-data', default=None)
    parser.add_argument('--allow-fixture', action='store_true')
    arguments = parser.parse_args(argv)
    try:
        run_driver(arguments.runs_root, _entries(arguments.entries),
                   _integers(arguments.fold_rows, 'fold rows'),
                   _integers(arguments.seeds, 'seeds'),
                   arguments.fold_product, arguments.truth, arguments.corpus,
                   arguments.output, arguments.figure_data,
                   allow_fixture=arguments.allow_fixture)
    except DriverError as exc:
        print(f'error: {exc}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
