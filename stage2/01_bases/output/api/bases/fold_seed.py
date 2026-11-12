"""Accepted cold-drug assignment and the fold/seed index tables (ADR 0063).

``output/fold_n_seed/`` publishes the accepted tenv cold-drug assignment
(``cohorts.csv``, byte-identical to the accepted split) plus two index tables
that freeze the order of the three outer folds and the three training seeds.
An experiment declares how many folds and seeds it uses and takes that many
rows from the top of each table, so appending an identity never reshuffles an
earlier run.

The assignment is the single source of the cold-drug folds: nothing here reads
the HVG product (ADR 0063).
"""
import hashlib
import json
from numbers import Integral
from pathlib import Path

import pandas as pd

SCHEMA = 'stage2-fold-seed-index-v1'
# Accepted bases-cold-parent-v1 assignment (tenv-hvg manifest, #156/#157 records).
ACCEPTED_COHORTS_SHA256 = 'c1dcce7954b3fd1757dd68fdef4a9e5bf80a0015be6354d29c7be37628a843cc'
COHORTS_FILE = 'cohorts.csv'
FOLD_COLUMNS = ['row', 'fold', 'n_wells', 'n_train', 'n_validation', 'n_test',
                'n_test_drugs']
SEED_COLUMNS = ['row', 'seed']
SEED_ORDER = (17, 42, 2026)
ROLES = ('train', 'validation', 'test')


class FoldSeedError(ValueError):
    """Missing, malformed or digest-mismatched fold/seed product."""


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def _table(path, columns, kind):
    try:
        frame = pd.read_csv(path)
    except FileNotFoundError as exc:
        raise FoldSeedError(f'missing {kind} index table: {path}') from exc
    if list(frame.columns) != columns:
        raise FoldSeedError(f'{Path(path).name} must have columns {columns}')
    if frame.empty:
        raise FoldSeedError(f'{Path(path).name} has no rows')
    if frame.isna().any().any():
        raise FoldSeedError(f'{Path(path).name} has missing values')
    if frame['row'].tolist() != list(range(len(frame))):
        raise FoldSeedError(f'{Path(path).name} rows must be 0..{len(frame) - 1} in order')
    return frame


class FoldSeedIndex:
    """Reader for the accepted fold/seed index tables; selection is by row prefix.

    ``load`` verifies the manifest's declared digest for every file, so a
    tampered or half-published product fails closed.  ``select_folds`` and
    ``select_seeds`` return the first ``n`` table rows unchanged and raise
    :class:`FoldSeedError` for a request beyond the table.
    """

    def __init__(self, root, manifest, folds, seeds):
        self.root = Path(root)
        self.manifest = manifest
        self.folds = folds
        self.seeds = seeds

    @classmethod
    def load(cls, root, verify=True):
        root = Path(root)
        manifest_path = root / 'manifest.json'
        try:
            manifest = json.loads(manifest_path.read_text())
        except FileNotFoundError as exc:
            raise FoldSeedError(f'missing product manifest: {manifest_path}') from exc
        except json.JSONDecodeError as exc:
            raise FoldSeedError(f'manifest.json is not valid JSON: {exc}') from exc
        if manifest.get('schema') != SCHEMA:
            raise FoldSeedError(f'manifest schema is not {SCHEMA}')
        declared = {COHORTS_FILE: (manifest.get('cohorts') or {}).get('sha256')}
        declared.update({name: meta.get('sha256') for name, meta
                         in (manifest.get('files') or {}).items()})
        if verify:
            for name in (COHORTS_FILE, 'folds.csv', 'seeds.csv'):
                path = root / name
                if not path.is_file():
                    raise FoldSeedError(f'missing product file: {name}')
                if declared.get(name) != sha256(path):
                    raise FoldSeedError(f'{name}: manifest digest mismatch')
        folds = _table(root / 'folds.csv', FOLD_COLUMNS, 'fold')
        seeds = _table(root / 'seeds.csv', SEED_COLUMNS, 'seed')
        fold_order = list((manifest.get('folds') or {}).get('order') or [])
        if folds['fold'].tolist() != fold_order:
            raise FoldSeedError('folds.csv order differs from the manifest fold order')
        if fold_order != list(range(len(folds))):
            raise FoldSeedError('folds.csv must order fold identities 0..n-1')
        seed_order = list((manifest.get('seeds') or {}).get('order') or [])
        if seeds['seed'].tolist() != seed_order:
            raise FoldSeedError('seeds.csv order differs from the manifest seed order')
        # The accepted seeds are the first rows; a later append keeps every
        # earlier identity in place, so only the prefix is frozen.
        if tuple(seed_order[:len(SEED_ORDER)]) != SEED_ORDER:
            raise FoldSeedError(f'seeds.csv must start with {list(SEED_ORDER)}')
        if len(seeds) < len(SEED_ORDER):
            raise FoldSeedError(f'seeds.csv must index at least {len(SEED_ORDER)} seeds')
        if (manifest.get('folds') or {}).get('count') != len(folds):
            raise FoldSeedError('manifest fold count differs from folds.csv')
        if (manifest.get('seeds') or {}).get('count') != len(seeds):
            raise FoldSeedError('manifest seed count differs from seeds.csv')
        return cls(root, manifest, folds, seeds)

    def _prefix(self, table, n, kind):
        if isinstance(n, bool) or not isinstance(n, Integral) or n < 0:
            raise FoldSeedError(f'invalid {kind} count: {n!r}')
        if n > len(table):
            raise FoldSeedError(
                f'requested {n} {kind} but the {kind} index table holds {len(table)} rows')
        return table.iloc[:int(n)].copy().reset_index(drop=True)

    def select_folds(self, n):
        """First ``n`` fold rows, in accepted order; reject ``n`` beyond the table."""
        return self._prefix(self.folds, n, 'fold')

    def select_seeds(self, n):
        """First ``n`` seed rows, in accepted order; reject ``n`` beyond the table."""
        return self._prefix(self.seeds, n, 'seed')


def _validate_assignment(frame):
    """Check that ``frame`` is a complete cold-drug assignment; return its summary.

    The accepted assignment splits every drug into exactly one test fold while
    every fold partition still carries all wells and all three roles.  A
    parent-dose group never crosses roles, so the fold split is a drug split.
    """
    required = {'sample_id', 'drug_parent', 'dose_um', 'group_id', 'role', 'fold'}
    if not required <= set(frame.columns):
        raise FoldSeedError(f'assignment requires columns {sorted(required)}')
    frame = frame[list(sorted(required))].copy()
    for column in ('sample_id', 'drug_parent', 'group_id', 'role'):
        frame[column] = frame[column].astype(str)
    if frame.eq('').any().any():
        raise FoldSeedError('assignment has an empty required value')
    if not frame.sample_id.nunique():
        raise FoldSeedError('assignment has no complete well population')
    if not set(frame.role) <= set(ROLES):
        raise FoldSeedError('assignment has an unknown role')
    doses = pd.to_numeric(frame.dose_um, errors='raise')
    if not (doses > 0).all():
        raise FoldSeedError('assignment has a non-positive dose')
    frame['dose_key'] = doses.map(lambda dose: format(float(dose), '.12g'))
    fold_values = pd.to_numeric(frame.fold, errors='raise')
    if not (fold_values == fold_values.astype(int)).all():
        raise FoldSeedError('assignment fold must be an integer')
    frame['fold_num'] = fold_values.astype(int)
    group_keys = frame.groupby('group_id')[['drug_parent', 'dose_key']].nunique()
    if (group_keys > 1).any().any():
        raise FoldSeedError('a group identity spans several parent-dose pairs')
    pair_groups = frame.groupby(['drug_parent', 'dose_key']).group_id.nunique()
    if (pair_groups > 1).any():
        raise FoldSeedError('a parent-dose pair has several group identities')
    fold_ids = sorted(int(value) for value in set(frame.fold_num))
    if fold_ids != list(range(len(fold_ids))):
        raise FoldSeedError('fold identities must be 0..n-1')
    wells = set(frame.sample_id)
    all_drugs = set(frame.drug_parent)
    rows, test_sets = [], {}
    for fold in fold_ids:
        part = frame[frame.fold_num.eq(fold)]
        if part.sample_id.duplicated().any() or set(part.sample_id) != wells:
            raise FoldSeedError(f'fold {fold} must partition every well exactly once')
        if set(part.role) != set(ROLES):
            raise FoldSeedError(f'fold {fold} requires populated {ROLES} roles')
        cross_roles = part.groupby(['drug_parent', 'dose_key']).role.nunique().gt(1).any()
        if cross_roles:
            raise FoldSeedError('parent-dose group crosses roles')
        counts = part.role.value_counts().to_dict()
        test_drugs = set(part.loc[part.role.eq('test'), 'drug_parent'])
        if not test_drugs:
            raise FoldSeedError(f'fold {fold} has no test drug')
        test_sets[fold] = test_drugs
        rows.append({'row': fold, 'fold': fold, 'n_wells': len(part),
                     'n_train': int(counts.get('train', 0)),
                     'n_validation': int(counts.get('validation', 0)),
                     'n_test': int(counts.get('test', 0)),
                     'n_test_drugs': len(test_drugs)})
    covered = set()
    for fold, test_drugs in test_sets.items():
        overlap = covered & test_drugs
        if overlap:
            raise FoldSeedError(
                f'test drugs of fold {fold} also test elsewhere: {sorted(overlap)[:3]}')
        covered |= test_drugs
    if covered != all_drugs:
        raise FoldSeedError('test folds must cover every drug exactly once')
    if len(frame) != len(wells) * len(fold_ids):
        raise FoldSeedError('assignment must repeat every well exactly once per fold')
    return {'folds': rows, 'n_wells': len(wells), 'n_drugs': len(all_drugs),
            'n_groups': int(frame.group_id.nunique()), 'n_folds': len(fold_ids),
            'n_rows': len(frame)}


def build_fold_seed(cohorts, output, accepted_sha256=ACCEPTED_COHORTS_SHA256):
    """Re-home the accepted assignment and publish the two index tables.

    The assignment is copied byte-identically; the index tables and the
    manifest are regenerated deterministically from it.  The build refuses a
    source whose digest is not the accepted one, so an accidental re-freeze
    cannot silently replace the frozen folds.
    """
    import io

    source, output = Path(cohorts), Path(output)
    raw = source.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != accepted_sha256:
        raise FoldSeedError(
            f'source assignment digest {digest} is not the accepted digest {accepted_sha256}')
    try:
        frame = pd.read_csv(io.BytesIO(raw), keep_default_na=False,
                            dtype={'sample_id': str, 'drug_parent': str,
                                   'group_id': str})
    except pd.errors.ParserError as exc:
        raise FoldSeedError(f'assignment is not a valid CSV: {exc}') from exc
    summary = _validate_assignment(frame)
    output.mkdir(parents=True, exist_ok=True)
    (output / COHORTS_FILE).write_bytes(raw)
    fold_frame = pd.DataFrame(summary['folds'], columns=FOLD_COLUMNS)
    fold_frame.to_csv(output / 'folds.csv', index=False)
    seed_frame = pd.DataFrame({'row': range(len(SEED_ORDER)), 'seed': SEED_ORDER})
    seed_frame.to_csv(output / 'seeds.csv', index=False)
    files = {name: {'sha256': sha256(output / name),
                    'bytes': (output / name).stat().st_size,
                    'rows': len(pd.read_csv(output / name))}
             for name in ('folds.csv', 'seeds.csv')}
    manifest = {
        'schema': SCHEMA,
        'cohorts': {'file': COHORTS_FILE, 'sha256': digest, 'bytes': len(raw),
                    'n_wells': summary['n_wells'], 'n_drugs': summary['n_drugs'],
                    'n_groups': summary['n_groups'], 'n_folds': summary['n_folds'],
                    'algorithm': 'SHA256(seed|parent) sorted buckets; 10% remaining parents validation',
                    'seed': 'bases-cold-parent-v1'},
        'files': files,
        'folds': {'file': 'folds.csv', 'count': summary['n_folds'],
                  'order': [row['fold'] for row in summary['folds']]},
        'seeds': {'file': 'seeds.csv', 'count': len(SEED_ORDER),
                  'order': list(SEED_ORDER)},
        'selection': 'take the first n rows of each index table; n beyond the table is an error',
        'recipe': 'frozen cold-drug folds re-homed from the accepted bases-cold-parent-v1 assignment; index tables order identities by row (ADR 0063)',
    }
    (output / 'manifest.json').write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    return manifest


def _fold_counts(frame):
    """Recompute each fold's role partition from an assignment frame."""
    numbered = frame.assign(
        fold_num=pd.to_numeric(frame.fold, errors='raise').astype(int),
        role=frame.role.astype(str), drug_parent=frame.drug_parent.astype(str))
    counts = {}
    for fold, part in numbered.groupby('fold_num'):
        roles = part.role.value_counts().to_dict()
        counts[int(fold)] = {
            'row': int(fold), 'fold': int(fold), 'n_wells': int(len(part)),
            'n_train': int(roles.get('train', 0)),
            'n_validation': int(roles.get('validation', 0)),
            'n_test': int(roles.get('test', 0)),
            'n_test_drugs': int(part.loc[part.role.eq('test'), 'drug_parent'].nunique()),
        }
    return counts


def verify_fold_seed(root, expected_cohorts_sha256=ACCEPTED_COHORTS_SHA256):
    """Re-verify the published product offline; return a PASS/FAIL report.

    Checks the accepted assignment digest and its complete fold partition, the
    two index tables' row order, the frozen cold-drug split (disjoint test sets
    whose union is every drug exactly once), each fold's role counts against its
    well count, the manifest's declared digests and counts, and the reader's
    prefix/reject contract.  Every check reports PASS or FAIL; nothing raises.
    """
    import io

    root = Path(root)
    checks, files = [], {}

    def add(name, ok, detail=None):
        entry = {'check': name, 'status': 'PASS' if ok else 'FAIL'}
        if detail is not None:
            entry['detail'] = detail
        checks.append(entry)

    manifest = None
    try:
        manifest = json.loads((root / 'manifest.json').read_text())
    except (OSError, json.JSONDecodeError) as exc:
        add('manifest: schema and declared digests', False, {'error': str(exc)})
    if manifest is not None:
        add('manifest: schema is the fold/seed index schema',
            manifest.get('schema') == SCHEMA, {'schema': manifest.get('schema')})

    for name in (COHORTS_FILE, 'folds.csv', 'seeds.csv'):
        try:
            raw = (root / name).read_bytes()
            files[name] = {'sha256': hashlib.sha256(raw).hexdigest(), 'bytes': len(raw)}
        except OSError as exc:
            files[name] = {'sha256': None, 'bytes': 0, 'error': str(exc)}
    declared = {}
    if manifest is not None:
        declared[COHORTS_FILE] = (manifest.get('cohorts') or {}).get('sha256')
        for name, meta in (manifest.get('files') or {}).items():
            declared[name] = (meta or {}).get('sha256')
    digests_ok = (manifest is not None and files[COHORTS_FILE]['sha256'] is not None
                  and set(declared) == {COHORTS_FILE, 'folds.csv', 'seeds.csv'}
                  and all(declared.get(name) == files[name]['sha256']
                          for name in files))
    add('manifest: declared file digests match the product bytes', digests_ok,
        {'declared': sorted(declared), 'actual': {name: value['sha256']
                                                  for name, value in files.items()}})

    cohorts_sha = files[COHORTS_FILE]['sha256']
    add('cohorts.csv: accepted assignment digest',
        cohorts_sha is not None and cohorts_sha == expected_cohorts_sha256,
        {'expected_sha256': expected_cohorts_sha256, 'actual_sha256': cohorts_sha})

    frame = None
    cohorts_summary = {'sha256': cohorts_sha,
                       'expected_sha256': expected_cohorts_sha256}
    try:
        frame = pd.read_csv(io.BytesIO((root / COHORTS_FILE).read_bytes()),
                            keep_default_na=False,
                            dtype={'sample_id': str, 'drug_parent': str,
                                   'group_id': str})
    except Exception as exc:  # noqa: BLE001 - report the failure, never raise
        add('cohorts.csv: parses as the accepted assignment table', False,
            {'error': str(exc)})
    if frame is not None:
        cohorts_summary.update({
            'n_rows': int(len(frame)),
            'n_wells': int(frame.sample_id.astype(str).nunique()),
            'n_drugs': int(frame.drug_parent.astype(str).nunique()),
            'n_groups': int(frame.group_id.astype(str).nunique()),
        })

    assignment = None
    if frame is not None:
        try:
            assignment = _validate_assignment(frame)
            add('cohorts.csv: complete well, role and group assignment', True,
                {'n_wells': assignment['n_wells'], 'n_folds': assignment['n_folds'],
                 'n_groups': assignment['n_groups']})
        except Exception as exc:  # noqa: BLE001 - report the failure, never raise
            add('cohorts.csv: complete well, role and group assignment', False,
                {'error': str(exc)})

    if frame is not None:
        try:
            numbered = frame.assign(
                fold_num=pd.to_numeric(frame.fold, errors='raise').astype(int),
                role=frame.role.astype(str), drug_parent=frame.drug_parent.astype(str))
            test_sets = {int(fold): set(part.loc[part.role.eq('test'), 'drug_parent'])
                         for fold, part in numbered.groupby('fold_num')}
            all_drugs = set(numbered.drug_parent)
            covered, overlap = set(), set()
            for test_drugs in test_sets.values():
                overlap |= covered & test_drugs
                covered |= test_drugs
            once = sum(len(value) for value in test_sets.values()) == len(all_drugs)
            add('folds: test drugs are disjoint across folds and cover every drug exactly once',
                not overlap and once and covered == all_drugs,
                {'test_drugs_per_fold': {str(fold): len(value)
                                         for fold, value in sorted(test_sets.items())},
                 'n_drugs': len(all_drugs), 'overlap': sorted(overlap)[:3],
                 'covered': len(covered)})
        except Exception as exc:  # noqa: BLE001
            add('folds: test drugs are disjoint across folds and cover every drug exactly once',
                False, {'error': str(exc)})

    actual = {}
    if frame is not None:
        try:
            actual = _fold_counts(frame)
        except Exception:  # noqa: BLE001
            actual = {}
    folds_rows, folds_frame = [], None
    try:
        folds_frame = pd.read_csv(io.BytesIO((root / 'folds.csv').read_bytes()))
    except Exception as exc:  # noqa: BLE001
        add('folds.csv: parses with the index columns', False, {'error': str(exc)})
    if folds_frame is not None:
        columns_ok = list(folds_frame.columns) == FOLD_COLUMNS
        rows_ok = folds_frame.get('row', pd.Series(dtype=int)).tolist() == list(range(len(folds_frame)))
        add('folds.csv: rows 0..n-1 order the frozen folds', columns_ok and rows_ok,
            {'columns': list(folds_frame.columns), 'rows': rows_ok})
        mismatches = []
        sum_ok = True
        if columns_ok:
            for record in folds_frame.to_dict('records'):
                fold = int(record['fold'])
                expected = actual.get(fold)
                if expected is None or any(int(record[key]) != expected[key]
                                           for key in FOLD_COLUMNS if key != 'row'):
                    mismatches.append({'fold': fold, 'declared': record,
                                       'assignment': expected})
                if int(record['n_train']) + int(record['n_validation']) + int(record['n_test']) \
                        != int(record['n_wells']):
                    sum_ok = False
        add('folds.csv: role counts match the assignment and sum to the well count',
            not mismatches and sum_ok and set(int(row['fold']) for row in
                                              folds_frame.to_dict('records')) == set(actual),
            {'mismatches': mismatches[:3], 'role_sums_equal_well_count': sum_ok})
        fold_order = list((manifest or {}).get('folds', {}).get('order') or [])
        fold_count = (manifest or {}).get('folds', {}).get('count')
        add('folds.csv: manifest order and count match the table',
            fold_order == [int(value) for value in folds_frame['fold']]
            and fold_count == len(folds_frame),
            {'declared_order': fold_order, 'declared_count': fold_count,
             'table_order': [int(value) for value in folds_frame['fold']]})
        folds_rows = [{key: int(value) for key, value in record.items()}
                      for record in folds_frame.to_dict('records')]

    seeds_rows, seeds_frame = [], None
    try:
        seeds_frame = pd.read_csv(io.BytesIO((root / 'seeds.csv').read_bytes()))
    except Exception as exc:  # noqa: BLE001
        add('seeds.csv: parses with the index columns', False, {'error': str(exc)})
    if seeds_frame is not None:
        seed_order = list((manifest or {}).get('seeds', {}).get('order') or [])
        ok = (list(seeds_frame.columns) == SEED_COLUMNS
              and seeds_frame.get('row', pd.Series(dtype=int)).tolist() == list(range(len(seeds_frame)))
              and seeds_frame['seed'].tolist() == list(SEED_ORDER)
              and seed_order == list(SEED_ORDER))
        add('seeds.csv: rows 0-2 order the accepted seeds 17, 42, 2026', ok,
            {'order': [int(value) for value in seeds_frame['seed']],
             'declared_order': seed_order})
        seeds_rows = [{key: int(value) for key, value in record.items()}
                      for record in seeds_frame.to_dict('records')]

    try:
        index = FoldSeedIndex.load(root)
        prefix_ok = (index.select_folds(len(index.folds))['row'].tolist()
                     == index.folds['row'].tolist()
                     and index.select_seeds(len(index.seeds))['row'].tolist()
                     == index.seeds['row'].tolist())
        rejected = []
        for call, kind in ((index.select_folds, 'fold'), (index.select_seeds, 'seed')):
            try:
                call(1 + (len(index.folds) if kind == 'fold' else len(index.seeds)))
            except FoldSeedError:
                rejected.append(kind)
        add('selection: row prefix returns the first n rows and rejects beyond the table',
            bool(prefix_ok) and set(rejected) == {'fold', 'seed'},
            {'rejected_beyond_table': sorted(rejected)})
    except Exception as exc:  # noqa: BLE001
        add('selection: row prefix returns the first n rows and rejects beyond the table',
            False, {'error': str(exc)})

    if manifest is not None and assignment is not None:
        recorded = manifest.get('cohorts') or {}
        totals_ok = all(recorded.get(key) == assignment[key]
                        for key in ('n_wells', 'n_drugs', 'n_groups', 'n_folds'))
        index_counts_ok = ((manifest.get('folds') or {}).get('count') == assignment['n_folds']
                           and (manifest.get('seeds') or {}).get('count') == len(SEED_ORDER))
        add('manifest: recorded assignment totals and index counts match',
            totals_ok and index_counts_ok,
            {'totals': {key: recorded.get(key) for key in
                        ('n_wells', 'n_drugs', 'n_groups', 'n_folds')},
             'assignment': {key: assignment[key] for key in
                            ('n_wells', 'n_drugs', 'n_groups', 'n_folds')}})

    failures = [entry['check'] for entry in checks if entry['status'] != 'PASS']
    return {
        'result': 'PASS' if checks and not failures else 'FAIL',
        'schema': SCHEMA,
        'files': files,
        'cohorts': cohorts_summary,
        'folds': folds_rows,
        'seeds': seeds_rows,
        'checks': checks,
        'failures': failures,
    }
