"""The one live derivation for in-vitro bulk cardiomyocyte deposits (ticket #188).

ADR 0066 / spec D9 fixes one DEG implementation, ``stage2_util/derivation/disease_deg.R`` (limma:
voom for integer counts, ``lmFit`` for normalised log-scale expression), and
ADR 0062 fixes one place for the science: the step's public API boundary. This
module is the single *source adapter* the six isogenic ADR 0030 items share, so
the derivation rules live once and are unit-tested while each item's
``script/build_deg.py`` only binds its deposited bytes and declares its contrast
design (which arm is the isogenic control, which strata are pooled, why):

* ``read_series_matrix`` — parse a GEO series matrix's ``!Sample_*`` block into a
  per-sample annotation frame (title, accession, description, and every
  ``key: value`` characteristic). The series matrix is the deposit's own sample
  annotation; no genotype/arm is ever hand-entered.
* ``match_samples`` — bind the deposited expression matrix's column names to
  those samples by name (never by position), with a recorded basis and a
  fail-closed bijection check.
* ``read_deposited_matrix`` — the deposited feature x sample matrix, as float64,
  with the BOM/separator/versioned-Ensembl quirks handled once.
* ``integerise_counts`` — **the one shared helper** that turns a non-integer
  count-scale deposit (RSEM ``expected_count``, kallisto/Salmon estimated counts)
  into integers for the strict-count voom path. It records the rule, the maximum
  absolute deviation, the library-size scale and every gene whose rounded row
  sum is zero, because the count model drops all-zero rows and that loss must be
  reported, never silently filtered. The shared ``disease_deg.R`` count contract
  is not relaxed: #186/#187/#189 recorded products depend on it.
* ``derive_item_deg`` — orchestrates the above and runs the one live DEG worker,
  writing the item's ``local/deg/`` matrix, sample sheet, full DEG table and the
  tracked summary JSON. No science below this line in an item script.

Contrast direction is always ``numerator - denominator`` (positive log2FC = up in
the case/isogenic-mutant arm named first).
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ITEM_DEG_SCHEMA = 'stage2-item-deg-summary-v1'
# A bulk RNA-seq count-scale library sums to millions of counts. A deposit whose
# median column sum is far below that is not on a count scale (it is already
# normalised) and must take the ``lmFit`` expression path instead of being
# integerised: this guard makes the wrong declaration fail instead of rounding a
# normalised matrix into a meaningless count matrix.
MIN_COUNT_SCALE_LIBRARY_SUM = 1e5
MATRIX_COLUMNS = ('source_gene_id', 'gene_symbol')
_CHARACTERISTIC = re.compile(r'^!Sample_characteristics_ch1\t(.+)$')
_SAMPLE_FIELD = re.compile(r'^!Sample_([a-zA-Z0-9_]+)\t(.+)$')
_NAME_PREFIX = re.compile(r'^(library|label|sample)\s+name:\s*', re.IGNORECASE)
_BRACKET = re.compile(r'\[([^\[\]]+)\]')


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _sanitise(name):
    return re.sub(r'[^0-9a-zA-Z]+', '_', name).strip('_').lower()


def _unquote(field):
    field = field.strip()
    if len(field) >= 2 and field[0] == '"' and field[-1] == '"':
        return field[1:-1]
    return field


def _open_text(path):
    path = Path(path)
    if path.name.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8-sig', newline='')
    return path.open('r', encoding='utf-8-sig', newline='')


def read_series_matrix(path):
    """Parse a GEO series matrix into a per-sample annotation frame.

    The returned frame is indexed by GEO accession and carries ``title``,
    ``source_name_ch1``, ``description`` plus one column per characteristic key
    (e.g. ``genotype``, ``cell_type``, ``batch``). A repeated characteristic key
    with different values for one sample is an error, not a last-one-wins.
    """
    fields = {}
    order = []
    with _open_text(path) as handle:
        for line in handle:
            line = line.rstrip('\r\n')
            if line.startswith('!series_matrix_table_begin'):
                break
            if line.startswith('!series_matrix_table'):
                break
            if line.startswith('!Sample_characteristics_ch1\t'):
                values = [_unquote(item) for item in line.split('\t')[1:]]
                if not order:
                    raise ValueError(f'{path}: sample fields appear after characteristics')
                per_sample = []
                for value in values:
                    if ':' not in value:
                        raise ValueError(
                            f'{path}: characteristic without "key: value": {value!r}')
                    key, text = value.split(':', 1)
                    per_sample.append((_sanitise(key), text.strip()))
                for position, (key, text) in enumerate(per_sample):
                    target = fields[order[position]]
                    if key in target and target[key] != text:
                        # A repeated characteristic key with a different value is a
                        # conflicting design fact: fail closed, never last-one-wins.
                        raise ValueError(
                            f'{path}: conflicting characteristic {key!r} for '
                            f'{order[position]}: {target[key]!r} vs {text!r}')
                    target[key] = text
                continue
            match = _SAMPLE_FIELD.match(line)
            if not match:
                continue
            key = _sanitise(match.group(1))
            values = [_unquote(item) for item in match.group(2).split('\t')]
            if not order:
                order = values  # first !Sample_ row defines the column order
                fields = {accession: {} for accession in order}
            if len(values) != len(order):
                raise ValueError(f'{path}: inconsistent !Sample_ row {key!r}')
            for accession, value in zip(order, values):
                if key in fields[accession]:
                    # Extraction/data-processing protocol rows legitimately repeat;
                    # the first occurrence is kept and none of them carries a design
                    # fact. Characteristic keys are handled above and may not conflict.
                    continue
                fields[accession][key] = value
    if not fields:
        raise ValueError(f'{path}: no !Sample_ rows found')
    accessions = [fields[accession].get('geo_accession') for accession in order]
    if any(accession is None for accession in accessions):
        raise ValueError(f'{path}: samples without a geo_accession')
    series = pd.DataFrame.from_dict(fields, orient='index')
    series.index = pd.Index(accessions, name='gsm')
    if 'title' not in series:
        raise ValueError(f'{path}: series matrix lacks !Sample_title')
    for column in series.columns:
        series[column] = series[column].astype(object)
    return series


def normalise_token(text):
    """Comparison key for a sample name: lowercase alphanumerics only.

    A leading ``s`` followed by a digit is dropped so a GSM title token
    ``s20020813E7`` matches the deposited column ``2-0020813E7`` (the same
    sample name with separator punctuation and a one-letter prefix).
    """
    token = re.sub(r'[^0-9a-zA-Z]+', '', str(text)).lower()
    if len(token) > 1 and token[0] == 's' and token[1].isdigit():
        token = token[1:]
    return token


def _sample_keys(series, basis):
    keys = {}
    for gsm, row in series.iterrows():
        if basis == 'title_bracket':
            found = _BRACKET.search(str(row['title']))
            if not found:
                raise ValueError(f'sample {gsm}: title has no [...] token: {row["title"]!r}')
            raw = found.group(1)
        elif basis == 'description':
            if 'description' not in row or pd.isna(row['description']):
                raise ValueError(f'sample {gsm}: match basis "description" but the '
                                 'series matrix has no !Sample_description')
            raw = _NAME_PREFIX.sub('', str(row['description']))
        elif basis == 'title':
            raw = str(row['title'])
        else:
            raise ValueError(f'unknown match basis {basis!r}')
        keys[gsm] = normalise_token(raw)
    return keys


def match_samples(matrix_columns, series, basis):
    """Bind matrix columns to series-matrix samples by name; fail closed.

    Returns a frame indexed by the matrix column name in matrix column order,
    carrying the series annotation columns. Every matrix column must match
    exactly one series sample and every series sample must be matched, so a
    silently dropped or duplicated library is impossible.
    """
    keys = _sample_keys(series, basis)
    by_key = {}
    for gsm, key in keys.items():
        if key in by_key:
            raise ValueError(f'series samples {by_key[key]} and {gsm} share the name {key!r}')
        by_key[key] = gsm
    rows, seen = [], set()
    for column in matrix_columns:
        key = normalise_token(column)
        gsm = by_key.get(key)
        if gsm is None:
            raise ValueError(f'matrix column {column!r} matches no series sample '
                             f'(basis {basis!r})')
        if gsm in seen:
            raise ValueError(f'series sample {gsm} matched by two matrix columns')
        seen.add(gsm)
        payload = {'sample': column, 'gsm': gsm}
        payload.update(series.loc[gsm].to_dict())
        rows.append(pd.Series(payload))
    unmatched = sorted(set(series.index) - seen)
    if unmatched:
        raise ValueError(f'series samples with no matrix column: {unmatched}')
    return pd.DataFrame(rows).set_index('sample', drop=False)


def _derived_values(series, specification):
    """Evaluate a ``derived`` covariate declaration against the annotation frame."""
    values = []
    for gsm, row in series.iterrows():
        text = str(row[specification['from']])
        found = re.search(specification['regex'], text)
        if not found:
            raise ValueError(f'sample {gsm}: covariate regex found nothing in {text!r}')
        values.append(found.group(1))
    return pd.Series(values, index=series.index)


def build_sample_sheet(matrix_columns, series, spec):
    """Build the fit design sheet: sample, group, and declared blocking covariates."""
    matched = match_samples(matrix_columns, series, spec['sample_match'])
    group_source = spec['group_source']
    if isinstance(group_source, (list, tuple)):
        # A composite stratum (e.g. format x cell type x genotype) is one declared
        # group source whose raw value joins the named characteristics in order.
        raw = matched[list(group_source)].astype(str).agg('|'.join, axis=1)
    else:
        raw = matched[group_source].astype(str)
    groups = raw.map(spec['group_map'])
    missing = sorted(set(groups[groups.isna()].index.tolist()))
    if missing:
        raise ValueError(f'unmapped {group_source} value(s): {sorted(set(raw[groups.isna()]))}')
    sheet = pd.DataFrame({'sample': matched['sample'].to_numpy(),
                          'group': groups.to_numpy(dtype=object)})
    for key, column in spec.get('covariates', {}).items():
        sheet[column] = matched[key].astype(str).to_numpy()
    for column, specification in spec.get('derived_covariates', {}).items():
        derived = _derived_values(matched, specification)
        sheet[column] = derived.to_numpy(dtype=object)
    sheet['gsm'] = matched['gsm'].to_numpy()
    sheet['title'] = matched['title'].to_numpy()
    for label in ('numerator', 'denominator'):
        for contrast in spec['contrasts']:
            level = contrast[label]
            if not (sheet.group == level).any():
                raise ValueError(f'contrast arm {level!r} has no sample')
    return sheet, matched


def read_deposited_matrix(path, sep):
    """Read a deposited expression matrix: first column feature id, rest samples."""
    frame = pd.read_csv(path, sep=sep, index_col=0, encoding='utf-8-sig')
    if frame.shape[1] < 2:
        raise ValueError(f'{path}: matrix needs at least two sample columns')
    frame = frame.astype(np.float64)
    if frame.index.has_duplicates:
        raise ValueError(f'{path}: duplicate feature identifiers')
    if not np.isfinite(frame.to_numpy()).all():
        raise ValueError(f'{path}: matrix carries non-finite values')
    if (frame.to_numpy() < 0).any():
        raise ValueError(f'{path}: matrix carries negative values')
    return frame


def integerise_counts(frame):
    """Round a non-integer count-scale matrix to integers, recording every loss.

    The rule is ``numpy.rint`` (round to nearest, ties to even) applied to every
    value of a deposit whose values are on a count scale but not integers (RSEM
    ``expected_count``, kallisto/Salmon estimated counts). Returns
    ``(integer_frame, stats)`` where ``stats`` carries the rule, the maximum and
    mean absolute deviation, the library-size scale before/after, and the number
    of genes whose rounded row is all-zero (the count model drops all-zero rows,
    so those genes leave the test and that loss is reported, never filtered
    silently).
    """
    values = np.asarray(frame, dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError('count-scale integerisation requires finite values')
    if (values < 0).any():
        raise ValueError('count-scale integerisation requires non-negative values')
    rounded = np.rint(values)
    deviation = np.abs(values - rounded)
    integer = pd.DataFrame(rounded.astype(np.int64), index=frame.index,
                           columns=frame.columns)
    before = values.sum(axis=0)
    after = rounded.sum(axis=0)
    row_before = values.sum(axis=1) > 0
    row_after = rounded.sum(axis=1) > 0
    stats = {
        'rule': 'round to nearest integer, ties to even (numpy.rint)',
        'n_values': int(values.size),
        'n_non_integer_values': int((values != rounded).sum()),
        'max_abs_deviation': float(deviation.max()),
        'mean_abs_deviation': float(deviation.mean()),
        'library_size_before': {'min': float(before.min()), 'median': float(np.median(before)),
                                'max': float(before.max())},
        'library_size_after': {'min': int(after.min()), 'median': int(np.median(after)),
                               'max': int(after.max())},
        'n_genes_all_zero_before': int((~row_before).sum()),
        'n_genes_all_zero_after': int((~row_after).sum()),
        'n_genes_dropped_by_rounding': int((row_before & ~row_after).sum()),
        'dropped_gene_ids': frame.index[row_before & ~row_after].astype(str).tolist(),
    }
    if stats['library_size_before']['median'] < MIN_COUNT_SCALE_LIBRARY_SUM:
        raise ValueError(
            'declared count-scale deposit has a median library sum of '
            f"{stats['library_size_before']['median']:.0f}: not a count-scale matrix "
            '(use the lmFit expression path for a normalised deposit)')
    if stats['max_abs_deviation'] >= 1.0:  # unreachable for np.rint; kept as a contract
        raise ValueError('integerisation deviation >= 1 count: not a count-scale deposit')
    return integer, stats


def validate_integral_counts(frame):
    """Assert an as-deposited matrix really is integer counts; return its scale."""
    values = np.asarray(frame, dtype=np.float64)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('count deposit requires finite non-negative values')
    if not np.equal(values, np.rint(values)).all():
        raise ValueError('deposit declared integer counts but carries non-integer values')
    totals = values.sum(axis=0)
    if float(np.median(totals)) < MIN_COUNT_SCALE_LIBRARY_SUM:
        raise ValueError(
            'declared count deposit has a median library sum of '
            f'{np.median(totals):.0f}: not a count-scale matrix '
            '(use the lmFit expression path for a normalised deposit)')
    return {'mode': 'integer-as-deposited',
            'n_values': int(values.size),
            'library_size': {'min': int(totals.min()), 'median': int(np.median(totals)),
                             'max': int(totals.max())}}


def _write_tsv_gz(frame, path):
    """Write a TSV gzip with a fixed mtime so the bytes are reproducible."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = frame.to_csv(sep='\t', index=False, lineterminator='\n').encode()
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode='wb', mtime=0) as handle:
        handle.write(payload)
    path.write_bytes(buffer.getvalue())


def _relative_to_item(item_dir, path):
    try:
        return Path(path).resolve().relative_to(Path(item_dir).resolve()).as_posix()
    except ValueError:
        return str(path)


def derive_item_deg(item_dir, spec):
    """Derive one item's local expression matrix, DEG table and tracked summary."""
    from . import disease_deg

    item_dir = Path(item_dir)
    item = item_dir.name
    local = item_dir / 'local' / 'deg'
    local.mkdir(parents=True, exist_ok=True)
    matrix_path = item_dir / spec['matrix']
    series_path = item_dir / spec['series_matrix']
    frame = read_deposited_matrix(matrix_path, spec['sep'])
    series = read_series_matrix(series_path)
    sheet, matched = build_sample_sheet(list(frame.columns), series, spec)
    samples_path = local / f'{item}_samples.csv'
    sheet.to_csv(samples_path, index=False)
    kind = spec['kind']
    if kind == 'counts':
        if spec.get('integerisation') == 'round_half_to_even':
            matrix, values_record = integerise_counts(frame)
        else:
            values_record = validate_integral_counts(frame)
            matrix = frame.astype(np.int64)
    elif kind == 'expression':
        if spec.get('integerisation'):
            raise ValueError('a normalised deposit is not integerised')
        matrix = frame
        values_record = {'mode': 'normalised-expression-as-deposited'}
    else:
        raise ValueError(f'unknown deposition kind {kind!r}')
    expression = pd.DataFrame({'source_gene_id': matrix.index.astype(str),
                               'gene_symbol': matrix.index.astype(str)})
    for column in matrix.columns:
        expression[column] = matrix[column].to_numpy()
    expr_path = local / f'{item}_expr.tsv.gz'
    _write_tsv_gz(expression, expr_path)
    design = {'kind': 'counts' if kind == 'counts' else 'expression',
              'formula': spec['formula'], 'factor': 'group',
              'contrasts': spec['contrasts']}
    deg_path = local / f'{item}_deg.csv'
    fit = disease_deg.run_deg_worker(expr_path, samples_path, design, deg_path)
    table = pd.read_csv(deg_path)
    missing = set(disease_deg.LOCAL_DEG_COLUMNS) - set(table.columns)
    if missing:
        raise ValueError(f'{item}: DEG worker output lacks {sorted(missing)}')
    produced = sorted(table.contrast.astype(str).unique().tolist())
    declared = sorted(contrast['label'] for contrast in spec['contrasts'])
    if produced != declared:
        raise ValueError(f'{item}: worker produced {produced}, declared {declared}')
    if table[['log2fc', 'stat', 'p_value', 'base_mean']].isna().all().any():
        raise ValueError(f'{item}: a DEG statistic column is entirely missing')
    case_arms = {contrast['numerator'] for contrast in spec['contrasts']}
    control_arms = {contrast['denominator'] for contrast in spec['contrasts']}
    if case_arms & control_arms:
        raise ValueError(f'{item}: an arm is both case and control')
    group_counts = {str(key): int(value)
                    for key, value in sheet.group.value_counts().sort_index().items()}
    contrasts = {}
    for contrast in spec['contrasts']:
        contrasts[contrast['label']] = {
            'n_case': int((sheet.group == contrast['numerator']).sum()),
            'n_control': int((sheet.group == contrast['denominator']).sum()),
            'numerator': contrast['numerator'],
            'denominator': contrast['denominator'],
        }
    record = {
        'schema': ITEM_DEG_SCHEMA,
        'item': item,
        'category': 'genetic_cardiomyopathy',
        'deposition': {
            'matrix': {'path': spec['matrix'], 'sha256': sha256(matrix_path),
                       'bytes': int(matrix_path.stat().st_size),
                       'id_kind': spec.get('id_kind', 'symbol'),
                       'feature_count': int(frame.shape[0])},
            'series_matrix': {'path': spec['series_matrix'],
                              'sha256': sha256(series_path),
                              'bytes': int(series_path.stat().st_size)},
            'values': values_record,
        },
        'samples': {
            'path': _relative_to_item(item_dir, samples_path),
            'sha256': sha256(samples_path),
            'n': int(len(sheet)),
            'match_basis': spec['sample_match'],
            'group_source': spec['group_source'],
            'group_counts': group_counts,
            'covariates': sorted(list(spec.get('covariates', {}).values())
                                 + list(spec.get('derived_covariates', {}).keys())),
            'isogenic_control': spec['isogenic_control'],
        },
        'design': design,
        'contrasts': contrasts,
        'n_case': int(sheet.group.isin(case_arms).sum()),
        'n_control': int(sheet.group.isin(control_arms).sum()),
        'strata': spec['strata'],
        'method': {
            'deposition': ('integer counts' if kind == 'counts'
                           else 'normalised expression (log-scale as deposited)'),
            'model': ('limma-voom (voom library-size normalisation, moderated t)'
                      if kind == 'counts' else
                      'limma on the deposited log-scale matrix (moderated t)'),
            'implementation': 'stage2/util/output/api/stage2_util/derivation/disease_deg.R',
            'implementation_sha256': sha256(Path(disease_deg.DEG_WORKER)),
            'derivation': 'stage2/util/output/api/stage2_util/derivation/bulk_deg.py',
            'derivation_sha256': sha256(Path(__file__).resolve()),
            'test': 'moderated t, BH-adjusted p-value',
            'contrast_direction': spec['contrast_direction'],
            'limma_version': fit['limma_version'],
            'r_version': fit['r_version'],
        },
        'fit_summary': fit,
        'deg': {'path': _relative_to_item(item_dir, deg_path),
                'sha256': sha256(deg_path), 'rows': int(len(table)),
                'n_source_genes': int(table.source_gene_id.nunique()),
                'bytes': int(deg_path.stat().st_size)},
        'expr': {'path': _relative_to_item(item_dir, expr_path),
                 'sha256': sha256(expr_path)},
    }
    summary_path = local / f'{item}_deg.json'
    summary_path.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    return {'result': 'PASS', 'item': item, 'rows': int(len(table)),
            'contrasts': produced, 'n_case': record['n_case'],
            'n_control': record['n_control'],
            'deg_sha256': record['deg']['sha256'],
            'summary': _relative_to_item(item_dir, summary_path)}
