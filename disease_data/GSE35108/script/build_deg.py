#!/usr/bin/env python
"""GSE35108 — patient DCM iPSC-CM vs normal iPSC-CM (one-colour Affymetrix array).

The item-local derivation: the GEO series matrix (the deposited RMA log2
expression table plus the ``!Sample_*`` design block) and the GPL6244 platform
annotation (transcript cluster -> Entrez gene symbol) are bound here, the
probes are collapsed to one probe per gene, and the one live DEG implementation
``bases.disease_deg.R`` (limma on the deposited log2 matrix) is run through
``bases.disease_deg.run_deg_worker``. Nothing below that line re-implements the
fit.

Verified design (from the deposit, not from the ticket prose):

* 17 libraries total, of which **8 are cardiomyocytes**: a 2x2 design of
  genotype (Normal / DCM iPSC) x adenoviral Serca2a gene therapy
  (control / Serca2a), 2 libraries per cell. The remaining 9 libraries are
  fibroblast, H7 hESC and undifferentiated iPSC lines and are not part of the
  cardiomyocyte contrast.
* GPL6244 = [HuGene-1_0-st] Affymetrix Human Gene 1.0 ST array, one-colour
  (``!Sample_channel_count == 1``); the deposit is RMA-normalised log2
  expression (``!Sample_data_processing`` says "RMA as normalization method").
* Patient-derived DCM iPSC lines vs normal iPSC lines: a clinical,
  non-genetic-model comparison (no mutation was introduced, no isogenic
  control is deposited). The lines are different individuals, so genotype and
  line identity are confounded — recorded, not corrected.

Declared contrasts (both are DCM-vs-normal within one treatment stratum; the
strata are not averaged, they collapse per gene by maximum effect under the
ticket-01 multi-contrast rule)::

    DCM_control-Normal_control   (2 vs 2, untreated)
    DCM_Serca2a-Normal_Serca2a   (2 vs 2, after Serca2a gene therapy)

Usage (arc python):
    /home/chang/miniconda3/envs/arc/bin/python script/build_deg.py
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ITEM = Path(__file__).resolve().parents[1]
REPO = ITEM
while not (REPO / '.git').exists():
    if REPO.parent == REPO:
        raise SystemExit('cannot locate the repository root')
    REPO = REPO.parent
sys.path.insert(0, str(REPO / 'stage2/util/output/api'))
from stage2_util.derivation import deg_rank, disease_deg  # noqa: E402

SERIES_MATRIX = 'source/GSE35108_series_matrix.txt.gz'
PLATFORM_ANNOT = 'source/GPL6244.annot.gz'
PROBE_COLUMN = 'ID_REF'
PADJ_THRESHOLD = 0.05
DEG_LOG2FC_THRESHOLD = 2.0

SPEC = {
    'item': 'GSE35108',
    'category': 'nHF',
    'tier': 2,
    'evidence_type': 'bulk',
    'series_matrix': SERIES_MATRIX,
    'platform_annot': PLATFORM_ANNOT,
    'platform': 'GPL6244',
    'platform_title': ('[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array '
                       '[transcript (gene) version]'),
    'expected_sample_count': 17,
    'expected_platform_id': 'GPL6244',
    'cm_cell_type_prefix': 'Cardiomyocytes',
    'sample_match': 'geo_accession',
    'group_source': 'characteristics (cell type, genotype/variation, gene therapy)',
    'group_map': {
        ('Normal', 'control'): 'Normal_control',
        ('Normal', 'Serca2a'): 'Normal_Serca2a',
        ('DCM', 'control'): 'DCM_control',
        ('DCM', 'Serca2a'): 'DCM_Serca2a',
    },
    'kind': 'expression',
    'formula': '~0 + group',
    'contrasts': [
        {'label': 'DCM_control-Normal_control',
         'numerator': 'DCM_control', 'denominator': 'Normal_control'},
        {'label': 'DCM_Serca2a-Normal_Serca2a',
         'numerator': 'DCM_Serca2a', 'denominator': 'Normal_Serca2a'},
    ],
    'contrast_direction': ('positive log2FC = up in the patient DCM iPSC-CM arm '
                           'vs the normal iPSC-CM arm within the same treatment stratum'),
    'strata': {
        'decision': ('the deposit is a 2x2 genotype x Serca2a gene-therapy design; both '
                     'genotype contrasts are declared and collapse per gene by maximum '
                     'effect (ticket-01 rule), not averaged into one main effect'),
        'pooled': True,
        'blocking_covariates': [],
        'why': ('the Serca2a treatment is applied symmetrically to both genotypes, so each '
                'treatment stratum carries a valid DCM-vs-normal contrast; the untreated '
                'stratum is the primary disease contrast and its DEG count is recorded as '
                'a sensitivity in the tracked summary'),
    },
    'non_genetic_model_basis': (
        'patient-derived DCM iPSC-CM lines compared with normal iPSC-CM lines; no mutation '
        'was introduced and no isogenic control is deposited, so this is the clinical '
        '(non-mutation-modelled) iPSC-CM form that ADR 0071 admits to nHF tier 2'),
    'sensitivity': {
        'contrast': 'DCM_control-Normal_control',
        'why': 'the untreated 2 vs 2 stratum is the primary disease contrast',
    },
    'probe_collapse': {
        'rule': ('one probe per gene symbol: the probe with the maximum mean expression '
                 'across the 8 selected CM libraries; probes with no gene symbol or a '
                 'multi-symbol annotation (containing "///") are not used'),
        'symbol_column': 'Gene symbol',
        'id_column': 'ID',
    },
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_series_matrix(path):
    """Return ``(series, sample_fields, characteristic_rows, matrix)`` from a series matrix."""
    series = {}
    fields = {}
    characteristics = []
    with gzip.open(path, 'rt', encoding='utf-8-sig') as handle:
        for line in handle:
            line = line.rstrip('\r\n')
            if line.startswith('!series_matrix_table_begin'):
                break
            if line.startswith('!Series_'):
                key, _, value = line.partition('\t')
                series.setdefault(key[len('!Series_'):].lower(), []).append(
                    value.strip().strip('"'))
                continue
            if line.startswith('!Sample_characteristics_ch1\t'):
                characteristics.append([value.strip().strip('"')
                                        for value in line.split('\t')[1:]])
                continue
            match = re.match(r'^!Sample_([A-Za-z0-9_]+)\t(.+)$', line)
            if match:
                fields.setdefault(match.group(1).lower(), []).append(
                    [value.strip().strip('"') for value in match.group(2).split('\t')])
    with gzip.open(path, 'rt', encoding='utf-8-sig') as handle:
        for line in handle:
            if line.startswith('!series_matrix_table_begin'):
                break
        matrix = pd.read_csv(handle, sep='\t', index_col=0, dtype=str, comment='!')
    # The probe identifiers are numeric-looking strings; keep them strings so the
    # platform-annotation join is exact (no int/str coercion).
    matrix.index = matrix.index.astype(str)
    matrix.index.name = 'ID_REF'
    matrix = matrix.apply(pd.to_numeric, errors='raise')
    return series, fields, characteristics, matrix


def read_platform_annotation(path, spec):
    """Return the platform table's probe id -> gene symbol mapping."""
    with gzip.open(path, 'rt', encoding='utf-8-sig') as handle:
        for line in handle:
            if line.startswith('!platform_table_begin'):
                break
        table = pd.read_csv(io.StringIO(handle.read()), sep='\t', dtype=str,
                            usecols=[spec['probe_collapse']['id_column'],
                                     spec['probe_collapse']['symbol_column']])
    ids = table[spec['probe_collapse']['id_column']].astype(str).str.strip()
    symbols = table[spec['probe_collapse']['symbol_column']]
    return dict(zip(ids, symbols)), int(len(table))


def single_symbol(value):
    """A probe's gene symbol when it is unambiguous; ``''`` for missing/multi-mapped."""
    if not isinstance(value, str):
        return ''
    value = value.strip()
    if not value or '///' in value:
        return ''
    return value


def build_sample_sheet(fields, characteristics, spec):
    """Bind the deposited design block to one fit design for the 8 CM libraries."""
    titles = fields['title'][0]
    accessions = fields['geo_accession'][0]
    if len(titles) != spec['expected_sample_count'] or len(accessions) != len(titles):
        raise ValueError(f'expected {spec["expected_sample_count"]} samples, '
                         f'found titles={len(titles)} accessions={len(accessions)}')
    cell_type = [value.split(':', 1)[1].strip() for value in characteristics[0]]
    genotype = [value.split(':', 1)[1].strip() for value in characteristics[1]]
    therapy = [value.split(':', 1)[1].strip() if value else ''
               for value in characteristics[2]]
    selected = [index for index, value in enumerate(cell_type)
                if value.startswith(spec['cm_cell_type_prefix'])]
    if not selected:
        raise ValueError('no cardiomyocyte library selected')
    rows = []
    unmapped = set()
    for index in selected:
        key = (genotype[index], therapy[index])
        arm = spec['group_map'].get(key)
        if arm is None:
            unmapped.add(key)
            continue
        rows.append({'sample': accessions[index], 'group': arm,
                     'title': titles[index], 'cell_type': cell_type[index],
                     'genotype': genotype[index], 'gene_therapy': therapy[index]})
    if unmapped:
        raise ValueError(f'unmapped genotype/therapy value(s): {sorted(unmapped)}')
    sheet = pd.DataFrame(rows)
    if sheet['sample'].duplicated().any():
        raise ValueError('a GSM is selected twice')
    return sheet, {'total_samples': len(titles),
                   'cm_libraries': int(len(selected))}


def collapse_probes(matrix, symbols, sample_columns, spec):
    """Collapse transcript clusters to one probe per gene symbol (max mean expression)."""
    probe_symbol = pd.Series(
        {str(probe): single_symbol(symbols.get(str(probe), '')) for probe in matrix.index},
        dtype=object)
    mapped = probe_symbol[probe_symbol != '']
    means = matrix.loc[:, sample_columns].mean(axis=1)
    selected = []
    for _symbol, probes in mapped.groupby(mapped).groups.items():
        selected.append(means.loc[probes].astype(float).idxmax())
    selected = sorted(selected)
    stats = {
        'rule': spec['probe_collapse']['rule'],
        'n_probes': int(len(matrix)),
        'n_probes_with_single_symbol': int(len(mapped)),
        'n_probes_multi_or_missing_symbol': int(len(matrix) - len(mapped)),
        'n_genes': int(len(selected)),
        'n_probes_dropped_by_collapse': int(len(mapped) - len(selected)),
    }
    if not selected:
        raise ValueError('probe collapse selected no gene')
    return selected, probe_symbol, stats


def write_tsv_gz(frame, path):
    """Write a reproducible TSV gzip: fixed mtime and the declared column order."""
    payload = frame.to_csv(sep='\t', index=False, lineterminator='\n').encode()
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode='wb', mtime=0) as handle:
        handle.write(payload)
    Path(path).write_bytes(buffer.getvalue())


def deg_counts(table):
    """v2 DEG counts per contrast and per gene (max-effect collapse over contrasts)."""
    per_contrast = {}
    for label, rows in table.groupby('contrast'):
        flag = (rows.padj < PADJ_THRESHOLD) & (rows.log2fc.abs() >= DEG_LOG2FC_THRESHOLD)
        per_contrast[str(label)] = {
            'n_deg': int(flag.sum()),
            'n_up': int((flag & (rows.log2fc > 0)).sum()),
            'n_down': int((flag & (rows.log2fc < 0)).sum()),
            'n_padj_only': int((rows.padj < PADJ_THRESHOLD).sum()),
            'n_rows': int(len(rows)),
        }
    work = table.assign(
        effect=deg_rank.deg_effect(table.log2fc, table.padj),
        deg=(table.padj < PADJ_THRESHOLD) & (table.log2fc.abs() >= DEG_LOG2FC_THRESHOLD))
    # Multi-contrast collapse (ticket-01 rule): a gene's winning contrast is its
    # maximum-effect one, and its direction is read from that contrast.
    winner = work.sort_values(['gene_symbol', 'effect'], ascending=[True, False],
                              kind='stable').drop_duplicates('gene_symbol', keep='first')
    collapsed = winner[winner.deg]
    return {
        'definition': (f'padj < {PADJ_THRESHOLD} AND |log2FC| >= '
                       f'{DEG_LOG2FC_THRESHOLD} in the contrast (ADR 0071); a gene is a '
                       f'DEG when any declared contrast passes, collapsible per gene'),
        'per_contrast': per_contrast,
        'collapsed': {
            'n_genes': int(len(winner)),
            'n_deg': int(len(collapsed)),
            'n_up': int(collapsed.log2fc.gt(0).sum()),
            'n_down': int(collapsed.log2fc.lt(0).sum()),
            'n_padj_only_genes': int(table.loc[table.padj < PADJ_THRESHOLD,
                                               'gene_symbol'].nunique()),
        },
    }


def fit(expr, sheet, contrasts, spec, local, tag):
    """Write the fit inputs, run the one DEG worker, and read its table back."""
    expr_path = local / f'{tag}_expr.tsv.gz'
    samples_path = local / f'{tag}_samples.csv'
    write_tsv_gz(expr, expr_path)
    sheet.to_csv(samples_path, index=False)
    design = {'kind': spec['kind'], 'formula': spec['formula'], 'factor': 'group',
              'contrasts': contrasts}
    deg_path = local / f'{tag}_deg.csv'
    summary = disease_deg.run_deg_worker(expr_path, samples_path, design, deg_path)
    table = pd.read_csv(deg_path)
    missing = set(disease_deg.LOCAL_DEG_COLUMNS) - set(table.columns)
    if missing:
        raise ValueError(f'{tag}: DEG worker output lacks {sorted(missing)}')
    produced = sorted(table.contrast.astype(str).unique().tolist())
    declared = sorted(contrast['label'] for contrast in contrasts)
    if produced != declared:
        raise ValueError(f'{tag}: worker produced {produced}, declared {declared}')
    return expr_path, samples_path, deg_path, summary, table


def main():
    for relative in (SPEC['series_matrix'], SPEC['platform_annot']):
        if not (ITEM / relative).is_file():
            raise SystemExit(f'run script/fetch.sh first ({relative} missing)')
    item = SPEC['item']
    local = ITEM / 'local' / 'deg'
    local.mkdir(parents=True, exist_ok=True)
    series_path = ITEM / SPEC['series_matrix']
    annot_path = ITEM / SPEC['platform_annot']
    series, fields, characteristics, matrix = read_series_matrix(series_path)
    platform_ids = {value for row in fields.get('platform_id', []) for value in row}
    if platform_ids != {SPEC['expected_platform_id']}:
        raise ValueError(f'platform differs from the verified GPL6244: {sorted(platform_ids)}')
    channel = {value for row in fields.get('channel_count', []) for value in row}
    if channel != {'1'}:
        raise ValueError(f'expected a one-colour array, channel_count={sorted(channel)}')
    processing = ' '.join(value for row in fields.get('data_processing', [])
                          for value in row)
    if 'RMA' not in processing:
        raise ValueError('the deposit is not RMA-normalised as verified in the GEO record')
    sheet, design_counts = build_sample_sheet(fields, characteristics, SPEC)
    group_counts = {str(key): int(value)
                    for key, value in sheet.group.value_counts().sort_index().items()}
    for label, expected in (('Normal_control', 2), ('Normal_Serca2a', 2),
                            ('DCM_control', 2), ('DCM_Serca2a', 2)):
        if group_counts.get(label, 0) != expected:
            raise ValueError(f'arm {label}: {group_counts.get(label, 0)} libraries, '
                             f'expected {expected}')
    sample_columns = sheet['sample'].tolist()
    missing = set(sample_columns) - set(matrix.columns)
    if missing:
        raise ValueError(f'selected samples absent from the deposited matrix: {sorted(missing)}')
    symbols, annotation_rows = read_platform_annotation(annot_path, SPEC)
    selected, _probe_symbol, collapse_stats = collapse_probes(
        matrix, symbols, sample_columns, SPEC)
    expression = pd.DataFrame({
        'source_gene_id': [str(probe) for probe in selected],
        'gene_symbol': [single_symbol(symbols.get(str(probe), '')) for probe in selected],
    })
    for column in sample_columns:
        expression[column] = matrix.loc[selected, column].to_numpy(dtype=np.float64)
    values = matrix.loc[selected, sample_columns].to_numpy(dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError('the deposited RMA matrix carries non-finite values')
    if values.max() < 10:
        raise ValueError('the deposited matrix does not look like log2 expression')
    values_record = {
        'mode': 'normalised-expression-as-deposited (RMA log2, one-colour array)',
        'n_values': int(values.size),
        'min': float(values.min()),
        'median': float(np.median(values)),
        'max': float(values.max()),
        'n_non_positive': int((values <= 0).sum()),
        'one_colour': True,
        'normalisation': 'RMA (Affymetrix Expression Console, per !Sample_data_processing)',
    }
    _expr_path, samples_path, deg_path, fit_summary, table = fit(
        expression, sheet, SPEC['contrasts'], SPEC, local, item)
    counts = deg_counts(table)
    sensitivity = SPEC['sensitivity']
    sensitivity_contrast = [contrast for contrast in SPEC['contrasts']
                            if contrast['label'] == sensitivity['contrast']][0]
    sensitivity_sheet = sheet[sheet.group.isin([sensitivity_contrast['numerator'],
                                                sensitivity_contrast['denominator']]
                                               )].reset_index(drop=True)
    sensitivity_expr = expression[['source_gene_id', 'gene_symbol']
                                  + sensitivity_sheet['sample'].tolist()]
    with tempfile.TemporaryDirectory() as scratch:
        _e, _s, _d, _fit, sensitivity_table = fit(
            sensitivity_expr, sensitivity_sheet, [sensitivity_contrast], SPEC,
            Path(scratch), f'{item}_untreated')
        sensitivity_counts = deg_counts(sensitivity_table)
    contrasts = {}
    for contrast in SPEC['contrasts']:
        contrasts[contrast['label']] = {
            'n_case': int((sheet.group == contrast['numerator']).sum()),
            'n_control': int((sheet.group == contrast['denominator']).sum()),
            'numerator': contrast['numerator'],
            'denominator': contrast['denominator'],
        }
    case_arms = {contrast['numerator'] for contrast in SPEC['contrasts']}
    control_arms = {contrast['denominator'] for contrast in SPEC['contrasts']}
    if case_arms & control_arms:
        raise ValueError('an arm is both case and control')
    record = {
        'schema': 'stage2-item-deg-summary-v1',
        'item': item,
        'category': SPEC['category'],
        'tier': SPEC['tier'],
        'evidence_type': SPEC['evidence_type'],
        'non_genetic_model_basis': SPEC['non_genetic_model_basis'],
        'series': {
            'title': series.get('title', [''])[0],
            'platform': SPEC['platform'],
            'platform_title': SPEC['platform_title'],
            'assay': 'one-colour Affymetrix Gene 1.0 ST array (RMA log2 expression)',
            'total_libraries': design_counts['total_samples'],
            'cm_libraries': design_counts['cm_libraries'],
            'excluded_libraries': (design_counts['total_samples']
                                   - design_counts['cm_libraries']),
            'exclusion_reason': ('the deposit also profiles fibroblast, H7 hESC and '
                                 'undifferentiated iPSC lines; only the cardiomyocyte '
                                 'libraries carry the disease contrast'),
            'patient_vs_control': ('patient-derived familial-DCM iPSC-CM lines vs normal '
                                   'iPSC-CM lines; not isogenic and no causal variant is '
                                   'declared in the deposit, so genotype and line identity '
                                   'are confounded (recorded, not corrected)'),
        },
        'deposition': {
            'matrix': {'path': SPEC['series_matrix'], 'sha256': sha256(series_path),
                       'bytes': int(series_path.stat().st_size),
                       'probe_count': int(matrix.shape[0])},
            'platform_annot': {'path': SPEC['platform_annot'],
                               'sha256': sha256(annot_path),
                               'bytes': int(annot_path.stat().st_size),
                               'annotation_rows': int(len(symbols))},
            'values': values_record,
        },
        'probe_collapse': collapse_stats,
        'samples': {
            'path': str(samples_path.relative_to(ITEM)),
            'sha256': sha256(samples_path),
            'n': int(len(sheet)),
            'match_basis': SPEC['sample_match'],
            'group_source': SPEC['group_source'],
            'group_counts': group_counts,
        },
        'design': {'kind': SPEC['kind'], 'formula': SPEC['formula'], 'factor': 'group',
                   'contrasts': SPEC['contrasts']},
        'contrasts': contrasts,
        'n_case': int(sheet.group.isin(case_arms).sum()),
        'n_control': int(sheet.group.isin(control_arms).sum()),
        'strata': SPEC['strata'],
        'deg_counts': counts,
        'sensitivity': {
            'contrast': sensitivity['contrast'],
            'why': sensitivity['why'],
            'n_case': 2,
            'n_control': 2,
            'n_deg': sensitivity_counts['per_contrast'][sensitivity['contrast']]['n_deg'],
            'n_padj_only': sensitivity_counts['per_contrast'][
                sensitivity['contrast']]['n_padj_only'],
            'note': ('an item-local sensitivity fit of the untreated 2 vs 2 stratum only; '
                     'it is not a declared contrast of the published dataset'),
        },
        'method': {
            'deposition': 'deposited RMA log2 expression (one-colour array)',
            'collapse': collapse_stats['rule'],
            'model': 'limma on the deposited log2 matrix (moderated t)',
            'implementation': 'stage2/util/output/api/stage2_util/derivation/disease_deg.R',
            'implementation_sha256': sha256(Path(disease_deg.DEG_WORKER)),
            'derivation': str(Path(__file__).resolve().relative_to(REPO)),
            'derivation_sha256': sha256(Path(__file__).resolve()),
            'test': 'moderated t, BH-adjusted p-value',
            'contrast_direction': SPEC['contrast_direction'],
            'limma_version': fit_summary['limma_version'],
            'r_version': fit_summary['r_version'],
        },
        'fit_summary': fit_summary,
        'deg': {'path': str(deg_path.relative_to(ITEM)), 'sha256': sha256(deg_path),
                'rows': int(len(table)),
                'n_source_genes': int(table.gene_symbol.nunique()),
                'bytes': int(deg_path.stat().st_size)},
        'expr': {'path': str(_expr_path.relative_to(ITEM)), 'sha256': sha256(_expr_path)},
    }
    summary_path = local / f'{item}_deg.json'
    summary_path.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS', 'item': item,
                      'n_libraries': record['samples']['n'],
                      'n_deg': counts['collapsed']['n_deg'],
                      'deg_sha256': record['deg']['sha256'],
                      'summary': str(summary_path.relative_to(ITEM))}, indent=2,
                     sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
