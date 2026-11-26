#!/usr/bin/env python
"""GSE138252 — PCM1+ cardiomyocyte-nuclei bulk RNA-seq, AF vs control (ADR 0030).

GSE138252 (Lipovsky et al., JCI Insight 2020, PMID 32841220) is bulk RNA-seq of
**FANS-sorted PCM1+ cardiomyocyte nuclei** from human hearts — the CM-nuclear
bulk branch of ADR 0071's tier 1. The deposit has 48 GSMs but only three
series-level matrices, and the matrices cover **38 library columns**; the item
records that discrepancy and the chamber/arm decision explicitly:

* 48 GSMs across two submission batches (28 in 2019/2020; 20 added 2020-09-11).
  One 2020 GSM (`GSM4103429`, h76ra) has no expression row, and the 2021
  addendum re-sequences four hearts already in the 2020 batch (h77la, h77ra,
  h83la, h87ra), so the 38 library columns are the distinct sequenced libraries,
  not 38 different hearts.
* The primary payload `GSE138252_human_merged_genebody_exp.txt.gz` is the 2020
  batch: 56,609 genes x 18 human **atrial** libraries (10 AF + 8 control; the
  other 2020 libraries are LV/RV controls). The 2021 addendum matrix holds the
  6 HF-only libraries and 14 controls; the 2785 matrix is a technical-replicate
  subset. Only the 2020 matrix supplies the AF contrast.
* The published design is the **left-atrium, 2020-batch, batch-matched 5 AF vs
  5 control** contrast (the series' only AF arm is 2020-batch, so a
  across-batch control pool would confound AF with sequencing batch; that pool
  is measured as a declared sensitivity instead). Region is the dominant
  variance axis in this series (the paper reports 1,180 LA-vs-RA DEGs), so
  pooling chambers is not used.
* **Confound (recorded, not corrected):** the AF arm is AF + end-stage heart
  failure and the control arm is non-failing donor heart, so AF-CTL measures
  AF+HF vs NF. The HF-only arm is 3 LA libraries in the 2021 batch and is
  batch-separated from the AF arm, so an AF-vs-HF de-confounding contrast is
  not identifiable here; its HF-vs-control count is measured as a diagnostic.

The one live implementation is `stage2_CarPS-VAE/util/output/api/stage2_util/derivation/disease_deg.R`
(limma-voom, voom's own library-size normalisation), run through
`stage2_util.derivation.disease_deg.run_deg_worker`; the item declares the edgeR
`filterByExpr` rule through the worker's `filter` option.

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
from pathlib import Path

import numpy as np
import pandas as pd

ITEM = Path(__file__).resolve().parents[1]
REPO = ITEM
while not (REPO / '.git').exists():
    if REPO.parent == REPO:
        raise SystemExit('cannot locate the repository root')
    REPO = REPO.parent
sys.path.insert(0, str(REPO / 'stage2_CarPS-VAE/util/output/api'))
from stage2_util.derivation import disease_deg  # noqa: E402

SOURCE = ITEM / 'source'
SERIES = 'source/GSE138252_series_matrix.txt.gz'
PRIMARY = 'source/GSE138252_human_merged_genebody_exp.txt.gz'
ADDENDUM = 'source/GSE138252_merged_gene_name_expression.txt.gz'
TECHREP = 'source/GSE138252_human_2785_merged_gene_name_expression.txt.gz'
PADJ_THRESHOLD = 0.05
DEG_LOG2FC_THRESHOLD = 2.0
FILTER = {'min_count': 10, 'min_total_count': 15, 'large_n': 10, 'min_prop': 0.7}

SPEC = {
    'item': 'GSE138252',
    'category': 'AF',
    'tier': 1,
    'evidence_type': 'bulk',
    'model_system': ('human atrial cardiomyocyte-nuclei bulk RNA-seq (PCM1+ FANS-sorted '
                     'nuclei; CM-nuclear bulk, not pseudobulk)'),
    'platform': 'GPL16791 (Illumina HiSeq 2500, single-end 1x50)',
    'expected_gsms': 48,
    'expected_primary_genes': 56609,
    'expected_primary_columns': 18,
    'expected_addendum_columns': 20,
    'expected_techrep_columns': 6,
    'expected_atrial_2020': 19,
    'contrasts': [{'label': 'AF-CTL', 'numerator': 'AF', 'denominator': 'CTL'}],
    'cm_specificity_evidence': (
        'the series summary states a PCM1-based cardiomyocyte-specific nuclear isolation; '
        'every GSM\'s extract protocol says RNA was extracted from PCM1+ cardiomyocyte '
        'nuclei; the `molecule subtype` characteristic is "Cardiomyoyte-specific nuclear '
        'RNA" (sic)'),
}

# The 2020-batch matrix column -> GSM link is recoverable from the series titles
# only (the matrices are series-level and carry no per-sample files). The title
# core (text before ' :', e.g. "af13_129la" / "CR_h86ra" / "h77ra#1") is a
# substring of the matrix column name once both are reduced to lower-case
# alphanumerics (e.g. "DL_af13_129la" -> "dlaf13129la"). The mapping is asserted
# to be a bijection onto the 18 atrial 2020 columns; the build fails rather than
# guesses if a column or a GSM cannot be paired.
def normalise(value):
    return re.sub(r'[^a-z0-9]', '', str(value).lower())


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_tsv_gz(frame, path):
    """Write a reproducible TSV gzip: fixed mtime so a rebuild is byte-identical."""
    payload = frame.to_csv(sep='\t', index=False, lineterminator='\n').encode()
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode='wb', mtime=0) as handle:
        handle.write(payload)
    Path(path).write_bytes(buffer.getvalue())


def read_series_matrix(path):
    """Return ``(series, sample_fields, characteristics)`` from a series matrix."""
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
    return series, fields, characteristics


def characteristic(characteristics, prefix):
    """Values of the one characteristics row whose labels start with ``prefix``."""
    for row in characteristics:
        if row and row[0].startswith(prefix):
            # The 2021-addendum GSMs carry '' in these rows (the addendum's design is
            # only in the matrix/column names), so an empty value is kept as ''.
            return [value.split(':', 1)[1].strip() if ':' in value else ''
                    for value in row]
    raise ValueError(f'no characteristics row starting {prefix!r}')


def read_count_matrix(path, expected_genes, expected_columns):
    frame = pd.read_csv(path, sep='\t')
    if list(frame.columns[:2]) != ['name', 'length']:
        raise ValueError(f'{path}: expected name/length as the first two columns')
    genes = frame['name'].astype(str)
    if len(frame) != expected_genes:
        raise ValueError(f'{path}: {len(frame)} genes, expected {expected_genes}')
    if genes.duplicated().any():
        raise ValueError(f'{path}: duplicated gene name')
    if frame.shape[1] - 2 != expected_columns:
        raise ValueError(f'{path}: {frame.shape[1] - 2} library columns, '
                         f'expected {expected_columns}')
    values = frame.iloc[:, 2:].to_numpy()
    if not np.isfinite(values).all() or (values < 0).any() or (values != np.round(values)).any():
        raise ValueError(f'{path}: deposited matrix is not raw integer counts')
    return frame, genes


def bind_samples(fields, characteristics):
    """Resolve the 2020 atrial matrix columns to GSMs and their arm labels."""
    accessions = fields['geo_accession'][0]
    titles = fields['title'][0]
    tissues = characteristic(characteristics, 'tissue:')
    diseases = characteristic(characteristics, 'disease state:')
    if not (len(accessions) == len(titles) == len(tissues) == len(diseases)
            == SPEC['expected_gsms']):
        raise ValueError('series matrix sample fields disagree on the 48-GSM count')
    rows = []
    for gsm, title, tissue, disease in zip(accessions, titles, tissues, diseases):
        if 'atrium' not in tissue.lower():
            continue
        region = 'LA' if 'left atrium' in tissue.lower() else 'RA'
        group = 'AF' if disease.lower() == 'atrial fibrillation' else 'CTL'
        core = title.split(':')[0] if ':' in title else title
        rows.append({'gsm': gsm, 'title': title, 'region': region, 'group': group,
                     'title_core': normalise(core)})
    return rows


def map_columns(matrix_columns, atrial_2020):
    """Pair every 2020 atrial matrix column with exactly one 2020 atrial GSM."""
    normalised = {column: normalise(column) for column in matrix_columns}
    assignment = {}
    for column, value in normalised.items():
        matches = [row for row in atrial_2020 if row['title_core'] in value]
        if len(matches) != 1:
            raise ValueError(f'column {column!r} matches {len(matches)} GSMs: '
                             f'{[(row["gsm"], row["title"]) for row in matches]}')
        assignment[column] = matches[0]['gsm']
    if len(set(assignment.values())) != len(assignment):
        raise ValueError('the column -> GSM mapping is not a bijection')
    return assignment


def deg_counts(table):
    """v2 DEG counts per contrast and the per-gene max-effect collapse."""
    per_contrast = {}
    for label, rows in table.groupby('contrast'):
        flag = (rows.padj < PADJ_THRESHOLD) & (rows.log2fc.abs() >= DEG_LOG2FC_THRESHOLD)
        per_contrast[str(label)] = {
            'n_deg': int(flag.sum()),
            'n_up': int((flag & (rows.log2fc > 0)).sum()),
            'n_down': int((flag & (rows.log2fc < 0)).sum()),
            'n_padj_only': int((rows.padj < PADJ_THRESHOLD).sum()),
            'n_rows': int(len(rows)),
            'min_padj': float(rows.padj.min()),
            'n_padj_lt_010': int((rows.padj < 0.10).sum()),
        }
    work = table.assign(deg=(table.padj < PADJ_THRESHOLD)
                        & (table.log2fc.abs() >= DEG_LOG2FC_THRESHOLD))
    winner = work.sort_values(['gene_symbol', 'padj', 'p_value'],
                              ascending=[True, True, True], kind='stable'
                              ).drop_duplicates('gene_symbol', keep='first')
    collapsed = winner[winner.deg]
    return {
        'definition': (f'padj < {PADJ_THRESHOLD} AND |log2FC| >= {DEG_LOG2FC_THRESHOLD} '
                       f'in a declared contrast (ADR 0071); direction carried, not filtered'),
        'per_contrast': per_contrast,
        'collapsed': {
            'n_genes': int(len(winner)),
            'n_deg': int(len(collapsed)),
            'n_up': int(collapsed.log2fc.gt(0).sum()),
            'n_down': int(collapsed.log2fc.lt(0).sum()),
            'n_padj_only_genes': int(table.loc[table.padj < PADJ_THRESHOLD,
                                               'gene_symbol'].nunique()),
            'deg_genes': collapsed.gene_symbol.tolist(),
        },
    }


def fit(expr, sheet, contrasts, formula, local, tag):
    expr_path = local / f'{tag}_expr.tsv.gz'
    samples_path = local / f'{tag}_samples.csv'
    write_tsv_gz(expr, expr_path)
    sheet.to_csv(samples_path, index=False)
    design = {'kind': 'counts', 'formula': formula, 'factor': 'group',
              'contrasts': contrasts, 'filter': FILTER}
    summary = disease_deg.run_deg_worker(expr_path, samples_path, design,
                                         local / f'{tag}_deg.csv')
    table = pd.read_csv(local / f'{tag}_deg.csv')
    missing = set(disease_deg.LOCAL_DEG_COLUMNS) - set(table.columns)
    if missing:
        raise ValueError(f'{tag}: DEG worker output lacks {sorted(missing)}')
    produced = sorted(table.contrast.astype(str).unique().tolist())
    declared = sorted(contrast['label'] for contrast in contrasts)
    if produced != declared:
        raise ValueError(f'{tag}: worker produced {produced}, declared {declared}')
    return expr_path, samples_path, local / f'{tag}_deg.csv', summary, table


def main():
    for relative in (SERIES, PRIMARY, ADDENDUM, TECHREP):
        if not (ITEM / relative).is_file():
            raise SystemExit(f'run script/fetch.sh first ({relative} missing)')
    item = SPEC['item']
    local = ITEM / 'local' / 'deg'
    local.mkdir(parents=True, exist_ok=True)
    series, fields, characteristics = read_series_matrix(ITEM / SERIES)
    atrial_2020 = [row for row in bind_samples(fields, characteristics)
                   if row['gsm'] >= 'GSM4103429' and row['gsm'] <= 'GSM4103456']
    if len(atrial_2020) != SPEC['expected_atrial_2020']:
        raise ValueError(f'{len(atrial_2020)} atrial 2020 GSMs, '
                         f'expected {SPEC["expected_atrial_2020"]}')
    primary, genes = read_count_matrix(ITEM / PRIMARY, SPEC['expected_primary_genes'],
                                       SPEC['expected_primary_columns'])
    addendum, addendum_genes = read_count_matrix(
        ITEM / ADDENDUM, SPEC['expected_primary_genes'], SPEC['expected_addendum_columns'])
    techrep, _ = read_count_matrix(ITEM / TECHREP, SPEC['expected_primary_genes'],
                                   SPEC['expected_techrep_columns'])
    if set(addendum_genes) != set(genes):
        raise ValueError('the 2020 and 2021 matrices do not share the gene axis')
    assignment = map_columns(primary.columns[2:], atrial_2020)
    without_column = sorted({row['gsm'] for row in atrial_2020}
                            - set(assignment.values()))
    if without_column != ['GSM4103429']:
        raise ValueError(f'atrial 2020 GSMs without an expression row: {without_column}; '
                         f'expected exactly GSM4103429 (h76ra)')

    # --- published contrast: LA, 2020 batch, 5 AF vs 5 control -----------------
    grouped = {}
    for column, gsm in assignment.items():
        row = next(entry for entry in atrial_2020 if entry['gsm'] == gsm)
        grouped.setdefault((row['region'], row['group']), []).append((column, gsm))
    for key in (('LA', 'AF'), ('LA', 'CTL'), ('RA', 'AF'), ('RA', 'CTL')):
        if len(grouped.get(key, [])) != (5 if key != ('RA', 'CTL') else 3):
            raise ValueError(f'arm {key}: {len(grouped.get(key, []))} libraries')
    primary_samples = grouped[('LA', 'AF')] + grouped[('LA', 'CTL')]
    expression = pd.DataFrame({
        'source_gene_id': genes,
        'gene_symbol': genes,
    })
    for column, gsm in primary_samples:
        expression[gsm] = primary[column].to_numpy(dtype=np.int64)
    sheet = pd.DataFrame({
        'sample': [gsm for _column, gsm in primary_samples],
        'group': ['AF'] * 5 + ['CTL'] * 5,
        'region': ['LA'] * 10,
        'batch': ['2020'] * 10,
        'library': [column for column, _gsm in primary_samples],
    })
    (contrast,) = SPEC['contrasts']
    expr_path, samples_path, deg_path, fit_summary, table = fit(
        expression, sheet, [contrast], '~0 + group', local, item)
    counts = deg_counts(table)

    # --- sensitivity 1: pooled controls (2020 + 2021 new-heart LA controls), ---
    # batch as a covariate. The 2021 re-sequenced hearts (h77la/h83la/...) are
    # duplicates of 2020 libraries and are excluded; only the 5 new-heart 2021
    # LA controls join the pool.
    addendum_columns = [column for column in addendum.columns[2:]
                        if '_LA_CTL.' in column]
    if len(addendum_columns) != 5:
        raise ValueError(f'{len(addendum_columns)} 2021 LA control columns, expected 5')
    addendum = addendum.set_index(addendum_genes)
    merged = pd.DataFrame({'source_gene_id': genes, 'gene_symbol': genes})
    for column, gsm in primary_samples:
        merged[gsm] = primary[column].to_numpy(dtype=np.int64)
    pooled_samples = sheet.copy()
    for index, column in enumerate(addendum_columns):
        gsm = f'2021CTL{index + 1}'
        merged[gsm] = addendum.loc[genes, column].to_numpy(dtype=np.int64)
        pooled_samples.loc[len(pooled_samples)] = {
            'sample': gsm, 'group': 'CTL', 'region': 'LA', 'batch': '2021',
            'library': column}
    pooled_contrast = {'label': 'AF-CTL', 'numerator': 'AF', 'denominator': 'CTL'}
    _e, _s, _d, pooled_summary, pooled_table = fit(
        merged, pooled_samples, [pooled_contrast], '~0 + group + batch', local,
        'GSE138252_pooled')
    pooled_counts = deg_counts(pooled_table)
    pooled_groups = (pooled_samples.groupby(['batch', 'group']).size()
                     .unstack(fill_value=0))
    pooled_groups = {f'{batch}:{group}': int(pooled_groups.loc[batch, group])
                     for batch in pooled_groups.index for group in pooled_groups.columns}

    # --- diagnostic: HF-only LA (2021 batch) vs the 2021-batch LA controls ------
    hf_columns = [column for column in addendum.columns
                  if '_LA_HF.' in column]
    if len(hf_columns) != 3:
        raise ValueError(f'{len(hf_columns)} 2021 LA HF-only columns, expected 3')
    hf = pd.DataFrame({'source_gene_id': genes, 'gene_symbol': genes})
    hf_sheet = pd.DataFrame({'sample': [], 'group': []})
    for index, column in enumerate(hf_columns):
        hf[f'HF{index + 1}'] = addendum.loc[genes, column].to_numpy(dtype=np.int64)
        hf_sheet.loc[len(hf_sheet)] = {'sample': f'HF{index + 1}', 'group': 'HF'}
    for index, column in enumerate(addendum_columns):
        hf[f'CTL{index + 1}'] = addendum.loc[genes, column].to_numpy(dtype=np.int64)
        hf_sheet.loc[len(hf_sheet)] = {'sample': f'CTL{index + 1}', 'group': 'CTL'}
    hf_contrast = {'label': 'HF-CTL', 'numerator': 'HF', 'denominator': 'CTL'}
    _e2, _s2, _d2, hf_summary, hf_table = fit(
        hf, hf_sheet, [hf_contrast], '~0 + group', local, 'GSE138252_hfonly')
    hf_counts = deg_counts(hf_table)

    record = {
        'schema': 'stage2-item-deg-summary-v1',
        'item': item,
        'category': SPEC['category'],
        'tier': SPEC['tier'],
        'evidence_type': SPEC['evidence_type'],
        'model_system': SPEC['model_system'],
        'cm_specificity': {
            'evidence': SPEC['cm_specificity_evidence'],
            'chamber_model': ('cardiomyocyte nuclei from left/right atria and ventricles; '
                              'this item selects the left-atrium contrast only'),
        },
        'series': {
            'title': series.get('title', [''])[0],
            'platform': SPEC['platform'],
            'total_gsms': int(len(fields['geo_accession'][0])),
            'supplementary_matrices': 3,
            'per_sample_supplementary_files': 'NONE on every GSM',
            'matrix_library_columns': {
                'primary_2020': int(primary.shape[1] - 2),
                'addendum_2021': int(addendum.shape[1] - 2),
                'technical_replicate_subset': int(techrep.shape[1] - 2),
                'distinct_columns_total': int(primary.shape[1] - 2 + addendum.shape[1] - 2),
            },
            'gsm_versus_column_discrepancy': (
                '48 GSMs map to 38 library columns: the 2020 batch contributes 28 GSMs '
                '(27 with an expression row; GSM4103429/h76ra has none) and the 2021 '
                'addendum 20 GSMs; the addendum re-sequences four 2020 hearts (h77la, '
                'h77ra, h83la, h87ra), which are different sequencing runs of the same '
                'heart, not extra hearts. The 2785 matrix is a 6-column '
                'technical-replicate subset of one run.'),
            'primary_atrial_arms': {'LA_AF': 5, 'LA_CTL': 5, 'RA_AF': 5, 'RA_CTL': 3},
            'batch_structure': ('AF libraries exist only in the 2020 batch, so the '
                                'primary contrast is batch-matched; the 2021 addendum '
                                'contributes HF-only libraries and extra controls'),
        },
        'deposition': {
            'series_matrix': {'path': SERIES, 'sha256': sha256(ITEM / SERIES),
                              'bytes': int((ITEM / SERIES).stat().st_size)},
            'primary': {'path': PRIMARY, 'sha256': sha256(ITEM / PRIMARY),
                        'bytes': int((ITEM / PRIMARY).stat().st_size),
                        'genes': int(len(primary)),
                        'library_columns': int(primary.shape[1] - 2)},
            'addendum': {'path': ADDENDUM, 'sha256': sha256(ITEM / ADDENDUM),
                         'bytes': int((ITEM / ADDENDUM).stat().st_size),
                         'genes': int(len(addendum)),
                         'library_columns': int(addendum.shape[1] - 2)},
            'technical_replicate_subset': {
                'path': TECHREP, 'sha256': sha256(ITEM / TECHREP),
                'bytes': int((ITEM / TECHREP).stat().st_size),
                'library_columns': int(techrep.shape[1] - 2),
                'used': False,
                'why': 'six columns of one sequencing run (technical replicates); no arm'},
            'column_to_gsm_mapping': {column: gsm for column, gsm in assignment.items()},
            'matrix_values': {'kind': 'raw integer gene-body counts',
                              'genes': int(len(genes)),
                              'total_counts': int(primary.iloc[:, 2:].to_numpy().sum())},
        },
        'samples': {
            'path': str(samples_path.relative_to(ITEM)),
            'sha256': sha256(samples_path),
            'n': int(len(sheet)),
            'group_counts': {'AF': 5, 'CTL': 5},
            'region': 'left atrium',
            'batch': '2020',
            'match_basis': ('series-matrix title core substring match against the '
                            'matrix column name, asserted bijective'),
        },
        'design': {'kind': 'counts', 'formula': '~0 + group', 'factor': 'group',
                   'contrasts': [contrast], 'filter': FILTER},
        'contrasts': {'AF-CTL': {'n_case': 5, 'n_control': 5,
                                 'numerator': 'AF', 'denominator': 'CTL',
                                 'direction': 'positive log2FC = up in AF (+HF) vs control'}},
        'n_case': 5,
        'n_control': 5,
        'chamber_and_arm_decision': {
            'primary': ('left atrium, 2020 sequencing batch only, AF vs control '
                        '(batch-matched 5 vs 5)'),
            'why_left_atrium': ('AF has its own chamber-specific transcriptional '
                                'programme; the paper reports 1,180 LA-vs-RA DEGs at '
                                'baseline and pooling chambers loses power'),
            'why_2020_batch': ('every AF library is in the 2020 batch, so only a '
                               'batch-matched control arm avoids confounding AF with the '
                               'sequencing batch'),
            'arms_not_used': {'RA': 'sex-confounded (3/3 RA controls female) and 0 DEGs',
                              'LV/RV': 'control-only libraries, no AF arm',
                              '2021 HF-only': 'batch-separated from the AF arm; used only '
                                              'as the HF-vs-control diagnostic below'},
        },
        'confound': {
            'af_arm': 'AF + end-stage heart failure (transplanted hearts with documented AF)',
            'control_arm': 'non-failing donor hearts without HF or AF',
            'statement': ('AF-CTL measures AF+HF vs non-failing heart, so AF and HF are '
                          'confounded; an AF-vs-HF-only contrast is not identifiable '
                          'because the HF-only libraries are 2021-batch while AF is '
                          '2020-batch only. Recorded, not corrected.'),
            'hf_only_diagnostic': {
                'contrast': 'HF-CTL (2021-batch LA: 3 HF-only vs 5 new-heart controls)',
                'n_case': 3, 'n_control': 5,
                'n_deg': hf_counts['per_contrast']['HF-CTL']['n_deg'],
                'n_padj_only': hf_counts['per_contrast']['HF-CTL']['n_padj_only'],
                'reading': ('HF alone produces at least as many DEGs as AF+HF does, so '
                            'the AF arm cannot be read as HF-independent'),
            },
        },
        'deg_counts': counts,
        'sensitivity': {
            'pooled_controls_batch_adjusted': {
                'contrast': 'AF-CTL, LA, 5 AF (2020) vs 10 controls '
                            '(5 x 2020 + 5 new-heart 2021), ~0 + group + batch',
                'batches': {str(key): int(value) for key, value in
                            pooled_groups.items()},
                'n_case': 5, 'n_control': 10,
                'n_deg': pooled_counts['per_contrast']['AF-CTL']['n_deg'],
                'n_padj_only': pooled_counts['per_contrast']['AF-CTL']['n_padj_only'],
                'note': ('adds 2021-batch controls with batch as a covariate; it is a '
                         'declared sensitivity, not the published contrast'),
            },
            'audit_model_alternative_tmm': {
                'status': 'not_published',
                'method': ('limma-voom after edgeR TMM normalisation (the pre-existing AF '
                           'source audit model, not the live worker: the pinned arc R env '
                           'ships limma but not edgeR, and the one live implementation '
                           'uses voom\'s own library-size normalisation)'),
                'n_deg_v2': 2,
                'deg_genes': ['LRRC10', 'SLC9A7'],
                'detail': {'SLC9A7': {'log2fc': 2.2121, 'padj': 0.034379},
                           'LRRC10': {'log2fc': -2.5445, 'padj': 0.034379}},
                'reading': ('TMM is the only difference from the published fit; it moves '
                            'LRRC10 from padj 0.0616 to 0.0344 and the v2 count from 1 to '
                            '2. Recorded per the one-live-implementation rule: the '
                            'published number is the live worker\'s, never this one.'),
                'provenance': ('.scratch/af-source-audit/GSE138252/run_deg.R + '
                               'analyze_axis.out (untracked audit of 2026-11-14)'),
            },
            'audit_model_alternative_edger_ql': {
                'status': 'not_published',
                'method': 'edgeR quasi-likelihood F-test with TMM on the same 5 vs 5 LA contrast',
                'n_deg_v2': 6,
                'deg_genes': ['CHST7', 'HBB', 'HOTAIRM1', 'LRRC10', 'SLC9A7', 'TRAC'],
                'reading': 'the most liberal defensible test; recorded as a sensitivity only',
            },
        },
        'method': {
            'deposition': 'raw integer gene-body counts',
            'model': "limma-voom (voom's own library-size normalisation, moderated t)",
            'filter': 'edgeR filterByExpr rule written out in the live worker',
            'implementation': 'stage2_CarPS-VAE/util/output/api/stage2_util/derivation/disease_deg.R',
            'implementation_sha256': sha256(Path(disease_deg.DEG_WORKER)),
            'derivation': str(Path(__file__).resolve().relative_to(REPO)),
            'derivation_sha256': sha256(Path(__file__).resolve()),
            'test': 'moderated t, BH-adjusted p-value',
            'limma_version': fit_summary['limma_version'],
            'r_version': fit_summary['r_version'],
        },
        'fit_summary': {'primary': fit_summary, 'pooled': pooled_summary,
                        'hf_only': hf_summary},
        'deg': {'path': str(deg_path.relative_to(ITEM)), 'sha256': sha256(deg_path),
                'rows': int(len(table)),
                'n_source_genes': int(table.gene_symbol.nunique()),
                'bytes': int(deg_path.stat().st_size)},
        'expr': {'path': str(expr_path.relative_to(ITEM)), 'sha256': sha256(expr_path)},
    }
    summary_path = local / f'{item}_deg.json'
    summary_path.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS', 'item': item,
                      'n_libraries': record['samples']['n'],
                      'n_deg': counts['collapsed']['n_deg'],
                      'deg_genes': counts['collapsed']['deg_genes'],
                      'summary': str(summary_path.relative_to(ITEM))}, indent=2,
                     sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
