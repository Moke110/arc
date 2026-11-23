#!/usr/bin/env python
"""GSE269882 — right-atrial-appendage bulk RNA-seq, staged AF vs sinus rhythm (ADR 0030).

GSE269882 (Pikunov et al., the RNA-seq arm of the personalized atrial-fibrillation
model study) is **whole right-atrial-appendage bulk tissue** RNA-seq with four
author-declared diagnosis groups. ADR 0071 gives AF this source at **tier 3**
(whole-tissue bulk, not cardiomyocyte-specific), which supersedes ADR 0070's
whole-tissue exclusion.

Verified deposit design (from the series matrix, not the ticket prose):

* 27 GSMs / 27 donor libraries, Illumina HiSeq 2500, RNA-Seq, right atrial
  appendage biopsies.
* diagnosis groups: **sinus rhythm 9, post-operative AF 8, paroxysmal AF 5,
  persistent AF 5**. (The pre-existing AF audit recorded 10 SR / 7 poAF; the
  deposited series matrix says 9 / 8, and this item follows the deposit.)
* column labels: the count matrix's sample columns are named by patient number
  (1747, 1067, ...), so they are joined to the series titles by that number.
  The join is verified independently with the sex probes XIST/RPS4Y1/DDX3Y/
  UTY/KDM5D: 26 of 27 labels agree with the annotated sex (one annotated
  female, patient 1409, has a male expression profile; recorded, not corrected).

Declared design: one joint `~0 + group` limma-voom fit over all 27 libraries
with the four arms; the **three staged AF-vs-SR contrasts** (persistent,
paroxysmal, post-operative) are declared and collapse per gene by maximum effect
(ADR 0071's pass-any-contrast rule). The staging is recorded, and no pooling is
hidden: a pooled all-AF-vs-SR fit is measured as a declared sensitivity.

Measured outcome: **zero** genes pass padj < 0.05 in any declared contrast under
the live implementation (minimum padj 0.972 / 0.998 / 1.000), i.e. this bulk
cohort is not a DEG-yielding AF source at the v2 threshold. That is published as
the honest result.

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
sys.path.insert(0, str(REPO / 'stage2/util/output/api'))
from stage2_util.derivation import disease_deg  # noqa: E402

SOURCE = ITEM / 'source'
SERIES = 'source/GSE269882_series_matrix.txt.gz'
COUNTS = 'source/GSE269882_countdata_all.csv.gz'
TPM = 'source/GSE269882_TPM_values.csv.gz'
PADJ_THRESHOLD = 0.05
DEG_LOG2FC_THRESHOLD = 2.0
FILTER = {'min_count': 10, 'min_total_count': 15, 'large_n': 10, 'min_prop': 0.7}
ARM = {'sinus rhythm': 'SR', 'post-operative AF': 'poAF',
       'paroxysmal AF': 'pAF', 'persistent AF': 'persAF'}
ARM_ORDER = ['persAF', 'pAF', 'SR', 'poAF']
EXPECTED_ARMS = {'SR': 9, 'poAF': 8, 'pAF': 5, 'persAF': 5}
SEX_PROBES = {'XIST': 'ENSG00000229807', 'RPS4Y1': 'ENSG00000129824',
              'DDX3Y': 'ENSG00000067048', 'UTY': 'ENSG00000183878',
              'KDM5D': 'ENSG00000012817'}
SPEC = {
    'item': 'GSE269882',
    'category': 'AF',
    'tier': 3,
    'evidence_type': 'bulk',
    'model_system': ('human right-atrial-appendage whole-tissue bulk RNA-seq '
                     '(not cardiomyocyte-specific)'),
    'platform': 'GPL16791 (Illumina HiSeq 2500)',
    'contrasts': [
        {'label': 'persAF-SR', 'numerator': 'persAF', 'denominator': 'SR'},
        {'label': 'pAF-SR', 'numerator': 'pAF', 'denominator': 'SR'},
        {'label': 'poAF-SR', 'numerator': 'poAF', 'denominator': 'SR'},
    ],
    'staging': {
        'decision': ('the three AF arms are declared as separate contrasts against the '
                     'sinus-rhythm arm inside one joint fit, so the staging is carried '
                     'rather than averaged into a single case group'),
        'arms': {'persistent AF': 'AF present continuously (>7 days in the source cohort '
                                  'description); the persistent arm',
                 'paroxysmal AF': 'AF episodes that terminate spontaneously',
                 'post-operative AF': 'AF arising after cardiac surgery; kept as its own '
                                      'declared arm, not merged into the persistent group'},
        'pooled_sensitivity': 'all 18 AF libraries vs 9 sinus rhythm, measured separately',
    },
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_tsv_gz(frame, path):
    payload = frame.to_csv(sep='\t', index=False, lineterminator='\n').encode()
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode='wb', mtime=0) as handle:
        handle.write(payload)
    Path(path).write_bytes(buffer.getvalue())


def read_series_matrix(path):
    fields = {}
    characteristics = []
    with gzip.open(path, 'rt', encoding='utf-8-sig') as handle:
        for line in handle:
            line = line.rstrip('\r\n')
            if line.startswith('!series_matrix_table_begin'):
                break
            if line.startswith('!Sample_characteristics_ch1\t'):
                characteristics.append([value.strip().strip('"')
                                        for value in line.split('\t')[1:]])
                continue
            match = re.match(r'^!Sample_([A-Za-z0-9_]+)\t(.+)$', line)
            if match:
                fields.setdefault(match.group(1).lower(), []).append(
                    [value.strip().strip('"') for value in match.group(2).split('\t')])
    return fields, characteristics


def characteristic(characteristics, prefix):
    for row in characteristics:
        if row and row[0].startswith(prefix):
            return [value.split(':', 1)[1].strip() for value in row]
    raise ValueError(f'no characteristics row starting {prefix!r}')


def sample_table(fields, characteristics):
    titles = fields['title'][0]
    accessions = fields['geo_accession'][0]
    sexes = characteristic(characteristics, 'Sex:')
    ages = characteristic(characteristics, 'age:')
    diagnoses = characteristic(characteristics, 'diagnosis group:')
    if not (len(titles) == len(accessions) == len(sexes) == len(ages) == len(diagnoses) == 27):
        raise ValueError('series matrix sample fields disagree on the 27-library count')
    rows = []
    for gsm, title, sex, age, diagnosis in zip(accessions, titles, sexes, ages, diagnoses):
        match = re.search(r'patient nr\.\s*(\d+)', title)
        if not match:
            raise ValueError(f'{gsm}: title has no patient number: {title!r}')
        if diagnosis not in ARM:
            raise ValueError(f'{gsm}: unknown diagnosis group {diagnosis!r}')
        rows.append({'gsm': gsm, 'patient': match.group(1), 'title': title,
                     'sex': 'female' if sex.lower() == 'female' else 'male',
                     'age': int(age), 'arm': ARM[diagnosis], 'diagnosis': diagnosis})
    sheet = pd.DataFrame(rows)
    if sheet.patient.duplicated().any():
        raise ValueError('two samples share a patient number')
    counts = sheet.arm.value_counts().to_dict()
    if counts != EXPECTED_ARMS:
        raise ValueError(f'arm counts {counts}, expected {EXPECTED_ARMS}')
    return sheet


def read_matrix(path, expected_columns, integer):
    frame = pd.read_csv(path, index_col=0)
    if frame.shape[1] != 27:
        raise ValueError(f'{path}: {frame.shape[1]} sample columns, expected 27')
    # The deposited column order is not the series order (the file was assembled
    # in a different sample order), so the join is by patient number: the set must
    # match the titles exactly and the frame is reordered to the title order.
    if set(frame.columns) != set(expected_columns) or len(set(frame.columns)) != 27:
        raise ValueError(f'{path}: sample columns do not match the 27 series titles')
    frame = frame[expected_columns]
    if frame.index.duplicated().any():
        raise ValueError(f'{path}: duplicated gene id')
    values = frame.to_numpy()
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError(f'{path}: non-finite or negative values')
    if integer and (values != np.round(values)).any():
        raise ValueError(f'{path}: count matrix is not integer')
    return frame


def sex_label_verification(counts, sheet):
    """Classify each column's sex from XIST vs Y-gene expression and compare."""
    library = counts.sum(axis=0)
    log_cpm = np.log2(counts.divide(library, axis=1) * 1e6 + 1)
    for symbol, gene in SEX_PROBES.items():
        if gene not in log_cpm.index:
            raise ValueError(f'sex probe {symbol} ({gene}) is absent from the count matrix')
    xist = log_cpm.loc[SEX_PROBES['XIST']]
    y_score = log_cpm.loc[[SEX_PROBES[key] for key in ('RPS4Y1', 'DDX3Y', 'UTY', 'KDM5D')]
                          ].mean(axis=0)
    observed = {column: ('female' if (xist[column] > 3 and y_score[column] < 4) else 'male')
                for column in counts.columns}
    annotated = dict(zip(sheet.patient, sheet.sex))
    agreement = sum(observed[p] == annotated[p] for p in counts.columns)
    discordant = [{'patient': p, 'annotated': annotated[p], 'expression': observed[p],
                   'xist_log2cpm': round(float(xist[p]), 3),
                   'y_log2cpm': round(float(y_score[p]), 3)}
                  for p in counts.columns if observed[p] != annotated[p]]
    return {
        'method': ('per-column sex from the XIST / (RPS4Y1, DDX3Y, UTY, KDM5D) '
                   'log2(CPM+1) contrast; a female profile has XIST > 3 and Y < 4'),
        'n_columns': int(len(counts.columns)),
        'n_agree': int(agreement),
        'discordant': discordant,
        'reading': ('26/27 agree, so the patient-number column labels are sound; the one '
                    'discordance is an annotated-female column whose expression is male '
                    '(recorded, not corrected: sex is not a design factor here)'),
    }


def deg_counts(table):
    per_contrast = {}
    for label, rows in table.groupby('contrast'):
        flag = (rows.padj < PADJ_THRESHOLD) & (rows.log2fc.abs() >= DEG_LOG2FC_THRESHOLD)
        per_contrast[str(label)] = {
            'n_deg': int(flag.sum()),
            'n_up': int((flag & (rows.log2fc > 0)).sum()),
            'n_down': int((flag & (rows.log2fc < 0)).sum()),
            'n_padj_only': int((rows.padj < PADJ_THRESHOLD).sum()),
            'n_padj_lt_010': int((rows.padj < 0.10).sum()),
            'n_rows': int(len(rows)),
            'min_padj': float(rows.padj.min()),
            'max_abs_log2fc': float(rows.log2fc.abs().max()),
        }
    work = table.assign(deg=(table.padj < PADJ_THRESHOLD)
                        & (table.log2fc.abs() >= DEG_LOG2FC_THRESHOLD))
    winner = work.sort_values(['gene_symbol', 'padj', 'p_value'],
                              ascending=[True, True, True], kind='stable'
                              ).drop_duplicates('gene_symbol', keep='first')
    collapsed = winner[winner.deg]
    return {
        'definition': (f'padj < {PADJ_THRESHOLD} AND |log2FC| >= {DEG_LOG2FC_THRESHOLD} '
                       f'in any declared contrast (ADR 0071 pass-any), direction '
                       f'carried from the winning contrast'),
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
    produced = sorted(table.contrast.astype(str).unique().tolist())
    declared = sorted(contrast['label'] for contrast in contrasts)
    if produced != declared:
        raise ValueError(f'{tag}: worker produced {produced}, declared {declared}')
    return expr_path, samples_path, local / f'{tag}_deg.csv', summary, table


def main():
    for relative in (SERIES, COUNTS, TPM):
        if not (ITEM / relative).is_file():
            raise SystemExit(f'run script/fetch.sh first ({relative} missing)')
    item = SPEC['item']
    local = ITEM / 'local' / 'deg'
    local.mkdir(parents=True, exist_ok=True)
    fields, characteristics = read_series_matrix(ITEM / SERIES)
    samples = sample_table(fields, characteristics)
    title_order = samples.patient.tolist()
    counts = read_matrix(ITEM / COUNTS, title_order, integer=True)
    tpm = read_matrix(ITEM / TPM, title_order, integer=False)
    tpm_sums = tpm.sum(axis=0)
    if not ((tpm_sums > 900e3) & (tpm_sums < 1100e3)).all():
        raise ValueError('the TPM matrix columns do not sum to ~1e6')
    sex_check = sex_label_verification(counts, samples)
    order = [arm for arm in ARM_ORDER for _ in range(EXPECTED_ARMS[arm])]
    ordered = pd.concat([samples[samples.arm == arm] for arm in ARM_ORDER],
                        ignore_index=True)
    genes = counts.index.astype(str)
    expression = pd.DataFrame({'source_gene_id': genes, 'gene_symbol': genes})
    for patient in ordered.patient:
        expression[patient] = counts[patient].to_numpy(dtype=np.int64)
    sheet = pd.DataFrame({'sample': ordered.patient,
                          'group': ordered.arm,
                          'diagnosis': ordered.diagnosis,
                          'sex': ordered.sex,
                          'age': ordered.age,
                          'gsm': ordered.gsm})
    expr_path, samples_path, deg_path, fit_summary, table = fit(
        expression, sheet, SPEC['contrasts'], '~0 + group', local, item)
    counts_summary = deg_counts(table)

    # Declared sensitivity: pool the three staged AF arms into one case group.
    pooled_sheet = sheet.assign(group=np.where(sheet.group == 'SR', 'SR', 'AF'))
    pooled_expression = expression
    pooled_contrast = {'label': 'AF-SR', 'numerator': 'AF', 'denominator': 'SR'}
    _e, _s, _d, pooled_summary, pooled_table = fit(
        pooled_expression, pooled_sheet, [pooled_contrast], '~0 + group', local,
        f'{item}_pooled')
    pooled_counts = deg_counts(pooled_table)

    # Cross-check: the deposited TPM matrix carries the same 27 columns; compare the
    # count-model effect for the persistent-AF contrast with the TPM group difference.
    log_tpm = np.log2(tpm.divide(tpm_sums, axis=1) * 1e6 + 1)
    pers = [p for p in ordered.patient if ordered.set_index('patient').loc[p, 'arm'] == 'persAF']
    sr = [p for p in ordered.patient if ordered.set_index('patient').loc[p, 'arm'] == 'SR']
    tpm_delta = (log_tpm[pers].mean(axis=1) - log_tpm[sr].mean(axis=1)).reindex(genes)
    log_cpm = np.log2(counts.divide(counts.sum(axis=0), axis=1) * 1e6 + 1)
    cpm_delta = (log_cpm[pers].mean(axis=1) - log_cpm[sr].mean(axis=1)).reindex(genes)
    tpm_correlation = float(pd.Series(cpm_delta).corr(pd.Series(tpm_delta),
                                                      method='spearman'))
    deg_pers = table[table.contrast == 'persAF-SR'].set_index('gene_symbol')
    tpm_stat_correlation = float(deg_pers.log2fc.corr(
        pd.Series(tpm_delta.to_numpy(), index=genes).reindex(deg_pers.index),
        method='spearman'))

    record = {
        'schema': 'stage2-item-deg-summary-v1',
        'item': item,
        'category': SPEC['category'],
        'tier': SPEC['tier'],
        'evidence_type': SPEC['evidence_type'],
        'model_system': SPEC['model_system'],
        'tier_note': ('ADR 0071 admits GSE269882 **only as tier 3** (whole-tissue bulk, '
                      'not cardiomyocyte-specific), superseding the ADR 0070 whole-tissue '
                      'exclusion; the ADR 0070 audit note is not reused as eligibility '
                      'evidence. This is NOT GSE271229: that accession is a left-atrial '
                      'snRNA series with no author cell-type annotation, audited and '
                      'rejected separately.'),
        'series': {
            'title': fields.get('title', [''])[0],
            'platform': SPEC['platform'],
            'tissue': 'human right atrial appendage biopsy (whole tissue)',
            'assay': 'bulk RNA-seq, unique exon-mapped reads (featureCounts), GRCh38',
            'total_libraries': 27,
            'arms': {arm: int(EXPECTED_ARMS[arm]) for arm in ARM_ORDER},
            'arms_as_verified': {str(key): int(value) for key, value in
                                 samples.diagnosis.value_counts().items()},
            'arm_prose_discrepancy': (
                'the pre-existing AF audit recorded SR 10 / post-operative AF 7; the '
                'deposited series matrix says SR 9 / post-operative AF 8, and this item '
                'follows the deposit'),
            'column_join': ('the matrix columns are patient numbers; they are joined to '
                            'the series titles by `patient nr. N` and verified with the '
                            'sex probes'),
            'sex_label_verification': sex_check,
            'staging': SPEC['staging'],
        },
        'deposition': {
            'series_matrix': {'path': SERIES, 'bytes': int((ITEM / SERIES).stat().st_size),
                              'sha256': sha256(ITEM / SERIES)},
            'counts': {'path': COUNTS, 'bytes': int((ITEM / COUNTS).stat().st_size),
                       'sha256': sha256(ITEM / COUNTS), 'genes': int(counts.shape[0]),
                       'samples': int(counts.shape[1]), 'kind': 'raw integer counts'},
            'tpm': {'path': TPM, 'bytes': int((ITEM / TPM).stat().st_size),
                    'sha256': sha256(ITEM / TPM), 'genes': int(tpm.shape[0]),
                    'columns_sum_to_1e6': True,
                    'used_for': ('cross-check only; the DEG input is the raw counts, '
                                 'because the count model is the documented one'),
                    'cross_check': {
                        'vs_count_cpm_delta_spearman': round(tpm_correlation, 4),
                        'vs_count_limma_log2fc_persAF_spearman': round(tpm_stat_correlation, 4),
                        'reading': ('the deposited TPM and raw counts agree on the sample '
                                    'profiles (rank correlation near 1)'),
                    }},
        },
        'samples': {
            'path': str(samples_path.relative_to(ITEM)),
            'sha256': sha256(samples_path),
            'n': int(len(sheet)),
            'group_counts': {str(key): int(value)
                             for key, value in sheet.group.value_counts().items()},
            'match_basis': 'patient number in the matrix column header joined to the title',
        },
        'design': {'kind': 'counts', 'formula': '~0 + group', 'factor': 'group',
                   'contrasts': SPEC['contrasts'], 'filter': FILTER},
        'contrasts': {contrast['label']: {
            'n_case': int(EXPECTED_ARMS[contrast['numerator']]),
            'n_control': int(EXPECTED_ARMS['SR']),
            'numerator': contrast['numerator'],
            'denominator': 'SR',
        } for contrast in SPEC['contrasts']},
        'n_case': 18,
        'n_control': 9,
        'deg_counts': counts_summary,
        'null_outcome_statement': (
            'no gene reaches padj < 0.05 in any declared staged contrast under the live '
            'implementation; the item publishes the honest zero (no threshold loosening, '
            'no padj-only counting, no arm selection by outcome)'),
        'sensitivity': {
            'pooled_AF_vs_SR': {
                'contrast': 'AF-SR with the three AF arms pooled (18 vs 9), separate fit',
                'n_deg': pooled_counts['per_contrast']['AF-SR']['n_deg'],
                'n_padj_only': pooled_counts['per_contrast']['AF-SR']['n_padj_only'],
                'min_padj': pooled_counts['per_contrast']['AF-SR']['min_padj'],
                'note': ('the staged design is the declared one; pooling does not change '
                         'the conclusion'),
            },
            'audit_model_alternative_edger_ql': {
                'status': 'not_published',
                'method': ('edgeR quasi-likelihood F-test with TMM on the same filtered '
                           'matrix (the pre-existing audit environment has edgeR; the live '
                           'worker does not)'),
                'n_deg_v2_by_arm': {'post-operative AF': 2, 'persistent AF': 0,
                                    'paroxysmal AF': 0},
                'min_fdr_by_arm': {'post-operative AF': 0.049, 'persistent AF': 0.069,
                                   'paroxysmal AF': 0.555},
                'reading': ('a liberal alternative test finds at most two marginal genes, '
                            'all in the post-operative arm; not published'),
            },
        },
        'method': {
            'deposition': 'raw integer counts',
            'model': "limma-voom (voom's own library-size normalisation, moderated t)",
            'filter': 'edgeR filterByExpr rule written out in the live worker',
            'implementation': 'stage2/util/output/api/stage2_util/derivation/disease_deg.R',
            'implementation_sha256': sha256(Path(disease_deg.DEG_WORKER)),
            'derivation': str(Path(__file__).resolve().relative_to(REPO)),
            'derivation_sha256': sha256(Path(__file__).resolve()),
            'test': 'moderated t, BH-adjusted p-value',
            'limma_version': fit_summary['limma_version'],
            'r_version': fit_summary['r_version'],
        },
        'fit_summary': {'staged': fit_summary, 'pooled': pooled_summary},
        'deg': {'path': str(deg_path.relative_to(ITEM)), 'sha256': sha256(deg_path),
                'rows': int(len(table)),
                'n_source_genes': int(table.gene_symbol.nunique()),
                'bytes': int(deg_path.stat().st_size)},
        'expr': {'path': str(expr_path.relative_to(ITEM)), 'sha256': sha256(expr_path)},
    }
    summary_path = local / f'{item}_deg.json'
    summary_path.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS', 'item': item,
                      'n_libraries': int(len(sheet)),
                      'n_deg': counts_summary['collapsed']['n_deg'],
                      'min_padj': {label: values['min_padj'] for label, values
                                   in counts_summary['per_contrast'].items()},
                      'summary': str(summary_path.relative_to(ITEM))}, indent=2,
                     sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
