#!/usr/bin/env python
"""Reichart-CM — one joint-design limma-voom fit; only the nHF contrast publishes.

The genotype-negative DCM donors are contrasted against the object's normal
donors. Genotype-positive DCM and ARVC donors stay in the same joint design
(``~ 0 + status``) so dispersion is estimated on all 78 exported donors, but
this ticket publishes no DEG for them: ticket #188 owns the
``genetic_cardiomyopathy`` contrast and consumes the tracked pseudobulk +
genotype export produced by ``build_pseudobulk.py``.

Inputs (declared ``local/`` additions; the pseudobulk/genotype pair is the
tracked compact export for ticket #188):

    local/pseudobulk/reichart-cm_cm_pseudobulk.csv.gz
    local/genotype/reichart-cm_donor_genotypes.csv

Outputs:

    local/deg/reichart-cm_deg.csv    long table, nHF contrast only (gitignored)
    local/deg/reichart-cm_deg.json   tracked summary (design/method/digests)

Usage (arc python):
    /home/chang/miniconda3/envs/arc/bin/python script/build_deg.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ITEM = Path(__file__).resolve().parents[1]
REPO = ITEM
while not (REPO / '.git').exists():
    if REPO.parent == REPO:
        raise SystemExit('cannot locate the repository root')
    REPO = REPO.parent
sys.path.insert(0, str(REPO / 'stage2_CarPS-VAE/util/output/api'))
from stage2_util.derivation import disease_deg  # noqa: E402

LOCAL = ITEM / 'local'
PSEUDOBULK = LOCAL / 'pseudobulk' / 'reichart-cm_cm_pseudobulk.csv.gz'
GENOTYPES = LOCAL / 'genotype' / 'reichart-cm_donor_genotypes.csv'
DEG = LOCAL / 'deg' / 'reichart-cm_deg.csv'
SAMPLES = LOCAL / 'deg' / 'reichart-cm_samples.csv'
SUMMARY = LOCAL / 'deg' / 'reichart-cm_deg.json'
SCHEMA = 'stage2-item-deg-summary-v1'
CONTRAST = 'dcm_genotype_negative-normal'
DESIGN = {
    'kind': 'counts',
    'formula': '~0 + group',
    'factor': 'group',
    'contrasts': [{'label': CONTRAST, 'numerator': 'dcm_genotype_negative',
                   'denominator': 'normal'}],
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    if not PSEUDOBULK.is_file() or not GENOTYPES.is_file():
        raise SystemExit('run script/build_pseudobulk.py first')
    genotypes = pd.read_csv(GENOTYPES)
    if genotypes.donor_id.duplicated().any():
        raise SystemExit('genotype table repeats a donor')
    samples = pd.DataFrame({'sample': genotypes.donor_id,
                            'group': genotypes.fit_group,
                            'donor': genotypes.donor_id})
    LOCAL.joinpath('deg').mkdir(parents=True, exist_ok=True)
    samples.to_csv(SAMPLES, index=False)
    pseudobulk = pd.read_csv(PSEUDOBULK)
    id_columns = ['gene_index', 'gene_symbol']
    matrix = pd.DataFrame({'source_gene_id': pseudobulk.gene_symbol,
                           'gene_symbol': pseudobulk.gene_symbol})
    for donor in genotypes.donor_id:
        matrix[donor] = pseudobulk[donor]
    matrix.to_csv(LOCAL / 'deg' / 'reichart-cm_pseudobulk_matrix.tsv.gz',
                  index=False, sep='\t', compression={'method': 'gzip', 'mtime': 0})
    matrix_path = LOCAL / 'deg' / 'reichart-cm_pseudobulk_matrix.tsv.gz'
    summary = disease_deg.run_deg_worker(matrix_path, SAMPLES, DESIGN, DEG)
    table = pd.read_csv(DEG)
    missing = set(disease_deg.LOCAL_DEG_COLUMNS) - set(table.columns)
    if missing:
        raise SystemExit(f'R output lacks {sorted(missing)}')
    if sorted(table.contrast.unique()) != [CONTRAST]:
        raise SystemExit('the Reichart fit published an unexpected contrast')
    n_case = int((samples.group == 'dcm_genotype_negative').sum())
    n_control = int((samples.group == 'normal').sum())
    record = {
        'schema': SCHEMA,
        'item': 'Reichart-CM',
        'category': 'nHF',
        'gate': 'cell_type == "cardiac muscle cell" (CM-only object)',
        'pseudobulk': {
            'path': 'local/pseudobulk/reichart-cm_cm_pseudobulk.csv.gz',
            'sha256': sha256(PSEUDOBULK),
            'dtype': 'integer counts summed per donor across that donor\'s CM cells',
            'axis_aligned': True,
            'donors': int(len(genotypes)),
            'excluded_donors': genotypes.loc[genotypes.nHF_role == 'excluded',
                                            'donor_id'].tolist(),
            'cells_per_donor': {row.donor_id: int(row.cells)
                                for row in genotypes.itertuples()},
        },
        'genotype_table': {'path': 'local/genotype/reichart-cm_donor_genotypes.csv',
                           'sha256': sha256(GENOTYPES),
                           'note': 'tracked compact export for ticket #188'},
        'design': DESIGN,
        'design_ambiguities': {
            'joint_design': ('one joint fit ~0 + status over normal / '
                             'genotype-negative DCM / genotype-positive DCM / ARVC; '
                             'only the genotype-negative DCM vs normal contrast is '
                             'published here'),
            'genotype_first': ('donors are split by Primary.Genetic.Diagnosis, not by '
                               'the study disease label; genotype-positive DCM and the '
                               '8 ARVC donors move out of nHF'),
        },
        'method': {
            'model': 'limma-voom (voom library-size normalisation, moderated t)',
            'implementation': 'stage2_CarPS-VAE/util/output/api/stage2_util/derivation/disease_deg.R',
            'implementation_sha256': sha256(Path(disease_deg.DEG_WORKER)),
            'test': 'moderated t, BH-adjusted p-value',
            'contrast_direction': ('positive log2FC = up in genotype-negative DCM vs '
                                   'normal'),
            'limma_version': summary['limma_version'],
            'r_version': summary['r_version'],
        },
        'contrasts': {CONTRAST: {'n_case': n_case, 'n_control': n_control}},
        'fit_summary': summary,
        'deg': {'path': 'local/deg/reichart-cm_deg.csv', 'sha256': sha256(DEG),
                'rows': int(len(table)),
                'n_source_genes': int(table.source_gene_id.nunique())},
    }
    SUMMARY.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS', 'rows': len(table),
                      'n_case': n_case, 'n_control': n_control,
                      'deg_sha256': record['deg']['sha256']}, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
