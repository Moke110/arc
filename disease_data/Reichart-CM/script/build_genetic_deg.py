#!/usr/bin/env python
"""Reichart-CM — the ``genetic_cardiomyopathy`` contrasts (ticket #188, Option A).

Ticket #187 fitted one joint design over the exported donors but published only
the ``nHF`` contrast (genotype-negative DCM vs normal). Under ADR 0066's
genotype-first split the genotype-positive DCM donors and the 8 ARVC donors are
``genetic_cardiomyopathy`` datasets, so this script re-runs the same single live
DEG implementation (``bases/disease_deg.R``, limma-voom) over the **tracked
export #187 produced** and publishes the two genetic contrasts:

    local/deg/reichart-cm_genetic_deg.csv    long table (gitignored, regenerable)
    local/deg/reichart-cm_genetic_deg.json   tracked summary (design/method/digests)

Inputs are consumed read-only and **not re-derived** here:

    local/pseudobulk/reichart-cm_cm_pseudobulk.csv.gz  (tracked per-donor CM pseudobulk)
    local/genotype/reichart-cm_donor_genotypes.csv     (tracked per-donor genotype split)

The 78 exported donors (18 normal / 8 genotype-negative DCM / 44
genotype-positive DCM / 8 ARVC; the single NCCM donor is excluded from the
export) are fitted in one design ``~0 + group`` so the contrasts share
dispersion; only the genotype-positive-DCM vs normal and ARVC vs normal
contrasts are selected into this item-local table. The genotype-negative
contrast stays #187's and is not touched.

Usage (arc python, from the item directory):
    /home/chang/miniconda3/envs/arc/bin/python script/build_genetic_deg.py
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
sys.path.insert(0, str(REPO / 'stage2/01_bases/execution'))
import public_api  # noqa: F401  (binds the single public implementation)
from bases import disease_deg  # noqa: E402

LOCAL = ITEM / 'local'
PSEUDOBULK = LOCAL / 'pseudobulk' / 'reichart-cm_cm_pseudobulk.csv.gz'
GENOTYPES = LOCAL / 'genotype' / 'reichart-cm_donor_genotypes.csv'
DEG = LOCAL / 'deg' / 'reichart-cm_genetic_deg.csv'
SAMPLES = LOCAL / 'deg' / 'reichart-cm_genetic_samples.csv'
MATRIX = LOCAL / 'deg' / 'reichart-cm_genetic_matrix.tsv.gz'
SUMMARY = LOCAL / 'deg' / 'reichart-cm_genetic_deg.json'
SCHEMA = 'stage2-item-deg-summary-v1'
CATEGORY = 'genetic_cardiomyopathy'
CONTRASTS = [
    {'label': 'dcm_genotype_positive-normal', 'numerator': 'dcm_genotype_positive',
     'denominator': 'normal'},
    {'label': 'arvc-normal', 'numerator': 'arvc', 'denominator': 'normal'},
]
DESIGN = {'kind': 'counts', 'formula': '~0 + group', 'factor': 'group',
          'contrasts': CONTRASTS}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    if not PSEUDOBULK.is_file() or not GENOTYPES.is_file():
        raise SystemExit('the tracked #187 pseudobulk/genotype export is missing')
    genotypes = pd.read_csv(GENOTYPES)
    if genotypes.donor_id.duplicated().any():
        raise SystemExit('genotype table repeats a donor')
    expected_groups = {'normal', 'dcm_genotype_negative', 'dcm_genotype_positive', 'arvc'}
    if set(genotypes.fit_group) != expected_groups:
        raise SystemExit(f'unexpected exported genotype groups: {sorted(set(genotypes.fit_group))}')
    excluded = genotypes.loc[genotypes.nHF_role == 'excluded', 'donor_id'].tolist()
    if excluded:
        raise SystemExit(f'export must not contain NCCM/excluded donors: {excluded}')
    samples = pd.DataFrame({'sample': genotypes.donor_id,
                            'group': genotypes.fit_group,
                            'donor': genotypes.donor_id})
    LOCAL.joinpath('deg').mkdir(parents=True, exist_ok=True)
    samples.to_csv(SAMPLES, index=False)
    pseudobulk = pd.read_csv(PSEUDOBULK)
    matrix = pd.DataFrame({'source_gene_id': pseudobulk.gene_symbol,
                           'gene_symbol': pseudobulk.gene_symbol})
    for donor in genotypes.donor_id:
        matrix[donor] = pseudobulk[donor]
    matrix.to_csv(MATRIX, index=False, sep='\t', compression={'method': 'gzip', 'mtime': 0})
    fit = disease_deg.run_deg_worker(MATRIX, SAMPLES, DESIGN, DEG)
    table = pd.read_csv(DEG)
    missing = set(disease_deg.LOCAL_DEG_COLUMNS) - set(table.columns)
    if missing:
        raise SystemExit(f'R output lacks {sorted(missing)}')
    produced = sorted(table.contrast.astype(str).unique().tolist())
    declared = sorted(contrast['label'] for contrast in CONTRASTS)
    if produced != declared:
        raise SystemExit(f'fit produced {produced}, expected {declared}')
    contrasts = {}
    for contrast in CONTRASTS:
        contrasts[contrast['label']] = {
            'n_case': int((samples.group == contrast['numerator']).sum()),
            'n_control': int((samples.group == contrast['denominator']).sum()),
            'numerator': contrast['numerator'],
            'denominator': contrast['denominator'],
        }
    record = {
        'schema': SCHEMA,
        'item': 'Reichart-CM',
        'category': CATEGORY,
        'gate': 'cell_type == "cardiac muscle cell" (CM-only cellxgene object)',
        'inputs': {
            'pseudobulk': {'path': 'local/pseudobulk/reichart-cm_cm_pseudobulk.csv.gz',
                           'sha256': sha256(PSEUDOBULK),
                           'note': ('tracked #187 export consumed read-only; the pseudobulk '
                                    'is NOT re-derived here')},
            'genotype': {'path': 'local/genotype/reichart-cm_donor_genotypes.csv',
                         'sha256': sha256(GENOTYPES),
                         'note': ('tracked #187 genotype split; genotype-first per '
                                  'Primary.Genetic.Diagnosis, NCCM donor excluded from the '
                                  'export')},
        },
        'design': DESIGN,
        'design_ambiguities': {
            'genotype_first': ('donors are split by Primary.Genetic.Diagnosis, not the study '
                               'disease label; the genotype-positive DCM and ARVC donors are '
                               'the genetic_cardiomyopathy datasets'),
            'joint_design': ('one joint fit ~0 + group over the 78 exported donors '
                             '(18 normal / 8 genotype-negative DCM / 44 genotype-positive DCM '
                             '/ 8 ARVC); this item-local table publishes the two genetic '
                             'contrasts only and leaves #187\'s genotype-negative contrast '
                             'and its tracked summary untouched'),
            'shared_controls': ('both genetic contrasts use the same 18 normal donors, which '
                                'are also nHF\'s Reichart controls; the normal signal is '
                                'therefore counted in both categories (declared, not '
                                'corrected)'),
        },
        'method': {
            'model': 'limma-voom (voom library-size normalisation, moderated t)',
            'implementation': 'stage2/01_bases/output/api/bases/disease_deg.R',
            'implementation_sha256': sha256(Path(disease_deg.DEG_WORKER)),
            'derivation': 'disease_data/Reichart-CM/script/build_genetic_deg.py',
            'derivation_sha256': sha256(Path(__file__).resolve()),
            'test': 'moderated t, BH-adjusted p-value',
            'contrast_direction': ('positive log2FC = up in the genotype-positive DCM (or '
                                   'ARVC) donors vs the object\'s normal donors'),
            'limma_version': fit['limma_version'],
            'r_version': fit['r_version'],
        },
        'contrasts': contrasts,
        'n_case': int(genotypes.fit_group.isin(
            [c['numerator'] for c in CONTRASTS]).sum()),
        'n_control': int((genotypes.fit_group == 'normal').sum()),
        'fit_summary': fit,
        'deg': {'path': 'local/deg/reichart-cm_genetic_deg.csv',
                'sha256': sha256(DEG), 'rows': int(len(table)),
                'n_source_genes': int(table.source_gene_id.nunique())},
    }
    SUMMARY.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS', 'rows': int(len(table)),
                      'contrasts': produced, 'n_case': record['n_case'],
                      'n_control': record['n_control'],
                      'deg_sha256': record['deg']['sha256']}, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
