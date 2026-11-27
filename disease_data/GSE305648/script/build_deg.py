#!/usr/bin/env python
"""GSE305648 — RBM20 DCM hiPSC cardioids: KO / PL / RQ against isogenic WT (WTC11).

Four genotype arms on the WTC11 background (wildtype, RBM20 KO, RBM20 PL, RBM20
RQ), 3 cardioid libraries each, with the three library batches rotated across the
arms (each arm sees A/B/C once). One joint fit `~0 + group + batch` blocks the
batch and yields the three mutant-vs-WT contrasts; the published dataset
collapses them per gene by maximum effect (ticket-01 multi-contrast rule).

Usage (arc python):
    /home/chang/miniconda3/envs/arc/bin/python script/build_deg.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ITEM = Path(__file__).resolve().parents[1]
REPO = ITEM
while not (REPO / '.git').exists():
    if REPO.parent == REPO:
        raise SystemExit('cannot locate the repository root')
    REPO = REPO.parent
sys.path.insert(0, str(REPO / 'cpst/util/output/api'))
from cpst_util.derivation import bulk_deg  # noqa: E402

SPEC = {
    'matrix': 'source/GSE305648_counts.tsv.gz',
    'sep': '\t',
    'series_matrix': 'source/GSE305648_series_matrix.txt.gz',
    'id_kind': 'ensembl',
    'sample_match': 'title',
    'group_source': 'genotype',
    'group_map': {'wildtype': 'wt', 'RBM20_KO': 'ko', 'RBM20_PL': 'pl', 'RBM20_RQ': 'rq'},
    'covariates': {'batch': 'batch'},
    'kind': 'counts',
    'formula': '~0 + group + batch',
    'contrasts': [
        {'label': 'RBM20_KO-WT', 'numerator': 'ko', 'denominator': 'wt'},
        {'label': 'RBM20_PL-WT', 'numerator': 'pl', 'denominator': 'wt'},
        {'label': 'RBM20_RQ-WT', 'numerator': 'rq', 'denominator': 'wt'},
    ],
    'contrast_direction': ('positive log2FC = up in the RBM20 mutant arm (KO / PL / RQ) '
                           'vs the isogenic wildtype cardioid'),
    'isogenic_control': {
        'arm': 'RBM20 wildtype cardioids, 3 libraries',
        'basis': ('all four arms are the WTC11 iPSC background; KO / PL / RQ are engineered '
                  'RBM20 alleles and the wildtype cardioids are the isogenic control'),
        'deposit_evidence': ('series-matrix !Sample_characteristics_ch1 "genotype: '
                             'RBM20_KO/RBM20_PL/RBM20_RQ/wildtype" (3 each)'),
    },
    'strata': {
        'decision': 'three mutant genotype arms (KO, PL, RQ) vs one isogenic wildtype arm',
        'pooled': True,
        'blocking_covariates': ['batch'],
        'why': ('the three mutant alleles are distinct RBM20 lesions on one WTC11 background, '
                'so each keeps its own contrast inside one joint fit and the three are '
                'collapsed per gene by maximum effect; the library batch (A/B/C) is rotated '
                'across the arms and enters the fit as a blocking covariate'),
    },
}


def main():
    if not (ITEM / SPEC['matrix']).is_file():
        raise SystemExit('run script/fetch.sh first')
    result = bulk_deg.derive_item_deg(ITEM, SPEC)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == '__main__':
    sys.exit(main())
