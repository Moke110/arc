#!/usr/bin/env python
"""GSE253223 — acquired PKP2 loss in wild-type hiPSC-CM (siPKP2 vs siNeg, days 2/4/6/8).

This is a **knockdown, not a patient variant**: wild-type iCell Cardiomyocytes2
were transfected with PKP2 siRNA (siPKP2, 4 libraries per harvest day) or a
negative-control siRNA (siNeg, 1 library per day) and harvested on days 2, 4, 6
and 8. The four harvest days are pooled inside one joint fit with `day` as a
blocking factor — one siNeg library per day means a per-day contrast would rest
on a single control, and the day effect is a strong technical confounder — and
the fit yields the single declared siPKP2-vs-siNeg contrast.

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
sys.path.insert(0, str(REPO / 'stage2/util/output/api'))
from stage2_util.derivation import bulk_deg  # noqa: E402

SPEC = {
    'matrix': 'source/GSE253223_TNYA0016.counts.txt.gz',
    'sep': '\t',
    'series_matrix': 'source/GSE253223_series_matrix.txt.gz',
    'id_kind': 'symbol',
    'sample_match': 'title_bracket',
    'group_source': 'treatment',
    'group_map': {'siNeg': 'control', 'siPKP2': 'case'},
    'derived_covariates': {'day': {'from': 'title', 'regex': r'_D(\d+)'}},
    'kind': 'counts',
    'integerisation': 'round_half_to_even',
    'formula': '~0 + group + day',
    'contrasts': [
        {'label': 'siPKP2-siNeg', 'numerator': 'case', 'denominator': 'control'},
    ],
    'contrast_direction': ('positive log2FC = up in the PKP2-knockdown arm (siPKP2) vs the '
                           'negative-control siRNA (siNeg)'),
    'isogenic_control': {
        'arm': 'siNeg (non-targeting siRNA), 4 libraries (one per harvest day)',
        'basis': ('a single wild-type hiPSC-CM line (iCell Cardiomyocytes2) in both arms; '
                  'the negative-control siRNA is the matched control. This is an acquired '
                  'knockdown, not a patient variant and not an engineered isogenic line, '
                  'so the isogenic exception is applied to the model system, not to a genotype'),
        'deposit_evidence': ('series-matrix !Sample_characteristics_ch1 "treatment: siNeg" / '
                             '"treatment: siPKP2", "genotype: WT" on all 20 libraries'),
    },
    'strata': {
        'decision': 'four harvest days (2/4/6/8) pooled in one fit with day as a blocking covariate',
        'pooled': True,
        'blocking_covariates': ['day'],
        'why': ('each day contributes one siNeg library and four siPKP2 libraries, so a '
                'per-day case/control contrast would rest on a single control; pooling the '
                'days with day as a covariate uses all 20 libraries for one knockdown effect '
                'and removes the harvest-time effect'),
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
