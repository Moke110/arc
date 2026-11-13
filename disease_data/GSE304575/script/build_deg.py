#!/usr/bin/env python
"""GSE304575 — LMNA p.His222Pro DCM hiPSC-CM vs CRISPR-corrected isogenic control.

The cleanest isogenic deposit of the six: 6 LMNA H222P heterozygous libraries
against 6 CRISPR/Cas9-corrected control libraries, a single batch value
(`6/15/2023`) on every sample, so one `~0 + group` fit yields the declared
H222P-vs-corrected contrast with no stratum decision to make.

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
sys.path.insert(0, str(REPO / 'stage2/01_bases/execution'))
import public_api  # noqa: F401  (binds the single public implementation)
from bases import bulk_deg  # noqa: E402

SPEC = {
    'matrix': 'source/GSE304575_Raw_Counts_table.csv.gz',
    'sep': ';',
    'series_matrix': 'source/GSE304575_series_matrix.txt.gz',
    'id_kind': 'symbol',
    'sample_match': 'description',
    'group_source': 'genotype',
    'group_map': {'Corrected Isogenic control': 'corrected', 'LMNAH222P': 'lmna_h222p'},
    'kind': 'counts',
    'formula': '~0 + group',
    'contrasts': [
        {'label': 'LMNA_H222P-corrected', 'numerator': 'lmna_h222p',
         'denominator': 'corrected'},
    ],
    'contrast_direction': ('positive log2FC = up in the LMNA p.His222Pro heterozygous '
                           'arm vs the CRISPR-corrected isogenic control'),
    'isogenic_control': {
        'arm': 'Corrected Isogenic control, 6 libraries',
        'basis': ('patient hiPSC-CM carrying LMNA c.665A>C p.His222Pro heterozygous, '
                  'corrected with CRISPR/Cas9 to the isogenic control line'),
        'deposit_evidence': ('series-matrix !Sample_characteristics_ch1 "genotype: Corrected '
                             'Isogenic control" (6) and "genotype: LMNAH222P" (6); '
                             '"batch: 6/15/2023" on all 12'),
    },
    'strata': {
        'decision': 'single case/control pair; no stratum to pool',
        'pooled': False,
        'blocking_covariates': [],
        'why': ('one genotype arm against one isogenic control arm; the only batch value '
                'is constant across all 12 libraries, so no blocking term is estimable '
                'or needed'),
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
