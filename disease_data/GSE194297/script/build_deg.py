#!/usr/bin/env python
"""GSE194297 — DMD point-mutation cardiac organoids vs the CRISPR-corrected isogenic line.

6 day-56 cardiac organoid libraries: 3 from the DMD patient-derived hiPSC line
carrying the point mutation and 3 from the same line after CRISPR/Cas9 correction.
One `~0 + group` fit yields the declared DMD-vs-corrected contrast; there is no
stratum to pool.

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
sys.path.insert(0, str(REPO / 'stage2_CarPS-VAE/util/output/api'))
from stage2_util.derivation import bulk_deg  # noqa: E402

SPEC = {
    'matrix': 'source/GSE194297_DMD_vs_DMD-Iso.absolute_counts.csv.gz',
    'sep': ';',
    'series_matrix': 'source/GSE194297_series_matrix.txt.gz',
    'id_kind': 'symbol',
    'sample_match': 'title_bracket',
    'group_source': 'genotype_variation',
    'group_map': {
        'Contains specific DMD point mutation': 'dmd',
        'DMD point mutation was corrected via CRISPR/Cas technology': 'corrected',
    },
    'kind': 'counts',
    'formula': '~0 + group',
    'contrasts': [
        {'label': 'DMD-DMD_corrected', 'numerator': 'dmd', 'denominator': 'corrected'},
    ],
    'contrast_direction': ('positive log2FC = up in the DMD point-mutation organoids vs the '
                           'CRISPR/Cas-corrected isogenic line'),
    'isogenic_control': {
        'arm': 'CRISPR/Cas9-corrected isogenic line, 3 day-56 cardiac organoids',
        'basis': ('the DMD point mutation was corrected in the same patient-derived hiPSC '
                  'line, so the corrected organoids are the isogenic control'),
        'deposit_evidence': ('series-matrix !Sample_characteristics_ch1 "genotype/variation: '
                             'Contains specific DMD point mutation" (3) and "DMD point '
                             'mutation was corrected via CRISPR/Cas technology" (3)'),
    },
    'strata': {
        'decision': 'single case/control pair (DMD vs corrected); no stratum to pool',
        'pooled': False,
        'blocking_covariates': [],
        'why': 'one mutant line against its corrected isogenic line at one time point (day 56)',
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
