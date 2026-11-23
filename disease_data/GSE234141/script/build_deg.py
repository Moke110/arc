#!/usr/bin/env python
"""GSE234141 — TTN truncating-variant DCM hiPSC-CM, 2D monolayer and EHT, atrial- and
ventricular-like.

The 24 libraries span genotype (WT/TTNtv isogenic pair) x cell type
(atrial-like / ventricular-like) x format (2D / EHT). One joint fit over all 24
yields **four paired TTNtv-vs-WT contrasts**, one per format x cell-type stratum;
the published dataset collapses the four per gene by maximum effect (ticket-01
multi-contrast rule). The strata are pooled into one dataset rather than
published as four because they are two culture formats and two differentiation
outcomes of the same isogenic pair, not independent cohorts, and the deposit
itself is 3 libraries per cell (no stratum can be treated as its own study).

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

_TTN = 'WT/TTNtv (p.Ser21134Trpfs*25)'
_TWO_D = 'iPSC-derived cardiomyocytes (2D/monolayer)'
_EHT = 'iPSC-derived engineered heart tissue (3D)'


def _group(format_name, cell_type, genotype):
    return '|'.join((format_name, cell_type, genotype))


SPEC = {
    'matrix': 'source/GSE234141_counts.txt.gz',
    'sep': '\t',
    'series_matrix': 'source/GSE234141_series_matrix.txt.gz',
    'id_kind': 'ensembl',
    'sample_match': 'description',
    'group_source': ['tissue', 'cell_type', 'genotype'],
    'group_map': {
        _group(_TWO_D, 'ventricular-like', _TTN): '2d_ventricular_ttn',
        _group(_TWO_D, 'ventricular-like', 'WT/WT'): '2d_ventricular_wt',
        _group(_TWO_D, 'atrial-like', _TTN): '2d_atrial_ttn',
        _group(_TWO_D, 'atrial-like', 'WT/WT'): '2d_atrial_wt',
        _group(_EHT, 'ventricular-like', _TTN): 'eht_ventricular_ttn',
        _group(_EHT, 'ventricular-like', 'WT/WT'): 'eht_ventricular_wt',
        _group(_EHT, 'atrial-like', _TTN): 'eht_atrial_ttn',
        _group(_EHT, 'atrial-like', 'WT/WT'): 'eht_atrial_wt',
    },
    'kind': 'counts',
    'formula': '~0 + group',
    'contrasts': [
        {'label': 'TTNtv-WT.2d_ventricular', 'numerator': '2d_ventricular_ttn',
         'denominator': '2d_ventricular_wt'},
        {'label': 'TTNtv-WT.2d_atrial', 'numerator': '2d_atrial_ttn',
         'denominator': '2d_atrial_wt'},
        {'label': 'TTNtv-WT.eht_ventricular', 'numerator': 'eht_ventricular_ttn',
         'denominator': 'eht_ventricular_wt'},
        {'label': 'TTNtv-WT.eht_atrial', 'numerator': 'eht_atrial_ttn',
         'denominator': 'eht_atrial_wt'},
    ],
    'contrast_direction': ('positive log2FC = up in the TTN truncating-variant arm vs the '
                           'CRISPR-corrected isogenic WT/WT arm of the same format and cell type'),
    'isogenic_control': {
        'arm': ('WT/WT isogenic line, 12 libraries (3 per format x cell-type stratum)'),
        'basis': ('the deposit pairs a patient TTNtv iPSC line with its CRISPR/Cas9-corrected '
                  'isogenic WT/WT line; the WT/WT libraries are the control in every stratum'),
        'deposit_evidence': ('series-matrix !Sample_characteristics_ch1 "genotype: WT/TTNtv '
                             '(p.Ser21134Trpfs*25)" and "genotype: WT/WT"'),
    },
    'strata': {
        'decision': ('four paired strata (2D vs EHT) x (ventricular-like vs atrial-like), '
                     'each TTNtv vs WT, collapsed per gene by maximum effect'),
        'pooled': True,
        'blocking_covariates': [],
        'why': ('format and cell type are culture/differentiation strata of one isogenic '
                'pair, not independent cohorts; each stratum keeps its own paired contrast '
                'inside one joint fit (shared dispersion) and the four are collapsed per '
                'gene by maximum effect, the ticket-01 rule for a multi-arm source'),
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
