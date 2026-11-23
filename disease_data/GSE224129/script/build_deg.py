#!/usr/bin/env python
"""GSE224129 — cMyBP-C (MYBPC3) HCM engineered cardiac tissue on one isogenic background.

The item-local derivation: the deposit's RSEM gene matrix (non-integer
``expected_count`` values on a count scale) and the GEO series matrix are bound
by the shared adapter ``bases.bulk_deg``, which integerises the count-scale matrix
with the shared rounding helper and runs the one live DEG implementation
``stage2_util/derivation/disease_deg.R`` (limma-voom). One joint fit over all 12 ECT libraries
yields two contrasts — MYBPC3 +/- and -/- against the isogenic +/+ arm — which
become one published dataset through the ticket-01 per-gene maximum-effect
collapse.

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
    'matrix': 'source/GSE224129_AHNNF3DMXX_rsem_GENE.txt.gz',
    'sep': '\t',
    'series_matrix': 'source/GSE224129_series_matrix.txt.gz',
    'id_kind': 'ensembl',
    'sample_match': 'title_bracket',
    'group_source': 'genotype',
    'group_map': {'cMyBP-C+/+': 'wt', 'cMyBP-C+/-': 'het', 'cMyBP-C-/-': 'hom'},
    'kind': 'counts',
    'integerisation': 'round_half_to_even',
    'formula': '~0 + group',
    'contrasts': [
        {'label': 'MYBPC3_het-MYBPC3_wt', 'numerator': 'het', 'denominator': 'wt'},
        {'label': 'MYBPC3_hom-MYBPC3_wt', 'numerator': 'hom', 'denominator': 'wt'},
    ],
    'contrast_direction': ('positive log2FC = up in the MYBPC3 mutant arm (+/- or -/-) '
                           'vs the isogenic cMyBP-C+/+ control'),
    'isogenic_control': {
        'arm': 'cMyBP-C+/+ ECT, 4 libraries',
        'basis': ('all three arms are the same DF19-9-11T.H iPSC background; CRISPR-Cas9 '
                  'introduced heterozygous and homozygous MYBPC3 frame-shift alleles, and '
                  'the +/+ arm is the isogenic control'),
        'deposit_evidence': ('series-matrix !Sample_characteristics_ch1 "genotype: '
                             'cMyBP-C+/+" (4 samples), "+/-" (4), "-/-" (4)'),
    },
    'strata': {
        'decision': 'two mutant genotype arms (heterozygous, homozygous) vs one isogenic +/+ control arm',
        'pooled': True,
        'blocking_covariates': [],
        'why': ('+/- and -/- are two doses of the same MYBPC3 lesion on one isogenic '
                'background: both contrasts come from one joint fit and collapse per gene '
                'by maximum effect (ticket-01 rule), not as two independent cohorts'),
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
