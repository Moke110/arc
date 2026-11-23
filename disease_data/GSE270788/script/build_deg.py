#!/usr/bin/env python
"""GSE270788 — one joint-design limma-voom DEG fit for the three HF arms.

All three arms (NICM 7, ICM 5, AMI 5) are contrasted against the same six Donor
controls, so they are **not** fitted as three case-control tests: one joint
design ``~ 0 + arm`` is fitted once over the 23 donor pseudobulk profiles and the
three contrasts are read from that single fit (ADR 0066 / spec D9). The AMI arm
carries no per-sample time-since-infarct field, so its acute status is
`author_declared` and is recorded here, not re-derived.

Inputs (built by build_pseudobulk.py, declared ``local/`` additions):

    local/pseudobulk/gse270788_cm_pseudobulk.csv.gz
    local/pseudobulk/gse270788_donors.csv

Outputs:

    local/deg/gse270788_deg.csv    long per-contrast table (gitignored, regenerable)
    local/deg/gse270788_deg.json   tracked summary (design/method/digests)

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
sys.path.insert(0, str(REPO / 'stage2/util/output/api'))
from stage2_util.derivation import disease_deg  # noqa: E402

LOCAL = ITEM / 'local'
PSEUDOBULK = LOCAL / 'pseudobulk' / 'gse270788_cm_pseudobulk.csv.gz'
DONORS = LOCAL / 'pseudobulk' / 'gse270788_donors.csv'
DEG = LOCAL / 'deg' / 'gse270788_deg.csv'
SAMPLES = LOCAL / 'deg' / 'gse270788_samples.csv'
SUMMARY = LOCAL / 'deg' / 'gse270788_deg.json'
SCHEMA = 'stage2-item-deg-summary-v1'
CONTRASTS = [
    {'label': 'NICM-Donor', 'numerator': 'NICM', 'denominator': 'Donor'},
    {'label': 'ICM-Donor', 'numerator': 'ICM', 'denominator': 'Donor'},
    {'label': 'AMI-Donor', 'numerator': 'AMI', 'denominator': 'Donor'},
]
DESIGN = {'kind': 'counts', 'formula': '~0 + group', 'factor': 'group',
          'contrasts': CONTRASTS}
AMBIGUITY = {
    'design': 'one joint fit ~ 0 + arm over all 23 donor pseudobulks',
    'why_joint': ('the NICM/ICM/AMI arms all use the same six Donor controls; '
                  'three separate case-control fits would count the Donor signal '
                  'three times and are not independent'),
    'donor_term': ('no donor covariate: each donor belongs to exactly one arm, so a '
                   'donor effect is nested in arm and not identifiable; the shared '
                   'control structure is handled by the single joint fit'),
    'AMI_time_window_basis': 'author_declared (no per-sample time-since-infarct field)',
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    if not PSEUDOBULK.is_file() or not DONORS.is_file():
        raise SystemExit('run script/build_pseudobulk.py first')
    donors = pd.read_csv(DONORS)
    samples = donors[['sample', 'arm']].rename(columns={'arm': 'group'})
    samples = samples.assign(donor=donors['sample'])
    LOCAL.joinpath('deg').mkdir(parents=True, exist_ok=True)
    samples.to_csv(SAMPLES, index=False)
    summary = disease_deg.run_deg_worker(PSEUDOBULK, SAMPLES, DESIGN, DEG)
    table = pd.read_csv(DEG)
    missing = set(disease_deg.LOCAL_DEG_COLUMNS) - set(table.columns)
    if missing:
        raise SystemExit(f'R output lacks {sorted(missing)}')
    counts = {}
    for contrast in CONTRASTS:
        counts[contrast['label']] = {
            'n_case': int((samples.group == contrast['numerator']).sum()),
            'n_control': int((samples.group == contrast['denominator']).sum()),
        }
    record = {
        'schema': SCHEMA,
        'item': 'GSE270788',
        'category': 'nHF (NICM) / ICM (ICM, AMI)',
        'gate': 'cell.type.l1 == "Cardiomyocyte"',
        'pseudobulk': {
            'path': 'local/pseudobulk/gse270788_cm_pseudobulk.csv.gz',
            'sha256': sha256(PSEUDOBULK),
            'dtype': 'integer counts summed per donor across that donor\'s CM cells',
            'cells_per_donor': {row.sample: int(row.cells_cm)
                                for row in donors.itertuples()},
        },
        'design': DESIGN,
        'design_ambiguities': {
            **AMBIGUITY,
            'arm_labels': ('the deposited labels are Donor 7 / NICM 6 / ICM 5 / AMI 5; '
                           'the survey/ticket statement Donor 6 / NICM 7 is not what the '
                           'immutable metadata or series matrix says, so the deposited '
                           'labels are used unchanged and the discrepancy is recorded, '
                           'not relabelled'),
        },
        'method': {
            'model': 'limma-voom (TMM size factors, voom precision weights, moderated t)',
            'implementation': 'stage2/util/output/api/stage2_util/derivation/disease_deg.R',
            'implementation_sha256': sha256(Path(disease_deg.DEG_WORKER)),
            'test': 'moderated t, BH-adjusted p-value',
            'contrast_direction': 'positive log2FC = up in the disease arm vs Donor',
            'limma_version': summary['limma_version'],
            'r_version': summary['r_version'],
        },
        'contrasts': counts,
        'fit_summary': summary,
        'deg': {'path': 'local/deg/gse270788_deg.csv', 'sha256': sha256(DEG),
                'rows': int(len(table)), 'n_source_genes': int(table.source_gene_id.nunique())},
    }
    SUMMARY.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS', 'rows': len(table),
                      'contrasts': sorted(table.contrast.unique().tolist()),
                      'deg_sha256': record['deg']['sha256']}, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
