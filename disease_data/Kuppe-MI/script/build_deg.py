#!/usr/bin/env python
"""Kuppe-MI — the ICM gene-axis DEG: IZ vs CTRL cardiomyocyte pseudobulk (ADR 0071).

Design (declared in ``SPEC`` below, all of it recorded in the tracked summary):

* **Profiles.** One CM pseudobulk per (donor, zone) from
  ``script/build_pseudobulk.py``; the zone is the authors' ``major_labl``. A
  donor whose CM nuclei sit in several zones therefore contributes several
  profiles — that is intended, because the contrast of interest is the
  ischaemic zone of the infarcted hearts against the control hearts, not a
  donor-level average (pooling the four zones per donor washes the ischaemic
  signal out: the #202 audit measured 3 DEGs for the donor-pooled contrast and
  5,675 for the zone contrast).
* **CM floor.** A profile enters the fit only when it holds at least
  ``MIN_CM_NUCLEI = 500`` CM nuclei. That drops the sparse donor-zone profiles
  (P13/P16/P4/P18/P19/P20/P3|IZ/P4). After the floor the arms are CTRL 4 (the
  four normal donors) and IZ 4 (P2, P9, P10, P15), plus BZ 3, RZ 5 and FZ 2.
* **One joint design.** ``~0 + zone`` over all 18 profiles, so dispersion is
  estimated from every zone; the IZ-vs-CTRL contrast is the published one and
  the other three zone contrasts stay recorded artifacts.
* **Donor blocking.** The same donor contributes up to three zone profiles, so
  the repeated profiles are not independent: ``duplicateCorrelation(block =
  donor)`` estimates one consensus intra-donor correlation and both ``voom`` and
  ``lmFit`` are given it.
* **The fit** goes through the one live implementation,
  ``cpst/util/output/api/cpst_util/derivation/disease_deg.R``; the summary it prints (limma/R versions,
  consensus correlation) is carried into the tracked JSON.

The output is the gene-axis DEG interface: a full-transcriptome per-contrast
table (Ensembl id, deposit symbol, signed log2FC, moderated t, p, BH padj, the
voom average log-expression) for the published contrast, plus the tracked JSON.
The frozen-axis counts in the JSON are a cross-check computed with the same
``cpst_util.derivation.disease_deg`` resolver the republish ticket uses (symbol first, then
Ensembl); they do not filter the published table.

Writes (ADR 0030; the tables are regenerable ``local/`` additions):

    local/deg/kuppe-mi_zone_deg.csv            all four zone contrasts (artifact)
    local/deg/kuppe-mi_deg.csv                 the published IZ-CTRL contrast
    local/deg/kuppe-mi_samples.csv             the fit's design sheet
    local/deg/kuppe-mi_pseudobulk_matrix.tsv.gz the fit's matrix input
    local/deg/kuppe-mi_deg.json                tracked summary

Usage (arc python, after script/build_pseudobulk.py):

    /home/chang/miniconda3/envs/arc/bin/python script/build_deg.py
"""
from __future__ import annotations

import gzip
import hashlib
import io
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ITEM = Path(__file__).resolve().parents[1]
# ADR 0030: discover the repository root by the .git directory, never by counting parents.
REPO = ITEM
while not (REPO / '.git').exists():
    if REPO.parent == REPO:
        raise SystemExit('cannot locate the repository root')
    REPO = REPO.parent
sys.path.insert(0, str(REPO / 'cpst/util/output/api'))
from cpst_util.derivation import disease_deg  # noqa: E402

LOCAL = ITEM / 'local'
PSEUDOBULK = LOCAL / 'pseudobulk' / 'kuppe-mi_cm_pseudobulk.csv.gz'
PROFILES = LOCAL / 'pseudobulk' / 'kuppe-mi_profiles.csv'
MATRIX = LOCAL / 'deg' / 'kuppe-mi_pseudobulk_matrix.tsv.gz'
SAMPLES = LOCAL / 'deg' / 'kuppe-mi_samples.csv'
ZONE_DEG = LOCAL / 'deg' / 'kuppe-mi_zone_deg.csv'
DEG = LOCAL / 'deg' / 'kuppe-mi_deg.csv'
SUMMARY = LOCAL / 'deg' / 'kuppe-mi_deg.json'
AXIS = REPO / 'cpst/util/output/gene-identity/genes.csv'
CROSSWALK = REPO / 'cpst/util/output/gene-identity/tables/gene_crosswalk.csv'
SCHEMA = 'stage2-item-deg-summary-v1'
CATEGORY = 'ICM'
TIER = 1
MIN_CM_NUCLEI = 500
PUBLISHED_CONTRAST = 'IZ-CTRL'
ZONES = ('CTRL', 'IZ', 'BZ', 'RZ', 'FZ')
EXPECTED_ZONE_PROFILES = {'CTRL': 4, 'IZ': 4, 'BZ': 3, 'RZ': 5, 'FZ': 2}
SPEC = {
    'item': 'Kuppe-MI',
    'category': CATEGORY,
    'tier': TIER,
    'kind': 'counts',
    'formula': '~0 + zone',
    'factor': 'zone',
    'block': 'donor',
    'normalisation': 'voom',
    'filter': {'min_count': 10, 'min_total_count': 15, 'large_n': 10, 'min_prop': 0.7},
    'contrasts': [{'label': f'{zone}-CTRL', 'numerator': zone, 'denominator': 'CTRL'}
                  for zone in ZONES if zone != 'CTRL'],
}
# The design object handed to disease_deg.R is the SPEC minus the item bookkeeping.
DESIGN = {key: SPEC[key] for key in
          ('kind', 'formula', 'factor', 'block', 'normalisation', 'filter', 'contrasts')}
MODEL_DESCRIPTION = (
    'limma-voom (voom library-size normalisation, moderated t) with donor '
    'blocking: duplicateCorrelation(block=donor) on the joint ~0 + zone design, '
    'applied to both voom and lmFit')


def write_gzip_text(path, text):
    """Write gzip bytes with a fixed mtime so a rebuild is byte-identical (ADR 0030)."""
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", mtime=0) as handle:
        handle.write(text.encode())
    Path(path).write_bytes(buffer.getvalue())


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def contrasts_from_spec():
    return [contrast['label'] for contrast in SPEC['contrasts']]


def count_deg(rows, padj_column='padj', effect_column='log2fc'):
    significant = rows[padj_column] < 0.05
    passing = significant & (rows[effect_column].abs() >= 2)
    return {
        'n_tested': int(len(rows)),
        'n_padj_lt_0.05': int(significant.sum()),
        'n_v2_deg': int(passing.sum()),
        'n_up': int((passing & (rows[effect_column] > 0)).sum()),
        'n_down': int((passing & (rows[effect_column] < 0)).sum()),
    }


def main():
    if not PSEUDOBULK.is_file() or not PROFILES.is_file():
        raise SystemExit('run script/build_pseudobulk.py first')
    pseudobulk = pd.read_csv(PSEUDOBULK)
    profiles = pd.read_csv(PROFILES)
    if profiles.profile.duplicated().any():
        raise SystemExit('the profile sheet repeats a profile')
    if sorted(pseudobulk.columns[2:]) != sorted(profiles.profile):
        raise SystemExit('the pseudobulk columns and the profile sheet disagree')
    # The 500-CM-nuclei floor is the donor-profile inclusion rule.
    kept = profiles[profiles.n_cm >= MIN_CM_NUCLEI].reset_index(drop=True)
    zone_counts = {zone: int((kept.zone == zone).sum()) for zone in ZONES}
    if zone_counts != EXPECTED_ZONE_PROFILES:
        raise SystemExit(f'zone profile counts after the {MIN_CM_NUCLEI}-CM floor are '
                         f'{zone_counts}, expected {EXPECTED_ZONE_PROFILES}')
    # The floor may only drop sparse profiles; what it drops is recorded, not hidden.
    dropped = profiles[profiles.n_cm < MIN_CM_NUCLEI].reset_index(drop=True)
    if len(dropped) + len(kept) != len(profiles):
        raise SystemExit('the CM floor and the profile sheet disagree')
    if dropped.empty:
        raise SystemExit('the CM floor dropped nothing; the design would be the unfiltered one')
    for donor, group in kept.groupby('donor'):
        if group.zone.duplicated().any():
            raise SystemExit(f'donor {donor} has two profiles in one zone')
    samples = pd.DataFrame({'sample': kept.profile, 'zone': kept.zone,
                            'donor': kept.donor, 'n_cm': kept.n_cm})
    if sorted(samples['zone'].unique()) != sorted(ZONES):
        raise SystemExit('a zone disappeared after the CM floor')
    LOCAL.joinpath('deg').mkdir(parents=True, exist_ok=True)
    samples.to_csv(SAMPLES, index=False)
    # The fit matrix must be in sample-sheet order (the worker asserts it).
    matrix = pseudobulk[['source_gene_id', 'gene_symbol']].copy()
    for profile in samples['sample']:
        matrix[profile] = pseudobulk[profile]
    write_gzip_text(MATRIX, matrix.to_csv(index=False, sep='\t'))
    summary = disease_deg.run_deg_worker(MATRIX, SAMPLES, DESIGN, ZONE_DEG)
    table = pd.read_csv(ZONE_DEG)
    missing = set(disease_deg.LOCAL_DEG_COLUMNS) - set(table.columns)
    if missing:
        raise SystemExit(f'R output lacks {sorted(missing)}')
    if sorted(table.contrast.unique()) != sorted(contrasts_from_spec()):
        raise SystemExit('the joint fit published an unexpected contrast set')
    published = table[table.contrast == PUBLISHED_CONTRAST].copy()
    if published.empty:
        raise SystemExit(f'the joint fit produced no {PUBLISHED_CONTRAST} rows')
    published = published.reset_index(drop=True)
    published.to_csv(DEG, index=False)
    # Per-contrast counts, the published one first; the other zones are artifacts.
    deg_counts = {}
    for contrast in contrasts_from_spec():
        rows = table[table.contrast == contrast]
        counts = count_deg(rows)
        counts['published'] = contrast == PUBLISHED_CONTRAST
        deg_counts[contrast] = counts
    published_counts = deg_counts[PUBLISHED_CONTRAST]
    if published_counts['n_v2_deg'] == 0:
        raise SystemExit(f'{PUBLISHED_CONTRAST} produced no v2 DEG')
    # Frozen-axis cross-check through the same resolver the republish ticket uses.
    axis = disease_deg.read_axis(AXIS)
    resolver = disease_deg.build_resolver(axis, disease_deg.read_crosswalk(CROSSWALK))
    aligned, alignment = disease_deg.align_source_table(published, resolver)
    axis_counts = count_deg(aligned)
    axis_counts['n_axis_genes'] = int(len(axis))
    n_case = int((samples.zone == 'IZ').sum())
    n_control = int((samples.zone == 'CTRL').sum())
    record = {
        'schema': SCHEMA,
        'item': SPEC['item'],
        'category': CATEGORY,
        'tier': TIER,
        'gate': 'cell_type_original == "Cardiomyocyte" (author CM label; 64,510 nuclei)',
        'zone_column': 'major_labl',
        'published_contrast': PUBLISHED_CONTRAST,
        'artifact_contrasts': [contrast for contrast in contrasts_from_spec()
                               if contrast != PUBLISHED_CONTRAST],
        'profile_floor': {
            'rule': f'profile n_cm >= {MIN_CM_NUCLEI}',
            'n_profiles_total': int(len(profiles)),
            'n_profiles_fitted': int(len(kept)),
            'dropped_profiles': [{'profile': row.profile, 'donor': row.donor,
                                  'zone': row.zone, 'n_cm': int(row.n_cm)}
                                 for row in dropped.itertuples()],
            'zone_profile_counts': zone_counts,
        },
        'samples': {
            'n': int(len(samples)),
            'group_counts': {zone: int((samples.zone == zone).sum()) for zone in ZONES},
            'cm_nuclei_per_profile': {row.sample: int(row.n_cm) for row in samples.itertuples()},
            'donors_per_arm': {
                'IZ': sorted(samples.loc[samples.zone == 'IZ', 'donor'].tolist()),
                'CTRL': sorted(samples.loc[samples.zone == 'CTRL', 'donor'].tolist()),
            },
            'path': 'local/deg/kuppe-mi_samples.csv',
            'sha256': sha256(SAMPLES),
        },
        'pseudobulk': {
            'path': 'local/pseudobulk/kuppe-mi_cm_pseudobulk.csv.gz',
            'sha256': sha256(PSEUDOBULK),
            'dtype': "integer counts summed per (donor, zone) profile over that profile's CM nuclei",
            'n_genes': int(len(pseudobulk)),
            'n_profiles': int(len(profiles)),
        },
        'design': DESIGN,
        'spec': SPEC,
        'design_ambiguities': {
            'zone_not_donor': (
                'profiles are per (donor, zone), not per donor: an infarct heart '
                'contributes its IZ/BZ/RZ/FZ CM nuclei as separate profiles, which is '
                'what makes the IZ-vs-CTRL contrast possible; a donor-level pooled '
                'contrast was measured and kept out of the axis (3 DEGs)'),
            'repeated_donors': (
                'a donor contributes up to three zone profiles, so the profiles are not '
                'independent; duplicateCorrelation(block=donor) is estimated on the joint '
                'design and applied to voom and lmFit'),
            'published_contrast': (
                'only IZ-CTRL is published: the ischaemic zone is the authors\' infarct-zone '
                'annotation and the contrast the ADR 0071 ICM evidence rests on; BZ/RZ/FZ '
                'are fitted for dispersion and recorded as artifacts'),
            'time_window': (
                'the deposit carries no per-donor time-since-infarct field, so acute-vs-chronic '
                'is not declared; the record is time_window_basis: not_declared'),
        },
        'contrasts': {
            PUBLISHED_CONTRAST: {
                'n_case': n_case, 'n_control': n_control,
                'case_donors': sorted(samples.loc[samples.zone == 'IZ', 'donor'].tolist()),
                'control_donors': sorted(samples.loc[samples.zone == 'CTRL', 'donor'].tolist()),
            },
        },
        'deg_counts': {
            'published_contrast': published_counts,
            'all_contrasts': deg_counts,
            'axis': axis_counts,
            'axis_definition': 'frozen 18,164-gene axis, resolved by cpst_util.derivation.disease_deg',
        },
        'alignment': alignment,
        'method': {
            'model': (MODEL_DESCRIPTION),
            'implementation': 'cpst/util/output/api/cpst_util/derivation/disease_deg.R',
            'implementation_sha256': sha256(Path(disease_deg.DEG_WORKER)),
            'count_filter': SPEC['filter'],
            'normalisation': SPEC['normalisation'],
            'block': SPEC['block'],
            'test': 'moderated t, BH-adjusted p-value',
            'contrast_direction': ('positive log2FC = up in the IZ (ischaemic zone) CM '
                                   'profiles vs the control-heart CM profiles'),
            'limma_version': summary['limma_version'],
            'r_version': summary['r_version'],
        },
        'fit_summary': summary,
        'deg': {'path': 'local/deg/kuppe-mi_deg.csv', 'sha256': sha256(DEG),
                'rows': int(len(published)),
                'n_source_genes': int(published.source_gene_id.nunique())},
        'zone_artifact': {'path': 'local/deg/kuppe-mi_zone_deg.csv', 'sha256': sha256(ZONE_DEG),
                          'rows': int(len(table)),
                          'contrasts': sorted(table.contrast.unique())},
        'matrix': {'path': 'local/deg/kuppe-mi_pseudobulk_matrix.tsv.gz',
                   'sha256': sha256(MATRIX)},
    }
    SUMMARY.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS',
                      'published_contrast': PUBLISHED_CONTRAST,
                      'n_tested': published_counts['n_tested'],
                      'padj_lt_0.05': published_counts['n_padj_lt_0.05'],
                      'v2_deg': published_counts['n_v2_deg'],
                      'v2_up': published_counts['n_up'], 'v2_down': published_counts['n_down'],
                      'axis_padj': axis_counts['n_padj_lt_0.05'],
                      'axis_v2': axis_counts['n_v2_deg'],
                      'limma': summary['limma_version'],
                      'deg_sha256': record['deg']['sha256']}, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
