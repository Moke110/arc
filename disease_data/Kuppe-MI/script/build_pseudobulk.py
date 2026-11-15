#!/usr/bin/env python
"""Kuppe-MI — CM gate and per-(donor, zone) cardiomyocyte pseudobulk (ADR 0030).

The cellxgene object `All-snRNA-Spatial multi-omic map of human myocardial
infarction` (Kuppe et al., Nature 2022) is a 191,795-nucleus, 20-donor, 28,380-
feature left-ventricle snRNA object with 11 author cell classes. The disease
evidence this item carries is cardiomyocyte-specific (ADR 0066/0071), so the
object is gated on its own author label ``cell_type_original == "Cardiomyocyte"``
(64,510 nuclei) and nothing is re-clustered.

The object is also **zoned**: the authors' ``major_labl`` column assigns every
nucleus to one physiological zone of the infarcted left ventricle — ``CTRL``
(control hearts), ``IZ`` (ischaemic/infarct zone), ``BZ`` (border zone), ``RZ``
(remote zone) and ``FZ`` (fibrotic zone). The zone is a property of the nucleus,
so one donor can contribute several zone profiles; the point of this source is
that the same heart's infarct zone is compared with the control hearts'
myocardium, which is why the profiles are built per (donor, zone) rather than
pooled per donor (pooling a donor's four zones washes the ischaemic signal out;
see the item's author_log notes and the #202 audit).

Counts are the deposited integer counts in ``raw/X`` (cellxgene
``raw_data_location: raw.X``, CSR over cells). Each profile is the sum of its
nuclei's counts **before** any testing: ``raw/X`` is streamed in 4,000-cell
chunks and accumulated in the linear space, so no dense cell x gene matrix is
ever materialised and non-CM entries are skipped at the entry level.

Writes (ADR 0030; the matrix and the profile sheet are regenerable ``local/``
additions, the JSON summary is tracked):

    local/pseudobulk/kuppe-mi_cm_pseudobulk.csv.gz  28,380 genes x 25 profiles
    local/pseudobulk/kuppe-mi_profiles.csv           profile -> donor/zone/CM counts
    local/pseudobulk/kuppe-mi_cm_pseudobulk.json     tracked summary (digests)

Usage (arc python, after script/fetch.sh && script/verify.sh):

    /home/chang/miniconda3/envs/arc/bin/python script/build_pseudobulk.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ITEM = Path(__file__).resolve().parents[1]
# ADR 0030: discover the repository root by the .git directory, never by counting parents.
REPO = ITEM
while not (REPO / '.git').exists():
    if REPO.parent == REPO:
        raise SystemExit('cannot locate the repository root')
    REPO = REPO.parent

SOURCE = ITEM / 'source'
OBJECT = SOURCE / '4780eb7d-4c7c-4630-8705-471f0de8e27a.h5ad'
LOCAL = ITEM / 'local'
SCHEMA = 'stage2-item-cm-pseudobulk-v1'
CELL_TYPE_COLUMN = 'cell_type_original'
CELL_TYPE_LABEL = 'Cardiomyocyte'
ZONE_COLUMN = 'major_labl'
ZONES = ('CTRL', 'IZ', 'BZ', 'RZ', 'FZ')
CELL_CHUNK = 4000
EXPECTED_NUCLEI = 191795
EXPECTED_CM_NUCLEI = 64510
EXPECTED_FEATURES = 28380
EXPECTED_DONORS = 20


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def categorical(handle, column):
    """Return ``(categories, codes)`` of one categorical obs/var column."""
    group = handle[column]
    if isinstance(group, h5py.Dataset):
        values = [value.decode() if isinstance(value, bytes) else str(value)
                  for value in group[:]]
        return values, np.arange(len(values), dtype=np.int64)
    categories = [value.decode() if isinstance(value, bytes) else str(value)
                  for value in group['categories'][:]]
    return categories, np.asarray(group['codes'][:], dtype=np.int64)


def donor_and_zone(handle):
    """Resolve the CM gate, the donor/zone labels and the (donor, zone) profiles."""
    cell_categories, cell_codes = categorical(handle, f'obs/{CELL_TYPE_COLUMN}')
    if CELL_TYPE_LABEL not in cell_categories:
        raise SystemExit(f'{CELL_TYPE_COLUMN} has no {CELL_TYPE_LABEL!r}: {cell_categories}')
    is_cm = cell_codes == cell_categories.index(CELL_TYPE_LABEL)
    if int(is_cm.sum()) != EXPECTED_CM_NUCLEI:
        raise SystemExit(f'CM gate selects {int(is_cm.sum())} nuclei, '
                         f'expected {EXPECTED_CM_NUCLEI}')
    donor_categories, donor_codes = categorical(handle, 'obs/donor_id')
    disease_categories, disease_codes = categorical(handle, 'obs/disease')
    zone_categories, zone_codes = categorical(handle, f'obs/{ZONE_COLUMN}')
    if sorted(zone_categories) != sorted(ZONES):
        raise SystemExit(f'unexpected {ZONE_COLUMN} values: {zone_categories}')
    if len(donor_categories) != EXPECTED_DONORS:
        raise SystemExit(f'expected {EXPECTED_DONORS} donors, found {len(donor_categories)}')
    donor_disease = {}
    for index, donor in enumerate(donor_categories):
        labels = set(disease_categories[code]
                     for code in disease_codes[donor_codes == index])
        if len(labels) != 1:
            raise SystemExit(f'donor {donor} carries more than one disease label: {labels}')
        donor_disease[donor] = labels.pop()
    # One profile per (donor, zone) that actually holds CM nuclei, in the
    # lexicographic profile order ("P1|CTRL", "P10|IZ", ...).
    profiles = []
    for donor_index, donor in enumerate(donor_categories):
        for zone in ZONES:
            zone_index = zone_categories.index(zone)
            n_cm = int(((donor_codes == donor_index) & (zone_codes == zone_index) & is_cm).sum())
            if n_cm:
                profiles.append({'profile': f'{donor}|{zone}', 'donor': donor,
                                 'zone': zone, 'disease': donor_disease[donor],
                                 'n_cm': n_cm})
    profiles.sort(key=lambda row: row['profile'])
    return is_cm, donor_categories, donor_codes, zone_categories, zone_codes, profiles


def pseudobulk(handle, is_cm, donor_categories, donor_codes, zone_categories,
               zone_codes, profiles):
    """Accumulate ``raw/X`` over each profile's CM nuclei; never densify cells x genes."""
    position = {row['profile']: index for index, row in enumerate(profiles)}
    matrix = handle['raw/X']
    indptr = np.asarray(matrix['indptr'][:], dtype=np.int64)
    n_cells = len(indptr) - 1
    if n_cells != EXPECTED_NUCLEI:
        raise SystemExit(f'raw/X holds {n_cells} cells, expected {EXPECTED_NUCLEI}')
    identifiers, identifier_codes = categorical(handle, 'var/_index')
    feature_categories, feature_codes = categorical(handle, 'var/feature_name')
    if len(identifiers) != EXPECTED_FEATURES:
        raise SystemExit(f'expected {EXPECTED_FEATURES} features, found {len(identifiers)}')
    n_genes = len(identifiers)
    # Per-nucleus profile assignment (donor x zone); -1 for a non-CM nucleus.
    profile_of_cell = np.full(n_cells, -1, dtype=np.int64)
    donor_index_of = {donor: index for index, donor in enumerate(donor_categories)}
    zone_index_of = {zone: index for index, zone in enumerate(zone_categories)}
    for row in profiles:
        mask = ((donor_codes == donor_index_of[row['donor']])
                & (zone_codes == zone_index_of[row['zone']]) & is_cm)
        if int(mask.sum()) != row['n_cm']:
            raise SystemExit(f'profile {row["profile"]}: CM nuclei disagree with the gate')
        profile_of_cell[mask] = position[row['profile']]
    if int((profile_of_cell >= 0).sum()) != int(is_cm.sum()):
        raise SystemExit('profile assignment lost CM nuclei')
    accumulator = np.zeros((len(profiles), n_genes), dtype=np.int64)
    n_entries = 0
    for start in range(0, n_cells, CELL_CHUNK):
        end = min(start + CELL_CHUNK, n_cells)
        first, last = int(indptr[start]), int(indptr[end])
        values = np.asarray(matrix['data'][first:last], dtype=np.float64)
        indices = np.asarray(matrix['indices'][first:last], dtype=np.int64)
        n_entries += int(len(values))
        lengths = np.diff(indptr[start:end + 1])
        cells = np.repeat(np.arange(start, end, dtype=np.int64), lengths)
        profile = profile_of_cell[cells]
        keep = profile >= 0                       # skip every non-CM entry
        if not keep.any():
            continue
        flat = profile[keep] * n_genes + indices[keep]
        block = np.bincount(flat, weights=values[keep],
                            minlength=len(profiles) * n_genes).reshape(len(profiles), n_genes)
        if np.abs(block - np.rint(block)).max() > 1e-6:
            raise SystemExit('raw/X holds non-integer values; expected integer counts')
        accumulator += np.rint(block).astype(np.int64)
    if n_entries != int(indptr[-1]):
        raise SystemExit(f'raw/X entry count {n_entries} != indptr[-1] {int(indptr[-1])}')
    identifiers = [identifiers[code] for code in identifier_codes]
    symbols = [feature_categories[code] for code in feature_codes]
    return identifiers, symbols, accumulator


def main():
    if not OBJECT.is_file():
        raise SystemExit('run script/fetch.sh first')
    with h5py.File(OBJECT, 'r') as handle:
        gated = donor_and_zone(handle)
        is_cm, donor_categories, donor_codes, zone_categories, zone_codes, profiles = gated
        identifiers, symbols, counts = pseudobulk(
            handle, is_cm, donor_categories, donor_codes, zone_categories, zone_codes,
            profiles)
    if not profiles:
        raise SystemExit('no CM profile')
    per_zone = {zone: int(sum(row['n_cm'] for row in profiles if row['zone'] == zone))
                for zone in ZONES}
    if sum(per_zone.values()) != int(is_cm.sum()):
        raise SystemExit('zone totals disagree with the CM gate')
    frame = pd.DataFrame({'source_gene_id': identifiers, 'gene_symbol': symbols})
    for index, row in enumerate(profiles):
        frame[row['profile']] = counts[index, :]
    profile_frame = pd.DataFrame(
        [{'profile': row['profile'], 'donor': row['donor'], 'zone': row['zone'],
          'disease': row['disease'], 'n_cm': row['n_cm'],
          'total_counts': int(counts[index, :].sum()),
          'n_genes_detected': int((counts[index, :] > 0).sum())}
         for index, row in enumerate(profiles)])
    LOCAL.joinpath('pseudobulk').mkdir(parents=True, exist_ok=True)
    pseudobulk_path = LOCAL / 'pseudobulk' / 'kuppe-mi_cm_pseudobulk.csv.gz'
    profiles_path = LOCAL / 'pseudobulk' / 'kuppe-mi_profiles.csv'
    frame.to_csv(pseudobulk_path, index=False, compression='gzip')
    profile_frame.to_csv(profiles_path, index=False)
    summary = {
        'schema': SCHEMA,
        'item': 'Kuppe-MI',
        'collection_id': '8191c283-0816-424b-9b61-c3e1d6258a77',
        'dataset_id': '1c739a3e-c3f5-49d5-98e0-73975e751201',
        'dataset_version_id': '4780eb7d-4c7c-4630-8705-471f0de8e27a',
        'gate': f'{CELL_TYPE_COLUMN} == "{CELL_TYPE_LABEL}"',
        'zone_column': ZONE_COLUMN,
        'n_cells': int(EXPECTED_NUCLEI),
        'n_cm_cells': int(is_cm.sum()),
        'n_features': int(len(identifiers)),
        'n_donors': int(EXPECTED_DONORS),
        'cm_cells_per_zone': per_zone,
        'n_profiles': int(len(profiles)),
        'profiles': [{'profile': row['profile'], 'donor': row['donor'],
                      'zone': row['zone'], 'disease': row['disease'],
                      'n_cm': int(row['n_cm'])} for row in profiles],
        'counts': ('integer counts summed per (donor, zone) profile from raw/X over '
                   "that profile's CM nuclei, before any testing"),
        'pseudobulk_sha256': sha256(pseudobulk_path),
        'profiles_sha256': sha256(profiles_path),
    }
    summary_path = LOCAL / 'pseudobulk' / 'kuppe-mi_cm_pseudobulk.json'
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS', 'cm_cells': summary['n_cm_cells'],
                      'profiles': summary['n_profiles'], 'per_zone': per_zone,
                      'genes': summary['n_features'],
                      'pseudobulk_sha256': summary['pseudobulk_sha256']}, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
