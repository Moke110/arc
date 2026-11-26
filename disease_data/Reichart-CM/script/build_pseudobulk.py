#!/usr/bin/env python
"""Reichart-CM — CM gate, per-donor pseudobulk and the genotype split (ADR 0066).

The cellxgene object is a CM-only h5ad (311,418 cells, 32,383 genes; 79 donors)
from the Reichart et al. cardiomyopathy collection. The object's own disease
labels are not used to define the ``nHF`` evidence: donors are split by
``Primary.Genetic.Diagnosis`` (genotype first, as ADR 0066 requires), so a DCM
donor carrying a declared pathogenic variant leaves ``nHF`` for ticket #188,
the 8 ARVC donors go to the #188 export, and the single NCCM donor is not used.

The raw counts are ``raw/X`` (CSR over cells). Each donor's CM counts are summed
into one pseudobulk profile **before** testing; the matrix is streamed in cell
chunks and accumulated in the linear space, never materialised as a dense
cell x gene matrix. The pseudobulk is then aligned to the frozen 18,164-gene
axis (Ensembl first, then symbol) and summed over aliases.

Writes (ADR 0030; tracked compact export is the ticket #188 input):

    local/pseudobulk/reichart-cm_cm_pseudobulk.csv.gz   axis genes x 78 donors, int
    local/genotype/reichart-cm_donor_genotypes.csv      per-donor split + cells
    local/pseudobulk/reichart-cm_cm_pseudobulk.json     tracked summary (digests)

Usage (arc python):
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
REPO = ITEM
while not (REPO / '.git').exists():
    if REPO.parent == REPO:
        raise SystemExit('cannot locate the repository root')
    REPO = REPO.parent
sys.path.insert(0, str(REPO / 'stage2_CarPS-VAE/util/output/api'))
from stage2_util.derivation import disease_deg  # noqa: E402

SOURCE = ITEM / 'source'
OBJECT = SOURCE / '27fdbdf9-deae-4aed-bf7f-21964ee62137.h5ad'
LOCAL = ITEM / 'local'
AXIS = REPO / 'stage2_CarPS-VAE/util/output/gene-identity/genes.csv'
CROSSWALK = REPO / 'stage2_CarPS-VAE/util/output/gene-identity/tables/gene_crosswalk.csv'
GATE = ('cell_type', 'cardiac muscle cell')
CELL_CHUNK = 4000
SCHEMA = 'stage2-item-cm-pseudobulk-v1'
DONOR_STATUS = {
    'normal': ('negative', 'control', 'normal'),
    'arvc': ('pathogenic', 'export_only', 'arvc'),
    'nccm': ('pathogenic', 'excluded', 'nccm'),
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def categorical(group):
    categories = [value.decode() if isinstance(value, bytes) else str(value)
                  for value in group['categories'][:]]
    return categories, group['codes'][:]


def donor_labels(handle):
    obs = handle['obs']
    donor_categories, donor_codes = categorical(obs['donor_id'])
    disease_categories, disease_codes = categorical(obs['disease'])
    genotype_categories, genotype_codes = categorical(obs['Primary.Genetic.Diagnosis'])
    cell_categories, cell_codes = categorical(obs['cell_type'])
    if [value for value in cell_categories] != [GATE[1]]:
        raise SystemExit(f'unexpected cell_type values: {cell_categories}')
    records = []
    for donor in donor_categories:
        positions = np.flatnonzero(donor_codes == donor_categories.index(donor))
        disease = disease_categories[disease_codes[positions[0]]]
        genotype = genotype_categories[genotype_codes[positions[0]]]
        if not np.all(disease_codes[positions] == disease_codes[positions[0]]):
            raise SystemExit(f'donor {donor} has more than one disease label')
        if not np.all(genotype_codes[positions] == genotype_codes[positions[0]]):
            raise SystemExit(f'donor {donor} has more than one genotype label')
        if disease == 'normal':
            label = 'normal'
        elif disease == 'arrhythmogenic right ventricular cardiomyopathy':
            label = 'arvc'
        elif disease == 'non-compaction cardiomyopathy':
            label = 'nccm'
        elif disease == 'dilated cardiomyopathy':
            label = 'dcm_genotype_negative' if genotype == 'PVneg' else 'dcm_genotype_positive'
        else:
            raise SystemExit(f'unrecognised disease label: {disease}')
        records.append({'donor_id': donor, 'disease': disease,
                        'primary_genetic_diagnosis': genotype,
                        'genotype': label, 'cells': int(len(positions))})
    frame = pd.DataFrame(records)
    split, summary = disease_deg.genotype_split(frame[['donor_id', 'genotype']])
    frame = frame.merge(split[['donor_id', 'genotype_status', 'nHF_role']], on='donor_id')
    return frame, summary


def pseudobulk(handle, donor_ids):
    matrix = handle['raw/X']
    indptr = matrix['indptr'][:]
    n_cells = len(indptr) - 1
    position = {donor: index for index, donor in enumerate(donor_ids)}
    donor_categories, donor_codes = categorical(handle['obs/donor_id'])
    cell_donor = np.asarray([position[donor_categories[code]] for code in donor_codes],
                            dtype=np.int64)
    if len(cell_donor) != n_cells:
        raise SystemExit('donor labels and raw/X cells disagree')
    var = handle['raw/var/_index'][:]
    genes = np.array([value.decode() if isinstance(value, bytes) else str(value)
                      for value in var])
    n_genes = len(genes)
    accumulator = np.zeros((len(donor_ids), n_genes), dtype=np.int64)
    n_entries = 0
    for start in range(0, n_cells, CELL_CHUNK):
        end = min(start + CELL_CHUNK, n_cells)
        first, last = int(indptr[start]), int(indptr[end])
        values = matrix['data'][first:last]
        indices = matrix['indices'][first:last]
        n_entries += int(len(values))
        lengths = np.diff(indptr[start:end + 1])
        cells = np.repeat(np.arange(start, end), lengths)
        donor = cell_donor[cells]
        flat = donor * n_genes + indices.astype(np.int64)
        block = np.bincount(flat, weights=values, minlength=len(donor_ids) * n_genes)
        accumulator += block.astype(np.int64).reshape(len(donor_ids), n_genes)
    expected = int(indptr[-1])
    if n_entries != expected:
        raise SystemExit(f'raw/X entry count {n_entries} != indptr[-1] {expected}')
    return genes, accumulator


def main():
    if not OBJECT.is_file():
        raise SystemExit('run script/fetch.sh first')
    with h5py.File(OBJECT, 'r') as handle:
        donors, split_summary = donor_labels(handle)
        donors = donors.sort_values('donor_id', kind='stable').reset_index(drop=True)
        exported = donors[donors.nHF_role != 'excluded'].reset_index(drop=True)
        genes, all_counts = pseudobulk(handle, donors.donor_id.tolist())
    donor_position = {donor: index for index, donor in enumerate(donors.donor_id)}
    counts = all_counts[[donor_position[donor] for donor in exported.donor_id], :]
    if int((exported.nHF_role == 'case').sum()) == 0:
        raise SystemExit('no genotype-negative DCM donor for nHF')
    if int((exported.nHF_role == 'control').sum()) == 0:
        raise SystemExit('no normal donor for nHF')
    axis = disease_deg.read_axis(AXIS)
    crosswalk = disease_deg.read_crosswalk(CROSSWALK)
    resolver = disease_deg.build_resolver(axis, crosswalk)
    indices, categories = resolver.resolve_many(genes.tolist())
    report = resolver.loss_summary(indices, categories)
    aligned = np.zeros((len(exported), len(axis)), dtype=np.int64)
    groups = {}
    for position, gene_index in enumerate(indices):
        if gene_index >= 0:
            groups.setdefault(int(gene_index), []).append(position)
    for gene_index, positions in groups.items():
        aligned[:, gene_index] += counts[:, positions].sum(axis=1)
    frame = pd.DataFrame({'gene_index': axis.gene_index.to_numpy(dtype=np.int64),
                          'gene_symbol': axis.symbol.to_numpy(dtype=object)})
    for position, donor in enumerate(exported.donor_id):
        frame[donor] = aligned[position, :]
    LOCAL.joinpath('pseudobulk').mkdir(parents=True, exist_ok=True)
    LOCAL.joinpath('genotype').mkdir(parents=True, exist_ok=True)
    pseudobulk_path = LOCAL / 'pseudobulk' / 'reichart-cm_cm_pseudobulk.csv.gz'
    genotype_path = LOCAL / 'genotype' / 'reichart-cm_donor_genotypes.csv'
    frame.to_csv(pseudobulk_path, index=False, compression='gzip')
    genotype_table = exported[['donor_id', 'disease', 'primary_genetic_diagnosis',
                               'genotype', 'genotype_status', 'nHF_role', 'cells']]
    genotype_table = genotype_table.assign(
        fit_group=genotype_table.genotype,
        exported_as=genotype_table.nHF_role.map(
            {'control': 'nHF control', 'case': 'nHF case',
             'export_only': 'ticket-188 export'}))
    genotype_table.to_csv(genotype_path, index=False)
    summary = {
        'schema': SCHEMA,
        'item': 'Reichart-CM',
        'collection_id': 'e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb',
        'dataset_version_id': '27fdbdf9-deae-4aed-bf7f-21964ee62137',
        'gate': f'{GATE[0]} == "{GATE[1]}"',
        'n_cells': int(exported.cells.sum()),
        'n_cells_all_donors': int(donors.cells.sum()),
        'n_donors': int(len(donors)),
        'n_donors_exported': int(len(exported)),
        'excluded_donors': donors.loc[donors.nHF_role == 'excluded', 'donor_id'].tolist(),
        'cells_per_donor': {row.donor_id: int(row.cells) for row in donors.itertuples()},
        'genotype_status_counts': split_summary,
        'nHF_case_donors': int((exported.nHF_role == 'case').sum()),
        'nHF_control_donors': int((exported.nHF_role == 'control').sum()),
        'genotype_export_donors': int((exported.nHF_role == 'export_only').sum()),
        'n_genes_native': int(len(genes)),
        'n_axis_genes_nonzero': int((frame[exported.donor_id.tolist()].sum(axis=1) > 0).sum()),
        'pseudobulk_dtype': 'integer counts summed per donor before testing',
        'alignment': report,
        'pseudobulk_sha256': sha256(pseudobulk_path),
        'genotype_sha256': sha256(genotype_path),
    }
    summary_path = LOCAL / 'pseudobulk' / 'reichart-cm_cm_pseudobulk.json'
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'result': 'PASS', 'donors': summary['n_donors_exported'],
                      'cells': summary['n_cells'],
                      'nHF_case': summary['nHF_case_donors'],
                      'nHF_control': summary['nHF_control_donors'],
                      'export': summary['genotype_export_donors'],
                      'alignment': report['n_resolved'],
                      'pseudobulk_sha256': summary['pseudobulk_sha256']}, indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
