"""Pinned source and eligibility records for ticket #187's disease categories.

ADR 0066 D10/D11: one category is one fusion source, and each source records its
species, model system, ``evidence_type``, priority ``tier``,
``shared_control_group``, ``genotype_status``, ``time_window_basis`` and
``exploratory`` flag. This module is the single tracked declaration the release
driver turns into each category's ``catalog.json``; eligibility is data, not
prose, so the published catalog carries it verbatim.

The ticket-187 categories are ``nHF``, ``ICM`` and ``AF`` only. The six isogenic
models (#188), the drug-class sources (#189) and the Reichart genotype-positive /
ARVC DEG (also #188) are other tickets and are not declared here.
"""
from __future__ import annotations

CATEGORIES = {
    'nHF': {
        'directory': 'nhf',
        'tier': 1,
        'evidence_type': 'cm_pseudobulk',
        'shared_control_group': ['gse270788:donor', 'reichart-cm:normal'],
        'rejected': [
            {'candidate': 'GSE183852',
             'reason': 'eligible human LV CM source, but CM labels live only inside '
                       '2.1-12 GB Seurat Robj objects; ad-hoc re-clustering is forbidden'},
            {'candidate': 'cellxgene myocarditis CM object',
             'reason': 'shares its normal donors with the Reichart object, so the two '
                       'cannot both count as independent sources'},
            {'candidate': 'GSE141910 / GSE116250 / GSE128188 / GSE282618',
             'reason': 'whole-heart/whole-ventricle bulk tissue, not cardiomyocyte-derived'},
        ],
    },
    'ICM': {
        'directory': 'icm',
        'tier': 1,
        'evidence_type': 'cm_pseudobulk',
        'shared_control_group': ['gse270788:donor'],
        'rejected': [
            {'candidate': 'GSE298023',
             'reason': 'labels only inside a 4.0 GB .rds; ad-hoc re-clustering is forbidden'},
            {'candidate': 'GSE145154',
             'reason': 'no author cell labels and region-heterogeneous; re-clustering '
                       'not permitted'},
            {'candidate': 'GSE141910 / GSE116250',
             'reason': 'whole-tissue bulk, not cardiomyocyte-derived'},
        ],
    },
    'AF': {
        'directory': 'af',
        'tier': 1,
        'evidence_type': 'cm_pseudobulk',
        'shared_control_group': ['gse255612:donor'],
        'rejected': [
            {'candidate': 'GSE128188',
             'reason': 'whole atrial appendage bulk tissue, not cardiomyocyte-derived'},
        ],
    },
}

# One record per published dataset; `contrasts` are labels produced by the item's
# single joint fit. `n_case`/`n_control` are asserted against the item summaries
# by the driver's tests.
DATASETS = [
    {
        'dataset_id': 'gse270788__nicm',
        'category': 'nHF',
        'source_item': 'disease_data/GSE270788',
        'local_deg': 'disease_data/GSE270788/local/deg/gse270788_deg.csv',
        'local_summary': 'disease_data/GSE270788/local/deg/gse270788_deg.json',
        'cells_per_donor': 'disease_data/GSE270788/local/pseudobulk/gse270788_donors.csv',
        'arm_or_class': 'NICM arm (6 donors) vs the seven shared Donor controls',
        'species': 'Homo sapiens',
        'model_system': 'human left-ventricle 10x Multiome CM pseudobulk',
        'evidence_type': 'cm_pseudobulk',
        'tier': 1,
        'genotype_status': 'not_declared',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse270788:donor',
        'n_case': 6,
        'n_control': 7,
        'contrasts': ['NICM-Donor'],
        'method': {
            'deposition': 'integer counts',
            'model': 'limma-voom (voom library-size normalisation, moderated t)',
            'design': 'one joint fit ~0 + arm over Donor/NICM/ICM/AMI (shared Donor controls)',
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'gse270788__icm',
        'category': 'ICM',
        'source_item': 'disease_data/GSE270788',
        'local_deg': 'disease_data/GSE270788/local/deg/gse270788_deg.csv',
        'local_summary': 'disease_data/GSE270788/local/deg/gse270788_deg.json',
        'cells_per_donor': 'disease_data/GSE270788/local/pseudobulk/gse270788_donors.csv',
        'arm_or_class': 'ICM (chronic ischemic) arm (5 donors) vs the seven shared Donor controls',
        'species': 'Homo sapiens',
        'model_system': 'human left-ventricle 10x Multiome CM pseudobulk',
        'evidence_type': 'cm_pseudobulk',
        'tier': 1,
        'genotype_status': 'not_declared',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse270788:donor',
        'n_case': 5,
        'n_control': 7,
        'contrasts': ['ICM-Donor'],
        'method': {
            'deposition': 'integer counts',
            'model': 'limma-voom (voom library-size normalisation, moderated t)',
            'design': 'one joint fit ~0 + arm over Donor/NICM/ICM/AMI (shared Donor controls)',
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'gse270788__ami',
        'category': 'ICM',
        'source_item': 'disease_data/GSE270788',
        'local_deg': 'disease_data/GSE270788/local/deg/gse270788_deg.csv',
        'local_summary': 'disease_data/GSE270788/local/deg/gse270788_deg.json',
        'cells_per_donor': 'disease_data/GSE270788/local/pseudobulk/gse270788_donors.csv',
        'arm_or_class': 'AMI (acute MI) arm (5 donors) vs the seven shared Donor controls',
        'species': 'Homo sapiens',
        'model_system': 'human left-ventricle 10x Multiome CM pseudobulk',
        'evidence_type': 'cm_pseudobulk',
        'tier': 1,
        'genotype_status': 'not_declared',
        'time_window_basis': 'author_declared',
        'exploratory': False,
        'shared_control_group': 'gse270788:donor',
        'n_case': 5,
        'n_control': 7,
        'contrasts': ['AMI-Donor'],
        'method': {
            'deposition': 'integer counts',
            'model': 'limma-voom (voom library-size normalisation, moderated t)',
            'design': ('one joint fit ~0 + arm over Donor/NICM/ICM/AMI; the AMI arm is the '
                       'acute end of the ischemic continuum, no per-sample time-since-infarct '
                       'field exists so acute status is author-declared'),
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'reichart-cm__dcm_genotype_negative',
        'category': 'nHF',
        'source_item': 'disease_data/Reichart-CM',
        'local_deg': 'disease_data/Reichart-CM/local/deg/reichart-cm_deg.csv',
        'local_summary': 'disease_data/Reichart-CM/local/deg/reichart-cm_deg.json',
        'cells_per_donor': 'disease_data/Reichart-CM/local/genotype/reichart-cm_donor_genotypes.csv',
        'arm_or_class': 'genotype-negative DCM donors vs the object\'s normal donors',
        'species': 'Homo sapiens',
        'model_system': 'human LV/RV/apex/IVS sc+snRNA CM-only cellxgene object pseudobulk',
        'evidence_type': 'cm_pseudobulk',
        'tier': 1,
        'genotype_status': 'negative',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'reichart-cm:normal',
        'n_case': 0,
        'n_control': 0,
        'contrasts': ['dcm_genotype_negative-normal'],
        'method': {
            'deposition': 'integer counts',
            'model': 'limma-voom (voom library-size normalisation, moderated t)',
            'design': ('one joint fit ~0 + status over normal / genotype-negative DCM / '
                       'genotype-positive DCM / ARVC; positive and ARVC donors share the fit '
                       'for dispersion but publish no DEG in this ticket'),
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'gse255612__af',
        'category': 'AF',
        'source_item': 'disease_data/GSE255612',
        'local_deg': 'disease_data/GSE255612/local/deg/gse255612_deg.csv',
        'local_summary': 'disease_data/GSE255612/local/deg/gse255612_deg.json',
        'cells_per_donor': 'disease_data/GSE255612/local/pseudobulk/gse255612_donors.csv',
        'arm_or_class': 'AF donors vs non-AF control donors (left-atrial CM)',
        'species': 'Homo sapiens',
        'model_system': 'human left-atrium 10x snRNA CM pseudobulk',
        'evidence_type': 'cm_pseudobulk',
        'tier': 1,
        'genotype_status': 'not_declared',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse255612:donor',
        'n_case': 0,
        'n_control': 0,
        'contrasts': ['AF-control'],
        'method': {
            'deposition': 'normalised expression (log1p(CP10K))',
            'model': 'limma on log2(mean_CP10K + 1) pseudobulk (moderated t)',
            'design': 'one joint fit ~0 + group over 34 donors',
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
]


def datasets(category):
    """The declared datasets of one category, in catalog order."""
    selected = [record for record in DATASETS if record['category'] == category]
    if not selected:
        raise KeyError(f'no ticket-187 datasets declared for category {category!r}')
    return selected


def category_record(category):
    if category not in CATEGORIES:
        raise KeyError(f'unknown ticket-187 category {category!r}')
    return CATEGORIES[category]
