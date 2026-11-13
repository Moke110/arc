"""Pinned source and eligibility records for ticket #187's disease categories.

ADR 0066 D10/D11: one category is one fusion source, and each source records its
species, model system, ``evidence_type``, priority ``tier``,
``shared_control_group``, ``genotype_status``, ``time_window_basis`` and
``exploratory`` flag. This module is the single tracked declaration the release
driver turns into each category's ``catalog.json``; eligibility is data, not
prose, so the published catalog carries it verbatim.

Local input paths are **item-relative** here (``source_item`` names the ADR
0030 directory, ``local_deg``/``cells_per_donor``/``local_summary`` are relative
to it). The driver composes the repo-relative path at runtime, so a merged
checkout does not carry a repo-relative literal for an item's gitignored
``local/`` payload (``script/check_paths.sh`` must stay green without the
per-worktree derived files). Every derived input is still pinned by sha256 in
the published catalog/manifest.

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
        'local_deg': 'local/deg/gse270788_deg.csv',
        'local_summary': 'local/deg/gse270788_deg.json',
        'cells_per_donor': 'local/pseudobulk/gse270788_donors.csv',
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
        'local_deg': 'local/deg/gse270788_deg.csv',
        'local_summary': 'local/deg/gse270788_deg.json',
        'cells_per_donor': 'local/pseudobulk/gse270788_donors.csv',
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
        'local_deg': 'local/deg/gse270788_deg.csv',
        'local_summary': 'local/deg/gse270788_deg.json',
        'cells_per_donor': 'local/pseudobulk/gse270788_donors.csv',
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
        'local_deg': 'local/deg/reichart-cm_deg.csv',
        'local_summary': 'local/deg/reichart-cm_deg.json',
        'cells_per_donor': 'local/genotype/reichart-cm_donor_genotypes.csv',
        'arm_or_class': 'genotype-negative DCM donors vs the object\'s normal donors',
        'species': 'Homo sapiens',
        'model_system': 'human LV/RV/apex/IVS sc+snRNA CM-only cellxgene object pseudobulk',
        'evidence_type': 'cm_pseudobulk',
        'tier': 1,
        'genotype_status': 'negative',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'reichart-cm:normal',
        'n_case': 8,
        'n_control': 18,
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
        'local_deg': 'local/deg/gse255612_deg.csv',
        'local_summary': 'local/deg/gse255612_deg.json',
        'cells_per_donor': 'local/pseudobulk/gse255612_donors.csv',
        'arm_or_class': 'AF donors vs non-AF control donors (left-atrial CM)',
        'species': 'Homo sapiens',
        'model_system': 'human left-atrium 10x snRNA CM pseudobulk',
        'evidence_type': 'cm_pseudobulk',
        'tier': 1,
        'genotype_status': 'not_declared',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse255612:donor',
        'n_case': 18,
        'n_control': 16,
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

# --- ticket #188: genetic_cardiomyopathy -------------------------------------
# One category, eight datasets (ADR 0066 *Category roster*): the two Reichart
# genotype-split human CM-pseudobulk exports (tier 1) plus the six isogenic
# hiPSC-CM bulk models (tier 2) that the ADR's *Isogenic exception* keeps as
# their own evidence type. Every local path here is item-relative (composed from
# ``source_item`` by the driver), so no repo-relative literal points at a
# gitignored ``local/`` payload and ``script/check_paths.sh`` stays green in a
# merged checkout.

CATEGORIES['genetic_cardiomyopathy'] = {
    'directory': 'genetic_cardiomyopathy',
    'tier': 1,
    'evidence_type': 'cm_pseudobulk',  # the winning tier; datasets carry their own
    'shared_control_group': ['reichart-cm:normal'],
    'isogenic_exception': {
        'applies': True,
        'statement': (
            'ADR 0066 Isogenic exception: the six in-vitro hiPSC-CM models carry a patient '
            'genotype with an isogenic control and are kept as their own evidence type inside '
            'this category even though a human CM-pseudobulk cohort (Reichart-CM) exists for '
            'the same category: a mixed-etiology clinical cohort and a single-genotype '
            'isogenic model answer different questions, and collapsing them loses the causal '
            'signature. The category therefore has one tier-1 human CM-pseudobulk pair '
            '(Reichart genotype-positive DCM + ARVC) and six tier-2 isogenic bulk datasets; '
            '`tier` records the winning tier and each dataset record carries its own '
            '`evidence_type` and `tier`.'),
    },
    'rejected': [
        {'candidate': 'GSE140965',
         'reason': 'eligible isogenic HCM MYBPC3 deposit, but its control sets are batch '
                   'confounded and the six isogenic models already cover MYBPC3 HCM'},
        {'candidate': 'GSE251993 / GSE262876 / GSE274124 / GSE221945 / GSE329393',
         'reason': 'eligible isogenic bulk models, but the declared roster fixes the six '
                   'isogenic accessions; admitting more would exceed the survey fetch list'},
        {'candidate': 'GSE289628 / GSE248187',
         'reason': 'ARVC DSP models: GSE289628 is confounded by AAV6-vs-siRNA vehicle and '
                   'GSE248187 deposits a DEG table only, so no counts could be re-fit'},
        {'candidate': 'GSE141910 / GSE116250 / GSE128188 / GSE282618',
         'reason': 'whole-heart/whole-ventricle bulk tissue, not cardiomyocyte-derived'},
        {'candidate': 'Reichart-CM single NCCM donor',
         'reason': 'ADR 0066 does not use the NCCM donor; it is excluded from the tracked '
                   'export and appears in no dataset of this category'},
    ],
}

DATASETS.extend([
    {
        'dataset_id': 'reichart-cm__dcm_genotype_positive',
        'category': 'genetic_cardiomyopathy',
        'source_item': 'disease_data/Reichart-CM',
        'local_deg': 'local/deg/reichart-cm_genetic_deg.csv',
        'local_summary': 'local/deg/reichart-cm_genetic_deg.json',
        'arm_or_class': 'genotype-positive DCM donors (44) vs the 18 normal donors',
        'species': 'Homo sapiens',
        'model_system': 'human LV/RV/apex/IVS sc+snRNA CM-only cellxgene object pseudobulk',
        'evidence_type': 'cm_pseudobulk',
        'tier': 1,
        'genotype_status': 'pathogenic',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'reichart-cm:normal',
        'n_case': 44,
        'n_control': 18,
        'contrasts': ['dcm_genotype_positive-normal'],
        'method': {
            'deposition': 'integer counts (per-donor CM pseudobulk)',
            'model': 'limma-voom (voom library-size normalisation, moderated t)',
            'design': ('one joint fit ~0 + group over the 78 donors exported by ticket #187 '
                       '(18 normal / 8 genotype-negative DCM / 44 genotype-positive DCM / '
                       '8 ARVC); this dataset selects the genotype-positive DCM vs normal '
                       'contrast'),
            'derivation': ('disease_data/Reichart-CM/script/build_genetic_deg.py runs the one '
                           'live implementation over ticket #187\'s tracked pseudobulk + '
                           'genotype export, read-only; the pseudobulk is not re-derived'),
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'reichart-cm__arvc',
        'category': 'genetic_cardiomyopathy',
        'source_item': 'disease_data/Reichart-CM',
        'local_deg': 'local/deg/reichart-cm_genetic_deg.csv',
        'local_summary': 'local/deg/reichart-cm_genetic_deg.json',
        'arm_or_class': 'ARVC donors (8) vs the 18 normal donors',
        'species': 'Homo sapiens',
        'model_system': 'human LV/RV/apex/IVS sc+snRNA CM-only cellxgene object pseudobulk',
        'evidence_type': 'cm_pseudobulk',
        'tier': 1,
        'genotype_status': 'pathogenic',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'reichart-cm:normal',
        'n_case': 8,
        'n_control': 18,
        'contrasts': ['arvc-normal'],
        'method': {
            'deposition': 'integer counts (per-donor CM pseudobulk)',
            'model': 'limma-voom (voom library-size normalisation, moderated t)',
            'design': ('one joint fit ~0 + group over the 78 donors exported by ticket #187 '
                       '(18 normal / 8 genotype-negative DCM / 44 genotype-positive DCM / '
                       '8 ARVC); this dataset selects the ARVC vs normal contrast'),
            'derivation': ('disease_data/Reichart-CM/script/build_genetic_deg.py runs the one '
                           'live implementation over ticket #187\'s tracked pseudobulk + '
                           'genotype export, read-only; the pseudobulk is not re-derived'),
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'gse224129__mybpc3_hcm',
        'category': 'genetic_cardiomyopathy',
        'source_item': 'disease_data/GSE224129',
        'local_deg': 'local/deg/GSE224129_deg.csv',
        'local_summary': 'local/deg/GSE224129_deg.json',
        'arm_or_class': 'MYBPC3 +/- and -/- engineered cardiac tissue vs the isogenic +/+ arm',
        'species': 'Homo sapiens',
        'model_system': 'human hiPSC-CM engineered cardiac tissue (isogenic MYBPC3 +/- and -/-)',
        'evidence_type': 'bulk',
        'tier': 2,
        'genotype_status': 'pathogenic',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse224129:isogenic',
        'n_case': 8,
        'n_control': 4,
        'contrasts': ['MYBPC3_het-MYBPC3_wt', 'MYBPC3_hom-MYBPC3_wt'],
        'isogenic_exception': {
            'applies': True,
            'control_arm': 'cMyBP-C+/+ ECT, 4 libraries (isogenic +/+ line, DF19-9-11T.H)',
            'why': ('ADR 0066 Isogenic exception: a single-genotype isogenic HCM model is '
                    'kept as its own evidence type even though the human CM-pseudobulk '
                    'Reichart cohort supplies this category, because a mixed-etiology clinical '
                    'cohort and a causal single-lesion model answer different questions'),
        },
        'strata': {
            'decision': ('two mutant genotype arms (heterozygous, homozygous) vs one isogenic '
                         '+/+ control arm; both contrasts from one joint fit, collapsed per '
                         'gene by maximum effect'),
            'pooled': True,
            'blocking_covariates': [],
        },
        'method': {
            'deposition': 'RSEM gene-level expected counts (non-integer, count scale)',
            'integerisation': ('bases.bulk_deg.integerise_counts: round to nearest, ties to '
                               'even; max |x - round(x)| 0.5; 52,572 non-integer values; '
                               'library sums 24.1-30.7M before rounding; 119 genes whose '
                               'rounded row is all-zero recorded and dropped by the count '
                               'model, never silently filtered'),
            'model': 'limma-voom (voom library-size normalisation, moderated t)',
            'design': ('one joint fit ~0 + group over all 12 ECT libraries; the +\u2044- and '
                       '-/- contrasts share the isogenic +/+ control'),
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'gse234141__ttntv_dcm',
        'category': 'genetic_cardiomyopathy',
        'source_item': 'disease_data/GSE234141',
        'local_deg': 'local/deg/GSE234141_deg.csv',
        'local_summary': 'local/deg/GSE234141_deg.json',
        'arm_or_class': ('TTN truncating-variant hiPSC-CM vs the isogenic WT/WT line across '
                         '2D/EHT x ventricular-like/atrial-like strata'),
        'species': 'Homo sapiens',
        'model_system': 'human hiPSC-CM 2D monolayer and engineered heart tissue (isogenic TTNtv)',
        'evidence_type': 'bulk',
        'tier': 2,
        'genotype_status': 'pathogenic',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse234141:isogenic',
        'n_case': 12,
        'n_control': 12,
        'contrasts': ['TTNtv-WT.2d_ventricular', 'TTNtv-WT.2d_atrial',
                      'TTNtv-WT.eht_ventricular', 'TTNtv-WT.eht_atrial'],
        'isogenic_exception': {
            'applies': True,
            'control_arm': 'WT/WT isogenic line, 12 libraries (3 per format x cell-type stratum)',
            'why': ('ADR 0066 Isogenic exception: the TTNtv model is a causal single-lesion '
                    'isogenic pair, kept as its own evidence type beside the human '
                    'CM-pseudobulk cohort'),
        },
        'strata': {
            'decision': ('four paired strata (format 2D vs EHT) x (cell type ventricular-like '
                         'vs atrial-like), each TTNtv vs its isogenic WT/WT, collapsed per '
                         'gene by maximum effect into one published dataset'),
            'pooled': True,
            'blocking_covariates': [],
            'why': ('format and cell type are culture/differentiation strata of one isogenic '
                    'pair, not independent cohorts; each stratum keeps its own paired contrast '
                    'inside one joint fit (shared dispersion) and the four collapse by maximum '
                    'effect (ticket-01 rule for a multi-arm source)'),
        },
        'method': {
            'deposition': 'integer counts',
            'model': 'limma-voom (voom library-size normalisation, moderated t)',
            'design': ('one joint fit ~0 + group over all 24 libraries with 8 composite '
                       'format x cell-type x genotype arms; four paired TTNtv-vs-WT contrasts'),
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'gse304575__lmna_h222p_dcm',
        'category': 'genetic_cardiomyopathy',
        'source_item': 'disease_data/GSE304575',
        'local_deg': 'local/deg/GSE304575_deg.csv',
        'local_summary': 'local/deg/GSE304575_deg.json',
        'arm_or_class': 'LMNA p.His222Pro heterozygous hiPSC-CM vs the CRISPR-corrected isogenic control',
        'species': 'Homo sapiens',
        'model_system': 'human hiPSC-CM (isogenic LMNA p.His222Pro vs CRISPR-corrected)',
        'evidence_type': 'bulk',
        'tier': 2,
        'genotype_status': 'pathogenic',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse304575:isogenic',
        'n_case': 6,
        'n_control': 6,
        'contrasts': ['LMNA_H222P-corrected'],
        'isogenic_exception': {
            'applies': True,
            'control_arm': 'Corrected Isogenic control, 6 libraries (CRISPR/Cas9-corrected line)',
            'why': ('ADR 0066 Isogenic exception: a causal single-variant isogenic DCM model '
                    'kept as its own evidence type beside the mixed-etiology human cohort'),
        },
        'strata': {
            'decision': 'single case/control pair; no stratum to pool',
            'pooled': False,
            'blocking_covariates': [],
            'why': ('one genotype arm against one isogenic control arm; the only batch value is '
                    'constant (6/15/2023) on all 12 libraries'),
        },
        'method': {
            'deposition': 'integer counts',
            'model': 'limma-voom (voom library-size normalisation, moderated t)',
            'design': 'one joint fit ~0 + group over all 12 libraries; single contrast',
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'gse305648__rbm20_dcm',
        'category': 'genetic_cardiomyopathy',
        'source_item': 'disease_data/GSE305648',
        'local_deg': 'local/deg/GSE305648_deg.csv',
        'local_summary': 'local/deg/GSE305648_deg.json',
        'arm_or_class': 'RBM20 KO / PL / RQ cardioids vs the isogenic WTC11 wildtype cardioids',
        'species': 'Homo sapiens',
        'model_system': 'human hiPSC cardioids (isogenic RBM20 KO / PL / RQ vs WTC11 wildtype)',
        'evidence_type': 'bulk',
        'tier': 2,
        'genotype_status': 'pathogenic',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse305648:isogenic',
        'n_case': 9,
        'n_control': 3,
        'contrasts': ['RBM20_KO-WT', 'RBM20_PL-WT', 'RBM20_RQ-WT'],
        'isogenic_exception': {
            'applies': True,
            'control_arm': 'RBM20 wildtype cardioids, 3 libraries (WTC11 background)',
            'why': ('ADR 0066 Isogenic exception: the three RBM20 alleles are causal '
                    'single-lesion isogenic models on one WTC11 background, kept as their own '
                    'evidence type'),
        },
        'strata': {
            'decision': ('three mutant genotype arms (KO, PL, RQ) vs one isogenic wildtype arm; '
                         'collapsed per gene by maximum effect'),
            'pooled': True,
            'blocking_covariates': ['batch'],
            'why': ('the mutant alleles are distinct RBM20 lesions, so each keeps its own '
                    'contrast inside one fit; the library batch (A/B/C) is rotated across the '
                    'arms and enters the design as a blocking covariate'),
        },
        'method': {
            'deposition': 'integer counts (Ensembl ids with version suffix)',
            'model': 'limma-voom (voom library-size normalisation, moderated t)',
            'design': ('one joint fit ~0 + group + batch over all 12 cardioid libraries; three '
                       'mutant-vs-WT contrasts'),
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'gse253223__pkp2_knockdown_arvc',
        'category': 'genetic_cardiomyopathy',
        'source_item': 'disease_data/GSE253223',
        'local_deg': 'local/deg/GSE253223_deg.csv',
        'local_summary': 'local/deg/GSE253223_deg.json',
        'arm_or_class': 'siPKP2 (PKP2 knockdown) vs siNeg in wild-type hiPSC-CM, days 2/4/6/8',
        'species': 'Homo sapiens',
        'model_system': 'human wild-type hiPSC-CM (acquired PKP2 siRNA knockdown)',
        'evidence_type': 'bulk',
        'tier': 2,
        'genotype_status': 'not_applicable',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse253223:isogenic',
        'n_case': 16,
        'n_control': 4,
        'contrasts': ['siPKP2-siNeg'],
        'isogenic_exception': {
            'applies': True,
            'control_arm': 'siNeg (non-targeting siRNA), 4 libraries (one per harvest day)',
            'why': ('ADR 0066 Isogenic exception is applied to the model system rather than to '
                    'a patient genotype: this is an acquired PKP2 knockdown on a single '
                    'wild-type line, so it is kept as its own in-vitro ARVC-loss evidence type '
                    'and is explicitly NOT a patient variant'),
        },
        'strata': {
            'decision': ('four harvest days (2/4/6/8) pooled in one fit with day as a blocking '
                         'covariate'),
            'pooled': True,
            'blocking_covariates': ['day'],
            'why': ('each day has one siNeg library and four siPKP2 libraries, so a per-day '
                    'contrast would rest on a single control; pooling with day as a covariate '
                    'uses all 20 libraries for one knockdown effect'),
        },
        'method': {
            'deposition': 'estimated counts (non-integer, count scale)',
            'integerisation': ('bases.bulk_deg.integerise_counts: round to nearest, ties to '
                               'even; max |x - round(x)| 0.5; 282,071 non-integer values; '
                               'library sums 25.2-41.9M before rounding; 4 genes whose rounded '
                               'row is all-zero recorded and dropped by the count model, never '
                               'silently filtered'),
            'model': 'limma-voom (voom library-size normalisation, moderated t)',
            'design': ('one joint fit ~0 + group + day over all 20 libraries; single '
                       'siPKP2-vs-siNeg contrast'),
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'gse194297__dmd_cardiomyopathy',
        'category': 'genetic_cardiomyopathy',
        'source_item': 'disease_data/GSE194297',
        'local_deg': 'local/deg/GSE194297_deg.csv',
        'local_summary': 'local/deg/GSE194297_deg.json',
        'arm_or_class': 'DMD point-mutation cardiac organoids vs the CRISPR-corrected isogenic line',
        'species': 'Homo sapiens',
        'model_system': 'human hiPSC cardiac organoids (DMD point mutation vs corrected isogenic)',
        'evidence_type': 'bulk',
        'tier': 2,
        'genotype_status': 'pathogenic',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse194297:isogenic',
        'n_case': 3,
        'n_control': 3,
        'contrasts': ['DMD-DMD_corrected'],
        'isogenic_exception': {
            'applies': True,
            'control_arm': 'CRISPR/Cas9-corrected isogenic line, 3 day-56 cardiac organoids',
            'why': ('ADR 0066 Isogenic exception: the DMD mutation was corrected in the same '
                    'patient-derived line, so the model is its own causal isogenic evidence '
                    'type'),
        },
        'strata': {
            'decision': 'single case/control pair at one time point; no stratum to pool',
            'pooled': False,
            'blocking_covariates': [],
            'why': 'one mutant line against its corrected isogenic line (day 56)',
        },
        'method': {
            'deposition': 'integer counts',
            'model': 'limma-voom (voom library-size normalisation, moderated t)',
            'design': 'one joint fit ~0 + group over all 6 organoid libraries; single contrast',
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
])
