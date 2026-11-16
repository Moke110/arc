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

The registry declares nHF, ICM, AF and genetic_cardiomyopathy; the unified
TOX category is declared by its independent pre-DEG admission ledger. ADR 0073
collects every qualifying declared source, with dated finite reassessments below.
"""
from __future__ import annotations

CATEGORIES = {
    'nHF': {
        'directory': 'nhf',
        'tier': 1,
        'evidence_type': 'cm_pseudobulk',
        'shared_control_group': ['gse270788:donor', 'reichart-cm:normal'],
        'rejected': [],  # Current finite reassessment is declared below.
    },
    'ICM': {
        'directory': 'icm',
        'tier': 1,
        'evidence_type': 'cm_pseudobulk',
        'shared_control_group': ['gse270788:donor'],
        'rejected': [
            {'candidate': 'GSE298023',
             'reason': ('eligible human LV CM source, but its author cell labels exist '
                        'only inside an unopened 4.0 GB .rds, so no author CM gate was '
                        'available and a gate invented from markers is not permitted; '
                        'the 1.8 GB per-sample 10x RAW.tar carries counts only')},
            {'candidate': 'GSE145154',
             'reason': 'no author cell labels and region-heterogeneous; re-clustering '
                       'not permitted'},
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
            'ADR 0066 Isogenic exception: the declared in-vitro hiPSC-CM models carry a disease '
            'genotype with an isogenic control and are kept as their own evidence type inside '
            'this category even though a human CM-pseudobulk cohort (Reichart-CM) exists for '
            'the same category: a mixed-etiology clinical cohort and a single-genotype '
            'isogenic model answer different questions, and collapsing them loses the causal '
            'signature. The category therefore has one tier-1 human CM-pseudobulk pair '
            '(Reichart genotype-positive DCM + ARVC) and tier-2 isogenic bulk datasets; '
            '`tier` records the winning tier and each dataset record carries its own '
            '`evidence_type` and `tier`.'),
    },
    'rejected': [],  # Current finite reassessment is declared below.
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

# --- ticket #203: nHF tier-2 / tier-3 sources --------------------------------
# GSE35108 is the tier-2 (patient-derived iPSC-CM) source and GSE165303 the tier-3
# (whole left-ventricle bulk) source ADR 0071 gives nHF. Both are admitted only
# after their series design was verified against the deposit (see each item's
# author_log.json): GSE35108 is 8 cardiomyocyte libraries in a 2x2 genotype x
# Serca2a design, not the series' 17 samples; GSE165303 is 101 RNA-seq libraries
# = 50 DCM + 51 non-failing, not "101 DCM samples".
DATASETS.extend([
    {
        'dataset_id': 'gse35108__dcm_ipsc_cm',
        'category': 'nHF',
        'source_item': 'disease_data/GSE35108',
        'local_deg': 'local/deg/GSE35108_deg.csv',
        'local_summary': 'local/deg/GSE35108_deg.json',
        'arm_or_class': ('patient-derived DCM iPSC-CM (4 libraries) vs normal iPSC-CM '
                         '(4 libraries) in a 2x2 genotype x Serca2a design'),
        'species': 'Homo sapiens',
        'model_system': ('human patient-derived iPSC-CM monolayer (DCM vs normal), '
                         'one-colour Affymetrix GPL6244 RMA log2 array'),
        'evidence_type': 'bulk',
        'tier': 2,
        'genotype_status': 'not_declared',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse35108:normal_ipsc_cm',
        'n_case': 4,
        'n_control': 4,
        'contrasts': ['DCM_control-Normal_control', 'DCM_Serca2a-Normal_Serca2a'],
        'non_genetic_model_basis': (
            'patient-derived familial-DCM iPSC-CM lines compared with normal iPSC-CM '
            'lines; no mutation was introduced and no isogenic control is deposited, so '
            'this is the clinical (non-mutation-modelled) iPSC-CM form ADR 0071 admits '
            'to a clinical category at tier 2'),
        'strata': {
            'decision': ('two treatment strata of the same genotype contrast (untreated, '
                         'and after adenoviral Serca2a gene therapy), each 2 vs 2; '
                         'collapsed per gene by maximum effect'),
            'pooled': True,
            'blocking_covariates': [],
            'why': ('the gene therapy is applied symmetrically to both genotypes, so each '
                    'stratum carries a valid DCM-vs-normal contrast; the untreated stratum '
                    'is the primary disease contrast and its DEG count is recorded as a '
                    'sensitivity in the item summary'),
        },
        'method': {
            'deposition': 'RMA log2 expression as deposited (one-colour array)',
            'collapse': ('the max-mean-expression GPL6244 probe represents each gene '
                         '(19,036 genes from 20,052 single-symbol probes; multi-symbol '
                         'and unsymbolised probes unused)'),
            'model': 'limma on the deposited log2 matrix (moderated t)',
            'design': ('one joint fit ~0 + group over the 8 cardiomyocyte libraries; two '
                       'DCM-vs-normal contrasts, one per treatment stratum'),
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'gse165303__dcm_lv_bulk',
        'category': 'nHF',
        'source_item': 'disease_data/GSE165303',
        'local_deg': 'local/deg/GSE165303_deg.csv',
        'local_summary': 'local/deg/GSE165303_deg.json',
        'arm_or_class': 'DCM left-ventricle tissue (50 donors) vs non-failing LV (51 donors)',
        'species': 'Homo sapiens',
        'model_system': ('human left-ventricle bulk tissue (whole tissue, not '
                         'cardiomyocyte-specific), bulk RNA-seq'),
        'evidence_type': 'bulk',
        'tier': 3,
        'genotype_status': 'not_declared',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse165303:non_failing_lv',
        'n_case': 50,
        'n_control': 51,
        'contrasts': ['DCM-NF'],
        'model_basis': (
            'whole left-ventricle tissue bulk: not cardiomyocyte-specific, so it is a '
            'tier-3 source; ADR 0073 collects it regardless of higher-tier yield'),
        'method': {
            'deposition': ('log2(normalised counts) as deposited, with a fixed floor at '
                           '3.821687527 (35.2% of cells); not raw counts'),
            'model': 'limma on the deposited log2 matrix (moderated t)',
            'design': ('one joint fit ~0 + group over all 101 RNA-seq libraries; no '
                       'clinical covariate is deposited'),
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
])

# --- ticket #202: ICM tier-1 Kuppe-MI (IZ vs CTRL) ---------------------------
# The cellxgene Kuppe MI all-snRNA object supplies ICM's human cardiomyocyte
# pseudobulk evidence at tier 1. ADR 0073 also collects qualifying lower tiers.
# The object is gated on its own CM label and
# the profiles are per (donor, zone), not per donor: an infarct heart contributes
# its ischaemic-zone CM nuclei as one profile, which is what makes the
# ischaemic-zone-vs-control-heart contrast possible. A 500-CM-nuclei profile floor
# removes the sparse donor-zone profiles, one joint ~0 + zone fit estimates
# dispersion over all five zones, and the repeated donor profiles are blocked with
# duplicateCorrelation(block = donor). Only the IZ-vs-CTRL contrast publishes a
# DEG; the border/remote/fibrotic-zone contrasts are fitted for dispersion and
# recorded as artifacts in the item summary. The deposit carries no per-donor
# time-since-infarct field, so acute-vs-chronic stays undeclared (see
# `time_window_note`), not asserted.
CATEGORIES['ICM']['shared_control_group'] = (
    list(CATEGORIES['ICM']['shared_control_group']) + ['kuppe-mi:normal'])
DATASETS.append({
    'dataset_id': 'kuppe-mi__iz',
    'category': 'ICM',
    'source_item': 'disease_data/Kuppe-MI',
    'local_deg': 'local/deg/kuppe-mi_deg.csv',
    'local_summary': 'local/deg/kuppe-mi_deg.json',
    'cells_per_donor': 'local/pseudobulk/kuppe-mi_profiles.csv',
    'arm_or_class': ('IZ (ischaemic zone) cardiomyocyte profiles of 4 infarct donors '
                     '(P2, P9, P10, P15) vs the 4 normal-heart CM profiles '
                     '(P1, P7, P8, P17), after the 500-CM-nuclei profile floor'),
    'species': 'Homo sapiens',
    'model_system': ("human left-ventricle 10x 3' v3 snRNA CM pseudobulk, one profile "
                     'per (donor, zone) of the infarcted heart'),
    'evidence_type': 'cm_pseudobulk',
    'tier': 1,
    'genotype_status': 'not_declared',
    'time_window_basis': None,
    'time_window_note': ('the deposit carries no per-donor time-since-infarct field, so the '
                         'deposit does not separate acute from chronic MI; the record is '
                         'undeclared, never asserted'),
    'exploratory': False,
    'shared_control_group': 'kuppe-mi:normal',
    'n_case': 4,
    'n_control': 4,
    'contrasts': ['IZ-CTRL'],
    'method': {
        'deposition': 'integer counts (cellxgene raw.X)',
        'model': ('limma-voom (voom library-size normalisation, moderated t) with donor '
                  'blocking: duplicateCorrelation(block = donor) on the joint design'),
        'design': ('one joint fit ~0 + zone over the 18 (donor, zone) profiles that pass the '
                   '500-CM-nuclei floor (CTRL 4 / IZ 4 / BZ 3 / RZ 5 / FZ 2); this dataset '
                   'selects the IZ vs CTRL contrast, the other three zone contrasts are '
                   'artifact rows of the same fit'),
        'test': 'moderated t',
        'padj': 'Benjamini-Hochberg',
    },
})

# --- ticket #204: AF returns (tier-1 CM-specific sources + tier-3 bulk) ------
# ADR 0071 decision 7 re-admits AF after ADR 0069/0070 removed it: its CM-specific
# sources (GSE255612, GSE138252, GSE238242) are tier 1; the whole-atrial bulk
# sources enter at tier 3 independently of their yield under ADR 0073.
# Every path below is item-relative; the
# driver composes the repo-relative path, so no gitignored `local/` literal points
# at the checkout. GSE255612's record is restored unchanged from before ADR 0069:
# its DEG item is not rebuilt by this ticket.
CATEGORIES['AF'] = {
    'directory': 'af',
    'tier': 1,
    'evidence_type': 'cm_pseudobulk',
    'shared_control_group': ['gse255612:donor', 'gse138252:control_la',
                             'gse238242:sr', 'gse269882:sinus_rhythm'],
    'rejected': [
        {'candidate': 'GSE271229',
         'accession': 'GSE271229',
         'status': 'not_admitted',
         'reason': ('left-atrial snRNA-seq (9 donors: 3 no-AF / 3 paroxysmal / 3 '
                    'persistent) whose deposit ships NO author cell-type annotation — '
                    'only per-sample CellRanger barcodes/features/matrix triples and an '
                    'empty series-matrix data table. A CM gate would have to be invented '
                    'from markers, which is forbidden, so the DEG step is not testable. '
                    'This is a different accession from GSE269882, not an alias of it.')},
        {'candidate': 'dbGaP phs003644',
         'status': 'not_downloadable',
         'reason': 'controlled access; not fetched'},
    ],
    'historical_tier3_owned': [
        {'candidate': 'GSE128188',
         'item': 'disease_data/GSE128188',
         'reason': ('already owned human atrial bulk item (AF 5 vs 5 SR, paired LA+RA). '
                    'ADR 0071 admits whole-tissue bulk only at tier 3, so this is a '
                    'candidate for the #206 roster; it is not a rejection and it is not '
                    'declared here because #204 builds only the three new items.')},
    ],
    'historical_tiering': {
        'tier_1': {
            'sources': ['gse255612__af', 'gse138252__af', 'gse238242__af'],
            'measured_union_v2_deg': 1,
            'measured_union_genes': ['SLC9A7'],
            'note': ('measured with the one live DEG implementation under ADR 0071\'s '
                     'definition (padj < 0.05 and |log2FC| >= 2): GSE255612 0 (its '
                     'existing v1-contract table), GSE138252 1 (SLC9A7), GSE238242 0. '
                     'Reconciled by #211: ADR 0071 decision 7 now states this '
                     'live-worker single-gene union (SLC9A7), and the pre-existing '
                     'audit\'s 2-gene figure -- its limma-voom + edgeR TMM result, TMM '
                     'being the only difference, since the pinned arc R env ships limma '
                     'but not edgeR and the live worker therefore uses voom\'s own '
                     'library-size normalisation -- is recorded there and in the '
                     'GSE138252 item summary as a sensitivity, never as the published '
                     'number.'),
        },
        'tier_3': {
            'sources': ['gse269882__af'],
            'measured_union_v2_deg': 0,
            'note': ('GSE269882 is whole-tissue right-atrial-appendage bulk and is '
                     'admitted only at tier 3; its three staged AF-vs-SR contrasts and '
                     'the pooled sensitivity all give 0 v2 DEGs under the live worker. '
                     'GSE128188 (owned, undeclared here) is the other tier-3 candidate.'),
        },
        'historical_stop_decision': ('tier 1 union 1 < 2,000, so the collection descends to tier 3; '
                          'the tier-3 union is 0, so AF\'s total axis-resolved DEG union '
                          'is the tier-1 gene set (1 gene under the live worker). The '
                          'category is published with the tier-3 evidence type recorded '
                          'on GSE269882.'),
    },
}

DATASETS.extend([
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
    {
        'dataset_id': 'gse138252__af',
        'category': 'AF',
        'source_item': 'disease_data/GSE138252',
        'local_deg': 'local/deg/GSE138252_deg.csv',
        'local_summary': 'local/deg/GSE138252_deg.json',
        'arm_or_class': ('left-atrium AF+HF (5 libraries) vs non-failing control '
                         '(5 libraries), 2020 sequencing batch only'),
        'species': 'Homo sapiens',
        'model_system': ('human atrial PCM1+ cardiomyocyte-nuclei bulk RNA-seq '
                         '(CM-nuclear bulk, not whole tissue)'),
        'evidence_type': 'bulk',
        'tier': 1,
        'genotype_status': 'not_declared',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse138252:control_la',
        'n_case': 5,
        'n_control': 5,
        'contrasts': ['AF-CTL'],
        'model_basis': ('cardiomyocyte-derived (PCM1+ FANS-sorted nuclei), so tier 1 '
                        'under ADR 0071 even though the model is a bulk library; the '
                        'AF+HF-vs-non-failing confound is recorded in the item summary'),
        'strata': {
            'decision': ('left atrium, 2020 batch only: the AF arm exists only in that '
                         'batch, so a batch-matched control arm is the primary contrast '
                         'and a pooled-controls fit is a declared sensitivity'),
            'pooled': False,
            'blocking_covariates': [],
            'why': ('region is the dominant variance axis in this series and pooling '
                    'chambers loses power; the RA arm is sex-confounded'),
        },
        'method': {
            'deposition': 'raw integer gene-body counts',
            'model': "limma-voom (voom's own library-size normalisation, moderated t)",
            'design': ('one joint fit ~0 + group over the 5 AF + 5 control LA libraries '
                       'of the 2020 batch (edgeR filterByExpr rule)'),
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'gse238242__af',
        'category': 'AF',
        'source_item': 'disease_data/GSE238242',
        'local_deg': 'local/deg/GSE238242_deg.csv',
        'local_summary': 'local/deg/GSE238242_deg.json',
        'cells_per_donor': 'local/pseudobulk/GSE238242_profiles.csv',
        'arm_or_class': ('persistent-AF donors (3) vs sinus-rhythm donors (4), '
                         'left-atrial-appendage CM pseudobulk'),
        'species': 'Homo sapiens',
        'model_system': ('human left-atrial-appendage 10x multiome snRNA arm, author '
                         'CM label, per-donor CM count pseudobulk'),
        'evidence_type': 'cm_pseudobulk',
        'tier': 1,
        'genotype_status': 'not_declared',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse238242:sr',
        'n_case': 3,
        'n_control': 4,
        'contrasts': ['AF-SR'],
        'method': {
            'deposition': 'raw integer counts (dense TSV, genes x nuclei)',
            'model': ("limma-voom (voom's own library-size normalisation, moderated t) "
                      'on the per-donor sum-of-CM-counts pseudobulk'),
            'design': ('one joint fit ~0 + group over the 7 donor pseudobulk profiles '
                       '(4 SR / 3 persistent AF)'),
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'gse269882__af',
        'category': 'AF',
        'source_item': 'disease_data/GSE269882',
        'local_deg': 'local/deg/GSE269882_deg.csv',
        'local_summary': 'local/deg/GSE269882_deg.json',
        'arm_or_class': ('the three staged AF arms (persistent 5, paroxysmal 5, '
                         'post-operative 8) vs the 9 sinus-rhythm libraries'),
        'species': 'Homo sapiens',
        'model_system': ('human right-atrial-appendage whole-tissue bulk RNA-seq '
                         '(not cardiomyocyte-specific)'),
        'evidence_type': 'bulk',
        'tier': 3,
        'genotype_status': 'not_declared',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse269882:sinus_rhythm',
        'n_case': 18,
        'n_control': 9,
        'contrasts': ['persAF-SR', 'pAF-SR', 'poAF-SR'],
        'model_basis': ('whole right-atrial-appendage tissue: not cardiomyocyte-specific, '
                        'so it is a tier-3 source; ADR 0073 collects every qualifying tier, '
                        'superseding the old count-based stopping condition'),
        'strata': {
            'decision': ('the three diagnosis arms are declared as separate contrasts '
                         'against sinus rhythm inside one joint fit; the staging is '
                         'carried, never averaged into one case group, and the '
                         'post-operative arm is its own declared arm'),
            'pooled': True,
            'blocking_covariates': [],
            'why': ('a pooled all-AF-vs-SR fit is measured as a declared sensitivity '
                    'and does not change the conclusion (0 v2 DEGs)'),
        },
        'method': {
            'deposition': 'raw integer counts (featureCounts, GRCh38)',
            'model': "limma-voom (voom's own library-size normalisation, moderated t)",
            'design': ('one joint fit ~0 + group over all 27 libraries with the four '
                       'diagnosis arms; three staged AF-vs-SR contrasts'),
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
])

# The ladder is derived from each record's tier; the public tier_ledger records
# the complete cumulative axis-resolved union. Source admission does not depend
# on DEG yield. Historical measured ADR0071 records above remain labelled history.
DATASETS.extend([
    {
        'dataset_id': 'gse116250__hf_lv_bulk',
        'category': 'nHF',
        'source_item': 'disease_data/GSE116250',
        'local_deg': 'local/deg/GSE116250_deg.csv',
        'local_summary': 'local/deg/GSE116250_deg.json',
        # The #179-era HF item deposits the shared DEG table with its own column
        # names; the canonical contract is applied by `bases.disease_deg.read_local_deg`.
        'local_deg_columns': {'log2FC': 'log2fc', 'mean_expr': 'base_mean'},
        'arm_or_class': ('heart failure (50: 37 dilated + 13 ischaemic cardiomyopathy) vs '
                         'non-failing left ventricle (14 donors), three fitted arm contrasts'),
        'species': 'Homo sapiens',
        'model_system': ('human left-ventricular bulk tissue (whole tissue, not '
                         'cardiomyocyte-specific), bulk RNA-seq'),
        'evidence_type': 'bulk',
        'tier': 3,
        'genotype_status': 'not_declared',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse116250:nf',
        'n_case': 50,
        'n_control': 14,
        'contrasts': ['HF_vs_NF', 'DCM_vs_NF', 'ICM_vs_NF'],
        'model_basis': ('whole left-ventricular tissue bulk: not cardiomyocyte-specific, '
                        'so it is a tier-3 source; ADR 0073 collects every qualifying tier '),
        'strata': {
            'decision': ('one joint fit ~0 + group over the 64 libraries (14 NF / 37 DCM / '
                         '13 ICM) publishes the primary HF_vs_NF contrast and the two arm '
                         'contrasts as contrast rows of the same fit; the multi-contrast '
                         'collapse keeps the maximum-effect contrast per gene and the DEG '
                         'flag is pass-any-contrast'),
            'pooled': True,
            'blocking_covariates': [],
            'why': ('the three labels are arms of one model, not independent sources; the '
                    'tier-3 contribution is their pass-any-contrast v2 DEG union (43 '
                    'axis DEGs measured, 144 genome-wide). The 13-donor ICM arm is an '
                    'end-stage patient ICM comparison; ADR 0073 also records it once in '
                    'the ICM category at tier 3, with the shared NF controls explicit'),
        },
        'method': {
            'deposition': ('RPKM as deposited (Cufflinks / GRCh37.82), analysed as '
                           'log2(RPKM + 1); the submitters removed mitochondrial genes'),
            'model': ('limma on the deposited log2(RPKM + 1) matrix; Ensembl annotations '
                      'collapsed to HGNC symbols by maximum mean expression'),
            'design': ('one joint fit ~0 + group over all 64 libraries with group in '
                       '{NF, DCM, ICM}; three contrasts (HF_vs_NF primary)'),
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
    {
        'dataset_id': 'gse141910__hf_lv_bulk',
        'category': 'nHF',
        'source_item': 'disease_data/GSE141910',
        'local_deg': 'local/deg/GSE141910_deg.csv',
        'local_summary': 'local/deg/GSE141910_deg.json',
        'local_deg_columns': {'log2FC': 'log2fc', 'mean_expr': 'base_mean'},
        'arm_or_class': 'heart failure (200: 166 DCM / 28 HCM / 6 PPCM) vs non-failing LV (166 donors)',
        'species': 'Homo sapiens',
        'model_system': ('human left-ventricular free-wall bulk tissue (whole tissue, not '
                         'cardiomyocyte-specific), bulk RNA-seq, MAGNet consortium'),
        'evidence_type': 'bulk',
        'tier': 3,
        'genotype_status': 'not_declared',
        'time_window_basis': None,
        'exploratory': False,
        'shared_control_group': 'gse141910:nf',
        'n_case': 200,
        'n_control': 166,
        'contrasts': ['HF_vs_NF'],
        'model_basis': ('whole left-ventricular free-wall bulk: not cardiomyocyte-specific, '
                        'so it is a tier-3 source; all qualifying tiers are collected (ADR 0073)'),
        'strata': {
            'decision': ('single pooled HF-vs-NF contrast; the deposited etiology labels '
                         '(DCM/HCM/PPCM) are recorded but not modelled'),
            'pooled': True,
            'blocking_covariates': [],
            'why': ('the item is the programme\'s HF primary anchor; the 28 hereditary HCM '
                    'donors inside the 200 HF samples are a recorded limitation of the '
                    'pooled contrast, not a genetic-cardiomyopathy source (that category '
                    'uses the genotype-declared cohorts)'),
        },
        'method': {
            'deposition': ('the submitters\' log2-scale voom-transformed abundance as '
                           'deposited (no raw counts are deposited)'),
            'model': ('limma on the deposited log2 abundance; HGNC symbols collapsed by '
                      'maximum mean expression with the retired Ensembl v67 ids kept as '
                      'an identifier fallback'),
            'design': 'one joint fit ~0 + group over all 366 libraries; one contrast',
            'test': 'moderated t',
            'padj': 'Benjamini-Hochberg',
        },
    },
])


# ADR 0073 complete collection reopens the owned, formerly deferred AF source.
DATASETS.append({
    'dataset_id': 'gse128188__af', 'category': 'AF',
    'source_item': 'disease_data/GSE128188',
    'local_deg': 'local/deg/GSE128188_deg.csv',
    'local_summary': 'local/deg/GSE128188_deg.json',
    'local_deg_columns': {'log2FC': 'log2fc', 'baseMean': 'base_mean'},
    'arm_or_class': 'AF vs sinus rhythm; paired LA+RA counts pooled per donor',
    'species': 'Homo sapiens', 'model_system': 'human paired atrial appendage bulk RNA-seq',
    'evidence_type': 'bulk', 'tier': 3, 'genotype_status': 'not_declared',
    'time_window_basis': None, 'exploratory': False,
    'shared_control_group': 'gse128188:sr', 'n_case': 5, 'n_control': 5,
    'contrasts': ['AF_vs_SR'],
    'model_basis': 'human atrial whole-tissue bulk; all qualifying tier-3 evidence collected (ADR 0073)',
    'method': {'deposition': 'integer counts', 'model': 'donor-pooled DESeq2 ~ group',
               'test': 'Wald; unshrunken log2FC', 'padj': 'Benjamini-Hochberg'},
    'strata': {'decision': 'sum each donor paired LA+RA counts before disease comparison',
               'pooled': True, 'blocking_covariates': [],
               'why': 'AF status is donor-level; donor fixed effects are confounded with disease'},
})

CATEGORIES['AF']['source_review'] = {
    'reviewed_on': '2026-11-16', 'policy': 'all_qualifying_declared_sources',
    'scope': 'Existing pinned source survey, rejected/screened records and owned external items; not a claim of exhaustive public discovery',
    'records': [
        {'candidate': 'GSE128188', 'disposition': 'admitted', 'tier': 3,
         'reason': 'Previously deferred owned human atrial bulk; 5 AF/5 SR donors, paired appendages pooled before DESeq2; no quality exclusion applies',
         'source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128188',
         'provenance': 'disease_data/GSE128188/author_log.json'},
        {'candidate': 'GSE271229', 'disposition': 'not_admitted',
         'reason': 'Reviewed author-label requirement retained; unannotated cell matrices do not supply a CM gate'},
        {'candidate': 'dbGaP phs003644', 'disposition': 'not_downloadable',
         'reason': 'Controlled access remains outside fetched evidence'},
    ]}

# The already pinned patient ICM arm qualifies independently within ICM at tier 3.
# Its shared NF controls with the nHF source remain explicit; cross-disease
# independence is not claimed or corrected by the fusion.
_icm_bulk = dict(next(row for row in DATASETS if row['dataset_id'] == 'gse116250__hf_lv_bulk'))
_icm_bulk.update(dataset_id='gse116250__icm_lv_bulk', category='ICM', n_case=13,
                 contrasts=['ICM_vs_NF'],
                 arm_or_class='13 patient ICM vs 14 non-failing LV donors',
                 model_basis='human ICM left-ventricular whole-tissue bulk, tier 3; complete qualifying collection under ADR 0073')
_icm_bulk['strata'] = dict(decision='ICM_vs_NF arm of the existing joint three-group limma fit',
                          pooled=False, blocking_covariates=[],
                          why='same dataset contributes once within ICM; NF controls shared with its nHF record')
DATASETS.append(_icm_bulk)
CATEGORIES['ICM']['source_review'] = {
    'reviewed_on': '2026-11-16', 'policy': 'all_qualifying_declared_sources',
    'records': [{'candidate': 'GSE116250 ICM arm', 'disposition': 'admitted', 'tier': 3,
                 'reason': 'Actual patient ICM vs NF human LV bulk; historical tier1 sufficiency is not an admission exclusion',
                 'source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116250',
                 'provenance': 'disease_data/GSE116250/author_log.json'}]}

DATASETS.append({
    'dataset_id': 'gse183852__dcm_cm', 'category': 'nHF',
    'source_item': 'disease_data/GSE183852',
    'local_deg': 'local/deg/GSE183852_deg.csv',
    'local_summary': 'local/deg/GSE183852_deg.json',
    'arm_or_class': 'author Cardiomyocytes DCM vs Donor DESeq2 table',
    'species': 'Homo sapiens', 'model_system': 'human LV cardiomyocyte snRNA-seq',
    'evidence_type': 'cm_pseudobulk', 'tier': 1, 'genotype_status': 'not_declared',
    'time_window_basis': None, 'exploratory': False,
    'shared_control_group': 'gse183852:donor', 'n_case': 13, 'n_control': 25,
    'contrasts': ['DCM-Donor'],
    'model_basis': 'Full author cell-type-specific DE table retained; 13/25 are snRNA libraries, not asserted unique donors; no raw-count refit or per-cell annotation',
    'method': {'deposition': 'full 26946-row author Cardiomyocytes workbook sheet',
               'model': 'author DESeq2; statistics preserved',
               'test': 'author differential expression statistics',
               'padj': 'author adjusted p values; missing values remain missing'},
})
CATEGORIES['nHF']['rejected'] = [r for r in CATEGORIES['nHF']['rejected'] if r['candidate'] != 'GSE183852']
CATEGORIES['nHF']['source_review'] = {
    'reviewed_on': '2026-11-16', 'policy': 'all_qualifying_declared_sources',
    'records': [{'candidate': 'GSE183852', 'disposition': 'admitted', 'tier': 1,
                 'reason': 'Full author per-gene CM DE table permits exact threshold/rank replay without Robj extraction; low yield is not exclusion',
                 'source': 'https://doi.org/10.1038/s44161-022-00028-6',
                 'provenance': 'disease_data/GSE183852/author_log.json'}]}

# Finite ADR 0073 re-review: count sufficiency and age alone do not exclude a source.
CATEGORIES['nHF']['rejected'] = [
    {'candidate': 'cellxgene myocarditis CM object',
     'reason': 'Shares normal donors with Reichart; not an independent additional source'},
    {'candidate': 'GSE282618',
     'reason': 'HHD exploratory 3-vs-2 subject design, not a declared nHF disease cohort'},
    {'candidate': 'GSE185100', 'status': 'author_cm_input_unavailable',
     'reason': 'GEO supplies count matrices but no author per-cell CM assignment. The publication-linked analysis repository references local integrated objects and full CM test tables that are not deposited there. The portal selected/rounded export does not provide a complete unrounded test table. No marker-derived gate or reclustering is substituted; age and low estimated yield are not exclusion grounds.',
     'reentry': 'Versioned author CM barcode assignments with counts, or complete author CM DCM-vs-ND logFC/padj statistics'}]
CATEGORIES['nHF']['source_review']['records'].append({
    'candidate': 'GSE185100', 'disposition': 'not_admitted_author_cm_input_unavailable',
    'reason': CATEGORIES['nHF']['rejected'][-1]['reason'],
    'reentry': CATEGORIES['nHF']['rejected'][-1]['reentry'],
    'sources': ['https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185100',
                'https://doi.org/10.1161/CIRCULATIONAHA.121.057763',
                'https://github.com/neda-mehdiabadi/Fetal-Gene-Program-snRNAseq',
                'https://github.com/neda-mehdiabadi/HeartExplorer'],
    'uncertainty': 'Author objects may be obtainable elsewhere; current finite public inventory does not expose the required inputs. Portal certificate chain fails default verification.'})

# Reopened isogenic sources after finite primary-input review (ADR 0073).
DATASETS.extend([{'dataset_id': 'gse248187__dsp_eht',
  'category': 'genetic_cardiomyopathy',
  'source_item': 'disease_data/GSE248187',
  'local_deg': 'local/deg/GSE248187_deg.csv',
  'local_summary': 'local/deg/GSE248187_deg.json',
  'arm_or_class': 'Human DSP knockout EHT versus isogenic control',
  'species': 'Homo sapiens',
  'model_system': 'Human DSP knockout EHT versus isogenic control',
  'evidence_type': 'bulk',
  'tier': 2,
  'genotype_status': 'pathogenic',
  'time_window_basis': None,
  'exploratory': False,
  'shared_control_group': 'gse248187:isogenic',
  'n_case': 3,
  'n_control': 3,
  'contrasts': ['DSP_KO-isogenic'],
  'isogenic_exception': {'applies': True,
                         'control_arm': 'isogenic control',
                         'why': 'Unchanged ADR0066 causal isogenic-model exception; reopened under '
                                'ADR0073 complete qualifying-source collection.'},
  'method': {'deposition': 'See pinned item summary and author_log for full input/test provenance',
             'model': 'Complete author edgeR test',
             'design': 'Human DSP knockout EHT versus isogenic control',
             'test': 'Author likelihood ratio',
             'padj': 'Benjamini-Hochberg'}},
 {'dataset_id': 'gse221945__tbx5_brs',
  'category': 'genetic_cardiomyopathy',
  'source_item': 'disease_data/GSE221945',
  'local_deg': 'local/deg/GSE221945_deg.csv',
  'local_summary': 'local/deg/GSE221945_deg.json',
  'arm_or_class': 'TBX5G145R Brugada hiPSC-CM versus corrected isogenic control',
  'species': 'Homo sapiens',
  'model_system': 'TBX5G145R Brugada hiPSC-CM versus corrected isogenic control',
  'evidence_type': 'bulk',
  'tier': 2,
  'genotype_status': 'pathogenic',
  'time_window_basis': None,
  'exploratory': False,
  'shared_control_group': 'gse221945:isogenic',
  'n_case': 4,
  'n_control': 4,
  'contrasts': ['TBX5_G145R-corrected'],
  'isogenic_exception': {'applies': True,
                         'control_arm': 'corrected isogenic control',
                         'why': 'Unchanged ADR0066 causal isogenic-model exception; reopened under '
                                'ADR0073 complete qualifying-source collection.'},
  'method': {'deposition': 'See pinned item summary and author_log for full input/test provenance',
             'model': 'Complete author DESeq2 Wald test, sign reversed to disease-minus-corrected',
             'design': 'TBX5G145R Brugada hiPSC-CM versus corrected isogenic control',
             'test': 'Author Wald',
             'padj': 'Benjamini-Hochberg'}},
 {'dataset_id': 'gse329393__rbm20',
  'category': 'genetic_cardiomyopathy',
  'source_item': 'disease_data/GSE329393',
  'local_deg': 'local/deg/GSE329393_deg.csv',
  'local_summary': 'local/deg/GSE329393_deg.json',
  'arm_or_class': 'RBM20 mutant hiPSC-CM versus family-matched corrected controls',
  'species': 'Homo sapiens',
  'model_system': 'RBM20 mutant hiPSC-CM versus family-matched corrected controls',
  'evidence_type': 'bulk',
  'tier': 2,
  'genotype_status': 'pathogenic',
  'time_window_basis': None,
  'exploratory': False,
  'shared_control_group': 'gse329393:isogenic',
  'n_case': 6,
  'n_control': 10,
  'contrasts': ['RBM20_LVNC-corrected', 'RBM20_DCM-corrected'],
  'isogenic_exception': {'applies': True,
                         'control_arm': 'family-matched corrected controls',
                         'why': 'Unchanged ADR0066 causal isogenic-model exception; reopened under '
                                'ADR0073 complete qualifying-source collection.'},
  'method': {'deposition': 'See pinned item summary and author_log for full input/test provenance',
             'model': 'limma-voom',
             'design': 'RBM20 mutant hiPSC-CM versus family-matched corrected controls',
             'test': 'moderated t',
             'padj': 'Benjamini-Hochberg'}}])
CATEGORIES['genetic_cardiomyopathy']['rejected'] = [{'candidate': 'GSE262876',
  'reason': 'Primary paper limitations identify one noncarrier and one D389V subject, with '
            'isogenic CRISPR comparisons still needed; donor background is confounded with '
            'genotype despite GEO abstract wording.',
  'source': 'https://doi.org/10.3390/cells13221913',
  'reentry': 'Verified genotype-matched isogenic comparison'},
 {'candidate': 'GSE274124',
  'reason': 'Patient iPSC31 versus unrelated WTC11 control; no isogenic correction, so '
            'donor-background confounding fails the causal isogenic exception.',
  'source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274124',
  'reentry': 'Isogenic corrected/introduced genotype or independently admissible patient-CM '
             'cohort'},
 {'candidate': 'GSE251993',
  'reason': 'Deposited raw_counts table has fractional values on an undeclared scale (column '
            'sums107.36–384.74, maximum18.52). It cannot be treated as raw counts or assumed log2 '
            'abundance without changing the effect threshold meaning.',
  'source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251993',
  'reentry': 'Author-defined expression units, valid raw counts, or full per-gene differential '
             'tests'},
 {'candidate': 'GSE140965',
  'reason': 'Deposited matrix has35 columns versus30 GSM, duplicated control names and extra '
            'Ex27Het groups; primary sample/control/batch mapping does not bind all columns. No '
            'positional matching or unsupported batch correction is invented.',
  'source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140965',
  'reentry': 'Author sample/control/batch mapping or complete differential test statistics'},
 {'candidate': 'GSE289628',
  'reason': 'AAV6-versus-siRNA vehicle is confounded with DSP intervention; unchanged '
            'causal-control restriction retained.',
  'source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289628'},
 {'candidate': 'Reichart-CM single NCCM donor',
  'reason': 'Single NCCM donor excluded by existing cohort admission; no new independent '
            'case-control evidence.'}]
CATEGORIES['genetic_cardiomyopathy']['source_review'] = {'reviewed_on': '2026-11-16',
 'policy': 'all_qualifying_declared_sources',
 'scope': 'Existing screened/rejected genetic candidate inventory; primary metadata and actual '
          'deposited tables checked, no sufficiency or fixed-roster exclusion',
 'records': [{'candidate': 'GSE248187',
              'disposition': 'admitted',
              'tier': 2,
              'provenance': 'disease_data/GSE248187/author_log.json',
              'source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248187',
              'reason': 'Human DSP knockout EHT versus isogenic control'},
             {'candidate': 'GSE221945',
              'disposition': 'admitted',
              'tier': 2,
              'provenance': 'disease_data/GSE221945/author_log.json',
              'source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221945',
              'reason': 'TBX5G145R Brugada hiPSC-CM versus corrected isogenic control'},
             {'candidate': 'GSE329393',
              'disposition': 'admitted',
              'tier': 2,
              'provenance': 'disease_data/GSE329393/author_log.json',
              'source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329393',
              'reason': 'RBM20 mutant hiPSC-CM versus family-matched corrected controls'},
             {'candidate': 'GSE262876',
              'reason': 'Primary paper limitations identify one noncarrier and one D389V subject, '
                        'with isogenic CRISPR comparisons still needed; donor background is '
                        'confounded with genotype despite GEO abstract wording.',
              'source': 'https://doi.org/10.3390/cells13221913',
              'reentry': 'Verified genotype-matched isogenic comparison',
              'disposition': 'not_admitted'},
             {'candidate': 'GSE274124',
              'reason': 'Patient iPSC31 versus unrelated WTC11 control; no isogenic correction, so '
                        'donor-background confounding fails the causal isogenic exception.',
              'source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274124',
              'reentry': 'Isogenic corrected/introduced genotype or independently admissible '
                         'patient-CM cohort',
              'disposition': 'not_admitted'},
             {'candidate': 'GSE251993',
              'reason': 'Deposited raw_counts table has fractional values on an undeclared scale '
                        '(column sums107.36–384.74, maximum18.52). It cannot be treated as raw '
                        'counts or assumed log2 abundance without changing the effect threshold '
                        'meaning.',
              'source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE251993',
              'reentry': 'Author-defined expression units, valid raw counts, or full per-gene '
                         'differential tests',
              'disposition': 'not_admitted'},
             {'candidate': 'GSE140965',
              'reason': 'Deposited matrix has35 columns versus30 GSM, duplicated control names and '
                        'extra Ex27Het groups; primary sample/control/batch mapping does not bind '
                        'all columns. No positional matching or unsupported batch correction is '
                        'invented.',
              'source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140965',
              'reentry': 'Author sample/control/batch mapping or complete differential test '
                         'statistics',
              'disposition': 'not_admitted'},
             {'candidate': 'GSE289628',
              'reason': 'AAV6-versus-siRNA vehicle is confounded with DSP intervention; unchanged '
                        'causal-control restriction retained.',
              'source': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289628',
              'disposition': 'not_admitted'},
             {'candidate': 'Reichart-CM single NCCM donor',
              'reason': 'Single NCCM donor excluded by existing cohort admission; no new '
                        'independent case-control evidence.',
              'disposition': 'not_admitted'}]}

# Every retained finite rejection also appears in the dated public review.
# Admission records above preserve richer method/provenance evidence.
for _category in CATEGORIES.values():
    _review = _category.setdefault('source_review', {
        'reviewed_on': '2026-11-16', 'policy': 'all_qualifying_declared_sources',
        'records': []})
    _review.setdefault('scope', 'Previously screened, rejected and owned source inventory; no claim of unbounded public discovery')
    _reviewed = {row['candidate'] for row in _review['records']}
    for _rejection in _category.get('rejected', []):
        if _rejection['candidate'] not in _reviewed:
            _review['records'].append(dict(_rejection, disposition='not_admitted'))
