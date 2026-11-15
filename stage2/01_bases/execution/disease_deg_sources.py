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

The ticket-187 categories are ``nHF`` and ``ICM`` only (the ``AF`` category was
removed by ADR 0069 because its only dataset carried no ordering evidence). The
six isogenic models (#188), the drug-class sources (#189) and the Reichart
genotype-positive / ARVC DEG (also #188) are other tickets and are not declared
here.
"""
from __future__ import annotations

CATEGORIES = {
    'nHF': {
        'directory': 'nhf',
        'tier': 1,
        'evidence_type': 'cm_pseudobulk',
        'shared_control_group': ['gse270788:donor', 'reichart-cm:normal'],
        # NOTE (docs sweep #211, 2026-11-15). The strings below are DECLARATION
        # DATA: `build_disease_deg.py` copies them verbatim into the published
        # nhf/ and icm/ `catalog.json`, whose digests the accepted disease-fusion,
        # gene-axis and view manifests pin. They are therefore NOT edited here
        # (that would desynchronise the frozen v2 release from its catalogs and
        # fail `build_disease_deg.py verify`); this comment is the correction of
        # record until the next disease-DEG republish rewrites them.
        #
        # GSE183852's real obstacle is measured, not the label-extraction claim
        # below: its authors' 6.7 MB `MOESM4` `Cardiomyocytes` sheet carries the
        # CM labels, and the v2 definition yields 70 DEGs (padj < 0.05 and
        # |log2FC| >= 2), 50 of them on the frozen axis (26,946 genes tested; the
        # padj-only anchors 6,280 at |log2FC| >= 0.5 and 3,597 at >= 1 reproduce).
        # Reading an already-annotated `meta.data` would not be re-clustering, so
        # the "labels live only inside Robj objects / ad-hoc re-clustering is
        # forbidden" reason misdescribes the trade-off: the 8.1-12 GB Robj label
        # extraction is simply not worth 2.5% of a 2,000-gene target. Its
        # padj-only yield (4,927 on-axis) would only matter under a softer fold
        # tier, which the v2 rule does not admit -- so nHF closes at 1,661
        # axis-resolved DEGs by decision, not by missing capability.
        #
        # The third entry's "not cardiomyocyte-derived" ground is superseded by
        # ADR 0071 tier 3: GSE141910, GSE116250 and GSE128188 are admitted (or
        # available) as whole-tissue tier-3 sources, and GSE165303 joined them;
        # only GSE282618 remains un-admitted, for its 3-vs-2 design.
        'rejected': [
            {'candidate': 'GSE183852',
             'reason': 'eligible human LV CM source, but CM labels live only inside '
                       '2.1-12 GB Seurat Robj objects; ad-hoc re-clustering is forbidden'},
            {'candidate': 'cellxgene myocarditis CM object',
             'reason': 'shares its normal donors with the Reichart object, so the two '
                       'cannot both count as independent sources'},
            {'candidate': 'GSE141910 / GSE116250 / GSE128188 / GSE282618',
             'reason': 'whole-heart/whole-ventricle bulk tissue, not cardiomyocyte-derived'},
            {'candidate': 'GSE185100',
             'accession': 'GSE185100',
             'status': 'unverified_not_admitted',
             'size': ('757 MB RAW.tar (722 MB of per-sample processed tars; the two '
                      'recordings disagree)'),
             'donor_counts': ('4 early-onset DCM vs 3 non-diseased LV-apex donors '
                              '(pediatric, 4-14 y)'),
             'reason': ('recorded unverified and not admitted (#203): the NCBI FTP '
                        'throttled the 757 MB RAW.tar to roughly 50 KB/s and the '
                        'download was abandoned, so the deposit was never opened. The '
                        'series is also 4 vs 3 donors of a pediatric early-onset '
                        'cohort, which could not meet the evidence bar even if the '
                        'archive were readable. Nothing was downloaded and no gene-axis '
                        'evidence is derived from it.')},
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
            'tier-3 source only and enters nHF only after tiers 1-2 fall short of 2,000 '
            'DEGs (ADR 0071)'),
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
# pseudobulk evidence at tier 1 — ADR 0071's "ICM gains the Kuppe MI IZ vs CTRL
# contrast ... and stops at tier 1". The object is gated on its own CM label and
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
# sources (GSE255612, GSE138252, GSE238242) are tier 1 and their measured DEG union
# is far short of 2,000, so the category descends to tier-3 whole-atrial bulk
# (GSE269882; GSE128188 is already owned). Every path below is item-relative; the
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
    'tier3_owned': [
        {'candidate': 'GSE128188',
         'item': 'disease_data/GSE128188',
         'reason': ('already owned human atrial bulk item (AF 5 vs 5 SR, paired LA+RA). '
                    'ADR 0071 admits whole-tissue bulk only at tier 3, so this is a '
                    'candidate for the #206 roster; it is not a rejection and it is not '
                    'declared here because #204 builds only the three new items.')},
    ],
    'tiering': {
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
        'stop_decision': ('tier 1 union 1 < 2,000, so the collection descends to tier 3; '
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
                        'so it is a tier-3 source only and enters AF only after tiers 1-2 '
                        'fall short of 2,000 DEGs (ADR 0071), which supersedes the ADR '
                        '0070 whole-tissue exclusion'),
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

# --- ticket #206: nHF tier-3 bulk roster + the six-category tier ladders ------
# ADR 0071 decision 2 collects each category's sources in ascending tier order.
# The ladder is derived from each record's own `tier` (no second declaration):
#   nHF   1 = gse270788__nicm + reichart genotype-negative
#         2 = gse35108__dcm_ipsc_cm
#         3 = gse165303__dcm_lv_bulk + GSE116250 + GSE141910
#   ICM   1 = gse270788__icm + gse270788__ami + kuppe-mi__iz
#   AF    1 = gse255612__af + gse138252__af + gse238242__af   (no tier-2 source)
#         3 = gse269882__af
#   genetic_cardiomyopathy  1 = the two Reichart exports
#                           2 = the six isogenic hiPSC-CM models
#   cardiotox_TOP2 / cardiotox_other  2 = their pooled class sources (the
#         category begins at tier 2; `bases.cardiotox` declares those sources)
# The stop rule and its axis-resolved criterion live once, in
# `bases.disease_deg.TIER_STOP_RULE` / `tier_ledger`; `build_disease_deg.py ledger`
# computes a category's cumulative axis-resolved union and stop point from these
# declarations, and the publish path writes the resulting `tiering` record.
# nHF's ceiling is accepted (ADR 0071 decision-2 note, 2026-11-15): GSE183852 is
# closed on both routes (measured 70 v2 DEGs / 50 on the axis; an 8.1-12 GB Seurat
# Robj label extraction is not warranted), so nHF tops out at ~1,661 axis DEGs and
# `DEG2000_nHF` is the whole union rather than a 2,000-gene cut.
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
                        'so it is a tier-3 source only and enters nHF only after tiers '
                        '1-2 fall short of 2,000 axis-resolved DEGs (ADR 0071)'),
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
                    'end-stage remodelling proxy, recorded here and never re-used as the '
                    'ICM category source, whose evidence is tier-1 CM pseudobulk'),
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
                        'so it is a tier-3 source only and enters nHF only after tiers 1-2 '
                        'fall short of 2,000 axis-resolved DEGs (ADR 0071)'),
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
