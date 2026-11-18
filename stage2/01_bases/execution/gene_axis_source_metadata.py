"""Verified primary-source metadata for the gene-axis source table.

One record per published gene-axis source: the 41 ``dataset_id`` values the five
released ``stage2/01_bases/output/disease-deg/<category>/catalog.json`` files
declare, plus the separate ``tenv_hvg`` non-disease row. This module is pure data keyed by
``dataset_id`` and is the only place the rule document's accession and sequencing
technology columns read from (ticket #198 renders the table). It performs no I/O
and imports nothing outside the standard library.

Every value below was looked up at the primary source on ``VERIFIED_ON`` and is
never copied from a catalog's ``model_system`` prose. ``technology`` is the short
table phrase the primary record supports, ``platform`` carries the GPL id and
instrument behind that claim, and ``platform_note`` states the platform evidence
itself (the series record plus the GPL/instrument) or, where the primary record
publishes no GEO platform, why the absence is real rather than a harvesting gap.
``source_url`` is the record the row was checked against. ``accession_note``
records a series hierarchy where one exists. ``sample_size_note`` records how the
repo's ``n_case``/``n_control`` and the primary record's sample table reconcile
wherever they are not directly comparable (a different counting basis, a subset
of the series, or a split only derivable from the paper); it is ``None`` where
the primary count equals ``n_case + n_control``. The release's numbers are never
changed here.

``tenv_hvg`` is accession-less by design: the tenv-HVG rank is the published
descriptive per-well dispersion rank of the repo corpus built from ``GSE262419``
(ADR 0066), so ``GSE262419`` is the corpus source, not an axis accession. The
table cell carries only the short ``accession_note``; the long rationale lives
here and in ``platform_note``. ``TENV_TREATED_WELLS``/``TENV_VEHICLE_WELLS`` are
the published tenv corpus populations pinned in ``output/view-5k/manifest.json``
``sources.tenv`` (4,478 rows = 4,206 treated + 272 vehicle wells), so the one
non-catalog sample-size cell has a single source.

The historical DToxS class rows, the three GSE270788 rows and the three Reichart-CM rows
share byte-identical accession/technology/source/platform blocks through the
``_DTOXS_CLASS`` / ``_GSE270788`` / ``_REICHART`` base records: one accession,
several arms. Every other row is one accession and one arm.

Two platform labels in the items' own tracked DEG summaries disagreed with the
GPL record they cite, and ``platform_note`` records the history rather than
silently repeating either reading: ``GSE165303``'s DEG summary labels
``GPL20795`` "Illumina NovaSeq 6000" while the ``GPL20795`` record is *HiSeq X
Ten*, and ``GSE116250``'s summary labels ``GPL16791`` "Illumina HiSeq 2000"
while the ``GPL16791`` record is *Illumina HiSeq 2500*. The GPL record is the
platform claim, so ``platform`` carries it; #211 corrected the two items'
``author_log.json`` and their DEG-build scripts, and the published
``local/deg/*_deg.json`` bytes keep the pre-correction string because they are
digest-pinned by the accepted disease-DEG catalogs (see the notes below).

Harvest method (recorded so a reviewer can re-run it):

* NCBI E-utilities ``esummary`` (``db=gds``) for every GSE and GPL uid.
* GEO FTP ``matrix/<GSE>_series_matrix.txt.gz``, read only up to the
  ``!series_matrix_table_begin`` line, for the sample count and overall design.
* GEO FTP ``suppl/`` listings and ``filelist.txt``, plus the deposits' own
  published supplementary metadata (``GSE270788_metadata.csv.gz`` for the
  ``HF.etiology`` donor arms; the DToxS ``Conv-DEG-Calc`` Calc-table headers for
  the post-QC well columns).
* The cellxgene Curation API v1 collection endpoint for the Reichart-CM object
  and for the Kuppe-MI object.
* NCBI PMC full-text XML (PMC9528698, Reichart et al., Science 2022) for the
  donor-cohort split.
* For the nine rows #209 added, the platform half is re-read from the same
  E-utilities ``esummary`` (``db=gds``) call on the dataset: the GSE record's
  ``gpl`` and ``n_samples`` fields plus the ``GPL[ACCN] AND gpl[ETYP]`` record
  title, and the sample-size half from the item's tracked
  ``local/deg/*_deg.json`` / ``author_log.json`` (the 2020-batch slice of
  GSE138252, the snRNA arm of GSE238242, the 500-CM-nuclei profile floor of
  Kuppe-MI). No value below is copied from a ticket or a catalog's prose.

Requests were limited to at most 3 per second with exponential backoff, and every
response was cached under ``.scratch/source-verify/raw/``. No expression payload
was read — metadata only. That raw harvest is untracked and will not exist in
another worktree; the durable human-readable record of what was verified is
issue #197's posted comment.
"""
from __future__ import annotations

from collections.abc import Iterable

SCHEMA = 'stage2-gene-axis-source-metadata-v1'
VERIFIED_ON = '2026-11-15'
TENV_HVG = 'tenv_hvg'
#: Published tenv corpus populations (``output/view-5k/manifest.json``
#: ``sources.tenv``): 4,478 rows = ``TENV_TREATED_WELLS + TENV_VEHICLE_WELLS``.
TENV_TREATED_WELLS = 4206
TENV_VEHICLE_WELLS = 272

_REICHART_COLLECTION = 'cellxgene collection e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb'
_REICHART_URL = 'https://cellxgene.cziscience.com/collections/e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb'

_ENTRY_FIELDS = (
    'accession', 'accession_note', 'technology', 'source_url',
    'verified_on', 'platform', 'platform_note', 'sample_size_note',
)
_PLATFORM_FIELDS = ('gpl', 'instrument')

# --- shared base records: one accession, several released arms ---------------
_GSE270788 = {
    'accession': 'GSE270788',
    'accession_note': 'SubSeries of SuperSeries GSE218392',
    'technology': '10x Multiome snRNA-seq (GEX arm, CellRanger ARC)',
    'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270788',
    'verified_on': VERIFIED_ON,
    'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
    'platform_note': ('GEO series record and GPL24676 platform record '
                      '(Illumina NovaSeq 6000)'),
}
_REICHART = {
    'accession': _REICHART_COLLECTION,
    'technology': "10x 3' snRNA-seq (CM-only cellxgene object)",
    'source_url': _REICHART_URL,
    'verified_on': VERIFIED_ON,
    'platform': {'gpl': None, 'instrument': None},
    'platform_note': ('no GEO platform: cellxgene h5ad object; the cellxgene record '
                      'and the paper do not name the instrument'),
}
_DTOXS_CLASS = {
    'accession': 'GSE217421',
    'accession_note': 'Conv arm of SuperSeries GSE217424',
    'technology': 'bulk RNA-seq (hiPSC-CM plate wells, 48 h)',
    'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217421',
    'verified_on': VERIFIED_ON,
    'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
    'platform_note': ('GEO series record and GPL24676 platform record '
                      '(Illumina NovaSeq 6000)'),
}

SOURCES = {
    # --- nHF: GSE270788, a 10x Multiome series whose GEX arm the repo uses -----
    'gse270788__nicm': {
        **_GSE270788,
        'sample_size_note': (
            'The series reports 46 GSMs = 23 GEX (snRNA) + 23 ATAC libraries; the repo arms '
            'count the 23 donors instead (deposit metadata HF.etiology: Donor 7, NICM 6, '
            'ICM 5, AMI 5), which is the 6/7 donor basis of this row.'
        ),
    },
    'reichart-cm__dcm_genotype_negative': {
        **_REICHART,
        'accession_note': (
            'cellxgene dataset f7995301-7551-4e1d-8396-ffe3c9497ace, '
            'version 27fdbdf9-deae-4aed-bf7f-21964ee62137'
        ),
        'sample_size_note': (
            'The collection carries 79 donors (61 failing + 18 controls); the repo arms are '
            'the 8 PV-negative DCM donors and the 18 controls stated in the paper, i.e. '
            '79 = 78 exported (8 + 44 + 8 + 18) + 1 NCCM.'
        ),
    },
    'gse35108__dcm_ipsc_cm': {
        'accession': 'GSE35108',
        'accession_note': None,
        'technology': 'one-colour microarray (RMA log2; human patient-derived iPSC-CM)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35108',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL6244',
                     'instrument': '[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array'},
        'platform_note': (
            'GEO series record and GPL6244 platform record ([HuGene-1_0-st] Affymetrix '
            'Human Gene 1.0 ST Array [transcript (gene) version]); this is an array, so '
            'the platform record names the array rather than a sequencer'
        ),
        'sample_size_note': (
            'The series carries 17 GSMs, of which 8 are cardiomyocyte libraries in a 2x2 '
            'genotype x Serca2a design (2 per cell); this repo arm is the 4 DCM vs 4 normal '
            'CM libraries across both treatment strata, not all 17.'
        ),
    },
    'gse165303__dcm_lv_bulk': {
        'accession': 'GSE165303',
        'accession_note': ('RNA-seq arm of a multi-assay series '
                           '(HiChIP / Hi-C / ChIP-seq / ATAC-seq / RNA-seq)'),
        'technology': 'bulk RNA-seq (human left-ventricle whole tissue, log2 normalised counts)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165303',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL20795', 'instrument': 'HiSeq X Ten'},
        'platform_note': (
            'GEO series record and GPL20795 platform record (HiSeq X Ten); the item\'s '
            'published DEG summary still labels GPL20795 "Illumina NovaSeq 6000", which '
            'the GPL record does not support, so the GPL record is the platform claim '
            'here (#211 corrected the item\'s author_log.json and script/build_deg.py)'
        ),
        'sample_size_note': (
            'The series carries 177 samples across the HiChIP / Hi-C / ChIP-seq / ATAC-seq / '
            'RNA-seq arms; this repo arm is the 101-library bulk RNA-seq arm (50 DCM + 51 '
            'non-failing), not the whole series.'
        ),
    },
    'gse116250__hf_lv_bulk': {
        'accession': 'GSE116250',
        'accession_note': None,
        'technology': 'bulk RNA-seq (human left-ventricular tissue, RPKM as deposited)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116250',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL16791', 'instrument': 'Illumina HiSeq 2500'},
        'platform_note': (
            'GEO series record and GPL16791 platform record (Illumina HiSeq 2500); the '
            'item\'s published DEG summary still calls the same GPL16791 "Illumina HiSeq '
            '2000", so the GPL record is the platform claim here (#211 corrected the '
            'item\'s author_log.json and script/build_deg.R)'
        ),
        'sample_size_note': None,
    },
    'gse141910__hf_lv_bulk': {
        'accession': 'GSE141910',
        'accession_note': None,
        'technology': ('bulk RNA-seq (human left-ventricular free-wall tissue, '
                       'submitted voom log2)'),
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141910',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL16791', 'instrument': 'Illumina HiSeq 2500'},
        'platform_note': (
            'GEO series record and GPL16791 platform record (Illumina HiSeq 2500)'
        ),
        'sample_size_note': None,
    },
    # --- ICM: the same GSE270788 GEX arm, split by aetiology -------------------
    'gse270788__icm': {
        **_GSE270788,
        'sample_size_note': (
            'The series reports 46 GSMs = 23 GEX + 23 ATAC; this repo arm is the 5 ICM '
            'donors against the same 7 donor controls (HF.etiology), not the 46 libraries.'
        ),
    },
    'gse270788__ami': {
        **_GSE270788,
        'sample_size_note': (
            'The series reports 46 GSMs = 23 GEX + 23 ATAC; this repo arm is the 5 AMI '
            'donors against the same 7 donor controls (HF.etiology), not the 46 libraries.'
        ),
    },
    'kuppe-mi__iz': {
        'accession': 'cellxgene collection 8191c283-0816-424b-9b61-c3e1d6258a77',
        'accession_note': (
            'cellxgene dataset 1c739a3e-c3f5-49d5-98e0-73975e751201, '
            'version 4780eb7d-4c7c-4630-8705-471f0de8e27a'
        ),
        'technology': "10x 3' v3 snRNA-seq (CM pseudobulk, per donor-zone profile)",
        'source_url': ('https://cellxgene.cziscience.com/collections/'
                       '8191c283-0816-424b-9b61-c3e1d6258a77'),
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': None, 'instrument': None},
        'platform_note': (
            'no GEO platform: cellxgene h5ad object; the cellxgene record and the paper '
            'name the 10x 3\' v3 chemistry but not the instrument'
        ),
        'sample_size_note': (
            'The object carries 191,795 nuclei from 20 donors; this repo arm is the 4 '
            'ischaemic-zone profiles (donors P2, P9, P10, P15) vs the 4 normal-heart profiles '
            '(P1, P7, P8, P17) left by the 500-CM-nuclei profile floor, not the 20 donors.'
        ),
    },
    # --- genetic_cardiomyopathy ------------------------------------------------
    'reichart-cm__dcm_genotype_positive': {
        **_REICHART,
        'accession_note': 'cellxgene dataset f7995301-7551-4e1d-8396-ffe3c9497ace',
        'sample_size_note': (
            'Partial match: the repo n_case=44 is derivable (61 failing - 8 PV-negative '
            '- 8 ARVC - 1 NCCM) but no fetched primary record states 44 directly; the '
            'per-patient DCM-vs-ACM labels live only in the paper\'s Tables S1/S2.'
        ),
    },
    'reichart-cm__arvc': {
        **_REICHART,
        'accession_note': 'cellxgene dataset f7995301-7551-4e1d-8396-ffe3c9497ace',
        'sample_size_note': (
            'Partial match: the repo n_case=8 is the object\'s implied ARVC donor count '
            '(61 - 8 PV-negative - 44 DCM-PV+ - 1 NCCM); the paper states PKP2 n=6 '
            'explicitly but never states "8 ARVC".'
        ),
    },
    'gse224129__mybpc3_hcm': {
        'accession': 'GSE224129',
        'accession_note': None,
        'technology': 'bulk RNA-seq (hiPSC-CM engineered cardiac tissue)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224129',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'platform_note': ('GEO series record and GPL24676 platform record '
                          '(Illumina NovaSeq 6000)'),
        'sample_size_note': None,
    },
    'gse234141__ttntv_dcm': {
        'accession': 'GSE234141',
        'accession_note': None,
        'technology': 'bulk RNA-seq (hiPSC-CM 2D and engineered heart tissue)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234141',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL28337', 'instrument': 'DNBSEQ-G50'},
        'platform_note': 'GEO series record and GPL28337 platform record (DNBSEQ-G50)',
        'sample_size_note': None,
    },
    'gse304575__lmna_h222p_dcm': {
        'accession': 'GSE304575',
        'accession_note': None,
        'technology': 'bulk RNA-seq (hiPSC-CM, isogenic)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304575',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'platform_note': ('GEO series record and GPL24676 platform record '
                          '(Illumina NovaSeq 6000)'),
        'sample_size_note': None,
    },
    'gse305648__rbm20_dcm': {
        'accession': 'GSE305648',
        'accession_note': None,
        'technology': 'bulk RNA-seq (total RNA, rRNA-depleted; hiPSC cardioids)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305648',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL30173', 'instrument': 'NextSeq 2000'},
        'platform_note': 'GEO series record and GPL30173 platform record (NextSeq 2000)',
        'sample_size_note': None,
    },
    'gse253223__pkp2_knockdown_arvc': {
        'accession': 'GSE253223',
        'accession_note': 'SubSeries of GSE253226',
        'technology': 'bulk RNA-seq (hiPSC-CM, siRNA knockdown)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253223',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'platform_note': ('GEO series record and GPL24676 platform record '
                          '(Illumina NovaSeq 6000)'),
        'sample_size_note': None,
    },
    'gse194297__dmd_cardiomyopathy': {
        'accession': 'GSE194297',
        'accession_note': None,
        'technology': 'bulk RNA-seq (hiPSC cardiac organoids)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194297',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL20301', 'instrument': 'Illumina HiSeq 4000'},
        'platform_note': 'GEO series record and GPL20301 platform record (Illumina HiSeq 4000)',
        'sample_size_note': None,
    },
    # --- AF: three tier-1 CM-specific sources, one tier-3 bulk atrial source ----
    'gse255612__af': {
        'accession': 'GSE255612',
        'accession_note': None,
        'technology': "10x 3' v3 snRNA-seq (left-atrium CM pseudobulk)",
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255612',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'platform_note': (
            'GEO series record and GPL24676 platform record (Illumina NovaSeq 6000)'
        ),
        'sample_size_note': None,
    },
    'gse138252__af': {
        'accession': 'GSE138252',
        'accession_note': None,
        'technology': 'bulk RNA-seq (PCM1+ cardiomyocyte nuclei, 2020 batch)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138252',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL16791', 'instrument': 'Illumina HiSeq 2500'},
        'platform_note': (
            'GEO series record and GPL16791 platform record (Illumina HiSeq 2500)'
        ),
        'sample_size_note': (
            'The series carries 48 GSMs / 38 distinct library columns across the 2020 matrix '
            'and the 2021 addendum; this repo arm is the batch-matched 2020 left-atrial slice '
            '(5 AF+HF vs 5 non-failing) and does not use the addendum\'s 6 HF-only + 14 '
            'control libraries.'
        ),
    },
    'gse238242__af': {
        'accession': 'GSE238242',
        'accession_note': None,
        'technology': '10x multiome snRNA-seq (snRNA arm, author CM label)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE238242',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'platform_note': (
            'GEO series record and GPL24676 platform record (Illumina NovaSeq 6000)'
        ),
        'sample_size_note': (
            'The series carries 14 GSMs = 7 snRNA + 7 snATAC libraries; this repo arm is the '
            '3 persistent-AF vs 4 sinus-rhythm donors of the snRNA arm, not the 14 libraries.'
        ),
    },
    'gse269882__af': {
        'accession': 'GSE269882',
        'accession_note': None,
        'technology': 'bulk RNA-seq (human atrial appendage whole tissue)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269882',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL16791', 'instrument': 'Illumina HiSeq 2500'},
        'platform_note': (
            'GEO series record and GPL16791 platform record (Illumina HiSeq 2500)'
        ),
        'sample_size_note': None,
    },
    # --- unified TOX --------------------------------------------------------
    'gse243674': {
        'accession': 'GSE243674',
        'accession_note': None,
        'technology': 'bulk RNA-seq (hiPSC-CM monolayer, 3 h / 24 h)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243674',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL21697', 'instrument': 'NextSeq 550'},
        'platform_note': 'GEO series record and GPL21697 platform record (NextSeq 550)',
        'sample_size_note': (
            'The series reports 72 GSMs = 6 donors x 6 conditions x 2 time points; the repo '
            'contrast is the 3 h slice (4 anthracyclines x 6 donors = 24 case vs 6 vehicle), '
            'not the whole series.'
        ),
    },
    'dtoxs_ipsc': {
        **_DTOXS_CLASS,
        'accession_note': 'Conv arm of SuperSeries GSE217424; sibling DGE arm GSE217423',
        'sample_size_note': 'Pooled admitted drugs within matched cell-line and plate strata; n_case_basis and n_control_basis are per-contrast well counts, not GSM totals.',
    },
    'gse244740': {
        'accession': 'GSE244740',
        'accession_note': None,
        'technology': 'bulk RNA-seq (384-well ScreenSeq hiPSC-CM plate, 24 h)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244740',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'platform_note': ('GEO series record and GPL24676 platform record '
                          '(Illumina NovaSeq 6000)'),
        'sample_size_note': (
            'The series reports 3,072 GSMs = 8 x 384-well plates (all 24 h). '
            'The unified TOX source pools 32 admitted drugs within each of eight disjoint '
            'plate strata versus same-plate DMSO: 1,920 treated wells and 536 control wells '
            'in total. Per-contrast counts and the summed-well basis are retained in the '
            'catalog; these are wells, not independent donors.'
        ),
    },
    # --- tenv-HVG fusion row (not a DEG source) -------------------------------
    TENV_HVG: {
        'accession': None,
        'accession_note': '—（无 accession；语料源 GSE262419）',
        'technology': 'TempO-Seq targeted transcriptomic panel (96-well, HiSeq 2500)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262419',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL16791', 'instrument': 'Illumina HiSeq 2500'},
        'platform_note': (
            'GEO series record and GPL16791 platform record (Illumina HiSeq 2500); the row '
            'has no axis accession because GSE262419 is the corpus source behind the '
            'published tenv-HVG rank (ADR 0066), not an axis dataset'
        ),
        'sample_size_note': (
            f'Not a contrast, so there is no n_case/n_control: the row ranks the '
            f'{TENV_TREATED_WELLS:,} treated wells of the GSE262419 corpus (4,558 GSMs over '
            f'16 plates; the corpus keeps 4,510 native wells x 19,613 genes = '
            f'{TENV_TREATED_WELLS:,} treated + {TENV_VEHICLE_WELLS:,} vehicle + 32 blanks).'
        ),
    },
}


# ADR 0072 dataset-level TOX metadata, checked against the pinned GEO series headers.
SOURCES.update({'gse106297': {'accession': 'GSE106297',
               'accession_note': None,
               'technology': 'bulk RNA-seq (Illumina HiSeq 2000)',
               'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106297',
               'verified_on': '2026-11-16',
               'platform': {'gpl': 'GPL11154', 'instrument': 'Illumina HiSeq 2000'},
               'platform_note': 'Pinned deposited series matrix Sample_platform_id and '
                                'Sample_instrument_model headers.',
               'sample_size_note': 'Pooled admitted drug versus matched vehicle contrasts; '
                                   'published n_case_basis and n_control_basis retain the '
                                   'per-contrast counting basis.'},
 'gse157282': {'accession': 'GSE157282',
               'accession_note': None,
               'technology': 'bulk RNA-seq (Illumina NovaSeq 6000)',
               'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157282',
               'verified_on': '2026-11-16',
               'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
               'platform_note': 'Pinned deposited series matrix Sample_platform_id and '
                                'Sample_instrument_model headers.',
               'sample_size_note': 'Pooled admitted drug versus matched vehicle contrasts; '
                                   'published n_case_basis and n_control_basis retain the '
                                   'per-contrast counting basis.'},
 'gse232331': {'accession': 'GSE232331',
               'accession_note': None,
               'technology': 'bulk RNA-seq (NextSeq 500)',
               'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232331',
               'verified_on': '2026-11-16',
               'platform': {'gpl': 'GPL18573', 'instrument': 'NextSeq 500'},
               'platform_note': 'Pinned deposited series matrix Sample_platform_id and '
                                'Sample_instrument_model headers.',
               'sample_size_note': 'Pooled admitted drug versus matched vehicle contrasts; '
                                   'published n_case_basis and n_control_basis retain the '
                                   'per-contrast counting basis.'},
 'gse235470': {'accession': 'GSE235470',
               'accession_note': None,
               'technology': 'bulk RNA-seq (NextSeq 500)',
               'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235470',
               'verified_on': '2026-11-16',
               'platform': {'gpl': 'GPL18573', 'instrument': 'NextSeq 500'},
               'platform_note': 'Pinned deposited series matrix Sample_platform_id and '
                                'Sample_instrument_model headers.',
               'sample_size_note': 'Pooled admitted drug versus matched vehicle contrasts; '
                                   'published n_case_basis and n_control_basis retain the '
                                   'per-contrast counting basis.'},
 'gse242692': {'accession': 'GSE242692',
               'accession_note': None,
               'technology': 'bulk RNA-seq (Illumina NovaSeq 6000)',
               'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242692',
               'verified_on': '2026-11-16',
               'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
               'platform_note': 'Pinned deposited series matrix Sample_platform_id and '
                                'Sample_instrument_model headers.',
               'sample_size_note': 'Pooled admitted drug versus matched vehicle contrasts; '
                                   'published n_case_basis and n_control_basis retain the '
                                   'per-contrast counting basis.'},
 'gse247541': {'accession': 'GSE247541',
               'accession_note': None,
               'technology': 'bulk RNA-seq (Illumina NovaSeq 6000)',
               'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247541',
               'verified_on': '2026-11-16',
               'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
               'platform_note': 'Pinned deposited series matrix Sample_platform_id and '
                                'Sample_instrument_model headers.',
               'sample_size_note': 'Pooled admitted drug versus matched vehicle contrasts; '
                                   'published n_case_basis and n_control_basis retain the '
                                   'per-contrast counting basis.'},
 'gse255505': {'accession': 'GSE255505',
               'accession_note': None,
               'technology': 'bulk RNA-seq (Illumina NovaSeq 6000)',
               'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255505',
               'verified_on': '2026-11-16',
               'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
               'platform_note': 'Pinned deposited series matrix Sample_platform_id and '
                                'Sample_instrument_model headers.',
               'sample_size_note': 'Pooled admitted drug versus matched vehicle contrasts; '
                                   'published n_case_basis and n_control_basis retain the '
                                   'per-contrast counting basis.'},
 'gse262419': {'accession': 'GSE262419',
               'accession_note': None,
               'technology': 'TempO-Seq targeted transcriptomics (Illumina HiSeq 2500)',
               'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262419',
               'verified_on': '2026-11-16',
               'platform': {'gpl': 'GPL16791', 'instrument': 'Illumina HiSeq 2500'},
               'platform_note': 'Pinned deposited series matrix Sample_platform_id and '
                                'Sample_instrument_model headers.',
               'sample_size_note': 'Pooled admitted drug versus matched vehicle contrasts; '
                                   'published n_case_basis and n_control_basis retain the '
                                   'per-contrast counting basis.'},
 'gse313746': {'accession': 'GSE313746',
               'accession_note': None,
               'technology': 'bulk RNA-seq (Element AVITI)',
               'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313746',
               'verified_on': '2026-11-16',
               'platform': {'gpl': 'GPL34295', 'instrument': 'Element AVITI'},
               'platform_note': 'Pinned deposited series matrix Sample_platform_id and '
                                'Sample_instrument_model headers.',
               'sample_size_note': 'Pooled admitted drug versus matched vehicle contrasts; '
                                   'published n_case_basis and n_control_basis retain the '
                                   'per-contrast counting basis.'}})

for _entry in SOURCES.values():
    if _entry['verified_on'] == '2026-11-16':
        _entry['platform_note'] += ' ' + _entry['platform']['gpl'] + ' ' + _entry['platform']['instrument']


# ADR0073 reopened sources. Primary series/GPL records and author tables reviewed
# 2026-11-16; counts below describe the declared analysis unit, not whole series.
SOURCES['gse116250__icm_lv_bulk'] = dict(SOURCES['gse116250__hf_lv_bulk'],
    verified_on='2026-11-16', sample_size_note='ICM arm:13 patients versus14 non-failing; shared NF controls with the independent DCM arm are declared, not counted as additional donors.')


def _successor_source(accession, technology, gpl, instrument, sample_note):
    return {'accession': accession, 'accession_note': None, 'technology': technology,
            'source_url': f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}',
            'verified_on': '2026-11-16', 'platform': {'gpl': gpl, 'instrument': instrument},
            'platform_note': f'Official GEO series platform and instrument records: {gpl} {instrument}; matched human RNA analysis arm only.',
            'sample_size_note': sample_note}


SOURCES.update({
    'gse128188__af': _successor_source('GSE128188', 'bulk RNA-seq (human atrial appendage)',
        'GPL18573', 'Illumina NextSeq 500',
        '5 AF versus5 sinus-rhythm donors; paired LA/RA libraries pooled per donor,20 libraries total.'),
    'gse183852__dcm_cm': _successor_source('GSE183852', 'single-nucleus RNA-seq (author cardiomyocyte full statistics)',
        'GPL24676', 'Illumina NovaSeq 6000',
        'Official snRNA arm13 DCM/25 donor libraries; full author Cardiomyocytes table retained. Library counts are not asserted to be unique donor counts; no local refit claimed.'),
    'gse298023__dcm_cm': _successor_source('GSE298023', '10x single-nucleus RNA-seq (author CM donor pseudobulk)',
        'GPL24676', 'Illumina NovaSeq 6000',
        '6 DCM hearts versus4 NF; author Cardiomyocytes labels, donor-pooled counts and run-pattern-adjusted joint model.'),
    'gse298023__icm_cm': _successor_source('GSE298023', '10x single-nucleus RNA-seq (author CM donor pseudobulk)',
        'GPL24676', 'Illumina NovaSeq 6000',
        '7 unique ICM hearts versus4 NF; repeated scar/non-scar regions pooled once per heart, shared NF controls declared.'),
    'gse329393__rbm20': _successor_source('GSE329393', 'bulk RNA-seq (RBM20 hiPSC-CM)',
        'GPL24676', 'Illumina NovaSeq 6000',
        '6 mutant versus10 corrected libraries across2 matched families; joint family-specific contrasts, not16 independent donors.'),
    'gse248187__dsp_eht': _successor_source('GSE248187', 'bulk RNA-seq (human DSP knockout EHT author edgeR statistics)',
        'GPL24676', 'Illumina NovaSeq 6000',
        '3 DSP knockout versus3 isogenic-control human EHT libraries; other rat arms excluded; complete author test retained without count refit.'),
    'gse221945__tbx5_brs': _successor_source('GSE221945', 'bulk RNA-seq (TBX5 hiPSC-CM author DESeq2 statistics)',
        'GPL21290', 'Illumina HiSeq 3000',
        '4 Brugada versus4 corrected isogenic libraries explicitly named in the complete author test; remaining population controls and knockout libraries excluded.'),
})


def validate(dataset_ids: Iterable[str]) -> None:
    """Raise ``ValueError`` unless ``dataset_ids`` exactly matches the published coverage.

    ``dataset_ids`` is the caller's list of published dataset ids plus ``tenv_hvg``
    (the five catalogs' concatenation, in any order). It must contain no duplicate;
    its set must equal ``SOURCES``' key set — nothing missing, nothing extra; and
    every entry must be complete: an accession or a stated accession-less reason,
    a non-empty technology, source URL and verification date, a non-empty
    ``platform_note`` (which, when both ``gpl`` and ``instrument`` are null, must
    explain the absence), and a platform dict with both ``gpl`` and ``instrument``.
    """
    identifiers = list(dataset_ids)
    duplicates = sorted({dataset_id for dataset_id in identifiers
                         if identifiers.count(dataset_id) > 1})
    if duplicates:
        raise ValueError(f'duplicate dataset_id(s): {duplicates}')
    missing = sorted(set(SOURCES) - set(identifiers))
    extra = sorted(set(identifiers) - set(SOURCES))
    if missing or extra:
        raise ValueError(
            f'{SCHEMA} coverage mismatch: {len(missing)} missing {missing}, '
            f'{len(extra)} extra {extra}'
        )
    for dataset_id, entry in SOURCES.items():
        if set(entry) != set(_ENTRY_FIELDS):
            raise ValueError(
                f'{dataset_id}: entry fields {sorted(entry)} != {sorted(_ENTRY_FIELDS)}'
            )
        accession = (entry['accession'] or '').strip()
        accession_note = (entry['accession_note'] or '').strip()
        if not (accession or accession_note):
            raise ValueError(f'{dataset_id}: needs a non-empty accession or accession_note')
        for field in ('technology', 'source_url'):
            value = entry[field]
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f'{dataset_id}: {field} must be a non-empty string')
        if entry['verified_on'] not in (VERIFIED_ON, '2026-11-16'):
            raise ValueError(
                f'{dataset_id}: verified_on {entry["verified_on"]!r} != {VERIFIED_ON!r}'
            )
        platform = entry['platform']
        if set(platform) != set(_PLATFORM_FIELDS):
            raise ValueError(
                f'{dataset_id}: platform keys {sorted(platform)} != {sorted(_PLATFORM_FIELDS)}'
            )
        for key, value in platform.items():
            if value is not None and not isinstance(value, str):
                raise ValueError(f'{dataset_id}: platform[{key!r}] must be str or None')
        # Every row needs the platform evidence, or an explicit statement that the
        # primary record publishes none (the three Reichart-CM rows).
        platform_note = entry['platform_note']
        if not isinstance(platform_note, str) or not platform_note.strip():
            raise ValueError(f'{dataset_id}: platform_note must be a non-empty string')
        sample_size_note = entry['sample_size_note']
        if sample_size_note is not None and not isinstance(sample_size_note, str):
            raise ValueError(f'{dataset_id}: sample_size_note must be str or None')
