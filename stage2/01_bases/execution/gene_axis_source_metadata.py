"""Verified primary-source metadata for the gene-axis source table.

One record per published gene-axis source: the 18 ``dataset_id`` values the five
released ``stage2/01_bases/output/disease-deg/<category>/catalog.json`` files
declare, plus the ``tenv_hvg`` fusion row. This module is pure data keyed by
``dataset_id`` and is the only place the rule document's accession and sequencing
technology columns read from (ticket #198 renders the table). It performs no I/O
and imports nothing outside the standard library.

Every value below was looked up at the primary source on ``VERIFIED_ON`` and is
never copied from a catalog's ``model_system`` prose. ``technology`` is the short
table phrase the primary record supports, ``platform`` carries the GPL id and
instrument behind that claim, and ``source_url`` is the record the row was
checked against. ``accession_note`` records a series hierarchy where one exists
and, for ``tenv_hvg``, why the row is accession-less. ``sample_size_note``
records how the repo's ``n_case``/``n_control`` and the primary record's sample
table reconcile wherever they are not directly comparable (a different counting
basis, a subset of the series, or a split only derivable from the paper); it is
``None`` where the primary count equals ``n_case + n_control``. The release's
numbers are never changed here.

Harvest method (recorded so a reviewer can re-run it):

* NCBI E-utilities ``esummary`` (``db=gds``) for every GSE and GPL uid.
* GEO FTP ``matrix/<GSE>_series_matrix.txt.gz``, read only up to the
  ``!series_matrix_table_begin`` line, for the sample count and overall design.
* GEO FTP ``suppl/`` listings and ``filelist.txt``, plus the deposits' own
  published supplementary metadata (``GSE270788_metadata.csv.gz`` for the
  ``HF.etiology`` donor arms; the DToxS ``Conv-DEG-Calc`` Calc-table headers for
  the post-QC well columns).
* The cellxgene Curation API v1 collection endpoint for the Reichart-CM object.
* NCBI PMC full-text XML (PMC9528698, Reichart et al., Science 2022) for the
  donor-cohort split.

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

_REICHART_COLLECTION = 'cellxgene collection e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb'
_REICHART_URL = 'https://cellxgene.cziscience.com/collections/e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb'

_ENTRY_FIELDS = (
    'accession', 'accession_note', 'technology', 'source_url',
    'verified_on', 'platform', 'sample_size_note',
)
_PLATFORM_FIELDS = ('gpl', 'instrument')

SOURCES = {
    # --- nHF: GSE270788, a 10x Multiome series whose GEX arm the repo uses -----
    'gse270788__nicm': {
        'accession': 'GSE270788',
        'accession_note': 'SubSeries of SuperSeries GSE218392',
        'technology': '10x Multiome snRNA-seq (GEX arm, CellRanger ARC)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270788',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'sample_size_note': (
            'The series reports 46 GSMs = 23 GEX (snRNA) + 23 ATAC libraries; the repo arms '
            'count the 23 donors instead (deposit metadata HF.etiology: Donor 7, NICM 6, '
            'ICM 5, AMI 5), which is the 6/7 donor basis of this row.'
        ),
    },
    'reichart-cm__dcm_genotype_negative': {
        'accession': _REICHART_COLLECTION,
        'accession_note': (
            'cellxgene dataset f7995301-7551-4e1d-8396-ffe3c9497ace, '
            'version 27fdbdf9-deae-4aed-bf7f-21964ee62137'
        ),
        'technology': "10x 3' snRNA-seq (CM-only cellxgene object)",
        'source_url': _REICHART_URL,
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': None, 'instrument': None},
        'sample_size_note': (
            'The collection carries 79 donors (61 failing + 18 controls); the repo arms are '
            'the 8 PV-negative DCM donors and the 18 controls stated in the paper, i.e. '
            '79 = 78 exported (8 + 44 + 8 + 18) + 1 NCCM.'
        ),
    },
    # --- ICM: the same GSE270788 GEX arm, split by aetiology -------------------
    'gse270788__icm': {
        'accession': 'GSE270788',
        'accession_note': 'SubSeries of SuperSeries GSE218392',
        'technology': '10x Multiome snRNA-seq (GEX arm, CellRanger ARC)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270788',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'sample_size_note': (
            'The series reports 46 GSMs = 23 GEX + 23 ATAC; this repo arm is the 5 ICM '
            'donors against the same 7 donor controls (HF.etiology), not the 46 libraries.'
        ),
    },
    'gse270788__ami': {
        'accession': 'GSE270788',
        'accession_note': 'SubSeries of SuperSeries GSE218392',
        'technology': '10x Multiome snRNA-seq (GEX arm, CellRanger ARC)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270788',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'sample_size_note': (
            'The series reports 46 GSMs = 23 GEX + 23 ATAC; this repo arm is the 5 AMI '
            'donors against the same 7 donor controls (HF.etiology), not the 46 libraries.'
        ),
    },
    # --- genetic_cardiomyopathy ------------------------------------------------
    'reichart-cm__dcm_genotype_positive': {
        'accession': _REICHART_COLLECTION,
        'accession_note': 'cellxgene dataset f7995301-7551-4e1d-8396-ffe3c9497ace',
        'technology': "10x 3' snRNA-seq (CM-only cellxgene object)",
        'source_url': _REICHART_URL,
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': None, 'instrument': None},
        'sample_size_note': (
            'Partial match: the repo n_case=44 is derivable (61 failing - 8 PV-negative '
            '- 8 ARVC - 1 NCCM) but no fetched primary record states 44 directly; the '
            'per-patient DCM-vs-ACM labels live only in the paper\'s Tables S1/S2.'
        ),
    },
    'reichart-cm__arvc': {
        'accession': _REICHART_COLLECTION,
        'accession_note': 'cellxgene dataset f7995301-7551-4e1d-8396-ffe3c9497ace',
        'technology': "10x 3' snRNA-seq (CM-only cellxgene object)",
        'source_url': _REICHART_URL,
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': None, 'instrument': None},
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
        'sample_size_note': None,
    },
    'gse234141__ttntv_dcm': {
        'accession': 'GSE234141',
        'accession_note': None,
        'technology': 'bulk RNA-seq (hiPSC-CM 2D and engineered heart tissue)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234141',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL28337', 'instrument': 'DNBSEQ-G50'},
        'sample_size_note': None,
    },
    'gse304575__lmna_h222p_dcm': {
        'accession': 'GSE304575',
        'accession_note': None,
        'technology': 'bulk RNA-seq (hiPSC-CM, isogenic)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE304575',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'sample_size_note': None,
    },
    'gse305648__rbm20_dcm': {
        'accession': 'GSE305648',
        'accession_note': None,
        'technology': 'bulk RNA-seq (total RNA, rRNA-depleted; hiPSC cardioids)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE305648',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL30173', 'instrument': 'NextSeq 2000'},
        'sample_size_note': None,
    },
    'gse253223__pkp2_knockdown_arvc': {
        'accession': 'GSE253223',
        'accession_note': 'SubSeries of GSE253226',
        'technology': 'bulk RNA-seq (hiPSC-CM, siRNA knockdown)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253223',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'sample_size_note': None,
    },
    'gse194297__dmd_cardiomyopathy': {
        'accession': 'GSE194297',
        'accession_note': None,
        'technology': 'bulk RNA-seq (hiPSC cardiac organoids)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194297',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL20301', 'instrument': 'Illumina HiSeq 4000'},
        'sample_size_note': None,
    },
    # --- cardiotox_TOP2 --------------------------------------------------------
    'gse243674': {
        'accession': 'GSE243674',
        'accession_note': None,
        'technology': 'bulk RNA-seq (hiPSC-CM monolayer, 3 h / 24 h)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243674',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL21697', 'instrument': 'NextSeq 550'},
        'sample_size_note': (
            'The series reports 72 GSMs = 6 donors x 6 conditions x 2 time points; the repo '
            'contrast is the 3 h slice (4 anthracyclines x 6 donors = 24 case vs 6 vehicle), '
            'not the whole series.'
        ),
    },
    'dtoxs_ipsc': {
        'accession': 'GSE217421',
        'accession_note': 'Conv arm of SuperSeries GSE217424; sibling DGE arm GSE217423',
        'technology': 'bulk RNA-seq (hiPSC-CM plate wells, 48 h)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217421',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'sample_size_note': (
            'Well-vs-GSM basis: the series carries 1,171 GSMs (4 wells per drug x 6 iPSC-CM '
            'lines, 48 h); the repo 42/41 are the post-QC well columns of the deposit\'s own '
            'Conv-DEG-Calc tables, reproduced exactly (7+2+11+9+6+7 = 42 case, '
            '7+3+10+6+8+7 = 41 control).'
        ),
    },
    'gse244740': {
        'accession': 'GSE244740',
        'accession_note': None,
        'technology': 'bulk RNA-seq (384-well ScreenSeq hiPSC-CM plate, 24 h)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244740',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'sample_size_note': (
            'The series reports 3,072 GSMs = 8 x 384-well plates (all 24 h); the repo uses '
            'Screening plate 3 (S4881Nr8) only: DOX 24 + IDA 24 + MTX 24 = 72 case and '
            'DMSO 22 control.'
        ),
    },
    # --- cardiotox_other: the three DToxS drug classes ------------------------
    'dtoxs_ipsc__atc_l4-L01EA': {
        'accession': 'GSE217421',
        'accession_note': 'Conv arm of SuperSeries GSE217424',
        'technology': 'bulk RNA-seq (hiPSC-CM plate wells, 48 h)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217421',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'sample_size_note': (
            'Well-vs-GSM basis: the series carries 1,171 GSMs; the repo 64/41 are the '
            'Calc-table well columns (12+2+16+11+8+15 = 64 case, 7+3+10+6+8+7 = 41 control) '
            'across the six lines, where MSN02-04R has only PON.'
        ),
    },
    'dtoxs_ipsc__atc_l4-L01EB': {
        'accession': 'GSE217421',
        'accession_note': 'Conv arm of SuperSeries GSE217424',
        'technology': 'bulk RNA-seq (hiPSC-CM plate wells, 48 h)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217421',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'sample_size_note': (
            'Well-vs-GSM basis: the series carries 1,171 GSMs; the repo 36/38 are the '
            'Calc-table well columns (6+0+7+6+7+10 = 36 case, 7+0+10+6+8+7 = 38 control), '
            'with n_control 38 because the MSN02-04R stratum has no L01EB table.'
        ),
    },
    'dtoxs_ipsc__atc_l4-L01EX': {
        'accession': 'GSE217421',
        'accession_note': 'Conv arm of SuperSeries GSE217424',
        'technology': 'bulk RNA-seq (hiPSC-CM plate wells, 48 h)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217421',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'},
        'sample_size_note': (
            'Well-vs-GSM basis: the series carries 1,171 GSMs; the repo 72/41 are the '
            'Calc-table well columns (12+4+20+10+14+12 = 72 case, 7+3+10+6+8+7 = 41 control).'
        ),
    },
    # --- tenv-HVG fusion row (not a DEG source) -------------------------------
    TENV_HVG: {
        'accession': None,
        'accession_note': (
            'No accession is declared in the gene axis: the tenv-HVG rank is the published '
            'descriptive per-well dispersion rank of the repo corpus built from GSE262419 '
            '(ADR 0066), so GSE262419 is the corpus source, not an axis accession.'
        ),
        'technology': 'TempO-Seq targeted transcriptomic panel (96-well, HiSeq 2500)',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262419',
        'verified_on': VERIFIED_ON,
        'platform': {'gpl': 'GPL16791', 'instrument': 'Illumina HiSeq 2500'},
        'sample_size_note': (
            'Not a contrast, so there is no n_case/n_control: the row ranks the 4,206 treated '
            'wells of the GSE262419 corpus (4,558 GSMs over 16 plates; the corpus keeps 4,510 '
            'native wells x 19,613 genes = 4,206 treated + 272 vehicle + 32 blanks).'
        ),
    },
}


def validate(dataset_ids: Iterable[str]) -> None:
    """Raise ``ValueError`` unless ``dataset_ids`` exactly matches the published coverage.

    ``dataset_ids`` is the caller's list of published dataset ids plus ``tenv_hvg``
    (the five catalogs' concatenation, in any order). It must contain no duplicate;
    its set must equal ``SOURCES``' key set — nothing missing, nothing extra; and
    every entry must be complete: an accession or a stated accession-less reason,
    a non-empty technology, source URL and verification date, and a platform dict
    with both ``gpl`` and ``instrument``.
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
        if entry['verified_on'] != VERIFIED_ON:
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
        sample_size_note = entry['sample_size_note']
        if sample_size_note is not None and not isinstance(sample_size_note, str):
            raise ValueError(f'{dataset_id}: sample_size_note must be str or None')
