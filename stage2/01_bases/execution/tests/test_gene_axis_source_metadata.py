"""External-contract tests for the verified gene-axis source metadata (#197).

Only the mapping's published constants and its ``validate`` helper are exercised:
the coverage contract against the five released disease-deg catalogs, the
required fields of every entry, and a handful of primary-source values a
reviewer can re-verify against issue #197's comment. Nothing internal is
asserted, and the module is pure data (it does no I/O), so the catalogs are read
here, from disk, the same way the docx generator will read them.
"""
import copy
import json
import sys
import unittest
from pathlib import Path
from unittest import mock

EXECUTION = Path(__file__).resolve().parents[1]
STEP = EXECUTION.parent
sys.path.insert(0, str(EXECUTION))

import gene_axis_source_metadata as metadata

CATEGORIES = ('nhf', 'icm', 'genetic_cardiomyopathy', 'cardiotox_TOP2', 'cardiotox_other')
DISEASE_DEG = STEP / 'output' / 'disease-deg'
TENV = 'tenv_hvg'
ENTRY_FIELDS = {
    'accession', 'accession_note', 'technology', 'source_url',
    'verified_on', 'platform', 'platform_note', 'sample_size_note',
}
REICHART_COLLECTION_ID = 'e75342a8-0f3b-4ec5-8ee1-245a23e0f7cb'


def published_dataset_ids():
    """The dataset ids the five released catalogs publish, in catalog order."""
    dataset_ids = []
    for category in CATEGORIES:
        catalog = json.loads((DISEASE_DEG / category / 'catalog.json').read_text(encoding='utf-8'))
        dataset_ids.extend(record['dataset_id'] for record in catalog['datasets'])
    return dataset_ids


class CoverageContractTest(unittest.TestCase):
    def test_validate_passes_for_every_published_dataset_plus_tenv(self):
        metadata.validate(published_dataset_ids() + [TENV])

    def test_key_set_is_exactly_the_published_datasets_plus_tenv(self):
        published = published_dataset_ids()
        self.assertEqual(len(published), 18)
        self.assertEqual(set(metadata.SOURCES), set(published) | {TENV})
        self.assertEqual(len(metadata.SOURCES), len(set(published)) + 1)

    def test_validate_rejects_a_missing_dataset_id(self):
        identifiers = published_dataset_ids() + [TENV]
        for missing in (published_dataset_ids()[0], TENV):
            with self.subTest(missing=missing):
                with self.assertRaises(ValueError):
                    metadata.validate([i for i in identifiers if i != missing])

    def test_validate_rejects_an_extra_dataset_id(self):
        with self.assertRaises(ValueError):
            metadata.validate(published_dataset_ids() + [TENV, 'gse999999__not_released'])

    def test_validate_rejects_a_duplicate_dataset_id(self):
        identifiers = published_dataset_ids() + [TENV]
        with self.assertRaises(ValueError):
            metadata.validate(identifiers + [published_dataset_ids()[0]])

    def test_validate_rejects_an_entry_with_a_blank_required_field(self):
        identifiers = published_dataset_ids() + [TENV]
        broken = copy.deepcopy(metadata.SOURCES)
        broken['gse243674']['technology'] = '   '
        with mock.patch.dict(metadata.SOURCES, broken, clear=True):
            with self.assertRaises(ValueError):
                metadata.validate(identifiers)


class EntryContractTest(unittest.TestCase):
    def test_schema_and_verification_date_are_the_declared_constants(self):
        self.assertEqual(metadata.SCHEMA, 'stage2-gene-axis-source-metadata-v1')
        self.assertEqual(metadata.VERIFIED_ON, '2026-11-15')

    def test_every_entry_carries_the_required_fields(self):
        for dataset_id, entry in metadata.SOURCES.items():
            with self.subTest(dataset_id=dataset_id):
                self.assertEqual(set(entry), ENTRY_FIELDS)
                self.assertEqual(set(entry['platform']), {'gpl', 'instrument'})
                self.assertEqual(entry['verified_on'], metadata.VERIFIED_ON)
                accession = (entry['accession'] or '').strip()
                accession_note = (entry['accession_note'] or '').strip()
                self.assertTrue(accession or accession_note)
                for field in ('technology', 'source_url', 'platform_note'):
                    self.assertIsInstance(entry[field], str)
                    self.assertTrue(entry[field].strip())
                self.assertTrue(entry['source_url'].startswith('https://'))
                self.assertTrue(entry['sample_size_note'] is None
                                or (isinstance(entry['sample_size_note'], str)
                                    and entry['sample_size_note'].strip()))

    def test_platform_note_backs_or_explains_the_platform_evidence(self):
        for dataset_id, entry in metadata.SOURCES.items():
            with self.subTest(dataset_id=dataset_id):
                platform, note = entry['platform'], entry['platform_note']
                self.assertEqual(bool(platform['gpl']), bool(platform['instrument']))
                if platform['gpl'] is None:
                    self.assertIn('no GEO platform', note)
                else:
                    self.assertIn(platform['gpl'], note)
                    self.assertIn(platform['instrument'], note)

    def test_sample_size_notes_are_non_empty_when_present(self):
        notes = [entry['sample_size_note'] for entry in metadata.SOURCES.values()
                 if entry['sample_size_note'] is not None]
        self.assertEqual(len(notes), 13, 'expected thirteen reconciled datasets')
        for note in notes:
            self.assertIsInstance(note, str)
            self.assertTrue(note.strip())

    def test_tenv_sample_size_constants_match_the_published_view_manifest(self):
        manifest = json.loads((STEP / 'output/view-2k/manifest.json').read_text(
            encoding='utf-8'))
        tenv = manifest['sources']['tenv']
        self.assertEqual(metadata.TENV_TREATED_WELLS, tenv['treated'])
        self.assertEqual(metadata.TENV_VEHICLE_WELLS, tenv['vehicle'])
        self.assertEqual(metadata.TENV_TREATED_WELLS + metadata.TENV_VEHICLE_WELLS,
                         tenv['rows'])

    def test_tenv_is_explicitly_accession_less_with_a_stated_reason(self):
        entry = metadata.SOURCES[TENV]
        self.assertIsNone(entry['accession'])
        self.assertIsInstance(entry['accession_note'], str)
        self.assertIn('GSE262419', entry['accession_note'])
        self.assertLessEqual(len(entry['accession_note']), 60,
                             'the table cell stays short; the rationale lives in the docstring')

    def test_tenv_is_the_only_entry_without_an_accession(self):
        accession_less = [dataset_id for dataset_id, entry in metadata.SOURCES.items()
                          if not entry['accession']]
        self.assertEqual(accession_less, [TENV])


class PrimarySourceSpotCheckTest(unittest.TestCase):
    """A handful of published values, re-checkable against issue #197's comment."""

    def test_gse234141_is_dnbseq_g50_not_illumina(self):
        entry = metadata.SOURCES['gse234141__ttntv_dcm']
        self.assertEqual(entry['accession'], 'GSE234141')
        self.assertEqual(entry['platform'], {'gpl': 'GPL28337', 'instrument': 'DNBSEQ-G50'})
        self.assertNotIn('Illumina', entry['platform']['instrument'])

    def test_gse243674_is_nextseq_550(self):
        entry = metadata.SOURCES['gse243674']
        self.assertEqual(entry['accession'], 'GSE243674')
        self.assertEqual(entry['platform'], {'gpl': 'GPL21697', 'instrument': 'NextSeq 550'})

    def test_gse270788_is_10x_multiome_on_novaseq_6000(self):
        for dataset_id in ('gse270788__nicm', 'gse270788__icm', 'gse270788__ami'):
            with self.subTest(dataset_id=dataset_id):
                entry = metadata.SOURCES[dataset_id]
                self.assertEqual(entry['accession'], 'GSE270788')
                self.assertIn('10x Multiome', entry['technology'])
                self.assertEqual(entry['platform'],
                                 {'gpl': 'GPL24676', 'instrument': 'Illumina NovaSeq 6000'})

    def test_gse305648_is_nextseq_2000_total_rna(self):
        entry = metadata.SOURCES['gse305648__rbm20_dcm']
        self.assertEqual(entry['platform'], {'gpl': 'GPL30173', 'instrument': 'NextSeq 2000'})
        self.assertIn('rRNA', entry['technology'])

    def test_reichart_is_a_cellxgene_collection_not_a_geo_accession(self):
        for dataset_id in ('reichart-cm__dcm_genotype_negative',
                           'reichart-cm__dcm_genotype_positive',
                           'reichart-cm__arvc'):
            with self.subTest(dataset_id=dataset_id):
                entry = metadata.SOURCES[dataset_id]
                self.assertTrue(entry['accession'].startswith('cellxgene collection '))
                self.assertNotIn('GSE', entry['accession'])
                self.assertIn(REICHART_COLLECTION_ID, entry['accession'])
                self.assertIn(REICHART_COLLECTION_ID, entry['source_url'])
                self.assertIsNone(entry['platform']['gpl'])
                self.assertIsNone(entry['platform']['instrument'])
                self.assertIn('no GEO platform', entry['platform_note'])
                self.assertIn('do not name the instrument', entry['platform_note'])

    def test_tenv_source_is_the_gse262419_corpus(self):
        entry = metadata.SOURCES[TENV]
        self.assertIn('GSE262419', entry['source_url'])
        self.assertEqual(entry['platform'],
                         {'gpl': 'GPL16791', 'instrument': 'Illumina HiSeq 2500'})


class ValidateFailureTest(unittest.TestCase):
    """The validate failures a reviewer can re-trigger, including their wording."""

    def test_coverage_mismatch_error_names_the_schema(self):
        with self.assertRaises(ValueError) as caught:
            metadata.validate(published_dataset_ids() + [TENV, 'gse999999__not_released'])
        self.assertIn(metadata.SCHEMA, str(caught.exception))

    def test_validate_rejects_a_blank_platform_note(self):
        broken = copy.deepcopy(metadata.SOURCES)
        broken['gse243674']['platform_note'] = '   '
        with mock.patch.dict(metadata.SOURCES, broken, clear=True):
            with self.assertRaises(ValueError):
                metadata.validate(published_dataset_ids() + [TENV])

    def test_validate_rejects_a_null_platform_without_an_explanation(self):
        broken = copy.deepcopy(metadata.SOURCES)
        broken['gse243674']['platform'] = {'gpl': None, 'instrument': None}
        broken['gse243674']['platform_note'] = ''
        with mock.patch.dict(metadata.SOURCES, broken, clear=True):
            with self.assertRaises(ValueError):
                metadata.validate(published_dataset_ids() + [TENV])


if __name__ == '__main__':
    unittest.main()
