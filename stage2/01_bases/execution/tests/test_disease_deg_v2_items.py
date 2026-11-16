#!/usr/bin/env python
"""New nHF tier-2 / tier-3 external items for the gene-axis v2 rebuild (#203).

Hermetic: the fixtures here are the items' *tracked* artefacts only — the
``local/deg/<GSE>_deg.json`` summaries (committed), each item's
``script/build_deg.py`` declarations (imported, never executed) and
``script/fetch.sh`` pins, plus the ``author_log.json`` records. The gitignored
``source/`` and ``local/`` payloads are not required, so the suite runs on a
fresh clone.

What is asserted is the externally visible contract of the two new sources:

* ``disease_deg_sources`` declares them in ``nHF`` with the right tier,
  evidence type, arms and contrasts, and records ``GSE185100`` as
  unverified-not-admitted (with accession, size, donor counts and reason), with
  no source record and no item directory anywhere;
* the tracked item summaries carry the *verified* series design (arms, n,
  platform, deposited scale) and agree with the declaration and with
  ``fetch.sh``'s byte pins;
* the tier-2 declaration records a non-genetic model basis (ADR 0071's clinical
  category rule), and the tier-3 declaration records that the source is
  whole-tissue, not cardiomyocyte-specific.
"""
import importlib.util
import json
import re
import sys
import unittest
from importlib.machinery import SourceFileLoader
from pathlib import Path

TESTS = Path(__file__).resolve().parent
EXECUTION = TESTS.parent
# ADR 0030: discover the root by the directory containing .git, never by counting
# parents; the sibling product tests use the same walk.
REPO = next(path for path in Path(__file__).resolve().parents if (path / ".git").exists())
sys.path.insert(0, str(EXECUTION))
import public_api  # noqa: F401  (binds the single public implementation)
import disease_deg_sources  # noqa: E402

ITEMS = ('GSE35108', 'GSE165303')
NEW_DATASETS = {
    'gse35108__dcm_ipsc_cm': ('GSE35108', 2, 4, 4,
                             ['DCM_control-Normal_control',
                              'DCM_Serca2a-Normal_Serca2a']),
    'gse165303__dcm_lv_bulk': ('GSE165303', 3, 50, 51, ['DCM-NF']),
}


def load_item_module(accession, kind='build_deg'):
    path = REPO / 'disease_data' / accession / 'script' / f'{kind}.py'
    loader = SourceFileLoader(f'v2_{accession}_{kind}', str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def read_json(path):
    return json.loads(Path(path).read_text())


def summary(accession):
    return read_json(REPO / 'disease_data' / accession / 'local/deg' /
                     f'{accession}_deg.json')


def item_dir(accession):
    return REPO / 'disease_data' / accession


def declared(accession):
    records = [record for record in disease_deg_sources.datasets('nHF')
               if record.get('source_item') == f'disease_data/{accession}']
    return records


class SourceDeclarationTest(unittest.TestCase):
    """The nHF roster gains exactly the two new sources, at the right tiers."""

    def test_new_sources_are_declared_in_nhf(self):
        for dataset_id, (_item, tier, n_case, n_control, contrasts) in NEW_DATASETS.items():
            record = next((row for row in disease_deg_sources.DATASETS
                           if row['dataset_id'] == dataset_id), None)
            self.assertIsNotNone(record, dataset_id)
            self.assertEqual(record['category'], 'nHF')
            self.assertEqual(record['tier'], tier)
            self.assertEqual(record['evidence_type'], 'bulk')
            self.assertEqual(record['n_case'], n_case)
            self.assertEqual(record['n_control'], n_control)
            self.assertEqual(record['contrasts'], contrasts)
            self.assertEqual(record['source_item'], f'disease_data/{_item}')
            self.assertEqual(record['local_deg'], f'local/deg/{_item}_deg.csv')
            self.assertEqual(record['local_summary'], f'local/deg/{_item}_deg.json')
            self.assertFalse(record['exploratory'])
            self.assertEqual(record['genotype_status'], 'not_declared')

    def test_dataset_ids_are_unique(self):
        ids = [record['dataset_id'] for record in disease_deg_sources.DATASETS]
        self.assertEqual(len(ids), len(set(ids)))

    def test_tier_two_record_declares_a_non_genetic_model_basis(self):
        (record,) = declared('GSE35108')
        basis = record['non_genetic_model_basis']
        self.assertIn('patient-derived', basis)
        self.assertIn('no mutation', basis)
        self.assertIn('no isogenic control', basis)

    def test_tier_three_record_declares_whole_tissue_bulk(self):
        (record,) = declared('GSE165303')
        self.assertIn('not cardiomyocyte-specific', record['model_basis'])
        self.assertIn('tier-3', record['model_basis'])
        self.assertNotIn('non_genetic_model_basis', record)


class Gse185100RecordingTest(unittest.TestCase):
    """The abandoned candidate is a recorded rejection, never a source."""

    def test_gse185100_has_a_measured_input_gap_and_reentry_condition(self):
        category = disease_deg_sources.category_record('nHF')
        rejected = category['rejected']
        entry = next((row for row in rejected if row.get('candidate') == 'GSE185100'), None)
        self.assertIsNotNone(entry, 'GSE185100 missing from the nHF rejected list')
        self.assertEqual(entry['status'], 'author_cm_input_unavailable')
        self.assertIn('no author per-cell CM assignment', entry['reason'])
        self.assertIn('complete unrounded test table', entry['reason'])
        self.assertIn('age and low estimated yield are not exclusion grounds', entry['reason'])
        self.assertIn('CM barcode assignments with counts', entry['reentry'])
        self.assertIn('complete author CM DCM-vs-ND', entry['reentry'])
        review = next(row for row in category['source_review']['records']
                      if row['candidate'] == 'GSE185100')
        self.assertEqual(review['disposition'], 'not_admitted_author_cm_input_unavailable')
        self.assertEqual(review['reason'], entry['reason'])
        self.assertEqual(review['reentry'], entry['reentry'])
        self.assertEqual(len(review['sources']), 4)
        self.assertIn('https://doi.org/10.1161/CIRCULATIONAHA.121.057763', review['sources'])
        self.assertIn('Author objects may be obtainable elsewhere', review['uncertainty'])

    def test_gse185100_is_not_a_source_or_an_item(self):
        self.assertFalse(any('GSE185100' in record.get('source_item', '')
                             for record in disease_deg_sources.DATASETS))
        self.assertFalse(any('GSE185100' in record['dataset_id']
                             for record in disease_deg_sources.DATASETS))
        self.assertFalse((REPO / 'disease_data' / 'GSE185100').exists())


class ItemSummaryTest(unittest.TestCase):
    """Each tracked summary states the design verified against the deposit."""

    def test_common_summary_contract(self):
        for accession in ITEMS:
            record = summary(accession)
            self.assertEqual(record['schema'], 'stage2-item-deg-summary-v1')
            self.assertEqual(record['item'], accession)
            self.assertEqual(record['category'], 'nHF')
            self.assertEqual(record['evidence_type'], 'bulk')
            self.assertEqual(record['design']['kind'], 'expression')
            self.assertEqual(record['design']['formula'], '~0 + group')
            self.assertEqual(sorted(record['contrasts']),
                             sorted(declared(accession)[0]['contrasts']))
            self.assertEqual(record['method']['implementation'],
                             'stage2/01_bases/output/api/bases/disease_deg.R')
            self.assertFalse(record['deg']['path'].startswith('/'))
            self.assertFalse(record['expr']['path'].startswith('/'))
            self.assertFalse(record['samples']['path'].startswith('/'))
            for key in ('deg', 'expr', 'samples'):
                self.assertTrue(record[key]['path'].startswith('local/deg/'))
                self.assertEqual(len(record[key]['sha256']), 64)

    def test_gse35108_verified_design(self):
        record = summary('GSE35108')
        series = record['series']
        self.assertEqual(record['tier'], 2)
        self.assertEqual(series['total_libraries'], 17)
        self.assertEqual(series['cm_libraries'], 8)
        self.assertEqual(series['excluded_libraries'], 9)
        self.assertEqual(series['platform'], 'GPL6244')
        self.assertTrue(record['deposition']['values']['one_colour'])
        self.assertIn('RMA', record['deposition']['values']['normalisation'])
        self.assertEqual(record['samples']['n'], 8)
        self.assertEqual(record['samples']['group_counts'],
                         {'DCM_Serca2a': 2, 'DCM_control': 2,
                          'Normal_Serca2a': 2, 'Normal_control': 2})
        self.assertEqual(record['n_case'], 4)
        self.assertEqual(record['n_control'], 4)
        for label, arm in record['contrasts'].items():
            self.assertEqual((arm['n_case'], arm['n_control']), (2, 2), label)
        collapse = record['probe_collapse']
        self.assertEqual(collapse['n_probes_with_single_symbol']
                         + collapse['n_probes_multi_or_missing_symbol'],
                         collapse['n_probes'])
        self.assertEqual(collapse['n_probes_with_single_symbol']
                         - collapse['n_probes_dropped_by_collapse'],
                         collapse['n_genes'])
        counts = record['deg_counts']['collapsed']
        self.assertGreater(counts['n_deg'], 0)
        self.assertEqual(counts['n_up'] + counts['n_down'], counts['n_deg'])
        self.assertEqual(record['deg']['rows'],
                         record['deg']['n_source_genes'] * 2)
        self.assertEqual(record['sensitivity']['contrast'], 'DCM_control-Normal_control')
        self.assertEqual((record['sensitivity']['n_case'],
                          record['sensitivity']['n_control']), (2, 2))

    def test_gse165303_verified_design(self):
        record = summary('GSE165303')
        series = record['series']
        self.assertEqual(record['tier'], 3)
        self.assertEqual(series['rna_libraries'], 101)
        self.assertEqual(series['non_rna_libraries'], 46)
        self.assertEqual(series['total_libraries_in_series'], 177)
        self.assertEqual(series['total_libraries_on_gpl20795'], 147)
        self.assertIn('left-ventricle', series['tissue'])
        self.assertEqual(record['samples']['n'], 101)
        self.assertEqual(record['samples']['group_counts'], {'DCM': 50, 'NF': 51})
        self.assertEqual(record['n_case'], 50)
        self.assertEqual(record['n_control'], 51)
        values = record['deposition']['values']
        self.assertEqual(values['floor_value'], 3.821687527)
        self.assertGreater(values['high_expression_correlation'], 0.99)
        self.assertGreater(values['floor_fraction'], 0.3)
        self.assertEqual(record['deposition']['gene_ids']['n_rows'],
                         record['deg']['n_source_genes'])
        counts = record['deg_counts']['collapsed']
        self.assertGreater(counts['n_deg'], 0)
        per_contrast = record['deg_counts']['per_contrast']['DCM-NF']
        self.assertEqual(per_contrast['n_up'] + per_contrast['n_down'],
                         per_contrast['n_deg'])
        self.assertEqual(record['deg']['rows'], record['deg']['n_source_genes'])
        # The authors' own DESeq2 columns are carried only as a cross-check.
        check = record['author_table_check']
        self.assertEqual(check['n_rows_joined'], record['deg']['rows'])
        self.assertGreater(check['log2fc_spearman'], 0.9)
        self.assertGreater(check['author_n_v2'], check['limma_n_v2'])


class BuildSpecAndPinTest(unittest.TestCase):
    """The build declaration and fetch pins agree with the tracked summary."""

    def test_spec_matches_summary(self):
        for accession in ITEMS:
            spec = load_item_module(accession).SPEC
            record = summary(accession)
            self.assertEqual(spec['item'], accession)
            self.assertEqual(spec['category'], 'nHF')
            self.assertEqual(spec['tier'], record['tier'])
            self.assertEqual(spec['kind'], 'expression')
            self.assertEqual(spec['formula'], '~0 + group')
            self.assertEqual(sorted(contrast['label'] for contrast in spec['contrasts']),
                             sorted(record['contrasts']))
            deposition = record['deposition']
            # GSE35108's deposited matrix *is* the series matrix; GSE165303 keeps
            # the master table separate from the series matrix's sample block.
            expected = (deposition['series_matrix']['path']
                        if 'series_matrix' in deposition
                        else deposition['matrix']['path'])
            self.assertEqual(spec['series_matrix'], expected)

    def test_fetch_pins_match_the_deposited_bytes(self):
        for accession in ITEMS:
            record = summary(accession)
            script = (item_dir(accession) / 'script/fetch.sh').read_text()
            pins = {}
            for line in script.splitlines():
                match = re.search(r'/([A-Za-z0-9._-]+)\t\1\t(\d+)', line)
                if match:
                    pins[match.group(1)] = int(match.group(2))
            deposition = record['deposition']
            for key in ('series_matrix', 'master_table', 'platform_annot'):
                if key not in deposition:
                    continue
                entry = deposition[key]
                name = Path(entry['path']).name
                self.assertIn(name, pins, f'{accession}: {key} not pinned')
                self.assertEqual(pins[name], entry['bytes'], f'{accession}: {key} bytes')
            self.assertGreaterEqual(len(pins), 2)

    def test_author_log_pins_the_source_tree(self):
        for accession in ITEMS:
            log = read_json(item_dir(accession) / 'author_log.json')
            self.assertEqual(log['category'], 'disease_data')
            self.assertEqual(log['manifest'], 'manifest.sha256')
            self.assertEqual(log['verify_result'], 'pass')
            (source,) = log['sources']
            self.assertEqual(source['id'], accession)
            self.assertEqual(source['files'], 2)
            self.assertEqual(len(source['tree_sha256']), 64)
            declared_paths = {entry['path'] for entry in log['local_additions']}
            self.assertIn('.gitignore', declared_paths)
            record = summary(accession)
            for key in ('deg', 'expr', 'samples'):
                self.assertIn(record[key]['path'], declared_paths)
            # No absolute path in the tracked provenance record.
            self.assertNotRegex(json.dumps(log), r'"/(home|mnt|tmp)/')


if __name__ == '__main__':
    unittest.main()
