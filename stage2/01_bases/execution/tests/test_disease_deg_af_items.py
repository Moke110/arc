#!/usr/bin/env python
"""New AF external items for the gene-axis v2 rebuild (#204).

Hermetic: the fixtures are the items' *tracked* artefacts only — the
``local/deg/<GSE>_deg.json`` summaries and the GSE238242 pseudobulk summary
(committed), each item's ``script/build_deg.py`` / ``build_pseudobulk.py``
declarations (imported, never executed), ``script/fetch.sh`` byte pins and the
``author_log.json`` records. The gitignored ``source/`` and ``local/`` payloads
are not required, so the suite runs on a fresh clone.

What is asserted is the externally visible contract of the three new AF sources
and their declaration under ADR 0071:

* ``disease_deg_sources`` declares the returning ``AF`` category with the four
  sources at their tiers (GSE255612/GSE138252/GSE238242 tier 1, GSE269882 tier 3),
  records GSE271229 as audited-and-not-admitted and GSE128188 as an owned
  tier-3 candidate (not a dataset of this ticket);
* each tracked summary carries the *verified* deposit design (accession arms, n,
  platform, assay), agrees with the declaration, and reports the measured v2 DEG
  count under the one live implementation — including the two honest zeros;
* GSE138252 records the 48-GSM / 38-column discrepancy, the LA + 2020-batch
  chamber model, the AF+HF confound and the TMM alternative as a sensitivity;
* GSE238242 records the snRNA-only arm, the author CM gate and the count
  pseudobulk; GSE269882 records the tier-3 admission, the staged arms and the
  patient-number join verified with sex probes.
"""
import hashlib
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

ITEMS = ('GSE138252', 'GSE238242', 'GSE269882')
NEW_DATASETS = {
    'gse138252__af': ('GSE138252', 1, 'bulk', 5, 5, ['AF-CTL']),
    'gse238242__af': ('GSE238242', 1, 'cm_pseudobulk', 3, 4, ['AF-SR']),
    'gse269882__af': ('GSE269882', 3, 'bulk', 18, 9,
                      ['persAF-SR', 'pAF-SR', 'poAF-SR']),
}


def load_item_module(accession, kind='build_deg'):
    path = REPO / 'disease_data' / accession / 'script' / f'{kind}.py'
    loader = SourceFileLoader(f'af_{accession}_{kind}', str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def read_json(path):
    return json.loads(Path(path).read_text())


def summary(accession):
    return read_json(REPO / 'disease_data' / accession / 'local/deg' /
                     f'{accession}_deg.json')


def pseudobulk_summary():
    return read_json(REPO / 'disease_data/GSE238242/local/pseudobulk' /
                     'GSE238242_cm_pseudobulk.json')


def item_dir(accession):
    return REPO / 'disease_data' / accession


def declared(accession):
    return [record for record in disease_deg_sources.datasets('AF')
            if record.get('source_item') == f'disease_data/{accession}']


class AfDeclarationTest(unittest.TestCase):
    """AF returns as a declared category with its four tiered sources."""

    def test_af_category_is_declared_again(self):
        record = disease_deg_sources.category_record('AF')
        self.assertEqual(record['directory'], 'af')
        self.assertEqual(record['tier'], 1)
        self.assertEqual(record['evidence_type'], 'cm_pseudobulk')
        for group in ('gse255612:donor', 'gse138252:control_la', 'gse238242:sr',
                      'gse269882:sinus_rhythm'):
            self.assertIn(group, record['shared_control_group'])

    def test_new_sources_are_declared_with_tier_and_contrasts(self):
        for dataset_id, (item, tier, evidence, n_case, n_control, contrasts) \
                in NEW_DATASETS.items():
            record = next((row for row in disease_deg_sources.DATASETS
                           if row['dataset_id'] == dataset_id), None)
            self.assertIsNotNone(record, dataset_id)
            self.assertEqual(record['category'], 'AF')
            self.assertEqual(record['tier'], tier, dataset_id)
            self.assertEqual(record['evidence_type'], evidence, dataset_id)
            self.assertEqual(record['n_case'], n_case, dataset_id)
            self.assertEqual(record['n_control'], n_control, dataset_id)
            self.assertEqual(record['contrasts'], contrasts, dataset_id)
            self.assertEqual(record['source_item'], f'disease_data/{item}')
            self.assertEqual(record['local_deg'], f'local/deg/{item}_deg.csv')
            self.assertEqual(record['local_summary'], f'local/deg/{item}_deg.json')
            self.assertFalse(record['exploratory'])
            self.assertEqual(record['genotype_status'], 'not_declared')

    def test_gse255612_record_is_restored_unchanged(self):
        (record,) = declared('GSE255612')
        self.assertEqual(record['category'], 'AF')
        self.assertEqual(record['tier'], 1)
        self.assertEqual(record['evidence_type'], 'cm_pseudobulk')
        self.assertEqual((record['n_case'], record['n_control']), (18, 16))
        self.assertEqual(record['contrasts'], ['AF-control'])
        self.assertEqual(record['shared_control_group'], 'gse255612:donor')

    def test_gse271229_is_recorded_not_admitted(self):
        rejected = disease_deg_sources.category_record('AF')['rejected']
        entry = next((row for row in rejected
                      if row.get('candidate') == 'GSE271229'), None)
        self.assertIsNotNone(entry, 'GSE271229 missing from the AF rejected list')
        self.assertEqual(entry['accession'], 'GSE271229')
        self.assertEqual(entry['status'], 'not_admitted')
        self.assertIn('author cell-type annotation', entry['reason'])
        self.assertIn('invented from markers', entry['reason'])
        self.assertIn('different accession from GSE269882', entry['reason'])
        # It is not an item (no directory) and not a dataset of this ticket.
        self.assertNotIn('GSE271229',
                         {path.name for path in (REPO / 'disease_data').iterdir()
                          if path.is_dir()})
        self.assertFalse(any('GSE271229' in row['source_item']
                             for row in disease_deg_sources.DATASETS))
        self.assertFalse(any('GSE271229' in row['dataset_id']
                             for row in disease_deg_sources.DATASETS))

    def test_gse128188_is_recorded_as_an_owned_tier3_candidate(self):
        record = disease_deg_sources.category_record('AF')
        entry = next((row for row in record['tier3_owned']
                      if row['candidate'] == 'GSE128188'), None)
        self.assertIsNotNone(entry, 'GSE128188 missing from AF tier3_owned')
        self.assertEqual(entry['item'], 'disease_data/GSE128188')
        self.assertIn('candidate', entry['reason'])
        self.assertFalse(any(row['source_item'] == 'disease_data/GSE128188'
                             for row in disease_deg_sources.DATASETS
                             if row['category'] == 'AF'))

    def test_tiering_record_states_the_measured_union(self):
        tiering = disease_deg_sources.category_record('AF')['tiering']
        self.assertEqual(tiering['tier_1']['measured_union_v2_deg'], 1)
        self.assertEqual(tiering['tier_1']['measured_union_genes'], ['SLC9A7'])
        self.assertEqual(tiering['tier_3']['measured_union_v2_deg'], 0)
        # The ADR's 2-gene figure is flagged as the TMM audit number, not ours.
        self.assertIn('TMM', tiering['tier_1']['note'])
        self.assertIn('live worker', tiering['stop_decision'])

    def test_tier_three_record_declares_whole_tissue_bulk(self):
        (record,) = declared('GSE269882')
        self.assertIn('not cardiomyocyte-specific', record['model_basis'])
        self.assertIn('tier-3', record['model_basis'])
        self.assertEqual(record['tier'], 3)

    def test_dataset_ids_are_unique(self):
        ids = [record['dataset_id'] for record in disease_deg_sources.DATASETS]
        self.assertEqual(len(ids), len(set(ids)))

    def test_category_matches_every_item_summary(self):
        for accession in ITEMS:
            record = summary(accession)
            (declaration,) = declared(accession)
            self.assertEqual(record['tier'], declaration['tier'])
            self.assertEqual(record['evidence_type'], declaration['evidence_type'])
            self.assertEqual(sorted(record['contrasts']), sorted(declaration['contrasts']))


class ItemSummaryTest(unittest.TestCase):
    """Each tracked summary states the design verified against the deposit."""

    def test_common_summary_contract(self):
        for accession in ITEMS:
            record = summary(accession)
            self.assertEqual(record['schema'], 'stage2-item-deg-summary-v1')
            self.assertEqual(record['item'], accession)
            self.assertEqual(record['category'], 'AF')
            self.assertEqual(record['design']['kind'], 'counts')
            self.assertEqual(record['design']['formula'], '~0 + group')
            self.assertEqual(record['design']['factor'], 'group')
            labels = sorted(row['label'] for row in record['design']['contrasts'])
            self.assertEqual(sorted(record['contrasts']), labels, accession)
            self.assertEqual(record['method']['implementation'],
                             'stage2/01_bases/output/api/bases/disease_deg.R')
            self.assertEqual(len(record['method']['implementation_sha256']), 64)
            for key in ('deg', 'expr'):
                self.assertFalse(record[key]['path'].startswith('/'))
                self.assertTrue(record[key]['path'].startswith('local/'))
                self.assertEqual(len(record[key]['sha256']), 64)
            self.assertFalse(record['samples']['path'].startswith('/'))
            self.assertEqual(len(record['samples']['sha256']), 64)

    def test_summary_pins_the_live_implementation_it_ran(self):
        worker = REPO / 'stage2/01_bases/output/api/bases/disease_deg.R'
        digest = hashlib.sha256(worker.read_bytes()).hexdigest()
        for accession in ITEMS:
            self.assertEqual(summary(accession)['method']['implementation_sha256'],
                             digest, accession)

    def test_every_local_output_is_declared(self):
        for accession in ITEMS:
            declared_paths = {entry['path'] for entry in
                              read_json(item_dir(accession) / 'author_log.json')[
                                  'local_additions']}
            record = summary(accession)
            self.assertIn(record['deg']['path'], declared_paths, accession)
            self.assertIn(record['expr']['path'], declared_paths, accession)
            self.assertIn(record['samples']['path'], declared_paths, accession)

    def test_gse138252_verified_design(self):
        record = summary('GSE138252')
        series = record['series']
        self.assertEqual(record['tier'], 1)
        self.assertEqual(series['total_gsms'], 48)
        self.assertEqual(series['matrix_library_columns'],
                         {'primary_2020': 18, 'addendum_2021': 20,
                          'technical_replicate_subset': 6, 'distinct_columns_total': 38})
        self.assertEqual(series['primary_atrial_arms'],
                         {'LA_AF': 5, 'LA_CTL': 5, 'RA_AF': 5, 'RA_CTL': 3})
        self.assertIn('GSM4103429', series['gsm_versus_column_discrepancy'])
        self.assertIn('re-sequences four 2020 hearts',
                      series['gsm_versus_column_discrepancy'])
        self.assertEqual(record['samples']['n'], 10)
        self.assertEqual(record['samples']['region'], 'left atrium')
        self.assertEqual(record['samples']['batch'], '2020')
        self.assertEqual((record['n_case'], record['n_control']), (5, 5))
        counts = record['deg_counts']
        self.assertEqual(counts['collapsed']['n_deg'], 1)
        self.assertEqual(counts['collapsed']['deg_genes'], ['SLC9A7'])
        self.assertEqual(counts['per_contrast']['AF-CTL']['n_deg'], 1)
        self.assertEqual(counts['per_contrast']['AF-CTL']['n_up'], 1)
        self.assertEqual(counts['per_contrast']['AF-CTL']['n_down'], 0)
        # The chamber/arm decision and the AF+HF confound are recorded.
        decision = record['chamber_and_arm_decision']
        self.assertIn('left atrium', decision['primary'])
        self.assertIn('2020', decision['primary'])
        self.assertIn('confounded', record['confound']['statement'])
        self.assertEqual(record['confound']['af_arm'].split(' + ')[0], 'AF')
        # HF alone is at least as DEG-rich as AF+HF, which is the confound check.
        self.assertGreater(record['confound']['hf_only_diagnostic']['n_deg'],
                           counts['collapsed']['n_deg'])

    def test_gse138252_column_mapping_is_a_bijection(self):
        mapping = summary('GSE138252')['deposition']['column_to_gsm_mapping']
        self.assertEqual(len(mapping), 18)
        self.assertEqual(len(set(mapping.values())), 18)
        for gsm in mapping.values():
            self.assertRegex(gsm, r'^GSM41034(2[9]|3\d|4\d|5[0-6])$')
        self.assertNotIn('GSM4103429', mapping.values())

    def test_gse138252_records_the_tmm_alternative_as_a_sensitivity(self):
        sensitivity = summary('GSE138252')['sensitivity']
        self.assertEqual(sensitivity['audit_model_alternative_tmm']['status'],
                         'not_published')
        self.assertEqual(sensitivity['audit_model_alternative_tmm']['n_deg_v2'], 2)
        self.assertIn('LRRC10',
                      sensitivity['audit_model_alternative_tmm']['deg_genes'])
        self.assertEqual(summary('GSE138252')['deg_counts']['collapsed']['n_deg'], 1)
        self.assertEqual(sensitivity['pooled_controls_batch_adjusted']['n_control'], 10)

    def test_gse238242_verified_design_and_honest_zero(self):
        record = summary('GSE238242')
        series = record['series']
        self.assertEqual(record['tier'], 1)
        self.assertEqual(series['total_libraries'], 14)
        self.assertEqual(series['rna_libraries'], 7)
        self.assertEqual(series['atac_libraries_not_fetched'], 7)
        self.assertIn('never downloaded', series['excluded_libraries_reason'])
        self.assertEqual(record['samples']['n'], 7)
        self.assertEqual(record['samples']['group_counts'], {'AF': 3, 'SR': 4})
        self.assertEqual((record['n_case'], record['n_control']), (3, 4))
        annotation = series['author_cell_annotation']
        self.assertEqual(pseudobulk_summary()['cell_type_label'], 'CM')
        self.assertEqual(annotation['nuclei'], 11986)
        self.assertEqual(annotation['cm_nuclei'], 3049)
        self.assertEqual(record['deg_counts']['per_contrast']['AF-SR']['n_deg'], 0)
        self.assertEqual(record['deg_counts']['collapsed']['n_deg'], 0)
        self.assertGreater(record['deg_counts']['per_contrast']['AF-SR']['min_padj'],
                           0.05)
        self.assertIn('zero', record['null_outcome_statement'])

    def test_gse238242_pseudobulk_summary_agrees_with_the_deg_summary(self):
        pseudobulk = pseudobulk_summary()
        record = summary('GSE238242')
        self.assertEqual(pseudobulk['schema'], 'stage2-item-cm-pseudobulk-v1')
        self.assertEqual(pseudobulk['cell_type_label'], 'CM')
        self.assertEqual(pseudobulk['pseudobulk']['genes'], 38231)
        self.assertEqual(pseudobulk['pseudobulk']['profiles'], 7)
        self.assertEqual(pseudobulk['rhythm_counts'], {'AF': 3, 'SR': 4})
        self.assertEqual(record['pseudobulk']['sha256'],
                         pseudobulk['pseudobulk']['matrix_sha256'])
        self.assertEqual(record['pseudobulk']['profiles_sha256'],
                         pseudobulk['pseudobulk']['profiles_sha256'])

    def test_gse269882_verified_design_and_staging(self):
        record = summary('GSE269882')
        series = record['series']
        self.assertEqual(record['tier'], 3)
        self.assertEqual(series['total_libraries'], 27)
        self.assertEqual(series['arms'], {'SR': 9, 'poAF': 8, 'pAF': 5, 'persAF': 5})
        self.assertEqual(series['arms_as_verified'],
                         {'sinus rhythm': 9, 'post-operative AF': 8,
                          'paroxysmal AF': 5, 'persistent AF': 5})
        self.assertIn('10', series['arm_prose_discrepancy'])
        self.assertEqual(len(record['design']['contrasts']), 3)
        self.assertEqual(sorted(record['contrasts']),
                         ['pAF-SR', 'persAF-SR', 'poAF-SR'])
        self.assertEqual((record['n_case'], record['n_control']), (18, 9))
        for label, contrast in record['contrasts'].items():
            self.assertEqual(contrast['n_control'], 9, label)
        self.assertEqual(record['deg_counts']['collapsed']['n_deg'], 0)
        for label, arm in record['deg_counts']['per_contrast'].items():
            self.assertEqual(arm['n_deg'], 0, label)
            self.assertGreater(arm['min_padj'], 0.05, label)
        self.assertIn('tier 3', record['tier_note'])
        self.assertIn('GSE271229', record['tier_note'])

    def test_gse269882_column_join_is_verified_with_sex_probes(self):
        verification = summary('GSE269882')['series']['sex_label_verification']
        self.assertEqual(verification['n_columns'], 27)
        self.assertEqual(verification['n_agree'], 26)
        (discordant,) = verification['discordant']
        self.assertEqual(discordant['patient'], '1409')
        self.assertEqual(discordant['annotated'], 'female')
        self.assertEqual(discordant['expression'], 'male')
        cross = summary('GSE269882')['deposition']['tpm']['cross_check']
        self.assertGreater(cross['vs_count_limma_log2fc_persAF_spearman'], 0.8)
        self.assertGreater(cross['vs_count_cpm_delta_spearman'], 0.8)
        pooled = summary('GSE269882')['sensitivity']['pooled_AF_vs_SR']
        self.assertEqual(pooled['n_deg'], 0)


class BuildSpecAndPinTest(unittest.TestCase):
    """The build declaration and fetch pins agree with the tracked summary."""

    def test_spec_matches_summary(self):
        for accession in ITEMS:
            spec = load_item_module(accession).SPEC
            record = summary(accession)
            self.assertEqual(spec['item'], accession)
            self.assertEqual(spec['category'], 'AF')
            self.assertEqual((spec['tier'], spec['evidence_type']),
                             (record['tier'], record['evidence_type']))
            self.assertEqual(sorted(contrast['label'] for contrast in spec['contrasts']),
                             sorted(row['label']
                                    for row in record['design']['contrasts']))

    def test_gse238242_pseudobulk_spec(self):
        module = load_item_module('GSE238242', 'build_pseudobulk')
        self.assertEqual(module.CM_LABEL, 'CM')
        self.assertEqual(module.EXPECTED_CM_NUCLEI, 3049)
        self.assertEqual(module.EXPECTED_NUCLEI, 11986)
        self.assertEqual(module.EXPECTED_GENES, 38231)
        self.assertEqual(len(module.LIBRARIES), 7)

    def test_fetch_pins_match_the_deposited_bytes(self):
        for accession in ITEMS:
            record = summary(accession)
            script = (item_dir(accession) / 'script/fetch.sh').read_text()
            pins = {}
            for line in script.splitlines():
                match = re.search(r'/([A-Za-z0-9._-]+)\t\1\t(\d+)', line)
                if match:
                    pins[match.group(1)] = int(match.group(2))
            self.assertGreaterEqual(len(pins), 3, accession)
            deposition = record['deposition']
            entries = [value for key, value in deposition.items()
                       if isinstance(value, dict) and 'path' in value and 'bytes' in value]
            self.assertTrue(entries, accession)
            for entry in entries:
                name = Path(entry['path']).name
                self.assertIn(name, pins, f'{accession}: {name} not pinned')
                self.assertEqual(pins[name], entry['bytes'], f'{accession}: {name}')

    def test_author_log_pins_the_source_tree(self):
        expected_files = {'GSE138252': 4, 'GSE238242': 9, 'GSE269882': 3}
        for accession in ITEMS:
            log = read_json(item_dir(accession) / 'author_log.json')
            self.assertEqual(log['category'], 'disease_data')
            self.assertEqual(log['manifest'], 'manifest.sha256')
            self.assertEqual(log['verify_result'], 'pass')
            (source,) = log['sources']
            self.assertEqual(source['id'], accession)
            self.assertEqual(source['files'], expected_files[accession])
            self.assertEqual(len(source['tree_sha256']), 64)
            declared_paths = {entry['path'] for entry in log['local_additions']}
            self.assertIn('.gitignore', declared_paths)
            self.assertIn(f'local/deg/{accession}_deg.json', declared_paths)
            self.assertNotRegex(json.dumps(log), r'"/(home|mnt|tmp)/')


if __name__ == '__main__':
    unittest.main()
