"""Disease-DEG publication contract for the human CM sources (ticket #187).

Hermetic: the fixture writes a tiny repository-shaped tree (frozen 18,164-gene
axis, crosswalk with alias/ambiguous rows, fake item manifests and hand-worked
item-local DEG tables), then the real driver publishes/replays against it. Every
expected number is hand-derived in the test, not read back from the build.

Covers the ticket's acceptance items: DEG table schema, rank semantics and
gene_index tie-breaks, NA/unmeasured handling, genotype split, one-joint-design
arm collapse, determinism (publish twice -> identical digests) and the
no-absolute-paths rule for manifests.
"""
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

TESTS = Path(__file__).resolve().parent
EXECUTION = TESTS.parent
sys.path.insert(0, str(EXECUTION))
import public_api  # noqa: F401  (binds the single public implementation)
import build_disease_deg as driver  # noqa: E402
import disease_deg_fixtures as fixtures  # noqa: E402
from bases import disease_deg, deg_rank  # noqa: E402

SCRATCH = EXECUTION / 'checkpoints' / 'ticket187-tests'


def resolver_fixture():
    axis = fixtures._axis()
    crosswalk = fixtures._crosswalk()
    return disease_deg.build_resolver(axis, crosswalk), axis.symbol.to_numpy(dtype=object)


def products(local, dataset_id):
    resolver, symbols = resolver_fixture()
    return disease_deg.build_dataset_products(local, resolver, symbols, dataset_id)


class DegTableContractTest(unittest.TestCase):
    def setUp(self):
        self.local = disease_deg.contrast_selection(
            fixtures.gse270788_deg(), 'gse270788__nicm', ['NICM-Donor'])
        self.table, self.rank, self.stats, self.report = products(
            self.local, 'gse270788__nicm')

    def test_published_table_schema(self):
        self.assertEqual(tuple(self.table.columns), disease_deg.DEG_COLUMNS)
        self.assertEqual(len(self.table), 4)
        self.assertEqual(self.table.gene_index.tolist(), [1, 5, 2, 3])
        self.assertEqual(self.table['rank'].tolist(), [1, 2, 3, 4])
        np.testing.assert_allclose(self.table.rank_norm.to_numpy(), [0.25, 0.5, 0.75, 1.0])
        self.assertEqual(self.table.gene_symbol.tolist(),
                         ['G00001', 'G00005', 'G00002', 'G00003'])
        # |log2FC| * -log10(padj): gene 0 = 2 * 4 = 8; gene 4 = 1.5 * 3 = 4.5.
        np.testing.assert_allclose(self.table.effect.to_numpy(), [8.0, 4.5, 2.0, 0.0])

    def test_na_padj_is_measured_but_non_significant(self):
        row = self.table[self.table.gene_index == 3].iloc[0]
        self.assertFalse(bool(row.significant))
        self.assertEqual(int(row.direction), 0)
        self.assertEqual(float(row.effect), 0.0)
        self.assertTrue(np.isnan(row.padj))

    def test_unmeasured_gene_gets_n_measured_plus_one(self):
        self.assertEqual(self.stats['n_measured'], 4)
        self.assertEqual(self.stats['unmeasured_rank'], 5)
        self.assertEqual(int(self.rank.loc[0, 'rank']), 5)
        self.assertEqual(float(self.rank.loc[0, 'rank_norm']), 1.0)
        self.assertEqual(int(self.rank.loc[1, 'rank']), 1)
        self.assertEqual(len(self.rank), fixtures.AXIS_SIZE)
        self.assertEqual(self.rank.gene_index.tolist(), list(range(fixtures.AXIS_SIZE)))

    def test_alignment_reports_every_loss(self):
        report = self.report
        self.assertEqual(report['n_source_rows'], 7)
        self.assertEqual(report['n_resolved'], 5)
        self.assertEqual(report['alignment_categories']['ambiguous'], 1)
        self.assertEqual(report['alignment_categories']['unmapped'], 1)
        self.assertEqual(report['alignment_categories']['crosswalk_alias'], 2)
        self.assertEqual(report['n_axis_genes_resolved'], 4)

    def test_duplicate_alias_takes_max_effect_row(self):
        # PREV_G2 aliases to gene_index 1 but its effect is smaller than G00002's,
        # so the table reports the G00002 statistic, not the alias row's.
        row = self.table[self.table.gene_index == 2].iloc[0]
        self.assertAlmostEqual(float(row.log2fc), -1.0)
        self.assertAlmostEqual(float(row.stat), -3.0)

    def test_validate_rejects_broken_tables(self):
        broken = self.table.drop(columns=['base_mean'])
        with self.assertRaises(ValueError):
            disease_deg.validate_table(broken)
        unsorted = self.table.iloc[::-1].reset_index(drop=True)
        with self.assertRaises(ValueError):
            disease_deg.validate_table(unsorted)
        with self.assertRaises(ValueError):
            disease_deg.validate_rank(self.rank.iloc[:100])


class TieBreakTest(unittest.TestCase):
    def test_equal_effect_breaks_on_gene_index(self):
        table, rank, stats, _ = products(fixtures.gse255612_deg(), 'gse255612__af')
        first, second = table.iloc[0], table.iloc[1]
        self.assertEqual(int(first.gene_index), 10)
        self.assertEqual(int(second.gene_index), 11)
        self.assertAlmostEqual(float(first.effect), float(second.effect))
        self.assertLess(int(first['rank']), int(second['rank']))
        self.assertEqual(int(rank.loc[12, 'rank']), 3)
        self.assertTrue(np.isnan(table.iloc[2].padj))
        self.assertFalse(bool(table.iloc[2].significant))


class JointDesignTest(unittest.TestCase):
    def test_one_fit_supplies_every_arm_contrast(self):
        joint = fixtures.gse270788_deg()
        self.assertEqual(sorted(joint.contrast.unique()),
                         ['AMI-Donor', 'ICM-Donor', 'NICM-Donor'])
        icm = disease_deg.contrast_selection(joint, 'gse270788__icm', ['ICM-Donor'])
        ami = disease_deg.contrast_selection(joint, 'gse270788__ami', ['AMI-Donor'])
        self.assertEqual(set(icm.contrast), {'ICM-Donor'})
        self.assertEqual(set(ami.contrast), {'AMI-Donor'})
        self.assertEqual(len(icm), len(ami))

    def test_missing_contrast_and_empty_selection_are_errors(self):
        joint = fixtures.gse270788_deg()
        with self.assertRaises(ValueError):
            disease_deg.contrast_selection(joint, 'gse270788__icm', ['AMI-Donor-typo'])
        selected = disease_deg.contrast_selection(
            joint, 'gse270788__icm', ['ICM-Donor', 'AMI-Donor'])
        self.assertEqual(set(selected.contrast), {'ICM-Donor', 'AMI-Donor'})


class GenotypeSplitTest(unittest.TestCase):
    def setUp(self):
        self.frame = pd.DataFrame({
            'donor_id': ['D1', 'D2', 'D3', 'D4', 'D5', 'D6'],
            'genotype': ['normal', 'dcm_genotype_negative', 'dcm_genotype_positive',
                         'arvc', 'nccm', 'normal'],
        })
        self.split, self.summary = disease_deg.genotype_split(self.frame)

    def test_positive_genotype_leaves_nhf(self):
        roles = dict(zip(self.split.donor_id, self.split.nHF_role))
        self.assertEqual(roles['D1'], 'control')
        self.assertEqual(roles['D2'], 'case')
        self.assertEqual(roles['D3'], 'export_only')
        self.assertEqual(roles['D4'], 'export_only')
        self.assertEqual(roles['D5'], 'excluded')
        nHF_cases = set(self.split.loc[self.split.nHF_role == 'case', 'donor_id'])
        nHF_controls = set(self.split.loc[self.split.nHF_role == 'control', 'donor_id'])
        self.assertNotIn('D3', nHF_cases)       # genotype-positive DCM is not nHF
        self.assertEqual(nHF_controls, {'D1', 'D6'})
        self.assertEqual(self.split.loc[self.split.donor_id == 'D5'].genotype_status.iloc[0],
                         'pathogenic')

    def test_unknown_label_fails(self):
        with self.assertRaises(ValueError):
            disease_deg.genotype_split(pd.DataFrame({'donor_id': ['X'],
                                                     'genotype': ['sarcoidosis']}))
        with self.assertRaises(ValueError):
            disease_deg.genotype_split(pd.DataFrame({'donor_id': ['X', 'X'],
                                                     'genotype': ['normal', 'normal']}))


class ReleaseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        SCRATCH.mkdir(parents=True, exist_ok=True)
        cls.root = Path(tempfile.mkdtemp(prefix='fixture-', dir=SCRATCH))
        cls.input_root, cls.sources = fixtures.write_fixture(cls.root)
        cls.output = cls.root / 'out'
        cls.command = ['publish', '--all', '--root', str(cls.output),
                       '--input-root', str(cls.input_root),
                       '--sources', str(cls.sources)]
        assert driver.main(cls.command) == 0

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.root, ignore_errors=True)

    def _manifests(self):
        return {category: json.loads((self.output / category / 'manifest.json').read_text())
                for category in ('af', 'icm', 'nhf')}

    def test_published_layout(self):
        for category in ('af', 'icm', 'nhf'):
            self.assertTrue((self.output / category / 'catalog.json').is_file())
            self.assertTrue((self.output / category / 'manifest.json').is_file())
        catalog = json.loads((self.output / 'nhf' / 'catalog.json').read_text())
        self.assertEqual(catalog['schema'], disease_deg.CATALOG_SCHEMA)
        self.assertEqual([record['dataset_id'] for record in catalog['datasets']],
                         ['gse270788__nicm', 'reichart-cm__dcm_genotype_negative'])
        self.assertEqual(catalog['axis_size'], fixtures.AXIS_SIZE)
        record = next(record for record in catalog['datasets']
                      if record['dataset_id'] == 'gse270788__nicm')
        for key in ('species', 'model_system', 'evidence_type', 'tier',
                    'genotype_status', 'time_window_basis', 'exploratory',
                    'shared_control_group', 'n_case', 'n_control', 'method',
                    'deg', 'rank'):
            self.assertIn(key, record)

    def test_every_dataset_carries_case_control_counts_and_pinned_derived_inputs(self):
        for category in ('af', 'icm', 'nhf'):
            catalog = json.loads((self.output / category / 'catalog.json').read_text())
            self.assertTrue(catalog['datasets'])
            for record in catalog['datasets']:
                # A bulk/pseudobulk DEG always has both arms; zero is a metadata bug.
                self.assertGreater(record['n_case'], 0, record['dataset_id'])
                self.assertGreater(record['n_control'], 0, record['dataset_id'])
                self.assertTrue(record['derived_inputs'], record['dataset_id'])
                for entry in record['derived_inputs']:
                    self.assertEqual(len(entry['sha256']), 64)
                    self.assertFalse(entry['relative'].startswith('disease_data/'))
                for key in record['inputs']:
                    self.assertFalse(key.startswith('/'))
                    self.assertTrue((self.input_root / key).is_file(), key)

    def test_manifests_have_no_absolute_paths(self):
        for category, manifest in self._manifests().items():
            driver.release.check_relative_paths(manifest)
            for key in manifest['inputs']:
                self.assertFalse(key.startswith('/'), f'{category}: absolute input {key}')
            for key in manifest['outputs']:
                self.assertFalse(key.startswith('/'), f'{category}: absolute output {key}')

    def test_publish_twice_is_byte_identical(self):
        first = {category: (self.output / category / 'manifest.json').read_bytes()
                 for category in ('af', 'icm', 'nhf')}
        assert driver.main(self.command) == 0
        for category in ('af', 'icm', 'nhf'):
            self.assertEqual((self.output / category / 'manifest.json').read_bytes(),
                             first[category], category)

    def test_verify_replays_the_release(self):
        assert driver.main(['verify', '--all', '--root', str(self.output),
                            '--input-root', str(self.input_root),
                            '--sources', str(self.sources)]) == 0

    def test_reader_returns_tables_ranks_and_eligibility(self):
        reader = disease_deg.DiseaseDeg(self.output)
        self.assertEqual(reader.categories(), ['af', 'icm', 'nhf'])
        self.assertEqual(reader.datasets('ICM'), ['gse270788__icm', 'gse270788__ami'])
        table = reader.table('ICM', 'gse270788__ami')
        self.assertEqual(tuple(table.columns), disease_deg.DEG_COLUMNS)
        self.assertEqual(table.contrast.unique().tolist(), ['AMI-Donor'])
        rank = reader.rank('nHF', 'gse270788__nicm')
        self.assertEqual(len(rank), fixtures.AXIS_SIZE)
        eligibility = reader.eligibility('ICM', 'gse270788__ami')
        self.assertEqual(eligibility['time_window_basis'], 'author_declared')
        self.assertEqual(eligibility['evidence_type'], 'cm_pseudobulk')
        self.assertEqual(eligibility['tier'], 1)
        self.assertEqual(eligibility['n_case'], 5)

    def test_tampered_payload_is_rejected(self):
        reader = disease_deg.DiseaseDeg(self.output)
        path = self.output / 'af' / 'gse255612__af.deg.csv'
        original = path.read_bytes()
        path.write_bytes(original.replace(b'G00010', b'G00099', 1))
        try:
            with self.assertRaises(ValueError):
                reader.table('AF', 'gse255612__af')
        finally:
            path.write_bytes(original)


if __name__ == '__main__':
    unittest.main()
