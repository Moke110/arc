"""Disease-DEG contract for the isogenic CM models + genetic_cardiomyopathy (#188).

Hermetic wherever it can be: synthetic matrices/series matrices in a temp tree,
a tiny fixture repo for the driver replay, and the committed catalog/manifest/
DEG products read through the public reader. Every expected number is
hand-derived here; no test re-derives it by calling the implementation it is
checking.

Covers the ticket's acceptance items: the shared integerisation contract (rule,
deviation, library scale, genes dropped by rounding to zero, no NaN/negative
reaching voom), the as-deposited count check, GEO series-matrix sample binding,
the declared isogenic control arm of every source, the multi-arm strata/pooling
decision being recorded *and* applied, the published catalog schema, the
within-dataset rank semantics and multi-contrast collapse, determinism,
no absolute paths, and the genotype split (genotype-positive DCM / ARVC donors
belong only to genetic_cardiomyopathy, never to nHF).
"""
import importlib.util
import json
import shutil
import sys
import tempfile
import unittest
from importlib.machinery import SourceFileLoader
from pathlib import Path

import numpy as np
import pandas as pd

TESTS = Path(__file__).resolve().parent
EXECUTION = TESTS.parent
# ADR 0030: discover the root by the directory containing .git, never by counting
# parents; the sibling acceptance test uses the same walk.
REPO = next(path for path in Path(__file__).resolve().parents if (path / ".git").exists())
sys.path.insert(0, str(EXECUTION))
import public_api  # noqa: F401  (binds the single public implementation)
import disease_deg_fixtures as fixtures  # noqa: E402
import build_disease_deg as driver  # noqa: E402
import disease_deg_sources  # noqa: E402
from bases import bulk_deg, deg_rank, disease_deg  # noqa: E402

RESULT_ROOT = REPO / 'stage2/01_bases/output/disease-deg'
CATEGORY_DIR = RESULT_ROOT / 'genetic_cardiomyopathy'
ITEMS = ('GSE224129', 'GSE234141', 'GSE304575', 'GSE305648', 'GSE253223', 'GSE194297')
# The four deposits that split a case arm across strata/arms and must record how
# they are pooled: three collapse several contrasts, GSE253223 pools harvest days.
MULTI_ARM = ('gse224129__mybpc3_hcm', 'gse234141__ttntv_dcm', 'gse305648__rbm20_dcm',
             'gse253223__pkp2_knockdown_arvc')
MULTI_CONTRAST = ('gse224129__mybpc3_hcm', 'gse234141__ttntv_dcm', 'gse305648__rbm20_dcm')
SINGLE_ARM = ('gse304575__lmna_h222p_dcm', 'gse194297__dmd_cardiomyopathy')
EXPECTED_DATASETS = (
    'reichart-cm__dcm_genotype_positive', 'reichart-cm__arvc',
    'gse224129__mybpc3_hcm', 'gse234141__ttntv_dcm', 'gse304575__lmna_h222p_dcm',
    'gse305648__rbm20_dcm', 'gse253223__pkp2_knockdown_arvc',
    'gse194297__dmd_cardiomyopathy',
)


def load_item_spec(accession):
    """Import one item's ``script/build_deg.py`` and return its declared SPEC."""
    path = REPO / 'disease_data' / accession / 'script' / 'build_deg.py'
    loader = SourceFileLoader(f'isogenic_{accession}', str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module.SPEC


def write_series_matrix(path, samples):
    """Write a minimal GEO series matrix from ``[(gsm, title, description, chars)]``."""
    lines = []
    fields = [
        ('!Sample_title', lambda sample: sample['title']),
        ('!Sample_geo_accession', lambda sample: sample['gsm']),
        ('!Sample_description', lambda sample: sample['description']),
    ]
    for key, getter in fields:
        lines.append(key + '\t' + '\t'.join(f'"{getter(sample)}"' for sample in samples))
    keys = []
    for sample in samples:
        for characteristic in sample['characteristics']:
            keys.append(characteristic.split(':', 1)[0].strip())
    for index, label in enumerate(keys):
        values = []
        for sample in samples:
            matches = [item for item in sample['characteristics']
                       if item.split(':', 1)[0].strip() == label]
            values.append(matches[0])
        lines.append('!Sample_characteristics_ch1\t'
                     + '\t'.join(f'"{value}"' for value in values))
    lines.append('!series_matrix_table_begin')
    path.write_text('\n'.join(lines) + '\n')


def sample(gsm, title, description='x', characteristics=()):
    return {'gsm': gsm, 'title': title, 'description': description,
            'characteristics': list(characteristics)}


class IntegerisationContractTest(unittest.TestCase):
    """The one shared helper that puts a fractional count deposit on the voom path."""

    def frame(self, values):
        return pd.DataFrame(np.asarray(values, dtype=np.float64),
                            index=[f'G{index}' for index in range(len(values))],
                            columns=['S1', 'S2', 'S3'])

    def test_round_half_to_even_and_recorded_deviation(self):
        frame = pd.DataFrame(
            [[0.5, 1.5, 2.5], [100000.0, 100000.0, 100000.0]],
            index=['A', 'B'], columns=['S1', 'S2', 'S3'])
        integer, stats = bulk_deg.integerise_counts(frame)
        # round-half-to-even: 0.5 -> 0, 1.5 -> 2, 2.5 -> 2.
        self.assertEqual(integer.loc['A'].tolist(), [0, 2, 2])
        self.assertEqual(str(integer.dtypes.iloc[0]), 'int64')
        self.assertEqual(stats['rule'], 'round to nearest integer, ties to even (numpy.rint)')
        self.assertEqual(stats['n_non_integer_values'], 3)
        self.assertEqual(stats['max_abs_deviation'], 0.5)
        self.assertEqual(stats['n_genes_dropped_by_rounding'], 0)
        self.assertGreaterEqual(stats['library_size_before']['median'], 1e5)

    def test_genes_dropped_by_rounding_to_zero_are_reported(self):
        frame = pd.DataFrame(
            [[0.2, 0.1, 0.4], [100000.4, 100000.1, 100000.2], [0.0, 0.0, 0.0]],
            index=['tiny', 'solid', 'empty'], columns=['S1', 'S2', 'S3'])
        integer, stats = bulk_deg.integerise_counts(frame)
        self.assertEqual(integer.loc['tiny'].tolist(), [0, 0, 0])
        self.assertEqual(stats['n_genes_dropped_by_rounding'], 1)
        self.assertEqual(stats['dropped_gene_ids'], ['tiny'])
        self.assertEqual(stats['n_genes_all_zero_before'], 1)
        self.assertEqual(stats['n_genes_all_zero_after'], 2)

    def test_rejects_nan_and_negative_counts(self):
        frame = self.frame([[1.0, np.nan, 2.0]])
        with self.assertRaises(ValueError):
            bulk_deg.integerise_counts(frame)
        frame = self.frame([[-1.0, 2.0, 3.0]])
        with self.assertRaises(ValueError):
            bulk_deg.integerise_counts(frame)

    def test_rejects_a_normalised_deposit_declared_as_counts(self):
        # A log-scale/normalised matrix is not on a count scale: the shared guard
        # fails closed instead of rounding it into a meaningless count matrix.
        frame = self.frame([[1.2, 2.4, 0.7]])
        with self.assertRaises(ValueError):
            bulk_deg.integerise_counts(frame)
        with self.assertRaises(ValueError):
            bulk_deg.validate_integral_counts(self.frame([[1, 2, 3]]))


class CountDepositAsDepositedTest(unittest.TestCase):
    def test_accepts_integral_counts_and_records_scale(self):
        frame = pd.DataFrame([[10, 20, 30], [0, 5, 7]], index=['A', 'B'],
                             columns=['S1', 'S2', 'S3']).astype(float)
        frame = frame + 100000.0
        record = bulk_deg.validate_integral_counts(frame)
        self.assertEqual(record['mode'], 'integer-as-deposited')
        self.assertEqual(record['n_values'], 6)
        self.assertGreaterEqual(record['library_size']['min'], 1e5)

    def test_rejects_fractional_values_declared_integral(self):
        frame = pd.DataFrame([[10.5, 20, 30]], index=['A'], columns=['S1', 'S2', 'S3'])
        with self.assertRaises(ValueError):
            bulk_deg.validate_integral_counts(frame)


class SeriesMatrixBindingTest(unittest.TestCase):
    def test_parses_characteristics_and_matches_by_bracket_token(self):
        with tempfile.TemporaryDirectory() as scratch:
            path = Path(scratch) / 'series.txt'
            write_series_matrix(path, [
                sample('GSM1', 'clone A [s20020801E7]', characteristics=['genotype: +/+']),
                sample('GSM2', 'clone B [s20020802E7]', characteristics=['genotype: +/-']),
            ])
            series = bulk_deg.read_series_matrix(path)
            self.assertEqual(series.loc['GSM1', 'genotype'], '+/+')
            matched = bulk_deg.match_samples(['2-0020801E7', '2-0020802E7'], series,
                                             'title_bracket')
            self.assertEqual(matched.set_index('sample').loc['2-0020801E7', 'gsm'], 'GSM1')
            self.assertEqual(matched.set_index('sample').loc['2-0020802E7', 'gsm'], 'GSM2')

    def test_fails_closed_on_unmatched_and_duplicate_samples(self):
        with tempfile.TemporaryDirectory() as scratch:
            path = Path(scratch) / 'series.txt'
            write_series_matrix(path, [
                sample('GSM1', 'clone A [X1]', characteristics=['genotype: +/+']),
                sample('GSM2', 'clone B [X2]', characteristics=['genotype: +/-']),
            ])
            series = bulk_deg.read_series_matrix(path)
            with self.assertRaises(ValueError):
                bulk_deg.match_samples(['X1', 'X9'], series, 'title_bracket')
            with self.assertRaises(ValueError):
                bulk_deg.match_samples(['X1'], series, 'title_bracket')
            with self.assertRaises(ValueError):
                bulk_deg.match_samples(['X1', 'x1'], series, 'title_bracket')

    def test_conflicting_characteristic_fails_closed(self):
        with tempfile.TemporaryDirectory() as scratch:
            path = Path(scratch) / 'series.txt'
            write_series_matrix(path, [
                sample('GSM1', 'a [X1]', characteristics=['genotype: +/+']),
            ])
            text = path.read_text().replace(
                '!series_matrix_table_begin',
                '!Sample_characteristics_ch1\t"genotype: -/-"\n!series_matrix_table_begin')
            path.write_text(text)
            with self.assertRaises(ValueError):
                bulk_deg.read_series_matrix(path)

    def test_description_basis_requires_a_description(self):
        with tempfile.TemporaryDirectory() as scratch:
            path = Path(scratch) / 'series.txt'
            write_series_matrix(path, [sample('GSM1', 'a [X1]', description='lib1',
                                              characteristics=['genotype: +/+'])])
            series = bulk_deg.read_series_matrix(path)
            self.assertEqual(bulk_deg.match_samples(['lib1'], series, 'description')
                             .iloc[0]['gsm'], 'GSM1')


class ItemDeclarationTest(unittest.TestCase):
    """Every isogenic item declares its isogenic control and its stratum decision."""

    def test_all_six_items_declare_the_shared_derivation(self):
        for accession in ITEMS:
            spec = load_item_spec(accession)
            self.assertEqual(spec['kind'], 'counts', accession)
            self.assertIn('isogenic_control', spec)
            self.assertIn('strata', spec)
            self.assertTrue(spec['isogenic_control']['arm'], accession)
            self.assertTrue(spec['isogenic_control']['basis'], accession)

    def test_isogenic_control_is_the_declared_non_disease_arm(self):
        expected = {
            'GSE224129': ({'wt'}, 'cMyBP-C+/+'),
            'GSE234141': ({'2d_ventricular_wt', '2d_atrial_wt', 'eht_ventricular_wt',
                           'eht_atrial_wt'}, 'WT/WT'),
            'GSE304575': ({'corrected'}, 'Corrected Isogenic control'),
            'GSE305648': ({'wt'}, 'wildtype'),
            'GSE253223': ({'control'}, 'siNeg'),
            'GSE194297': ({'corrected'}, 'corrected'),
        }
        for accession, (arms, evidence) in expected.items():
            spec = load_item_spec(accession)
            controls = {contrast['denominator'] for contrast in spec['contrasts']}
            self.assertEqual(controls, arms, accession)
            deposit_evidence = spec['isogenic_control']['deposit_evidence']
            self.assertIn(evidence.split()[0], deposit_evidence, accession)

    def test_multi_arm_deposits_record_their_pooling_decision(self):
        for accession, dataset_id in zip(('GSE224129', 'GSE234141', 'GSE305648',
                                          'GSE253223'), MULTI_ARM):
            spec = load_item_spec(accession)
            strata = spec['strata']
            self.assertTrue(strata['pooled'], accession)
            self.assertTrue(strata['decision'], accession)
            self.assertTrue(strata['why'], accession)
            self.assertGreater(len(spec['contrasts']) + len(strata['blocking_covariates']), 1,
                               f'{accession}: a multi-arm deposit must pool arms or contrasts')
        self.assertIn('batch', load_item_spec('GSE305648')['strata']['blocking_covariates'])
        self.assertIn('day', load_item_spec('GSE253223')['strata']['blocking_covariates'])

    def test_single_arm_deposits_declare_no_pooling(self):
        for accession in ('GSE304575', 'GSE194297'):
            spec = load_item_spec(accession)
            self.assertFalse(spec['strata']['pooled'], accession)
            self.assertEqual(len(spec['contrasts']), 1, accession)

    def test_knockdown_caveat_is_recorded(self):
        spec = load_item_spec('GSE253223')
        self.assertIn('knockdown', json.dumps(spec['isogenic_control']).lower())
        self.assertEqual(spec['formula'], '~0 + group + day')


class ItemSummaryTest(unittest.TestCase):
    """The tracked item summaries: derivation record the published catalog relies on."""

    def summary(self, accession):
        path = REPO / 'disease_data' / accession / 'local/deg' / f'{accession}_deg.json'
        return json.loads(path.read_text())

    def test_rounding_record_is_complete_for_fractional_deposits(self):
        for accession in ('GSE224129', 'GSE253223'):
            values = self.summary(accession)['deposition']['values']
            self.assertIn('round to nearest', values['rule'])
            self.assertLess(values['max_abs_deviation'], 1.0)
            self.assertGreater(values['library_size_before']['median'], 1e6)
            self.assertGreaterEqual(values['n_genes_dropped_by_rounding'], 0)
            self.assertEqual(len(values['dropped_gene_ids']),
                             values['n_genes_dropped_by_rounding'])
        for accession in ('GSE234141', 'GSE304575', 'GSE305648', 'GSE194297'):
            values = self.summary(accession)['deposition']['values']
            self.assertEqual(values['mode'], 'integer-as-deposited', accession)

    def test_every_summary_is_self_consistent(self):
        for accession in ITEMS:
            record = self.summary(accession)
            self.assertEqual(record['schema'], 'stage2-item-deg-summary-v1')
            self.assertEqual(record['category'], 'genetic_cardiomyopathy')
            self.assertGreater(record['n_case'], 0)
            self.assertGreater(record['n_control'], 0)
            self.assertEqual(record['design']['formula'], load_item_spec(accession)['formula'])
            self.assertEqual(sorted(record['contrasts']),
                             sorted(c['label'] for c in load_item_spec(accession)['contrasts']))
            self.assertEqual(record['samples']['isogenic_control'],
                             load_item_spec(accession)['isogenic_control'])
            self.assertEqual(record['strata'], load_item_spec(accession)['strata'])
            self.assertEqual(len(record['deg']['sha256']), 64)
            self.assertEqual(len(record['method']['implementation_sha256']), 64)
            self.assertNotIn(str(REPO), json.dumps(record))


class PublishedCatalogTest(unittest.TestCase):
    def setUp(self):
        self.catalog = disease_deg.read_catalog(CATEGORY_DIR / 'catalog.json')
        self.manifest = json.loads((CATEGORY_DIR / 'manifest.json').read_text())
        self.records = {record['dataset_id']: record for record in self.catalog['datasets']}

    def test_eight_datasets_in_the_canonical_shape(self):
        self.assertEqual(self.catalog['schema'], disease_deg.CATALOG_SCHEMA)
        self.assertEqual(self.catalog['category'], 'genetic_cardiomyopathy')
        self.assertEqual(self.catalog['directory'], 'genetic_cardiomyopathy')
        self.assertEqual(tuple(self.records), tuple(EXPECTED_DATASETS))
        self.assertEqual(self.catalog['tier'], 1)
        for dataset_id, record in self.records.items():
            self.assertGreater(record['n_case'], 0, dataset_id)
            self.assertGreater(record['n_control'], 0, dataset_id)
            self.assertEqual(record['category'], 'genetic_cardiomyopathy')
            self.assertIn(record['evidence_type'], ('bulk', 'cm_pseudobulk'))
            self.assertIn(record['tier'], (1, 2))
            self.assertTrue(record['method'])
            self.assertTrue(record['inputs'])
            self.assertTrue(record['derived_inputs'])
            for derived in record['derived_inputs']:
                self.assertEqual(len(derived['sha256']), 64)
                self.assertFalse(derived['relative'].startswith('/'))
                self.assertTrue(derived['relative'].startswith('local/'))

    def test_isogenic_exception_recorded_per_source(self):
        self.assertTrue(self.catalog['isogenic_exception']['applies'])
        for dataset_id in EXPECTED_DATASETS:
            record = self.records[dataset_id]
            if dataset_id.startswith('reichart-cm__'):
                self.assertNotIn('isogenic_exception', record)
                continue
            exception = record['isogenic_exception']
            self.assertTrue(exception['applies'], dataset_id)
            self.assertIn('ADR 0066 Isogenic exception', exception['why'])
            self.assertTrue(exception['control_arm'], dataset_id)

    def test_multi_arm_strata_recorded_and_applied(self):
        for dataset_id in MULTI_ARM:
            record = self.records[dataset_id]
            self.assertTrue(record['strata']['pooled'], dataset_id)
            self.assertTrue(record['strata']['decision'], dataset_id)
        for dataset_id in MULTI_CONTRAST:
            self.assertGreater(len(self.records[dataset_id]['contrasts']), 1, dataset_id)
        # GSE253223 pools four harvest days with day as a blocking covariate
        # rather than producing four per-day contrasts.
        pkp2 = self.records['gse253223__pkp2_knockdown_arvc']
        self.assertEqual(pkp2['strata']['blocking_covariates'], ['day'])
        self.assertEqual(len(pkp2['contrasts']), 1)
        for dataset_id in SINGLE_ARM:
            self.assertFalse(self.records[dataset_id]['strata']['pooled'], dataset_id)
            self.assertEqual(len(self.records[dataset_id]['contrasts']), 1, dataset_id)

    def test_rounding_record_reaches_the_catalog_method(self):
        for dataset_id, genes in (('gse224129__mybpc3_hcm', '119'),
                                  ('gse253223__pkp2_knockdown_arvc', '4')):
            text = self.records[dataset_id]['method']['integerisation']
            self.assertIn('round to nearest, ties to even', text)
            self.assertIn('max |x - round(x)| 0.5', text)
            self.assertIn(f'{genes} genes whose rounded row is all-zero', text)
        for dataset_id in ('gse234141__ttntv_dcm', 'gse304575__lmna_h222p_dcm',
                           'gse305648__rbm20_dcm', 'gse194297__dmd_cardiomyopathy'):
            self.assertNotIn('integerisation', self.records[dataset_id]['method'])

    def test_reichart_datasets_record_the_read_only_derivation(self):
        for dataset_id in ('reichart-cm__dcm_genotype_positive', 'reichart-cm__arvc'):
            method = self.records[dataset_id]['method']
            self.assertIn('build_genetic_deg.py', method['derivation'])
            self.assertIn('not re-derived', method['derivation'])
            self.assertEqual(method['deposition'], 'integer counts (per-donor CM pseudobulk)')

    def test_manifest_replays_the_catalog_without_absolute_paths(self):
        self.assertEqual(self.manifest['schema'], disease_deg.RELEASE_SCHEMA)
        self.assertEqual(self.manifest['category'], 'genetic_cardiomyopathy')
        self.assertEqual(self.manifest['catalog_sha256'], self.catalog['catalog_sha256'])
        self.assertEqual(self.catalog['catalog_sha256'], disease_deg.catalog_digest(self.catalog))
        self.assertEqual(self.manifest['counts']['datasets'], 8)
        for name, digest in self.manifest['inputs'].items():
            self.assertFalse(name.startswith('/'))
            self.assertEqual(len(digest), 64)
        for name, value in self.manifest['outputs'].items():
            self.assertFalse(name.startswith('/'))
            self.assertEqual(len(value['sha256']), 64)
        self.assertNotIn(str(REPO), json.dumps(self.manifest))

    def test_reader_returns_tables_and_ranks(self):
        reader = disease_deg.DiseaseDeg(RESULT_ROOT)
        for dataset_id in EXPECTED_DATASETS:
            table = reader.table('genetic_cardiomyopathy', dataset_id)
            rank = reader.rank('genetic_cardiomyopathy', dataset_id)
            disease_deg.validate_table(table)
            disease_deg.validate_rank(rank)
            self.assertEqual(len(rank), deg_rank.AXIS_SIZE)
            self.assertEqual(rank.gene_index.tolist(), list(range(deg_rank.AXIS_SIZE)))


class RankSemanticsTest(unittest.TestCase):
    def setUp(self):
        self.catalog = disease_deg.read_catalog(CATEGORY_DIR / 'catalog.json')
        self.records = {record['dataset_id']: record for record in self.catalog['datasets']}

    def table(self, dataset_id):
        return pd.read_csv(CATEGORY_DIR / self.records[dataset_id]['deg'],
                           float_precision='round_trip')

    def test_effect_rank_and_deg_only_rank_norm(self):
        # ADR 0071 v2: only DEGs are ranked 1..n_DEG; every other gene sits at
        # n_DEG + 1 and rank_norm = rank / n_DEG (no n_measured denominator).
        table = pd.read_csv(CATEGORY_DIR / self.records['gse224129__mybpc3_hcm']['deg'],
                            float_precision='round_trip')
        padj = table.padj.to_numpy(dtype=np.float64)
        log2fc = table.log2fc.to_numpy(dtype=np.float64)
        expected = np.where(np.isnan(padj), 0.0, np.abs(log2fc) * -np.log10(padj))
        np.testing.assert_allclose(table.effect.to_numpy(dtype=np.float64), expected,
                                   rtol=0, atol=1e-12)
        deg = table.deg.to_numpy(dtype=bool)
        n_deg = int(deg.sum())
        self.assertGreater(n_deg, 0)
        self.assertEqual(int((table.n_contrasts_deg > 0).sum()), n_deg)
        ranks = table['rank'].to_numpy(dtype=np.int64)
        self.assertEqual(sorted(ranks.tolist()),
                         list(range(1, n_deg + 1)) + [n_deg + 1] * (len(table) - n_deg))
        np.testing.assert_allclose(table.rank_norm.to_numpy(dtype=np.float64),
                                   ranks / float(n_deg))
        # the DEGs are ordered by descending effect with a gene_index tie-break
        order = np.lexsort((table.gene_index.to_numpy(dtype=np.int64)[deg],
                            -table.effect.to_numpy(dtype=np.float64)[deg]))
        np.testing.assert_array_equal(ranks[deg], np.arange(1, n_deg + 1))
        np.testing.assert_array_equal(
            table.gene_index.to_numpy(dtype=np.int64)[deg],
            table.gene_index.to_numpy(dtype=np.int64)[deg][order])

    def test_non_deg_genes_get_n_deg_plus_one(self):
        record = self.records['gse234141__ttntv_dcm']
        rank = pd.read_csv(CATEGORY_DIR / record['rank'])
        n_deg = int(record['n_deg'])
        degs = rank[rank['rank'] <= n_deg]
        non_degs = rank[rank['rank'] > n_deg]
        self.assertEqual(len(degs), n_deg)
        self.assertEqual(set(non_degs['rank'].unique()), {n_deg + 1})
        self.assertEqual(len(degs) + len(non_degs), deg_rank.AXIS_SIZE)
        np.testing.assert_allclose(non_degs.rank_norm.to_numpy(dtype=np.float64),
                                   (n_deg + 1) / float(n_deg))

    def test_multi_contrast_collapse_keeps_the_maximum_effect_contrast(self):
        for dataset_id in MULTI_CONTRAST:
            record = self.records[dataset_id]
            table = self.table(dataset_id)
            self.assertGreater(len(record['contrasts']), 1)
            # every declared contrast won at least one gene (else declaration and
            # fit disagree) and the published labels are a subset of the declared.
            published = set(table.contrast.astype(str))
            self.assertTrue(published <= set(record['contrasts']), dataset_id)
            self.assertEqual(published, set(record['contrasts']), dataset_id)
            self.assertTrue(table.significant.any(), dataset_id)


class GenotypeSplitTest(unittest.TestCase):
    def setUp(self):
        self.genotypes = pd.read_csv(
            REPO / 'disease_data/Reichart-CM/local/genotype/reichart-cm_donor_genotypes.csv')

    def test_positive_and_arvc_donors_are_export_only(self):
        split, summary = disease_deg.genotype_split(self.genotypes[['donor_id', 'genotype']])
        self.assertEqual(summary.get('negative:case'), 8)
        self.assertEqual(summary.get('negative:control'), 18)
        self.assertEqual(summary.get('pathogenic:export_only'), 52)
        self.assertNotIn('pathogenic:excluded', summary)  # the NCCM donor is not exported
        self.assertEqual(len(self.genotypes), 78)

    def test_genotype_positive_and_arvc_only_in_this_category(self):
        nhf = disease_deg_sources.datasets('nHF')
        labels = [contrast for record in nhf for contrast in record['contrasts']]
        self.assertIn('dcm_genotype_negative-normal', labels)
        for label in labels:
            self.assertNotIn('positive', label)
            self.assertNotIn('arvc', label)
        genetic = {record['dataset_id']: record
                   for record in disease_deg_sources.datasets('genetic_cardiomyopathy')}
        self.assertEqual(genetic['reichart-cm__dcm_genotype_positive']['n_case'], 44)
        self.assertEqual(genetic['reichart-cm__dcm_genotype_positive']['n_control'], 18)
        self.assertEqual(genetic['reichart-cm__arvc']['n_case'], 8)
        self.assertEqual(genetic['reichart-cm__arvc']['n_control'], 18)
        self.assertEqual(genetic['reichart-cm__dcm_genotype_positive']['shared_control_group'],
                         'reichart-cm:normal')

    def test_exported_pseudobulk_has_no_nccm_donor(self):
        pseudobulk = pd.read_csv(
            REPO / 'disease_data/Reichart-CM/local/pseudobulk/reichart-cm_cm_pseudobulk.csv.gz',
            nrows=2)
        donors = [column for column in pseudobulk.columns
                  if column not in ('gene_index', 'gene_symbol')]
        self.assertEqual(len(donors), 78)
        kept = set(self.genotypes.donor_id)
        self.assertEqual(set(donors), kept)
        # The single NCCM donor (H84, TPM1) is excluded from the export entirely,
        # so it is neither a column nor an exported role.
        self.assertNotIn('H84', set(self.genotypes.donor_id))
        self.assertFalse((self.genotypes.nHF_role == 'excluded').any())


class DerivationPathTest(unittest.TestCase):
    """A synthetic fractional deposit goes through the shared helper into voom."""

    def test_fractional_counts_become_finite_integers_before_voom(self):
        with tempfile.TemporaryDirectory() as scratch:
            item = Path(scratch) / 'GSEFIXTURE'
            (item / 'source').mkdir(parents=True)
            genes = ['G00000', 'G00001', 'G00002', 'G00003']
            matrix = pd.DataFrame({
                '1-AAA': [1000000.4, 50000.6, 40.4, 30000.0],
                '2-BBB': [1200000.2, 40000.3, 60.6, 25000.0],
                '3-CCC': [1100000.6, 45000.1, 20.2, 28000.0],
                '4-DDD': [1050000.1, 48000.7, 30.3, 26000.0],
            }, index=genes)
            matrix.to_csv(item / 'source/matrix.csv.gz', sep=';',
                          compression='gzip', index_label='gene')
            write_series_matrix(item / 'source/series.txt', [
                sample('GSM1', 'wt 1 [1-AAA]', characteristics=['genotype: WT/WT']),
                sample('GSM2', 'wt 2 [2-BBB]', characteristics=['genotype: WT/WT']),
                sample('GSM3', 'het 1 [3-CCC]', characteristics=['genotype: WT/KO']),
                sample('GSM4', 'het 2 [4-DDD]', characteristics=['genotype: WT/KO']),
            ])
            spec = {
                'matrix': 'source/matrix.csv.gz', 'sep': ';',
                'series_matrix': 'source/series.txt', 'id_kind': 'symbol',
                'sample_match': 'title_bracket', 'group_source': 'genotype',
                'group_map': {'WT/WT': 'wt', 'WT/KO': 'het'},
                'kind': 'counts', 'integerisation': 'round_half_to_even',
                'formula': '~0 + group',
                'contrasts': [{'label': 'het-wt', 'numerator': 'het',
                               'denominator': 'wt'}],
                'contrast_direction': 'positive log2FC = up in het vs wt',
                'isogenic_control': {'arm': 'WT/WT', 'basis': 'synthetic',
                                     'deposit_evidence': 'synthetic'},
                'strata': {'decision': 'single pair', 'pooled': False,
                           'blocking_covariates': []},
            }
            result = bulk_deg.derive_item_deg(item, spec)
            self.assertEqual(result['result'], 'PASS')
            summary = json.loads((item / 'local/deg/GSEFIXTURE_deg.json').read_text())
            values = summary['deposition']['values']
            self.assertEqual(values['rule'], 'round to nearest integer, ties to even (numpy.rint)')
            self.assertEqual(summary['contrasts']['het-wt']['n_case'], 2)
            self.assertEqual(summary['contrasts']['het-wt']['n_control'], 2)
            self.assertIn('limma_version', summary['method'])
            expression = pd.read_csv(item / 'local/deg/GSEFIXTURE_expr.tsv.gz', sep='\t')
            numeric = expression.drop(columns=['source_gene_id', 'gene_symbol']).to_numpy()
            self.assertTrue(np.isfinite(numeric).all())
            self.assertTrue((numeric >= 0).all())
            np.testing.assert_array_equal(numeric, np.rint(numeric))
            # the deposited fractional extremes are recorded, not silently absorbed
            self.assertEqual(values['n_non_integer_values'], 12)


def _fixture_sources():
    categories = {
        'genetic_cardiomyopathy': {
            'directory': 'genetic_cardiomyopathy', 'tier': 1,
            'evidence_type': 'cm_pseudobulk',
            'shared_control_group': ['reichart-cm:normal'],
            'isogenic_exception': {'applies': True, 'statement': 'fixture exception'},
            'rejected': [],
        },
    }
    base = {
        'category': 'genetic_cardiomyopathy',
        'species': 'Homo sapiens', 'model_system': 'fixture', 'evidence_type': 'bulk',
        'tier': 2, 'genotype_status': 'pathogenic', 'time_window_basis': None,
        'exploratory': False, 'shared_control_group': 'fixture:isogenic',
        'method': {'deposition': 'integer counts', 'model': 'limma-voom',
                   'design': 'one joint fit', 'test': 'moderated t', 'padj': 'BH'},
    }
    datasets = [
        dict(base, dataset_id='iso__multi', source_item='fixture_items/ISO-FIX',
             local_deg='local/deg/iso_deg.csv', local_summary='local/deg/iso_deg.json',
             arm_or_class='two contrasts', n_case=4, n_control=2,
             contrasts=['A-wt', 'B-wt'],
             isogenic_exception={'applies': True, 'control_arm': 'wt', 'why': 'fixture'},
             strata={'decision': 'collapsed', 'pooled': True, 'blocking_covariates': []}),
        dict(base, dataset_id='reichart__positive', source_item='fixture_items/REICH-FIX',
             local_deg='local/deg/reichart_deg.csv',
             local_summary='local/deg/reichart_deg.json', arm_or_class='positive',
             n_case=3, n_control=1, contrasts=['positive-normal'],
             evidence_type='cm_pseudobulk', tier=1),
    ]
    return categories, datasets


def _write_driver_fixture(root):
    root = Path(root)
    identity = root / 'stage2/01_bases/output/gene-identity'
    (identity / 'tables').mkdir(parents=True, exist_ok=True)
    fixtures._axis().to_csv(identity / 'genes.csv', index=False)
    fixtures._crosswalk().to_csv(identity / 'tables/gene_crosswalk.csv', index=False)
    rows = []
    for gene, effect_rows in (
        ('G00000', [('A-wt', 2.0, 1e-4), ('B-wt', 1.0, 0.5)]),
        ('G00001', [('A-wt', 1.0, 0.5), ('B-wt', 3.0, 1e-3)]),
        ('G00002', [('A-wt', -1.5, None)]),
    ):
        for contrast, log2fc, padj in effect_rows:
            rows.append({'source_gene_id': gene, 'gene_symbol': gene,
                         'contrast': contrast, 'log2fc': log2fc, 'stat': log2fc,
                         'p_value': 1e-3 if padj is not None else np.nan,
                         'padj': padj, 'base_mean': 1.0})
    iso = root / 'fixture_items/ISO-FIX'
    (iso / 'local/deg').mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(iso / 'local/deg/iso_deg.csv', index=False)
    (iso / 'local/deg/iso_deg.json').write_text('{fixture}\n')
    (iso / 'manifest.sha256').write_text('00  fixture\n')
    reich = root / 'fixture_items/REICH-FIX'
    (reich / 'local/deg').mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{'source_gene_id': 'G00003', 'gene_symbol': 'G00003',
                   'contrast': 'positive-normal', 'log2fc': 1.0, 'stat': 2.0,
                   'p_value': 0.01, 'padj': 0.01, 'base_mean': 1.0}]
                 ).to_csv(reich / 'local/deg/reichart_deg.csv', index=False)
    (reich / 'local/deg/reichart_deg.json').write_text('{fixture}\n')
    (reich / 'manifest.sha256').write_text('00  fixture\n')
    categories, datasets = _fixture_sources()
    module = root / 'fixture_sources.py'
    module.write_text(
        '"""Fixture source declaration for test_disease_deg_isogenic.py."""\n'
        f'CATEGORIES = {categories!r}\n'
        f'DATASETS = {datasets!r}\n'
        'def datasets(category):\n'
        '    return [record for record in DATASETS if record["category"] == category]\n'
        'def category_record(category):\n'
        '    return CATEGORIES[category]\n')
    return root, module


class DriverReplayTest(unittest.TestCase):
    def setUp(self):
        self.scratch = Path(tempfile.mkdtemp(prefix='ticket188-'))
        self.root, self.sources = _write_driver_fixture(self.scratch / 'repo')
        self.output = self.scratch / 'repo/stage2/01_bases/output/disease-deg'

    def tearDown(self):
        shutil.rmtree(self.scratch, ignore_errors=True)

    def publish(self, root):
        from bases import release as release_module
        context = driver.Context(root, self.root, driver.load_sources(self.sources))
        output = root / 'genetic_cardiomyopathy'
        release_module.publish(
            output,
            lambda staging: driver.build_category(
                context, staging, 'genetic_cardiomyopathy'),
            driver.verify_category, root=self.scratch / 'checkpoints')
        return driver.verify_category(output)

    def test_publish_twice_is_byte_identical(self):
        first = self.publish(self.output)
        first_bytes = {name: (self.output / 'genetic_cardiomyopathy' / name).read_bytes()
                       for name in ('catalog.json', 'manifest.json')}
        self.output.joinpath('genetic_cardiomyopathy').rename(
            self.scratch / 'first-release')
        second = self.publish(self.output)
        self.assertEqual(first['manifest_sha256'], second['manifest_sha256'])
        for name, payload in first_bytes.items():
            self.assertEqual(payload,
                             (self.output / 'genetic_cardiomyopathy' / name).read_bytes())

    def test_fixture_catalog_collapses_contrasts_by_maximum_effect(self):
        context = driver.Context(self.output, self.root, driver.load_sources(self.sources))
        (self.scratch / 'candidate').mkdir()
        result = driver.build_category(context, self.scratch / 'candidate',
                                       'genetic_cardiomyopathy')
        self.assertEqual(result['counts']['datasets'], 2)
        candidate = self.scratch / 'candidate'
        catalog = json.loads((candidate / 'catalog.json').read_text())
        record = next(record for record in catalog['datasets']
                      if record['dataset_id'] == 'iso__multi')
        self.assertEqual(record['n_measured'], 3)
        self.assertEqual(record['n_axis_unmeasured'], fixtures.AXIS_SIZE - 3)
        table = pd.read_csv(candidate / record['deg'], float_precision='round_trip')
        # effects: G00001 via B = 3 * 3 = 9 (rank 1), G00000 via A = 2 * 4 = 8
        # (rank 2), G00002 via A = 0 (rank 3).
        self.assertEqual(table.gene_index.tolist(), [1, 0, 2])
        self.assertEqual(table.contrast.tolist(), ['B-wt', 'A-wt', 'A-wt'])
        np.testing.assert_allclose(table.effect.to_numpy(), [9.0, 8.0, 0.0])
        self.assertEqual(record['strata']['pooled'], True)
        self.assertEqual(record['isogenic_exception']['control_arm'], 'wt')
        self.assertTrue(catalog['isogenic_exception']['applies'])
        rank = pd.read_csv(candidate / record['rank'])
        # G00001 and G00000 are the DEGs (ranks 1 and 2); every other gene sits
        # at n_deg + 1 = 3 under the v2 DEG-only rank.
        self.assertEqual(record['n_deg'], 2)
        self.assertEqual(int(rank.loc[3, 'rank']), 3)
        np.testing.assert_allclose(float(rank.loc[3, 'rank_norm']), 1.5)

    def test_undeclared_and_never_winning_contrasts_fail_verification(self):
        from bases import release as release_module
        context = driver.Context(self.output, self.root, driver.load_sources(self.sources))
        (self.scratch / 'candidate').mkdir()
        driver.build_category(context, self.scratch / 'candidate', 'genetic_cardiomyopathy')
        candidate = self.scratch / 'candidate'
        catalog = json.loads((candidate / 'catalog.json').read_text())
        record = next(record for record in catalog['datasets']
                      if record['dataset_id'] == 'iso__multi')
        # The fit produced A-wt and B-wt; dropping B-wt from the declaration must
        # fail verification ("undeclared contrast in the published table").
        record['contrasts'] = ['A-wt']
        disease_deg.write_catalog(candidate / 'catalog.json', catalog)
        manifest = json.loads((candidate / 'manifest.json').read_text())
        manifest['catalog_sha256'] = catalog['catalog_sha256']
        manifest['outputs']['catalog.json'] = {
            'sha256': release_module.sha256(candidate / 'catalog.json')}
        (candidate / 'manifest.json').write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + '\n')
        with self.assertRaises(release_module.ReleaseError):
            driver.verify_category(candidate)

    def test_declared_but_never_winning_contrast_fails_verification(self):
        from bases import release as release_module
        context = driver.Context(self.output, self.root, driver.load_sources(self.sources))
        (self.scratch / 'candidate').mkdir()
        driver.build_category(context, self.scratch / 'candidate', 'genetic_cardiomyopathy')
        candidate = self.scratch / 'candidate'
        catalog = json.loads((candidate / 'catalog.json').read_text())
        record = next(record for record in catalog['datasets']
                      if record['dataset_id'] == 'iso__multi')
        # A declared arm that wins no gene means declaration and fit disagree.
        record['contrasts'] = ['A-wt', 'B-wt', 'C-wt']
        disease_deg.write_catalog(candidate / 'catalog.json', catalog)
        manifest = json.loads((candidate / 'manifest.json').read_text())
        manifest['catalog_sha256'] = catalog['catalog_sha256']
        manifest['outputs']['catalog.json'] = {
            'sha256': release_module.sha256(candidate / 'catalog.json')}
        (candidate / 'manifest.json').write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + '\n')
        with self.assertRaises(release_module.ReleaseError):
            driver.verify_category(candidate)

if __name__ == '__main__':
    unittest.main()
