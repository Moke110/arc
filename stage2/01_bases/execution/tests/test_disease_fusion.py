#!/usr/bin/env python
"""Per-disease level-1 fusion, ``DEG2000`` and the catalog shape (#205).

Covers the ADR 0071 decisions 3 and 4 from the outside:

* ``source_deg_rank`` turns one published compact DEG table into its v2 DEG-only
  rank, from both published table generations (the v2 pass-any ``deg`` flag and
  the v1 collapsed winning-contrast rows);
* ``level_one_rank``/``fuse_disease`` sum ``rank / n_DEG`` over the informative
  sources, rank with the ``gene_index`` tie-break, exclude zero-DEG sources with
  a recorded reason (never a NaN), and cut ``DEG2000`` at exactly the declared
  limit — all of it when the union is smaller;
* the AF extreme case (one informative source, ``n_DEG = 1``) and the empty
  union are representable, including through the published driver;
* the category catalog shape (``tier``/``evidence_type``/``n_deg``/``deg2000``
  and a consistent ``tiering`` record), the emitted ``output/disease-fusion/``
  product and its replay are validated, path-free and deterministic.

The fixtures are a synthetic 18,164-gene frozen axis and tiny per-source compact
tables; the only real-file test replays the tracked products, which carry no
external payload.
"""
import hashlib
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
# parents; the sibling product tests use the same walk.
REPO = next(path for path in Path(__file__).resolve().parents if (path / ".git").exists())
sys.path.insert(0, str(EXECUTION))
import public_api  # noqa: F401  (binds the single public implementation)
import disease_deg_sources  # noqa: E402
import build_disease_fusion as driver  # noqa: E402
from bases import disease_deg, gene_axis  # noqa: E402

AXIS_SIZE = 18164
RESULT_ROOT = REPO / 'stage2/01_bases/output/disease-fusion'


# ---------------------------------------------------------------- fixtures

def axis_frame():
    """The frozen 18,164-gene identity table (symbols only; no payload read)."""
    return pd.DataFrame({
        'gene_index': np.arange(AXIS_SIZE, dtype=np.int64),
        'symbol': [f'G{index:05d}' for index in range(AXIS_SIZE)],
        'hgnc_id': [f'HGNC:{index + 1}' for index in range(AXIS_SIZE)],
        'ensembl_id': [f'ENSG{index:011d}' for index in range(AXIS_SIZE)],
        'match_type': ['direct'] * AXIS_SIZE,
    })


def deg_table(rows):
    """A compact published DEG table from ``(gene_index, log2fc, padj)`` rows."""
    frame = pd.DataFrame(rows, columns=['gene_index', 'log2fc', 'padj'])
    frame['contrast'] = 'fixture'
    return frame


# Two ICM sources whose DEG endpoints are mirrored: the level-1 sums of genes 3
# and 4 tie exactly, so the gene_index tie-break decides (3 before 4).
ICM_B = deg_table([(3, 4.0, 1e-6), (4, 2.0, 1e-3)])
ICM_C = deg_table([(4, 4.0, 1e-6), (3, 2.0, 1e-3)])
# A zero-DEG source: no gene passes padj < 0.05 and |log2FC| >= 2.
ZERO = deg_table([(1, 0.4, 0.5), (2, -0.3, 0.8), (3, 0.2, 0.4)])
# The AF extreme case: one informative source with exactly one DEG.
SINGLE = deg_table([(7, 3.5, 1e-4), (8, 0.4, 0.6)])
# A v2 published table: its pass-any flag marks gene 5 a DEG in a contrast the
# collapsed winning row would not flag (small effect, non-significant padj).
V2_TABLE = deg_table([(5, 0.4, 0.6), (6, 3.0, 1e-5)]).assign(
    deg=lambda frame: frame.gene_index.isin([5, 6]))

FIXTURE_CATEGORIES = {
    'nhf': {
        'label': 'nHF', 'tier': 1, 'evidence_type': 'cm_pseudobulk',
        'datasets': [
            {'dataset_id': 'nhf_a', 'tier': 1, 'evidence_type': 'cm_pseudobulk',
             'table': deg_table([(1, 3.0, 1e-4), (2, -2.5, 1e-5)])},
            {'dataset_id': 'nhf_zero', 'tier': 1, 'evidence_type': 'cm_pseudobulk',
             'table': ZERO},
        ],
    },
    'icm': {
        'label': 'ICM', 'tier': 1, 'evidence_type': 'cm_pseudobulk',
        'datasets': [
            {'dataset_id': 'icm_b', 'tier': 1, 'evidence_type': 'cm_pseudobulk',
             'table': ICM_B},
            {'dataset_id': 'icm_c', 'tier': 1, 'evidence_type': 'cm_pseudobulk',
             'table': ICM_C},
        ],
    },
    'af': {
        'label': 'AF', 'tier': 1, 'evidence_type': 'cm_pseudobulk',
        'datasets': [
            {'dataset_id': 'af_single', 'tier': 1, 'evidence_type': 'cm_pseudobulk',
             'table': SINGLE},
        ],
    },
    'empty': {
        'label': 'empty', 'tier': 2, 'evidence_type': 'bulk',
        'datasets': [
            {'dataset_id': 'empty_zero', 'tier': 2, 'evidence_type': 'bulk',
             'table': ZERO},
        ],
    },
    'v2': {
        'label': 'v2', 'tier': 2, 'evidence_type': 'bulk',
        'datasets': [
            {'dataset_id': 'v2_published', 'tier': 2, 'evidence_type': 'bulk',
             'table': V2_TABLE},
        ],
    },
}

# A declared dataset with no accepted disease-deg product: the catalog must
# report it as not_published rather than dropping it.
NOT_PUBLISHED = {'nhf': ['nhf_declared_missing']}
REGISTRY_CATEGORIES = {
    'nHF': {'directory': 'nhf', 'tier': 1, 'evidence_type': 'cm_pseudobulk'},
    'ICM': {'directory': 'icm', 'tier': 1, 'evidence_type': 'cm_pseudobulk'},
    'AF': {'directory': 'af', 'tier': 1, 'evidence_type': 'cm_pseudobulk'},
    'empty': {'directory': 'empty', 'tier': 2, 'evidence_type': 'bulk'},
    'v2': {'directory': 'v2', 'tier': 2, 'evidence_type': 'bulk'},
    # Declared but never published in the fixture: the unpublished record path.
    'unpublished': {'directory': 'unpublished', 'tier': 3, 'evidence_type': 'bulk'},
}


def write_fixture(root):
    """Write a hermetic repo tree with the axis and the fixture disease products."""
    root = Path(root)
    identity = root / 'stage2/01_bases/output/gene-identity'
    identity.mkdir(parents=True, exist_ok=True)
    axis_frame().to_csv(identity / 'genes.csv', index=False, lineterminator='\n')
    disease = root / driver.DISEASE_DEG_RELATIVE
    for directory, spec in FIXTURE_CATEGORIES.items():
        target = disease / directory
        target.mkdir(parents=True, exist_ok=True)
        outputs, records = {}, []
        for dataset in spec['datasets']:
            name = f"{dataset['dataset_id']}.deg.csv"
            path = target / name
            dataset['table'].to_csv(path, index=False, lineterminator='\n',
                                    float_format='%.17g')
            outputs[name] = {'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                             'rows': int(len(dataset['table']))}
            records.append({
                'dataset_id': dataset['dataset_id'],
                'category': spec['label'],
                'tier': dataset['tier'],
                'evidence_type': dataset['evidence_type'],
                'deg': name,
                'informative': True,
                'method': {'model': 'fixture'},
                'inputs': {'fixture.csv': hashlib.sha256(b'fixture').hexdigest()},
            })
        catalog = {
            'schema': disease_deg.CATALOG_SCHEMA,
            'category': spec['label'],
            'directory': directory,
            'tier': spec['tier'],
            'evidence_type': spec['evidence_type'],
            'axis_size': AXIS_SIZE,
            'axis': {'path': 'stage2/01_bases/output/gene-identity/genes.csv',
                     'sha256': hashlib.sha256(
                         (identity / 'genes.csv').read_bytes()).hexdigest(),
                     'size': AXIS_SIZE},
            'rank_contract': disease_deg.RANK_SCHEMA,
            'datasets': records,
            'rejected': [],
        }
        digest = disease_deg.write_catalog(target / 'catalog.json', catalog)
        outputs['catalog.json'] = {
            'sha256': hashlib.sha256((target / 'catalog.json').read_bytes()).hexdigest()}
        manifest = {
            'schema': disease_deg.RELEASE_SCHEMA,
            'category': spec['label'],
            'recipe': 'fixture',
            'axis': {'path': 'stage2/01_bases/output/gene-identity/genes.csv'},
            'inputs': {},
            'outputs': dict(sorted(outputs.items())),
            'counts': {'datasets': len(records)},
            'catalog_sha256': digest,
        }
        (target / 'manifest.json').write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    sources = root / 'fixture_sources.py'
    sources.write_text(_sources_module(), encoding='utf-8')
    return root, sources


def _sources_module():
    datasets = []
    for directory, spec in FIXTURE_CATEGORIES.items():
        for dataset in spec['datasets']:
            datasets.append({'dataset_id': dataset['dataset_id'],
                             'category': spec['label'],
                             'directory': directory})
        for dataset_id in NOT_PUBLISHED.get(directory, ()):
            datasets.append({'dataset_id': dataset_id, 'category': spec['label'],
                             'directory': directory})
    datasets.append({'dataset_id': 'unpublished_source', 'category': 'unpublished',
                     'directory': 'unpublished'})
    return (
        '"""Fixture source declaration generated by test_disease_fusion.py."""\n'
        f'CATEGORIES = {REGISTRY_CATEGORIES!r}\n'
        f'DATASETS = {datasets!r}\n'
        'def datasets(category):\n'
        '    return [record for record in DATASETS if record["category"] == category]\n'
        'def category_record(category):\n'
        '    return CATEGORIES[category]\n'
    )


def _strings(value):
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield from _strings(key)
            yield from _strings(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            yield from _strings(item)


def context_for(root, sources):
    """The driver context: fixture input root, fusion output under the fixture."""
    return driver.Context(output=Path(root) / 'stage2/01_bases/output/disease-fusion',
                          input_root=Path(root), sources=driver.load_sources(sources),
                          staging_root=Path(root) / 'staging')


# ---------------------------------------------------------------- source rank

class SourceDegRankTests(unittest.TestCase):
    """One published table -> the v2 DEG-only rank (both table generations)."""

    def test_v1_collapsed_rows_derive_the_v2_flag(self):
        result = disease_deg.source_deg_rank(deg_table([(1, 3.0, 1e-4),
                                                       (2, 0.5, 0.5),
                                                       (3, -2.0, 1e-3)]))
        self.assertEqual(result['n_deg'], 2)
        self.assertEqual(result['flag_source'], 'collapsed_winning_contrast')
        self.assertTrue(result['informative'])
        self.assertEqual(result['rank'].shape, (AXIS_SIZE,))
        self.assertEqual(result['rank'][1], 1)  # larger effect ranks first
        self.assertEqual(result['rank'][3], 2)
        self.assertEqual(result['rank'][2], 3)  # non-DEG gets n_DEG + 1
        self.assertEqual(result['rank'][0], 3)

    def test_published_v2_flag_is_authoritative(self):
        # The collapsed row is a non-significant small effect, but the published
        # pass-any flag says the gene is a DEG; the flag must win.
        result = disease_deg.source_deg_rank(V2_TABLE)
        self.assertEqual(result['n_deg'], 2)
        self.assertEqual(result['flag_source'], 'published_deg_flag')
        self.assertLessEqual(result['rank'][5], 2)
        self.assertLessEqual(result['rank'][6], 2)

    def test_zero_deg_source_is_flagged_non_informative(self):
        result = disease_deg.source_deg_rank(ZERO)
        self.assertEqual(result['n_deg'], 0)
        self.assertFalse(result['informative'])
        self.assertEqual(int(result['rank'].min()), 1)
        self.assertEqual(int(result['rank'].max()), 1)
        self.assertIn('is a DEG', result['reason'])

    def test_table_without_measurements_fails_closed(self):
        with self.assertRaises(ValueError):
            disease_deg.source_deg_rank(pd.DataFrame({'gene_index': [], 'log2fc': [],
                                                      'padj': []}))
        with self.assertRaises(ValueError):
            disease_deg.source_deg_rank(pd.DataFrame({'gene_index': [1]}))
        with self.assertRaises(ValueError):
            disease_deg.source_deg_rank(deg_table([(AXIS_SIZE, 3.0, 1e-4)]))
        with self.assertRaises(ValueError):
            disease_deg.source_deg_rank(deg_table([(1, 3.0, 1e-4),
                                                  (1, 2.0, 1e-4)]))


# ---------------------------------------------------------------- level 1

class LevelOneFusionTests(unittest.TestCase):
    """sum(rank / n_DEG), ascending rank, gene_index tie-break."""

    @staticmethod
    def source(dataset_id, rows, tier=1):
        derived = disease_deg.source_deg_rank(deg_table(rows))
        return {'dataset_id': dataset_id, 'tier': tier, 'rank': derived['rank'],
                'n_deg': derived['n_deg'], 'informative': derived['informative'],
                'informativeness_reason': derived['reason']}

    def test_sum_and_gene_index_tie_break(self):
        fusion = disease_deg.fuse_disease(
            [self.source('b', [(3, 4.0, 1e-6), (4, 2.0, 1e-3)]),
             self.source('c', [(4, 4.0, 1e-6), (3, 2.0, 1e-3)])])
        # genes 3 and 4 tie at 1/2 + 2/2 = 1.5, so gene_index decides.
        self.assertEqual(int(fusion['rank'][3]), 1)
        self.assertEqual(int(fusion['rank'][4]), 2)
        # A gene that is no DEG in either source contributes 2 * 3/2 = 3.0.
        np.testing.assert_allclose(fusion['sum'][0], 3.0)
        self.assertEqual(int(fusion['n_deg']), 2)
        self.assertEqual(fusion['deg2000'], 2)
        self.assertTrue(np.isfinite(fusion['sum']).all())

    def test_non_deg_term_is_n_deg_plus_one_over_n_deg(self):
        fusion = disease_deg.fuse_disease([self.source('a', [(1, 3.0, 1e-4), (9, 0.5, 0.5)])])
        self.assertEqual(fusion['n_deg'], 1)
        np.testing.assert_allclose(fusion['sum'][0], 2.0)  # (1+1)/1
        np.testing.assert_allclose(fusion['sum'][1], 1.0)  # 1/1
        np.testing.assert_allclose(fusion['sum'][2], 2.0)

    def test_zero_deg_source_is_excluded_and_recorded(self):
        fusion = disease_deg.fuse_disease(
            [self.source('b', [(3, 4.0, 1e-6), (9, 0.5, 0.5)]), self.source('zero', [(1, 0.4, 0.5)])])
        self.assertEqual(fusion['n_sources_declared'], 2)
        self.assertEqual(fusion['n_sources_informative'], 1)
        self.assertEqual(fusion['n_deg'], 1)
        self.assertTrue(np.isfinite(fusion['sum']).all())
        excluded = fusion['excluded_sources']
        self.assertEqual([entry['dataset_id'] for entry in excluded], ['zero'])
        self.assertEqual(excluded[0]['n_deg'], 0)
        self.assertIn('is a DEG', excluded[0]['reason'])
        ledger = {row['dataset_id']: row for row in fusion['sources']}
        self.assertFalse(ledger['zero']['informative'])
        self.assertEqual(ledger['zero']['n_deg_in_deg2000'], 0)
        self.assertTrue(ledger['b']['informative'])

    def test_zero_deg_never_poisons_a_sum_with_a_real_source(self):
        # The excluded NaN term must leave the informative source's sums intact.
        alone = disease_deg.fuse_disease([self.source('b', [(3, 4.0, 1e-6), (9, 0.5, 0.5)])])
        mixed = disease_deg.fuse_disease(
            [self.source('b', [(3, 4.0, 1e-6), (9, 0.5, 0.5)]),
             self.source('zero', [(1, 0.4, 0.5)])])
        np.testing.assert_array_equal(alone['sum'], mixed['sum'])
        np.testing.assert_array_equal(alone['rank'], mixed['rank'])

    def test_empty_union_is_representable(self):
        fusion = disease_deg.fuse_disease([self.source('zero', [(1, 0.4, 0.5)])])
        self.assertFalse(fusion['informative'])
        self.assertEqual(fusion['n_deg'], 0)
        self.assertEqual(fusion['deg2000'], 0)
        self.assertIsNone(fusion['rank'])
        self.assertIsNone(fusion['sum'])
        self.assertFalse(fusion['deg_union'].any())
        self.assertFalse(fusion['deg2000_mask'].any())
        self.assertIn('no informative source', fusion['reason'])
        rows = disease_deg.deg2000_rows(fusion)
        self.assertEqual(list(rows.columns), list(disease_deg.DEG2000_COLUMNS))
        self.assertTrue(rows.empty)

    def test_singleton_disease_publishes_its_one_gene(self):
        # AF's extreme case: one informative source, n_DEG = 1.
        fusion = disease_deg.fuse_disease([self.source('af', [(7, 3.5, 1e-4), (9, 0.5, 0.6)])])
        self.assertTrue(fusion['informative'])
        self.assertEqual(fusion['n_sources_informative'], 1)
        self.assertEqual(fusion['n_deg'], 1)
        self.assertEqual(fusion['deg2000'], 1)
        self.assertEqual(int(fusion['rank'][7]), 1)
        rows = disease_deg.deg2000_rows(fusion)
        self.assertEqual(rows.gene_index.tolist(), [7])
        self.assertEqual(rows.disease_rank.tolist(), [1])
        self.assertEqual(rows.disease_rank_sum.tolist(), [1.0])
        self.assertTrue(np.isfinite(fusion['sum']).all())

    def test_duplicate_source_fails_closed(self):
        source = self.source('a', [(1, 3.0, 1e-4), (9, 0.5, 0.5)])
        with self.assertRaises(ValueError):
            disease_deg.fuse_disease([source, source])
        with self.assertRaises(ValueError):
            disease_deg.fuse_disease([])

    def test_union_occupies_the_first_ranks(self):
        fusion = disease_deg.fuse_disease(
            [self.source('a', [(1, 3.0, 1e-4), (2, 2.5, 1e-3)]),
             self.source('b', [(3, 4.0, 1e-6), (9, 0.5, 0.5)])])
        union_ranks = fusion['rank'][fusion['deg_union']]
        self.assertEqual(sorted(union_ranks.tolist()),
                         list(range(1, fusion['n_deg'] + 1)))
        self.assertEqual(fusion['deg2000'], min(fusion['n_deg'],
                                                disease_deg.DEG2000_LIMIT))


class Deg2000TruncationTests(unittest.TestCase):
    """The declared 2,000 cut is exact: 2,199 -> 2,000, 120 -> 120, subset always."""

    @staticmethod
    def synthetic_source(n_deg, offset=0, n_genes=AXIS_SIZE):
        """A source whose genes ``offset..offset+n_deg-1`` are its DEGs, in order."""
        rank = np.full(n_genes, n_deg + 1, dtype=np.int64)
        genes = np.arange(offset, offset + n_deg, dtype=np.int64)
        rank[genes] = np.arange(1, n_deg + 1, dtype=np.int64)
        return {'dataset_id': f's{offset}', 'tier': 1, 'rank': rank,
                'n_deg': int(n_deg), 'informative': True}

    def test_truncation_at_exactly_the_declared_limit(self):
        fusion = disease_deg.fuse_disease([self.synthetic_source(2199)])
        self.assertEqual(fusion['n_deg'], 2199)
        self.assertEqual(fusion['deg2000'], disease_deg.DEG2000_LIMIT)
        self.assertEqual(int(fusion['deg2000_mask'].sum()), 2000)
        self.assertTrue((fusion['deg2000_mask'] <= fusion['deg_union']).all())
        rows = disease_deg.deg2000_rows(fusion)
        self.assertEqual(len(rows), 2000)
        self.assertEqual(rows.disease_rank.tolist(), list(range(1, 2001)))
        self.assertEqual(rows.disease_rank_sum.tolist(),
                         [rank / 2199 for rank in range(1, 2001)])

    def test_smaller_union_publishes_all_of_it(self):
        fusion = disease_deg.fuse_disease([self.synthetic_source(120)])
        self.assertEqual(fusion['n_deg'], 120)
        self.assertEqual(fusion['deg2000'], 120)
        rows = disease_deg.deg2000_rows(fusion)
        self.assertEqual(len(rows), 120)
        self.assertTrue(set(rows.gene_index) <= set(np.flatnonzero(
            fusion['deg_union']).tolist()))

    def test_boundary_at_the_limit_keeps_everything(self):
        fusion = disease_deg.fuse_disease([self.synthetic_source(2000)])
        self.assertEqual(fusion['deg2000'], 2000)
        self.assertEqual(len(disease_deg.deg2000_rows(fusion)), 2000)


# ---------------------------------------------------------------- catalog shape

def valid_catalog(**overrides):
    """A minimal category catalog carrying the v2 fusion fields."""
    catalog = {
        'schema': disease_deg.CATALOG_SCHEMA,
        'category': 'nHF',
        'directory': 'nhf',
        'tier': 1,
        'evidence_type': 'cm_pseudobulk',
        'axis_size': AXIS_SIZE,
        'catalog_sha256': '0' * 64,
        'rank_contract': disease_deg.RANK_SCHEMA,
        'n_deg': 2199,
        'deg2000': 2000,
        'tiering': {
            'tiers_entered': [1, 2],
            'tier_union': [{'tier': 1, 'n_deg_union': 900},
                           {'tier': 2, 'n_deg_union': 2199}],
            'stopped_at_tier': 2,
            'stop_reached': True,
            'limit': disease_deg.DEG2000_LIMIT,
        },
        'datasets': [
            {'dataset_id': 'a', 'tier': 1, 'informative': True, 'n_deg': 900,
             'n_deg_in_deg2000': 900},
            # #206: an nHF tier-2 record must declare the clinical, non-genetic basis.
            {'dataset_id': 'b', 'tier': 2, 'informative': True, 'n_deg': 1800,
             'n_deg_in_deg2000': 1400, 'genotype_status': 'not_declared',
             'non_genetic_model_basis': 'fixture: patient-derived iPSC-CM, no mutation '
                                        'introduced and no isogenic control'},
        ],
    }
    catalog.update(overrides)
    return catalog


class CategoryCatalogShapeTests(unittest.TestCase):
    """The v2 catalog fields are present, typed and internally consistent."""

    def test_valid_catalog_passes(self):
        catalog = valid_catalog()
        self.assertIs(disease_deg.validate_catalog(catalog), catalog)

    def test_missing_field_fails_closed(self):
        for key in ('tier', 'evidence_type', 'n_deg', 'deg2000', 'tiering'):
            catalog = valid_catalog()
            catalog.pop(key)
            with self.assertRaises(ValueError, msg=key):
                disease_deg.validate_catalog(catalog)

    def test_tier_must_be_one_of_three(self):
        for value in (0, 4, '1', True):
            with self.assertRaises(ValueError):
                disease_deg.validate_catalog(valid_catalog(tier=value))
        disease_deg.validate_catalog(valid_catalog(tier=3))

    def test_deg2000_must_be_min_of_union_and_limit(self):
        with self.assertRaises(ValueError):
            disease_deg.validate_catalog(valid_catalog(deg2000=1999))
        with self.assertRaises(ValueError):
            disease_deg.validate_catalog(valid_catalog(n_deg=120, deg2000=200))
        catalog = valid_catalog(n_deg=120, deg2000=120, tiering={
            'tiers_entered': [1],
            'tier_union': [{'tier': 1, 'n_deg_union': 120}],
            'stopped_at_tier': 1,
            'stop_reached': False,
            'limit': disease_deg.DEG2000_LIMIT,
        })
        catalog['datasets'] = [
            {'dataset_id': 'a', 'tier': 1, 'informative': True, 'n_deg': 90,
             'n_deg_in_deg2000': 90},
            {'dataset_id': 'b', 'tier': 1, 'informative': True, 'n_deg': 50,
             'n_deg_in_deg2000': 30},
        ]
        disease_deg.validate_catalog(catalog)

    def test_ledger_is_bounded_by_the_source_and_the_cut(self):
        catalog = valid_catalog()
        catalog['datasets'][0]['n_deg_in_deg2000'] = 901
        with self.assertRaises(ValueError):
            disease_deg.validate_catalog(catalog)
        catalog = valid_catalog()
        catalog['datasets'][1]['n_deg_in_deg2000'] = 2001
        with self.assertRaises(ValueError):
            disease_deg.validate_catalog(catalog)

    def test_union_bounds_and_coverage(self):
        with self.assertRaises(ValueError):
            disease_deg.validate_catalog(valid_catalog(
                n_deg=800, deg2000=800, tiering={
                    'tiers_entered': [1],
                    'tier_union': [{'tier': 1, 'n_deg_union': 800}],
                    'stopped_at_tier': 1, 'stop_reached': False,
                    'limit': disease_deg.DEG2000_LIMIT}))
        with self.assertRaises(ValueError):
            disease_deg.validate_catalog(valid_catalog(
                n_deg=2800, deg2000=2000, tiering={
                    'tiers_entered': [1, 2],
                    'tier_union': [{'tier': 1, 'n_deg_union': 1000},
                                   {'tier': 2, 'n_deg_union': 2800}],
                    'stopped_at_tier': 2, 'stop_reached': True,
                    'limit': disease_deg.DEG2000_LIMIT}))

    def test_tiering_must_agree_with_the_published_union(self):
        def tiering(**changes):
            record = dict(valid_catalog()['tiering'])
            record.update(changes)
            return record
        with self.assertRaises(ValueError):  # last union != published n_deg
            disease_deg.validate_catalog(valid_catalog(tiering=tiering(
                tier_union=[{'tier': 1, 'n_deg_union': 900},
                            {'tier': 2, 'n_deg_union': 2100}])))
        with self.assertRaises(ValueError):  # stopped_at_tier not the last tier
            disease_deg.validate_catalog(valid_catalog(tiering=tiering(
                stopped_at_tier=1)))
        with self.assertRaises(ValueError):  # stop_reached disagrees with the union
            disease_deg.validate_catalog(valid_catalog(tiering=tiering(
                stop_reached=False)))
        with self.assertRaises(ValueError):  # a lower tier already met the limit
            disease_deg.validate_catalog(valid_catalog(tiering={
                'tiers_entered': [1, 2],
                'tier_union': [{'tier': 1, 'n_deg_union': 2100},
                               {'tier': 2, 'n_deg_union': 2199}],
                'stopped_at_tier': 2, 'stop_reached': True,
                'limit': disease_deg.DEG2000_LIMIT}))
        # #206: a category whose evidence begins at tier 2 (the cardiotox classes)
        # legitimately enters there; descending tiers are still refused.
        cardiotox = valid_catalog(category='cardiotox_TOP2', tier=2, tiering={
            'tiers_entered': [2],
            'tier_union': [{'tier': 2, 'n_deg_union': 2199}],
            'stopped_at_tier': 2, 'stop_reached': True,
            'limit': disease_deg.DEG2000_LIMIT})
        cardiotox['datasets'] = [
            {'dataset_id': 'a', 'tier': 2, 'informative': True, 'n_deg': 900,
             'n_deg_in_deg2000': 900},
            {'dataset_id': 'b', 'tier': 2, 'informative': True, 'n_deg': 1800,
             'n_deg_in_deg2000': 1400},
        ]
        disease_deg.validate_catalog(cardiotox)
        with self.assertRaises(ValueError):  # tiers not entered in ascending order
            disease_deg.validate_catalog(valid_catalog(tiering={
                'tiers_entered': [3, 2],
                'tier_union': [{'tier': 3, 'n_deg_union': 1000},
                               {'tier': 2, 'n_deg_union': 2199}],
                'stopped_at_tier': 2, 'stop_reached': True,
                'limit': disease_deg.DEG2000_LIMIT}))
        with self.assertRaises(ValueError):  # shrinking cumulative union
            disease_deg.validate_catalog(valid_catalog(tiering={
                'tiers_entered': [1, 2],
                'tier_union': [{'tier': 1, 'n_deg_union': 2200},
                               {'tier': 2, 'n_deg_union': 2199}],
                'stopped_at_tier': 2, 'stop_reached': True,
                'limit': disease_deg.DEG2000_LIMIT}))

    def test_non_informative_source_is_outside_the_union_bounds(self):
        # A tied-effects source keeps its DEG table but is excluded from the
        # fusion, so its larger DEG set does not bound the published union.
        catalog = valid_catalog(n_deg=100, deg2000=100, tiering={
            'tiers_entered': [1],
            'tier_union': [{'tier': 1, 'n_deg_union': 100}],
            'stopped_at_tier': 1, 'stop_reached': False,
            'limit': disease_deg.DEG2000_LIMIT,
        })
        catalog['datasets'] = [
            {'dataset_id': 'kept', 'tier': 1, 'informative': True, 'n_deg': 100,
             'n_deg_in_deg2000': 100},
            {'dataset_id': 'tied', 'tier': 1, 'informative': False, 'n_deg': 500,
             'n_deg_in_deg2000': 0},
        ]
        disease_deg.validate_catalog(catalog)
        catalog['datasets'][1]['n_deg_in_deg2000'] = 1
        with self.assertRaises(ValueError):
            disease_deg.validate_catalog(catalog)
        catalog['datasets'][1].pop('informative')
        with self.assertRaises(ValueError):
            disease_deg.validate_catalog(catalog)

    def test_no_absolute_paths(self):
        catalog = valid_catalog()
        catalog['rejected'] = [{'candidate': '/home/someone/secret'}]
        with self.assertRaises(ValueError):
            disease_deg.validate_catalog(catalog)

    def test_build_catalog_publishes_and_validates_the_v2_fields(self):
        records = [{'dataset_id': 'a', 'tier': 1, 'informative': True},
                   {'dataset_id': 'b', 'tier': 2, 'informative': True,
                    'genotype_status': 'not_declared',
                    'non_genetic_model_basis': 'fixture: patient-derived iPSC-CM, '
                                               'no mutation introduced'}]
        catalog = disease_deg.build_catalog(
            'nHF', 'nhf', 1, 'cm_pseudobulk', pd.DataFrame(), '1' * 64, records,
            n_deg=2199, deg2000=2000, tiering=valid_catalog()['tiering'],
            source_deg={'a': {'n_deg': 900, 'n_deg_in_deg2000': 900},
                        'b': {'n_deg': 1800, 'n_deg_in_deg2000': 1400}})
        self.assertEqual(catalog['n_deg'], 2199)
        self.assertEqual(catalog['deg2000'], 2000)
        self.assertEqual(catalog['datasets'][0]['n_deg'], 900)
        self.assertEqual(catalog['datasets'][1]['n_deg_in_deg2000'], 1400)
        disease_deg.validate_catalog(catalog)
        with self.assertRaises(ValueError):  # partial fusion fields
            disease_deg.build_catalog('nHF', 'nhf', 1, 'bulk', pd.DataFrame(),
                                      '1' * 64, records, n_deg=2199)
        with self.assertRaises(ValueError):  # a declared dataset without a ledger
            disease_deg.build_catalog(
                'nHF', 'nhf', 1, 'bulk', pd.DataFrame(), '1' * 64, records,
                n_deg=2199, deg2000=2000, tiering=valid_catalog()['tiering'],
                source_deg={'a': {'n_deg': 900, 'n_deg_in_deg2000': 900}})


# ---------------------------------------------------------------- driver

class FusionDriverTests(unittest.TestCase):
    """The emitted per-disease fusion product: shape, exclusions, replay, hygiene."""

    def publish(self, root):
        root, sources = write_fixture(root)
        context = context_for(root, sources)
        result = driver.publish(context)
        return context, result

    def test_catalog_records_the_fusion_numbers_and_ledger(self):
        with tempfile.TemporaryDirectory() as temporary:
            context, result = self.publish(Path(temporary))
            catalog = disease_deg.DiseaseFusion(context.output).catalog()
            records = {record['directory']: record for record in catalog['diseases']}
            # nHF: two published sources, one of them zero-DEG (excluded).
            nhf = records['nhf']
            self.assertTrue(nhf['published'])
            self.assertEqual(nhf['declared_source_count'], 3)  # two + one missing
            self.assertEqual(nhf['informative_source_count'], 1)
            self.assertEqual(nhf['n_deg'], 2)
            self.assertEqual(nhf['deg2000'], 2)
            self.assertEqual([row['dataset_id'] for row in nhf['excluded_sources']],
                             ['nhf_zero'])
            self.assertIn('is a DEG', nhf['excluded_sources'][0]['reason'])
            self.assertEqual([row['dataset_id'] for row in nhf['not_published']],
                             ['nhf_declared_missing'])
            ledger = {row['dataset_id']: row for row in nhf['sources']}
            self.assertEqual(ledger['nhf_a']['n_deg'], 2)
            self.assertEqual(ledger['nhf_a']['n_deg_in_deg2000'], 2)
            self.assertEqual(ledger['nhf_zero']['n_deg'], 0)
            self.assertEqual(ledger['nhf_zero']['n_deg_in_deg2000'], 0)
            # ICM: two informative sources; the tie breaks by gene_index.
            icm = records['icm']
            self.assertEqual(icm['informative_source_count'], 2)
            self.assertEqual(icm['n_deg'], 2)
            self.assertEqual(icm['deg2000'], 2)
            # AF: one source, one DEG -- the degenerate singleton is well formed.
            af = records['af']
            self.assertEqual(af['n_deg'], 1)
            self.assertEqual(af['deg2000'], 1)
            rows = disease_deg.DiseaseFusion(context.output).deg2000('af')
            self.assertEqual(rows.gene_index.tolist(), [7])
            self.assertEqual(rows.disease_rank.tolist(), [1])
            # Empty union: published, header-only DEG2000 file.
            empty = records['empty']
            self.assertTrue(empty['published'])
            self.assertEqual(empty['n_deg'], 0)
            self.assertEqual(empty['deg2000'], 0)
            self.assertTrue(disease_deg.DiseaseFusion(
                context.output).deg2000('empty').empty)
            # Declared but unpublished category stays visible.
            unpublished = records['unpublished']
            self.assertFalse(unpublished['published'])
            self.assertEqual(unpublished['declared_source_count'], 1)
            # The v2 published flag is honoured.
            self.assertEqual(records['v2']['sources'][0]['flag_source'],
                             'published_deg_flag')
            self.assertEqual(result['published'], 5)
            self.assertEqual(result['diseases'], 6)

    def test_published_catalog_is_path_free_and_replays(self):
        with tempfile.TemporaryDirectory() as temporary:
            context, result = self.publish(Path(temporary))
            reader = disease_deg.DiseaseFusion(context.output)
            catalog = reader.catalog()
            manifest = reader.manifest()
            for payload in (catalog, manifest):
                for value in _strings(payload):
                    self.assertFalse(value.startswith('/'), value)
            self.assertEqual(manifest['catalog_sha256'], catalog['catalog_sha256'])
            replay = driver.replay(context)
            self.assertEqual(replay['result'], 'PASS')
            for record in catalog['diseases']:
                if not record['published']:
                    continue
                rows = reader.deg2000(record['directory'])
                self.assertEqual(list(rows.columns),
                                 list(disease_deg.DEG2000_COLUMNS))
                self.assertEqual(len(rows), record['deg2000'])

    def test_publish_is_deterministic_across_roots(self):
        with tempfile.TemporaryDirectory() as first, \
                tempfile.TemporaryDirectory() as second:
            context_a, _ = self.publish(Path(first))
            context_b, _ = self.publish(Path(second))
            for name in ('catalog.json', 'manifest.json'):
                self.assertEqual((context_a.output / name).read_bytes(),
                                 (context_b.output / name).read_bytes())

    def test_tampering_is_detected(self):
        with tempfile.TemporaryDirectory() as temporary:
            context, _ = self.publish(Path(temporary))
            target = context.output / 'nhf.deg2000.csv'
            original = target.read_text()
            target.write_text(original + '\n')
            with self.assertRaises(ValueError):
                driver.verify_fusion(context.output)
            target.write_text(original)
            self.assertEqual(driver.verify_fusion(context.output)['result'], 'PASS')

    def test_committed_fusion_catalog_replays_from_the_tracked_products(self):
        self.assertTrue((RESULT_ROOT / 'catalog.json').is_file(),
                        'the disease-fusion product must be committed')
        context = driver.Context(output=RESULT_ROOT, input_root=REPO,
                                 sources=disease_deg_sources)
        reader = disease_deg.DiseaseFusion(RESULT_ROOT)
        catalog = reader.catalog()
        self.assertEqual(catalog['limit'], disease_deg.DEG2000_LIMIT)
        self.assertEqual(catalog['axis']['size'], AXIS_SIZE)
        directories = [record['directory'] for record in catalog['diseases']]
        published = [record['directory'] for record in catalog['diseases']
                     if record['published']]
        self.assertTrue(set(('nhf', 'icm', 'genetic_cardiomyopathy',
                             'cardiotox_TOP2', 'cardiotox_other')) <= set(published))
        for directory in published:
            rows = reader.deg2000(directory)
            record = reader.record(directory)
            self.assertEqual(len(rows), record['deg2000'])
            if record['deg2000']:
                self.assertEqual(rows.disease_rank.tolist(),
                                 list(range(1, record['deg2000'] + 1)))
        for directory in directories:
            record = reader.record(directory)
            self.assertEqual(record['deg2000'] or 0,
                             min(record['n_deg'] or 0, disease_deg.DEG2000_LIMIT))
        # #206 declares two nHF tier-3 sources (GSE116250, GSE141910) whose
        # disease-deg products are published by #208, so the committed fusion
        # catalog's declared/not_published roster deliberately lags the registry
        # until that republish re-runs this driver. Rebuild from the same accepted
        # products and require every fused number and every published ledger to be
        # unchanged; the only tolerated delta is the declared roster, pinned here.
        with tempfile.TemporaryDirectory() as temporary:
            driver.build(context, Path(temporary))
            rebuilt = disease_deg.read_catalog(Path(temporary) / 'catalog.json')
        accepted = {record['directory']: record for record in catalog['diseases']}
        replayed = {record['directory']: record for record in rebuilt['diseases']}
        self.assertEqual(set(accepted), set(replayed))
        for directory, record in accepted.items():
            other = replayed[directory]
            for key in ('published', 'tier', 'evidence_type', 'informative_source_count',
                        'n_deg', 'deg2000', 'deg2000_file', 'sources',
                        'excluded_sources', 'tiers_present', 'reason'):
                self.assertEqual(record[key], other[key], (directory, key))
            if directory != 'nhf':
                self.assertEqual(record['not_published'], other['not_published'],
                                 directory)
        lagged = ({entry['dataset_id'] for entry in replayed['nhf']['not_published']}
                  - {entry['dataset_id'] for entry in accepted['nhf']['not_published']})
        self.assertEqual(lagged, {'gse116250__hf_lv_bulk', 'gse141910__hf_lv_bulk'},
                         'the #206 roster additions must be the only fusion-catalog lag')
        self.assertEqual(replayed['nhf']['declared_source_count'],
                         accepted['nhf']['declared_source_count'] + 2)


if __name__ == '__main__':
    unittest.main()
