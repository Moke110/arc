#!/usr/bin/env python
"""Tier order, the axis-resolved stop rule and the six-disease roster (#206).

Covers ADR 0071 decision 2 from the outside:

* each category's ordered tier ladder comes from its own declared source records
  (``nHF`` 1 -> 2 -> 3, ``ICM`` 1, ``AF`` 1 -> 3 with no tier-2 source, the
  cardiotox classes starting at tier 2);
* ``tier_ledger`` measures the cumulative DEG union **resolved to the frozen
  axis** after each tier, stops at the first tier that reaches 2,000, records the
  stop point, and refuses a source in a deeper tier after a tier that already
  reached the limit;
* the clinical tier-2 gate: ``nHF``/``ICM``/``AF`` tier-2 evidence must declare a
  non-genetic model basis; a mutation-modelled record is rejected there but
  accepted under ``genetic_cardiomyopathy``;
* the catalog carries the computed ``tiering`` record and the declared
  ``tier_sources`` ladder, and the record agrees with the published union.

The fixtures are synthetic DEG sets over the frozen-axis index space; the
registry assertions read only the tracked source declaration, never an external
payload.
"""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

TESTS = Path(__file__).resolve().parent
EXECUTION = TESTS.parent
REPO = next(path for path in Path(__file__).resolve().parents if (path / ".git").exists())
sys.path.insert(0, str(EXECUTION))
import public_api  # noqa: F401  (binds the single public implementation)
import disease_deg_sources  # noqa: E402
from bases import disease_deg  # noqa: E402

LIMIT = disease_deg.DEG2000_LIMIT


def dataset(dataset_id, tier, **extra):
    record = {'dataset_id': dataset_id, 'tier': tier}
    record.update(extra)
    return record


def measured(*genes, informative=None):
    genes = sorted(int(gene) for gene in genes)
    entry = {'deg': genes, 'n_deg': len(genes)}
    if informative is not None:
        entry['informative'] = informative
    return entry


def ledger(category, datasets, sets, limit=LIMIT):
    return disease_deg.tier_ledger(category, datasets, sets, limit=limit)


class TierLadderTests(unittest.TestCase):
    """The ladder is the declaration: ascending tiers, declared order inside each."""

    def test_ladder_groups_sources_by_their_declared_tier(self):
        ladder = disease_deg.tier_ladder([
            dataset('a', 1), dataset('b', 3), dataset('c', 1), dataset('d', 2)])
        self.assertEqual(ladder, [{'tier': 1, 'sources': ['a', 'c']},
                                  {'tier': 2, 'sources': ['d']},
                                  {'tier': 3, 'sources': ['b']}])

    def test_a_missing_tier_is_absent_not_padded(self):
        ladder = disease_deg.tier_ladder([dataset('a', 1), dataset('b', 3)])
        self.assertEqual([rung['tier'] for rung in ladder], [1, 3])

    def test_a_category_may_begin_at_tier_two(self):
        ladder = disease_deg.tier_ladder([dataset('a', 2), dataset('b', 2)])
        self.assertEqual([rung['tier'] for rung in ladder], [2])

    def test_out_of_contract_declarations_fail_closed(self):
        for value in (0, 4, '1', True, None):
            with self.assertRaises(ValueError, msg=value):
                disease_deg.tier_ladder([dataset('a', value)])
        with self.assertRaises(ValueError):  # repeated source
            disease_deg.tier_ladder([dataset('a', 1), dataset('a', 2)])
        with self.assertRaises(ValueError):  # no source
            disease_deg.tier_ladder([])
        with self.assertRaises(ValueError):  # missing id
            disease_deg.tier_ladder([dict(tier=1)])

    def test_the_six_category_registry_declares_its_tiers(self):
        expected = {
            'nHF': [1, 2, 3],
            'ICM': [1, 3],
            'AF': [1, 3],
            'genetic_cardiomyopathy': [1, 2],
        }
        for category, tiers in expected.items():
            ladder = disease_deg.tier_ladder(disease_deg_sources.datasets(category))
            self.assertEqual([rung['tier'] for rung in ladder], tiers, category)
        by_tier = {rung['tier']: rung['sources']
                   for rung in disease_deg.tier_ladder(disease_deg_sources.datasets('nHF'))}
        self.assertEqual(len(by_tier[1]), 4)
        self.assertIn('gse298023__dcm_cm', by_tier[1])
        self.assertEqual(len(by_tier[2]), 1)
        self.assertEqual(len(by_tier[3]), 3)
        self.assertIn('gse165303__dcm_lv_bulk', by_tier[3])
        self.assertIn('gse116250__hf_lv_bulk', by_tier[3])
        self.assertIn('gse141910__hf_lv_bulk', by_tier[3])
        af = disease_deg.tier_ladder(disease_deg_sources.datasets('AF'))
        self.assertEqual([len(rung['sources']) for rung in af], [3, 2])


class TierStopRuleTests(unittest.TestCase):
    """The four acceptance cases plus the refusals the rule must make."""

    def test_tier_one_short_enters_tier_two(self):
        datasets = [dataset('a', 1), dataset('b', 2)]
        record = ledger('genetic_cardiomyopathy', datasets,
                        {'a': measured(1, 2), 'b': measured(3)})['tiering']
        self.assertEqual(record['tiers_entered'], [1, 2])
        self.assertEqual([row['n_deg_union'] for row in record['tier_union']], [2, 3])
        self.assertEqual(record['collection_policy'], disease_deg.COLLECTION_POLICY)

    def test_tier_one_at_or_above_the_limit_stops_the_descent(self):
        datasets = [dataset('a', 1)]
        record = ledger('genetic_cardiomyopathy', datasets,
                        {'a': measured(*range(LIMIT))})['tiering']
        self.assertEqual(record['tiers_entered'], [1])
        self.assertEqual(record['collection_policy'], disease_deg.COLLECTION_POLICY)
        self.assertEqual(record['tier_union'][-1]['n_deg_union'], LIMIT)

    def test_tier_two_short_enters_tier_three(self):
        datasets = [dataset('a', 1), dataset('b', 2), dataset('c', 3)]
        record = ledger('genetic_cardiomyopathy', datasets,
                        {'a': measured(1), 'b': measured(2), 'c': measured(3)})['tiering']
        self.assertEqual(record['tiers_entered'], [1, 2, 3])
        self.assertEqual([row['n_deg_union'] for row in record['tier_union']], [1, 2, 3])
        self.assertEqual(record['collection_policy'], disease_deg.COLLECTION_POLICY)

    def test_a_source_after_two_thousand_is_collected(self):
        result = ledger('genetic_cardiomyopathy',
                        [dataset('done', 1), dataset('deeper', 2)],
                        {'done': measured(*range(LIMIT)), 'deeper': measured(LIMIT + 1)})
        self.assertEqual(result['n_deg'], 2001)
        self.assertEqual(result['tiering']['tiers_entered'], [1, 2])

    def test_a_tier_between_two_entered_tiers_is_not_back_filled(self):
        datasets = [dataset('a', 1), dataset('c', 3)]
        record = ledger('AF', datasets,
                        {'a': measured(7), 'c': measured()})['tiering']
        self.assertEqual(record['tiers_entered'], [1, 3])
        self.assertEqual([row['n_deg_union'] for row in record['tier_union']], [1, 1])

    def test_non_informative_sources_are_recorded_but_add_no_gene(self):
        datasets = [dataset('a', 1), dataset('b', 3)]
        result = ledger('AF', datasets,
                        {'a': measured(11, 12), 'b': measured(13, 14, informative=False)})
        self.assertEqual([row['n_deg_union'] for row in result['tiering']['tier_union']],
                         [2, 2])
        self.assertEqual([row['informative'] for row in result['sources']],
                         [True, False])
        self.assertEqual(result['n_deg'], 2)
        self.assertEqual(result['deg2000'], 2)

    def test_the_union_is_axis_resolved_not_genome_wide(self):
        # A source with 3,000 genome-wide DEGs of which only 1,600 resolve to the
        # frozen axis contributes 1,600 here: the off-axis genes are structurally
        # unreachable (they can never enter DEG2000 or the axis), so they cannot
        # count toward the 2,000 stop test. Measured axis-resolved, this source
        # does not stop collection; the same signal on 2,000 axis genes would.
        short = ledger('genetic_cardiomyopathy', [dataset('a', 1)],
                       {'a': measured(*range(1600))})['tiering']
        self.assertEqual(short['collection_policy'], disease_deg.COLLECTION_POLICY)
        self.assertEqual(short['tier_union'][-1]['n_deg_union'], 1600)
        reached = ledger('genetic_cardiomyopathy', [dataset('a', 1)],
                         {'a': measured(*range(LIMIT))})['tiering']
        self.assertEqual(reached['collection_policy'], disease_deg.COLLECTION_POLICY)

    def test_boolean_axis_mask_and_gene_index_vector_agree(self):
        mask = np.zeros(disease_deg.AXIS_SIZE, dtype=bool)
        mask[[3, 5, 7]] = True
        by_mask = ledger('AF', [dataset('a', 1)], {'a': {
            'deg': mask, 'n_deg': 3, 'informative': True}})
        by_index = ledger('AF', [dataset('a', 1)], {'a': measured(3, 5, 7)})
        self.assertEqual(by_mask['n_deg'], by_index['n_deg'])

    def test_measured_sets_must_cover_exactly_the_declared_sources(self):
        with self.assertRaises(ValueError):  # declared without a measurement
            ledger('AF', [dataset('a', 1)], {})
        with self.assertRaises(ValueError):  # measured but not declared
            ledger('AF', [dataset('a', 1)], {'a': measured(1), 'b': measured(2)})
        with self.assertRaises(ValueError):  # n_deg disagrees with the set
            ledger('AF', [dataset('a', 1)], {'a': {'deg': [1, 2], 'n_deg': 1}})
        with self.assertRaises(ValueError):  # gene outside the frozen axis
            ledger('AF', [dataset('a', 1)], {'a': measured(disease_deg.AXIS_SIZE)})
        with self.assertRaises(ValueError):  # repeated gene
            ledger('AF', [dataset('a', 1)], {'a': {'deg': [1, 1], 'n_deg': 2}})

    def test_the_stop_record_publishes_the_axis_resolved_criterion(self):
        record = ledger('AF', [dataset('a', 1)], {'a': measured(1)})['tiering']
        self.assertEqual(record['criterion'], disease_deg.TIER_CRITERION)
        self.assertIn('axis-resolved, not genome-wide', record['rule'])
        self.assertNotIn('limit', record)


class ClinicalTierBasisTests(unittest.TestCase):
    """ADR 0071 decision 2: tier-2 evidence in a clinical category is non-genetic."""

    basis = ('patient-derived iPSC-CM lines compared with normal iPSC-CM lines; no '
             'mutation was introduced and no isogenic control is deposited')

    def test_a_declared_non_genetic_basis_is_accepted(self):
        record = dataset('ok', 2, genotype_status='not_declared',
                         non_genetic_model_basis=self.basis)
        for category in disease_deg.CLINICAL_CATEGORIES:
            result = ledger(category, [record], {'ok': measured(1)})
            self.assertEqual(result['tiering']['tiers_entered'], [2], category)

    def test_a_tier_two_record_without_a_basis_is_rejected(self):
        for category in disease_deg.CLINICAL_CATEGORIES:
            with self.assertRaisesRegex(ValueError, 'non-genetic model basis'):
                ledger(category, [dataset('bad', 2)], {'bad': measured(1)})

    def test_a_mutation_modelled_tier_two_record_is_rejected(self):
        pathological = dataset('bad', 2, genotype_status='pathogenic',
                               non_genetic_model_basis=self.basis)
        isogenic = dataset('bad', 2, isogenic_exception={'applies': True},
                           non_genetic_model_basis=self.basis)
        for category in disease_deg.CLINICAL_CATEGORIES:
            for record in (pathological, isogenic):
                with self.assertRaisesRegex(ValueError, 'genetic_cardiomyopathy'):
                    ledger(category, [record], {'bad': measured(1)})

    def test_the_same_mutation_modelled_record_belongs_to_genetic_cardiomyopathy(self):
        record = dataset('g', 2, genotype_status='pathogenic')
        result = ledger('genetic_cardiomyopathy', [record], {'g': measured(1)})
        self.assertEqual(result['tiering']['tiers_entered'], [2])

    def test_a_tier_one_record_is_not_gated(self):
        result = ledger('nHF', [dataset('a', 1)], {'a': measured(1)})
        self.assertEqual(result['tiering']['tiers_entered'], [1])

    def test_the_registry_tier_two_nhf_record_declares_its_basis(self):
        (record,) = [row for row in disease_deg_sources.datasets('nHF')
                     if row['tier'] == 2]
        self.assertEqual(record['dataset_id'], 'gse35108__dcm_ipsc_cm')
        self.assertIn('no mutation', record['non_genetic_model_basis'])
        disease_deg.validate_clinical_tier_basis('nHF', record)


class CatalogTieringTests(unittest.TestCase):
    """The category catalog publishes the computed record and its ladder."""

    def catalog(self, category, datasets, sets, tier, directory=None, **extra):
        result = ledger(category, datasets, sets)
        records = []
        for row in result['sources']:
            records.append({'dataset_id': row['dataset_id'], 'tier': row['tier'],
                            'informative': row['informative'],
                            **{key: value for key, value in extra.get('records', {})
                               .get(row['dataset_id'], {}).items()}})
        counts = {row['dataset_id']: {'n_deg': row['n_deg'], 'n_deg_in_deg2000':
                                      min(row['n_deg'], result['deg2000'])}
                  for row in result['sources']}
        return disease_deg.build_catalog(
            category, directory or category.lower(), tier, 'bulk',
            pd.DataFrame(), '1' * 64, records,
            n_deg=result['n_deg'], deg2000=result['deg2000'],
            tiering=result['tiering'], tier_sources=result['tier_sources'],
            source_deg=counts)

    def test_published_record_carries_tiering_and_the_ladder(self):
        catalog = self.catalog('AF', [dataset('a', 1), dataset('c', 3)],
                               {'a': measured(1), 'c': measured()}, tier=1)
        self.assertEqual(catalog['tiering']['tiers_entered'], [1, 3])
        self.assertEqual(catalog['tier_sources'],
                         [{'tier': 1, 'sources': ['a']},
                          {'tier': 3, 'sources': ['c']}])
        self.assertEqual(catalog['tiering']['criterion'], disease_deg.TIER_CRITERION)

    def test_a_ladder_that_disagrees_with_the_record_is_rejected(self):
        catalog = self.catalog('AF', [dataset('a', 1), dataset('c', 3)],
                               {'a': measured(1), 'c': measured()}, tier=1)
        catalog['tier_sources'] = [{'tier': 1, 'sources': ['a', 'c']}]
        with self.assertRaises(ValueError):
            disease_deg.validate_catalog(catalog)

    def test_a_ladder_source_at_the_wrong_tier_is_rejected(self):
        catalog = self.catalog('AF', [dataset('a', 1), dataset('c', 3)],
                               {'a': measured(1), 'c': measured()}, tier=1)
        catalog['tier_sources'] = [{'tier': 1, 'sources': ['a']},
                                   {'tier': 3, 'sources': ['c']}]
        catalog['datasets'][1]['tier'] = 2
        with self.assertRaises(ValueError):
            disease_deg.validate_catalog(catalog)

    def test_cardiotox_class_tiering_starts_at_tier_two(self):
        catalog = self.catalog('cardiotox_TOP2', [dataset('t', 2)],
                               {'t': measured(*range(LIMIT + 1))}, tier=2)
        self.assertEqual(catalog['tiering']['tiers_entered'], [2])
        self.assertEqual(catalog['tiering']['collection_policy'], disease_deg.COLLECTION_POLICY)
        self.assertEqual(catalog['deg2000'], LIMIT)

    def test_clinical_tier_two_records_are_gated_in_the_catalog(self):
        with self.assertRaises(ValueError):
            ledger('nHF', [dataset('a', 2, genotype_status='pathogenic')],
                   {'a': measured(1)})
        catalog = self.catalog('genetic_cardiomyopathy', [dataset('a', 2)],
                               {'a': measured(1)}, tier=2)
        catalog['datasets'][0]['genotype_status'] = 'pathogenic'
        self.assertEqual(disease_deg.validate_catalog(catalog), catalog)


class LegacyLocalDegAdapterTests(unittest.TestCase):
    """The pre-#187 HF items deposit their own column names (#206 declarations)."""

    def frame(self):
        return pd.DataFrame({
            'gene_symbol': ['A1BG', 'A1CF'],
            'contrast': ['HF_vs_NF', 'HF_vs_NF'],
            'log2FC': [2.5, -3.0],
            'stat': [4.0, -5.0],
            'p_value': [0.0001, 0.0002],
            'padj': [0.001, 0.002],
            'mean_expr': [6.0, 7.0],
        })

    def test_declared_aliases_fill_the_canonical_contract(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'legacy_deg.csv'
            self.frame().to_csv(path, index=False)
            frame = disease_deg.read_local_deg(
                path, aliases={'log2FC': 'log2fc', 'mean_expr': 'base_mean'})
        self.assertEqual(tuple(frame.columns[:1]), ('gene_symbol',))
        self.assertIn('source_gene_id', frame.columns)
        self.assertEqual(frame.source_gene_id.tolist(), frame.gene_symbol.tolist())
        self.assertEqual(frame.log2fc.tolist(), [2.5, -3.0])
        self.assertEqual(frame.base_mean.tolist(), [6.0, 7.0])

    def test_an_unaliased_legacy_table_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / 'legacy_deg.csv'
            self.frame().to_csv(path, index=False)
            with self.assertRaises(ValueError):
                disease_deg.read_local_deg(path)

    def test_the_two_new_nhf_tier_three_records_use_the_adapter(self):
        for dataset_id in ('gse116250__hf_lv_bulk', 'gse141910__hf_lv_bulk'):
            record = next(row for row in disease_deg_sources.DATASETS
                          if row['dataset_id'] == dataset_id)
            self.assertEqual(record['category'], 'nHF')
            self.assertEqual(record['tier'], 3)
            self.assertEqual(record['evidence_type'], 'bulk')
            self.assertEqual(record['local_deg_columns'],
                             {'log2FC': 'log2fc', 'mean_expr': 'base_mean'})
            self.assertIn('not cardiomyocyte-specific', record['model_basis'])
            self.assertNotIn('non_genetic_model_basis', record)


if __name__ == '__main__':
    unittest.main()
