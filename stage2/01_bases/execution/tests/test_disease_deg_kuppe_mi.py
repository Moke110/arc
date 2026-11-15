#!/usr/bin/env python
"""Kuppe-MI: the ICM tier-1 CM-pseudobulk gene-axis source (#202).

Hermetic: only the item's *tracked* artefacts are read — ``local/deg/kuppe-mi_deg.json``,
``local/pseudobulk/kuppe-mi_cm_pseudobulk.json``, ``author_log.json``,
``manifest.sha256``, the ``script/`` declarations (imported, never executed) and the
shared ``disease_deg_sources`` registry. The 1.86 GB ``source/`` object and the
gitignored ``local/`` payloads are not required, so the suite runs on a fresh clone.

What is asserted is the externally visible contract of the source:

* the registry declares it in ``ICM`` at tier 1 with the IZ-vs-CTRL contrast, the
  4-vs-4 donor arms, its own control group and an undeclared time window;
* the tracked summaries carry the *verified* deposit design (cellxgene object id and
  version, the CM gate, the authors' zone column, 25 donor-zone profiles), the
  500-CM-nuclei profile floor and its dropped profiles, the one joint
  ``~0 + zone`` design with donor blocking, and the published-vs-artifact
  contrast split;
* the DEG table's declared interface (path, digest, columns) and the frozen-axis
  cross-check stay self-consistent;
* ``fetch.sh``'s byte pin, the offline ``verify.sh``, the manifest and the
  author log agree with each other and with the item's summary.
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

ITEM = 'Kuppe-MI'
ITEM_DIR = REPO / 'disease_data' / ITEM
DATASET = 'kuppe-mi__iz'
CONTRAST = 'IZ-CTRL'
ARTIFACT_CONTRASTS = ['BZ-CTRL', 'RZ-CTRL', 'FZ-CTRL']
SOURCE_BYTES = 1862945764
AXIS_SIZE = 18164
IZ_DONORS = ['P10', 'P15', 'P2', 'P9']
CONTROL_DONORS = ['P1', 'P17', 'P7', 'P8']


def load_item_module(name):
    path = ITEM_DIR / 'script' / f'{name}.py'
    loader = SourceFileLoader(f'kuppe_{name}', str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def read_json(path):
    return json.loads(Path(path).read_text())


def summary():
    return read_json(ITEM_DIR / 'local/deg/kuppe-mi_deg.json')


def pseudobulk_summary():
    return read_json(ITEM_DIR / 'local/pseudobulk/kuppe-mi_cm_pseudobulk.json')


def author_log():
    return read_json(ITEM_DIR / 'author_log.json')


def declaration():
    return next(record for record in disease_deg_sources.DATASETS
                if record['dataset_id'] == DATASET)


class DeclarationTest(unittest.TestCase):
    """The shared ICM roster gains exactly the Kuppe IZ dataset."""

    def test_declared_in_icm_at_tier_one(self):
        record = declaration()
        self.assertEqual(record['category'], 'ICM')
        self.assertEqual(record['tier'], 1)
        self.assertEqual(record['evidence_type'], 'cm_pseudobulk')
        self.assertEqual(record['source_item'], f'disease_data/{ITEM}')
        self.assertEqual(record['local_deg'], 'local/deg/kuppe-mi_deg.csv')
        self.assertEqual(record['local_summary'], 'local/deg/kuppe-mi_deg.json')
        self.assertEqual(record['cells_per_donor'], 'local/pseudobulk/kuppe-mi_profiles.csv')
        self.assertEqual(record['contrasts'], [CONTRAST])
        self.assertEqual((record['n_case'], record['n_control']), (4, 4))
        self.assertEqual(record['shared_control_group'], 'kuppe-mi:normal')
        self.assertFalse(record['exploratory'])

    def test_time_window_is_undeclared_not_asserted(self):
        record = declaration()
        self.assertIsNone(record['time_window_basis'])
        self.assertIn('time-since-infarct', record['time_window_note'])
        self.assertIn('never asserted', record['time_window_note'])

    def test_own_control_group_is_declared_at_category_level(self):
        groups = disease_deg_sources.category_record('ICM')['shared_control_group']
        self.assertIn('kuppe-mi:normal', groups)
        # The #187/#203 ICM membership is untouched by the additive entry.
        self.assertIn('gse270788:donor', groups)
        ids = [record['dataset_id'] for record in disease_deg_sources.datasets('ICM')]
        self.assertIn('gse270788__icm', ids)
        self.assertIn('gse270788__ami', ids)
        self.assertIn(DATASET, ids)

    def test_dataset_id_is_unique(self):
        ids = [record['dataset_id'] for record in disease_deg_sources.DATASETS]
        self.assertEqual(len(ids), len(set(ids)))


class SummaryTest(unittest.TestCase):
    """The tracked DEG summary states the design verified against the deposit."""

    def test_common_contract(self):
        record = summary()
        self.assertEqual(record['schema'], 'stage2-item-deg-summary-v1')
        self.assertEqual(record['item'], ITEM)
        self.assertEqual(record['category'], 'ICM')
        self.assertEqual(record['tier'], 1)
        self.assertEqual(record['published_contrast'], CONTRAST)
        self.assertEqual(sorted(record['artifact_contrasts']), sorted(ARTIFACT_CONTRASTS))
        self.assertIn('Cardiomyocyte', record['gate'])
        self.assertEqual(record['zone_column'], 'major_labl')
        self.assertEqual(sorted(record['contrasts']), [CONTRAST])
        self.assertEqual(record['contrasts'][CONTRAST]['n_case'], 4)
        self.assertEqual(record['contrasts'][CONTRAST]['n_control'], 4)
        self.assertEqual(sorted(record['contrasts'][CONTRAST]['case_donors']), IZ_DONORS)
        self.assertEqual(sorted(record['contrasts'][CONTRAST]['control_donors']),
                         CONTROL_DONORS)

    def test_profile_floor_is_recorded(self):
        floor = summary()['profile_floor']
        self.assertEqual(floor['rule'], 'profile n_cm >= 500')
        self.assertEqual(floor['n_profiles_total'], 25)
        self.assertEqual(floor['n_profiles_fitted'], 18)
        dropped = floor['dropped_profiles']
        self.assertEqual(len(dropped), floor['n_profiles_total'] - floor['n_profiles_fitted'])
        self.assertTrue(all(entry['n_cm'] < 500 for entry in dropped))
        self.assertEqual(floor['zone_profile_counts'],
                         {'CTRL': 4, 'IZ': 4, 'BZ': 3, 'RZ': 5, 'FZ': 2})
        # The floor removes sparse donor-zone profiles; every fitted donor keeps
        # at least one profile and no donor appears twice in one zone.
        self.assertEqual(len({entry['donor'] for entry in dropped}),
                         len({entry['profile'].split('|')[0] for entry in dropped}))

    def test_samples_and_arms(self):
        samples = summary()['samples']
        self.assertEqual(samples['n'], 18)
        self.assertEqual(samples['group_counts'],
                         {'CTRL': 4, 'IZ': 4, 'BZ': 3, 'RZ': 5, 'FZ': 2})
        self.assertEqual(sorted(samples['donors_per_arm']['IZ']), IZ_DONORS)
        self.assertEqual(sorted(samples['donors_per_arm']['CTRL']), CONTROL_DONORS)
        per_profile = samples['cm_nuclei_per_profile']
        self.assertEqual(len(per_profile), 18)
        self.assertTrue(all(count >= 500 for count in per_profile.values()))
        self.assertEqual(sum(per_profile.values()), sum(samples['cm_nuclei_per_profile'].values()))
        # The four IZ donors' IZ profiles are the ones the audit measured.
        for donor, n_cm in {'P2': 1140, 'P9': 1952, 'P10': 999, 'P15': 582}.items():
            self.assertEqual(per_profile[f'{donor}|IZ'], n_cm, donor)

    def test_design_declares_donor_blocking(self):
        design = summary()['design']
        self.assertEqual(design['kind'], 'counts')
        self.assertEqual(design['formula'], '~0 + zone')
        self.assertEqual(design['factor'], 'zone')
        self.assertEqual(design['block'], 'donor')
        self.assertIn(design['normalisation'], ('voom', 'tmm'))
        self.assertEqual(sorted(contrast['label'] for contrast in design['contrasts']),
                         sorted([CONTRAST] + ARTIFACT_CONTRASTS))
        self.assertEqual(design['filter']['min_count'], 10)
        self.assertEqual(design['filter']['min_total_count'], 15)
        ambiguities = summary()['design_ambiguities']
        self.assertIn('duplicateCorrelation', ambiguities['repeated_donors'])
        self.assertIn('only IZ-CTRL is published', ambiguities['published_contrast'])
        self.assertIn('not_declared', ambiguities['time_window'])

    def test_deg_counts_are_consistent(self):
        counts = summary()['deg_counts']
        published = counts['published_contrast']
        self.assertEqual(published['n_up'] + published['n_down'], published['n_v2_deg'])
        self.assertGreaterEqual(published['n_padj_lt_0.05'], published['n_v2_deg'])
        self.assertGreater(published['n_v2_deg'], 0)
        self.assertEqual(sorted(counts['all_contrasts']),
                         sorted([CONTRAST] + ARTIFACT_CONTRASTS))
        self.assertTrue(counts['all_contrasts'][CONTRAST]['published'])
        for contrast in ARTIFACT_CONTRASTS:
            self.assertFalse(counts['all_contrasts'][contrast]['published'])
        axis = counts['axis']
        self.assertEqual(axis['n_axis_genes'], AXIS_SIZE)
        self.assertLessEqual(axis['n_tested'], AXIS_SIZE)
        self.assertLessEqual(axis['n_v2_deg'], published['n_v2_deg'])
        self.assertEqual(axis['n_up'] + axis['n_down'], axis['n_v2_deg'])
        self.assertGreater(axis['n_v2_deg'], 0)

    def test_method_records_the_live_implementation_and_the_fit(self):
        method = summary()['method']
        self.assertEqual(method['implementation'],
                         'stage2/01_bases/output/api/bases/disease_deg.R')
        self.assertEqual(len(method['implementation_sha256']), 64)
        self.assertIn('limma-voom', method['model'])
        self.assertIn('blocking', method['model'])
        self.assertRegex(method['limma_version'], r'^\d+\.\d+\.\d+$')
        self.assertRegex(method['r_version'], r'^R version ')
        fit = summary()['fit_summary']
        self.assertEqual(fit['block'], 'donor')
        self.assertEqual(fit['n_blocks'], 14)
        self.assertGreater(fit['n_genes_dropped'], 0)
        self.assertEqual(fit['n_genes_matrix'], 28380)
        self.assertEqual(fit['n_genes_tested'], fit['n_genes_input'])
        self.assertLess(abs(fit['consensus_correlation']), 1.0)

    def test_declared_paths_and_digests(self):
        record = summary()
        for key in ('deg', 'zone_artifact', 'matrix'):
            self.assertFalse(record[key]['path'].startswith('/'))
            self.assertTrue(record[key]['path'].startswith('local/'))
            self.assertEqual(len(record[key]['sha256']), 64)
        self.assertFalse(record['pseudobulk']['path'].startswith('/'))
        self.assertEqual(len(record['pseudobulk']['sha256']), 64)
        self.assertEqual(record['pseudobulk']['n_genes'], 28380)
        self.assertEqual(record['pseudobulk']['n_profiles'], 25)
        self.assertNotRegex(json.dumps(record), r'"/(home|mnt|tmp)/')

    def test_pseudobulk_summary_states_the_gate_and_the_zones(self):
        record = pseudobulk_summary()
        self.assertEqual(record['schema'], 'stage2-item-cm-pseudobulk-v1')
        self.assertEqual(record['gate'], 'cell_type_original == "Cardiomyocyte"')
        self.assertEqual(record['zone_column'], 'major_labl')
        self.assertEqual(record['n_cells'], 191795)
        self.assertEqual(record['n_cm_cells'], 64510)
        self.assertEqual(record['n_features'], 28380)
        self.assertEqual(record['n_donors'], 20)
        self.assertEqual(record['cm_cells_per_zone'],
                         {'CTRL': 18310, 'IZ': 4866, 'BZ': 10594, 'RZ': 23141, 'FZ': 7599})
        self.assertEqual(record['n_profiles'], 25)
        self.assertEqual(len({entry['profile'] for entry in record['profiles']}), 25)
        self.assertIn('before any testing', record['counts'])
        self.assertEqual(len(record['pseudobulk_sha256']), 64)
        self.assertEqual(len(record['profiles_sha256']), 64)


class BuildSpecTest(unittest.TestCase):
    """The build declarations, fetch pin and provenance record agree."""

    def test_deg_spec_matches_the_summary(self):
        module = load_item_module('build_deg')
        spec = module.SPEC
        record = summary()
        self.assertEqual(spec['item'], ITEM)
        self.assertEqual(spec['category'], 'ICM')
        self.assertEqual(spec['tier'], 1)
        self.assertEqual(spec['block'], 'donor')
        self.assertEqual(spec['normalisation'], record['design']['normalisation'])
        self.assertEqual(spec['formula'], record['design']['formula'])
        self.assertEqual(sorted(contrast['label'] for contrast in spec['contrasts']),
                         sorted([CONTRAST] + ARTIFACT_CONTRASTS))
        self.assertEqual(module.MIN_CM_NUCLEI, 500)
        self.assertEqual(module.PUBLISHED_CONTRAST, CONTRAST)
        self.assertEqual(module.EXPECTED_ZONE_PROFILES,
                         summary()['profile_floor']['zone_profile_counts'])
        # The declared design object is exactly what the shared worker received.
        self.assertEqual(sorted(contrast['label'] for contrast in module.DESIGN['contrasts']),
                         sorted(contrast['label'] for contrast in record['design']['contrasts']))
        self.assertEqual(module.DESIGN['block'], record['design']['block'])
        self.assertEqual(module.DESIGN['filter'], record['design']['filter'])

    def test_pseudobulk_builder_constants_match_the_deposit(self):
        module = load_item_module('build_pseudobulk')
        self.assertEqual(module.CELL_TYPE_COLUMN, 'cell_type_original')
        self.assertEqual(module.CELL_TYPE_LABEL, 'Cardiomyocyte')
        self.assertEqual(module.ZONE_COLUMN, 'major_labl')
        self.assertEqual(module.ZONES, ('CTRL', 'IZ', 'BZ', 'RZ', 'FZ'))
        self.assertEqual(module.EXPECTED_NUCLEI, pseudobulk_summary()['n_cells'])
        self.assertEqual(module.EXPECTED_CM_NUCLEI, pseudobulk_summary()['n_cm_cells'])
        self.assertEqual(module.EXPECTED_FEATURES, pseudobulk_summary()['n_features'])
        self.assertEqual(module.EXPECTED_DONORS, pseudobulk_summary()['n_donors'])

    def test_fetch_pins_the_deposited_bytes(self):
        script = (ITEM_DIR / 'script/fetch.sh').read_text()
        self.assertIn('8191c283-0816-424b-9b61-c3e1d6258a77', script)
        self.assertIn('1c739a3e-c3f5-49d5-98e0-73975e751201', script)
        self.assertIn('4780eb7d-4c7c-4630-8705-471f0de8e27a', script)
        self.assertRegex(script, rf'BYTES={SOURCE_BYTES}\b')
        self.assertIn(f'NAME="4780eb7d-4c7c-4630-8705-471f0de8e27a.h5ad"', script)

    def test_verify_is_offline_and_the_manifest_pins_the_object(self):
        verify = (ITEM_DIR / 'script/verify.sh').read_text()
        self.assertIn('manifest.sha256', verify)
        self.assertNotIn('curl', verify)
        lines = [line for line in (ITEM_DIR / 'manifest.sha256').read_text().splitlines()
                 if line.strip()]
        self.assertEqual(len(lines), 1)
        digest, name = lines[0].split()
        self.assertEqual(len(digest), 64)
        self.assertEqual(name, '4780eb7d-4c7c-4630-8705-471f0de8e27a.h5ad')

    def test_author_log_pins_the_source_tree(self):
        log = author_log()
        self.assertEqual(log['category'], 'disease_data')
        self.assertEqual(log['manifest'], 'manifest.sha256')
        self.assertEqual(log['verify_result'], 'pass')
        (source,) = log['sources']
        self.assertIn('cellxgene:', source['id'])
        self.assertIn('8191c283-0816-424b-9b61-c3e1d6258a77', source['id'])
        self.assertIn('1c739a3e-c3f5-49d5-98e0-73975e751201', source['id'])
        self.assertIn('4780eb7d-4c7c-4630-8705-471f0de8e27a', source['version'])
        self.assertEqual(source['bytes'], SOURCE_BYTES)
        self.assertEqual(source['files'], 1)
        manifest = (ITEM_DIR / 'manifest.sha256').read_bytes()
        self.assertEqual(source['tree_sha256'], hashlib.sha256(manifest).hexdigest())
        declared = {entry['path'] for entry in log['local_additions']}
        self.assertIn('.gitignore', declared)
        record = summary()
        for path in (record['deg']['path'], record['zone_artifact']['path'],
                     record['matrix']['path'], record['pseudobulk']['path'],
                     record['samples']['path'],
                     'local/deg/kuppe-mi_deg.json',
                     'local/pseudobulk/kuppe-mi_cm_pseudobulk.json'):
            self.assertIn(path, declared, path)
        self.assertNotRegex(json.dumps(log), r'"/(home|mnt|tmp)/')

    def test_item_gitignore_covers_the_regenerable_additions(self):
        patterns = [line.strip() for line in
                    (ITEM_DIR / '.gitignore').read_text().splitlines()
                    if line.strip() and not line.startswith('#')]
        for pattern in ('local/pseudobulk/*.csv.gz', 'local/pseudobulk/*.csv',
                        'local/deg/*.csv', 'local/deg/*.tsv.gz'):
            self.assertIn(pattern, patterns)
        # The tracked summaries must not be ignored.
        for pattern in patterns:
            self.assertNotIn('*.json', pattern)

    def test_item_layout(self):
        for name in ('script/fetch.sh', 'script/verify.sh', 'script/hash_source.py',
                     'script/build_pseudobulk.py', 'script/build_deg.py',
                     'author_log.json', 'manifest.sha256', '.gitignore'):
            self.assertTrue((ITEM_DIR / name).is_file(), name)
        for name in ('local/deg/kuppe-mi_deg.json',
                     'local/pseudobulk/kuppe-mi_cm_pseudobulk.json'):
            self.assertTrue((ITEM_DIR / name).is_file(), name)


if __name__ == '__main__':
    unittest.main()
