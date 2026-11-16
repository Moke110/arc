"""Successor human explanation at the existing DOCX public build seam (#218)."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

EXECUTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXECUTION))
import make_gene_axis_docx as builder
import gene_axis_source_metadata as metadata
STEP = Path(os.environ.get('ARC_DOCX_STEP', EXECUTION.parent))

class SuccessorDocumentTests(unittest.TestCase):
    def test_metadata_exactly_covers_declared_successor_sources(self):
        ids = []
        for category in builder.gene_axis.RELEASE_CATEGORIES:
            ids.extend(d['dataset_id'] for d in json.loads((STEP / 'output/disease-deg' / category / 'catalog.json').read_text())['datasets'])
        self.assertEqual(len(ids), 41)
        metadata.validate(ids + [metadata.TENV_HVG])

    def test_af_retention_is_explicitly_excluded_from_coverage(self):
        catalogs = {c: STEP / 'output/disease-deg' / c / 'catalog.json' for c in builder.gene_axis.RELEASE_CATEGORIES}
        fusion_path = STEP / 'output/disease-fusion/catalog.json'
        rows = builder.read_disease_table_rows(catalogs, json.loads(fusion_path.read_text()), fusion_path, builder.read_axis_symbols(STEP / 'output/gene-identity/genes.csv'))
        af = next(r for r in rows if r['disease'] == 'AF')
        self.assertEqual(af['DEG union'], '25')
        self.assertEqual(af['DEG1000'], '25')
        self.assertEqual(af['eligible'], 'no')
        self.assertIn('100', af['note'])
        coverage = builder.read_coverage(STEP / 'output/gene-axis/source_topn_coverage.csv')
        self.assertTrue(coverage.loc[coverage['group'] == 'af', ['numerator', 'denominator', 'rate']].isna().all().all())
