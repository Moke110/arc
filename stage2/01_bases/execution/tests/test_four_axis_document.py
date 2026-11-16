"""Four-axis reading aid at the public document build/read seams."""
from pathlib import Path
import sys
import unittest

EXECUTION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXECUTION))
import make_gene_axis_docx as builder


class FourAxisDocumentTests(unittest.TestCase):
    def test_default_destination_and_rule_describe_successor(self):
        self.assertEqual(builder.default_paths()['out'], EXECUTION.parent / 'gene_axis_fill_rule.docx')
        text = '\n'.join(builder.RULE_STEPS)
        self.assertIn('source_tier_highest', text)
        self.assertIn('1,000', text)
        self.assertIn('3,000', text)
        self.assertIn('4,000', text)
        self.assertNotIn('达到即停', text)
        self.assertNotIn('DEG2000', text)
        self.assertNotIn('f2000', text)

    def test_build_is_portrait_with_separate_category_tables_and_dynamic_eligibility(self):
        import json
        import tempfile
        import types
        import pandas as pd
        from docx import Document
        from docx.shared import Cm
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            catalogs, sources, records = {}, {}, []
            coverage = []
            for category in builder.gene_axis.RELEASE_CATEGORIES:
                label = builder.gene_axis.CATEGORY_LABELS[category]
                count = 25 if category == 'af' else 120
                dataset_id = category + '_fixture'
                sources[dataset_id] = {'accession': 'GSE000001', 'accession_note': None, 'technology': 'bulk RNA-seq'}
                dataset = dict(dataset_id=dataset_id, n_case=4, n_control=4, n_deg=count, arm_or_class='case versus control')
                catalog = dict(category=label, datasets=[dataset], tier_sources=[{'tier': 1}], n_deg=count, deg1000=count)
                path = root / (category + '.json'); path.write_text(json.dumps(catalog)); catalogs[category] = path
                cut = category + '.deg1000.csv'
                pd.DataFrame({'gene_index': range(count)}).to_csv(root / cut, index=False)
                records.append(dict(category=label, declared_source_count=1, n_deg=count, deg1000=count,
                                    deg1000_file=cut, eligible=count >= 100))
                for size in builder.gene_axis.AXIS_SIZES:
                    excluded = count < 100
                    coverage.append([size, size, category, label, None if excluded else count,
                                     None if excluded else count, None if excluded else 100.0, 'union below 100' if excluded else ''])
            sources['tenv_hvg'] = {'accession': None, 'accession_note': 'GSE262419 corpus', 'technology': 'TempO-Seq'}
            metadata = types.SimpleNamespace(SOURCES=sources, TENV_HVG='tenv_hvg', validate=lambda ids: None)
            for size in builder.gene_axis.AXIS_SIZES:
                coverage.append([size, size, 'tenv_hvg', 'tenv-HVG', size, size, 100.0, ''])
            coverage_path = root / 'coverage.csv'
            pd.DataFrame(coverage, columns=builder.COVERAGE_COLUMNS).to_csv(coverage_path, index=False)
            fusion = root / 'fusion.json'; fusion.write_text(json.dumps({'diseases': records}))
            axis = root / 'genes.csv'; pd.DataFrame({'gene_index': range(120), 'symbol': ['G'+str(i) for i in range(120)]}).to_csv(axis,index=False)
            manifest = root / 'manifest.json'; manifest.write_text(json.dumps({'release': {'algorithm': builder.gene_axis.ALGORITHM,
                'fallback_second_half_slots': {'2000': 0, '3000': 0, '4000': 500, '5000': 1000}}}))
            figure = EXECUTION.parent / 'figure/gene_axis_rank_sum_histogram.png'
            output = builder.build(coverage_path, figure, root / 'rule.docx', catalogs, metadata,
                                   fusion, axis, manifest)
            first_bytes = output.read_bytes()
            builder.build(coverage_path, figure, root / 'replay.docx', catalogs, metadata,
                          fusion, axis, manifest)
            self.assertEqual(first_bytes, (root / 'replay.docx').read_bytes())
            doc = Document(output)
            self.assertAlmostEqual(doc.sections[0].page_width.cm, 21, places=2)
            self.assertAlmostEqual(doc.sections[0].page_height.cm, 29.7, places=2)
            source_tables = [t for t in doc.tables if t.rows[0].cells[0].text == 'category']
            self.assertEqual(len(source_tables), 6)
            self.assertEqual([t.rows[1].cells[0].text for t in source_tables[:-1]],
                             [builder.gene_axis.CATEGORY_LABELS[c] for c in builder.gene_axis.RELEASE_CATEGORIES])
            text = '\n'.join(p.text for p in doc.paragraphs)
            self.assertIn('25', '\n'.join(c.text for t in doc.tables for r in t.rows for c in r.cells))
            self.assertIn('N=4000 用 500', text)
            self.assertNotIn('AF 的一个', text)

    def test_four_axis_histogram_remains_readable_in_portrait_document(self):
        import importlib.util
        import tempfile
        from PIL import Image
        path = EXECUTION.parent / 'figure/script/render_gene_axis_rank_sum_histogram.py'
        spec = importlib.util.spec_from_file_location('histogram224', path)
        renderer = importlib.util.module_from_spec(spec); spec.loader.exec_module(renderer)
        panels = [dict(panel=str(n), axis=n, panel_title='N = '+str(n)) for n in (2000,3000,4000,5000)]
        histograms = {str(n): [(0.0, 1.0, 18164, n)] for n in (2000,3000,4000,5000)}
        with tempfile.TemporaryDirectory() as root:
            out = Path(root) / 'four.png'
            renderer.render('DEG1000 f1000', panels, histograms, out)
            with Image.open(out) as image:
                self.assertLess(image.width / image.height, 1.5)
