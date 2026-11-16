"""Public adapter CLI: retain non-DEGs and preserve deposited statistics."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import openpyxl

class AdapterTest(unittest.TestCase):
    def test_full_table_preserves_missing_adjusted_p_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'source.xlsx'
            book = openpyxl.Workbook()
            sheet = book.active
            sheet.title = 'Cardiomyocytes'
            sheet.append(['Gene', 'log2FoldChange', 'padj'])
            sheet.append(['NPPA', 2.5, 0.001])
            sheet.append(['MYH6', -0.3, None])
            book.save(source)
            subprocess.run([sys.executable, str(Path(__file__).with_name('build_deg.py')),
                            '--source', str(source), '--out', str(root / 'out')], check=True)
            with (root / 'out/GSE183852_deg.csv').open() as f:
                rows = list(csv.DictReader(f))
            self.assertEqual([r['gene_symbol'] for r in rows], ['NPPA', 'MYH6'])
            self.assertEqual(rows[0]['log2fc'], '2.5')
            self.assertEqual(rows[1]['padj'], '')
            self.assertEqual(rows[0]['n_case'], '13')
            self.assertEqual(rows[0]['n_control'], '25')

if __name__ == '__main__':
    unittest.main()
