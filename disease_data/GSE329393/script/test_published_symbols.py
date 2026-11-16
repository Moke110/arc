"""Integration check of the public item output against its source annotation."""
import csv,unittest
from pathlib import Path
ITEM=Path(__file__).resolve().parents[1]
class PublishedTest(unittest.TestCase):
 @unittest.skipUnless((ITEM/'local/deg/GSE329393_deg.csv').exists(),'materialized item output required')
 def test_nppa_retains_source_symbol(self):
  with (ITEM/'local/deg/GSE329393_deg.csv').open() as f: rows=[r for r in csv.DictReader(f) if r['source_gene_id']=='4878']
  self.assertEqual(len(rows),2)
  self.assertEqual({r['gene_symbol'] for r in rows},{'NPPA'})
if __name__=='__main__':unittest.main()
