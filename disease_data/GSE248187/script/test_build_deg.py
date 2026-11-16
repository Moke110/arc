import csv,gzip,subprocess,sys,tempfile,unittest
from pathlib import Path
class AdapterTest(unittest.TestCase):
 def test_preserves_complete_author_statistics(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t);src=p/'input.gz'
   with gzip.open(src,'wt') as f:f.write('geneID\tgene_name\tUC1_cpm\tUC2_cpm\tUC3_cpm\tDSP1_cpm\tDSP2_cpm\tDSP3_cpm\tlogFC\tlogCPM\tLR\tp-val\tFDR-pval\nENSG1\tNPPA\t1\t2\t3\t4\t5\t6\t2.5\t3\t10\t0.001\t0.01\nENSG2\tMYH6\t1\t1\t1\t1\t1\t1\t0.1\t1\t0.1\t0.8\t0.9\n')
   subprocess.run([sys.executable,str(Path(__file__).with_name('build_deg.py')),'--source',str(src),'--out',str(p/'out')],check=True)
   with (p/'out/GSE248187_deg.csv').open() as f:rows=list(csv.DictReader(f))
   self.assertEqual(len(rows),2);self.assertEqual(rows[0]['log2fc'],'2.5');self.assertEqual(rows[0]['padj'],'0.01');self.assertEqual(rows[0]['n_case'],'3')
if __name__=='__main__':unittest.main()
