import csv,gzip,subprocess,sys,tempfile,unittest
from pathlib import Path
class AdapterTest(unittest.TestCase):
 def test_reverses_isogenic_minus_case_effect_and_stat_only(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t);src=p/'input.gz'
   with gzip.open(src,'wt') as f:f.write(',Feature_gene_name,baseMean,log2FoldChange,stat,pvalue,padj\nENSG1,NPPA,10,-2.5,-3,0.001,0.01\nENSG2,MYH6,4,0.2,0.1,0.7,NA\n')
   subprocess.run([sys.executable,str(Path(__file__).with_name('build_deg.py')),'--source',str(src),'--out',str(p/'out')],check=True)
   with (p/'out/GSE221945_deg.csv').open() as f:rows=list(csv.DictReader(f))
   self.assertEqual(len(rows),2);self.assertEqual(rows[0]['log2fc'],'2.5');self.assertEqual(rows[0]['stat'],'3.0');self.assertEqual(rows[0]['padj'],'0.01');self.assertEqual(rows[1]['padj'],'');self.assertEqual(rows[0]['n_case'],'4')
if __name__=='__main__':unittest.main()
