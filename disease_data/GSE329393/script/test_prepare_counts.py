import gzip,subprocess,sys,tempfile,unittest,csv
from pathlib import Path
class AdapterTest(unittest.TestCase):
 def test_keeps_numeric_entrez_identifiers_and_actual_sample_names(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t);src=p/'input.gz'
   with gzip.open(src,'wt') as f:f.write('ensembl_gene_id\thgnc_symbol\t12706\t14110\n100287102\tNA\t4\t5\n123\tNPPA\t6\t7\n')
   subprocess.run([sys.executable,str(Path(__file__).with_name('prepare_counts.py')),'--source',str(src),'--out',str(p/'counts.tsv.gz')],check=True)
   with gzip.open(p/'counts.tsv.gz','rt') as f:rows=list(csv.reader(f,delimiter='\t'))
   self.assertEqual(rows,[['source_gene_id','12706','14110'],['100287102','4','5'],['123','6','7']])
if __name__=='__main__':unittest.main()
