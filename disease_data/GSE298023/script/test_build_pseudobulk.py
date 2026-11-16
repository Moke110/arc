import csv,subprocess,tempfile,unittest
from pathlib import Path
R='/home/chang/miniconda3/envs/biomni_e1/bin/Rscript'
class PseudobulkCliTest(unittest.TestCase):
 def test_pools_author_cm_cells_by_heart_and_retains_run_pattern(self):
  with tempfile.TemporaryDirectory() as t:
   p=Path(t);fixture=p/'fixture.R';source=p/'author.rds';out=p/'out'
   fixture.write_text('''suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Matrix))
a<-commandArgs(TRUE)
options(Seurat.object.assay.version='v3')
m<-Matrix(matrix(c(1,2,3,4,100,200,5,6),nrow=2),sparse=TRUE)
rownames(m)<-c('NPPA','MYH6');colnames(m)<-paste0('c',1:4)
meta<-data.frame(ident=c('Cardiomyocytes','Cardiomyocytes','Fibroblasts','Cardiomyocytes'),heart_number=c('Heart_01','Heart_01','Heart_02','Heart_02'),condition=c('DCM','DCM','NF','NF'),region=c('DCM','DCM','NF','NF'),run_number=c('1','2','2','2'),row.names=colnames(m))
x<-CreateSeuratObject(counts=m,meta.data=meta);saveRDS(x,a[1])
''')
   subprocess.run([R,str(fixture),str(source)],check=True,capture_output=True)
   result=subprocess.run([R,str(Path(__file__).with_name('build_pseudobulk.R')),str(source),str(out)],capture_output=True,text=True)
   self.assertEqual(result.returncode,0,result.stderr)
   with (out/'gse298023_cm_pseudobulk.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
   self.assertEqual([(r['Heart_01'],r['Heart_02']) for r in rows],[('4','5'),('6','6')])
   with (out/'gse298023_donors.csv').open() as f:d=list(csv.DictReader(f))
   self.assertEqual([(r['sample'],r['group'],r['cells_cm'],r['batch']) for r in d],[('Heart_01','DCM','2','run_1_2'),('Heart_02','NF','1','run_2')])
if __name__=='__main__':unittest.main()
