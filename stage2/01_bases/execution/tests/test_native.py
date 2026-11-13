"""Native preparation is tested through its standalone build/replay boundary."""
import sys
from pathlib import Path
import tempfile
import unittest
import pandas as pd
import h5py
from anndata.io import read_elem
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.native import build_tenv, verify_tenv

class NativePreparationTest(unittest.TestCase):
    def test_duplicate_probes_sum_and_raw_wells_survive_without_other_steps(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); source=root/'source'; source.mkdir()
            pd.DataFrame({'Genes':['A_1','A_2','B_3'],'Plate1-A1':[1,2,4],'Plate1-A2':[2,3,5],'Plate1-A3':[0,1,1]}).to_csv(source/'GSE262419_Plate1.csv.gz',index=False)
            pd.DataFrame({'Plate_ID':[1]*3,'Well_ID':['A01','A02','A03'],'Chemical_name':['VEH','ethanol','Blank'],'Chemical_Concentration_uM':[0,1,0],'Chemical_Index':['0','1','2']}).to_csv(source/'GSE262419_hash.csv.gz',index=False)
            hgnc=root/'hgnc.tsv'
            pd.DataFrame({'hgnc_id':['HGNC:1','HGNC:2'],'symbol':['A','B'],'alias_symbol':['',''],'prev_symbol':['',''],'entrez_id':['1','3'],'ensembl_gene_id':['E1','E2']}).to_csv(hgnc,sep='\t',index=False)
            annotation=root/'drugs.csv'
            pd.DataFrame([{'lookup_name':'ethanol','dataset_inchikey':'','parent_inchikey':'LFQSCWFLJHTTHZ-UHFFFAOYSA-N','inchikey':'LFQSCWFLJHTTHZ-UHFFFAOYSA-N','smiles':'CCO'}]).to_csv(annotation,index=False)
            output=root/'native'; result=build_tenv(source,hgnc,annotation,output)
            self.assertEqual(result['shape'],[3,2])
            with h5py.File(output/'tenv.h5ad') as f:
                self.assertEqual(f['X'][:].tolist(),[[3,4],[5,5],[1,1]])
                self.assertEqual(read_elem(f['obs']).is_blank.tolist(),[False,False,True])
                self.assertEqual(read_elem(f['var']).hgnc_id.tolist(),['HGNC:1','HGNC:2'])
                self.assertEqual(set(f['layers']),set())
            self.assertEqual(verify_tenv(source,hgnc,annotation,output)['result'],'PASS')

if __name__ == "__main__":
    unittest.main()
