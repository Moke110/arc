"""Prepare deposited numeric Entrez counts then run the one live DEG worker."""
import csv,gzip,json,sys
from pathlib import Path
from prepare_counts import prepare
ITEM=Path(__file__).resolve().parents[1]
REPO=ITEM
while not (REPO/'.git').exists():REPO=REPO.parent
sys.path.insert(0,str(REPO/'stage2_CarPS-VAE/util/output/api'))
from stage2_util.derivation import bulk_deg
SPEC={'matrix': 'local/deg/GSE329393_counts.tsv.gz', 'sep': '\t', 'series_matrix': 'source/GSE329393_series_matrix.txt.gz', 'id_kind': 'entrez', 'sample_match': 'description', 'group_source': 'genotype', 'group_map': {'6crispr': 'ctrl6', '6RBM': 'mut6', '9crispr': 'ctrl9', '9RBM': 'mut9'}, 'kind': 'counts', 'formula': '~0 + group', 'contrasts': [{'label': 'RBM20_LVNC-corrected', 'numerator': 'mut6', 'denominator': 'ctrl6'}, {'label': 'RBM20_DCM-corrected', 'numerator': 'mut9', 'denominator': 'ctrl9'}], 'contrast_direction': 'positive log2FC = mutant minus matched-family CRISPR-corrected control', 'isogenic_control': {'arm': '5 corrected libraries per family', 'basis': 'GEO genotype6crispr/9crispr vs6RBM/9RBM, patient-specific CRISPR corrected controls'}, 'strata': {'decision': 'one fit four family/genotype groups; two matched-family contrasts collapse as one dataset', 'pooled': False, 'blocking_covariates': [], 'why': 'family and genotype represented jointly, no cross-family comparison or duplicate dataset votes'}}
if __name__=='__main__':
 prepare(ITEM/'source/GSE329393_raw_counts.tsv.gz',ITEM/SPEC['matrix'])
 result=bulk_deg.derive_item_deg(ITEM,SPEC)
 # The deposited symbol column is all NA. Resolve numeric Entrez IDs only
 # when exactly one approved HGNC record exists; never pick an ambiguous symbol.
 annotation=REPO/'pert_data/Annotation/source/hgnc_complete_set.txt'
 with annotation.open() as f:
  candidates={}
  for r in csv.DictReader(f,delimiter='\t'):
   if r['status']=='Approved' and r['entrez_id']:
    candidates.setdefault(r['entrez_id'],set()).add(r['symbol'])
 symbols={k:next(iter(v)) for k,v in candidates.items() if len(v)==1}
 target=ITEM/'local/deg/GSE329393_deg.csv'
 with target.open() as f:
  reader=csv.DictReader(f);fields=reader.fieldnames;rows=list(reader)
 for r in rows:r['gene_symbol']=symbols.get(r['source_gene_id'],'')
 with target.open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');writer.writeheader();writer.writerows(rows)
 summary=ITEM/'local/deg/GSE329393_deg.json';record=json.loads(summary.read_text())
 record['deg'].update(sha256=bulk_deg.sha256(target),bytes=target.stat().st_size)
 record['method']['item_adapter_sha256']=bulk_deg.sha256(Path(__file__))
 record['symbol_annotation']={'path':'pert_data/Annotation/source/hgnc_complete_set.txt','sha256':bulk_deg.sha256(annotation),'binding':'exact numeric Entrez ID to unique approved HGNC symbol; missing/ambiguous empty','n_source_ids_unmapped_or_ambiguous':len({r['source_gene_id'] for r in rows if not r['gene_symbol']})}
 summary.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
 result['deg_sha256']=record['deg']['sha256']
 print(json.dumps(result,indent=2))
