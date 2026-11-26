#!/usr/bin/env python3
"""One donor-pooled, run-pattern-adjusted fit of the author CM count source."""
import hashlib,json,sys
from pathlib import Path
import pandas as pd
ITEM=Path(__file__).resolve().parents[1]
REPO=ITEM
while not (REPO/'.git').exists():REPO=REPO.parent
sys.path.insert(0,str(REPO/'stage2_CarPS-VAE/util/output/api'))
from stage2_util.derivation import disease_deg
CONTRASTS=[dict(label='DCM-NF',numerator='DCM',denominator='NF'),dict(label='ICM-NF',numerator='ICM',denominator='NF')]
DESIGN=dict(kind='counts',formula='~0 + group + batch',factor='group',contrasts=CONTRASTS)
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def main():
 pb=ITEM/'local/pseudobulk/gse298023_cm_pseudobulk.tsv';donors=ITEM/'local/pseudobulk/gse298023_donors.csv';out=ITEM/'local/deg/gse298023_deg.csv';out.parent.mkdir(parents=True,exist_ok=True)
 samples=pd.read_csv(donors)
 if samples['donor'].duplicated().any() or samples.groupby('donor')['group'].nunique().max()!=1:raise ValueError('one profile per author heart is required')
 fit=disease_deg.run_deg_worker(pb,donors,DESIGN,out)
 table=pd.read_csv(out)
 source=ITEM/'source/GSE298023_HH_dataset_seurat_object_2026.rds'
 record=dict(schema='stage2-item-deg-summary-v1',item='GSE298023',category='nHF / ICM',gate='author ident == Cardiomyocytes',
  source=dict(path='source/'+source.name,sha256=digest(source),bytes=source.stat().st_size),
  pseudobulk=dict(path='local/pseudobulk/'+pb.name,sha256=digest(pb),dtype='integer counts summed per author heart_number over CM cells and repeated scar/nonscar regions',cells_per_donor={r.sample:int(r.cells_cm) for r in samples.itertuples()}),
  donors=dict(path='local/pseudobulk/'+donors.name,sha256=digest(donors)),design=DESIGN,
  method=dict(model='limma-voom; moderated t; BH adjusted p values',implementation='stage2_CarPS-VAE/util/output/api/stage2_util/derivation/disease_deg.R',implementation_sha256=digest(disease_deg.DEG_WORKER),adapter_sha256=digest(__file__),pseudobulk_adapter_sha256=digest(ITEM/'script/build_pseudobulk.R'),r_version=fit['r_version'],limma_version=fit['limma_version']),
  contrasts={c['label']:dict(n_case=int((samples.group==c['numerator']).sum()),n_control=int((samples.group==c['denominator']).sum())) for c in CONTRASTS},
  fit_summary=fit,
  limitations=['Author heart_number defines the independent donor, including pooled male/female NF assignments already in the deposited object; no new sex-demultiplexing or marker gate.',
  'Sum repeated scar/nonscar profiles before fitting; no region-level pseudoreplication. Region composition remains a recorded ICM limitation.',
  'NF is run2 only. Disease cohorts also include run2, so disease is estimable with the run-pattern nuisance factor; mixed-run Heart03 gets its own run_1_2 pattern. No unsupported clinical covariates or genotype exclusion.'],
  deg=dict(path='local/deg/'+out.name,sha256=digest(out),bytes=out.stat().st_size,rows=len(table),n_source_genes=int(table.source_gene_id.nunique())))
 (out.parent/'gse298023_deg.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
 print(json.dumps(dict(result='PASS',contrasts=record['contrasts'],deg=record['deg']),indent=2))
if __name__=='__main__':main()
