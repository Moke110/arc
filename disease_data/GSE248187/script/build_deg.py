"""Preserve the complete deposited human DSP EHT edgeR test, without refitting."""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
ITEM = Path(__file__).resolve().parents[1]

def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def build(source, out):
    out.mkdir(parents=True, exist_ok=True)
    target = out / 'GSE248187_deg.csv'
    fields = ['source_gene_id','gene_symbol','contrast','log2fc','stat','p_value','padj','base_mean','n_case','n_control']
    seen = set()
    with gzip.open(source,'rt') as f, target.open('w',newline='') as g:
        reader = csv.DictReader(f, delimiter='\t')
        required = {'geneID','gene_name','logFC','logCPM','LR','p-val','FDR-pval'}
        if not required <= set(reader.fieldnames):
            raise ValueError('missing deposited human EHT test columns')
        writer = csv.DictWriter(g,fieldnames=fields,lineterminator='\n');writer.writeheader()
        for r in reader:
            if not r['geneID'] or r['geneID'] in seen:
                raise ValueError('missing/duplicate source gene')
            seen.add(r['geneID'])
            for k in ['logFC','logCPM','LR','p-val','FDR-pval']:
                if not math.isfinite(float(r[k])):
                    raise ValueError('nonfinite author statistic')
            if not 0 <= float(r['FDR-pval']) <= 1 or not 0 <= float(r['p-val']) <= 1:
                raise ValueError('invalid author p value')
            writer.writerow(dict(source_gene_id=r['geneID'],gene_symbol=r['gene_name'],contrast='DSP_KO-isogenic',log2fc=r['logFC'],stat=r['LR'],p_value=r['p-val'],padj=r['FDR-pval'],base_mean='',n_case=3,n_control=3))
    summary = {'series':'GSE248187','category':'genetic_cardiomyopathy','species':'Homo sapiens','method':'Complete author edgeR human EHT statistics preserved; no local count refit','contrasts':{'DSP_KO-isogenic':{'n_case':3,'n_control':3}},'deposition':{'path':'source/'+source.name,'sha256':sha(source),'bytes':source.stat().st_size},'deg':{'path':'local/deg/'+target.name,'rows':len(seen),'n_source_genes':len(seen),'sha256':sha(target),'bytes':target.stat().st_size},'limitations':['Three EHT RNA libraries per arm, not three independent patients.','Human isogenic DSP knockout EHT includes supporting fibroblasts.','Author logCPM is not invented as a linear base_mean; that field is empty.','Rat NRVM DSP/JUP/PKP2 tables are excluded.']}
    (out/'GSE248187_deg.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=ITEM/'source/GSE248187_DSP_EHTs_Differential_Expression_Table.txt.gz');p.add_argument('--out',type=Path,default=ITEM/'local/deg');a=p.parse_args();build(a.source,a.out)
