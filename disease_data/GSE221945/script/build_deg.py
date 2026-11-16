"""Preserve complete author DESeq2 test, reversing control-minus-case direction."""
import argparse,csv,gzip,hashlib,json,math
from pathlib import Path
ITEM=Path(__file__).resolve().parents[1]

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def build(source,out):
    out.mkdir(parents=True,exist_ok=True);target=out/'GSE221945_deg.csv';seen=set()
    fields=['source_gene_id','gene_symbol','contrast','log2fc','stat','p_value','padj','base_mean','n_case','n_control']
    with gzip.open(source,'rt') as f,target.open('w',newline='') as g:
        reader=csv.DictReader(f)
        if not {'','Feature_gene_name','baseMean','log2FoldChange','stat','pvalue','padj'}<=set(reader.fieldnames):raise ValueError('unexpected author test columns')
        writer=csv.DictWriter(g,fieldnames=fields,lineterminator='\n');writer.writeheader()
        for r in reader:
            if not r[''] or r[''] in seen:raise ValueError('missing/duplicate source gene')
            seen.add(r['']);effect=float(r['log2FoldChange']);stat=float(r['stat'])
            if not math.isfinite(effect) or not math.isfinite(stat):raise ValueError('nonfinite author effect/stat')
            for k in ['padj','pvalue']:
                if r[k] not in ['NA',''] and not 0<=float(r[k])<=1:raise ValueError('invalid author p value')
            writer.writerow(dict(source_gene_id=r[''],gene_symbol=r['Feature_gene_name'],contrast='TBX5_G145R-corrected',log2fc=-effect,stat=-stat,p_value='' if r['pvalue']=='NA' else r['pvalue'],padj='' if r['padj']=='NA' else r['padj'],base_mean=r['baseMean'],n_case=4,n_control=4))
    summary={'series':'GSE221945','species':'Homo sapiens','category':'genetic_cardiomyopathy','method':'Complete deposited DESeq2 Wald test; log2FoldChange and stat sign reversed from isogenic-minus-Brugada to Brugada-minus-corrected; p values unchanged','contrasts':{'TBX5_G145R-corrected':{'n_case':4,'n_control':4}},'deposition':{'path':'source/'+source.name,'sha256':sha(source),'bytes':source.stat().st_size},'deg':{'path':'local/deg/'+target.name,'rows':len(seen),'n_source_genes':len(seen),'sha256':sha(target),'bytes':target.stat().st_size},'limitations':['Four mutant and four corrected RNA libraries per exact table columns; not eight independent patients.','Series overall_design says3 corrected lines but sample matrix and differential table contain4 libraries; exact measured libraries retained.','Author min5 tested-gene filter preserved; no local count fit or population-control/knockout contrast added.']}
    (out/'GSE221945_deg.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=ITEM/'source/GSE221945_Isogenic_vs_Brugada_min5_DESeq2.csv.gz');p.add_argument('--out',type=Path,default=ITEM/'local/deg');a=p.parse_args();build(a.source,a.out)
