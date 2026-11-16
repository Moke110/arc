"""Drop deposited symbol annotation, retaining actual numeric Entrez source keys."""
import argparse
import csv
import gzip
import io
from pathlib import Path
ITEM=Path(__file__).resolve().parents[1]

def prepare(source,out):
    out.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(source,'rt') as f, out.open('wb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0,filename='') as zipped:
            with io.TextIOWrapper(zipped,newline='') as g:
                reader=csv.reader(f,delimiter='\t');header=next(reader)
                if header[:2]!=['ensembl_gene_id','hgnc_symbol']:
                    raise ValueError('unexpected deposit columns')
                writer=csv.writer(g,delimiter='\t',lineterminator='\n');writer.writerow(['source_gene_id']+header[2:]);seen=set()
                for row in reader:
                    if not row[0].isdigit() or row[0] in seen or len(row)!=len(header):
                        raise ValueError('missing/duplicate numeric Entrez key or malformed row')
                    seen.add(row[0])
                    if any(not x.isdigit() for x in row[2:]):raise ValueError('noninteger count')
                    writer.writerow([row[0]]+row[2:])
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=ITEM/'source/GSE329393_raw_counts.tsv.gz');p.add_argument('--out',type=Path,default=ITEM/'local/deg/GSE329393_counts.tsv.gz');a=p.parse_args();prepare(a.source,a.out)
