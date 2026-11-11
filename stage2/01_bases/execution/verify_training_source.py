"""Verify one surviving indexed training payload from bases-owned raw counts."""
import argparse
import json
from pathlib import Path
from bases.corpora import verify_training_source
root=next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source',choices=['tenv','tahoe'],required=True)
p.add_argument('--native',type=Path,default=root/'stage2/01_bases/output/native-v1')
p.add_argument('--training',type=Path,default=root/'stage2/01_bases/output/training-v1')
p.add_argument('--tables',type=Path,default=Path(__file__).resolve().parents[1]/'output/native-v1/tables')
p.add_argument('--drugs',type=Path,default=root/'pert_data/Tahoe-100M/source/metadata/drug_metadata.parquet')
p.add_argument('--record',type=Path)
a=p.parse_args();result=verify_training_source(a.native,a.training,a.source,a.tables,a.drugs)
s=json.dumps(result,indent=2,sort_keys=True)+'\n'
if a.record:a.record.parent.mkdir(parents=True,exist_ok=True);a.record.write_text(s)
print(s)
