"""Verify one surviving indexed training payload from bases-owned raw counts."""
import argparse
import json
from pathlib import Path
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.corpora import verify_training_source
from bases import raw
root=next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source',choices=['tenv','tahoe'],required=True)
raw.add_binding_arguments(p)
p.add_argument('--training',type=Path,default=root/'stage2/01_bases/output/training-corpora')
p.add_argument('--drugs',type=Path,default=root/'pert_data/Tahoe-100M/source/metadata/drug_metadata.parquet')
p.add_argument('--record',type=Path)
a=p.parse_args();bindings=raw.cli_bindings(a,root)
result=verify_training_source(bindings,a.training,a.source,bindings["genes.csv"].parent,a.drugs)
s=json.dumps(result,indent=2,sort_keys=True)+'\n'
if a.record:a.record.parent.mkdir(parents=True,exist_ok=True);a.record.write_text(s)
print(s)
