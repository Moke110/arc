"""Build/replay bases-owned raw tenv counts without restoring deleted steps."""
import argparse
import json
from pathlib import Path
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.native import build_tenv,verify_tenv
root=next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--data-root',type=Path,default=root/'pert_data')
p.add_argument('--output',type=Path,default=root/'stage2/01_bases/output/tenv-raw-counts')
p.add_argument('--verify',action='store_true');p.add_argument('--record',type=Path)
a=p.parse_args()
result=(verify_tenv if a.verify else build_tenv)(a.data_root/'GSE262419/source',a.data_root/'Annotation/source/hgnc_complete_set.txt',a.data_root/'Annotation/local/cm_compound_annotation.csv',a.output)
s=json.dumps(result,indent=2,sort_keys=True)+'\n'
if a.record:a.record.parent.mkdir(parents=True,exist_ok=True);a.record.write_text(s)
print(s)
