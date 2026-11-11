"""Build/replay bounded shared expression entirely within bases."""
import argparse,json
from pathlib import Path
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases import raw
root=next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())
p=argparse.ArgumentParser(description=__doc__)
raw.add_binding_arguments(p)
p.add_argument('--output',type=Path,default=root/'stage2/01_bases/output/shared-expression')
p.add_argument('--sources',nargs='+',choices=['tenv','tahoe'],default=['tahoe','tenv'])
p.add_argument('--verify',action='store_true');p.add_argument('--record',type=Path)
a=p.parse_args();bindings=raw.cli_bindings(a,root)
result=raw.verify(bindings,a.output) if a.verify else raw.build(bindings,a.output,sources=tuple(a.sources))
s=json.dumps(result,indent=2,sort_keys=True)+'\n'
if a.record:a.record.parent.mkdir(parents=True,exist_ok=True);a.record.write_text(s)
print(s)
