"""Active round1 learner/replay entrypoint; deployment training is retired."""
from pathlib import Path
import sys

STEP=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(STEP/'output/api'))
if __name__=='__main__':
    # Check before parsing or touching paths; the frozen learner remains available
    # for audit, while this active entrypoint accepts only round1 scenarios.
    import argparse
    boundary=argparse.ArgumentParser(add_help=False)
    boundary.add_argument('--scenario',choices=['cold-drug','held-dose'],default='cold-drug')
    boundary.add_argument('--candidate',type=Path)
    request,_=boundary.parse_known_args()
    if request.candidate is not None and (request.candidate/'config.json').is_file():
        import json
        scenario=json.loads((request.candidate/'config.json').read_text()).get('scenario')
        if scenario=='final':boundary.error('final origins are historical; active round1 excludes them')
    source=STEP/'execution/historical/transigen_v1.py'
    # Preserve __file__ for the accepted CLI's fresh-process self-replay path.
    exec(compile(source.read_bytes(),str(Path(__file__).resolve()),'exec'),globals())
