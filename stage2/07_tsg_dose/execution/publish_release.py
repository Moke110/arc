"""Retired model-publication entrypoint; round1 publishes report evidence only.

Frozen implementation is retained in historical/publish_release_v1.py for the
original 147-fit acceptance. It is excluded from active downstream contracts.
"""
import importlib.util
from pathlib import Path

STEP=Path(__file__).resolve().parents[1]
_spec=importlib.util.spec_from_file_location('_historical_dose_publication',STEP/'execution/historical/publish_release_v1.py')
_archive=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_archive)
_archive.STEP=STEP
# Integrity helpers are shared by report-only acceptance; none writes weights.
for _name in ['repository_root','read','write','digest','file_snapshot',
              'verify_file_snapshot','capture_preservation','csv_rows','csv_write']:
    if hasattr(_archive,_name):globals()[_name]=getattr(_archive,_name)

def publish(*args,**kwargs):
    raise ValueError('deployment weight publication is retired; use round1 report acceptance')

assemble=install=downstream=publish

def main():
    raise SystemExit('deployment weight publication and downstream weight replay are retired; use round1 reports')

if __name__=='__main__':main()
