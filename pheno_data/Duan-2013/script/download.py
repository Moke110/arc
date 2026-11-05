#!/usr/bin/env python3
"""Acquire or verify Duan 2013 original MTT article and Figure2."""
from pathlib import Path
import argparse
import sys
ITEM=Path(__file__).resolve().parent.parent
ROOT=next(p for p in ITEM.parents if (p/'.git').exists())
sys.path.insert(0,str(ROOT/'script'))
from phenotype_archive import acquire,verify
ASSETS={'PMC4002788.html':'https://pmc.ncbi.nlm.nih.gov/articles/PMC4002788/',
        'aps2012185f2.jpg':'https://cdn.ncbi.nlm.nih.gov/pmc/blobs/55ec/4002788/62dd84a93b5a/aps2012185f2.jpg'}
PROVENANCE={'upstream_identifier':'10.1038/aps.2012.185','version':'Published2013;PMC4002788 original article; source hashes frozen',
    'license':{'status':'Copyright2013 CPS and SIMM; no open reuse license established','evidence':'PMC4002788.html copyright'},
    'validation':{'status':'validated_original_formats_and_sha256','expected_html_markers':{'PMC4002788.html':['10.1038/aps.2012.185','23396376','Safflor yellow A protects neonatal rat cardiomyocytes']}}}
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    if args.verify: verify(ITEM,'duan_2013','PMC4002788',ASSETS)
    else: acquire(ITEM,'duan_2013','PMC4002788',ASSETS,PROVENANCE)
