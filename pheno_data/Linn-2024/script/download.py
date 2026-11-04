#!/usr/bin/env python3
"""Acquire or verify Linn 2024 article and exact cytotoxicity figure."""
from pathlib import Path
import argparse
import sys
ITEM=Path(__file__).resolve().parent.parent
ROOT=next(p for p in ITEM.parents if (p/'.git').exists())
sys.path.insert(0,str(ROOT/'script'))
from phenotype_archive import acquire,verify
BASE='https://pmc-oa-opendata.s3.amazonaws.com/PMC10855687.1/'
ASSETS={'article.xml':BASE+'PMC10855687.1.xml','metadata.json':BASE+'PMC10855687.1.json','ijms-25-01427-g003.jpg':BASE+'ijms-25-01427-g003.jpg'}
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    if args.verify: verify(ITEM,'linn_2024','PMC10855687',ASSETS)
    else: acquire(ITEM,'linn_2024','PMC10855687',ASSETS,{'upstream_identifier':'10.3390/ijms25031427','version':'PMC10855687.1','license':{'status':'CC BY 4.0','evidence':'article.xml permissions; metadata.json'},'validation':{'status':'validated_original_formats_and_sha256'}})
