#!/usr/bin/env python3
"""Acquire or verify Liu 2021 icariin article and XML-confirmed Figure2."""
from pathlib import Path
import argparse
import sys
ITEM=Path(__file__).resolve().parent.parent
ROOT=next(p for p in ITEM.parents if (p/'.git').exists())
sys.path.insert(0,str(ROOT/'script'))
from phenotype_archive import acquire,verify
BASE='https://pmc-oa-opendata.s3.amazonaws.com/PMC8564343.1/'
ASSETS={'article.xml':BASE+'PMC8564343.1.xml','FEB4-11-2966-g004.jpg':BASE+'FEB4-11-2966-g004.jpg'}
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    if args.verify: verify(ITEM,'liu_2021','PMC8564343',ASSETS)
    else: acquire(ITEM,'liu_2021','PMC8564343',ASSETS,{'upstream_identifier':'10.1002/2211-5463.13276','version':'PMC8564343.1','license':{'status':'CC BY 4.0','evidence':'article.xml permissions'},'validation':{'status':'validated_original_formats_and_sha256'}})
