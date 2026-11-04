#!/usr/bin/env python3
"""Acquire or verify the pinned original Sirenko archive."""
from pathlib import Path
import argparse
import sys
ITEM = Path(__file__).resolve().parent.parent
ROOT = next(p for p in ITEM.parents if (p/'.git').exists())
sys.path.insert(0, str(ROOT/'script'))
from phenotype_archive import acquire, verify
ASSETS = {
    'publisher-mmc3.xls': 'https://ars.els-cdn.com/content/image/1-s2.0-S0041008X13004171-mmc3.xls',
    'article.xml': 'https://pmc-oa-opendata.s3.amazonaws.com/PMC3900303.1/PMC3900303.1.xml',
    'metadata.json': 'https://pmc-oa-opendata.s3.amazonaws.com/PMC3900303.1/PMC3900303.1.json',
}
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', help='offline verification only')
    args = parser.parse_args()
    if args.verify:
        verify(ITEM, 'sirenko_2013', 'PMC3900303', ASSETS)
    else:
        acquire(ITEM, 'sirenko_2013', 'PMC3900303', ASSETS, {
            'upstream_identifier':'10.1016/j.taap.2013.09.017',
            'version':'PMC3900303.1; publisher Table 3 bytes pinned by manifest',
            'license':{'evidence':'article.xml copyright and metadata.json license_code=TDM; is_pmc_openaccess=false','status':'rights reserved; public accessibility does not establish redistribution permission'},
            'validation':{'status':'validated_original_formats_and_sha256','table':'Supplemental Table 3, 24hr_beat_count'}})
