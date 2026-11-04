from pathlib import Path
import sys
ITEM = Path(__file__).resolve().parent.parent
ROOT = next(p for p in ITEM.parents if (p/'.git').exists())
sys.path.insert(0, str(ROOT/'script'))
from phenotype_archive import acquire, verify
import argparse
p = ITEM
base = 'https://pmc-oa-opendata.s3.amazonaws.com/PMC8371130.1/'
assets = {name: base+name for name in ['PMC8371130.1.xml','PMC8371130.1.json','41598_2021_96033_MOESM1_ESM.pdf','41598_2021_96033_MOESM2_ESM.xlsx','41598_2021_96033_MOESM4_ESM.xlsx']}
parser = argparse.ArgumentParser(description='Acquire or verify Heliste original drug-screen assets')
parser.add_argument('--verify', action='store_true')
args = parser.parse_args()
if args.verify:
    verify(p, 'heliste_2021', 'PMC8371130', assets)
else:
    acquire(p, 'heliste_2021', 'PMC8371130', assets, {'upstream_identifier':'10.1038/s41598-021-96033-z','version':'PMC8371130.1; supplementary datasets S1/S2/S3','license':{'status':'CC BY 4.0, third-party exceptions as article specifies','evidence':'PMC metadata license_code=CC BY and article license text'},'validation':{'status':'formats_and_hashes_validated','scope':'Dataset S1 drug CTG; S2 shRNA is not included; S3 selected hit names','ancillary_not_archived':{'asset':'41598_2021_96033_MOESM3_ESM.xlsx','meaning':'Dataset S2 shRNA, not drug measurements','reason':'two bounded official refetch attempts timed out at 30 seconds; not required for CTG drug extraction'}}})
