"""Acquire exact-name PubChem identity evidence; corpus build never calls the network."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import time
import threading
from urllib.parse import quote
from urllib.request import urlopen, Request
import openpyxl
from rdkit import Chem, rdBase

ITEM = Path(__file__).resolve().parent.parent
_REQUEST_LOCK = threading.Lock()
_LAST_REQUEST = 0.0
API = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/'

def normalize(s):
    return re.sub(r'[^a-z0-9]', '', s.lower().replace('hcl', 'hydrochloride').replace('(+)', 'plus').replace('(-)', 'minus'))

def fetch(url):
    global _LAST_REQUEST
    for attempt in range(3):
        try:
            with _REQUEST_LOCK:
                time.sleep(max(0, .23-(time.monotonic()-_LAST_REQUEST)))
                _LAST_REQUEST=time.monotonic()
            with urlopen(Request(url, headers={'User-Agent':'ARC-scientific-identity-curation/1.0'}), timeout=25) as response:
                return json.load(response)
        except Exception as exc:
            if getattr(exc, 'code', None) == 404:
                return {'error': 'not_found'}
            if attempt == 2:
                return {'error': type(exc).__name__+': '+str(exc)}
            time.sleep(2+attempt)

def resolve(label):
    names=supplied_queries(label)
    queries=names[1:]+names[:1] if len(names)>1 else names
    attempts = []
    for query in dict.fromkeys(queries):
        url = API+'name/'+quote(query, safe='')+'/property/IsomericSMILES/JSON'
        props = fetch(url)
        evidence = {'query':query, 'url':url, 'response':props,
                    'retrieved_at':datetime.now(timezone.utc).isoformat()}
        attempts.append(evidence)
        entries = props.get('PropertyTable',{}).get('Properties',[])
        if len(entries) != 1:
            continue
        entry = entries[0]
        synurl = API+'cid/'+str(entry['CID'])+'/synonyms/JSON'
        syn = fetch(synurl)
        evidence.update(synonym_url=synurl, synonym_response=syn)
        synonyms = syn.get('InformationList',{}).get('Information',[{}])[0].get('Synonym',[])
        if not any(normalize(query)==normalize(s) for s in synonyms):
            continue
        smiles = entry.get('SMILES')
        mol = Chem.MolFromSmiles(smiles or '')
        if mol is None or not smiles:
            continue
        candidate = {'source_label':label, 'status':'resolved', 'canonical_smiles':Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True),
                'evidence':{'query':query,'url':url,'cid':entry['CID'],'smiles':smiles,'synonym_url':synurl,'retrieved_at':evidence['retrieved_at']},
                'reason':None, 'attempts':attempts}
        try:
            validate_record(candidate)
            return candidate
        except ValueError:
            continue
    ambiguous=any(len(a['response'].get('PropertyTable',{}).get('Properties',[]))>1 for a in attempts)
    return {'source_label':label,'status':'ambiguous' if ambiguous else 'unresolved','canonical_smiles':None,'evidence':None,
            'reason':'No single PubChem structure with exact supplied-name/alias synonym confirmation; see attempts', 'attempts':attempts}

def supplied_queries(label):
    label=' '.join(label.split())
    match=re.match(r'^(.*?)\s+\(([^()]+)\)(.*?)$',label)
    return list(dict.fromkeys([label.strip()] + ([' '.join((match[1]+' '+match[3]).split()),match[2].strip()] if match else [])))


def validate_record(record):
    """Validate frozen raw identity evidence offline, without trusting its digest alone."""
    label=record['source_label']
    if record['status'] not in {'resolved','unresolved','ambiguous'}:
        raise ValueError('Unknown resolution status')
    if record['status'] != 'resolved':
        if record['canonical_smiles'] is not None or record['evidence'] is not None:
            raise ValueError('Pending identity must not carry an admitted structure')
        return
    evidence=record['evidence']
    query=evidence['query']
    if query not in supplied_queries(label):
        raise ValueError('Identity query is not a supplied name')
    url=API+'name/'+quote(query,safe='')+'/property/IsomericSMILES/JSON'
    if evidence['url'] != url:
        raise ValueError('Identity URL differs from supplied query')
    attempts=[a for a in record['attempts'] if a['query']==query and a['url']==url]
    if len(attempts)!=1:
        raise ValueError('Identity raw property evidence missing or duplicated')
    attempt=attempts[0]
    props=attempt['response'].get('PropertyTable',{}).get('Properties',[])
    if len(props)!=1 or props[0].get('CID')!=evidence['cid'] or props[0].get('SMILES')!=evidence['smiles']:
        raise ValueError('Identity property evidence is not unique or differs')
    synurl=API+'cid/'+str(evidence['cid'])+'/synonyms/JSON'
    info=attempt.get('synonym_response',{}).get('InformationList',{}).get('Information',[])
    if attempt.get('synonym_url')!=synurl or evidence.get('synonym_url')!=synurl or len(info)!=1 or info[0].get('CID')!=evidence['cid']:
        raise ValueError('Identity synonym evidence differs')
    synonyms={normalize(s) for s in info[0].get('Synonym',[])}
    if normalize(query) not in synonyms:
        raise ValueError('Identity query lacks exact synonym confirmation')
    # A trade alias alone cannot erase a salt explicitly supplied by the source.
    salt_pattern=r'\b(hcl|hydrochloride|hydrobromide|chloride|bromide|sodium|disodium|potassium|calcium|sulfate|sulphate|phosphate|mesylate|besylate|tartrate|bitartrate|maleate|fumarate|succinate|salicylate|acetate|citrate|nitrate)\b'
    chemical_names=supplied_queries(label)[1:] or [label.strip()]
    if normalize(chemical_names[0]) not in synonyms:
        raise ValueError('Primary supplied chemical name lacks exact synonym confirmation')
    for chemical_name in chemical_names:
        if (re.search(salt_pattern,chemical_name,re.I) or re.search(r'\([RSrs+\-]\)|^(?:R|S)-|\b(?:cis|trans)-',chemical_name)) and normalize(chemical_name) not in synonyms:
            raise ValueError('Supplied salt/stereochemistry lacks exact synonym confirmation')
    mol=Chem.MolFromSmiles(evidence['smiles'])
    if mol is None or Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True)!=record['canonical_smiles']:
        raise ValueError('Canonical structure differs from evidence')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--retry-unresolved',action='store_true')
    args=parser.parse_args()
    out=ITEM/'local/chemical_resolution_v1.json'
    old=json.loads(out.read_text()) if out.exists() else {'records':[]}
    cache={r['source_label']:r for r in old['records']}
    for label, record in list(cache.items()):
        try:
            validate_record(record)
        except (ValueError, KeyError):
            del cache[label]
    workbook=openpyxl.load_workbook(ITEM/'source/41598_2021_96033_MOESM2_ESM.xlsx',read_only=True,data_only=True)
    labels=[r[0] for r in list(workbook.active.values)[4:]]
    todo=[s for s in labels if s not in cache or (args.retry_unresolved and cache[s]['status']!='resolved')]
    doc={'schema_version':1,'canonicalization':{'tool':'RDKit','version':rdBase.rdkitVersion,'options':'MolToSmiles(canonical=True,isomericSmiles=True); no salt stripping'},
         'policy':'One PubChem structure and exact synonym confirmation after case/punctuation/HCl normalization; only source-supplied parenthetical aliases may be queried; removing an alias retains any following salt suffix; primary chemical name must also match an exact synonym. Explicit salt/stereochemical portions require matching synonyms. Supplied salt and stereochemistry retained; no active-parent inference.', 'records':[]}
    with ThreadPoolExecutor(max_workers=8) as pool:
        for index,record in enumerate(pool.map(resolve,todo),1):
            cache[record['source_label']]=record
            doc['records']=[cache[s] for s in labels if s in cache]
            out.parent.mkdir(parents=True,exist_ok=True)
            out.write_text(json.dumps(doc,indent=2,sort_keys=True,ensure_ascii=False)+'\n')
            print(index,len(todo),record['status'],record['source_label'],flush=True)
    doc['records']=[cache[s] for s in labels]
    out.write_text(json.dumps(doc,indent=2,sort_keys=True,ensure_ascii=False)+'\n')
    logpath=ITEM/'author_log.json'
    log=json.loads(logpath.read_text())
    log['local_additions']=[a for a in log['local_additions'] if a['path']!='local/chemical_resolution_v1.json']
    log['local_additions'].append({'path':'local/chemical_resolution_v1.json','origin':'PubChem PUG REST exact-name/alias identity evidence, complete raw property and synonym responses; script/resolve_chemicals.py','version':'v1','bytes':out.stat().st_size,'sha256':hashlib.sha256(out.read_bytes()).hexdigest()})
    logpath.write_text(json.dumps(log,indent=2,ensure_ascii=False)+'\n')
    print('resolved',sum(r['status']=='resolved' for r in doc['records']),'of',len(labels))

if __name__=='__main__': main()
