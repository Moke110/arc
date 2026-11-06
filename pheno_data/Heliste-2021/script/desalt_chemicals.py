"""Freeze and validate desalted Heliste parent identities (spec #54, ticket #55)."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import quote

from rdkit import Chem, rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize
from resolve_chemicals import API, ITEM, fetch, normalize, supplied_queries, validate_record as validate_original

POLICY = 'RDKit LargestFragmentChooser(preferOrganic=True), Uncharger(canonicalOrder=True), canonical isomeric SMILES; active coordination-metal drugs retain all coordination components and charges, removing only water and alkali counterions; Bacitracin zinc requires matching Bacitracin parent evidence before zinc removal; preserve permanent charges, isotopes and stereochemistry; no tautomer or stereoisomer selection'
METALS = set(range(21,31)) | set(range(39,49)) | set(range(57,81)) | {13,31,49,50,81,82,83}


def parent_smiles(smiles, zinc_salt=False):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or not molecule.GetNumAtoms():
        raise ValueError('Invalid reference structure')
    if not zinc_salt and any(a.GetAtomicNum() in METALS for a in molecule.GetAtoms()):
        retained = []
        for fragment in Chem.GetMolFrags(molecule,asMols=True,sanitizeFrags=True):
            if fragment.GetNumAtoms()==1 and fragment.GetAtomWithIdx(0).GetAtomicNum() in {3,11,19,37,55}:
                continue
            if Chem.MolToSmiles(fragment)=='O':
                continue
            retained.append(fragment)
        if not retained:
            raise ValueError('No coordination parent remains')
        parent = retained[0]
        for fragment in retained[1:]:
            parent = Chem.CombineMols(parent,fragment)
        return Chem.MolToSmiles(parent,canonical=True,isomericSmiles=True)
    parent = rdMolStandardize.LargestFragmentChooser(preferOrganic=True).choose(molecule)
    parent = rdMolStandardize.Uncharger(canonicalOrder=True).uncharge(parent)
    return Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)


def primary_name(label):
    names = supplied_queries(label)
    return names[1] if len(names)>1 else names[0]


def decide(original, extra_attempts):
    """Replay the resolution decision from frozen property and synonym responses."""
    validate_original(original)
    if original['status']=='resolved' and original['source_label']!='Bacitracin zinc':
        return 'resolved', parent_smiles(original['evidence']['smiles']), None
    primary = primary_name(original['source_label'])
    candidates = []
    # Conventional Vitamin D2 identity distinguishes the mislabeled isotope reference.
    # Bacitracin's parent query must match the zinc salt's complete stereochemistry.
    target = {'Vitamin D2':'Ergocalciferol','Bacitracin zinc':'Bacitracin'}.get(primary)
    attempts = [a for a in extra_attempts if 'query' in a] if target else original['attempts']
    for attempt in attempts:
        query = attempt['query']
        allowed = {target} if target else set(supplied_queries(original['source_label']))
        if query not in allowed or attempt['url'] != API+'name/'+quote(query,safe='')+'/property/IsomericSMILES/JSON':
            raise ValueError('Parent query differs from supported source-name interpretation')
        entries = attempt['response'].get('PropertyTable',{}).get('Properties',[])
        if target and len(entries)!=1:
            continue
        for entry in entries:
            cid = entry['CID']
            synonym = next((s for s in extra_attempts if s.get('cid')==cid),None)
            if synonym is None:
                raise ValueError('Missing candidate synonym evidence; cannot infer convergence')
            info = synonym['response'].get('InformationList',{}).get('Information',[])
            if synonym['url']!=API+'cid/'+str(cid)+'/synonyms/JSON' or len(info)!=1 or info[0].get('CID')!=cid:
                raise ValueError('Parent synonym evidence differs')
            names = {normalize(s) for s in info[0].get('Synonym',[])}
            if normalize(target or primary) not in names or normalize(query) not in names:
                continue
            candidate = parent_smiles(entry['SMILES'])
            if primary=='Bacitracin zinc' and candidate!=parent_smiles(original['evidence']['smiles'],zinc_salt=True):
                raise ValueError('Bacitracin parent differs from supplied zinc reference')
            candidates.append(candidate)
    unique = set(candidates)
    if len(unique)==1:
        return 'resolved', unique.pop(), None
    if original['source_label']=='Sevelamer HCl':
        reason = 'Sevelamer is a polymer; exact source and parent-name lookups provide no supported single finite molecular structure.'
    elif original['source_label']=='Verteporfin (Visudyne)':
        reason = 'Visudyne references differ in stereochemical specification; salt removal cannot establish one source-supported stereoisomer.'
    else:
        reason = 'Source-supported candidate structures do not establish a unique desalted parent.'
    return 'ambiguous' if original['status']=='ambiguous' else 'unresolved', None, reason


def validate_record(record):
    expected = decide(record['original_reference'],record['parent_attempts'])
    if (record['status'],record['canonical_smiles'],record['reason']) != expected:
        raise ValueError('Desalted parent differs from frozen identity evidence')
    if record['source_label'] != record['original_reference']['source_label']:
        raise ValueError('Parent label differs from original reference')


def acquire():
    original_path = ITEM/'local/chemical_resolution_v1.json'
    original = json.loads(original_path.read_text())
    records = []
    for old in original['records']:
        extra = []
        if old['status']!='resolved' or old['source_label']=='Bacitracin zinc':
            primary = primary_name(old['source_label'])
            target = {'Bacitracin zinc':'Bacitracin','Vitamin D2':'Ergocalciferol','Sevelamer HCl':'Sevelamer'}.get(primary)
            if target:
                url = API+'name/'+quote(target,safe='')+'/property/IsomericSMILES/JSON'
                extra.append({'query':target,'url':url,'response':fetch(url),'retrieved_at':datetime.now(timezone.utc).isoformat()})
            cids = {p['CID'] for a in old['attempts']+extra for p in a.get('response',{}).get('PropertyTable',{}).get('Properties',[])}
            for cid in sorted(cids):
                url = API+'cid/'+str(cid)+'/synonyms/JSON'
                extra.append({'cid':cid,'url':url,'response':fetch(url),'retrieved_at':datetime.now(timezone.utc).isoformat()})
        status,smiles,reason = decide(old,extra)
        records.append({'source_label':old['source_label'],'status':status,'canonical_smiles':smiles,'reason':reason,'original_reference':old,'parent_attempts':extra})
    document = {'schema_version':2,'original_snapshot_sha256':hashlib.sha256(original_path.read_bytes()).hexdigest(),'canonicalization':{'tool':'RDKit','version':rdBase.rdkitVersion,'options':POLICY},'policy':POLICY,'records':records}
    out = ITEM/'local/chemical_resolution_v2.json'
    out.write_text(json.dumps(document,indent=2,sort_keys=True,ensure_ascii=False)+'\n')
    log_path = ITEM/'author_log.json'
    log = json.loads(log_path.read_text())
    log['local_additions'] = [a for a in log['local_additions'] if a['path']!='local/chemical_resolution_v2.json']
    log['local_additions'].append({'path':'local/chemical_resolution_v2.json','origin':'User-authorized desalted parents, deterministic RDKit transformation of immutable v1 evidence plus targeted raw PubChem property/synonym responses; script/desalt_chemicals.py','version':'v2','bytes':out.stat().st_size,'sha256':hashlib.sha256(out.read_bytes()).hexdigest()})
    log_path.write_text(json.dumps(log,indent=2,ensure_ascii=False)+'\n')
    print('resolved',sum(r['status']=='resolved' for r in records),'of',len(records))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--acquire',action='store_true',required=True)
    parser.parse_args()
    acquire()
