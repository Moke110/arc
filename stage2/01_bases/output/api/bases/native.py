"""Standalone bounded raw-count preparation from immutable GSE262419 plates.

Only one source plate is materialized. Counts retain native HGNC order, including
controls and technical blanks; projection/normalization belongs to training build.
"""
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem
from scipy.sparse import csr_matrix
from .metadata.annotation import HgncIndex, CompoundAnnotation
from .metadata.genes import map_probe_space
from .metadata.text import normalize_name

SCHEMA='stage2-bases-native-tenv-v1'

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def well_name(value):
    return re.sub(r'^([A-Z]+)0*(\d+)$',r'\1\2',str(value))

def source_plan(source,hgnc,annotation):
    source=Path(source); plates=sorted(source.glob('GSE262419_Plate*.csv.gz'))
    if not plates:raise ValueError('no immutable tenv plate files')
    index=HgncIndex.load(Path(hgnc)); registry=CompoundAnnotation.load(Path(annotation))
    probes={}
    for plate in plates:
        for gene in pd.read_csv(plate,usecols=['Genes']).Genes.astype(str):probes.setdefault(gene,len(probes))
    mapping=map_probe_space('temp_env',list(probes),index)
    identities={}; matches={}
    for identity,match in zip(mapping.hgnc_ids,mapping.match_types,strict=True):
        if identity is not None and identity not in identities:
            identities[identity]=len(identities);matches[identity]=match
    ids=list(identities)
    var=pd.DataFrame({'hgnc_id':ids,'entrez_id':[index.ref(h).entrez_id for h in ids],
        'ensembl_id':[index.ref(h).ensembl_id for h in ids],'match_type':[matches[h] for h in ids]},
        index=pd.Index([index.ref(h).symbol for h in ids],name='symbol'))
    probe_columns={raw:identities.get(h,-1) for raw,h in zip(mapping.raw_ids,mapping.hgnc_ids,strict=True)}
    hashes=pd.read_csv(source/'GSE262419_hash.csv.gz')
    info={(int(r.Plate_ID),well_name(r.Well_ID)):r for r in hashes.itertuples(index=False)}
    records=[]; selections=[]
    for plate in plates:
        columns=list(pd.read_csv(plate,nrows=0).columns[1:]);selected=[]
        for position,column in enumerate(columns):
            plate_id,_,well=column.partition('-'); entry=info.get((int(plate_id.replace('Plate','')),well_name(well)))
            if entry is None:continue
            name=normalize_name(entry.Chemical_name);control=name in ('VEH','Media');blank=name=='Blank'
            name='DMSO' if name=='VEH' else name
            ann=registry.get(name) if not(control or blank) else None
            # Native publication retained only resolved small molecules plus controls/blanks.
            if not(control or blank) and (name.upper()=='ET-1' or not ann or not ann.get('parent_inchikey')):continue
            atc=(ann or {}).get('atc','');moa=(ann or {}).get('moa','');mesh=(ann or {}).get('mesh_classes','')
            smiles=(ann or {}).get('smiles','')
            if smiles:
                from rdkit import Chem
                mol=Chem.MolFromSmiles(smiles)
                if mol is None:raise ValueError(f'invalid registered structure: {name}')
                smiles=Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True)
            records.append(dict(well=f'temp_env::{column}',dataset_id='tenv',compound_name=name,
                drug_key=(ann or {}).get('parent_inchikey',''),inchikey=(ann or {}).get('inchikey','') or (ann or {}).get('parent_inchikey',''),smiles=smiles,
                dose_um=0. if control or blank else float(entry.Chemical_Concentration_uM),dose_source='' if control or blank else 'measured',
                time_h=24.,is_control=control,vehicle='Media' if name=='Media' else 'DMSO',plate=plate_id,control_group='tenv::'+plate_id,
                cell_type='ipsc_cm',donor='NA',atc=atc,moa=moa,mesh_classes=mesh,pharm_atc_moa=bool(atc or moa),pharm_mesh=bool(atc or moa or mesh),
                annotation_source='chembl' if atc or moa else ('pubchem' if mesh else ''),cell_line='NA',n_raw_rows=1,is_blank=blank))
            selected.append(position)
        selections.append((plate,selected))
    obs=pd.DataFrame(records)
    groups=np.where(obs.is_control,obs.control_group+'|control|'+obs.vehicle,obs.control_group+'|'+obs.drug_key+'|'+obs.dose_um.astype(str))
    order=np.lexsort((obs.well.to_numpy(),groups));rep=np.zeros(len(obs),dtype='int32');previous=None;count=0
    for pos in order:
        if groups[pos]!=previous:count=0;previous=groups[pos]
        rep[pos]=count;count+=1
    obs['replicate']=rep
    # Preserve the native public observation schema and column order.
    fields=['well','dataset_id','compound_name','drug_key','inchikey','smiles','dose_um','dose_source','time_h','is_control','vehicle','plate','control_group','cell_type','donor','replicate','atc','moa','mesh_classes','pharm_atc_moa','pharm_mesh','annotation_source','cell_line','n_raw_rows','is_blank']
    obs=obs[fields];obs.index=pd.Index(obs.well,name='well')
    if obs.index.has_duplicates:raise ValueError('duplicate source well')
    evidence=pd.DataFrame({'source_probe':mapping.raw_ids,'hgnc_id':[h or '' for h in mapping.hgnc_ids],'match_type':mapping.match_types,'native_column':[probe_columns[p] for p in mapping.raw_ids]})
    return obs,var,probe_columns,selections,evidence

def plate_batches(selections,probe_columns,n_genes):
    for path,selected in selections:
        frame=pd.read_csv(path)
        target=np.array([probe_columns[p] for p in frame.Genes.astype(str)])
        keep=target>=0
        project=csr_matrix((np.ones(keep.sum(),dtype='int64'),(np.arange(keep.sum()),target[keep])),shape=(keep.sum(),n_genes))
        for start in range(0,len(selected),64):
            positions=selected[start:start+64]
            raw=frame.iloc[:,np.array(positions)+1].to_numpy().T
            if not np.isfinite(raw).all() or (raw<0).any() or (raw!=np.floor(raw)).any():raise ValueError(f'noninteger raw counts: {path}')
            counts=np.asarray(raw[:,keep].astype('int64')@project)
            if counts.max(initial=0)>np.iinfo(np.int32).max:raise ValueError('native counts overflow')
            yield counts.astype('int32')

def input_seals(source,hgnc,annotation,selections):
    files=[Path(hgnc),Path(annotation),Path(source)/'GSE262419_hash.csv.gz',*[p for p,_ in selections]]
    return {p.name:{'sha256':digest(p),'bytes':p.stat().st_size} for p in files}

def build_tenv(source,hgnc,annotation,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    destination=output/'tenv.h5ad'
    if destination.exists():raise ValueError(f'native corpus exists: {destination}')
    obs,var,columns,selections,evidence=source_plan(source,hgnc,annotation)
    seals=input_seals(source,hgnc,annotation,selections)
    fd,tmp=tempfile.mkstemp(prefix='.tenv-',suffix='.h5ad',dir=output);os.close(fd)
    try:
        with h5py.File(tmp,'w') as f:
            f.attrs.update({'encoding-type':'anndata','encoding-version':'0.1.0'})
            write_elem(f,'var',var);write_elem(f,'uns',{'schema':SCHEMA,'counts':'raw counts summed by resolved native HGNC; no normalization'})
            write_elem(f,'layers',{})
            x=f.create_dataset('X',shape=(len(obs),len(var)),dtype='int32',chunks=(min(64,len(obs)),len(var)),compression='gzip',compression_opts=1,shuffle=True)
            x.attrs.update({'encoding-type':'array','encoding-version':'0.2.0'})
            totals=[];detected=[];offset=0
            for block in plate_batches(selections,columns,len(var)):
                x[offset:offset+len(block)]=block;offset+=len(block)
                totals.extend(block.sum(axis=1,dtype='int64'));detected.extend((block>0).sum(axis=1))
            if offset!=len(obs) or min(totals)<=0:raise ValueError('invalid native row accounting/library')
            obs.insert(22,'lib_size',np.array(totals,dtype='int64'));obs.insert(23,'n_genes_detected',np.array(detected,dtype='int32'))
            write_elem(f,'obs',obs)
        os.replace(tmp,destination)
    finally:
        if Path(tmp).exists():Path(tmp).unlink()
    evidence.to_csv(output/'tenv_probe_mapping.csv',index=False)
    result={'schema':SCHEMA,'shape':[len(obs),len(var)],'inputs':seals,'files':{'tenv.h5ad':digest(destination),'tenv_probe_mapping.csv':digest(output/'tenv_probe_mapping.csv')},'populations':{'treated':int((~obs.is_control&~obs.is_blank).sum()),'control':int(obs.is_control.sum()),'blank':int(obs.is_blank.sum())}}
    (output/'tenv_manifest.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    return result

def verify_tenv(source,hgnc,annotation,output):
    output=Path(output);manifest=json.loads((output/'tenv_manifest.json').read_text())
    obs,var,columns,selections,evidence=source_plan(source,hgnc,annotation)
    if manifest['inputs']!=input_seals(source,hgnc,annotation,selections):raise ValueError('immutable input seals differ')
    for name,value in manifest['files'].items():
        if digest(output/name)!=value:raise ValueError(f'native payload hash differs: {name}')
    # Independent source replay sums each probe directly, rather than matrix multiplication.
    with h5py.File(output/'tenv.h5ad') as f:
        actual_obs=read_elem(f['obs']);actual_var=read_elem(f['var'])
        pd.testing.assert_frame_equal(actual_obs[obs.columns],obs,check_dtype=False,check_categorical=False)
        pd.testing.assert_frame_equal(actual_var,var,check_dtype=False,check_categorical=False)
        offset=0
        for path,selected in selections:
            frame=pd.read_csv(path);targets=np.array([columns[p] for p in frame.Genes.astype(str)])
            for start in range(0,len(selected),64):
                positions=selected[start:start+64];raw=frame.iloc[:,np.array(positions)+1].to_numpy().T
                expected=np.zeros((len(positions),len(var)),dtype='int64')
                for probe,target in enumerate(targets):
                    if target>=0:expected[:,target]+=raw[:,probe].astype('int64')
                np.testing.assert_array_equal(f['X'][offset:offset+len(positions)],expected)
                np.testing.assert_array_equal(actual_obs.lib_size.iloc[offset:offset+len(positions)],expected.sum(axis=1))
                offset+=len(positions)
    return {'result':'PASS','source_rows_replayed':offset,'shape':manifest['shape'],'manifest_sha256':digest(output/'tenv_manifest.json')}
