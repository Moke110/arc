"""Freeze the public 5k single-dose population, Morgan inputs and shared splits."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np
import pandas as pd


def repo_root(path):
    for root in (Path(path).resolve(),*Path(path).resolve().parents):
        if (root/'.git').exists():return root
    raise ValueError('repository root not found')

ROOT=repo_root(__file__)
STEP=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'stage2/01_bases/output/api'))
sys.path.insert(0,str(STEP/'output/api'))
from bases.view import open_view
from tenv5k.protocol import make_splits, PROTOCOL_SCHEMA, INPUT_SCHEMA, FOLDS


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n')


def prepare(view_dir,destination):
    """Candidate-only preparation; stable publication is owned by the matrix driver."""
    from rdkit import Chem,DataStructs,rdBase
    from rdkit.Chem import rdFingerprintGenerator
    destination=Path(destination)
    if destination.exists() and any(destination.iterdir()):raise ValueError('candidate must be empty')
    destination.mkdir(parents=True,exist_ok=True)
    view_dir=Path(view_dir)
    with open_view(view_dir) as view:
        source=view.source('tenv');obs=source.obs
        selected=np.flatnonzero((obs.role=='treated')&(obs.dose_um==1))
        population=obs.iloc[selected].copy()
        if len(population)!=1402 or population.drug_parent_inchikey.nunique()!=459 or source.axis_size!=5000:
            raise ValueError('formal population differs from 1402 wells / 459 parents / 5000 genes')
        if set(population.cell_line)!={'hiPSC-CM'}:raise ValueError('unexpected cell context')
        controls={i:[int(r) for r in source.matched_control_rows(int(row))] for i,row in enumerate(selected)}
        original_controls=sorted({r for rows in controls.values() for r in rows})
        all_rows=np.r_[selected,original_controls]
        # Memmaps and bounded public windows keep the whole view out of RAM.
        arrays={name:np.lib.format.open_memmap(destination/(name+'.npy'),mode='w+',dtype='float32',shape=(len(all_rows),5000)) for name in ['expression','control','delta']}
        for target,row in enumerate(all_rows):
            for name,reader in [('expression',source.X),('control',source.control),('delta',source.delta)]:arrays[name][target]=reader(int(row),int(row)+1)[0]
        for array in arrays.values():array.flush()
        del arrays
        control_lookup={r:len(selected)+i for i,r in enumerate(original_controls)}
        mapped={str(i):[control_lookup[r] for r in rows] for i,rows in controls.items()}
        population['view_row']=selected
        population['profile_id']=list(population.index)
        population.to_csv(destination/'population.csv',index=False)
        full_obs=obs.iloc[all_rows].copy();full_obs['profile_id']=list(full_obs.index)
        full_obs.to_csv(destination/'all_obs.csv',index=False)
        pd.DataFrame({'gene_index':source.gene_index}).to_csv(destination/'genes.csv',index=False)
        write_json(destination/'matched_controls.json',mapped)
        make_splits(population).to_csv(destination/'splits.csv',index=False)
        parents=population[['drug_parent_inchikey','drug_parent_smiles']].drop_duplicates().sort_values('drug_parent_inchikey')
        if parents.drug_parent_inchikey.duplicated().any():raise ValueError('parent has conflicting structures')
        generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
        fingerprints=[]
        for smiles in parents.drug_parent_smiles:
            mol=Chem.MolFromSmiles(smiles)
            if mol is None:raise ValueError('invalid parent SMILES')
            bits=np.zeros(2048,dtype=np.uint8);DataStructs.ConvertToNumpyArray(generator.GetFingerprint(mol),bits);fingerprints.append(bits)
        np.save(destination/'morgan.npy',np.asarray(fingerprints,dtype=np.uint8))
        parents.to_csv(destination/'parents.csv',index=False)
        manifest=view.manifest()
        if digest(view_dir/'tenv.h5ad')!=manifest['outputs']['tenv.h5ad']['sha256']:raise ValueError('view payload digest mismatch')
    config={'schema':PROTOCOL_SCHEMA,'split_seed':226,'statistical_seed':2262026,'training_seeds':[17,42,2026],'folds':FOLDS,'validation_fraction':0.1,'rounding':'ceil outer-training units; parents cold, wells warm','dose_um':1,'axis_size':5000,'early_stopping':{'monitor':'complete validation treated-expression MSE','strict_lower':True,'patience':50,'max_epochs':1000,'restore_best':True,'validation_seed':22617},'morgan':{'radius':2,'bits':2048,'includeChirality':True,'binary':True,'rdkit':rdBase.rdkitVersion},'n_treated':1402,'n_parents':459,'n_unique_controls':len(original_controls),'view_manifest_sha256':digest(view_dir/'manifest.json'),'view_tenv_sha256':digest(view_dir/'tenv.h5ad')}
    write_json(destination/'config.json',config)
    write_json(destination/'manifest.json',{'schema':INPUT_SCHEMA,'inputs':{'view_manifest':config['view_manifest_sha256'],'view_tenv':config['view_tenv_sha256']},'outputs':{p.name:digest(p) for p in sorted(destination.iterdir()) if p.is_file()}})
    return config

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--view',type=Path,default=ROOT/'stage2/01_bases/output/view-5k')
    parser.add_argument('--candidate',type=Path,required=True)
    args=parser.parse_args();print(json.dumps(prepare(args.view,args.candidate),indent=2))
