"""Digest-verified all-dose input freeze with bounded memmapped matrices."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from . import shared
from tenv5k.inputs import sha256
from .protocol import INPUT_SCHEMA, PROTOCOL_SCHEMA, make_splits
REQUIRED={'config.json','population.csv','all_obs.csv','parents.csv','genes.csv','splits.csv','matched_controls.json','expression.npy','control.npy','delta.npy','morgan.npy'}

class Inputs:
    def __init__(self,directory):
        self.directory=Path(directory)
        self.manifest=json.loads((self.directory/'manifest.json').read_text())
        if self.manifest.get('schema')!=INPUT_SCHEMA or set(self.manifest.get('outputs',{}))!=REQUIRED:raise ValueError('input manifest schema or coverage mismatch')
        for name,digest in self.manifest['outputs'].items():
            if Path(name).name!=name or sha256(self.directory/name)!=digest:raise ValueError('input digest mismatch: '+name)
        self.config=json.loads((self.directory/'config.json').read_text())
        if self.config.get('schema')!=PROTOCOL_SCHEMA or self.config.get('split_seed')!=226 or self.config.get('folds')!=3:raise ValueError('input protocol mismatch')
        self.obs=pd.read_csv(self.directory/'population.csv').set_index('profile_id',drop=False)
        self.all_obs=pd.read_csv(self.directory/'all_obs.csv').set_index('profile_id',drop=False)
        self.parents=pd.read_csv(self.directory/'parents.csv');self.genes=pd.read_csv(self.directory/'genes.csv');self.splits=pd.read_csv(self.directory/'splits.csv')
        self.controls={int(k):v for k,v in json.loads((self.directory/'matched_controls.json').read_text()).items()}
        for name in ['expression','control','delta','morgan']:setattr(self,name,np.load(self.directory/(name+'.npy'),mmap_mode='r'))
        if self.obs.index.duplicated().any() or list(self.all_obs.index[:len(self.obs)])!=list(self.obs.index):raise ValueError('input identity mismatch')
        if self.expression.shape!=self.control.shape or self.delta.shape!=self.expression.shape or self.expression.shape!=(len(self.all_obs),len(self.genes)):raise ValueError('expression shape mismatch')
        if self.morgan.shape!=(len(self.parents),2048) or self.parents.drug_parent_inchikey.duplicated().any():raise ValueError('Morgan identity mismatch')
        if not self.splits.equals(make_splits(self.obs)):raise ValueError('partition identity or coverage mismatch')
        if len(self.obs)!=4206 or self.obs.drug_parent_inchikey.nunique()!=459 or len(self.genes)!=5000:raise ValueError('formal all-dose population mismatch')
        for dose in [.2,1.,10.]:
            selected=self.obs[self.obs.dose_um==dose]
            if len(selected)!=1402 or selected.drug_parent_inchikey.nunique()!=459:raise ValueError('per-dose population mismatch')
    def rows(self,scenario,fold,partition):
        s=self.splits
        return s.loc[(s.scenario==scenario)&(s.fold==fold)&(s.partition==partition),'row'].to_numpy(dtype=int)
    def fingerprints(self,rows):
        lookup=dict(zip(self.parents.drug_parent_inchikey,range(len(self.parents))))
        return np.asarray(self.morgan[[lookup[p] for p in self.obs.iloc[rows].drug_parent_inchikey]],dtype=np.float32)
    def doses(self,rows):return self.obs.iloc[rows].dose_um.to_numpy(dtype=np.float32)
