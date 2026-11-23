"""Consume an accepted three-checkpoint release without training matrices."""
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from . import model as models
from .inputs import sha256
from .learner import runtime_identity,environment,predict
from .protocol import METHODS,SEEDS

ROLE='internal-validation-selected restored checkpoints; no new independent test'


class Release:
    def __init__(self,directory,source):
        self.directory=Path(directory)
        self.manifest=json.loads((self.directory/'manifest.json').read_text())
        if self.manifest.get('n_formal_fits')!=147:
            raise ValueError('release requires complete147-fit acceptance')
        self.source=Path(source)
        self._verify_files()
        self.contract=json.loads((self.directory/'input_contract.json').read_text())
        self.genes=pd.read_csv(self.directory/'genes.csv',dtype=str).gene_index.tolist()
        if len(self.genes)!=5000 or len(set(self.genes))!=5000 or self.contract.get('axis_size')!=5000:
            raise ValueError('release gene axis must contain5000 unique ordered genes')
        morgan=self.contract.get('morgan',{})
        if any(morgan.get(k)!=v for k,v in dict(radius=2,bits=2048,includeChirality=True,binary=True).items()):
            raise ValueError('release Morgan contract mismatch')
        if self.contract.get('input_manifest_sha256')!=self.manifest.get('input_manifest_sha256'):
            raise ValueError('release input contract mismatch')
        receipt=json.loads((self.directory/'selection_receipt.json').read_text())
        if any(receipt.get(k)!=v for k,v in dict(schema='tsg-dose-selection-v1',scenario='final',fold=0,method=self.manifest['method'],input_manifest_sha256=self.manifest['input_manifest_sha256']).items()):
            raise ValueError('release final selection mismatch')
        # Validate each constructor/state/transform independently; never retain
        # three model copies or require their original training matrices.
        for seed in SEEDS:self._load(seed)

    def _verify_files(self):
        manifest=self.manifest
        if manifest.get('schema')!='tsg-dose-release-v1' or manifest.get('scenario_counts')!={'cold-drug':63,'held-dose':63,'final':21} or manifest.get('seeds')!=list(SEEDS) or manifest.get('method') not in METHODS or manifest.get('role')!=ROLE:
            raise ValueError('release identity or complete matrix coverage mismatch')
        required={'genes.csv','input_contract.json','selection_receipt.json'}|{f'checkpoints/{s}.pt' for s in SEEDS}
        actual={str(p.relative_to(self.directory)) for p in self.directory.rglob('*') if p.is_file() and p!=self.directory/'manifest.json'}
        if set(manifest.get('outputs',{}))!=required or actual!=required:
            raise ValueError('release payload coverage mismatch')
        for name,digest in manifest['outputs'].items():
            path=self.directory/name
            if path.is_symlink() or sha256(path)!=digest:raise ValueError('release payload digest mismatch: '+name)

    def _load(self,seed):
        saved=torch.load(self.directory/'checkpoints'/f'{seed}.pt',map_location='cpu',weights_only=False)
        config=saved['config']
        expected=models.configuration(self.manifest['method'],config.get('transform'))
        identity=dict(scenario='final',fold=0,seed=seed,method=self.manifest['method'],smoke_rows=None,max_epochs=1000,patience=50,validation_seed=22617,inference_seed=seed+100000,native_source_sha256=models.SOURCE_DIGEST,runtime_identity=runtime_identity())
        if any(config.get(k)!=v for k,v in {**expected,**identity}.items()) or saved.get('input_manifest_sha256')!=self.manifest['input_manifest_sha256']:
            raise ValueError('release checkpoint configuration/input mismatch')
        if self.manifest.get('environments',{}).get(str(seed))!=environment(config['device']):raise ValueError('release checkpoint environment mismatch')
        os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8');torch.use_deterministic_algorithms(True)
        if config.get('backend')!=models.configure_backend(config['device']):raise ValueError('release checkpoint backend mismatch')
        state=saved['state_dict']
        if any(v.is_floating_point() and (v.dtype!=torch.float32 or not torch.isfinite(v).all()) for v in state.values()):raise ValueError('release checkpoint precision/nonfinite state mismatch')
        native=models.adapter(self.source,config['method'],config['transform']);network=native.build(seed,config['device'])
        if sum(p.numel() for p in network.parameters())!=self.manifest.get('parameter_counts',{}).get(str(seed)):raise ValueError('release parameter count mismatch')
        network.load_state_dict(state,strict=True)
        return native,network,config

    def predict(self,seed,*,controls,morgan,doses,gene_ids,profile_ids,sampling_seed=None):
        """Return ordered float32 treated expression; controls are Bases-normalized.

        Morgan2048 must use radius2/includeChirality binary fingerprints. Doses
        are positive µM. Native batch64 and sampled latent inference are retained;
        repeat exact rows/order/seed to reproduce a prediction. No labels enter.
        """
        if seed not in SEEDS:raise ValueError('unknown release training seed')
        if [str(v) for v in gene_ids]!=self.genes:raise ValueError('gene axis order mismatch')
        identities=pd.Index(profile_ids)
        if not len(identities) or identities.has_duplicates or identities.isna().any():raise ValueError('profile identities must be nonempty and unique')
        controls=np.asarray(controls,dtype=np.float32);morgan=np.asarray(morgan,dtype=np.float32)
        doses=models.positive_doses(doses)
        if controls.shape!=(len(identities),5000) or morgan.shape!=(len(identities),2048) or len(doses)!=len(identities):raise ValueError('prediction input shape mismatch')
        if not np.isfinite(controls).all() or not np.isfinite(morgan).all() or not np.isin(morgan,[0,1]).all():raise ValueError('prediction controls/Morgan must be finite; Morgan must be binary')
        self._verify_files();native,network,config=self._load(seed)
        class Data:
            control=controls
            def fingerprints(self,rows):return morgan[rows]
            def doses(self,rows):return doses[rows]
        result=predict(native,network,Data(),np.arange(len(identities)),config['inference_seed'] if sampling_seed is None else sampling_seed)
        return pd.DataFrame(result,index=identities,columns=self.genes)
