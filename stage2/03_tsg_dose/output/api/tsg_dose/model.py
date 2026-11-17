"""All-dose no-dose adapter: unchanged native TranSiGen arithmetic and RNG."""
import numpy as np
import torch
from torch import nn
from . import shared
from tenv5k.model_transigen import Adapter, CONFIG as NATIVE_CONFIG, SOURCE_DIGEST, configure_backend


def adapter(source,method='no-dose',transform=None):
    configuration(method,transform)
    if method=='no-dose':return Adapter(source)
    return DoseAdapter(source,method,transform)


CONCAT = ('input-numeric','input-embedding','latent-numeric','latent-embedding')


def positive_doses(doses):
    values=np.asarray(doses,dtype=np.float64)
    if values.ndim!=1 or not len(values) or not np.isfinite(values).all() or (values<=0).any():
        raise ValueError('doses must be finite positive one-dimensional values')
    return values


def fit_transform(method,doses):
    values=positive_doses(doses)
    if method=='no-dose':return None
    if method not in CONCAT:raise ValueError('unsupported dose method')
    logs=np.log10(values);scale=float(logs.std(ddof=0))
    if not np.isfinite(scale) or scale<=0:raise ValueError('training dose scale must be positive')
    return dict(kind='training-log10-dose-population-standardization-v1',mean=float(logs.mean()),std=scale)


def transform_doses(doses,transform):
    values=positive_doses(doses)
    if not isinstance(transform,dict) or set(transform)!={'kind','mean','std'} or transform['kind']!='training-log10-dose-population-standardization-v1':
        raise ValueError('invalid dose transform identity')
    mean,scale=transform['mean'],transform['std']
    if not np.isfinite([mean,scale]).all() or scale<=0:raise ValueError('invalid dose transform values')
    result=np.asarray((np.log10(values)-mean)/scale,dtype=np.float32)
    if not np.isfinite(result).all():raise ValueError('nonfinite standardized dose')
    return result


def configuration(method,transform):
    if method=='no-dose':
        if transform is not None:raise ValueError('no-dose requires no transform')
        return dict(NATIVE_CONFIG,method=method,method_identity='tsg-dose-no-dose-native-f32-v1',transform=transform,
                    initialization='native selective Xavier; unchanged no-dose constructor order')
    if method not in CONCAT:raise ValueError('unsupported dose method')
    transform_doses([1.],transform)
    location='input' if method.startswith('input-') else 'latent'
    embedding=method.endswith('embedding')
    order=['native constructor and selective Xavier']
    if location=='latent':order.append('replace both predictive heads Linear600→100 with PyTorch defaults')
    if embedding:order.append('dose Linear1→100 ReLU Linear100→100 with PyTorch defaults')
    return dict(NATIVE_CONFIG,method=method,method_identity='tsg-dose-'+method+'-100-f32-v1',
                transform=dict(transform),dose_input=True,dose_vector_dim=100,
                features_dim=2148 if location=='input' else 2048,
                predictive_head_dim=500 if location=='input' else 600,
                dose_embedding=[1,100,100] if embedding else None,
                constructor_order=order,initialization='; '.join(order),
                adaptations=[a for a in NATIVE_CONFIG['adaptations'] if 'no pretraining/dose' not in a]+['random initialization; no pretraining/KPGT',method+' dose conditioning'])


class DoseView:
    """Pack dose beside Morgan for the adapter; never expose targets to forward."""
    def __init__(self,data,transform):self.data,self.transform=data,transform
    def __getattr__(self,name):return getattr(self.data,name)
    def fingerprints(self,rows):
        scalar=transform_doses(self.data.doses(rows),self.transform)
        return np.concatenate([self.data.fingerprints(rows),scalar[:,None]],axis=1).astype(np.float32)


class DoseAdapter(Adapter):
    def __init__(self,source,method,transform):
        super().__init__(source)
        self.method,self.transform=method,dict(transform)

    def build(self,seed,device):
        torch.manual_seed(seed);self.device=torch.device(device)
        method=self.method
        native=self.native
        class Conditioned(native):
            def __init__(self,**kwargs):
                super().__init__(**kwargs)
                if method.startswith('latent-'):
                    self.mu_z2Fz1=nn.Sequential(nn.Linear(600,100))
                    self.logvar_z2Fz1=nn.Sequential(nn.Linear(600,100))
                if method.endswith('embedding'):
                    self.dose_embedding=nn.Sequential(nn.Linear(1,100),nn.ReLU(),nn.Linear(100,100))

            def forward(self,x1,features):
                if features.ndim!=2 or features.shape[1]!=2049:raise ValueError('dose transport must be Morgan2048 plus scalar')
                morgan,scalar=features[:,:2048],features[:,2048:]
                dose=self.dose_embedding(scalar) if method.endswith('embedding') else scalar.repeat(1,100)
                if method.startswith('input-'):
                    return super().forward(x1,torch.cat([morgan,dose],1))
                # Preserve native control/reconstruction/chemical/head/sample/decode order.
                z1,mu1,logvar1=self.encode_x1(x1)
                x1_rec=self.decode_x1(z1)
                chemical=self.feat_embeddings(morgan)
                conditioned=torch.cat([z1,chemical,dose],1)
                mu_pred,logvar_pred=self.mu_z2Fz1(conditioned),self.logvar_z2Fz1(conditioned)
                z2_pred=self.sample_latent(mu_pred,logvar_pred)
                return x1_rec,mu1,logvar1,self.decode_x2(z2_pred),mu_pred,logvar_pred,z2_pred

        network=Conditioned(n_genes=5000,n_latent=100,n_en_hidden=[1200],n_de_hidden=[800],
                            features_dim=2148 if method.startswith('input-') else 2048,
                            features_embed_dim=[400],init_w=True,beta=.1,dropout=.1,
                            random_seed=seed,device=self.device)
        return network.to(device=self.device,dtype=torch.float32)

    def predict(self,model,inputs,rows,batch_size=64):
        return super().predict(model,DoseView(inputs,self.transform),rows,batch_size)

    def train_epoch(self,model,inputs,rows,optimizer,rng):
        return super().train_epoch(model,DoseView(inputs,self.transform),rows,optimizer,rng)
