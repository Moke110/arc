"""Compact dose concatenations for the dual-source Round 1 generation."""
import torch
from torch import nn

from tsg_dose import model as predecessor
from tsg_dose.model import CONCAT, NATIVE_CONFIG, SOURCE_DIGEST, configure_backend, fit_transform, transform_doses


def adapter(source, method='no-dose', transform=None):
    if method not in CONCAT:
        return predecessor.adapter(source, method, transform)
    configuration(method, transform)
    return DoseAdapter(source, method, transform)


def configuration(method, transform):
    if method not in CONCAT:
        return predecessor.configuration(method, transform)
    transform_doses([1.],transform)
    location='input' if method.startswith('input-') else 'latent'
    embedding=method.endswith('embedding')
    order=['native constructor and selective Xavier']
    if location=='latent':order.append('replace both predictive heads Linear532→100 with PyTorch defaults')
    if embedding:order.append('dose Linear1→32 ReLU Linear32→32 with PyTorch defaults')
    return dict(NATIVE_CONFIG,method=method,method_identity='tsg-dose-'+method+'-32-f32-v2',
                transform=dict(transform),dose_input=True,dose_vector_dim=32,
                features_dim=2080 if location=='input' else 2048,
                predictive_head_dim=500 if location=='input' else 532,
                dose_embedding=[1,32,32] if embedding else None,
                constructor_order=order,initialization='; '.join(order),
                adaptations=[a for a in NATIVE_CONFIG['adaptations'] if 'no pretraining/dose' not in a]+['random initialization; no pretraining/KPGT',method+' dose conditioning'])


class DoseAdapter(predecessor.DoseAdapter):
    def build(self,seed,device):
        torch.manual_seed(seed);self.device=torch.device(device)
        method=self.method
        native=self.native
        class Conditioned(native):
            def __init__(self,**kwargs):
                super().__init__(**kwargs)
                if method.startswith('latent-'):
                    self.mu_z2Fz1=nn.Sequential(nn.Linear(532,100))
                    self.logvar_z2Fz1=nn.Sequential(nn.Linear(532,100))
                if method.endswith('embedding'):
                    self.dose_embedding=nn.Sequential(nn.Linear(1,32),nn.ReLU(),nn.Linear(32,32))

            def forward(self,x1,features):
                if features.ndim!=2 or features.shape[1]!=2049:raise ValueError('dose transport must be Morgan2048 plus scalar')
                morgan,scalar=features[:,:2048],features[:,2048:]
                dose=self.dose_embedding(scalar) if method.endswith('embedding') else scalar.repeat(1,32)
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
                            features_dim=2080 if method.startswith('input-') else 2048,
                            features_embed_dim=[400],init_w=True,beta=.1,dropout=.1,
                            random_seed=seed,device=self.device)
        return network.to(device=self.device,dtype=torch.float32)
