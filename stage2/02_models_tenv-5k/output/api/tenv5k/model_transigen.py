"""Native TranSiGen with declared full-float64 arithmetic variant.

Only the interface gene dimension changes to 5000. Prediction uses forward(x1,
features)[3], never the target encoder; native float32 CPU draws precede casting.
"""
import hashlib
import importlib.util
from pathlib import Path
import sys
from types import MethodType
import numpy as np
import torch
from .native_import import no_source_bytecode

SCRIPT_HASH = '4446df722b2dcd9633ae30253ee077d6836ee86862ddd72ecfa1460060c19eb5'
SOURCE_HASH = 'c328864a040a3c1a1e69b040f31a27572513b1974a93d722e3ff00d4cdd940c7'


def configure_backend(device):
    """Match ordinary author setup_seed CUDA policy; leave CPU flags untouched."""
    if torch.device(device).type != 'cuda':
        return {'policy': 'CPU backend unchanged'}
    torch.backends.cudnn.enabled = False
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    return dict(policy='author train_TranSiGen setup_seed CUDA',
                cudnn_enabled=torch.backends.cudnn.enabled,
                cudnn_benchmark=torch.backends.cudnn.benchmark,
                cudnn_deterministic=torch.backends.cudnn.deterministic)

def sample_latent_f64(model, mu, logvar):
    """Preserve native CPU float32 normal draws, then compute in model dtype."""
    std = logvar.mul(0.5).exp_()
    eps = torch.FloatTensor(std.size()).normal_().to(device=mu.device, dtype=mu.dtype)
    return eps.mul_(std).add_(mu)


class Adapter:
    name = 'TranSiGen'
    batch_size = 64

    def __init__(self, source):
        self.source = Path(source).resolve()
        path = self.source / 'src/model.py'
        if hashlib.sha256((self.source / 'src/train_TranSiGen.py').read_bytes()).hexdigest() != SCRIPT_HASH:
            raise ValueError('TranSiGen config script digest mismatch')
        if hashlib.sha256(path.read_bytes()).hexdigest() != SOURCE_HASH:
            raise ValueError('TranSiGen source digest mismatch')
        # Upstream imports its sibling utils by an unqualified name. Restore the
        # process module namespace after loading so other model sources coexist.
        previous = sys.modules.pop('utils', None)
        sys.path.insert(0, str(path.parent))
        try:
            spec = importlib.util.spec_from_file_location('_tenv5k_native_transigen', path)
            module = importlib.util.module_from_spec(spec)
            with no_source_bytecode():
                spec.loader.exec_module(module)
            self.native = module.TranSiGen
        finally:
            sys.path.pop(0)
            sys.modules.pop('utils', None)
            if previous is not None:
                sys.modules['utils'] = previous

    def build(self, seed, device):
        torch.manual_seed(seed)
        self.device = torch.device(device)
        model = self.native(n_genes=5000, n_latent=100, n_en_hidden=[1200],
                           n_de_hidden=[800], features_dim=2048,
                           features_embed_dim=[400], init_w=True, beta=0.1,
                           dropout=0.1, random_seed=seed, device=self.device)
        # Native selective initialization and its RNG draws occur in float32.
        model = model.double().to(self.device)
        model.sample_latent = MethodType(sample_latent_f64, model)
        return model

    def optimizer(self, model):
        return torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)

    def _tensor(self, value, model):
        return torch.tensor(np.asarray(value), dtype=next(model.parameters()).dtype, device=next(model.parameters()).device)

    def predict(self, model, inputs, rows, batch_size=64):
        model.eval()
        result = []
        with torch.no_grad():
            for start in range(0, len(rows), batch_size):
                batch = rows[start:start + batch_size]
                x1 = self._tensor(inputs.control[batch], model)
                features = self._tensor(inputs.fingerprints(batch), model)
                result.append(model(x1, features)[3].cpu().numpy())
        return np.concatenate(result) if result else np.empty((0, 5000), np.float64)

    def train_epoch(self, model, inputs, rows, optimizer, rng):
        model.train()
        order = rng.permutation(rows)
        total, count = 0.0, 0
        for start in range(0, len(order), self.batch_size):
            batch = order[start:start + self.batch_size]
            if len(batch) == 1:
                continue
            x1 = self._tensor(inputs.control[batch], model)
            x2 = self._tensor(inputs.expression[batch], model)
            features = self._tensor(inputs.fingerprints(batch), model)
            optimizer.zero_grad()
            x1_rec, mu1, logvar1, pred, mu_pred, logvar_pred, _ = model(x1, features)
            z2, mu2, logvar2 = model.encode_x2(x2)
            x2_rec = model.decode_x2(z2)
            loss = model.loss(x1, x1_rec, mu1, logvar1, x2, x2_rec,
                              mu2, logvar2, pred, mu_pred, logvar_pred)[0]
            if not torch.isfinite(loss):
                raise ValueError('non-finite native TranSiGen loss')
            loss.backward()
            optimizer.step()
            total += loss.item()
            count += len(batch)
        if not count:
            raise ValueError('no non-singleton TranSiGen training batch')
        return total / count

# Shared learner's functional adapter contract.
CONFIG = dict(method='TranSiGen', method_identity='transigen-native-f64-v1',
              compute_dtype='float64', initialization_dtype='float32',
              prediction_dtype='float64', noise_dtype='float32 CPU then cast', batch_size=64, n_genes=5000, n_latent=100,
              encoder=[1200], decoder=[800], features_dim=2048,
              features_embed_dim=[400], init_w=True, beta=0.1, dropout=0.1,
              molecule_feature='ECFP4', initialization_model='random',
              learning_rate=1e-3, weight_decay=1e-5, scheduler=None,
              dose_input=False, prediction='native forward(x1,features)[3]',
              singleton_train='skip singleton as upstream',
              upstream_commit='e7582e7794f13d4dfddfb15e95e3ba7cd24bc07a',
              config_script='src/train_TranSiGen.py', config_script_sha256=SCRIPT_HASH,
              singleton_train_batch='skip singleton as upstream',
              adaptations=['5000 gene interface', 'shared Morgan2048 ECFP4',
                           'random initialization; no pretraining/dose/KPGT',
                           'shared splits and fixed validation treated-expression MSE',
                           'x1-only native sampled forward for inference',
                           'native float32 initialization then full float64 computation',
                           'native float32 CPU normal draw then cast before scale/add'])
SOURCE_DIGEST = SOURCE_HASH

def build_model(source, seed, device):
    adapter = Adapter(source)
    model = adapter.build(seed, device)
    model._tenv_adapter = adapter
    return model

def make_optimizer(model):
    return model._tenv_adapter.optimizer(model)

def make_scheduler(optimizer):
    return None

def train_epoch(model, data, rows, optimizer, shuffle, device):
    return model._tenv_adapter.train_epoch(model, data, rows, optimizer, shuffle)

def sampled_prediction(model, data, rows, seed=42, device='cpu', batch_size=64):
    from .learning import fixed_rng
    with fixed_rng(seed, device):
        return model._tenv_adapter.predict(model, data, rows, batch_size)
