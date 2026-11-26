"""Native float32 TranSiGen with the matched-control study’s KPGT and dose inputs."""
import hashlib
import importlib.util
from operator import index
from pathlib import Path
import sys
import numpy as np
import torch
from .native_import import no_source_bytecode
from stage2_util.representations import WIDTHS

SCRIPT_HASH = '4446df722b2dcd9633ae30253ee077d6836ee86862ddd72ecfa1460060c19eb5'
SOURCE_DIGEST = 'c328864a040a3c1a1e69b040f31a27572513b1974a93d722e3ff00d4cdd940c7'


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


class Adapter:
    name = 'TranSiGen'
    batch_size = 64

    def __init__(self, source, transform):
        self.transform = dict(transform)
        self.source = Path(source).resolve()
        path = self.source / 'src/model.py'
        if hashlib.sha256((self.source / 'src/train_TranSiGen.py').read_bytes()).hexdigest() != SCRIPT_HASH:
            raise ValueError('TranSiGen config script digest mismatch')
        if hashlib.sha256(path.read_bytes()).hexdigest() != SOURCE_DIGEST:
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

    def build(self, seed, device, *, n_genes=5000):
        n_genes = index(n_genes)
        if n_genes < 1:
            raise ValueError("TranSiGen needs a positive gene-axis width")
        torch.manual_seed(seed)
        self.device = torch.device(device)
        model = self.network()(n_genes=n_genes, n_latent=100, n_en_hidden=[1200],
                               n_de_hidden=[800], features_dim=WIDTH,
                               features_embed_dim=[400], init_w=True, beta=0.1,
                               dropout=0.1, random_seed=seed, device=self.device)
        self.n_genes = n_genes
        # Native selective initialization, latent sampling and arithmetic stay float32.
        return model.to(device=self.device, dtype=torch.float32)

    def network(self):
        class Conditioned(self.native):
            def __init__(self, **kwargs):
                super().__init__(**kwargs)
                self.mu_z2Fz1 = torch.nn.Sequential(torch.nn.Linear(532, 100))
                self.logvar_z2Fz1 = torch.nn.Sequential(torch.nn.Linear(532, 100))

            def forward(self, x1, features):
                if features.ndim != 2 or features.shape[1] != WIDTH + 1:
                    raise ValueError("dose transport requires KPGT2304 plus one scalar")
                kpgt, scalar = features[:, :WIDTH], features[:, WIDTH:]
                dose = scalar.repeat(1, 32)
                z1, mu1, logvar1 = self.encode_x1(x1)
                x1_rec = self.decode_x1(z1)
                chemical = self.feat_embeddings(kpgt)
                conditioned = torch.cat([z1, chemical, dose], dim=1)
                mu_pred = self.mu_z2Fz1(conditioned)
                logvar_pred = self.logvar_z2Fz1(conditioned)
                z2_pred = self.sample_latent(mu_pred, logvar_pred)
                return x1_rec, mu1, logvar1, self.decode_x2(z2_pred), mu_pred, logvar_pred, z2_pred

        return Conditioned

    def optimizer(self, model):
        return torch.optim.Adam(model.parameters(), lr=1e-4, weight_decay=1e-5)

    def _tensor(self, value, model):
        return torch.tensor(np.asarray(value), dtype=next(model.parameters()).dtype, device=next(model.parameters()).device)

    def predict(self, model, inputs, rows, batch_size=64):
        inputs = DoseView(inputs, self.transform)
        model.eval()
        result = []
        with torch.no_grad():
            for start in range(0, len(rows), batch_size):
                batch = rows[start:start + batch_size]
                x1 = self._tensor(inputs.control[batch], model)
                features = self._tensor(inputs.fingerprints(batch), model)
                result.append(model(x1, features)[3].cpu().numpy())
        return np.concatenate(result) if result else np.empty((0, self.n_genes), np.float32)

    def train_epoch(self, model, inputs, rows, optimizer, rng):
        inputs = DoseView(inputs, self.transform)
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


WIDTH = WIDTHS["KPGT"]
TRANSFORM_KIND = "training-log10-dose-population-standardization-v1"


def positive_doses(doses):
    values = np.asarray(doses, dtype=np.float64)
    if values.ndim != 1 or not len(values) or not np.isfinite(values).all() or (values <= 0).any():
        raise ValueError("doses must be finite positive one-dimensional values")
    return values


def fit_transform(doses):
    logs = np.log10(positive_doses(doses))
    scale = float(logs.std(ddof=0))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("training dose scale must be positive")
    return dict(kind=TRANSFORM_KIND, mean=float(logs.mean()), std=scale)


def transform_doses(doses, transform):
    values = positive_doses(doses)
    if (not isinstance(transform, dict) or set(transform) != {"kind", "mean", "std"}
            or transform["kind"] != TRANSFORM_KIND):
        raise ValueError("invalid dose transform identity")
    mean, scale = transform["mean"], transform["std"]
    if not np.isfinite([mean, scale]).all() or scale <= 0:
        raise ValueError("invalid dose transform values")
    result = np.asarray((np.log10(values) - mean) / scale, dtype=np.float32)
    if not np.isfinite(result).all():
        raise ValueError("nonfinite standardized dose")
    return result


class DoseView:
    """Expose matched model inputs with one standardized dose per KPGT vector."""

    def __init__(self, data, transform):
        self.data = data
        self.transform = transform
        self.control = data.control
        self.expression = data.expression

    def fingerprints(self, rows):
        scalar = transform_doses(self.data.doses(rows), self.transform)
        return np.concatenate([self.data.fingerprints(rows), scalar[:, None]], axis=1).astype(np.float32)


def configuration(transform):
    transform_doses([1.], transform)
    order = ["native constructor and selective Xavier",
             "replace both predictive heads Linear532→100 with PyTorch defaults"]
    return dict(
        method="latent-numeric", method_identity="tsg-tenv-kpgt-latent-numeric-f32-v1",
        compute_dtype="float32", initialization_dtype="float32", prediction_dtype="float32",
        noise_dtype="native float32 CPU", batch_size=64, n_genes=5000, n_latent=100,
        encoder=[1200], decoder=[800], features_dim=WIDTH, features_embed_dim=[400],
        init_w=True, beta=0.1, dropout=0.1, molecule_feature="KPGT", initialization_model="random",
        learning_rate=1e-4, weight_decay=1e-5, scheduler=None, dose_input=True,
        prediction="native forward(x1,features)[3]", singleton_train="skip singleton as upstream",
        upstream_commit="e7582e7794f13d4dfddfb15e95e3ba7cd24bc07a",
        config_script="src/train_TranSiGen.py", config_script_sha256=SCRIPT_HASH,
        singleton_train_batch="skip singleton as upstream",
        adaptations=["5000 gene interface", "shared splits and fixed validation treated-expression MSE",
                     "x1-only native sampled forward for inference", "native float32 computation and CPU normal sampling",
                     "initial Adam learning rate 1e-4", "latent-numeric dose conditioning", "frozen KPGT2304 input",
                     "random response-model initialization"],
        transform=dict(transform), dose_vector_dim=32, predictive_head_dim=532, dose_embedding=None,
        constructor_order=order, initialization="; ".join(order), representation_width=WIDTH,
    )


def adapter(source, transform):
    configuration(transform)
    return Adapter(source, transform)
