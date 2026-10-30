#!/usr/bin/env python
"""PerturbNet adapter (ticket 14).

Author recipe (`pert_model/PerturbNet/source`, commit e1375b2c), as executed by
`notebooks/Tutorial_PerturbNet_Chemicals.ipynb`:

1. **ChemicalVAE** (`chemicalvae/chemicalVAE.py`): SMILES -> 120x35 one-hot ->
   196-d latent.  ZINC pretraining in the paper; here it is trained from
   scratch on the fold's training-drug representative SMILES only (525 epochs,
   batch 128, Adam lr 1e-4).  The 196-d latent is z-scored with mu/std
   recomputed on those training structures (`Standardize` semantics).
2. **scVI cell VAE** (scvi-tools 0.7.1 `scvi.model.SCVI`, ZINB, 700 epochs,
   `n_latent=10`), fitted on the fold's training treated wells (raw counts).
3. **cINN** (`cinn/flow.py`, `ConditionalFlatCouplingFlow` wrapped by
   `Net2NetFlow_scVIChemStdStatesFlow`): condition = standardized ChemicalVAE
   latent (196) + 3-bin dose one-hot (ticket 09), target = scVI latent (10).
   100 epochs, batch 128, Adam lr 4.5e-6 (betas 0.5/0.9, amsgrad).

Prediction follows `cinn/flow_generate.SCVIZ_CheckNet2Net`: sample the flow
conditioned on (drug, dose), decode the 10-d latent through the scVI decoder.
The benchmark needs one deterministic vector per test well, so the adapter
averages `flow_draws` flow samples and uses the ZINB rate (`px_rate`), i.e. the
per-draw predictive mean, instead of sampled counts; counts are expressed as
`log1p_cp10k` at the median training library size (the rate's cp10k form is
library-size invariant) and the matched plate DMSO mean is subtracted.

The TensorFlow 1.15 cell-VAE branch is not used: counts take the scVI/PyTorch
path (`cinn/flow.py:3289-3356`).  TF 1.15 still gets imported transitively by
`perturbnet.util`, and on this host its CUDA runtime segfaults, so the adapter
installs a minimal `tensorflow` stub before importing the vendored package.
Every deviation is recorded in `models/perturbnet/FIDELITY.md`.
"""

from __future__ import annotations

import os
import sys
import time
import types
from pathlib import Path

import numpy as np

# cuBLAS determinism for the 10-d linear algebra (ticket 10); must be set before
# the first CUDA call.
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")

# Model code is vendored read-only under `pert_model/`; never leave __pycache__.
sys.dont_write_bytecode = True

BENCH = Path(__file__).resolve().parents[3]
ROOT = BENCH.parents[2]
sys.path.insert(0, str(BENCH / "scripts"))
sys.path.insert(0, str(ROOT / "pert_model/PerturbNet/source"))


def _install_tensorflow_stub() -> None:
    """Satisfy `perturbnet.util`'s wildcard import of the TF cell-VAE module.

    `perturbnet/util.py:40` does ``from perturbnet.data_vae.vae import *``, and
    that module imports TensorFlow at module level even though the counts/scVI
    path never calls it.  TensorFlow 1.15 is installed in the faithful env but
    its CUDA 10 runtime cannot initialise next to torch 1.13/cu117 on this box
    (segfault), so the unused module is stubbed rather than imported.
    """

    if "tensorflow" in sys.modules:
        return
    tf = types.ModuleType("tensorflow")
    tf.distributions = types.ModuleType("tensorflow.distributions")
    tf.__version__ = "0.0.0-stub"
    sys.modules["tensorflow"] = tf
    sys.modules["tensorflow.distributions"] = tf.distributions


_install_tensorflow_stub()

import anndata as ad  # noqa: E402
import scvi  # noqa: E402
import torch  # noqa: E402
from rdkit import Chem  # noqa: E402

from perturbnet.chemicalvae.chemicalVAE import ChemicalVAE, vae_loss  # noqa: E402
from perturbnet.cinn.flow import (  # noqa: E402
    ConditionalFlatCouplingFlow,
    Net2NetFlow_scVIChemStdStatesFlow,
)
from perturbnet.pytorch_scvi.scvi_generate_z import scvi_predictive_z  # noqa: E402
from perturbnet.util import StandardizeLoad, smiles_to_hot  # noqa: E402

from adapter_api import Adapter, run_cli  # noqa: E402
from adapters.common import (  # noqa: E402
    N_GENES,
    Controls,
    case_summary,
    ensure_float32,
    load_json,
    log,
    write_json,
)

SMILES_ALPHABET = [
    "7", "6", "o", "]", "3", "s", "(", "-", "S", "/",
    "B", "4", "[", ")", "#", "I", "l", "O", "H", "c",
    "1", "@", "=", "n", "P", "8", "C", "2", "F", "5",
    "r", "N", "+", "\\", " ",
]
_ALPHABET = set(SMILES_ALPHABET)
MAX_LEN = 120
N_CHAR = 35
Z_DIM = 196


def _smoke() -> bool:
    return os.environ.get("BENCH_SMOKE", "").strip() not in ("", "0", "false", "False")


def _save_atomic(payload, path: Path) -> None:
    """Write a checkpoint through a temporary file, then rename it."""

    path = Path(path)
    tmp = Path(str(path) + ".tmp")
    torch.save(payload, tmp)
    os.replace(tmp, path)


def _seed_everything(seed: int) -> None:
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    try:
        torch.use_deterministic_algorithms(True, warn_only=True)
    except Exception:  # older torch builds without the kwarg
        pass


# ---------------------------------------------------------------------------
# Chemical structure


def sanitize_smiles(value: str) -> str | None:
    """Largest-fragment canonical SMILES (the benchmark's parent form).

    The corpus' representative SMILES still carry salts for 25/460 keys
    (`.Cl`, `.[Na+]`, ...).  The ChemicalVAE alphabet has no `.` or `a`, and
    the author's `smiles_to_hot` raises on unknown characters, so the salt
    counter-ions are dropped.  Returns ``None`` when RDKit cannot parse.
    """

    text = str(value).strip()
    if not text:
        return None
    molecule = Chem.MolFromSmiles(text)
    if molecule is None:
        return None
    fragments = Chem.GetMolFrags(molecule, asMols=True)
    if len(fragments) > 1:
        molecule = max(fragments, key=lambda mol: mol.GetNumHeavyAtoms())
    return Chem.MolToSmiles(molecule)


def encodable(smiles: str | None) -> bool:
    if not smiles or len(smiles) > MAX_LEN:
        return False
    return not (set(smiles) - _ALPHABET)


def onehot_smiles(smiles: list[str]) -> np.ndarray:
    """120x35 one-hot for already-sanitized, encodable SMILES."""

    return smiles_to_hot(smiles, max_len=MAX_LEN, padding="right", nchars=N_CHAR)


def load_case(path: Path, keys: tuple[str, ...]) -> dict[str, np.ndarray]:
    """Read only the bundle keys the adapter needs (peak-RSS discipline)."""

    with np.load(Path(path), allow_pickle=False) as archive:
        return {key: archive[key] for key in keys}


class PerturbNetAdapter(Adapter):
    name = "perturbnet"

    def __init__(self) -> None:
        self.cfg = load_json(BENCH / "models/perturbnet/config.json")

    # -- shared helpers -------------------------------------------------

    @staticmethod
    def _device() -> torch.device:
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")

    @staticmethod
    def _resolve(config: dict, cfg: dict, smoke: bool) -> dict:
        """Author config, with the `smoke` block folded into each stage."""

        merged = dict(cfg)
        merged["chemical_vae"] = dict(cfg["chemical_vae"])
        merged["cell_vae"] = dict(cfg["cell_vae"])
        merged["cinn"] = dict(cfg["cinn"])
        if smoke:
            block = cfg["smoke"]
            merged["chemical_vae"]["epochs"] = int(block["chemical_vae_epochs"])
            merged["cell_vae"]["epochs"] = int(block["cell_vae_epochs"])
            merged["cell_vae"]["frequency"] = int(block["cell_vae_frequency"])
            merged["cinn"]["epochs"] = int(block["cinn_epochs"])
            merged["chemical_vae_max_batches"] = int(block["chemical_vae_max_batches"])
            merged["cinn_max_batches"] = int(block["cinn_max_batches"])
            merged["valid_eval_max_rows"] = int(block["valid_eval_max_rows"])
            merged["flow_draws"] = int(block["flow_draws"])
        merged["seed"] = int(config.get("seed", cfg["seed"]))
        return merged

    @staticmethod
    def _drug_table(case: dict[str, np.ndarray]) -> dict[str, str]:
        table: dict[str, str] = {}
        for key, smiles in zip(case["drug_key"], case["smiles"]):
            table.setdefault(str(key), str(smiles))
        return table

    def _encode_chemical(self, drug_table: dict[str, str], chemvae, std_model):
        """Sanitized SMILES -> standardized 196-d latent (unencodable -> 0)."""

        keys = list(drug_table)
        clean = {key: sanitize_smiles(drug_table[key]) for key in keys}
        usable = [key for key in keys if encodable(clean[key])]
        latents = np.zeros((len(keys), Z_DIM), dtype=np.float32)
        if usable:
            onehot = onehot_smiles([clean[key] for key in usable])
            with torch.no_grad():
                _, _, _, sampled = chemvae(torch.from_numpy(onehot).to(std_model.device))
            latents[[keys.index(key) for key in usable]] = std_model.standardize_z(
                sampled.cpu().numpy()
            )
        index = {key: i for i, key in enumerate(keys)}
        dropped = sorted(set(keys) - set(usable))
        return latents, index, dropped

    @staticmethod
    def _dose_bins(doses: np.ndarray) -> list[float]:
        return [float(x) for x in np.unique(np.asarray(doses, dtype=np.float64))]

    @staticmethod
    def _dose_onehot(doses: np.ndarray, bins: list[float]) -> np.ndarray:
        values = np.asarray(doses, dtype=np.float64)
        out = np.zeros((len(values), len(bins)), dtype=np.float32)
        for row, value in enumerate(values):
            distances = np.abs(np.asarray(bins) - value)
            out[row, int(np.argmin(distances))] = 1.0
        return out

    def _condition(
        self,
        case: dict[str, np.ndarray],
        latents: np.ndarray,
        index: dict[str, int],
        bins: list[float],
    ) -> np.ndarray:
        rows = np.stack([latents[index[str(key)]] for key in case["drug_key"]])
        dose = self._dose_onehot(case["dose_um"], bins)
        return np.concatenate([rows, dose], axis=1).astype(np.float32)

    # -- stage 1: ChemicalVAE -------------------------------------------

    def _train_chemical_vae(self, structures, cfg, device, seed, max_batches):
        chemvae = ChemicalVAE(n_char=N_CHAR, max_len=MAX_LEN, z_dim=Z_DIM).to(device)
        optimizer = torch.optim.Adam(chemvae.parameters(), lr=float(cfg["lr"]))
        batch_size = int(cfg["batch_size"])
        epochs = int(cfg["epochs"])
        onehot = torch.from_numpy(onehot_smiles(structures))

        n = len(structures)
        n_train = int(n * float(cfg["train_ratio"]))
        rng = np.random.RandomState(seed)
        permutation = rng.permutation(n)
        train_idx = permutation[:n_train]
        valid_idx = permutation[n_train:]
        history: list[dict] = []
        started = time.time()

        for epoch in range(epochs):
            chemvae.train()
            order = rng.permutation(train_idx)
            total, seen = 0.0, 0
            for step, start in enumerate(range(0, len(order), batch_size)):
                if step >= max_batches:
                    break
                batch = order[start : start + batch_size]
                data = onehot[batch].to(device, dtype=torch.float)
                output, mean, logvar, _ = chemvae(data)
                loss = vae_loss(output, data, mean, logvar)
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                optimizer.step()
                total += float(loss.item())
                seen += len(batch)

            chemvae.eval()
            valid_total, valid_seen = 0.0, 0
            with torch.no_grad():
                for start in range(0, len(valid_idx), batch_size):
                    batch = valid_idx[start : start + batch_size]
                    data = onehot[batch].to(device, dtype=torch.float)
                    output, mean, logvar, _ = chemvae(data)
                    valid_total += float(vae_loss(output, data, mean, logvar).item())
                    valid_seen += len(batch)

            history.append(
                {
                    "epoch": epoch,
                    "train_loss": total / max(seen, 1),
                    "valid_loss": valid_total / max(valid_seen, 1),
                }
            )
            if epoch == 0 or (epoch + 1) % 50 == 0 or epoch + 1 == epochs:
                log(
                    f"chemvae epoch {epoch + 1}/{epochs} "
                    f"train={history[-1]['train_loss']:.3f} "
                    f"valid={history[-1]['valid_loss']:.3f}"
                )

        # Author `Standardize`: mu/std of the 196-d latent over the corpus.
        chemvae.eval()
        with torch.no_grad():
            chunks = []
            for start in range(0, n, 512):
                batch = onehot[start : start + 512].to(device, dtype=torch.float)
                _, _, _, sampled = chemvae(batch)
                chunks.append(sampled.cpu().numpy())
            latent = np.concatenate(chunks, axis=0)
        mu = latent.mean(axis=0)
        std = latent.std(axis=0)
        std = np.where(std < 1e-6, 1.0, std).astype(np.float32)
        return chemvae, mu.astype(np.float32), std, history, time.time() - started

    # -- stage 2: scVI ---------------------------------------------------

    def _train_cell_vae(self, counts, cfg, seed):
        torch.manual_seed(seed)
        model_cfg = cfg
        adata = ad.AnnData(X=np.asarray(counts, dtype=np.float32))
        scvi.data.setup_anndata(adata, layer=None)
        model = scvi.model.SCVI(
            adata,
            n_latent=int(model_cfg["n_latent"]),
            n_hidden=int(model_cfg["n_hidden"]),
            n_layers=int(model_cfg["n_layers"]),
            gene_likelihood=str(model_cfg["gene_likelihood"]),
            use_cuda=torch.cuda.is_available(),
        )
        started = time.time()
        model.train(
            n_epochs=int(cfg["epochs"]),
            frequency=int(cfg["frequency"]),
            lr=float(model_cfg.get("lr", 1e-3)),
        )
        return model, adata, time.time() - started

    # -- stage 3: cINN ---------------------------------------------------

    def _build_flow(self, cfg, device):
        flow_cfg = dict(cfg["flow"])
        return ConditionalFlatCouplingFlow(
            conditioning_dim=int(flow_cfg["conditioning_dim"]),
            embedding_dim=int(flow_cfg["embedding_dim"]),
            conditioning_depth=int(flow_cfg["conditioning_depth"]),
            n_flows=int(flow_cfg["n_flows"]),
            in_channels=int(flow_cfg["in_channels"]),
            hidden_dim=int(flow_cfg["hidden_dim"]),
            hidden_depth=int(flow_cfg["hidden_depth"]),
            activation=str(flow_cfg["activation"]),
            conditioner_use_bn=bool(flow_cfg["conditioner_use_bn"]),
        ).to(device)

    @staticmethod
    def _nll(wrapper, z_np, condition, device) -> float:
        z = torch.from_numpy(np.asarray(z_np, dtype=np.float32)).to(device)
        cond = torch.from_numpy(np.asarray(condition, dtype=np.float32)).to(device)
        loss, _ = wrapper.shared_step(z.unsqueeze(-1).unsqueeze(-1), cond, 0, split="valid")
        return float(loss.item())

    # -- Adapter API -----------------------------------------------------

    def fit(self, train, valid, config, run_dir, controls=None):
        if controls is None:
            raise RuntimeError("PerturbNet needs --controls (the shared DMSO means)")
        cfg = self._resolve(config, self.cfg, _smoke())
        seed = int(cfg["seed"])
        _seed_everything(seed)
        device = self._device()
        run_dir = Path(run_dir)
        checkpoint = run_dir / "checkpoint"
        checkpoint.mkdir(parents=True, exist_ok=True)

        train_case = load_case(
            train, ("X", "smiles", "dose_um", "drug_key", "well_id", "group_id", "control_group")
        )
        valid_case = load_case(
            valid, ("X", "smiles", "dose_um", "drug_key", "well_id", "group_id", "control_group")
        )
        log(f"device={device} train={case_summary(train_case)} valid={case_summary(valid_case)}")

        # ---- stage 1: ChemicalVAE on the training drugs ----------------
        drug_table = self._drug_table(train_case)
        clean = {key: sanitize_smiles(smiles) for key, smiles in drug_table.items()}
        structures = [clean[key] for key in sorted(drug_table) if encodable(clean[key])]
        dropped_train = sorted(key for key in drug_table if not encodable(clean[key]))
        if dropped_train:
            log(f"chemvae: {len(dropped_train)} training drugs are not encodable: {dropped_train}")
        log(f"chemvae: {len(structures)} unique training structures")
        chemvae, mu, std, chem_history, chem_wall = self._train_chemical_vae(
            structures,
            cfg["chemical_vae"],
            device,
            seed,
            int(cfg.get("chemical_vae_max_batches", 10**9)),
        )
        std_model = StandardizeLoad(mu, std, device)
        train_latents, train_index, _ = self._encode_chemical(drug_table, chemvae, std_model)
        _save_atomic(
            {
                "state_dict": chemvae.state_dict(),
                "mu": mu,
                "std": std,
                "structures": structures,
                "n_char": N_CHAR,
                "max_len": MAX_LEN,
                "z_dim": Z_DIM,
            },
            checkpoint / "chemicalvae.pt",
        )
        log(f"chemvae done in {chem_wall:.1f}s ({len(structures)} structures, {int(cfg['chemical_vae']['epochs'])} epochs)")

        # ---- stage 2: scVI on training counts --------------------------
        cell_model, adata_train, cell_wall = self._train_cell_vae(
            train_case["X"], cfg["cell_vae"], seed
        )
        log(f"scvi done in {cell_wall:.1f}s ({int(cfg['cell_vae']['epochs'])} epochs)")
        median_log_library = float(np.median(np.log(train_case["X"].sum(axis=1))))
        _save_atomic(
            {
                "state_dict": cell_model.model.state_dict(),
                "n_input": int(train_case["X"].shape[1]),
                "n_batch": int(cell_model.model.n_batch),
                "n_hidden": int(cfg["cell_vae"]["n_hidden"]),
                "n_latent": int(cfg["cell_vae"]["n_latent"]),
                "n_layers": int(cfg["cell_vae"]["n_layers"]),
                "gene_likelihood": str(cfg["cell_vae"]["gene_likelihood"]),
                "median_log_library": median_log_library,
            },
            checkpoint / "scvi_vae.pt",
        )

        # ---- stage 3: cINN ---------------------------------------------
        bins = self._dose_bins(train_case["dose_um"])
        if len(bins) != 3:
            raise RuntimeError(f"expected 3 dose bins, found {bins}")
        # Early metadata makes a partially trained checkpoint predictable:
        # predict only needs the seed, dose bins, draw count, and library.
        write_json(
            checkpoint / "meta.json",
            {
                "model": "perturbnet",
                "n_genes": N_GENES,
                "seed": seed,
                "dose_bins": bins,
                "flow_draws": int(cfg.get("flow_draws", 8)),
                "median_log_library": median_log_library,
                "partial": True,
                "smoke": _smoke(),
            },
        )
        _seed_everything(seed)  # cINN init independent of the scVI stage's RNG use
        flow = self._build_flow(cfg["cinn"], device)
        _save_atomic(
            {"state_dict": flow.state_dict(), "config": cfg["cinn"]["flow"]},
            checkpoint / "cinn.pt",
        )
        wrapper = Net2NetFlow_scVIChemStdStatesFlow(
            configured_flow=flow,
            first_stage_data=adata_train.X,
            cond_stage_data=None,
            model_con=chemvae,
            scvi_model=cell_model,
            std_model=std_model,
            cell_type=np.zeros((len(train_case["X"]), 0), dtype=np.float32),
            dose=self._dose_onehot(train_case["dose_um"], bins),
        )
        optimizer = wrapper.configure_flow_optimizers(lr=float(cfg["cinn"]["lr"]))
        cond_train = self._condition(train_case, train_latents, train_index, bins)
        valid_table = self._drug_table(valid_case)
        valid_latents, valid_index, dropped_valid = self._encode_chemical(
            valid_table, chemvae, std_model
        )
        if dropped_valid:
            log(f"chemvae: {len(dropped_valid)} valid drugs use the mean latent: {dropped_valid}")
        cond_valid = self._condition(valid_case, valid_latents, valid_index, bins)
        valid_adata = ad.AnnData(X=np.asarray(valid_case["X"], dtype=np.float32))
        scvi.data.setup_anndata(valid_adata, layer=None)

        batch_size = int(cfg["cinn"]["batch_size"])
        epochs = int(cfg["cinn"]["epochs"])
        max_batches = int(cfg.get("cinn_max_batches", 10**9))
        eval_rows = int(cfg.get("valid_eval_max_rows", len(cond_valid)))
        eval_rows = min(eval_rows, len(cond_valid))
        rng = np.random.RandomState(seed + 1)
        cinn_history: list[dict] = []
        best_valid = float("inf")
        best_state = None
        best_epoch = -1
        cinn_started = time.time()

        for epoch in range(epochs):
            flow.train()
            order = rng.permutation(len(cond_train))
            total, seen = 0.0, 0
            for step, start in enumerate(range(0, len(order), batch_size)):
                if step >= max_batches:
                    break
                batch = np.sort(order[start : start + batch_size])
                latents = cell_model.get_latent_representation(
                    indices=batch, give_mean=False
                )
                condition = torch.from_numpy(cond_train[batch]).to(device)
                target = (
                    torch.from_numpy(np.asarray(latents, dtype=np.float32))
                    .to(device)
                    .unsqueeze(-1)
                    .unsqueeze(-1)
                )
                loss, _ = wrapper.shared_step(target, condition, step)
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                optimizer.step()
                total += float(loss.item()) * len(batch)
                seen += len(batch)

            flow.eval()
            if eval_rows:
                latents_valid = cell_model.get_latent_representation(
                    adata=valid_adata, give_mean=True
                )[:eval_rows]
                valid_loss = self._nll(
                    wrapper, latents_valid, cond_valid[:eval_rows], device
                )
            else:
                valid_loss = float("nan")
            cinn_history.append(
                {
                    "epoch": epoch,
                    "train_loss": total / max(seen, 1),
                    "valid_loss": valid_loss,
                }
            )
            improved = valid_loss < best_valid
            if improved:
                best_valid = valid_loss
                best_epoch = epoch
                best_state = {
                    key: value.detach().cpu().clone()
                    for key, value in flow.state_dict().items()
                }
                if (epoch + 1) % 10 == 0 or epoch + 1 == epochs:
                    _save_atomic(
                        {"state_dict": best_state, "config": cfg["cinn"]["flow"]},
                        checkpoint / "cinn.pt",
                    )
            if epoch == 0 or (epoch + 1) % 10 == 0 or epoch + 1 == epochs:
                log(
                    f"cinn epoch {epoch + 1}/{epochs} train={cinn_history[-1]['train_loss']:.4f} "
                    f"valid={valid_loss:.4f} best={best_valid:.4f}@{best_epoch + 1}"
                )
        cinn_wall = time.time() - cinn_started

        if best_state is not None:
            flow.load_state_dict(best_state)

        # ---- checkpoints ----------------------------------------------
        _save_atomic(
            {
                "state_dict": chemvae.state_dict(),
                "mu": mu,
                "std": std,
                "structures": structures,
                "n_char": N_CHAR,
                "max_len": MAX_LEN,
                "z_dim": Z_DIM,
            },
            checkpoint / "chemicalvae.pt",
        )
        _save_atomic(
            {
                "state_dict": cell_model.model.state_dict(),
                "n_input": int(train_case["X"].shape[1]),
                "n_batch": int(cell_model.model.n_batch),
                "n_hidden": int(cfg["cell_vae"]["n_hidden"]),
                "n_latent": int(cfg["cell_vae"]["n_latent"]),
                "n_layers": int(cfg["cell_vae"]["n_layers"]),
                "gene_likelihood": str(cfg["cell_vae"]["gene_likelihood"]),
                "median_log_library": median_log_library,
            },
            checkpoint / "scvi_vae.pt",
        )
        _save_atomic(
            {"state_dict": flow.state_dict(), "config": cfg["cinn"]["flow"]},
            checkpoint / "cinn.pt",
        )
        meta = {
            "model": "perturbnet",
            "n_genes": N_GENES,
            "seed": seed,
            "dose_bins": bins,
            "flow_draws": int(cfg.get("flow_draws", 8)),
            "median_log_library": median_log_library,
            "best_valid_loss": best_valid,
            "best_epoch": best_epoch,
            "epochs": {"chemical_vae": int(cfg["chemical_vae"]["epochs"]),
                       "cell_vae": int(cfg["cell_vae"]["epochs"]),
                       "cinn": epochs},
            "wall_s": {"chemical_vae": chem_wall, "cell_vae": cell_wall, "cinn": cinn_wall},
            "dropped_train_drugs": dropped_train,
            "dropped_valid_drugs": dropped_valid,
            "smoke": _smoke(),
        }
        write_json(checkpoint / "meta.json", meta)
        write_json(
            run_dir / "fit_history.json",
            {
                "chemical_vae": chem_history,
                "cinn": cinn_history,
                "best_valid_loss": best_valid,
                "best_epoch": best_epoch,
                "smoke": _smoke(),
            },
        )
        size_mb = sum(p.stat().st_size for p in checkpoint.glob("*")) / 1e6
        log(
            f"checkpoint {checkpoint} ({size_mb:.1f} MB) best valid {best_valid:.4f} "
            f"@ epoch {best_epoch + 1}"
        )
        return checkpoint

    def predict(self, cases, checkpoint, run_dir, controls=None):
        if controls is None:
            raise RuntimeError("PerturbNet needs --controls (the shared DMSO means)")
        checkpoint = Path(checkpoint)
        meta = load_json(checkpoint / "meta.json")
        seed = int(meta["seed"])
        _seed_everything(seed)
        device = self._device()

        stored = torch.load(checkpoint / "chemicalvae.pt", map_location="cpu")
        chemvae = ChemicalVAE(n_char=N_CHAR, max_len=MAX_LEN, z_dim=Z_DIM)
        chemvae.load_state_dict(stored["state_dict"])
        chemvae.to(device).eval()
        std_model = StandardizeLoad(stored["mu"], stored["std"], device)

        cell = torch.load(checkpoint / "scvi_vae.pt", map_location="cpu")
        from scvi.core.modules import VAE as ScviVAE

        vae = ScviVAE(
            n_input=int(cell["n_input"]),
            n_batch=int(cell["n_batch"]),
            n_hidden=int(cell["n_hidden"]),
            n_latent=int(cell["n_latent"]),
            n_layers=int(cell["n_layers"]),
            gene_likelihood=str(cell["gene_likelihood"]),
        )
        vae.load_state_dict(cell["state_dict"])
        vae.to(device).eval()
        decoder = scvi_predictive_z(types.SimpleNamespace(model=vae))

        cinn = torch.load(checkpoint / "cinn.pt", map_location="cpu")
        flow = self._build_flow({"flow": cinn["config"]}, device)
        flow.load_state_dict(cinn["state_dict"])
        flow.eval()

        case = load_case(
            cases, ("smiles", "dose_um", "drug_key", "control_group", "well_id")
        )
        controls_obj = Controls(controls)
        table = self._drug_table(case)
        latents, index, dropped = self._encode_chemical(table, chemvae, std_model)
        if dropped:
            log(f"predict: {len(dropped)} unencodable drugs use the mean latent: {dropped}")
        condition = self._condition(case, latents, index, list(meta["dose_bins"]))
        control_mean = controls_obj.mean_for(case["control_group"])

        draws = int(meta.get("flow_draws", 8))
        batch_size = int(self.cfg["cinn"]["batch_size"])
        library = float(meta["median_log_library"])
        n = len(condition)
        cp10k = np.zeros((n, N_GENES), dtype=np.float32)
        with torch.no_grad():
            for draw in range(draws):
                for start in range(0, n, batch_size):
                    stop = min(start + batch_size, n)
                    cond = torch.from_numpy(condition[start:stop]).to(device)
                    z = flow.sample(cond).squeeze(-1).squeeze(-1)
                    lib = torch.full((stop - start, 1), library, device=device)
                    outputs = decoder.decoder_inference(
                        z,
                        lib,
                        batch_index=torch.zeros(stop - start, dtype=torch.long, device=device),
                    )
                    # px_rate is counts at exp(library); renormalise to cp10k.
                    counts = outputs["px_rate"] / float(np.exp(library))
                    cp10k[start:stop] += (counts * 1e4).cpu().numpy().astype(np.float32)
        cp10k /= float(draws)
        delta = np.log1p(cp10k) - control_mean
        log(
            f"predict {delta.shape} draws={draws} finite={np.isfinite(delta).all()} "
            f"mean={float(delta.mean()):.5f}"
        )
        return ensure_float32(delta)


if __name__ == "__main__":
    raise SystemExit(run_cli(PerturbNetAdapter()))
