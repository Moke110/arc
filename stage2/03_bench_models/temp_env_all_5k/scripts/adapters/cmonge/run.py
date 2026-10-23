#!/usr/bin/env python
"""CMonge adapter (ticket 15).

Author recipe (`.scratch/stage2-03-bench-models/research/03-squidiff-cmonge-decant.md`
§2, ticket 09):

- **Stage 1, autoencoder** (`configs/autoencoder-sciplex.yml`): 2x512 hidden,
  50-d latent, ``data_dim`` = gene count, GELU, MSE
  (``2 * optax.l2_loss(pred, batch).mean()``), AdamW lr 1e-4 / wd 1e-5,
  warmup-cosine 0 -> 1e-4 -> 1e-5 (warmup 100 steps, decay
  ``n_epochs * 100``), batch 512, 50 author "epochs" (one batch each,
  ``ae_trainer.py:149-152``).  Trained on controls + treated cells, then
  **frozen**.
- **Stage 2, conditional Monge map** (`configs/conditional-monge-sciplex.yml`
  with the demo config's context wiring): 4x64 MLP, ``dim_data`` 50,
  ``dim_cond_map`` 50, Sinkhorn fitting loss ``epsilon_fitting=1`` + Monge-gap
  regularizer ``epsilon_regularizer=1e-2`` (Euclidean), AdamW lr 1e-4 /
  wd 1e-5 with cosine decay (alpha 1e-2), batch 512, one condition sampled
  uniformly per iteration, **10,000 iterations** (the sciplex setting, not the
  500k all-drug setting).
- **Drug condition** (``models/embedding.py:39-144``): RDKit 2D descriptors
  (``cmonge.models.rdkit.RDKIT_PROPS``), first column dropped, zero-variance
  columns dropped, z-scored, then ``log(dose)`` appended (author's dose strings
  are integer nM).  The context is split into drug (desca) and dose (whole
  context) chunks by ``context_entity_bonds``.

Benchmark choices (all recorded in ``models/cmonge/FIDELITY.md``): the AE shows
the fold's training treated wells + all 192 shared DMSO cells; the descriptor
z-scoring is fitted on the fold's **training drug panel only** (the author fits
on the whole panel, mildly transductive); checkpoint selection is the author
loss on the cold valid drugs every 1,000 iterations; prediction maps one
matched DMSO cell per test well under the test drug+dose condition, decodes to
5,000 genes and subtracts the plate DMSO mean.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

# Model code is vendored read-only under `pert_model/`; never leave __pycache__.
sys.dont_write_bytecode = True
os.environ.setdefault("PYTHONDONTWRITEBYTECODE", "1")

BENCH = Path(__file__).resolve().parents[3]
ROOT = BENCH.parents[2]
sys.path.insert(0, str(BENCH / "scripts"))
sys.path.insert(0, str(ROOT / "pert_model/cmonge"))

from adapter_api import Adapter, run_cli  # noqa: E402
from adapters.common import (  # noqa: E402
    N_GENES,
    Controls,
    case_summary,
    ensure_float32,
    load_json,
    load_npz,
    log,
    write_json,
)

import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402
import optax  # noqa: E402
from flax import linen as nn  # noqa: E402
from flax import serialization  # noqa: E402
from flax.training.train_state import TrainState  # noqa: E402

from cmonge.metrics import fitting_loss, regularizer  # noqa: E402
from cmonge.models.nn import ConditionalPerturbationNetwork  # noqa: E402
from cmonge.models.rdkit import rdkit_feats  # noqa: E402
from cmonge.trainers.ae_trainer import AutoEncoder, mse_reconstruction_loss  # noqa: E402

PREDICT_BATCH = 128


def _smoke() -> bool:
    return os.environ.get("BENCH_SMOKE", "").strip() not in ("", "0", "false", "False")


def _activation(name: str):
    return {
        "gelu": nn.gelu,
        "relu": nn.relu,
        "leakyrelu": nn.leaky_relu,
    }[str(name)]


# ---------------------------------------------------------------------------
# RDKit descriptors and the author's standardisation


def _rdkit_rows(smiles) -> tuple[np.ndarray, int]:
    """RDKit2D rows for every SMILES, author transform without the flag column.

    ``cmonge.models.rdkit.rdkit_feats`` returns ``[True] + 200 values``; the
    author drops ``latent_0`` (the flag) and then variance-filters the panel.
    Non-finite descriptors become 0 exactly as ``embedding.py:78-86`` does.
    """

    unique = sorted({str(s) for s in smiles if str(s)})
    cache: dict[str, np.ndarray] = {}
    width: int | None = None
    for value in unique:
        features = rdkit_feats(value)
        if not features or len(features) < 2:
            if width is None:
                continue
            cache[value] = np.zeros(width, dtype=np.float64)
            continue
        vector = np.asarray(
            [0.0 if f is None else float(f) for f in features[1:]], dtype=np.float64
        )
        vector[~np.isfinite(vector)] = 0.0
        if width is None:
            width = int(vector.shape[0])
        cache[value] = vector
    if width is None:
        raise RuntimeError("RDKit2D produced no descriptor vector for any SMILES")
    rows = np.zeros((len(smiles), width), dtype=np.float64)
    for i, value in enumerate(smiles):
        rows[i] = cache.get(str(value), np.zeros(width, dtype=np.float64))
    return rows, width


def _fit_descriptor_standardizer(matrix: np.ndarray, threshold: float) -> dict:
    """Author ``embedding.py:99-113``: drop low-std columns, then z-score.

    ``pandas.DataFrame.std`` is ddof=1, so the standard deviation here is ddof=1.
    """

    mean = matrix.mean(axis=0)
    std = matrix.std(axis=0, ddof=1)
    keep = std > float(threshold)
    safe = np.where(keep, std, 1.0)
    return {
        "mean": mean.astype(np.float64),
        "std": safe.astype(np.float64),
        "keep": keep,
    }


def _apply_descriptor_standardizer(matrix: np.ndarray, stats: dict) -> np.ndarray:
    keep = np.asarray(stats["keep"], dtype=bool)
    mean = np.asarray(stats["mean"], dtype=np.float64)[keep]
    std = np.asarray(stats["std"], dtype=np.float64)[keep]
    return ((matrix[:, keep] - mean) / std).astype(np.float32)


def _drug_smiles(case: dict[str, np.ndarray]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for drug, smiles in zip(case["drug_key"], case["smiles"]):
        mapping.setdefault(str(drug), str(smiles))
    return mapping


def _condition_vectors(
    descriptor_rows: np.ndarray,
    order: list[str],
    drugs: list[str],
    doses: np.ndarray,
    stats: dict,
) -> np.ndarray:
    """``[z-scored descriptors, log(dose_nM)]`` per (drug, dose) condition.

    ``descriptor_rows`` is aligned with ``order`` (the sorted drug keys); the
    dose is the author's integer-nM ``log(dose)`` (our doses are µM).
    """

    normalized = _apply_descriptor_standardizer(descriptor_rows, stats)
    by_drug = {name: normalized[i] for i, name in enumerate(order)}
    dose_log = np.log(np.asarray(doses, dtype=np.float64) * 1000.0)
    descriptors = np.stack([by_drug[str(d)] for d in drugs]).astype(np.float32)
    return np.concatenate(
        [descriptors, dose_log.astype(np.float32)[:, None]], axis=1
    )


# ---------------------------------------------------------------------------
# Adapter


class CMongeAdapter(Adapter):
    name = "cmonge"

    def __init__(self) -> None:
        self.cfg = load_json(BENCH / "models/cmonge/config.json")

    # -- settings ------------------------------------------------------

    def _settings(self) -> dict:
        cfg = self.cfg
        settings = {
            "ae_epochs": int(cfg["ae"]["n_epochs"]),
            "ae_batch": int(cfg["ae"]["batch_size"]),
            "cm_iters": int(cfg["monge"]["num_train_iters"]),
            "cm_batch": int(cfg["monge"]["batch_size"]),
            "eval_every": int(cfg["monge"]["eval_every"]),
        }
        if _smoke():
            settings.update({k: int(v) for k, v in cfg["smoke"].items()})
        return settings

    # -- AE ------------------------------------------------------------

    def _ae_module(self, seed: int) -> AutoEncoder:
        ae = self.cfg["ae"]
        return AutoEncoder(
            hidden_dims=list(ae["hidden_dims"]),
            latent_dim=int(ae["latent_dim"]),
            data_dim=int(ae["data_dim"]),
            seed=int(seed),
            act_fn=_activation(ae["act_fn"]),
        )

    def _train_ae(
        self,
        pool: np.ndarray,
        module: AutoEncoder,
        settings: dict,
        seed: int,
    ) -> tuple[dict, list[dict]]:
        ae = self.cfg["ae"]
        rng = jax.random.PRNGKey(seed)
        rng, init_rng = jax.random.split(rng)
        params = module.init(init_rng, jnp.empty((1, int(ae["data_dim"]))))["params"]
        schedule = optax.warmup_cosine_decay_schedule(
            init_value=0.0,
            peak_value=float(ae["optimizer"]["lr"]),
            warmup_steps=int(ae["schedule"]["warmup_steps"]),
            decay_steps=int(settings["ae_epochs"]) * int(ae["schedule"]["decay_steps_per_epoch"]),
            end_value=float(ae["schedule"]["end_value"]),
        )
        optimizer = optax.adamw(
            learning_rate=schedule,
            weight_decay=float(ae["optimizer"]["weight_decay"]),
        )
        state = TrainState.create(apply_fn=module.apply, params=params, tx=optimizer)

        @jax.jit
        def train_step(state: TrainState, batch: jnp.ndarray):
            def loss_fn(params):
                return mse_reconstruction_loss(module, params, batch)

            loss, grads = jax.value_and_grad(loss_fn)(state.params)
            return state.apply_gradients(grads=grads), loss

        batch_size = int(settings["ae_batch"])
        n = pool.shape[0]
        sampler = np.random.default_rng(seed + 17)
        order = sampler.permutation(n)
        cursor = 0
        history: list[dict] = []
        for epoch in range(int(settings["ae_epochs"])):
            if cursor + batch_size > n:
                order = sampler.permutation(n)
                cursor = 0
            batch = pool[order[cursor : cursor + batch_size]]
            cursor += batch_size
            state, loss = train_step(state, jnp.asarray(batch))
            history.append({"epoch": epoch, "loss": float(loss)})
            if epoch % 10 == 0 or epoch + 1 == int(settings["ae_epochs"]):
                log(f"ae epoch {epoch + 1}/{settings['ae_epochs']} loss={float(loss):.6f}")
        return state.params, history

    @staticmethod
    def _encode(params: dict, module: AutoEncoder, matrix: np.ndarray) -> np.ndarray:
        @jax.jit
        def run(params: dict, batch: jnp.ndarray):
            return module.apply(
                {"params": params}, batch, method=lambda m, x: m.encoder(x)
            )

        out = np.empty((matrix.shape[0], int(module.latent_dim)), dtype=np.float32)
        for start in range(0, matrix.shape[0], 256):
            stop = min(start + 256, matrix.shape[0])
            batch = np.ascontiguousarray(matrix[start:stop], dtype=np.float32)
            out[start:stop] = np.asarray(run(params, jnp.asarray(batch)), dtype=np.float32)
        return out

    @staticmethod
    def _decode(params: dict, module: AutoEncoder, latents: np.ndarray) -> np.ndarray:
        @jax.jit
        def run(params: dict, batch: jnp.ndarray):
            return module.apply(
                {"params": params}, batch, method=lambda m, z: m.decoder(z)
            )

        out = np.empty((latents.shape[0], int(module.data_dim)), dtype=np.float32)
        for start in range(0, latents.shape[0], PREDICT_BATCH):
            stop = min(start + PREDICT_BATCH, latents.shape[0])
            batch = np.ascontiguousarray(latents[start:stop], dtype=np.float32)
            out[start:stop] = np.asarray(run(params, jnp.asarray(batch)), dtype=np.float32)
        return out

    # -- Monge map -----------------------------------------------------

    def _monge_module(self, dim_cond: int) -> ConditionalPerturbationNetwork:
        monge = self.cfg["monge"]
        dim_cond_map = int(monge["dim_cond_map"])
        return ConditionalPerturbationNetwork(
            dim_hidden=list(monge["hidden_dims"]),
            dim_data=int(monge["dim_data"]),
            dim_cond=int(dim_cond),
            dim_cond_map=dim_cond_map,
            act_fn=_activation(monge["act_fn"]),
            context_entity_bonds=[(0, int(dim_cond) - 1), (0, int(dim_cond))],
        )

    def _monge_loss(self, module, params, source, target, condition):
        mapped = module.apply(
            {"params": params}, source, condition, int(self.cfg["monge"]["num_contexts"])
        )
        fit = fitting_loss(
            target, mapped, epsilon_fitting=float(self.cfg["monge"]["fitting_loss"]["epsilon_fitting"])
        )
        reg = regularizer(
            source,
            mapped,
            epsilon_regularizer=float(
                self.cfg["monge"]["regularizer"]["epsilon_regularizer"]
            ),
            cost=str(self.cfg["monge"]["regularizer"]["cost"]),
        )
        return fit + reg, {"fitting": fit, "regularizer": reg}

    def _monge_state(self, module, seed: int, steps: int) -> TrainState:
        monge = self.cfg["monge"]
        rng = jax.random.PRNGKey(seed + 1)
        rng, init_rng = jax.random.split(rng)
        params = module.init(
            init_rng,
            jnp.ones((1, int(monge["dim_data"]))),
            jnp.ones((1, int(module.dim_cond))),
        )["params"]
        schedule = optax.cosine_decay_schedule(
            init_value=float(monge["optimizer"]["lr"]),
            decay_steps=int(steps),
            alpha=float(monge["schedule"]["alpha"]),
        )
        optimizer = optax.adamw(
            learning_rate=schedule,
            weight_decay=float(monge["optimizer"]["weight_decay"]),
        )
        return TrainState.create(apply_fn=module.apply, params=params, tx=optimizer)

    # -- Adapter API ---------------------------------------------------

    def fit(self, train, valid, config, run_dir, controls=None):
        if controls is None:
            raise RuntimeError("CMonge needs --controls (the shared DMSO controls)")
        settings = self._settings()
        seed = int(config.get("seed", self.cfg["seed"]))
        run_dir = Path(run_dir)
        checkpoint = run_dir / "checkpoint"
        checkpoint.mkdir(parents=True, exist_ok=True)
        tag = (
            f"{'smoke' if _smoke() else 'full'}:{settings['ae_epochs']}:"
            f"{settings['cm_iters']}"
        )

        train_case = load_npz(train)
        valid_case = load_npz(valid)
        control_sets = Controls(controls)
        log(
            f"backend={jax.default_backend()} train={case_summary(train_case)} "
            f"valid={case_summary(valid_case)} tag={tag}"
        )

        # ---- training data: treated wells + all shared DMSO cells ----
        treated = np.asarray(train_case["log1p_cp10k"], dtype=np.float32)
        controls_flat = control_sets.cells.reshape(-1, N_GENES).astype(np.float32)
        pool = np.vstack([treated, controls_flat])
        log(f"ae pool {pool.shape}")

        ae_module = self._ae_module(seed)
        ae_path = checkpoint / "ae.msgpack"
        ae_meta_path = checkpoint / "ae.json"
        ae_meta = {}
        if ae_meta_path.exists():
            ae_meta = load_json(ae_meta_path)
        reuse_ae = (
            ae_path.exists()
            and ae_meta.get("epochs") == settings["ae_epochs"]
            and ae_meta.get("batch") == settings["ae_batch"]
        )
        if reuse_ae:
            log(f"reusing frozen AE from {ae_path}")
            ae_params = serialization.from_bytes(
                ae_module.init(
                    jax.random.PRNGKey(seed + 2), jnp.empty((1, N_GENES))
                )["params"],
                ae_path.read_bytes(),
            )
            ae_history = ae_meta.get("history", [])
        else:
            ae_params, ae_history = self._train_ae(pool, ae_module, settings, seed)
            ae_path.write_bytes(serialization.to_bytes(ae_params))
            write_json(
                ae_meta_path,
                {
                    "epochs": settings["ae_epochs"],
                    "batch": settings["ae_batch"],
                    "params": int(sum(np.size(v) for v in jax.tree_util.tree_leaves(ae_params))),
                    "history": ae_history,
                    "tag": tag,
                },
            )
        del pool

        # ---- frozen latents ------------------------------------------
        control_latent = self._encode(ae_params, ae_module, controls_flat)
        del controls_flat
        treated_latent = self._encode(ae_params, ae_module, treated)
        valid_latent = self._encode(
            ae_params, ae_module, np.asarray(valid_case["log1p_cp10k"], dtype=np.float32)
        )
        log(
            f"latents control={control_latent.shape} train={treated_latent.shape} "
            f"valid={valid_latent.shape}"
        )

        # ---- conditions ----------------------------------------------
        train_smiles = _drug_smiles(train_case)
        valid_smiles = _drug_smiles(valid_case)
        test_smiles = _drug_smiles(load_npz(config["cases"]["test"]))
        all_smiles = {**train_smiles, **valid_smiles, **test_smiles}
        drug_order = sorted(all_smiles)
        descriptor_rows, width = _rdkit_rows([all_smiles[d] for d in drug_order])
        descriptor_index = {d: i for i, d in enumerate(drug_order)}
        train_rows = descriptor_rows[[descriptor_index[d] for d in sorted(train_smiles)]]
        stats = _fit_descriptor_standardizer(
            train_rows, float(self.cfg["embedding"]["variance_threshold"])
        )
        n_kept = int(np.asarray(stats["keep"]).sum())
        log(
            f"rdkit width={width} kept={n_kept} dropped={width - n_kept} "
            f"(train drugs={len(train_smiles)})"
        )

        train_group_rows: dict[str, np.ndarray] = {}
        for group in sorted(set(str(g) for g in train_case["group_id"])):
            train_group_rows[group] = np.flatnonzero(
                np.asarray(train_case["group_id"]).astype(str) == group
            )
        train_conditions = list(train_group_rows)
        train_drugs = [g.split("|")[0] for g in train_conditions]
        train_doses = np.asarray(
            [float(g.split("|")[1]) for g in train_conditions], dtype=np.float64
        )
        train_cond_vec = _condition_vectors(
            descriptor_rows, drug_order, train_drugs, train_doses, stats
        )
        train_targets = [treated_latent[train_group_rows[g]] for g in train_conditions]
        dim_cond = int(train_cond_vec.shape[1])

        valid_group_rows: dict[str, np.ndarray] = {}
        for group in sorted(set(str(g) for g in valid_case["group_id"])):
            valid_group_rows[group] = np.flatnonzero(
                np.asarray(valid_case["group_id"]).astype(str) == group
            )
        valid_conditions = list(valid_group_rows)
        valid_drugs = [g.split("|")[0] for g in valid_conditions]
        valid_doses = np.asarray(
            [float(g.split("|")[1]) for g in valid_conditions], dtype=np.float64
        )
        valid_cond_vec = _condition_vectors(
            descriptor_rows, drug_order, valid_drugs, valid_doses, stats
        )

        # ---- fixed validation batches (author loss, cold valid drugs) --
        batch = int(settings["cm_batch"])
        valid_batches = []
        for i, group in enumerate(valid_conditions):
            rng = np.random.default_rng(seed + 7919 * (i + 1))
            source_idx = rng.integers(0, control_latent.shape[0], size=batch)
            target_rows = valid_latent[valid_group_rows[group]]
            target_idx = rng.integers(0, target_rows.shape[0], size=batch)
            valid_batches.append(
                (
                    jnp.asarray(control_latent[source_idx]),
                    jnp.asarray(target_rows[target_idx]),
                    jnp.asarray(
                        np.broadcast_to(valid_cond_vec[i], (batch, dim_cond)),
                        dtype=jnp.float32,
                    ),
                )
            )

        # ---- Monge map ------------------------------------------------
        monge_module = self._monge_module(dim_cond)
        state = self._monge_state(monge_module, seed, int(settings["cm_iters"]))

        @jax.jit
        def train_step(state: TrainState, source, target, condition):
            def loss_fn(params):
                return self._monge_loss(
                    monge_module, params, source, target, condition
                )

            (loss, aux), grads = jax.value_and_grad(loss_fn, has_aux=True)(state.params)
            return state.apply_gradients(grads=grads), loss, aux

        @jax.jit
        def eval_loss(params, source, target, condition):
            return self._monge_loss(monge_module, params, source, target, condition)

        def validate(params) -> dict:
            totals = {"total": 0.0, "fitting": 0.0, "regularizer": 0.0}
            for source, target, condition in valid_batches:
                total, aux = eval_loss(params, source, target, condition)
                totals["total"] += float(total)
                totals["fitting"] += float(aux["fitting"])
                totals["regularizer"] += float(aux["regularizer"])
            n = max(len(valid_batches), 1)
            return {k: v / n for k, v in totals.items()}

        resume_path = checkpoint / "resume.msgpack"
        resume_meta_path = checkpoint / "resume.json"
        start_iter = 0
        best_valid = float("inf")
        best_params = None
        history: list[dict] = []
        best_path = checkpoint / "best.msgpack"
        if resume_path.exists() and resume_meta_path.exists():
            resume_meta = load_json(resume_meta_path)
            if resume_meta.get("tag") == tag and resume_meta.get("iter", 0) <= settings["cm_iters"]:
                state = state.replace(
                    params=serialization.from_bytes(state.params, resume_path.read_bytes())
                )
                start_iter = int(resume_meta["iter"])
                best_valid = float(resume_meta.get("best_valid", float("inf")))
                history = resume_meta.get("history", [])
                log(f"resuming CMonge at iteration {start_iter} (best valid {best_valid:.6f})")
            else:
                # A different budget/smoke tag: never let the old weights leak in.
                log("discarding checkpoints from another budget/smoke tag")
                for stale in (best_path, resume_path, resume_meta_path):
                    if stale.exists():
                        stale.unlink()

        if not resume_path.exists():
            # An untrained map is still a usable last checkpoint if the fold is
            # killed before the first grid point.
            resume_path.write_bytes(serialization.to_bytes(state.params))
            write_json(
                resume_meta_path,
                {"tag": tag, "iter": 0, "best_valid": best_valid, "history": []},
            )
        if not best_path.exists():
            # Predict prefers `best.msgpack`; seed it with the iteration-0 map so
            # a fold killed before the first improving eval still has a candidate.
            best_path.write_bytes(serialization.to_bytes(state.params))

        eval_every = int(settings["eval_every"])
        iterations = int(settings["cm_iters"])
        for step in range(start_iter, iterations):
            rng = np.random.default_rng(seed * 1_000_003 + step)
            ci = int(rng.integers(0, len(train_conditions)))
            source_idx = rng.integers(0, control_latent.shape[0], size=batch)
            target_rows = train_targets[ci]
            target_idx = rng.integers(0, target_rows.shape[0], size=batch)
            source = jnp.asarray(control_latent[source_idx])
            target = jnp.asarray(target_rows[target_idx])
            condition = jnp.asarray(
                np.broadcast_to(train_cond_vec[ci], (batch, dim_cond)), dtype=jnp.float32
            )
            state, loss, aux = train_step(state, source, target, condition)
            if (step + 1) % 100 == 0:
                log(
                    f"iter {step + 1}/{iterations} total={float(loss):.6f} "
                    f"fitting={float(aux['fitting']):.6f} "
                    f"regularizer={float(aux['regularizer']):.6f}"
                )
            if (step + 1) % eval_every == 0 or step + 1 == iterations:
                metrics = validate(state.params)
                record = {"iter": step + 1, **metrics}
                history.append(record)
                improved = metrics["total"] < best_valid
                if improved:
                    best_valid = metrics["total"]
                    best_params = jax.tree_util.tree_map(np.array, state.params)
                    best_path.write_bytes(serialization.to_bytes(best_params))
                log(
                    f"eval iter {step + 1} valid={metrics['total']:.6f} "
                    f"fitting={metrics['fitting']:.6f} "
                    f"regularizer={metrics['regularizer']:.6f} "
                    f"best={best_valid:.6f} improved={improved}"
                )
                resume_path.write_bytes(serialization.to_bytes(state.params))
                write_json(
                    resume_meta_path,
                    {
                        "tag": tag,
                        "iter": step + 1,
                        "best_valid": best_valid,
                        "history": history,
                    },
                )

        if best_params is None:
            if best_path.exists():
                best_params = serialization.from_bytes(
                    jax.tree_util.tree_map(np.zeros_like, state.params),
                    best_path.read_bytes(),
                )
            else:
                best_params = jax.tree_util.tree_map(np.array, state.params)
                best_valid = float(validate(state.params)["total"])
                best_path.write_bytes(serialization.to_bytes(best_params))

        meta = {
            "model": "cmonge",
            "tag": tag,
            "smoke": _smoke(),
            "seed": seed,
            "n_genes": N_GENES,
            "ae": self.cfg["ae"],
            "monge": self.cfg["monge"],
            "settings": settings,
            "descriptor": {
                "width": width,
                "keep": [bool(x) for x in np.asarray(stats["keep"])],
                "mean": [float(x) for x in stats["mean"]],
                "std": [float(x) for x in stats["std"]],
            },
            "dim_cond": dim_cond,
            "n_train_conditions": len(train_conditions),
            "n_valid_conditions": len(valid_conditions),
            "best_valid_loss": best_valid,
            "ae_history": ae_history,
            "history": history,
            "ae_params": int(sum(np.size(v) for v in jax.tree_util.tree_leaves(ae_params))),
            "map_params": int(sum(np.size(v) for v in jax.tree_util.tree_leaves(best_params))),
        }
        write_json(checkpoint / "meta.json", meta)
        write_json(
            run_dir / "fit_history.json",
            {
                "smoke": _smoke(),
                "best_valid_loss": best_valid,
                "ae_history": ae_history,
                "history": history,
            },
        )
        log(
            f"checkpoint {checkpoint} best valid {best_valid:.6f} "
            f"map params {meta['map_params']:,} ae params {meta['ae_params']:,}"
        )
        return checkpoint

    def predict(self, cases, checkpoint, run_dir, controls=None):
        if controls is None:
            raise RuntimeError("CMonge needs --controls (the shared DMSO controls)")
        checkpoint = Path(checkpoint)
        meta = load_json(checkpoint / "meta.json")
        seed = int(meta["seed"])
        ae_module = self._ae_module(seed)
        ae_params = serialization.from_bytes(
            ae_module.init(jax.random.PRNGKey(seed + 2), jnp.empty((1, N_GENES)))["params"],
            (checkpoint / "ae.msgpack").read_bytes(),
        )
        dim_cond = int(meta["dim_cond"])
        monge_module = self._monge_module(dim_cond)
        state = self._monge_state(monge_module, seed, int(meta["monge"]["num_train_iters"]))
        best = checkpoint / "best.msgpack"
        if not best.exists():
            best = checkpoint / "resume.msgpack"
        if not best.exists():
            raise RuntimeError(f"no trained Monge map under {checkpoint}")
        params = serialization.from_bytes(state.params, best.read_bytes())

        case = load_npz(cases)
        control_sets = Controls(controls)
        basal = control_sets.cell_for(case["control_group"], case["replicate"])
        plate_mean = control_sets.mean_for(case["control_group"])
        log(f"predict cases={case_summary(case)} basal={basal.shape}")

        stats = {
            "keep": np.asarray(meta["descriptor"]["keep"], dtype=bool),
            "mean": np.asarray(meta["descriptor"]["mean"], dtype=np.float64),
            "std": np.asarray(meta["descriptor"]["std"], dtype=np.float64),
        }
        smiles_by_drug = _drug_smiles(case)
        descriptor_rows, _ = _rdkit_rows([smiles_by_drug[d] for d in sorted(smiles_by_drug)])
        condition = _condition_vectors(
            descriptor_rows,
            sorted(smiles_by_drug),
            [str(d) for d in case["drug_key"]],
            np.asarray(case["dose_um"], dtype=np.float64),
            stats,
        )
        if condition.shape[1] != dim_cond:
            raise RuntimeError(
                f"condition width {condition.shape[1]} != trained {dim_cond}"
            )

        latent = self._encode(ae_params, ae_module, basal)
        mapped = np.empty_like(latent)

        @jax.jit
        def transport(x, c):
            return monge_module.apply(
                {"params": params}, x, c, int(self.cfg["monge"]["num_contexts"])
            )

        for start in range(0, latent.shape[0], PREDICT_BATCH):
            stop = min(start + PREDICT_BATCH, latent.shape[0])
            mapped[start:stop] = np.asarray(
                transport(
                    jnp.asarray(latent[start:stop]),
                    jnp.asarray(condition[start:stop]),
                ),
                dtype=np.float32,
            )
        decoded = self._decode(ae_params, ae_module, mapped)
        delta = decoded - plate_mean
        log(f"predict {delta.shape} finite={np.isfinite(delta).all()}")
        return ensure_float32(delta)


if __name__ == "__main__":
    raise SystemExit(run_cli(CMongeAdapter()))
