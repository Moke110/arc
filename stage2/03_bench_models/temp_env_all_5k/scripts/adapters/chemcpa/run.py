#!/usr/bin/env python
"""chemCPA adapter (ticket 13).

Author recipe (research/01-... §3, ticket 09):

- from-scratch ``ComPert`` at 5,000 genes (the L1000 977-gene pretrain is
  skipped: it is a 978-landmark artifact and its gain is deliberately not
  reproduced).
- ComPert: ``dim=32``, encoder/decoder 256x4, drug-embedding encoder 128x4,
  adversary 128x3, amortized doser 64x3, Gaussian NLL reconstruction plus
  adversarial drug/covariate terms (reg 9.10 / 16.17, penalty 0.455).
- Three Adam optimizers (AE 1.575e-3, adversary 1.192e-3, doser 1.575e-3),
  StepLR(6, 0.9), batch 256, 201 epochs, early stop patience 50.
- Molecule encoder G is frozen: the author's ``rdkit`` embedding is RDKit2D
  descriptors, standardised and variance-filtered.  We fit its mean/std and
  low-variance filter on the fold's **training drugs** only (the author fits on
  the whole dataset; train-only is stricter for cold drugs).

The shipping Hydra pipeline cannot run a new dataset as shipped (dict
val_dataloader, unconditional CUDA collate, DEG bookkeeping), so this adapter
drives ``ComPert`` directly with the author's optimizers, schedule, loss and
adversarial alternation.  Deviations are recorded in
``models/chemcpa/FIDELITY.md``.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import torch

BENCH = Path(__file__).resolve().parents[3]
ROOT = BENCH.parents[2]
sys.path.insert(0, str(BENCH / "scripts"))
sys.path.insert(0, str(ROOT / "pert_model/chemCPA/source"))

from adapter_api import Adapter, run_cli  # noqa: E402
from adapters.common import (  # noqa: E402
    N_GENES,
    Controls,
    apply_standardizer,
    case_summary,
    ensure_float32,
    fit_standardizer,
    load_json,
    load_npz,
    log,
    rdkit2d_descriptors,
    write_json,
)

from chemCPA.model import ComPert  # noqa: E402


def _smoke() -> bool:
    return os.environ.get("BENCH_SMOKE", "").strip() not in ("", "0", "false", "False")


def _normalise_dose(dose: np.ndarray, dmin: float, dmax: float) -> np.ndarray:
    if dmax <= dmin:
        return np.zeros_like(dose, dtype=np.float32)
    return ((dose - dmin) / (dmax - dmin)).astype(np.float32)


class ChemCPAAdapter(Adapter):
    name = "chemcpa"

    def __init__(self) -> None:
        self.cfg = load_json(BENCH / "models/chemcpa/config.json")

    # -- drug bookkeeping ---------------------------------------------

    def _drug_table(self, cases: dict[str, dict], train_drugs: set[str]):
        """sorted drug keys, control at index 0, plus frozen RDKit2D rows."""

        names = sorted({str(d) for case in cases.values() for d in case["drug_key"]})
        smiles_by_drug: dict[str, str] = {}
        for case in cases.values():
            for drug, smiles in zip(case["drug_key"], case["smiles"]):
                if str(drug) not in smiles_by_drug:
                    smiles_by_drug[str(drug)] = str(smiles)
        smiles = [smiles_by_drug[d] for d in names]
        raw = rdkit2d_descriptors(smiles)
        raw = raw[:, 1:]  # the author drops descriptor 0

        train_rows = np.array([names.index(d) for d in names if d in train_drugs])
        stats = fit_standardizer(raw[train_rows], min_std=float(self.cfg["embedding"]["min_std"]))
        embedding = apply_standardizer(raw, stats)
        embedding = np.vstack(
            [np.zeros((1, embedding.shape[1]), dtype=np.float32), embedding]
        )
        table = torch.nn.Embedding.from_pretrained(
            torch.from_numpy(embedding.astype(np.float32)), freeze=True
        )
        index = {name: i + 1 for i, name in enumerate(names)}
        return names, index, table, stats

    # -- tensors -------------------------------------------------------

    def _training_tensors(self, case, controls: Controls, drug_index, dmin, dmax):
        genes = np.asarray(case["log1p_cp10k"], dtype=np.float32)
        drugs = np.array([drug_index[str(d)] for d in case["drug_key"]], dtype=np.int64)
        dose = _normalise_dose(np.asarray(case["dose_um"], dtype=np.float32), dmin, dmax)
        return (
            torch.from_numpy(genes),
            torch.from_numpy(drugs),
            torch.from_numpy(dose),
        )

    def _control_tensors(self, controls: Controls):
        genes = controls.cells.reshape(-1, N_GENES).astype(np.float32)
        n = genes.shape[0]
        return (
            torch.from_numpy(genes),
            torch.zeros(n, dtype=torch.long),
            torch.zeros(n, dtype=torch.float32),
        )

    def _model(self, drug_table, num_covariates, seed, device):
        return ComPert(
            num_genes=N_GENES,
            num_drugs=int(drug_table.weight.shape[0]),
            num_covariates=num_covariates,
            device=str(device),
            seed=seed,
            hparams=dict(self.cfg["hparams"]),
            drug_embeddings=drug_table,
            use_drugs_idx=True,
            append_layer_width=None,
            decoder_activation=str(self.cfg["additional_params"]["decoder_activation"]),
            doser_type=str(self.cfg["additional_params"]["doser_type"]),
        )

    # -- Adapter API ---------------------------------------------------

    def fit(self, train, valid, config, run_dir, controls=None):
        if controls is None:
            raise RuntimeError("chemCPA needs --controls (the shared DMSO means)")
        cfg = dict(self.cfg)
        if _smoke():
            cfg["training"] = {
                **cfg["training"],
                "num_epochs": int(cfg["smoke"]["max_epochs"]),
                "max_batches": int(cfg["smoke"]["max_batches"]),
            }
        seed = int(config.get("seed", cfg["training"]["seed"]))
        torch.manual_seed(seed)
        np.random.seed(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        controls_obj = Controls(controls)
        train_case = load_npz(train)
        valid_case = load_npz(valid)
        test_case = load_npz(config["cases"]["test"])
        log(
            f"device={device} train={case_summary(train_case)} "
            f"valid={case_summary(valid_case)} test={case_summary(test_case)}"
        )

        names, drug_index, drug_table, stats = self._drug_table(
            {"train": train_case, "valid": valid_case, "test": test_case},
            {str(d) for d in train_case["drug_key"]},
        )
        dmin = float(np.min(train_case["dose_um"]))
        dmax = float(np.max(train_case["dose_um"]))
        log(f"drugs={len(names)} embedding_shape={tuple(drug_table.weight.shape)} dose=[{dmin},{dmax}]")

        train_genes, train_drugs, train_dose = self._training_tensors(
            train_case, controls_obj, drug_index, dmin, dmax
        )
        ctrl_genes, ctrl_drugs, ctrl_dose = self._control_tensors(controls_obj)
        train_genes = torch.cat([train_genes, ctrl_genes])
        train_drugs = torch.cat([train_drugs, ctrl_drugs])
        train_dose = torch.cat([train_dose, ctrl_dose])
        valid_genes, valid_drugs, valid_dose = self._training_tensors(
            valid_case, controls_obj, drug_index, dmin, dmax
        )
        valid_genes = torch.cat([valid_genes, ctrl_genes])
        valid_drugs = torch.cat([valid_drugs, ctrl_drugs])
        valid_dose = torch.cat([valid_dose, ctrl_dose])
        del train_case, valid_case, test_case

        model = self._model(drug_table, [1], seed, device).to(device)
        hparams = model.hparams
        has_covariates = model.num_covariates[0] > 0

        ae_params = (
            list(model.encoder.parameters())
            + list(model.decoder.parameters())
            + (list(model.drug_embedding_encoder.parameters()) if model.drug_embedding_encoder is not None else [])
        )
        if has_covariates:
            for embedding in model.covariates_embeddings:
                ae_params += list(embedding.parameters())
        adv_params = list(model.adversary_drugs.parameters())
        if has_covariates:
            for adversary in model.adversary_covariates:
                adv_params += list(adversary.parameters())
        dose_params = list(model.dosers.parameters())

        opt_ae = torch.optim.Adam(
            ae_params, lr=float(hparams["autoencoder_lr"]), weight_decay=float(hparams["autoencoder_wd"])
        )
        opt_adv = torch.optim.Adam(
            adv_params, lr=float(hparams["adversary_lr"]), weight_decay=float(hparams["adversary_wd"])
        )
        opt_dose = torch.optim.Adam(
            dose_params, lr=float(hparams["dosers_lr"]), weight_decay=float(hparams["dosers_wd"])
        )
        step_size = int(hparams["step_size_lr"])
        sch_ae = torch.optim.lr_scheduler.StepLR(opt_ae, step_size=step_size, gamma=0.9)
        sch_adv = torch.optim.lr_scheduler.StepLR(opt_adv, step_size=step_size, gamma=0.9)
        sch_dose = torch.optim.lr_scheduler.StepLR(opt_dose, step_size=step_size, gamma=0.9)

        batch_size = int(hparams["batch_size"])
        epochs = int(cfg["training"]["num_epochs"])
        patience = int(cfg["additional_params"]["patience"])
        adv_steps = int(hparams["adversary_steps"])
        reg_adv = float(hparams["reg_adversary"])
        reg_cov = float(hparams["reg_adversary_cov"])
        penalty = float(hparams["penalty_adversary"])
        covariate = torch.ones(1, 1, device=device)

        checkpoint = Path(run_dir) / "checkpoint"
        checkpoint.mkdir(parents=True, exist_ok=True)
        best_path = checkpoint / "model.pt"
        best_valid = float("inf")
        stale = 0
        history: list[dict] = []
        global_step = 0
        generator = torch.Generator().manual_seed(seed)
        max_batches = int(cfg["training"].get("max_batches", 10**9)) if _smoke() else 10**9

        def compute(model, genes, drugs, dose, training: bool):
            covariates = [covariate.expand(genes.shape[0], 1)]
            reconstruction, _, latents = model.predict(
                genes, drugs_idx=drugs, dosages=dose, covariates=covariates,
                return_latent_basal=True,
            )
            dim = reconstruction.shape[1] // 2
            mean, var = reconstruction[:, :dim], reconstruction[:, dim:]
            latent_basal = latents[0]
            recon = model.loss_autoencoder(input=mean, target=genes, var=var)
            adv_drug = model.loss_adversary_drugs(
                model.adversary_drugs(latent_basal),
                torch.nn.functional.one_hot(drugs, num_classes=model.num_drugs).float(),
            )
            adv_cov = torch.zeros((), device=genes.device)
            for i, adversary in enumerate(model.adversary_covariates):
                adv_cov = adv_cov + model.loss_adversary_covariates[i](
                    adversary(latent_basal), covariates[i].argmax(1)
                )
            return recon, adv_drug, adv_cov, latent_basal

        for epoch in range(epochs):
            model.train()
            order = torch.randperm(len(train_genes), generator=generator).numpy()
            totals = {"recon": 0.0, "adv_drug": 0.0, "adv_cov": 0.0}
            seen = 0
            for start in range(0, len(order), batch_size):
                if start // batch_size >= max_batches:
                    break
                idx = order[start : start + batch_size]
                genes = train_genes[idx].to(device)
                drugs = train_drugs[idx].to(device)
                dose = train_dose[idx].to(device)
                recon, adv_drug, adv_cov, latent_basal = compute(model, genes, drugs, dose, True)

                if global_step % adv_steps == 0:
                    opt_adv.zero_grad(set_to_none=True)
                    penalty_term = torch.zeros((), device=device)
                    if penalty > 0:
                        def gradient_penalty(output):
                            grad = torch.autograd.grad(output, latent_basal, create_graph=True)[0]
                            return (grad ** 2).mean()

                        penalty_term = penalty_term + gradient_penalty(
                            model.adversary_drugs(latent_basal).sum()
                        )
                        for adversary in model.adversary_covariates:
                            penalty_term = penalty_term + gradient_penalty(adversary(latent_basal).sum())
                    loss = adv_drug + adv_cov + penalty * penalty_term
                    loss.backward()
                    opt_adv.step()
                else:
                    opt_ae.zero_grad(set_to_none=True)
                    opt_dose.zero_grad(set_to_none=True)
                    loss = recon - reg_adv * adv_drug - reg_cov * adv_cov
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(ae_params, 1.0)
                    torch.nn.utils.clip_grad_norm_(dose_params, 1.0)
                    opt_ae.step()
                    opt_dose.step()

                totals["recon"] += float(recon.item()) * len(idx)
                totals["adv_drug"] += float(adv_drug.item()) * len(idx)
                totals["adv_cov"] += float(adv_cov.item()) * len(idx)
                seen += len(idx)
                global_step += 1

            sch_ae.step()
            sch_adv.step()
            sch_dose.step()

            model.eval()
            valid_total, valid_seen = 0.0, 0
            with torch.no_grad():
                for start in range(0, len(valid_genes), batch_size):
                    stop = min(start + batch_size, len(valid_genes))
                    genes = valid_genes[start:stop].to(device)
                    drugs = valid_drugs[start:stop].to(device)
                    dose = valid_dose[start:stop].to(device)
                    covariates = [covariate.expand(genes.shape[0], 1)]
                    reconstruction, _ = model.predict(
                        genes, drugs_idx=drugs, dosages=dose, covariates=covariates,
                        return_latent_basal=False,
                    )
                    dim = reconstruction.shape[1] // 2
                    mean, var = reconstruction[:, :dim], reconstruction[:, dim:]
                    loss = model.loss_autoencoder(input=mean, target=genes, var=var)
                    valid_total += float(loss.item()) * genes.shape[0]
                    valid_seen += genes.shape[0]
            valid_loss = valid_total / max(valid_seen, 1)
            record = {
                "epoch": epoch,
                "recon": totals["recon"] / max(seen, 1),
                "adv_drug": totals["adv_drug"] / max(seen, 1),
                "adv_cov": totals["adv_cov"] / max(seen, 1),
                "valid_loss": valid_loss,
                "global_step": global_step,
            }
            history.append(record)
            log(
                f"epoch {epoch + 1}/{epochs} recon={record['recon']:.6f} "
                f"adv_drug={record['adv_drug']:.4f} valid={valid_loss:.6f}"
            )

            if valid_loss < best_valid:
                best_valid = valid_loss
                stale = 0
                torch.save(
                    {
                        "model": "chemcpa",
                        "state_dict": model.state_dict(),
                        "drug_names": names,
                        "drug_index": drug_index,
                        "drug_embedding": drug_table.weight.detach().cpu(),
                        "standardizer": stats,
                        "dose_min": dmin,
                        "dose_max": dmax,
                        "hparams": dict(model.hparams),
                        "additional": dict(cfg["additional_params"]),
                        "seed": seed,
                        "n_genes": N_GENES,
                        "best_valid_loss": best_valid,
                        "epoch": epoch,
                    },
                    best_path,
                )
            else:
                stale += 1
                if not _smoke() and stale > patience:
                    log(f"early stop after {stale} epochs without improvement")
                    break

        if not best_path.exists():
            torch.save(
                {
                    "model": "chemcpa",
                    "state_dict": model.state_dict(),
                    "drug_names": names,
                    "drug_index": drug_index,
                    "drug_embedding": drug_table.weight.detach().cpu(),
                    "standardizer": stats,
                    "dose_min": dmin,
                    "dose_max": dmax,
                    "hparams": dict(model.hparams),
                    "additional": dict(cfg["additional_params"]),
                    "seed": seed,
                    "n_genes": N_GENES,
                    "best_valid_loss": best_valid,
                    "epoch": epochs - 1,
                },
                best_path,
            )

        write_json(
            Path(run_dir) / "fit_history.json",
            {"history": history, "best_valid_loss": best_valid, "smoke": _smoke()},
        )
        log(f"checkpoint {best_path} (best valid {best_valid:.6f})")
        return checkpoint

    def predict(self, cases, checkpoint, run_dir, controls=None):
        if controls is None:
            raise RuntimeError("chemCPA needs --controls (the shared DMSO means)")
        stored = torch.load(
            Path(checkpoint) / "model.pt", map_location="cpu", weights_only=False
        )
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        table = torch.nn.Embedding.from_pretrained(
            stored["drug_embedding"].float(), freeze=True
        )
        model = ComPert(
            num_genes=N_GENES,
            num_drugs=int(table.weight.shape[0]),
            num_covariates=[1],
            device=str(device),
            seed=int(stored["seed"]),
            hparams=dict(stored["hparams"]),
            drug_embeddings=table,
            use_drugs_idx=True,
            append_layer_width=None,
            decoder_activation=str(stored["additional"]["decoder_activation"]),
            doser_type=str(stored["additional"]["doser_type"]),
        ).to(device)
        model.load_state_dict(stored["state_dict"])
        model.eval()

        case = load_npz(cases)
        control_sets = Controls(controls)
        basal = control_sets.mean_for(case["control_group"])
        drugs = np.array(
            [stored["drug_index"][str(d)] for d in case["drug_key"]], dtype=np.int64
        )
        dose = _normalise_dose(
            np.asarray(case["dose_um"], dtype=np.float32),
            float(stored["dose_min"]),
            float(stored["dose_max"]),
        )
        out = np.empty((len(basal), N_GENES), dtype=np.float32)
        batch_size = 256
        covariate = torch.ones(1, 1, device=device)
        with torch.no_grad():
            for start in range(0, len(basal), batch_size):
                stop = min(start + batch_size, len(basal))
                genes = torch.from_numpy(basal[start:stop]).to(device)
                drug = torch.from_numpy(drugs[start:stop]).to(device)
                d = torch.from_numpy(dose[start:stop]).to(device)
                covariates = [covariate.expand(genes.shape[0], 1)]
                reconstruction, _ = model.predict(
                    genes, drugs_idx=drug, dosages=d, covariates=covariates,
                    return_latent_basal=False,
                )
                dim = reconstruction.shape[1] // 2
                mean = reconstruction[:, :dim]
                out[start:stop] = (mean - genes).cpu().numpy()
        log(f"predict {out.shape} finite={np.isfinite(out).all()}")
        return ensure_float32(out)


if __name__ == "__main__":
    raise SystemExit(run_cli(ChemCPAAdapter()))
