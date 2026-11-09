"""Probe the Latent XPert treatment path with a released frozen autoencoder."""
import csv
import hashlib
import json
from math import isfinite, log10
import time

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

from autoencoder import WIDTH, CrossAttention
from formal_fit import verify_formal_release
from inputs import sha256
from short_fit import ExpressionRelease, attention_runtime, load_checkpoint, measure, runtime


SEED = 20261109
REFERENCE_DOSE_UM = 0.05
MORGAN = {"radius": 2, "bits": 2048, "chirality": True, "binary": True, "invariants": "atom"}
INPUT_FILES = ("genes.csv", "wells.csv", "controls.csv", "controls.npy", "drugs.csv", "morgan.npy")


class AttentionBlock(nn.Module):
    def __init__(self):
        super().__init__()
        self.attention = CrossAttention()
        self.norm = nn.LayerNorm(WIDTH)

    def forward(self, queries, tokens):
        return self.norm(queries + self.attention(queries, tokens))


class TreatmentBranch(nn.Module):
    block_order = ("CA", "SA", "SA", "CA")

    def __init__(self):
        super().__init__()
        self.morgan_adapter = nn.Linear(MORGAN["bits"], WIDTH)
        self.dose_adapter = nn.Linear(1, WIDTH)
        self.blocks = nn.ModuleList(AttentionBlock() for _ in self.block_order)

    def conditions(self, morgan, log_dose):
        return torch.stack((self.morgan_adapter(morgan), self.dose_adapter(log_dose)), dim=1)

    def forward(self, latent, conditions):
        for kind, block in zip(self.block_order, self.blocks, strict=True):
            latent = block(latent, conditions if kind == "CA" else latent)
        return latent


def find_row(path, key, value):
    with path.open(newline="") as file:
        rows = [row for row in csv.DictReader(file) if row[key] == value]
    if len(rows) != 1:
        raise ValueError(f"{path.name}: expected one row with {key}={value}")
    return rows[0]


def load_case(directory, case_id, contract, verification_pin):
    # Only open predictor inputs. The release also contains evaluator-only truth.
    content = (directory / "manifest.json").read_bytes()
    digest = hashlib.sha256(content).hexdigest()
    if verification_pin is not None and digest != verification_pin:
        raise ValueError("tenv input manifest differs from the verified Step 03 release")
    manifest = json.loads(content)
    gene_count = contract["architecture"]["gene_count"]
    if (manifest["schema"] != "tenv-input-v1" or manifest["gene_count"] != gene_count
            or manifest["files"]["genes.csv"] != contract["axis_sha256"]):
        raise ValueError("tenv input gene axis differs from the autoencoder contract")
    if manifest["morgan"] != MORGAN:
        raise ValueError("tenv input requires the shared chirality-aware binary Morgan-2048")
    for name in INPUT_FILES:
        if sha256(directory / name) != manifest["files"][name]:
            raise ValueError(f"tenv input digest mismatch: {name}")

    case = find_row(directory / "wells.csv", "case_id", case_id)
    control = find_row(directory / "controls.csv", "control_index", case["control_index"])
    drug = find_row(directory / "drugs.csv", "drug_key", case["drug_key"])
    dose = float(case["dose_um"])
    if not isfinite(dose) or dose <= 0:
        raise ValueError("tenv treatment requires a positive measured dose")
    log_dose = log10(dose) - log10(REFERENCE_DOSE_UM)

    controls = np.load(directory / "controls.npy", mmap_mode="r", allow_pickle=False)
    morgan = np.load(directory / "morgan.npy", mmap_mode="r", allow_pickle=False)
    control_index, morgan_index = int(control["control_index"]), int(drug["morgan_index"])
    dmso_wells = int(control["n_wells"])
    if (controls.ndim != 2 or controls.shape[1] != gene_count or controls.dtype != np.float32
            or not 0 <= control_index < len(controls) or dmso_wells <= 0):
        raise ValueError("invalid same-plate DMSO baseline input")
    if (morgan.shape != (manifest["drug_count"], MORGAN["bits"]) or morgan.dtype != np.uint8
            or not 0 <= morgan_index < len(morgan)):
        raise ValueError("invalid Morgan-2048 input")
    baseline, fingerprint = controls[control_index].copy(), morgan[morgan_index].copy()
    if not np.isfinite(baseline).all() or (baseline < 0).any() or (fingerprint > 1).any():
        raise ValueError("invalid DMSO expression values or Morgan bits")
    return baseline, fingerprint, log_dose, {
        "tenv_manifest_sha256": digest,
        "tenv_files": {name: manifest["files"][name] for name in INPUT_FILES},
        "case": {
            "case_id": case_id, "drug_key": case["drug_key"], "plate": control["plate"],
            "dose_um": dose, "reference_dose_um": REFERENCE_DOSE_UM, "log10_relative_dose": log_dose,
            "control_index": control_index, "dmso_wells": dmso_wells, "morgan_index": morgan_index,
            "morgan_active_bits": int(fingerprint.sum()),
            "baseline_sha256": hashlib.sha256(baseline.tobytes()).hexdigest(),
        },
    }


def gradient_norm(module):
    gradients = [parameter.grad for parameter in module.parameters()]
    if any(grad is None or not torch.isfinite(grad).all() for grad in gradients):
        raise ValueError("treatment branch has missing or non-finite gradients")
    norm = float(torch.linalg.vector_norm(torch.cat([grad.flatten() for grad in gradients])))
    if not isfinite(norm) or norm <= 0:
        raise ValueError("treatment branch gradient norm must be finite and positive")
    return norm


def probe_treatment(autoencoder, device, baseline, fingerprint, log_dose):
    autoencoder.eval().requires_grad_(False)
    frozen = [parameter.detach().clone() for parameter in autoencoder.parameters()]
    torch.manual_seed(SEED)
    branch = TreatmentBranch().to(device)
    expression = torch.from_numpy(baseline[None, :]).to(device)
    morgan = torch.from_numpy(fingerprint[None, :]).to(device=device, dtype=torch.float32)
    dose = torch.tensor([[log_dose]], device=device, dtype=torch.float32)
    with attention_runtime(device):
        with torch.no_grad():
            latent = autoencoder.encode(expression)
        conditions = branch.conditions(morgan, dose)
        treated = branch(latent, conditions)
        post_treatment = autoencoder.decode(treated).float()
    delta = post_treatment - expression

    # Synthetic full-axis target for the backward path only. No observed treated target is loaded.
    target = torch.linspace(-1, 1, delta.shape[-1], device=device).expand_as(delta)
    centered = delta - delta.mean(dim=-1, keepdim=True)
    centered_target = target - target.mean(dim=-1, keepdim=True)
    loss = F.mse_loss(delta, target) + (1 - F.cosine_similarity(centered, centered_target)).mean()
    if any(not torch.isfinite(tensor).all()
           for tensor in (latent, conditions, treated, post_treatment, delta, loss)):
        raise ValueError("treatment smoke produced non-finite values")
    loss.backward()

    gradients = {name: gradient_norm(module) for name, module in branch.named_children()}
    unchanged = all(torch.equal(before, after) for before, after in
                    zip(frozen, autoencoder.parameters(), strict=True))
    gradients_absent = all(parameter.grad is None and not parameter.requires_grad
                           for parameter in autoencoder.parameters())
    if not unchanged or not gradients_absent:
        raise ValueError("autoencoder changed or received gradients during the treatment probe")
    post_array = post_treatment.detach().cpu().numpy()
    delta_array = delta.detach().cpu().numpy()
    delta_error = float(np.max(np.abs(delta_array - (post_array - baseline))))
    if delta_error != 0:
        raise ValueError("decoded treatment delta differs from post-treatment minus measured baseline")
    return {
        "seed": SEED, "block_order": list(branch.block_order), "latent_shape": list(latent.shape),
        "condition_shape": list(conditions.shape), "treated_latent_shape": list(treated.shape),
        "post_treatment_shape": list(post_treatment.shape), "delta_shape": list(delta.shape),
        "delta_max_abs_error": delta_error, "finite": True,
        "autoencoder_unchanged": unchanged, "autoencoder_gradients_absent": gradients_absent,
        "probe_objective": "full-axis MSE + (1 - Pearson) against a synthetic ramp",
        "probe_loss": float(loss.detach()), "gradient_l2": gradients,
    }


def smoke(args, source_pins, canonical_checkpoint, autoencoder_record, tenv_record):
    start = time.perf_counter()
    execution = runtime(args.device, args.fixture)
    checkpoint_sha = sha256(args.checkpoint)
    if args.checkpoint.resolve() == canonical_checkpoint.resolve():
        pinned_checkpoint = json.loads(autoencoder_record.read_text())["checkpoint_sha256"]
        if checkpoint_sha != pinned_checkpoint:
            raise ValueError("checkpoint differs from the verified Step 07 autoencoder release")
    release = ExpressionRelease(args.release, args.fixture, source_pins)
    model, payload = load_checkpoint(args.checkpoint, release.contract, args.device)
    metadata = payload["metadata"]
    if metadata.get("fit_stage") != "full_data_refit":
        raise ValueError("treatment smoke requires a full-data refit checkpoint")
    provenance = verify_formal_release(args.checkpoint.parent, metadata, checkpoint_sha, release)
    pin = None if args.fixture else json.loads(tenv_record.read_text())["manifest_sha256"]
    baseline, fingerprint, log_dose, case = load_case(args.tenv_input, args.case_id, release.contract, pin)
    report = probe_treatment(model, args.device, baseline, fingerprint, log_dose)
    return {**release.contract, **execution, **provenance, **case, **report,
            **measure(start, args.device), "checkpoint_sha256": checkpoint_sha}
