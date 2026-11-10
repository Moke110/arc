"""Model-side fixed file contract. Held-out tenv expression has no loader here."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch


DOSE_BOUNDARIES = np.array([.21, .41, .71, 1.01, 1.51, 3.1, 4.1, 7.1, 12.1], dtype=np.float32)
FIT_ROLES = ("train", "validation")
SOURCES = ("tahoe", "tenv")
AXIS_MODELS = {
    "landmark": "XPert_base",
    "G2000": "XPert_G2000",
    "G5000": "XPert_G5000",
    "G10000": "XPert_G10000",
    "G18164": "XPert_G18164",
}


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


class Predictor:
    def __init__(self, root, fold):
        if fold not in (0, 1, 2):
            raise ValueError("unknown tenv fold")
        self.root = Path(root).absolute()
        if self.root.resolve() != self.root or self.root.name != "predictor":
            raise ValueError("expected a real predictor directory")
        manifest = self.root / "manifest.json"
        if manifest.resolve() != manifest:
            raise ValueError("predictor manifest cannot be a symlink")
        self.manifest = json.loads(manifest.read_text())
        self.digest = sha256(manifest)
        self.fold = fold
        self.verified = set()
        if self.manifest["schema"] != "step09-predictor-v1":
            raise ValueError("unsupported Step 09 predictor")
        if fold not in self.manifest.get("ready_folds", (0, 1, 2)):
            raise ValueError("requested fold is not materialized for this gene axis")
        self.model = AXIS_MODELS.get(self.manifest["axis"])
        if self.model is None:
            raise ValueError("gene axis is evaluator-only or has no XPert model")
        self.genes = self.table("genes.csv")
        if (len(self.genes) != self.manifest["gene_count"] or self.genes.hgnc_id.duplicated().any()
                or self.genes.identity_index.duplicated().any()):
            raise ValueError("invalid predictor gene axis")

    def file(self, relative):
        # Paths are fixed by callers. A manifest can supply hashes only.
        path = self.root / relative
        if path.resolve() != path or not path.is_relative_to(self.root):
            raise ValueError(f"predictor path escaped or is a symlink: {relative}")
        if relative not in self.verified:
            if sha256(path) != self.manifest["files"][relative]:
                raise ValueError(f"predictor digest mismatch: {relative}")
            self.verified.add(relative)
        return path

    def table(self, relative):
        return pd.read_csv(self.file(relative), keep_default_na=False)

    def array(self, relative):
        return np.load(self.file(relative), mmap_mode="r", allow_pickle=False)

    def identity(self):
        cm = self.array("identity.npy")
        extra = self.array("identity_hdgfl3.npy")
        if cm.ndim != 2 or extra.shape != cm.shape[1:]:
            raise ValueError("invalid frozen identity lookup")
        combined = np.concatenate((cm, extra[None, :]))
        selected = self.genes.identity_index.to_numpy()
        if np.any(selected < 0) or np.any(selected >= len(combined)):
            raise ValueError("invalid predictor identity indices")
        return torch.from_numpy(combined[selected])

    def dataset(self, source, role, targets):
        if source not in SOURCES or role not in (*FIT_ROLES, "test"):
            raise ValueError("unknown predictor source or role")
        if role == "test" and (targets or source != "tenv"):
            raise ValueError("tenv test truth is evaluator-only")
        prefix = source if source == "tahoe" else f"tenv/fold_{self.fold}"
        cases = self.table(f"{prefix}/{role}.csv")
        treated = self.array(f"{prefix}/{role}_expression.npy") if targets else None
        return Profiles(cases, self.array(f"{source}/controls.npy"),
                        self.array(f"{source}/morgan.npy"), treated, len(self.genes))


class Profiles:
    def __init__(self, cases, controls, morgan, treated, genes):
        self.cases = cases
        self.controls = controls
        self.morgan = morgan
        self.treated = treated
        self.control_index = cases.control_index.to_numpy()
        self.morgan_index = cases.morgan_index.to_numpy()
        doses = cases.dose_um.to_numpy(dtype=np.float32)
        self.dose = np.searchsorted(DOSE_BOUNDARIES, doses, side="right")
        if (not len(cases) or cases.case_id.duplicated().any()
                or np.any(~np.isfinite(doses)) or np.any(doses <= 0)
                or np.any(self.control_index < 0) or np.any(self.control_index >= len(controls))
                or np.any(self.morgan_index < 0) or np.any(self.morgan_index >= len(morgan))
                or controls.ndim != 2 or controls.shape[1] != genes
                or morgan.ndim != 2 or morgan.shape[1] != 2048
                or (treated is not None and treated.shape != (len(cases), genes))):
            raise ValueError("invalid predictor cohort or matrix alignment")

    def __len__(self):
        return len(self.cases)

    def batch(self, rows, device):
        inputs = (self.controls[self.control_index[rows]], self.morgan[self.morgan_index[rows]].astype(np.float32),
                  self.dose[rows])
        tensors = tuple(torch.from_numpy(values).to(device) for values in inputs)
        treated = None if self.treated is None else torch.from_numpy(self.treated[rows]).to(device)
        return tensors, treated
