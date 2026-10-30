"""Shared plumbing for the ticket-13 real-model adapters (PRnet, CPA, chemCPA).

None of these helpers import model internals.  They only read the frozen
benchmark case bundles, the shared DMSO controls, and the adapter's own
`models/<model>/config.json`; every adapter keeps its model code behind its
own subprocess environment.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

N_GENES = 5000


# ---------------------------------------------------------------------------
# IO


def load_npz(path: Path) -> dict[str, np.ndarray]:
    """Load every array of an ``.npz`` bundle into memory (adapter process)."""

    with np.load(Path(path), allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def load_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def log(message: str) -> None:
    print(f"[adapter] {message}", flush=True)


# ---------------------------------------------------------------------------
# Controls


class Controls:
    """The 16 shared DMSO groups: ``control_group -> mean(cp10k log1p)``."""

    def __init__(self, path: Path):
        with np.load(Path(path), allow_pickle=False) as archive:
            self.groups = [str(g) for g in archive["control_groups"]]
            self.mean = np.asarray(archive["mean"], dtype=np.float32)
            self.cells = np.asarray(archive["cells"], dtype=np.float32)
        self.index = {group: i for i, group in enumerate(self.groups)}

    def mean_for(self, control_groups: np.ndarray) -> np.ndarray:
        groups = [str(g) for g in control_groups]
        try:
            rows = [self.index[g] for g in groups]
        except KeyError as exc:  # fail closed on an unknown plate
            raise KeyError(f"unknown control group in case bundle: {exc}") from exc
        return self.mean[rows]

    def cell_for(self, control_groups: np.ndarray, replicate: np.ndarray | None = None) -> np.ndarray:
        """One DMSO well per case (used by CPA's cell-level predict path).

        ``replicate`` is 0-based in the case bundle; it selects one of the 12
        DMSO cells on that plate.  The cell is only a basal input, never truth.
        """

        groups = [str(g) for g in control_groups]
        if replicate is None:
            replicate = np.zeros(len(groups), dtype=np.int32)
        replicate = np.mod(np.asarray(replicate, dtype=np.int32), self.cells.shape[1])
        return np.stack(
            [self.cells[self.index[g], int(r)] for g, r in zip(groups, replicate)]
        ).astype(np.float32)


# ---------------------------------------------------------------------------
# Chemical structure


def morgan_fingerprints(
    smiles: np.ndarray | list[str],
    n_bits: int = 1024,
    radius: int = 2,
    use_features: bool = True,
) -> np.ndarray:
    """Morgan (FCFP when ``use_features``) bit vectors, cached per SMILES."""

    from rdkit import Chem
    from rdkit.Chem import AllChem

    unique = sorted({str(s) for s in smiles if str(s)})
    cache: dict[str, np.ndarray] = {}
    for value in unique:
        molecule = Chem.MolFromSmiles(value)
        if molecule is None:
            cache[value] = np.zeros(n_bits, dtype=np.float32)
            continue
        fingerprint = AllChem.GetMorganFingerprintAsBitVect(
            molecule, radius, useFeatures=use_features, nBits=n_bits
        ).ToBitString()
        cache[value] = np.fromiter(
            (1.0 if bit == "1" else 0.0 for bit in fingerprint),
            dtype=np.float32,
            count=n_bits,
        )
    rows = np.zeros((len(smiles), n_bits), dtype=np.float32)
    for i, value in enumerate(smiles):
        rows[i] = cache.get(str(value), np.zeros(n_bits, dtype=np.float32))
    return rows


def rdkit2d_descriptors(smiles: np.ndarray | list[str]) -> np.ndarray:
    """RDKit2D descriptor matrix per SMILES (chemCPA's author molecule encoder).

    The author pipeline drops descriptor 0 and then standardises the remaining
    columns; standardisation is done by the caller because it must be fitted on
    the training drugs only.
    """

    from descriptastorus.descriptors.DescriptorGenerator import MakeGenerator

    generator = MakeGenerator(("RDKit2D",))
    unique = sorted({str(s) for s in smiles if str(s)})
    cache: dict[str, np.ndarray] = {}
    width: int | None = None
    for value in unique:
        result = generator.process(value)
        if result is None or len(result) < 2:
            if width is None:
                raise RuntimeError("RDKit2D generator returned no vector for any SMILES")
            cache[value] = np.zeros(width, dtype=np.float32)
            continue
        vector = np.asarray(result[1:], dtype=np.float32)
        vector[~np.isfinite(vector)] = 0.0
        width = vector.shape[0]
        cache[value] = vector
    if width is None:
        raise RuntimeError("no valid SMILES for the RDKit2D embedding")
    rows = np.zeros((len(smiles), width), dtype=np.float32)
    for i, value in enumerate(smiles):
        rows[i] = cache.get(str(value), np.zeros(width, dtype=np.float32))
    return rows


def fit_standardizer(matrix: np.ndarray, min_std: float = 0.01) -> dict:
    """Fit column mean/std on training drugs and drop near-constant columns."""

    mean = matrix.mean(axis=0)
    std = matrix.std(axis=0)
    keep = std > min_std
    std = np.where(keep, std, 1.0)
    return {
        "mean": mean.astype(np.float32),
        "std": std.astype(np.float32),
        "keep": keep,
    }


def apply_standardizer(matrix: np.ndarray, stats: dict) -> np.ndarray:
    keep = np.asarray(stats["keep"], dtype=bool)
    mean = np.asarray(stats["mean"], dtype=np.float32)
    std = np.asarray(stats["std"], dtype=np.float32)
    return ((matrix[:, keep] - mean[keep]) / std[keep]).astype(np.float32)


# ---------------------------------------------------------------------------
# Misc


def ensure_float32(array: np.ndarray) -> np.ndarray:
    out = np.asarray(array, dtype=np.float32)
    if not np.isfinite(out).all():
        raise RuntimeError("adapter produced a non-finite prediction")
    return out


def case_summary(case: dict[str, np.ndarray]) -> dict:
    return {
        "n_cases": int(len(case["well_id"])),
        "n_groups": int(len(np.unique(case["group_id"]))),
        "n_drugs": int(len(np.unique(case["drug_key"]))),
        "n_control_groups": int(len(np.unique(case["control_group"]))),
    }
