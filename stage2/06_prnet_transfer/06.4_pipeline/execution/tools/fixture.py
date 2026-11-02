"""Synthetic frozen-input fixture matching ticket #3's contract.

Test aid only, never a deliverable: it writes a tiny frozen input directory with
the exact file names / column semantics the pipeline loader expects, on the real
18,164-gene axis by default, so the CPU end-to-end test exercises the true
matrix shapes (encoder ``[18164+64] -> 128 -> 64``, decoder ``2 x 18164``) and
the frozen direction-AUPRC@200 metric instead of a toy stand-in.

Contents:

* ``cm_transfer.h5ad`` — dense ``int32`` X + dense ``float32``
  ``layers['log1p_cp10k']``; wells carry QC/eligibility/main-set flags, 5 folds
  with every main drug test once and valid once, plus eligibility-excluded
  (``main_set`` False) QC-pass wells for the #11 selector;
* ``folds.parquet`` — ``(drug_key, fold, role, plate)`` for main-set drugs;
* ``consensus_controls.npz`` / ``consensus_membership.parquet`` — per-fold c*
  built from the QC-pass DMSO wells of the plates touched by training drugs
  (per-plate mean, then equal weight across plates), with the exact membership;
* ``tahoe_transfer.h5ad`` — CSR int32 X + CSR float32 log1p;
  ``tahoe_control_groups.parquet`` (long ``group_index|group_key, row``) and
  ``tahoe_pairs.npz`` (``row``, ``group_index``); matched controls are
  ``(plate, cell_line)`` and plate13/14 are the frozen valid split;
* ``drug_smiles.csv`` — drug_key -> SMILES for the fingerprint encoder;
* ``locked_config.json`` — canonical frozen config (18,164 axis, k=200, 5 folds);
* ``manifest.json`` — sha256 of every artifact, which the loader verifies.

Run with the ``arc`` env (anndata 0.12 is the writer used by the frozen data
steps); ``tests/_support.py`` imports :func:`build_fixture` directly.
"""

from __future__ import annotations

import hashlib
import json
import zlib
import shutil
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

N_GENES_DEFAULT = 18164
SEED_DEFAULT = 20261021
TRAIN_PLATES = ("plate1", "plate2", "plate3")
VALID_PLATES = ("plate13", "plate14")
CELL_LINES = ("cellA", "cellB")
CM_PLATES = ("PLATE_A", "PLATE_B", "PLATE_C")
CM_DOSES = (0.2, 1.0, 10.0)
TAHOE_DOSES = (0.05, 0.5, 5.0)
N_FOLDS = 5


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _smiles(index: int) -> str:
    # methane, ethane, propane, ... — all valid RDKit SMILES, all distinct.
    return "C" * ((index % 40) + 1)


def _counts_log1p(counts: np.ndarray) -> np.ndarray:
    counts = counts.astype(np.float64)
    library = counts.sum(axis=1, keepdims=True)
    library[library == 0] = 1.0
    return np.log1p(counts / library * 10_000.0).astype(np.float32)


def _write_manifest(out_dir: Path, names: list[str]) -> dict[str, str]:
    hashes = {name: _sha256(out_dir / name) for name in names if name != "manifest.json"}
    (out_dir / "manifest.json").write_text(
        json.dumps(
            {
                "schema": "fixture.artifact_manifest/v1",
                "artifacts": {name: {"sha256": digest} for name, digest in hashes.items()},
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return hashes


def _cm_frame(
    n_genes: int,
    n_main_drugs: int,
    n_excluded: int,
    seed: int,
) -> tuple[ad.AnnData, pd.DataFrame, dict[str, Any]]:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []

    def add_well(**kwargs: Any) -> None:
        rows.append(kwargs)

    well_counter = 0

    def next_well(plate: str) -> str:
        """Plate-well id from a 12-column grid with spreadsheet column letters.

        ``A``..``Z`` then ``AA``..: a plain ``% 26`` wrap repeated the letter
        (and therefore the well id) once a fixture needed more than 312 wells
        per plate.  Ids stay unique for any fixture size.
        """
        nonlocal well_counter
        well_counter += 1
        index = well_counter - 1
        letters = ""
        column = index // 12
        while column >= 0:
            column, remainder = divmod(column, 26)
            letters = chr(ord("A") + remainder) + letters
            column -= 1
        return f"{plate}-{letters}{index % 12 + 1:02d}"

    # controls: per plate 4 QC-pass DMSO, 1 QC-fail DMSO, 1 QC-pass Media
    for plate in CM_PLATES:
        for i in range(6):
            is_media = i == 5
            qc = not (i == 4)
            add_well(
                well=next_well(plate),
                dataset_id="temp_env",
                compound_name="Media" if is_media else "DMSO",
                drug_key="",
                inchikey="",
                smiles="",
                dose_um=0.0,
                time_h=24.0,
                is_control=True,
                vehicle="Media" if is_media else "DMSO",
                plate=plate,
                control_group=f"temp_env::{plate}",
                cell_type="ipsc_cm",
                donor="d1",
                replicate=i,
                eligibility_decision="include",
                qc_pass=qc,
                qc_reason="" if qc else "low_library",
                main_set=True,
                replicate_raw=i,
                solvent="Media" if is_media else "DMSO",
            )

    # main-set drugs: round-robin over plates, all 3 doses x 2 replicates
    main_drugs = [f"DRUGC{i:03d}" for i in range(n_main_drugs)]
    drug_plate: dict[str, str] = {}
    for i, drug in enumerate(main_drugs):
        plate = CM_PLATES[i % len(CM_PLATES)]
        drug_plate[drug] = plate
        for dose_idx, dose in enumerate(CM_DOSES):
            for replicate in range(2):
                qc = not (i == 0 and dose_idx == 2 and replicate == 1)  # one QC-fail well
                add_well(
                    well=next_well(plate),
                    dataset_id="temp_env",
                    compound_name=drug,
                    drug_key=drug,
                    inchikey=f"INCHIKEY-{drug}",
                    smiles=_smiles(i),
                    dose_um=float(dose),
                    time_h=24.0,
                    is_control=False,
                    vehicle="DMSO",
                    plate=plate,
                    control_group=f"temp_env::{plate}",
                    cell_type="ipsc_cm",
                    donor="d1",
                    replicate=replicate,
                    eligibility_decision="include",
                    qc_pass=qc,
                    qc_reason="" if qc else "low_library",
                    main_set=True,
                    replicate_raw=replicate,
                    solvent="DMSO",
                )

    # eligibility-excluded but QC-pass wells (#11 expanded training candidates)
    excluded = [f"ENVC{i:03d}" for i in range(n_excluded)]
    excluded_plate: dict[str, str] = {}
    for i, drug in enumerate(excluded):
        plate = CM_PLATES[0]
        excluded_plate[drug] = plate
        for dose_idx, dose in enumerate((1.0, 10.0)):
            add_well(
                well=next_well(plate),
                dataset_id="temp_env",
                compound_name=drug,
                drug_key=drug,
                inchikey=f"INCHIKEY-{drug}",
                smiles=_smiles(60 + i),
                dose_um=float(dose),
                time_h=24.0,
                is_control=False,
                vehicle="DMSO",
                plate=plate,
                control_group=f"temp_env::{plate}",
                cell_type="ipsc_cm",
                donor="d1",
                replicate=dose_idx,
                eligibility_decision="exclude",
                qc_pass=True,
                qc_reason="",
                main_set=False,
                replicate_raw=dose_idx,
                solvent="DMSO",
            )

    obs = pd.DataFrame(rows)
    if obs["well"].duplicated().any():
        raise AssertionError("fixture CM well id is not unique (well grid wrapped)")
    counts = rng.poisson(lam=3.0, size=(len(obs), n_genes)).astype(np.int32)
    # Make the fixture's response signal detectable: treated wells on the first
    # few genes get a drug/dose-dependent shift, so AUPRC has something to rank.
    for i, row in obs.iterrows():
        if not row["is_control"]:
            shift = int(zlib.crc32(f"{row['drug_key']}|{float(row['dose_um']):g}".encode()) % 40) + 2
            counts[i, :shift] += (np.arange(shift) + 1)[::-1] * 2
    log1p = _counts_log1p(counts)
    obs = obs.set_index("well", drop=False)
    var = pd.DataFrame(
        {
            "entrez_id": np.arange(n_genes) + 1,
            "ensembl_id": [f"ENSG{i:011d}" for i in range(n_genes)],
            "match_type": "direct",
            "gene_measured": np.ones(n_genes, dtype=bool),
        },
        index=pd.Index([f"GENE{i:05d}" for i in range(n_genes)], name="symbol"),
    )
    adata = ad.AnnData(X=counts, obs=obs, var=var, layers={"log1p_cp10k": log1p})
    meta = {
        "main_drugs": main_drugs,
        "drug_plate": drug_plate,
        "excluded_drugs": excluded,
        "excluded_plate": excluded_plate,
        "wells": obs,
    }
    return adata, obs, meta


def _folds_frame(main_drugs: list[str], drug_plate: dict[str, str]) -> pd.DataFrame:
    records = []
    for drug_index, drug in enumerate(main_drugs):
        for fold in range(N_FOLDS):
            if drug_index % N_FOLDS == fold:
                role = "test"
            elif drug_index % N_FOLDS == (fold + 1) % N_FOLDS:
                role = "valid"
            else:
                role = "train"
            records.append(
                {"drug_key": drug, "fold": fold, "role": role, "plate": drug_plate[drug]}
            )
    return pd.DataFrame(records)


def _consensus(
    obs: pd.DataFrame, folds: pd.DataFrame, log1p: np.ndarray, n_genes: int
) -> tuple[dict[str, np.ndarray], pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    c_star = np.zeros((N_FOLDS, n_genes), dtype=np.float32)
    plate_info: dict[int, list[str]] = {}
    for fold in range(N_FOLDS):
        train_drugs = set(
            folds[(folds["fold"] == fold) & (folds["role"] == "train")]["drug_key"]
        )
        train_plates = sorted(
            {obs.loc[obs["drug_key"] == drug, "plate"].iloc[0] for drug in train_drugs}
        )
        plate_info[fold] = train_plates
        member_rows: dict[str, np.ndarray] = {}
        for plate in train_plates:
            mask = (
                obs["is_control"].to_numpy()
                & obs["qc_pass"].to_numpy()
                & (obs["vehicle"].to_numpy() == "DMSO")
                & (obs["plate"].to_numpy() == plate)
            )
            member_rows[plate] = np.flatnonzero(mask)
        n_plates = len(train_plates)
        plate_means = []
        for plate in train_plates:
            rows_idx = member_rows[plate]
            weight = 1.0 / (n_plates * len(rows_idx))
            for row in rows_idx:
                rows.append(
                    {
                        "fold": fold,
                        "plate": plate,
                        "well": obs.iloc[row]["well"],
                        "row": int(row),
                        "weight": float(weight),
                    }
                )
            plate_means.append(log1p[rows_idx].astype(np.float64).mean(axis=0))
        c_star[fold] = np.stack(plate_means).mean(axis=0).astype(np.float32)
    membership = pd.DataFrame(rows)
    return {"c_star": c_star, "folds": np.arange(N_FOLDS, dtype=np.int64)}, membership, {
        "train_plates": plate_info
    }


def _tahoe_frame(n_genes: int, seed: int) -> tuple[ad.AnnData, pd.DataFrame, dict[str, Any], dict[str, np.ndarray]]:
    rng = np.random.default_rng(seed + 1)
    rows: list[dict[str, Any]] = []
    controls: dict[tuple[str, str], list[int]] = {}

    def add_control(plate: str, line: str, physical_well: str) -> None:
        controls.setdefault((plate, line), []).append(len(rows))
        rows.append(
            {
                "observation_id": f"tahoe::{plate}:{physical_well}:{line}",
                "physical_well": physical_well,
                "plate": plate,
                "cell_line_id": line,
                "drug_key": "",
                "compound_name": "DMSO_TF",
                "smiles": "",
                "dose_um": 0.0,
                "is_control": True,
                "qc_pass": True,
                "has_structure": False,
                "tahoe_split": "train",
                "vehicle": "DMSO_TF",
            }
        )

    for plate in TRAIN_PLATES + VALID_PLATES:
        for line in CELL_LINES:
            for k in range(2):
                add_control(plate, line, f"smp_{plate}_{line}_ctl{k}")

    def add_treated(plate: str, line: str, drug: str, smiles: str, dose: float, split: str, idx: int) -> None:
        rows.append(
            {
                "observation_id": f"tahoe::{plate}:smp_{plate}_{line}_{idx}:{line}",
                "physical_well": f"smp_{plate}_{line}_{idx}",
                "plate": plate,
                "cell_line_id": line,
                "drug_key": drug,
                "compound_name": drug,
                "smiles": smiles,
                "dose_um": float(dose),
                "is_control": False,
                "qc_pass": True,
                "has_structure": True,
                "tahoe_split": split,
                "vehicle": "",
            }
        )

    counter = 0
    for p_idx, plate in enumerate(TRAIN_PLATES):
        for l_idx, line in enumerate(CELL_LINES):
            for d_idx in range(4):
                drug = f"DRUGT{(p_idx * 2 + d_idx) % 6:03d}"
                for dose in TAHOE_DOSES:
                    counter += 1
                    add_treated(plate, line, drug, _smiles(80 + (p_idx * 2 + d_idx)), dose, "train", counter)
    for p_idx, plate in enumerate(VALID_PLATES):
        for l_idx, line in enumerate(CELL_LINES):
            for d_idx in range(3):
                drug = f"DRUGTV{d_idx:03d}"
                drug_ix = 6 + d_idx
                for dose in TAHOE_DOSES:
                    counter += 1
                    add_treated(plate, line, drug, _smiles(80 + drug_ix), dose, "valid", counter)

    obs = pd.DataFrame(rows)
    if obs["observation_id"].duplicated().any():
        raise AssertionError("fixture Tahoe observation_id is not unique")
    counts = rng.poisson(lam=3.0, size=(len(obs), n_genes)).astype(np.int32)
    for i, row in obs.iterrows():
        if not row["is_control"]:
            shift = int(zlib.crc32(f"{row['drug_key']}|{row['cell_line_id']}".encode()) % 30) + 2
            counts[i, :shift] += (np.arange(shift) + 1)[::-1] * 2
    log1p = _counts_log1p(counts)
    obs = obs.set_index("observation_id", drop=False)
    var = pd.DataFrame(
        {
            "entrez_id": np.arange(n_genes) + 1,
            "ensembl_id": [f"ENSG{i:011d}" for i in range(n_genes)],
            "match_type": "direct",
            "gene_measured": np.ones(n_genes, dtype=bool),
        },
        index=pd.Index([f"GENE{i:05d}" for i in range(n_genes)], name="symbol"),
    )
    adata = ad.AnnData(
        X=sp.csr_matrix(counts), obs=obs, var=var,
        layers={"log1p_cp10k": sp.csr_matrix(log1p)},
    )

    group_rows: list[dict[str, Any]] = []
    group_index: dict[tuple[str, str], int] = {}
    for (plate, line), members in sorted(controls.items()):
        gid = len(group_index)
        group_index[(plate, line)] = gid
        for row in members:
            group_rows.append(
                {
                    "group_index": gid,
                    "group_key": f"{plate}|{line}",
                    "plate": plate,
                    "cell_line_id": line,
                    "row": int(row),
                }
            )
    groups = pd.DataFrame(group_rows)
    treated_rows = np.flatnonzero(~obs["is_control"].to_numpy())
    pair_group = np.array(
        [group_index[(obs.iloc[row]["plate"], obs.iloc[row]["cell_line_id"])] for row in treated_rows],
        dtype=np.int64,
    )
    pairs = {"row": treated_rows.astype(np.int64), "group_index": pair_group}
    smiles = {
        row["drug_key"]: row["smiles"]
        for _, row in obs[~obs["is_control"]].iterrows()
    }
    return adata, obs, {"control_groups": group_index, "smiles": smiles}, pairs


def build_fixture(
    out_dir: Path | str,
    n_genes: int = N_GENES_DEFAULT,
    n_main_drugs: int = 30,
    n_excluded: int = 3,
    seed: int = SEED_DEFAULT,
    overwrite: bool = True,
) -> dict[str, Any]:
    """Write the synthetic frozen input to *out_dir*; return its metadata."""
    out_dir = Path(out_dir)
    if out_dir.exists():
        if not overwrite:
            raise FileExistsError(f"fixture directory already exists: {out_dir}")
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)

    cm, obs, meta = _cm_frame(n_genes, n_main_drugs, n_excluded, seed)
    folds = _folds_frame(meta["main_drugs"], meta["drug_plate"])
    consensus, membership, consensus_meta = _consensus(
        obs, folds, cm.layers["log1p_cp10k"], n_genes
    )
    tahoe, tahoe_obs, tahoe_meta, pairs = _tahoe_frame(n_genes, seed)

    cm.write_h5ad(out_dir / "cm_transfer.h5ad")
    tahoe.write_h5ad(out_dir / "tahoe_transfer.h5ad")
    folds.to_parquet(out_dir / "folds.parquet", index=False)
    membership.to_parquet(out_dir / "consensus_membership.parquet", index=False)
    np.savez(out_dir / "consensus_controls.npz", **consensus)
    pd.read_parquet(out_dir / "consensus_membership.parquet")  # parquet magic check
    # write the long-form control groups straight from the constructed mapping
    control_rows = []
    for (plate, line), gid in sorted(tahoe_meta["control_groups"].items(), key=lambda kv: kv[1]):
        mask = (
            tahoe_obs["is_control"].to_numpy()
            & (tahoe_obs["plate"].to_numpy() == plate)
            & (tahoe_obs["cell_line_id"].to_numpy() == line)
        )
        for row in np.flatnonzero(mask):
            control_rows.append(
                {
                    "group_index": gid,
                    "group_key": f"{plate}|{line}",
                    "plate": plate,
                    "cell_line_id": line,
                    "row": int(row),
                }
            )
    pd.DataFrame(control_rows).to_parquet(out_dir / "tahoe_control_groups.parquet", index=False)
    np.savez(out_dir / "tahoe_pairs.npz", **pairs)

    smiles_map: dict[str, str] = {}
    for i, drug in enumerate(meta["main_drugs"]):
        smiles_map[drug] = _smiles(i)
    for i, drug in enumerate(meta["excluded_drugs"]):
        smiles_map[drug] = _smiles(60 + i)
    smiles_map.update(tahoe_meta["smiles"])
    pd.DataFrame(
        {"drug_key": sorted(smiles_map), "smiles": [smiles_map[k] for k in sorted(smiles_map)]}
    ).to_csv(out_dir / "drug_smiles.csv", index=False)

    locked = {
        "structure": {
            "n_genes": int(n_genes),
            "hidden": [128],
            "z_dim": 64,
            "c_dim": 64,
            "n_dim": 10,
            "drug_dim": 1024,
            "adaptor": [128],
            "dropout": 0.05,
        },
        "dose_encoding": {
            "fingerprint": "FCFP4", "n_bits": 1024, "radius": 2,
            "use_features": True, "dose_scale": "log10(dose_uM+1)",
        },
        "loss": "gaussian_nll",
        "optimizer": "adam",
        "recipe": {
            "batch_size": 32, "max_epochs": 2, "min_epochs": 1,
            "early_stop_patience": 2, "sched_factor": 0.5, "sched_patience": 1,
            "min_lr": 1e-8, "weight_decay": 1e-8,
        },
        "metric": {"name": "direction_auprc", "k": min(200, n_genes)},
        "n_folds": N_FOLDS,
        "seeds": [seed, seed + 1, seed + 2],
        "lrs": {"pretrain": 1e-3, "cm": [1e-4, 3e-4, 1e-3]},
        "prediction_noise": "zero",
    }
    (out_dir / "locked_config.json").write_text(
        json.dumps(locked, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    artifact_names = [
        "cm_transfer.h5ad", "folds.parquet", "consensus_controls.npz",
        "consensus_membership.parquet", "tahoe_transfer.h5ad",
        "tahoe_control_groups.parquet", "tahoe_pairs.npz",
        "drug_smiles.csv", "locked_config.json",
    ]
    hashes = _write_manifest(out_dir, artifact_names)
    return {
        "dir": str(out_dir),
        "n_genes": int(n_genes),
        "n_folds": N_FOLDS,
        "main_drugs": meta["main_drugs"],
        "excluded_drugs": meta["excluded_drugs"],
        "tahoe": {"control_groups": tahoe_meta["control_groups"]},
        "consensus": consensus_meta,
        "artifact_hashes": hashes,
    }


def _tahoe_controls(obs: pd.DataFrame) -> dict[tuple[str, str], list[int]]:
    out: dict[tuple[str, str], list[int]] = {}
    mask = obs["is_control"].to_numpy()
    for row in np.flatnonzero(mask):
        key = (obs.iloc[row]["plate"], obs.iloc[row]["cell_line_id"])
        out.setdefault(key, []).append(int(row))
    return out


def load_fixture_meta(data_dir: Path | str) -> dict[str, Any]:
    """Re-open a fixture and return the identities tests need."""
    import anndata as ad

    data_dir = Path(data_dir)
    cm = ad.read_h5ad(data_dir / "cm_transfer.h5ad")
    folds = pd.read_parquet(data_dir / "folds.parquet")
    membership = pd.read_parquet(data_dir / "consensus_membership.parquet")
    consensus = np.load(data_dir / "consensus_controls.npz")
    tahoe = ad.read_h5ad(data_dir / "tahoe_transfer.h5ad")
    return {
        "cm_obs": cm.obs,
        "folds": folds,
        "membership": membership,
        "consensus": consensus,
        "tahoe_obs": tahoe.obs,
    }
