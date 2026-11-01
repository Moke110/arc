#!/usr/bin/env python
"""Synthetic stand-in for the 04 pretrain view, used only to smoke-test the trainer.

The real pretraining input is
``stage2/04_tahoe_data/output/pretrain_view/tahoe_cm_schema.h5ad`` (67,018 x 18,164,
built by ``stage2/04_tahoe_data/execution/build_pretrain_view.py``).  That file is
produced by a peer step and does not exist while this step is being written, so this
script writes a *tiny* file with the exact same schema:

* ``X`` int32 CSR (300 x 18,164) and ``layers['log1p_cp10k']`` float32 CSR, same pattern;
* ``var`` = the frozen CM gene axis copied verbatim from
  ``stage2/02_cm_data/output/corpus/v1/genes.csv`` (``FAM153B`` stays all-zero, D2);
* ``obs`` = the 24 CM columns + the Tahoe-native extras + the three author shims,
  with 2 control groups, 40 drugs (2 of them structure-less) and plates ``plate1``
  (train) / ``plate13`` (valid);
* ``uns['control_mean']`` = the production encoding, a (2, 18,164) float64 pandas
  DataFrame indexed by ``control_group`` (spec §1.4);
* the root ``view``/``built_from``/``pseudobulk_unit`` attrs.

The fixture is a **test aid, not a deliverable**: ``*.h5ad`` is gitignored, and
``train_prnet.py`` never reads it unless ``--view`` points at it.

**Run this with the anndata 0.12 env** so the fixture carries the identical
``dataframe/0.2.0`` encoding that 04 writes.  anndata 0.10.8 (the pinned trainer
env) cannot *write* that DataFrame — HDF5 overflows the attribute holding the
18,164 column names — but it *reads* it fine, which is all the trainer needs.
The script refuses to silently downgrade; pass ``--control-mean-encoding
ndarray`` explicitly to exercise the trainer's fallback path instead.

Usage (from the repo root)::

    /home/chang/miniconda3/envs/arc/bin/python \
        stage2/05_prnet_pretrain/execution/_smoke/make_fixture.py

    # the smoke run itself uses the pinned GPU env:
    /home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python \
        stage2/05_prnet_pretrain/execution/train_prnet.py --smoke \
        --view stage2/05_prnet_pretrain/execution/_smoke/tahoe_cm_schema_smoke.h5ad
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

N_GENES = 18164
N_CONTROLS_PER_GROUP = 30
SEED = 20261021

# Plates chosen so the fixed train/valid split (plate13 + plate14 = valid) has rows
# on both sides without inventing a plate the split rule does not know about.
TRAIN_PLATE = "plate1"
VALID_PLATE = "plate13"
CONTROL_GROUPS = {
    TRAIN_PLATE: f"tahoe::{TRAIN_PLATE}::cellA",
    VALID_PLATE: f"tahoe::{VALID_PLATE}::cellB",
}
CELL_LINES = {TRAIN_PLATE: ("cellA", "Cell line A"), VALID_PLATE: ("cellB", "Cell line B")}
DOSES = (0.05, 0.5, 5.0)

# 38 structure-bearing drugs (simple, valid SMILES) + 2 structure-less ones, so both
# the fingerprint path and the `has_structure == False` drop path are exercised.
STRUCTURED = [
    ("ethanol", "CCO"),
    ("methanol", "CO"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("nicotine", "CN1CCC[C@H]1c1cccnc1"),
    ("urea", "NC(=O)N"),
    ("methane", "C"),
    ("benzene", "c1ccccc1"),
    ("toluene", "Cc1ccccc1"),
    ("phenol", "Oc1ccccc1"),
    ("aniline", "Nc1ccccc1"),
    ("acetic_acid", "CC(=O)O"),
    ("acetone", "CC(C)=O"),
    ("propane", "CCC"),
    ("butane", "CCCC"),
    ("cyclohexane", "C1CCCCC1"),
    ("pyridine", "c1ccncc1"),
    ("imidazole", "c1cnc[nH]1"),
    ("indole", "c1ccc2[nH]ccc2c1"),
    ("morpholine", "C1COCCN1"),
    ("piperidine", "C1CCNCC1"),
    ("pyrrole", "c1cc[nH]c1"),
    ("furan", "c1ccoc1"),
    ("thiophene", "c1ccsc1"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1"),
    ("glycine", "NCC(=O)O"),
    ("alanine", "C[C@@H](N)C(=O)O"),
    ("valine", "CC(C)[C@@H](N)C(=O)O"),
    ("leucine", "CC(C)C[C@@H](N)C(=O)O"),
    ("serine", "OC[C@@H](N)C(=O)O"),
    ("cysteine", "SC[C@@H](N)C(=O)O"),
    ("methionine", "CSCC[C@@H](N)C(=O)O"),
    ("lysine", "NCCCC[C@@H](N)C(=O)O"),
    ("arginine", "NC(=N)NCCCC[C@@H](N)C(=O)O"),
    ("histidine", "O=C(O)[C@@H](N)Cc1c[nH]cn1"),
]
NULL_STRUCTURE = ["nullstruct_alpha", "nullstruct_beta"]


def find_repo_root(start: Path) -> Path:
    """Walk up from *start* to the directory holding ``.git`` (never count levels)."""
    for candidate in (start, *start.parents):
        if (candidate / ".git").exists():
            return candidate
    raise SystemExit(f"cannot find the repository root: no .git at or above {start}")


def load_gene_axis(root: Path) -> pd.DataFrame:
    genes = pd.read_csv(root / "stage2/02_cm_data/output/corpus/v1/genes.csv")
    if len(genes) != N_GENES:
        raise SystemExit(f"genes.csv has {len(genes)} rows, expected {N_GENES}")
    return genes


def molecule_inchikey(smiles: str) -> str:
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchiKey(mol) if mol is not None else ""


def drug_order() -> list[str]:
    """40 drugs; the first structure-less one lands on the train plate, the other on valid."""
    names = [name for name, _ in STRUCTURED]
    return names[:6] + [NULL_STRUCTURE[0]] + names[6:] + [NULL_STRUCTURE[1]]


def build_obs(smiles_of: dict[str, str], drugs: list[str], rng: np.random.Generator) -> pd.DataFrame:
    """30 controls + 6 replicates per drug, spread over plate1 (train) and plate13 (valid)."""
    plate_of_drug = {drug: (TRAIN_PLATE if i < 25 else VALID_PLATE) for i, drug in enumerate(drugs)}
    records: list[dict] = []
    well_count = {TRAIN_PLATE: 0, VALID_PLATE: 0}

    def next_well(plate: str) -> str:
        well_count[plate] += 1
        return f"w{well_count[plate]:03d}"

    for plate in (TRAIN_PLATE, VALID_PLATE):
        line_id, line_name = CELL_LINES[plate]
        for _ in range(N_CONTROLS_PER_GROUP):
            records.append(
                dict(
                    plate=plate, cell_line_id=line_id, cell_name=line_name, well=next_well(plate),
                    compound_name="DMSO_TF", smiles="", dose_um=0.0, is_control=True,
                    control_group=CONTROL_GROUPS[plate], has_structure=True,
                    n_cells=int(rng.integers(200, 2000)),
                )
            )
        for drug in drugs:
            if plate_of_drug[drug] != plate:
                continue
            smiles = smiles_of[drug]
            for rep in range(6):
                records.append(
                    dict(
                        plate=plate, cell_line_id=line_id, cell_name=line_name,
                        well=next_well(plate), compound_name=drug, smiles=smiles,
                        dose_um=float(DOSES[rep % len(DOSES)]), is_control=False,
                        control_group=CONTROL_GROUPS[plate], has_structure=bool(smiles),
                        n_cells=int(rng.integers(200, 2000)),
                    )
                )

    obs = pd.DataFrame.from_records(records)
    obs["well"] = [f"tahoe::{r.plate}:{r.well}:{r.cell_line_id}" for r in obs.itertuples(index=False)]
    obs = obs.set_index("well")
    obs.index.name = "well"
    obs = obs.sort_index(kind="stable")

    obs["dataset_id"] = pd.Categorical(["tahoe"] * len(obs))
    obs["compound_name"] = pd.Categorical(obs["compound_name"])
    obs["has_structure"] = obs["has_structure"].astype(bool)
    obs["drug_key"] = [molecule_inchikey(s) if s else "" for s in obs["smiles"]]
    obs["inchikey"] = obs["drug_key"]  # no salts in the fixture, so parent == raw
    obs["dose_um"] = obs["dose_um"].astype(np.float64)
    obs["dose_source"] = "nominal_atlas"
    obs["time_h"] = np.nan
    obs["vehicle"] = pd.Categorical([("DMSO_TF" if c else "") for c in obs["is_control"]])
    obs["plate"] = pd.Categorical(obs["plate"])
    obs["control_group"] = obs["control_group"].astype(str)
    obs["cell_type"] = pd.Categorical([np.nan] * len(obs))
    obs["donor"] = np.nan
    obs["atc"] = ""
    obs["moa"] = ""
    obs["mesh_classes"] = ""
    obs["pharm_atc_moa"] = False
    obs["pharm_mesh"] = False
    obs["annotation_source"] = ""

    # `replicate` = within-group index, group = (control_group, drug_key, dose_um),
    # ordered by well (spec §1.3).
    obs["replicate"] = (
        obs.groupby(["control_group", "drug_key", "dose_um"], observed=True).cumcount().astype(np.int32)
    )

    obs["cell_line_id"] = pd.Categorical(obs["cell_line_id"])
    obs["cell_name"] = pd.Categorical(obs["cell_name"])
    obs["n_cells"] = obs["n_cells"].astype(np.int32)
    obs["moa_broad"] = pd.Categorical(["MOA_BROAD_A"] * len(obs))
    obs["moa_fine"] = pd.Categorical(["MOA_FINE_A"] * len(obs))
    obs["pubchem_cid"] = ["100000" if k else "" for k in obs["drug_key"]]
    obs["has_matched_control"] = True

    # Author-compatibility shims (spec §1.3), not part of CM's contract block.
    obs["SMILES"] = obs["smiles"]
    obs["dose"] = obs["dose_um"]

    first_control: dict[str, str] = {}
    for group, frame in obs[obs["is_control"]].groupby("control_group"):
        first_control[str(group)] = frame.index[0]
    obs["paired_control_index"] = [
        well if is_control else first_control[group]
        for well, is_control, group in zip(obs.index, obs["is_control"], obs["control_group"])
    ]
    return obs


def build_counts(obs: pd.DataFrame, genes: pd.Index, rng: np.random.Generator) -> np.ndarray:
    """Dense int32 counts; the CSR writer drops the zeros."""
    group_mean = {
        group: rng.gamma(shape=0.25, scale=1.8, size=N_GENES) for group in CONTROL_GROUPS.values()
    }
    drugs = list(dict.fromkeys(obs.loc[~obs["is_control"], "compound_name"].astype(str)))
    drug_factor = {drug: float(rng.uniform(0.7, 1.3)) for drug in drugs}
    gene_factor = np.exp(rng.normal(0.0, 0.3, size=N_GENES))

    counts = np.zeros((len(obs), N_GENES), dtype=np.int32)
    for i, (_, row) in enumerate(obs.iterrows()):
        mu = np.maximum(group_mean[row["control_group"]], 1e-3)
        if not bool(row["is_control"]):
            mu = mu * drug_factor[str(row["compound_name"])] * gene_factor
            mu = np.maximum(mu * (1.0 + 0.1 * np.log10(float(row["dose_um"]) + 1.0)), 1e-3)
        counts[i] = rng.poisson(mu).astype(np.int32)

    # D2: FAM153B has no Tahoe measurement, so its column is all-zero in the view.
    counts[:, genes.get_loc("FAM153B")] = 0
    return counts


def build_control_mean(counts: np.ndarray, obs: pd.DataFrame, genes: pd.Index) -> pd.DataFrame:
    groups = sorted(obs["control_group"].unique())
    values = np.zeros((len(groups), N_GENES), dtype=np.float64)
    for i, group in enumerate(groups):
        rows = np.flatnonzero(
            (obs["control_group"] == group).to_numpy() & obs["is_control"].to_numpy()
        )
        if rows.size == 0:
            raise SystemExit(f"control group {group} has no control row")
        values[i] = counts[rows].astype(np.float64).mean(axis=0)
    return pd.DataFrame(
        values, index=pd.Index(groups, name="control_group"), columns=genes.to_numpy()
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=None,
                        help="output .h5ad (default: next to this script)")
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument(
        "--control-mean-encoding",
        choices=("dataframe", "dataframe-t", "ndarray"),
        default="dataframe",
        help=(
            "how to persist uns['control_mean'].  'dataframe' is the production "
            "encoding of spec §1.4 and needs anndata >= 0.12 to be writable; "
            "'dataframe-t' is the transposed DataFrame; 'ndarray' keeps the fixture "
            "writable by the pinned 0.10.8 env and exercises the trainer fallback."
        ),
    )
    args = parser.parse_args(argv)

    root = find_repo_root(Path(__file__).resolve().parent)
    symbols = load_gene_axis(root)
    genes = pd.Index(symbols["symbol"].astype(str).to_numpy(), name="symbol")
    out = args.out or (Path(__file__).resolve().parent / "tahoe_cm_schema_smoke.h5ad")

    rng = np.random.default_rng(args.seed)
    smiles_of = dict(STRUCTURED) | {name: "" for name in NULL_STRUCTURE}
    drugs = drug_order()
    obs = build_obs(smiles_of, drugs, rng)
    counts = build_counts(obs, genes, rng)

    lib_size = counts.sum(axis=1, dtype=np.int64)
    obs["lib_size"] = lib_size
    obs["n_genes_detected"] = (counts > 0).sum(axis=1).astype(np.int32)

    X = sp.csr_matrix(counts)
    normalized = counts / np.maximum(lib_size, 1)[:, None] * 1e4
    log1p = sp.csr_matrix(np.log1p(normalized).astype(np.float32))

    control_mean = build_control_mean(counts, obs, genes)

    var = pd.DataFrame(
        {
            "entrez_id": symbols["entrez_id"].astype(str).to_numpy(),
            "ensembl_id": symbols["ensembl_id"].astype(str).to_numpy(),
            "match_type": pd.Categorical(
                symbols["match_type"], categories=["alias", "direct", "entrez"]
            ),
        },
        index=genes,
    )

    uns: dict[str, object] = {}
    if args.control_mean_encoding == "dataframe":
        # Production encoding (spec §1.4): anndata >= 0.12 writes dataframe/0.2.0.
        uns["control_mean"] = control_mean
    elif args.control_mean_encoding == "dataframe-t":
        # (18,164 x 700): the 700 column names fit in an HDF5 attribute; the 18,164
        # row index is a dataset, not an attribute.
        uns["control_mean"] = control_mean.T.copy()
    else:
        uns["control_mean"] = control_mean.to_numpy(dtype=np.float64)
        uns["control_mean_index"] = np.asarray(control_mean.index.astype(str))

    adata = ad.AnnData(X=X, obs=obs, var=var, uns=uns)
    adata.layers["log1p_cp10k"] = log1p
    try:
        adata.write_h5ad(out)
    except OSError as exc:
        if "object header" not in str(exc) or args.control_mean_encoding != "dataframe":
            raise
        raise SystemExit(
            f"anndata {getattr(ad, '__version__', '?')} cannot write the production "
            "DataFrame encoding of uns['control_mean'] (HDF5 attribute limit on the "
            "18,164 column names).\n"
            "Generate the fixture with the newer anndata env instead:\n"
            "  /home/chang/miniconda3/envs/arc/bin/python " + str(Path(__file__).resolve()) + "\n"
            "or pass --control-mean-encoding ndarray to exercise the trainer fallback."
        ) from exc

    with h5py.File(out, "r+") as handle:
        handle.attrs["view"] = "tahoe_cm_schema"
        handle.attrs["built_from"] = (
            "synthetic fixture: stage2/05_prnet_pretrain/execution/_smoke/make_fixture.py"
        )
        handle.attrs["pseudobulk_unit"] = "(plate, well, cell_line)"

    n_dropped = int(((~obs["is_control"]) & ~obs["has_structure"]).sum())
    print(f"wrote {out}")
    print(f"  anndata {getattr(ad, '__version__', '?')} | control_mean encoding "
          f"{args.control_mean_encoding} ({type(uns['control_mean']).__name__})")
    print(f"  shape {adata.shape}  X nnz {X.nnz}  log1p nnz {log1p.nnz}")
    print(f"  obs columns {len(obs.columns)}  drugs {len(drugs)} (2 structure-less)  "
          f"control groups {len(control_mean)}")
    print(f"  controls {int(obs['is_control'].sum())}  perturbed {int((~obs['is_control']).sum())}  "
          f"dropped(no structure) {n_dropped}")
    print(f"  split: {TRAIN_PLATE}={int((obs['plate'].astype(str) == TRAIN_PLATE).sum())} rows, "
          f"{VALID_PLATE}={int((obs['plate'].astype(str) == VALID_PLATE).sum())} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
