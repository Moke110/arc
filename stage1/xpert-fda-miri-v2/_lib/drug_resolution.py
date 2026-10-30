"""
Authoritative CM-compound -> L1000 drug-index resolution (workflow-owned).

Why this exists
---------------
Every sample that enters the CM finetune-data union must index a valid row of
the frozen drug-feature arrays (UniMol `input/all_drugs_unimol_arr.npy` and
the author HG reference), i.e. an L1000 drug index in `[0, 8981)`. Compounds
that resolve only to FDA-extension rows (>= 8981) or not at all are excluded
from the union and recorded in the step-03 manifest, so the finetune set is
auditable and no sample silently trains on a mismatched drug.

Resolution rule (author + dataset data ONLY, no external database)
------------------------------------------------------------------
  (a) exact canonical-SMILES match where the dataset or the DToxS small-molecule
      metadata provides a SMILES, else
  (b) drug-name match against the author L1000 drug-name table (`cmap_name` in
      the L1000 pretrain h5ad obs).

SMILES are canonicalized with `isomericSmiles=False` (connectivity only), which
is the v1-convention: DToxS's stored SMILES carry stereochemistry that the
L1000 table records in separate stereo/salt rows (e.g. sunitinib -> 2002, not
the Z-isomer row 9037). Name matching normalizes to lowercase alphanumerics
(e.g. `Acetylsalicylic_acid` -> `acetylsalicylicacid`).

Data sources
-----------
  1. author idx->SMILES mapping   workflow input/all_drugs_idx2smi.npy
  2. author L1000 drug names      cmap_name + pert_idx in the L1000 pretrain
     h5ad (pert_data/XPert-L1000/source/data/l1000_mdmt_full_336852.h5ad)
  3. DToxS small-molecule metadata pert_data/LINCS-LDS-1593/source/Metadata/
     Small_Molecule_Metadata.txt (sm_center_batch_id code, sm_name,
     sm_smiles_parent) — read-only. This is the released LINCS LDS-1593
     package (the GSE raw GEO deposit was deleted 2026-08-18; the v1-era
     package is frozen at dtoxs_primary/v1-frozen/).

Output
------
  03_cm_finetune_data/output/drug_manifest.csv
    columns: dataset, compound_name, compound_code (or ''), resolution
             (smiles | name | unmatched), drug_idx (int; -1 if unmatched),
             hg_class (in_hg | extension | unmatched), n_samples
  Only compounds with hg_class == in_hg enter the union h5ad.
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (WORKFLOW_ROOT, CM_DATA_ROOT, DTOXS_SMALL_MOL,  # noqa: E402
                   N_L1000_DRUGS, L1000_FULL_H5AD)

SMI_MAPPING = WORKFLOW_ROOT / "input" / "all_drugs_idx2smi.npy"
MANIFEST = (WORKFLOW_ROOT / "03_cm_finetune_data" / "output"
            / "drug_manifest.csv")


def _canon(smiles):
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def _norm_name(name):
    return "".join(ch for ch in str(name).lower() if ch.isalnum())


_SALT_SUFFIXES = [
    "hydrobromide", "hydrochloride", "sulfate", "sulphate", "phosphate",
    "maleate", "fumarate", "hydrogensulfate", "tosylate", "mesylate",
    "citrate", "oxalate", "sodium", "potassium", "dihydrate", "hydrate",
    "monohydrate",
]


def _name_variants(name):
    """Candidate normalized names for one compound (documented normalization).

    Tries the exact normalized name, then conservative same-compound
    normalizations: salt/solvate suffix stripping (e.g. Citalopram
    hydrobromide -> citalopram) and a leading racemic `D,L-` prefix strip
    (e.g. `D,1 Sotalol` -> sotalol). These are pure name normalization, not
    external alias data. Distinct compounds are never merged (e.g. the
    N-acetyl-procainamide metabolite is NOT mapped to procainamide).
    """
    base = _norm_name(name)
    variants = [base]
    low = str(name).lower()
    for suf in _SALT_SUFFIXES:
        if low.endswith(suf):
            stripped = _norm_name(low[: -len(suf)])
            if stripped and stripped not in variants:
                variants.append(stripped)
            break
    m = re.match(r"^d[,.]?1(?:[-\s])(.+)$", low)
    if m:
        stripped = _norm_name(m.group(1))
        if stripped and stripped not in variants:
            variants.append(stripped)
    return variants


def build_author_drug_index():
    """(smi2idx, name2idx) author drug indexes.

    smi2idx:  canonical SMILES (non-isomeric) -> smallest L1000 drug_idx.
    name2idx: normalized cmap_name -> list of L1000 drug_idx (smallest first).
    """
    idx2smi = np.load(SMI_MAPPING, allow_pickle=True).item()
    smi2idx = {}
    for idx, smi in idx2smi.items():
        c = _canon(smi)
        if c is not None:
            smi2idx.setdefault(c, int(idx))

    import scanpy as sc
    adata = sc.read_h5ad(L1000_FULL_H5AD, backed="r")
    obs = adata.obs[["cmap_name", "pert_idx"]].drop_duplicates("pert_idx")
    obs = obs[obs["cmap_name"].astype(str) != "nan"]
    name2idx = {}
    for c, i in zip(obs["cmap_name"], obs["pert_idx"]):
        name2idx.setdefault(_norm_name(c), []).append(int(i))
    for k in name2idx:
        name2idx[k] = sorted(set(name2idx[k]))
    return smi2idx, name2idx


def dtoxs_metadata_smiles_by_name():
    """DToxS small-molecule metadata: normalized sm_name -> canonical SMILES."""
    meta = pd.read_csv(DTOXS_SMALL_MOL, sep="\t")
    out = {}
    for name, smi in zip(meta["sm_name"], meta["sm_smiles_parent"]):
        if pd.isna(name) or pd.isna(smi):
            continue
        c = _canon(smi)
        if c is not None:
            out.setdefault(_norm_name(name), c)
    return out


def dtoxs_code_meta():
    """DToxS small-molecule metadata: sm_center_batch_id code -> (sm_name, canonical SMILES).

    Used by the DToxS arms (dtoxs_primary / dtoxs_ipsc), whose designs record
    the compound as the 3-letter batch code.
    """
    meta = pd.read_csv(DTOXS_SMALL_MOL, sep="\t")
    out = {}
    for _, r in meta.iterrows():
        code = r.get("sm_center_batch_id")
        if pd.isna(code) or pd.isna(r["sm_smiles_parent"]):
            continue
        c = _canon(r["sm_smiles_parent"])
        out.setdefault(str(code), (str(r["sm_name"]), c))
    return out


def resolve_compound(name, smiles, smi2idx, name2idx, meta_smi_by_name=None):
    """Resolve one compound to an L1000 drug index.

    Returns (drug_idx, method, note) where method is 'smiles' | 'name' and
    drug_idx is -1 when the compound cannot be resolved. `note` records
    ambiguities (multiple name hits) for the manifest.
    """
    if smiles:
        c = _canon(smiles)
        if c in smi2idx:
            return smi2idx[c], "smiles", ""
    # name path: prefer a DToxS-metadata SMILES for the same compound (author
    # + dataset provenance, exact canonical match), then the author name table.
    if meta_smi_by_name is not None:
        c = meta_smi_by_name.get(_norm_name(name))
        if c is not None and c in smi2idx:
            return smi2idx[c], "smiles", "dtxos-metadata-smiles"
    hits = []
    for variant in _name_variants(name):
        hits = name2idx.get(variant, [])
        if hits:
            break
    if hits:
        note = f"ambiguous-names={hits}" if len(hits) > 1 else ""
        return hits[0], "name", note
    return -1, "unmatched", ""


def classify(idx):
    if 0 <= idx < N_L1000_DRUGS:
        return "in_hg"
    if idx >= N_L1000_DRUGS:
        return "extension"
    return "unmatched"


def build_drug_manifest(compounds):
    """compounds: list of (dataset, compound_name, compound_code, smiles, n_samples).

    Returns the manifest DataFrame and writes it to the step-03 output.
    """
    smi2idx, name2idx = build_author_drug_index()
    meta_smi = dtoxs_metadata_smiles_by_name()
    rows = []
    for dataset, name, code, smiles, n in compounds:
        idx, method, note = resolve_compound(name, smiles, smi2idx, name2idx,
                                             meta_smi)
        rows.append({
            "dataset": dataset,
            "compound_name": name,
            "compound_code": code or "",
            "resolution": method,
            "drug_idx": idx,
            "hg_class": classify(idx),
            "note": note,
            "n_samples": n,
        })
    manifest = pd.DataFrame(rows).sort_values(
        ["dataset", "compound_name"]).reset_index(drop=True)
    return manifest


def write_manifest(compounds, out_path=None):
    manifest = build_drug_manifest(compounds)
    out_path = Path(out_path or MANIFEST)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(out_path, index=False)
    print(f"Drug manifest ({len(manifest)} compounds) -> {out_path}")
    return manifest
