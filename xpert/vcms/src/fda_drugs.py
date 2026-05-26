"""
Identify FDA-approved drugs in the 8981-drug embedding library by cross-referencing
with PubChem. Output: fda_drugs.csv (drug_name, embedding_index, canonical_smiles).

Strategy:
 1. Discover the HNID for "FDA Approved Drugs" via PubChem classification API.
 2. Fetch all FDA-approved CIDs using that HNID.
 3. Batch-fetch SMILES for those CIDs.
 4. RDKit-canonicalize and match against our 8981 SMILES.
"""
import sys
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from rdkit import Chem

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

OUTPUT_DIR = PROJECT_ROOT / "vcms" / "fda"
CACHE_DIR = OUTPUT_DIR / ".claude"
FDA_CSV = OUTPUT_DIR / "fda_drugs.csv"
CHECKPOINT = CACHE_DIR / "fda_checkpoint.json"

SMI_PATH = PROJECT_ROOT / "drug_embeddings" / "all_drugs_idx2smi_8981.npy"
META_PATH = PROJECT_ROOT / "metadata" / "l1000_sdst_drugs_info_8432.csv"

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
REQUEST_TIMEOUT = 30
BATCH_SIZE = 200
RETRY_SLEEP = 2

# Known: aspirin is FDA-approved, so we can discover the HNID from its classification
REFERENCE_DRUG = "aspirin"
FDA_HNID = 5614020  # Verified 2026-05-25: HNID for "FDA Approved Drugs" in PubChem


def load_checkpoint():
    if CHECKPOINT.exists():
        with open(CHECKPOINT) as f:
            return json.load(f)
    return {"fda_cids": [], "fda_smiles": {}, "matched_indices": {}, "stage": "init"}


def save_checkpoint(state):
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(CHECKPOINT, "w") as f:
        json.dump(state, f, default=str)


def pubchem_get(url, description=""):
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=REQUEST_TIMEOUT)
            if r.status_code == 200:
                return r
            if r.status_code in (429, 503):
                wait = RETRY_SLEEP * (2 ** attempt)
                print(f"  Rate limited ({description}), waiting {wait}s...")
                time.sleep(wait)
                continue
            print(f"  HTTP {r.status_code} for {description}")
            return None
        except requests.RequestException as e:
            print(f"  Request error ({description}), attempt {attempt + 1}: {e}")
            time.sleep(RETRY_SLEEP * (2 ** attempt))
    return None


# ---------------------------------------------------------------------------
# Stage 1: Fetch FDA-approved CIDs via HNID
# ---------------------------------------------------------------------------

def fetch_fda_cids(state):
    if state.get("fda_cids"):
        print(f"  {len(state['fda_cids'])} FDA CIDs already cached, skipping fetch.")
        return state

    print("Fetching FDA-approved compound CIDs from PubChem...")
    url = f"{PUBCHEM_BASE}/classification/hnid/{FDA_HNID}/compounds/JSON"
    r = pubchem_get(url, "FDA compounds via HNID")
    if r is None:
        print("  ERROR: Could not fetch FDA CIDs.")
        return state

    data = r.json()
    cids = data.get("IdentifierList", {}).get("CID", [])
    state["fda_cids"] = cids
    state["stage"] = "cids_fetched"
    save_checkpoint(state)
    print(f"  Total FDA-approved CIDs: {len(cids)}")
    return state


# ---------------------------------------------------------------------------
# Stage 2: Batch-fetch SMILES for FDA CIDs
# ---------------------------------------------------------------------------

def fetch_smiles_for_cids(state):
    fda_cids = state.get("fda_cids", [])
    if not fda_cids:
        return state

    already_fetched = set(int(c) for c in state["fda_smiles"].keys())
    remaining = [c for c in fda_cids if int(c) not in already_fetched]

    if not remaining:
        print(f"  All {len(fda_cids)} CIDs already have SMILES cached.")
        return state

    print(f"  Fetching SMILES for {len(remaining)} CIDs in batches of {BATCH_SIZE}...")

    for i in range(0, len(remaining), BATCH_SIZE):
        batch = remaining[i : i + BATCH_SIZE]
        cid_str = ",".join(str(c) for c in batch)
        url = f"{PUBCHEM_BASE}/compound/cid/{cid_str}/property/CanonicalSMILES/JSON"
        r = pubchem_get(url, f"SMILES batch {i // BATCH_SIZE + 1}")
        if r is None:
            continue

        for prop in r.json().get("PropertyTable", {}).get("Properties", []):
            cid = prop.get("CID")
            # PubChem may return ConnectivitySMILES or CanonicalSMILES key
            smi = prop.get("CanonicalSMILES") or prop.get("ConnectivitySMILES", "")
            if cid and smi:
                state["fda_smiles"][str(cid)] = smi

        state["stage"] = "smiles_fetching"
        save_checkpoint(state)
        done = min(i + BATCH_SIZE, len(remaining))
        print(f"    {done}/{len(remaining)}")
        time.sleep(0.15)

    state["stage"] = "smiles_fetched"
    save_checkpoint(state)
    print(f"  Total FDA SMILES fetched: {len(state['fda_smiles'])}")
    return state


# ---------------------------------------------------------------------------
# Stage 3: Match FDA SMILES against our embedding library
# ---------------------------------------------------------------------------

def load_our_drugs():
    print("Loading embedding library SMILES...")
    smi_dict = np.load(SMI_PATH, allow_pickle=True).item()
    our_canonical = {}
    for idx, smi in smi_dict.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            can = Chem.MolToSmiles(mol, canonical=True)
            our_canonical[can] = idx
    print(f"  Loaded {len(smi_dict)} SMILES, {len(our_canonical)} valid canonical")

    meta_df = pd.read_csv(META_PATH)
    idx_to_name = {}
    for _, row in meta_df.iterrows():
        node_idx = int(row["drug_node_idx"])
        if node_idx not in idx_to_name:
            name = row.get("cmap_name", "")
            if pd.isna(name) or name == "" or str(name).startswith("BRD-"):
                name = str(row.get("pert_id", f"Drug-{node_idx}"))
            idx_to_name[node_idx] = str(name)

    return smi_dict, our_canonical, idx_to_name


def match_fda_drugs(state):
    print("Matching FDA-approved SMILES against embedding library...")
    smi_dict, our_canonical, idx_to_name = load_our_drugs()

    fda_smiles = state.get("fda_smiles", {})
    print(f"  Cross-referencing {len(fda_smiles)} FDA SMILES...")

    matched = {}
    for cid, fda_smi in fda_smiles.items():
        mol = Chem.MolFromSmiles(fda_smi)
        if mol is None:
            continue
        can = Chem.MolToSmiles(mol, canonical=True)
        if can in our_canonical:
            emb_idx = int(our_canonical[can])
            if emb_idx not in matched:
                drug_name = idx_to_name.get(emb_idx, f"Drug-{emb_idx}")
                matched[str(emb_idx)] = {
                    "drug_name": drug_name,
                    "canonical_smiles": smi_dict[emb_idx],
                }

    state["matched_indices"] = matched
    state["stage"] = "matched"
    save_checkpoint(state)
    print(f"  Matched {len(matched)} FDA-approved drugs in embedding library")
    return state


# ---------------------------------------------------------------------------
# Stage 4: Write output CSV
# ---------------------------------------------------------------------------

def write_output(state):
    matched = state.get("matched_indices", {})
    if not matched:
        print("WARNING: No FDA-approved drugs found. Writing empty fda_drugs.csv.")

    rows = []
    for emb_idx, info in sorted(matched.items(), key=lambda x: int(x[0])):
        rows.append({
            "drug_name": info["drug_name"],
            "embedding_index": int(emb_idx),
            "canonical_smiles": info["canonical_smiles"],
        })

    df = pd.DataFrame(rows, columns=["drug_name", "embedding_index", "canonical_smiles"])
    df.to_csv(FDA_CSV, index=False)
    print(f"\nSaved {len(df)} FDA-approved drugs to {FDA_CSV}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    state = load_checkpoint()
    stage = state.get("stage", "init")
    print(f"Starting from stage: {stage}")

    if stage in ("init", "cids_fetched"):
        state = fetch_fda_cids(state)

    if stage in ("init", "cids_fetched", "smiles_fetching"):
        state = fetch_smiles_for_cids(state)

    if stage != "matched":
        state = match_fda_drugs(state)

    write_output(state)


if __name__ == "__main__":
    main()
