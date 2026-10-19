#!/usr/bin/env python3
"""
Recover `gse244740` compound SMILES into a bundled frozen input asset.

The GSE244740 series matrix records compound *names* only (no SMILES), and the
existing union build left this dataset's SMILES empty. The benchmark needs drug
structures for the fingerprint baselines (`ridge-fp`, `nearest-drug`) and for
chemistry-aware models, so the names are resolved once, here, and frozen as
`input/gse244740_smiles.tsv`.

Resolution order per compound:
  1. the DToxS small-molecule metadata (author + dataset provenance), by
     punctuation/case-insensitive name match;
  2. PubChem PUG-REST `compound/name/<name>/property/ConnectivitySMILES`.

This is an asset-generation step, not part of the benchmark run, so the
benchmark itself stays offline and deterministic. Unresolved compounds are
written with an empty SMILES and reported.

Usage:
    python _lib/fetch_gse244740_smiles.py [--out input/gse244740_smiles.tsv]
"""

import argparse
import gzip
import json
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import CM_DATA_ROOT, INPUT_DIR, dataset_path  # noqa: E402

PUB = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/ConnectivitySMILES/JSON"
DEFAULT_OUT = INPUT_DIR / "gse244740_smiles.tsv"

# Manual query overrides where the series-matrix token does not resolve directly.
NAME_QUERY = {
    "3_Azido_3_Deoxythymidine": "3'-Azido-3'-deoxythymidine",
    "5_Fluorouracil": "5-Fluorouracil",
    "Mitomycin_C": "Mitomycin C",
    "Amphotericin_B": "Amphotericin B",
    "Acetylsalicylic_acid": "Acetylsalicylic acid",
}


def norm(name):
    return "".join(ch for ch in str(name).lower() if ch.isalnum())


def series_drugs():
    path = dataset_path("gse244740") / "GSE244740_series_matrix.txt.gz"
    text = gzip.open(path, "rt").read()
    line = [l for l in text.splitlines()
            if l.startswith("!Sample_characteristics_ch1") and "treatment" in l][0]
    names = {x.split(",")[0].strip() for x in re.findall(r'"treatment: ([^"]*)"', line)}
    return sorted(n for n in names if n not in ("DMSO", "NoTreatment"))


def dtox_smiles_by_name():
    from rdkit import Chem

    meta = pd.read_csv(dataset_path("dtoxs_primary") / "Metadata"
                       / "Small_Molecule_Metadata.txt", sep="\t")
    out = {}
    for name, smi in zip(meta["sm_name"], meta["sm_smiles_parent"]):
        if pd.isna(name) or pd.isna(smi):
            continue
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            out.setdefault(norm(name), Chem.MolToSmiles(m, isomericSmiles=False))
    return out


def pubchem_smiles(name):
    import urllib.parse

    q = NAME_QUERY.get(name, name.replace("_", " "))
    url = PUB.format(urllib.parse.quote(q))
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 404:
                return None
            r.raise_for_status()
            props = r.json()["PropertyTable"]["Properties"][0]
            return props.get("ConnectivitySMILES") or props.get("CanonicalSMILES")
        except Exception:  # noqa: BLE001
            time.sleep(2 * (attempt + 1))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    args = ap.parse_args()
    dtox = dtox_smiles_by_name()
    rows = []
    for name in series_drugs():
        smi, source = None, ""
        if norm(name) in dtox:
            smi, source = dtox[norm(name)], "dtox-metadata"
        else:
            smi = pubchem_smiles(name)
            source = "pubchem" if smi else "unresolved"
        rows.append({"drug_name": name, "smiles": smi or "", "source": source})
        print(f"  {name}: {source}")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    n = sum(1 for r in rows if r["smiles"])
    print(f"Wrote {out} ({n}/{len(rows)} resolved)")


if __name__ == "__main__":
    main()
