"""
Query FDA-approved small molecule drugs from ChEMBL REST API.
Filters for small molecules with canonical SMILES and saves as fda_drugs.csv.
"""
import requests
import csv
import time

BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule"

def main():
    drugs = []
    offset = 0
    limit = 500
    total = None

    print("Fetching approved small molecules from ChEMBL (max_phase=4)...")
    while True:
        params = {
            "max_phase": 4,
            "molecule_type": "Small molecule",
            "limit": limit,
            "offset": offset,
            "format": "json",
        }
        try:
            resp = requests.get(BASE_URL, params=params, timeout=60)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            print(f"Error at offset {offset}: {e}, retrying...")
            time.sleep(5)
            continue

        if total is None:
            total = data.get("page_meta", {}).get("total_count", 0)
            print(f"Total approved small molecules in ChEMBL: {total}")

        molecules = data.get("molecules", [])
        if not molecules:
            break

        for mol in molecules:
            struct = mol.get("molecule_structures")
            smiles = struct.get("canonical_smiles") if struct else None
            if smiles:
                drugs.append({
                    "chembl_id": mol["molecule_chembl_id"],
                    "pref_name": mol.get("pref_name", ""),
                    "canonical_smiles": smiles,
                    "first_approval": mol.get("first_approval", ""),
                })

        offset += limit
        print(f"  Processed {min(offset, total)}/{total}, collected {len(drugs)} drugs with SMILES")
        if offset >= total:
            break
        time.sleep(0.1)

    print(f"\nTotal: {len(drugs)} FDA/globally approved small molecule drugs with canonical SMILES")

    outpath = "vcms/fda/fda_drugs.csv"
    with open(outpath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["chembl_id", "pref_name", "canonical_smiles", "first_approval"])
        writer.writeheader()
        writer.writerows(drugs)

    print(f"Saved to {outpath}")

if __name__ == "__main__":
    main()
