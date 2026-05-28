"""
Prepare FDA drug SMILES CSV for KPGT preprocessing.
Creates a minimal CSV with SMILES and a dummy target column.
"""
import pandas as pd
import os
import shutil

FDA_DIR = "vcms/fda"
KPGT_DATA = "/home/chang/arc/kpgt/feature_extract"

def main():
    df = pd.read_csv(f"{FDA_DIR}/fda_drugs.csv")
    smiles_list = df['canonical_smiles'].dropna().tolist()
    smiles_list = [s.strip() for s in smiles_list if s.strip()]

    # Create KPGT-compatible CSV (needs 'smiles' column + at least one task column)
    kpgt_df = pd.DataFrame({
        'smiles': smiles_list,
        'dummy_target': 0,
    })

    dataset_name = "FDA_drugs"
    os.makedirs(f"{KPGT_DATA}/{dataset_name}", exist_ok=True)
    kpgt_df.to_csv(f"{KPGT_DATA}/{dataset_name}/{dataset_name}.csv", index=False)
    print(f"Created {KPGT_DATA}/{dataset_name}/{dataset_name}.csv with {len(kpgt_df)} drugs")

    # Also copy to vcms/fda for reference
    kpgt_df.to_csv(f"{FDA_DIR}/fda_kpgt_input.csv", index=False)

if __name__ == "__main__":
    main()
