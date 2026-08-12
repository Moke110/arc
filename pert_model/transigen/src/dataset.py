from __future__ import annotations

from pathlib import Path
import pickle

import numpy as np
from torch.utils.data import Dataset

from .utils import load_from_HDF, subsetDict


class TranSiGenDataset(Dataset):
    def __init__(
        self,
        LINCS_index: np.ndarray,
        mol_feature_type: str,
        mid: np.ndarray,
        cid: np.ndarray,
        *,
        data_root: str | Path,
        processed_data_file: str = "processed_data.h5",
        mid2smi_file: str = "mid2smi.pickle",
        molecule_emb_file: str | None = None,
        molecule_emb_root: str | Path = "mol_emb",
    ) -> None:
        self.LINCS_index = LINCS_index
        self.mol_feature_type = mol_feature_type
        self.mid = mid
        self.cid = cid

        data_root = Path(data_root)
        self.LINCS_data = load_from_HDF(str(data_root / processed_data_file))
        with (data_root / mid2smi_file).open("rb") as f:
            self.mid2smi = pickle.load(f)

        emb_root = Path(molecule_emb_root)
        if molecule_emb_file is not None:
            emb_file = emb_root / molecule_emb_file
        elif self.mol_feature_type == "ECFP4":
            emb_file = emb_root / "example_ECFP4.pickle"
        elif self.mol_feature_type == "KPGT":
            emb_file = emb_root / "example_KPGT.pickle"
        else:
            raise ValueError(f"Unsupported mol_feature_type: {self.mol_feature_type}")

        with emb_file.open("rb") as f:
            self.smi2emb = pickle.load(f)

    def __getitem__(self, index: int):
        sub = subsetDict(self.LINCS_data, self.LINCS_index[index])
        mol_feature = self.smi2emb[self.mid2smi[self.mid[index]]]
        return sub["x1"], sub["x2"], mol_feature, self.mid[index], self.cid[index], sub["sig"]

    def __len__(self) -> int:
        return int(self.mid.shape[0])


class TranSiGenDatasetScreening(Dataset):
    def __init__(self, x1: np.ndarray, mol_feature: np.ndarray, mid: np.ndarray, cid: np.ndarray) -> None:
        self.x1 = x1
        self.mol_feature = mol_feature
        self.mid = mid
        self.cid = cid

    def __getitem__(self, index: int):
        return self.x1[index], self.mol_feature[index], self.mid[index], self.cid[index]

    def __len__(self) -> int:
        return int(self.mol_feature.shape[0])
