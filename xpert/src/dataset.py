import numpy as np
import torch
from scipy import sparse
from torch import tensor
from torch.utils.data import Dataset
from tqdm import tqdm


def _digitize(x: np.ndarray, bins: np.ndarray, side="both") -> np.ndarray:
    """
    Digitize the data into bins. This method spreads data uniformly when bins
    have same values.

    Args:

    x (:class:`np.ndarray`):
        The data to digitize.
    bins (:class:`np.ndarray`):
        The bins to use for digitization, in increasing order.
    side (:class:`str`, optional):
        The side to use for digitization. If "one", the left side is used. If
        "both", the left and right side are used. Default to "one".

    Returns:

    :class:`np.ndarray`:
        The digitized data.
    """

    left_digits = np.digitize(x, bins)
    if side == "one":
        return left_digits

    right_difits = np.digitize(x, bins, right=True)

    rands = np.random.rand(len(x))  # uniform random numbers

    digits = rands * (right_difits - left_digits) + left_digits
    digits = np.ceil(digits).astype(np.int64)
    return digits



def assign_pert_dose_corrected(pert_dose):
    if 0 <= pert_dose < 0.21:
        return 0
    elif 0.21 <= pert_dose < 0.41:
        return 1
    elif 0.41 <= pert_dose < 0.71:
        return 2
    elif 0.71 <= pert_dose < 1.01:
        return 3
    elif 1.01 <= pert_dose < 1.51:
        return 4
    elif 1.51 <= pert_dose < 3.1:
        return 5
    elif 3.1 <= pert_dose < 4.1:
        return 6
    elif 4.1 <= pert_dose < 7.1:
        return 7
    elif 7.1 <= pert_dose < 12.1:
        return 8
    elif pert_dose >= 12.1:
        return 9
    else:
        print(f'{pert_dose} is out of boundary!')
        return None


class MyDataset(Dataset):
    def __init__(self, raw_data, drug_feat, args, config, logger, max_value=None, min_value=None):

        self.raw_data = raw_data
        self.args = args
        self.config = config
        self.drug_feat = drug_feat
        self.logger = logger

        raw_data_items = self.raw_data.obs.copy()
        if 'pert_dose_idx' not in raw_data_items.columns:
            raw_data_items['pert_dose_idx'] = raw_data_items['pert_dose'].astype(np.float32).apply(assign_pert_dose_corrected)
        if 'pert_time_idx' not in raw_data_items.columns:
            pert_time2idx = {3: 0, 6: 1, 24: 2, 3.0: 0, 6.0: 1, 24.0: 2}
            raw_data_items['pert_time_idx'] = raw_data_items['pert_time'].astype(np.float32).map(pert_time2idx)

        if self.args.dataset == 'transigen_sdst':
            raw_data_items['pert_idx'] = raw_data_items['pert_id'].astype(np.int64)

        self.raw_data_items = raw_data_items

        self.trt_raw = self.raw_data.X
        if sparse.issparse(self.trt_raw):
            self.trt_raw = self.trt_raw.toarray()

        if self.args.dataset == 'transigen_sdst':
            self.ctl_raw = self.raw_data.uns['ctl_mode_avg']
        else:
            self.ctl_raw = self.raw_data.obsm['X_ctl']
        if sparse.issparse(self.ctl_raw):
            self.ctl_raw = self.ctl_raw.toarray()

        n_bins = config['dataset']['n_bins']
        self.logger.info(f'binning process... n_bins:{n_bins}')

        bins = np.quantile(np.array([min_value, max_value]), np.linspace(0, 1, n_bins - 1))
        self.trt_raw_binned = _digitize(self.trt_raw, bins, side="one").astype(np.uint8)

        if self.args.dataset == 'transigen_sdst':
            ctl_raw_binned_data = _digitize(np.array(list(self.ctl_raw.values())), bins, side="one").astype(np.uint8)
            self.ctl_raw_binned = {k: v for k, v in zip(self.ctl_raw.keys(), ctl_raw_binned_data)}
        else:
            self.ctl_raw_binned = _digitize(self.ctl_raw, bins, side="one").astype(np.uint8)

        # Pre-extract metadata columns as numpy arrays for fast __getitem__ access
        items = self.raw_data_items
        self._pert_id = items['pert_id'].values
        self._pert_idx = items['pert_idx'].values
        self._pert_dose = items['pert_dose'].values
        self._pert_dose_idx = items['pert_dose_idx'].values
        self._pert_time = items['pert_time'].values
        self._pert_time_idx = items['pert_time_idx'].values
        if 'cdsdb' in self.args.dataset:
            self._cell_id = items['Cancer subtype'].values
        else:
            self._cell_id = items['cell_iname'].values
        self._cell_idx = items['cell_idx'].values
        self._tissue_idx = items['tissue_idx'].values


    def __len__(self):
        return len(self.raw_data_items)


    def __getitem__(self, index):
        idx = index

        # get trt_raw_data and ctl_raw_data
        trt_raw_data = self.trt_raw[idx, :]
        if isinstance(trt_raw_data, sparse.csr_matrix):
            trt_raw_data = trt_raw_data.toarray().squeeze()

        trt_raw_data_binned = self.trt_raw_binned[idx, :]

        if self.args.dataset == 'transigen_sdst':
            ctl_id = self.raw_data_items['ctl_id'][idx].split(';')[0]
            ctl_raw_data = self.ctl_raw[ctl_id]
            ctl_raw_data_binned = self.ctl_raw_binned[ctl_id]
        else:
            ctl_raw_data = self.ctl_raw[idx, :]
            if isinstance(ctl_raw_data, sparse.csr_matrix):
                ctl_raw_data = ctl_raw_data.toarray().squeeze()
            ctl_raw_data_binned = self.ctl_raw_binned[idx, :]

        # get drug feat
        if self.args.dataset == 'transigen_sdst':
            pert_id = self._pert_id[idx]
            drug_feat = self.drug_feat[pert_id]
        else:
            pert_idx = self._pert_idx[idx]
            drug_feat = self.drug_feat[pert_idx]

        # transform to tensor
        drug_feat = tensor(drug_feat, dtype=torch.float32) if self.args.drug_feat != 'smi' else drug_feat
        trt_raw_data = tensor(trt_raw_data, dtype=torch.float32)
        ctl_raw_data = tensor(ctl_raw_data, dtype=torch.float32)
        trt_raw_data_binned = tensor(trt_raw_data_binned, dtype=torch.long)
        ctl_raw_data_binned = tensor(ctl_raw_data_binned, dtype=torch.long)

        pert_dose_idx = tensor(int(self._pert_dose_idx[idx]), dtype=torch.int64)
        pert_time_idx = tensor(int(self._pert_time_idx[idx]), dtype=torch.int64)

        pert_idx = tensor(int(self._pert_idx[idx]), dtype=torch.long)
        cell_idx = tensor(int(self._cell_idx[idx]), dtype=torch.long)
        tissue_idx = tensor(int(self._tissue_idx[idx]), dtype=torch.long)

        if self.args.mode != 'infer':
            combined_data = (
                trt_raw_data,
                ctl_raw_data,
                trt_raw_data_binned,
                ctl_raw_data_binned,
                drug_feat,
                pert_dose_idx,
                pert_time_idx,
                pert_idx,
                cell_idx,
                tissue_idx,
            )
        else:
            combined_data = (
                trt_raw_data,
                ctl_raw_data,
                trt_raw_data_binned,
                ctl_raw_data_binned,
                drug_feat,
                pert_dose_idx,
                pert_time_idx,
                pert_idx,
                cell_idx,
                tissue_idx,
                self._pert_id[idx],
                self._cell_id[idx],
                self._pert_dose[idx],
                self._pert_time[idx],
            )

        return combined_data
