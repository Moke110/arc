from __future__ import annotations

from collections import OrderedDict
import random
from typing import Any

import h5py
import numpy as np
import pandas as pd
import sklearn
import torch
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold


def setup_seed(random_seed: int) -> None:
    np.random.seed(random_seed)
    random.seed(random_seed)
    torch.manual_seed(random_seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(random_seed)
        torch.cuda.manual_seed_all(random_seed)
        torch.backends.cudnn.enabled = False
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True


def seed_worker(worker_id: int) -> None:
    worker_seed = torch.initial_seed() % (2**32)
    np.random.seed(worker_seed)
    random.seed(worker_seed)


def load_from_HDF(fname: str) -> dict[str, np.ndarray]:
    def _to_str_array(value: np.ndarray) -> np.ndarray:
        if value.size == 0:
            return value
        if np.issubdtype(value.dtype, np.bytes_):
            return value.astype(str)
        if value.dtype.kind == "O":
            flat = value.reshape(-1)
            if all(isinstance(x, (bytes, np.bytes_)) for x in flat):
                return np.array([x.decode("utf-8", errors="replace") for x in flat], dtype=str).reshape(value.shape)
            if all(isinstance(x, str) for x in flat):
                return np.array(flat, dtype=str).reshape(value.shape)
        return value

    data: dict[str, np.ndarray] = {}
    with h5py.File(fname, "r") as f:
        for key in f:
            value = np.asarray(f[key])
            value = _to_str_array(value)
            data[key] = value
    return data


def save_to_HDF(fname: str, data: dict[str, np.ndarray]) -> None:
    def _to_bytes_array(value: np.ndarray) -> np.ndarray:
        if value.size == 0:
            return value
        if np.issubdtype(value.dtype, np.str_):
            return value.astype(np.bytes_)
        if value.dtype.kind == "O":
            flat = value.reshape(-1)
            if all(isinstance(x, str) for x in flat):
                return np.array(flat, dtype=str).reshape(value.shape).astype(np.bytes_)
        return value

    with h5py.File(fname, "w") as f:
        for key, item in data.items():
            value = np.asarray(item)
            value = _to_bytes_array(value)
            f[key] = value


def subsetDict(d: dict[str, np.ndarray], ind: Any) -> dict[str, np.ndarray]:
    res: dict[str, np.ndarray] = {}
    if not isinstance(ind, np.ndarray):
        ind = np.asarray(ind)
    for k, v in d.items():
        if isinstance(v, np.ndarray):
            res[k] = v[ind]
    return res


def concat_dicts(a: OrderedDict, b: OrderedDict) -> OrderedDict:
    assert len(set(a.keys()) & set(b.keys())) == 0, "Can't unambiguously concatenate."
    return type(a)(list(a.items()) + list(b.items()))


def getSplitsByGroupKFold(
    groups: np.ndarray, n_splits: int, shuffle: bool, random_state: int
) -> tuple[list[int], list[int], list[int]]:
    assert n_splits >= 3
    kf = GroupKFold(n_splits=n_splits)

    if shuffle:
        unique_groups = np.unique(groups)
        rnd_renames = sklearn.utils.shuffle(np.arange(len(unique_groups)), random_state=random_state)
        groups_renamed = np.array([rnd_renames[np.argwhere(unique_groups == g)[0]] for g in groups])
        kfsplit = kf.split(X=np.zeros(groups.shape[0]), groups=groups_renamed)
    else:
        kfsplit = kf.split(X=np.zeros(groups.shape[0]), groups=groups)

    folds = [list(x[1]) for x in kfsplit]
    folds_nums = list(range(len(folds)))
    ind_tr = sum([folds[i] for i in folds_nums[:-2]], [])
    ind_va = folds[folds_nums[-2]]
    ind_te = folds[folds_nums[-1]]
    return ind_tr, ind_va, ind_te


def split_data(
    data: dict[str, np.ndarray], n_folds: int = 5, split_type: str = "random_split", rnds: int | None = None
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray]]:
    if split_type == "random_split":
        ind_tr, ind_va, ind_te = getSplitsByGroupKFold(data["sig"], n_folds, shuffle=True, random_state=rnds)
    elif split_type == "smiles_split":
        ind_tr, ind_va, ind_te = getSplitsByGroupKFold(
            data["canonical_smiles"], n_folds, shuffle=True, random_state=rnds
        )
    else:
        raise ValueError(f"Unsupported split_type: {split_type}")
    return subsetDict(data, ind_tr), subsetDict(data, ind_va), subsetDict(data, ind_te)


def rmse(targets: np.ndarray, preds: np.ndarray) -> float:
    return float(np.sqrt(mean_squared_error(targets, preds)))


def pearson(targets: np.ndarray, preds: np.ndarray) -> float:
    """Pearson r；当任一侧沿基因维方差为 0 或输入非有限值时 scipy 会给出 nan，此处给出可聚合的约定值。"""
    t = np.asarray(targets, dtype=np.float64).ravel()
    p = np.asarray(preds, dtype=np.float64).ravel()
    if t.size < 2 or p.size < 2 or t.size != p.size:
        return float("nan")
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(p)):
        return float("nan")
    st, sp = float(np.std(t)), float(np.std(p))
    if st == 0.0 or sp == 0.0:
        # 未定义：全常数谱与重建在向量相等时视为完全相关，否则视为无线性相关
        return 1.0 if np.allclose(t, p, rtol=0.0, atol=1e-6) else 0.0
    try:
        r = float(pearsonr(t, p)[0])
        return r if np.isfinite(r) else 0.0
    except ValueError:
        return float("nan")


def spearman(targets: np.ndarray, preds: np.ndarray) -> float:
    try:
        return float(spearmanr(targets, preds)[0])
    except ValueError:
        return float("nan")


def precision_100(label_test: np.ndarray, label_predict: np.ndarray) -> tuple[float, float]:
    """Negative / Positive Precision@100 (论文指标): 预测排序中最低/最高各 100 个基因与真值集合的重叠率。

    Returns:
        (negative_precision, positive_precision)，均在 [0, 1]。
    """
    k = 100
    num_pos = 100
    num_neg = 100
    label_test = np.argsort(label_test)
    label_predict = np.argsort(label_predict)
    neg_test = set(label_test[:num_neg])
    pos_test = set(label_test[-num_pos:])
    neg_predict = set(label_predict[:k])
    pos_predict = set(label_predict[-k:])
    return len(neg_test.intersection(neg_predict)) / k, len(pos_test.intersection(pos_predict)) / k


def get_metric_func(metric: str):
    if metric == "rmse":
        return rmse
    if metric == "mae":
        return mean_absolute_error
    if metric == "r2":
        return r2_score
    if metric == "pearson":
        return pearson
    if metric == "spearman":
        return spearman
    if metric == "precision100":
        return precision_100
    raise ValueError(f'Metric "{metric}" not supported.')
