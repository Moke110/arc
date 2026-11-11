"""Preserve continuous metrics at group, dose and run granularity."""
from functools import lru_cache
import numpy as np
import pandas as pd
from scipy.stats import rankdata

FIELDS = ('PCC', 'Spearman', 'RMSE', 'R2', 'PCC_native', 'Spearman_native', 'RMSE_native', 'R2_native')


@lru_cache(maxsize=8)
def group_truth(groups_path, matrix_path, digests):
    groups = pd.read_csv(groups_path)
    matrix = np.load(matrix_path, mmap_mode='r')
    if matrix.shape != (len(groups), 18164) or groups.group_id.duplicated().any():
        raise ValueError('invalid continuous truth axis or group identities')
    return {group: (matrix, i) for i, group in enumerate(groups.group_id)}


@lru_cache(maxsize=4)
def cached_groups(path, digest):
    frame = pd.read_csv(path, float_precision='round_trip').set_index(['role', 'model', 'group_id'])
    if not frame.index.is_unique:
        raise ValueError('duplicate preserved continuous group identity')
    return frame


def parent_identity(item):
    keys = ('prediction', 'cases', 'genes', 'train_cases', 'train_delta', 'train_controls')
    return {key: item[key]['sha256'] for key in keys if key in item}


class ContinuousMetrics:
    def __init__(self, item, resolve):
        self.item = item
        self.cache = None
        self.truth = {}
        if 'continuous_cache' in item:
            if item['continuous_cache_parent'] != parent_identity(item):
                raise ValueError('preserved continuous metrics differ from the frozen prediction parent')
            self.cache = cached_groups(str(resolve(item['continuous_cache'])), item['continuous_cache']['sha256'])
        elif 'continuous_truth' in item:
            axis = pd.read_csv(resolve(item['continuous_axis'])).set_index('hgnc_id').index
            genes = pd.read_csv(resolve(item['genes']))
            self.columns = axis.get_indexer(genes.hgnc_id)
            if np.any(self.columns < 0):
                raise ValueError('continuous prediction axis is outside measured CM truth')
            for release in item['continuous_truth']:
                values = group_truth(str(resolve(release['groups'])), str(resolve(release['delta'])),
                                     (release['groups']['sha256'], release['delta']['sha256']))
                if self.truth.keys() & values.keys():
                    raise ValueError('continuous truth group appears in multiple folds')
                self.truth.update(values)

    def row(self, group, prediction):
        if self.cache is not None:
            row = self.cache.loc[(self.item['role'], self.item['model'], group)]
            return {name: float(row[name]) if pd.notna(row[name]) else None for name in FIELDS}
        if not self.truth:
            return {}
        matrix, index = self.truth[group]
        truth = matrix[index, self.columns].astype(np.float64)
        prediction = np.asarray(prediction, dtype=np.float64)
        if truth.shape != prediction.shape or not np.isfinite(truth).all():
            raise ValueError('continuous prediction/truth axis differs')
        a, b = truth - truth.mean(), prediction - prediction.mean()
        variance, other = float(a @ a), float(b @ b)
        result = {'PCC': None, 'Spearman': None}
        if variance and other:
            result['PCC'] = float((a @ b) / np.sqrt(variance * other))
            result['Spearman'] = float(np.corrcoef(rankdata(truth), rankdata(prediction))[0, 1])
        error = float((truth - prediction) @ (truth - prediction))
        result['RMSE'] = float(np.sqrt(error / len(truth)))
        result['R2'] = 1 - error / variance if variance else None
        return result
