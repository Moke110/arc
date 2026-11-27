"""Chemical representation inputs aligned to model conditions."""

import numpy as np

from cpst_util.representations import WIDTHS


def aligned_features(vectors, arm, parents):
    """Bind the source vectors to each pooled condition, preserving duplicates."""
    parents_unique, inverse = np.unique(parents, return_inverse=True)
    features = vectors.vectors(arm, parents_unique)
    if features.shape != (len(parents_unique), WIDTHS[arm]) or features.dtype != np.float32:
        raise ValueError(f"unexpected frozen representation width or dtype: {arm}")
    return features[inverse]
