"""
Generalization-axis split masks.

**Cold-drug** (primary): 5-fold CV over the union of eligible drugs; each drug
is held out across all its doses and tested exactly once. Eligible = at least
`MIN_WELLS_PER_GROUP` treated wells in at least one dataset.

**Cold-dose** (secondary): leave-dose-level-out on `gse244740` only (the sole
multi-dose dataset), holding out one interior (interpolation) and one extreme
(extrapolation) dose level. No held-out dose level may appear in training.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (COLD_DOSE_EXTREME_LEVEL, COLD_DOSE_INTERIOR_LEVEL,  # noqa: E402
                   COLD_DRUG_FOLDS, COLD_DRUG_SEED, MIN_WELLS_PER_GROUP)


def eligible_drugs(groups, min_wells=MIN_WELLS_PER_GROUP):
    """Drugs with >= min_wells treated wells in at least one dataset."""
    per = groups.groupby(["dataset", "drug"])["n_wells"].sum().reset_index()
    ok = per[per["n_wells"] >= min_wells]
    return sorted(ok["drug"].unique())


def cold_drug_masks(groups, n_folds=COLD_DRUG_FOLDS, seed=COLD_DRUG_SEED,
                    min_wells=MIN_WELLS_PER_GROUP):
    """Assign every eligible drug to exactly one fold.

    Returns a DataFrame: drug, fold (int), eligible (bool).
    Ineligible drugs are recorded with fold = -1, eligible = False.
    """
    drugs = sorted(groups["drug"].unique())
    eligible = set(eligible_drugs(groups, min_wells))
    rng = np.random.default_rng(seed)
    order = list(drugs)
    rng.shuffle(order)
    elig_order = [d for d in order if d in eligible]
    fold_of = {}
    for i, d in enumerate(elig_order):
        fold_of[d] = i % n_folds
    rows = [{"drug": d, "fold": fold_of.get(d, -1), "eligible": d in eligible}
            for d in drugs]
    return pd.DataFrame(rows)


def cold_dose_masks(groups, interior=COLD_DOSE_INTERIOR_LEVEL,
                    extreme=COLD_DOSE_EXTREME_LEVEL):
    """Leave-dose-level-out on gse244740: two holdout configurations.

    Returns a DataFrame: dose_level (int), split ('interpolation'|
    'extrapolation'), held_out (bool).
    """
    g = groups[groups["dataset"] == "gse244740"]
    levels = sorted(int(x) for x in g["dose_level"].unique())
    if not levels:
        return pd.DataFrame(columns=["dose_level", "split", "held_out"])
    rows = []
    for level in levels:
        rows.append({"dose_level": level, "split": "interpolation",
                     "held_out": level == interior})
        rows.append({"dose_level": level, "split": "extrapolation",
                     "held_out": level == extreme})
    return pd.DataFrame(rows)
