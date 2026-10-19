"""
Per-plate control baseline and per-well Δ.

The baseline for a treated well is the mean of the control wells on the **same
plate** (Stage 1's ADR 0016 method, referenced only). Δ = normalized treated
well − its per-plate baseline. A plate with no controls falls back to the
dataset-wide control mean and is recorded (never silently dropped).
"""

import numpy as np


def per_plate_baselines(meta, norm, is_control):
    """plate -> baseline vector (mean of that plate's control wells).

    Returns (baselines dict, plates_without_controls list).
    """
    baselines, no_ctrl = {}, []
    is_control = np.asarray(is_control, dtype=bool)
    plate_col = meta["plate"].astype(str)
    for plate in sorted(plate_col.unique()):
        ctl = (plate_col == plate).to_numpy() & is_control
        if ctl.any():
            baselines[plate] = norm[ctl].mean(axis=0)
        else:
            no_ctrl.append(plate)
    if no_ctrl:  # dataset-wide fallback for control-less plates
        fallback = norm[is_control].mean(axis=0)
        for plate in no_ctrl:
            baselines[plate] = fallback
    return baselines, no_ctrl


def per_well_delta(meta, norm, baselines, is_control):
    """Δ for each treated well (controls get their baseline subtracted -> 0).

    Returns (treated_meta, delta, baseline_rows): the treated subset of `meta`,
    its Δ matrix, and the per-well baseline transcriptomes — the model input.
    """
    plate = meta["plate"].astype(str).values
    baseline_rows = np.vstack([baselines[p] for p in plate])
    delta = norm - baseline_rows
    treated = ~np.asarray(is_control, dtype=bool)
    return (meta.loc[treated].reset_index(drop=True), delta[treated],
            baseline_rows[treated])
