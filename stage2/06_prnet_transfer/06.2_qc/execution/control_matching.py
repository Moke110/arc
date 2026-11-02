"""Control matching: the join keys ticket #3 consumes.

Rules (locked in ``rules.RULES[*]['control_rule']``):

* temp_env — a treated well is matched only to **DMSO vehicle** wells on the **same
  physical plate**. Media is a different arm (not a vehicle control) and Blank is a
  technical blank, so neither can enter the pool.
* Tahoe — a treated profile is matched only to DMSO_TF profiles of the **same plate
  and the same cell line**. A Tahoe plate pools ~50 cell lines, so a cross-line mean
  is the Simpson's-paradox mixture ADR 0029 rejects.

Both functions work on a per-observation frame that carries the QC verdict
(``qc_status``) and a role column, and they count what is actually resolvable
*after* QC — a plate whose vehicle wells all failed QC leaves its treated wells
unmatched.
"""

from __future__ import annotations

import pandas as pd

CONTROL_KEYS = {
    "temp_env": ("plate", "vehicle"),
    "tahoe": ("plate", "cell_line_id"),
}

CONTROL_ROLES = {
    "temp_env": "vehicle_control_dmso",
    "tahoe": "vehicle_control_dmso_tf",
}

# Roles that need a matched control. Media / Blank are control-arm or blank-arm rows:
# they are never "treated", so they are never reported as unmatched.
ARMS = {
    "temp_env": ("treated_compound", "positive_control_et1"),
    "tahoe": ("treated_compound",),
}

# Roles counted as the training / scoring treatment arm in the availability table.
TREATED_ROLES = {
    "temp_env": "treated_compound",
    "tahoe": "treated_compound",
}


def control_key(frame: pd.DataFrame, source: str) -> pd.Series:
    """Stable string key for the matched-control stratum of each row."""
    cols = CONTROL_KEYS[source]
    parts = [frame[c].astype(str) for c in cols]
    key = parts[0]
    for c, p in zip(cols[1:], parts[1:]):
        key = key + "::" + c + "=" + p
    return (source + "::" + key).rename("matched_control_key")


def match_controls(frame: pd.DataFrame, source: str) -> pd.DataFrame:
    """Add matched-control counts/keys to a per-observation frame.

    Returns a frame indexed like ``frame`` with columns
    ``matched_control_key``, ``matched_control_n``, ``matched_control_n_total``,
    ``matched_control_ok`` and ``matched_control_status``.

    ``matched_control_n`` counts QC-passing controls in the row's stratum;
    ``matched_control_n_total`` counts them before the QC filter, so the loss caused
    by QC is visible. Control rows themselves carry ``n = NaN`` and
    ``status = "is_control"``.
    """
    key = control_key(frame, source)
    is_control = frame["role"].eq(CONTROL_ROLES[source])
    is_arm = frame["role"].isin(ARMS[source])
    eligible = is_control & frame["qc_status"].eq("pass")

    n_elig = eligible.groupby(key).sum().astype("int64")
    n_all = is_control.groupby(key).sum().astype("int64")

    out = pd.DataFrame(index=frame.index)
    out["matched_control_key"] = key
    out["matched_control_n"] = key.map(n_elig).fillna(0).astype("int64")
    out["matched_control_n_total"] = key.map(n_all).fillna(0).astype("int64")

    ok = pd.Series(pd.NA, index=frame.index, dtype="boolean")
    ok[is_arm] = out.loc[is_arm, "matched_control_n"] >= 1
    out["matched_control_ok"] = ok

    status = pd.Series("not_an_arm", index=frame.index, dtype=object)
    status[is_control] = "is_control"
    status[is_arm & (out["matched_control_n"] >= 1)] = "matched"
    status[is_arm & (out["matched_control_n"] == 0)] = "unmatched"
    out["matched_control_status"] = status

    # control rows and non-arm rows (blank / media) carry no matching verdict
    out.loc[~is_arm, "matched_control_n"] = pd.NA
    out.loc[~is_arm, "matched_control_n_total"] = pd.NA
    return out


def availability_table(frame: pd.DataFrame, source: str) -> pd.DataFrame:
    """Per stratum: eligible controls, treated kept, and treated left unmatched."""
    generated = ["matched_control_key", "matched_control_n", "matched_control_n_total",
                 "matched_control_ok", "matched_control_status"]
    frame = frame.drop(columns=[c for c in generated if c in frame.columns])
    matched = match_controls(frame, source)
    f = frame.join(matched)
    cols = list(CONTROL_KEYS[source])
    ctrl_role = CONTROL_ROLES[source]

    treated_role = TREATED_ROLES[source]
    rows = []
    for key, g in f.groupby("matched_control_key", observed=True):
        is_ctrl = g["role"].eq(ctrl_role)
        is_treated = g["role"].eq(treated_role)
        kept_treated = is_treated & g["qc_status"].eq("pass")
        rows.append({
            **{c: g[c].iloc[0] for c in cols},
            "matched_control_key": key,
            "control_observations_total": int(is_ctrl.sum()),
            "control_observations_qc_pass": int((is_ctrl & g["qc_status"].eq("pass")).sum()),
            "treated_observations_total": int(is_treated.sum()),
            "treated_observations_qc_pass": int(kept_treated.sum()),
            "treated_kept_unmatched": int((kept_treated & (g["matched_control_n"] == 0)).sum()),
        })
    out = pd.DataFrame(rows).sort_values(list(CONTROL_KEYS[source])).reset_index(drop=True)
    return out
