# FROZEN COPY — archived evidence, NOT executable code and NOT a step entry point.
# Origin: stage2/02_cm_data/execution/cmdata/paths.py (view/tier constants) @ commit 7e7d689 (2026-11-04).
# See ./README.md. Do not wire this into any build; do not edit it.
"""Frozen v1 view/tier constants (archival excerpt).

Extracted verbatim from the retired ``stage2/02_cm_data/execution/cmdata/paths.py``.
Only the constants that describe the five v1 views and their pharmacology tiers
are kept; the step's layout, raw-source roots and Stage 1 baseline belong to the
retired producer and are deliberately not frozen here.
"""

from __future__ import annotations

DEFAULT_BUNDLE = "datasets-1"
CORPUS_VERSION = "v1"

# Frozen view expectations (README / tickets 04-08).  `manifest.csv` remains
# authoritative; these are QA expectations and a guard against silent drift.
EXPECTED_VIEWS = {
    "dtoxs_cm.h5ad": {
        "perturbed": 957, "controls": 118, "drugs": 46,
    },
    "screen_ctox.h5ad": {
        "perturbed": 2520, "controls": 552, "drugs": 42,
    },
    "temp_env_atc-moa.h5ad": {
        "perturbed": 1047, "controls": 272, "drugs": 113,
    },
    "temp_env_atc-moa-mesh.h5ad": {
        "perturbed": 2310, "controls": 272, "drugs": 252,
    },
    "temp_env_all-smiles.h5ad": {
        "perturbed": 4206, "controls": 272, "drugs": 460,
    },
}

# Per-source pharmacology tier counts (drug keys).  These are the ticket-02
# audited counts and gate a rerun against a silent identity/pharmacology drift.
EXPECTED_TIER_COUNTS = {
    "dtoxs": {"atc-moa": 46, "atc-moa-mesh": 46, "all-smiles": 47},
    "screen_ctox": {"atc-moa": 42, "atc-moa-mesh": 42, "all-smiles": 42},
    "temp_env": {"atc-moa": 113, "atc-moa-mesh": 252, "all-smiles": 460},
}
TIER_DATASETS = {
    "dtoxs": ["dtoxs_primary", "dtoxs_ipsc"],
    "screen_ctox": ["screen_ctox"],
    "temp_env": ["temp_env"],
}

VIEWS = [
    ("dtoxs_cm.h5ad", ["dtoxs_primary", "dtoxs_ipsc"], "atc-moa"),
    ("screen_ctox.h5ad", ["screen_ctox"], "atc-moa"),
    ("temp_env_atc-moa.h5ad", ["temp_env"], "atc-moa"),
    ("temp_env_atc-moa-mesh.h5ad", ["temp_env"], "atc-moa-mesh"),
    ("temp_env_all-smiles.h5ad", ["temp_env"], "all-smiles"),
]

# The three temp_env views share one corpus projection; QA asserts it.
TEMP_ENV_VIEWS = [name for name, _, _ in VIEWS if name.startswith("temp_env_")]
