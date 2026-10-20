"""Filesystem layout and frozen expectations for the Stage 2 CM corpus."""

from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent          # stage2/02_cm_data/cmdata
CM_DATA_DIR = HERE.parent                       # stage2/02_cm_data
REPO_ROOT = CM_DATA_DIR.parents[1]
PERT_DATA = REPO_ROOT / "pert_data"
ANNOTATION_DIR = PERT_DATA / "annotation"

HGNC_PATH = ANNOTATION_DIR / "hgnc_complete_set.txt"
COMPOUND_ANNOTATION_PATH = ANNOTATION_DIR / "cm_compound_annotation.csv"

CORPUS_VERSION = "v1"
# Provenance of the vendored Stage 1 loaders (the data lake itself is not
# versioned by git). Kept as a constant so the manifest stays reproducible
# when the Stage 2 code is committed.
SOURCE_COMMIT = "c921f7d"
CORPUS_DIR = CM_DATA_DIR / "corpus" / CORPUS_VERSION
GENES_CSV = CORPUS_DIR / "genes.csv"
COMPOUNDS_CSV = CORPUS_DIR / "compounds.csv"
EXCLUSIONS_CSV = CORPUS_DIR / "excluded_compounds.csv"
CORPUS_README = CORPUS_DIR / "README.md"

DEFAULT_BUNDLE = "datasets-1"
BUNDLE_DIR = CM_DATA_DIR / DEFAULT_BUNDLE

# Raw sources (read-only).
DTOXS_ROOT = PERT_DATA / "dtoxs_primary"
DTOXS_LEVEL3 = DTOXS_ROOT / "Data" / "Level3"
DTOXS_SMALL_MOL = DTOXS_ROOT / "Metadata" / "Small_Molecule_Metadata.txt"
DTOXS_ENTREZ_MAP = DTOXS_ROOT / "Metadata" / "UCSC-hg38-Gene-HGNC-Entrez-Map.tsv"
DTOXS_DCIC_XLSX = (DTOXS_ROOT / "Data"
                   / "PromoCell-DGE-RNAseq-Experiment-Design.DCIC.xlsx")
DTOXS_IPSC = PERT_DATA / "dtoxs_ipsc"
GSE244740 = PERT_DATA / "gse244740"
GSE262419 = PERT_DATA / "gse262419"

SOURCE_DATASETS = ["dtoxs_primary", "dtoxs_ipsc", "gse244740", "gse262419"]

# Stage 1 regression baselines (read-only; outputs, never code).
STAGE1_OUTPUT = (REPO_ROOT / "stage1" / "xpert-fda-miri-v3"
                 / "03_cm_finetune_data" / "output")

# Frozen per-source expectations (ticket 03 / 08).  A loader drift is a hard
# failure in `sources.check_sources`.
EXPECTED_SOURCE_COUNTS = {
    "dtoxs_primary": {
        "wells": 842, "controls": 102, "perturbed": 740,
        "small_molecule_perturbed": 677, "distinct_compound_names": 54,
    },
    "dtoxs_ipsc": {
        "wells": 344, "controls": 16, "perturbed": 328,
        "small_molecule_perturbed": 296, "distinct_compound_names": 42,
    },
    "gse244740": {
        "wells": 3072, "controls": 552, "perturbed": 2520,
        "small_molecule_perturbed": 2520, "distinct_compound_names": 44,
    },
    "gse262419": {
        "wells": 4558, "controls": 272, "perturbed": 4254,
        "small_molecule_perturbed": 4206, "distinct_compound_names": 466,
        "blanks": 32,
    },
}

# Frozen view expectations (README / tickets 04-08).  `manifest.csv` remains
# authoritative; these are QA expectations and a guard against silent drift.
EXPECTED_VIEWS = {
    "dtoxs_cm.h5ad": {
        "perturbed": 957, "controls": 118, "drugs": 46,
    },
    "gse244740.h5ad": {
        "perturbed": 2520, "controls": 552, "drugs": 42,
    },
    "gse262419_atc-moa.h5ad": {
        "perturbed": 1047, "controls": 272, "drugs": 113,
    },
    "gse262419_atc-moa-mesh.h5ad": {
        "perturbed": 2310, "controls": 272, "drugs": 252,
    },
    "gse262419_all-smiles.h5ad": {
        "perturbed": 4206, "controls": 272, "drugs": 460,
    },
}

# Per-source pharmacology tier counts (drug keys).  These are the ticket-02
# audited counts and gate a rerun against a silent identity/pharmacology drift.
EXPECTED_TIER_COUNTS = {
    "dtoxs": {"atc-moa": 46, "atc-moa-mesh": 46, "all-smiles": 47},
    "gse244740": {"atc-moa": 42, "atc-moa-mesh": 42, "all-smiles": 42},
    "gse262419": {"atc-moa": 113, "atc-moa-mesh": 252, "all-smiles": 460},
}
TIER_DATASETS = {
    "dtoxs": ["dtoxs_primary", "dtoxs_ipsc"],
    "gse244740": ["gse244740"],
    "gse262419": ["gse262419"],
}

VIEWS = [
    ("dtoxs_cm.h5ad", ["dtoxs_primary", "dtoxs_ipsc"], "atc-moa"),
    ("gse244740.h5ad", ["gse244740"], "atc-moa"),
    ("gse262419_atc-moa.h5ad", ["gse262419"], "atc-moa"),
    ("gse262419_atc-moa-mesh.h5ad", ["gse262419"], "atc-moa-mesh"),
    ("gse262419_all-smiles.h5ad", ["gse262419"], "all-smiles"),
]

# The three gse262419 views share one corpus projection; QA asserts it.
GSE262419_VIEWS = [name for name, _, _ in VIEWS if name.startswith("gse262419_")]
