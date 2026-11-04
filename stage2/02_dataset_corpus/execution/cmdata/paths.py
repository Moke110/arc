"""Filesystem layout and frozen expectations for the Stage 2 CM corpus."""

from __future__ import annotations

from pathlib import Path


def find_repo_root(start: Path) -> Path:
    """Walk up from *start* to the directory that contains `.git`.

    Paths must not encode depth (ADR 0030 §4): the package now lives one level
    deeper, in `execution/cmdata/`, and a `parents[N]` count would silently
    resolve above the repository root.
    """
    for candidate in (start, *start.parents):
        if (candidate / ".git").exists():
            return candidate
    raise RuntimeError(f"cannot find the repository root above {start}")


HERE = Path(__file__).resolve().parent          # stage2/02_dataset_corpus/execution/cmdata
REPO_ROOT = find_repo_root(HERE)
# The step's two areas, derived from the repo root rather than by counting
# ancestors: `execution/` holds this package and the build, `output/` holds the
# corpus.  Ticket #14 folded the older CM step into `02_dataset_corpus`; the
# corpus and view/tier constants below still describe the retired v1 layout and
# are rewritten by ticket #15 (ADR 0033).
STEP_DIR = REPO_ROOT / "stage2" / "02_dataset_corpus"
EXECUTION_DIR = STEP_DIR / "execution"
OUTPUT_DIR = STEP_DIR / "output"
assert (EXECUTION_DIR / "cmdata" / "paths.py").resolve() == Path(__file__).resolve(), (
    f"02_dataset_corpus layout mismatch: {EXECUTION_DIR} does not hold this package"
)
PERT_DATA = REPO_ROOT / "pert_data"
ANNOTATION_DIR = PERT_DATA / "Annotation" / "source"   # the HGNC download
ANNOTATION_LOCAL = PERT_DATA / "Annotation" / "local"   # local additions

HGNC_PATH = ANNOTATION_DIR / "hgnc_complete_set.txt"
COMPOUND_ANNOTATION_PATH = ANNOTATION_LOCAL / "cm_compound_annotation.csv"

CORPUS_VERSION = "v1"
# Provenance of the vendored Stage 1 loaders (the data lake itself is not
# versioned by git). Kept as a constant so the manifest stays reproducible
# when the Stage 2 code is committed.
SOURCE_COMMIT = "c921f7d"
CORPUS_DIR = OUTPUT_DIR / "corpus" / CORPUS_VERSION
GENES_CSV = CORPUS_DIR / "genes.csv"
COMPOUNDS_CSV = CORPUS_DIR / "compounds.csv"
EXCLUSIONS_CSV = CORPUS_DIR / "excluded_compounds.csv"
CORPUS_README = CORPUS_DIR / "README.md"

DEFAULT_BUNDLE = "datasets-1"
BUNDLE_DIR = OUTPUT_DIR / DEFAULT_BUNDLE

# Raw sources (read-only).
DTOXS_ROOT = PERT_DATA / "LINCS-LDS-1593" / "source"
DTOXS_LEVEL3 = DTOXS_ROOT / "Data" / "Level3"
DTOXS_SMALL_MOL = DTOXS_ROOT / "Metadata" / "Small_Molecule_Metadata.txt"
DTOXS_ENTREZ_MAP = DTOXS_ROOT / "Metadata" / "UCSC-hg38-Gene-HGNC-Entrez-Map.tsv"
DTOXS_DCIC_XLSX = (DTOXS_ROOT / "Data"
                   / "PromoCell-DGE-RNAseq-Experiment-Design.DCIC.xlsx")
DTOXS_IPSC = PERT_DATA / "DToxS-iPSC" / "source"

# Raw GSE folders keep the GEO accession on disk; the Stage 2 dataset ids are
# the short names used everywhere the build writes data.
SCREEN_CTOX_DIR = PERT_DATA / "GSE244740" / "source"       # GSE244740
TEMP_ENV_DIR = PERT_DATA / "GSE262419" / "source"          # GSE262419

SOURCE_DATASETS = ["dtoxs_primary", "dtoxs_ipsc", "screen_ctox", "temp_env"]

# Stage 2 dataset id -> GEO accession.  This is the single place where the
# memorable name and the accession meet; the raw folder, the read-only Stage 1
# baseline and the papers keep the accession.  Documented in README.md,
# corpus/v1/README.md and ADR 0026.
SOURCE_ACCESSIONS = {
    "screen_ctox": "GSE244740",   # ScreenSeq hiPSC-CM, drug cardiotoxicity panel
    "temp_env": "GSE262419",      # TempO-Seq hiPSC-CM, environmental chemical panel
}

RAW_DIRS = {
    "dtoxs_primary": DTOXS_ROOT,
    "dtoxs_ipsc": DTOXS_IPSC,
    "screen_ctox": SCREEN_CTOX_DIR,
    "temp_env": TEMP_ENV_DIR,
}


def source_accession(dataset_id: str) -> str:
    """GEO accession for a Stage 2 dataset id, or "" for a DToxS panel."""
    return SOURCE_ACCESSIONS.get(dataset_id, "")


def stage1_dataset_name(dataset_id: str) -> str:
    """Name the read-only Stage 1 baseline uses for this Stage 2 dataset."""
    accession = SOURCE_ACCESSIONS.get(dataset_id)
    return accession.lower() if accession else dataset_id

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
    "screen_ctox": {
        "wells": 3072, "controls": 552, "perturbed": 2520,
        "small_molecule_perturbed": 2520, "distinct_compound_names": 44,
    },
    "temp_env": {
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
