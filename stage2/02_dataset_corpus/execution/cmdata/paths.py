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
# corpus.  Ticket #14 folded the older CM step into `02_dataset_corpus`; ticket
# #15 (ADR 0033) wrote the v2 corpus constants below and retired the view/tier
# layer (views are consumer-owned).
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

# Provenance of the vendored Stage 1 loaders (the data lake itself is not
# versioned by git). Kept as a constant so the record stays reproducible
# when the Stage 2 code is committed.
SOURCE_COMMIT = "c921f7d"

# Legacy identity-table helpers (`cmdata.corpus`): corpus/v1 held the gene set,
# the drug set, the exclusion report and the per-source gene sets, and no
# matrix.  Since ADR 0033 those tables live read-only in
# `stage2/_frozen/cm-data-v1/corpus/v1/`; ticket #16 re-exports their content
# under `output/corpus/v2/tables/`.  These paths are kept for the helper
# functions only; no build reads or writes them.
CORPUS_VERSION = "v1"
CORPUS_DIR = OUTPUT_DIR / "corpus" / CORPUS_VERSION
GENES_CSV = CORPUS_DIR / "genes.csv"
COMPOUNDS_CSV = CORPUS_DIR / "compounds.csv"
EXCLUSIONS_CSV = CORPUS_DIR / "excluded_compounds.csv"
CORPUS_README = CORPUS_DIR / "README.md"

# ---------------------------------------------------------------------------
# v2 corpus artifacts (ADR 0033; ticket #15 writes the four CM h5ads)
# ---------------------------------------------------------------------------
CORPUS_V2_DIR = OUTPUT_DIR / "corpus" / "v2"
TABLES_DIR = CORPUS_V2_DIR / "tables"           # ticket #16
TAHOE_H5AD = CORPUS_V2_DIR / "tahoe.h5ad"       # moved by ticket #14
BUILD_RECORD_DIR = EXECUTION_DIR / "build_record"
CORPUS_BUILD_RECORD = BUILD_RECORD_DIR / "corpus_build.json"
EXCLUDED_COMPOUNDS_RECORD = BUILD_RECORD_DIR / "excluded_compounds.csv"
REPRODUCIBILITY_RECORD = BUILD_RECORD_DIR / "corpus_reproducibility.json"

# Raw source id -> corpus dataset id (ADR 0033 decision 9, renames forward-only).
# `dtoxs_cm` splits back into its two arms; the GEO-sourced short names change.
SOURCE_TO_CORPUS = {
    "dtoxs_primary": "dtoxs_primary",
    "dtoxs_ipsc": "dtoxs_ipsc",
    "screen_ctox": "sctox",
    "temp_env": "tenv",
}
CORPUS_TO_SOURCE = {corpus: source for source, corpus in SOURCE_TO_CORPUS.items()}
CORPUS_IDS = ["dtoxs_primary", "dtoxs_ipsc", "sctox", "tenv"]


def corpus_h5ad(corpus_id: str) -> Path:
    """Path of one CM corpus h5ad under ``output/corpus/v2/``."""
    if corpus_id not in CORPUS_TO_SOURCE:
        raise KeyError(f"unknown corpus {corpus_id!r}")
    return CORPUS_V2_DIR / f"{corpus_id}.h5ad"


# The four CM corpora share one obs contract: the 24 ADR-0026 columns in their
# v1 order, then the three ticket-#15 additions (ADR 0033 decision 10).
CORPUS_OBS = [
    "well", "dataset_id", "compound_name", "drug_key", "inchikey", "smiles",
    "dose_um", "dose_source", "time_h", "is_control", "vehicle", "plate",
    "control_group", "cell_type", "donor", "replicate", "atc", "moa",
    "mesh_classes", "pharm_atc_moa", "pharm_mesh", "annotation_source",
    "lib_size", "n_genes_detected", "cell_line", "n_raw_rows", "is_blank",
]

# Row/axis/chemical acceptance counts (spec #12 contract table, ADR 0033 §2).
# A drift is a hard failure: the corpus is the information-maximal substrate
# and its shape is part of the published contract.
EXPECTED_CORPUS_COUNTS = {
    "dtoxs_primary": {
        "rows": 779, "perturbed": 677, "control": 102, "blank": 0,
        "chemicals": 47, "axis": 26328,
    },
    "dtoxs_ipsc": {
        "rows": 312, "perturbed": 296, "control": 16, "blank": 0,
        "chemicals": 37, "axis": 26328,
    },
    "sctox": {
        "rows": 3072, "perturbed": 2520, "control": 552, "blank": 0,
        "chemicals": 42, "axis": 31880,
    },
    "tenv": {
        "rows": 4510, "perturbed": 4206, "control": 272, "blank": 32,
        "chemicals": 460, "axis": 19613,
    },
}

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
