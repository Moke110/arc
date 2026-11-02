"""Frozen-input contract for Stage 2 step 06.3 (`feat/pt/t03-input`, issue #3).

This module is the single source of truth for every decision ticket #3 freezes:
the gene axis, the two-source artifact layout, the five drug-key folds, the
per-fold consensus control, the locked training config and the provenance rules.
`build_input.py` only orchestrates; `folds.py`, `consensus.py` and `writers.py`
consume the constants declared here.

Nothing in this module writes to a frozen artifact.  Frozen sources are opened
read-only and are addressed through `SourceRecord`s so the manifest can record
exactly which bytes were read (ADR 0030 §4: resolve the repo root by walking up
to `.git`; never count parents).
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------------- #
# Versioning / identities
# --------------------------------------------------------------------------- #

TOOL_VERSION = "1.0.0"
MANIFEST_SCHEMA = "06.3_input/manifest/1"
CONFIG_VERSION = "1.0.0"
CONTRACT_ID = "stage2/06_prnet_transfer/06.3_input/contract@1.0.0"

# --------------------------------------------------------------------------- #
# Fixed numeric contract (ADR 0032 + spec §统一表达输入契约/§切分与评价)
# --------------------------------------------------------------------------- #

N_GENES = 18164                    # frozen CM corpus axis, order preserved
N_FOLDS = 5                        # five drug-key folds
FOLD_SEED = 20261102               # derived_default: study start date, see config
MODEL_SEEDS = (20261021, 20261022, 20261023)   # derived_default (repo seed family)
TAHOE_VALID_PLATES = ("plate13", "plate14")     # existing 05 technical default

FAM153B_SYMBOL = "FAM153B"
FAM153B_ENSEMBL = "ENSG00000182230"

# Tahoe has no time column anywhere in the atlas (ADR 0029 update 2026-11-01,
# ADR 0031 D3).  The 24 h claim is literature provenance, never a metadata value.
TIME_SOURCE_TAHOE = "literature:10.1016/j.cell.2026.08.035"
TIME_SOURCE_NOTE = (
    "Atlas design states a single 24 h treatment. The quoted authors' sentence "
    "is in the v2 preprint (doi 10.1101/2025.02.20.639398), cited by ticket #2's "
    "locked qc_rules.json time_sources.tahoe_preprint_24h; the formal record is "
    "the Cell 2026 paper (doi 10.1016/j.cell.2026.08.035), the citation ADR 0029 "
    "already sanctions and the value written in obs.time_source. No atlas table "
    "carries a time column, so obs.time_h stays NaN and the value is provenance, "
    "not raw metadata, and must not be used to filter rows."
)

CM_REPLICATE_PROVENANCE = (
    "obs.replicate is a GENERATED within-group index (cmdata/views.py:50 "
    "_replicate_index, ordered by well over control_group|drug_key|dose_um); it "
    "agrees with the raw Replicate_Number on 99.665% of temp_env rows and is not "
    "biological.  The raw design value is obs.replicate_raw, re-joined from "
    "pert_data/GSE262419/source/GSE262419_hash.csv.gz by (Plate_ID, normalised "
    "Well_ID)."
)
CM_TIME_PROVENANCE = (
    "obs.time_h=24.0 on every temp_env row is a loader constant "
    "(cmdata/sources.py:539), not a field read from GEO. The 24 h design is "
    "sourced from the GSE262419 study paper: Tsai et al., Chem Res Toxicol 2024, "
    "doi:10.1021/acs.chemrestox.4c00193 (PMID 39046974), quoted in ticket #2's "
    "locked qc_rules.json time_sources.tempenv_paper_24h."
)

# Column name sets -----------------------------------------------------------

# The 24 inherited CM view columns, exact order (ADR 0026 / 02 README).
CM_OBS_INHERITED = [
    "well", "dataset_id", "compound_name", "drug_key", "inchikey", "smiles",
    "dose_um", "dose_source", "time_h", "is_control", "vehicle", "plate",
    "control_group", "cell_type", "donor", "replicate", "atc", "moa",
    "mesh_classes", "pharm_atc_moa", "pharm_mesh", "annotation_source",
    "lib_size", "n_genes_detected",
]

# The 35 inherited Tahoe view columns, exact order (ADR 0031 / 04 README).
TAHOE_OBS_INHERITED = CM_OBS_INHERITED + [
    "cell_line_id", "cell_name", "n_cells", "moa_broad", "moa_fine",
    "pubchem_cid", "has_structure", "has_matched_control",
    "SMILES", "dose", "paired_control_index",
]

CM_OBS_ADDED = [
    "observation_id", "physical_well", "plate_no", "replicate_raw", "solvent",
    "chemical_index", "eligibility_decision", "elig_class", "qc_pass",
    "qc_reason", "matched_control_ok", "matched_control_key", "main_set",
    "expansion_set",
]

TAHOE_OBS_ADDED = [
    "observation_id", "physical_well", "tahoe_split", "tahoe_valid_group",
    "qc_pass", "qc_reason", "matched_control_ok", "matched_control_key",
    "tahoe_usable", "time_source",
]

CM_OBS = CM_OBS_INHERITED + CM_OBS_ADDED
TAHOE_OBS = TAHOE_OBS_INHERITED + TAHOE_OBS_ADDED

# var carries the corpus identity columns plus the per-source measurement mask.
VAR_COLUMNS = ["hgnc_id", "entrez_id", "ensembl_id", "match_type",
               "gene_measured", "n_nonzero"]

ELIGIBILITY_DECISIONS = ("include", "exclude", "uncertain")
FOLD_ROLES = ("train", "valid", "test")
QC_SOURCES = ("cm", "tahoe")

# Artifact file names (the interface ticket #4 reads).
ARTIFACTS = {
    "cm": "cm_transfer.h5ad",
    "tahoe": "tahoe_transfer.h5ad",
    "folds": "folds.parquet",
    "consensus_controls": "consensus_controls.npz",
    "consensus_membership": "consensus_membership.parquet",
    "tahoe_control_groups": "tahoe_control_groups.parquet",
    "tahoe_pairs": "tahoe_pairs.npz",
    "locked_config": "locked_config.json",
    "manifest": "manifest.json",
    "counts": "counts.csv",
    "fold_sizes": "fold_sizes.csv",
    "consensus_summary": "consensus_summary.csv",
    "var_checks": "var_checks.csv",
    "contract": "contract.json",
}

# Tracked record files: small, force-added past the repository-wide `output/`
# and `*.csv` ignore rules (ADR 0030 / 04-05 precedent).  The two h5ad files are
# local-only (16 GB / 218 MB) and never committed.  Both .npz caches are part of
# the interface ticket #4 reads (c* and the treated-row -> control-group map), so
# they are tracked as well — they are kilobytes, and a frozen input whose
# required artifacts are absent from git is not a frozen input.
TRACKED_RECORD = [
    "manifest.json", "locked_config.json", "contract.json", "counts.csv",
    "fold_sizes.csv", "consensus_summary.csv", "var_checks.csv",
    "folds.parquet", "consensus_membership.parquet", "consensus_controls.npz",
    "tahoe_control_groups.parquet", "tahoe_pairs.npz",
]


# --------------------------------------------------------------------------- #
# Repository root (ADR 0030 §4: walk up to .git, never count parents)
# --------------------------------------------------------------------------- #

class ContractError(RuntimeError):
    """A frozen-input invariant is violated; the build must fail closed."""


def find_repo_root(start: Path | str | None = None) -> Path:
    path = Path(start).resolve() if start is not None else Path(__file__).resolve()
    if path.is_file():
        path = path.parent
    for candidate in (path, *path.parents):
        if (candidate / ".git").exists():
            return candidate
    raise ContractError(f"no .git directory at or above {path}")


ROOT = find_repo_root()
STEP = ROOT / "stage2" / "06_prnet_transfer"
SUBSTEP = STEP / "06.3_input"
EXECUTION = SUBSTEP / "execution"

# Data root.  The implementation lives in a git worktree that carries only
# tracked files, while the frozen 16 GB sources (and the `pert_data/` source
# trees) live in the main checkout.  `ARC_DATA_ROOT` says where to *read* data;
# it defaults to the code root, which is correct in the main checkout and after
# merge.  All recorded paths are relative to this root, so the merged repository
# reads the same strings.
DATA_ROOT = Path(os.environ.get("ARC_DATA_ROOT", str(ROOT))).expanduser().resolve()
STEP_DATA = DATA_ROOT / "stage2" / "06_prnet_transfer"

# Frozen inputs (read-only; never written by this ticket).
CORPUS_GENES = DATA_ROOT / "stage2/02_cm_data/output/corpus/v1/genes.csv"
CORPUS_SOURCE_GENE_SETS = DATA_ROOT / "stage2/02_cm_data/output/corpus/v1/source_gene_sets.csv"
CM_VIEW = DATA_ROOT / "stage2/02_cm_data/output/datasets-1/temp_env_all-smiles.h5ad"
CM_BUNDLE_MANIFEST = DATA_ROOT / "stage2/02_cm_data/output/datasets-1/manifest.csv"
TAHOE_VIEW = DATA_ROOT / "stage2/04_tahoe_data/output/pretrain_view/tahoe_cm_schema.h5ad"
TAHOE_VIEW_MANIFEST = DATA_ROOT / "stage2/04_tahoe_data/output/pretrain_view/manifest.csv"
TAHOE_CORPUS = DATA_ROOT / "stage2/04_tahoe_data/output/pseudobulk/tahoe_pseudobulk.h5ad"
GSE262419_HASH = DATA_ROOT / "pert_data/GSE262419/source/GSE262419_hash.csv.gz"


def relpath(path: Path | str) -> str:
    """Path string relative to the data root (the merged repository layout)."""
    resolved = Path(path).resolve()
    for base in (DATA_ROOT, ROOT):
        try:
            return str(resolved.relative_to(base))
        except ValueError:
            continue
    return str(resolved)

# Upstream handoffs (#1, #2).  Candidates are searched in order; every path can
# be overridden on the CLI.  See README "Upstream handoffs".
HANDOFF_CANDIDATES = {
    "eligibility": [
        # #1's designated durable artifact (one row per structure identity)
        "06.1_eligibility/output/eligibility_table.csv",
        "06.1_eligibility/output/eligibility.csv",
        "06.1_eligibility/output/eligibility.parquet",
        "06.1_eligibility/output/drug_eligibility.csv",
        "06.1_eligibility/output/eligibility_handoff.csv",
    ],
    "cm_qc": [
        "06.2_qc/output/cm_obs_qc.csv",
        "06.2_qc/output/cm_qc.csv",
        "06.2_qc/output/qc_cm.csv",
        "06.2_qc/output/temp_env_qc.csv",
        "06.2_qc/output/cm_observations.csv",
    ],
    "tahoe_qc": [
        "06.2_qc/output/tahoe_obs_qc.csv",
        "06.2_qc/output/tahoe_qc.csv",
        "06.2_qc/output/qc_tahoe.csv",
        "06.2_qc/output/tahoe_observations.csv",
    ],
    # A single combined per-observation table with a `source` column is accepted.
    "qc_combined": [
        # #2's designated durable handoff (71,576 rows: temp_env + tahoe)
        "06.2_qc/output/observations.csv.gz",
        "06.2_qc/output/obs_qc.csv",
        "06.2_qc/output/observation_qc.csv",
        "06.2_qc/output/qc_observations.csv",
    ],
}


def default_out_root() -> Path:
    """Repo-relative default: `<step>/06.3_input/output` under the data root."""
    return DATA_ROOT / "stage2/06_prnet_transfer"


def output_dir(out_root: Path | str | None = None) -> Path:
    root = Path(out_root).resolve() if out_root else (DATA_ROOT / "stage2/06_prnet_transfer")
    return root / "06.3_input" / "output"


# --------------------------------------------------------------------------- #
# Hashing / git
# --------------------------------------------------------------------------- #

def sha256_bytes(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def sha256_text(text: str) -> str:
    return sha256_bytes(text.encode("utf-8"))


def sha256_file(path: Path | str, chunk: int = 8 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def sha256_json(obj) -> str:
    """Canonical digest of a JSON-serialisable object (sorted keys, no spaces)."""
    return sha256_text(json.dumps(obj, sort_keys=True, separators=(",", ":"),
                                  ensure_ascii=False))


def canonical_json(obj) -> str:
    return json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=False) + "\n"


def git_state(start: Path | None = None) -> dict:
    cwd = str(start or ROOT)
    try:
        sha = subprocess.run(["git", "-C", cwd, "rev-parse", "HEAD"],
                             capture_output=True, text=True, check=True).stdout.strip()
        branch = subprocess.run(["git", "-C", cwd, "rev-parse", "--abbrev-ref", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
        status = subprocess.run(["git", "-C", cwd, "status", "--porcelain"],
                                capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError) as exc:      # pragma: no cover
        raise ContractError(f"cannot read git state in {cwd}: {exc}") from exc
    return {"sha": sha, "branch": branch, "dirty": bool(status.strip()),
            "dirty_paths": len(status.splitlines())}


def runtime_env() -> dict:
    import anndata
    import numpy
    import pandas
    import sys
    return {
        "python": sys.version.split()[0],
        "executable": sys.executable,
        "anndata": anndata.__version__,
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
    }



# Human-readable meaning of every field a downstream ticket may read.  The
# generated `contract.json` joins these with the dtypes measured on the frozen
# artifacts, so a rename or a dtype change cannot pass unnoticed.
FIELD_MEANINGS = {
    # shared CM/Tahoe inherited columns
    "well": "observation key: CM = `temp_env::Plate<N>-<well>`, Tahoe = `tahoe::<plate>:<smp>:<cell_line>`",
    "dataset_id": "source id (`temp_env` | `tahoe`)",
    "compound_name": "normalised display name of the perturbagen",
    "drug_key": "desalted parent InChIKey (stereo retained); '' on control rows",
    "inchikey": "exact structure key",
    "smiles": "canonical SMILES (Tahoe: derived from the atlas SMILES)",
    "dose_um": "dose in micromolar; 0.0 on control rows",
    "dose_source": "provenance of the dose value",
    "time_h": "exposure time in hours; 24.0 for CM (loader constant), NaN for Tahoe",
    "is_control": "vehicle-control row",
    "vehicle": "vehicle identity; CM DMSO/Media, Tahoe DMSO_TF on controls and '' on treated",
    "plate": "physical plate id",
    "control_group": "CM `temp_env::Plate<N>` (MIXES DMSO+Media — never a matched baseline); Tahoe `tahoe::<plate>::<cell_line>`",
    "cell_type": "cell type (CM ipsc_cm; Tahoe NaN — see cell_line_id)",
    "donor": "donor/line; CM literal 'NA', Tahoe NaN",
    "replicate": "GENERATED within-group index, not a raw replicate id",
    "atc": "ATC code (frozen annotation)",
    "moa": "mechanism of action (frozen annotation)",
    "mesh_classes": "PubChem MeSH classes",
    "pharm_atc_moa": "old automatic atc-moa tier membership (annotation, not eligibility)",
    "pharm_mesh": "old automatic atc-moa-mesh tier membership",
    "annotation_source": "annotation provenance",
    "lib_size": "X row sum",
    "n_genes_detected": "(X > 0) row sum",
    "cell_line_id": "Tahoe cell line (Cellosaurus id)",
    "cell_name": "Tahoe cell line name",
    "n_cells": "cells aggregated into this pseudobulk profile (not a technical replicate count)",
    "moa_broad": "atlas broad MoA",
    "moa_fine": "atlas fine MoA",
    "pubchem_cid": "PubChem CID",
    "has_structure": "structure resolvable; False = 2 Tahoe drugs / 300 treated rows kept for audit only",
    "has_matched_control": "inherited 04 flag (True on every row; use `matched_control_ok` instead)",
    "SMILES": "author-shim copy of smiles",
    "dose": "author-shim copy of dose_um",
    "paired_control_index": "author-shim pointer to the first control row of the row's control group",
    # CM added
    "observation_id": "unique observation key (= obs index `well`)",
    "physical_well": "physical well token (CM `B2`; Tahoe `smp_1495`)",
    "plate_no": "CM plate number 1-16",
    "replicate_raw": "raw design Replicate_Number re-joined from the GSE262419 hash table",
    "solvent": "raw design Solvent (CM)",
    "chemical_index": "raw design Chemical_Index (CM)",
    "eligibility_decision": "ticket #1 decision: include | exclude | uncertain (NaN on control rows)",
    "elig_class": "ticket #1 decision class (approved | clinical | mechanism_tool | other_use | unknown)",
    "qc_pass": "ticket #2 per-observation QC verdict",
    "qc_reason": "ticket #2 reason codes (';'-joined, first is primary)",
    "matched_control_ok": "ticket #2 matched-control availability; a control row is its own control",
    "matched_control_key": "ticket #2 matched-control group key (CM plate x vehicle; Tahoe plate x cell line)",
    "main_set": "treated AND eligibility include AND qc_pass — the fold/scored population",
    "expansion_set": "treated AND eligibility exclude AND qc_pass — ticket #11 expansion training only, never scored",
    # Tahoe added
    "tahoe_split": "train | valid by plate13/plate14 (technical default, NOT a cold-drug split)",
    "tahoe_valid_group": "validation grouping key: drug_key x dose x cell_line (controls: control x plate x cell_line)",
    "tahoe_usable": "~is_control AND has_structure AND qc_pass AND matched_control_ok",
    "time_source": "literature provenance for the 24 h claim; time_h itself stays NaN",
    # var
    "symbol": "HGNC symbol, the frozen gene-axis key/order",
    "hgnc_id": "HGNC id from corpus/v1 genes.csv",
    "entrez_id": "Entrez id from corpus/v1 genes.csv",
    "ensembl_id": "Ensembl id from corpus/v1 genes.csv",
    "match_type": "gene-identity match provenance from corpus/v1",
    "gene_measured": "per-source measurement mask; False only for the FAM153B Tahoe slot",
    "n_nonzero": "non-zero cells in this artifact's X column (audit only; NOT a measurement test)",
    # side tables
    "fold": "fold index 0-4",
    "role": "train | valid | test within that fold",
    "cluster_id": "identity cluster (drug_key + inchikey connectivity block); never split across folds",
    "test_fold": "the fold in which this drug is test",
    "valid_fold": "the fold in which this drug is valid",
    "n_wells": "main-set wells of this drug (on this plate, for folds.parquet)",
    "n_doses": "distinct doses of this drug (on this plate)",
    "n_plates": "plates carrying this drug's main-set wells (`plate` is ';'-joined)",
    "plate_weight": "consensus plate weight = 1 / n_plates",
    "n_wells_plate": "QC-pass DMSO candidate rows on that plate in that fold",
    "weight": "consensus well weight = plate_weight / n_wells_plate; sums to 1 per fold",
    "n_controls": "control observations in the group",
    "n_qc_pass": "QC-pass control observations in the group (0 => treated rows of this group are not usable)",
    "row": "0-based row position in the artifact (cm_transfer.h5ad for consensus membership; tahoe_transfer.h5ad for control groups and pairs)",
    "group_index": "integer control-group id, shared by tahoe_control_groups.parquet and tahoe_pairs.npz",
    "group_key": "canonical control group key `plate|cell_line_id` (never crosses cell lines)",
    "usable": "tahoe_usable of that treated row",
    "split": "tahoe_split of that treated row",
}

# --------------------------------------------------------------------------- #
# Reading the 16 GB h5ad files without materialising their matrices
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class H5adHeader:
    """The metadata of an h5ad artifact, decoded without reading X or layers."""

    obs: pd.DataFrame
    var: pd.DataFrame
    uns: dict
    shape: tuple[int, int]
    x_dtype: str
    x_nnz: int
    layer_names: list[str]


def _read_elem(payload):
    """`anndata.io.read_elem` on 0.11+, `anndata.experimental.read_elem` before."""
    try:                                            # anndata >= 0.11
        from anndata.io import read_elem
    except ImportError:                             # pragma: no cover - 0.10 reader
        from anndata.experimental import read_elem
    return read_elem(payload)


def read_h5ad_header(path: Path | str, *, with_uns: bool = True) -> H5adHeader:
    """Open a large h5ad with h5py and decode only obs/var/uns/shape.

    ``anndata.read_h5ad(..., backed="r")`` is **not** safe here: backed mode
    still materialises ``layers``, which is 8.1 GB for the frozen Tahoe view
    (67,018 x 18,164, 1.01e9 nnz) and 12.9 GB for the 25.8 GB corpus (1.61e9
    nnz).  That alone exceeds the WSL2 memory envelope and is what killed the
    2026-11-02 session.  Reading the element groups directly keeps the peak
    proportional to the metadata (a few MB), and `X` is never touched.
    """
    import h5py

    path = Path(path)
    with h5py.File(path, "r") as handle:
        parts: dict = {name: _read_elem(handle[name]) if name in handle else None
                       for name in ("obs", "var")}
        uns = _read_elem(handle["uns"]) if (with_uns and "uns" in handle) else {}
        matrix = handle["X"]
        if isinstance(matrix, h5py.Group):
            shape = tuple(int(v) for v in matrix.attrs["shape"])
            x_dtype = str(matrix["data"].dtype)
            x_nnz = int(matrix["indptr"][-1])
        else:                                       # dense encoding
            shape = tuple(int(v) for v in matrix.shape)
            x_dtype = str(matrix.dtype)
            x_nnz = int(shape[0] * shape[1])
        layer_names = sorted(handle["layers"].keys()) if "layers" in handle else []
    obs, var = parts["obs"], parts["var"]
    if obs is None or var is None:
        raise ContractError(f"{path}: h5ad has no obs/var group")
    return H5adHeader(obs=obs, var=var, uns=dict(uns) if uns else {}, shape=shape,
                      x_dtype=x_dtype, x_nnz=x_nnz, layer_names=layer_names)


# --------------------------------------------------------------------------- #
# Gene axis
# --------------------------------------------------------------------------- #

REQUIRED_GENE_COLUMNS = ["symbol", "hgnc_id", "entrez_id", "ensembl_id", "match_type"]


def load_gene_axis() -> pd.DataFrame:
    """The frozen 18,164-gene axis, in corpus order, index = symbol."""
    genes = pd.read_csv(CORPUS_GENES, dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED_GENE_COLUMNS if c not in genes.columns]
    if missing:
        raise ContractError(f"{CORPUS_GENES}: missing columns {missing}")
    genes = genes[REQUIRED_GENE_COLUMNS].copy()
    if len(genes) != N_GENES:
        raise ContractError(f"gene axis has {len(genes)} rows, expected {N_GENES}")
    if genes["symbol"].duplicated().any():
        raise ContractError("gene axis has duplicate symbols")
    if genes["symbol"].isna().any() or (genes["symbol"] == "").any():
        raise ContractError("gene axis has empty symbols")
    out = genes.set_index("symbol")
    out.index.name = "symbol"
    return out


def gene_measured_cm(axis: pd.DataFrame) -> pd.Series:
    """CM panel membership from the frozen per-source resolved gene set.

    The axis is the four-source intersection, so every axis gene is in the
    temp_env panel — asserted from `source_gene_sets.csv`, not assumed.  Two CM
    columns are all-zero (CDCA2, TEX38) yet *measured*; zero counts are not
    absence.
    """
    sets = pd.read_csv(CORPUS_SOURCE_GENE_SETS, dtype=str, keep_default_na=False)
    if "dataset_id" not in sets.columns or "symbol" not in sets.columns:
        raise ContractError(f"{CORPUS_SOURCE_GENE_SETS}: unexpected columns")
    panel = set(sets.loc[sets["dataset_id"] == "temp_env", "symbol"])
    if len(panel) < N_GENES:
        raise ContractError(f"temp_env resolved gene set has {len(panel)} symbols")
    return pd.Series([sym in panel for sym in axis.index], index=axis.index, name="gene_measured")


def tahoe_gene_space_ensembl() -> set[str]:
    """The Tahoe corpus's native 62,710-gene space, as Ensembl ids (var only)."""
    if not TAHOE_CORPUS.exists():
        raise ContractError(f"missing Tahoe corpus: {TAHOE_CORPUS}")
    var = read_h5ad_header(TAHOE_CORPUS, with_uns=False).var
    if "ensembl_id" not in var.columns:
        raise ContractError(f"{TAHOE_CORPUS}: var has no ensembl_id")
    return {str(x) for x in var["ensembl_id"].tolist()}


def gene_measured_tahoe(axis: pd.DataFrame, ensembl_space: set[str]) -> pd.Series:
    """Tahoe measurement mask = Ensembl membership in the atlas gene space.

    NOTE: 18 Tahoe columns are all-zero although measured (e.g. NKX2-1); only
    FAM153B is *absent* from the atlas.  Measurement is membership, not counts.
    """
    return pd.Series([str(e) in ensembl_space for e in axis["ensembl_id"]],
                     index=axis.index, name="gene_measured")


def gene_space_digest(ensembl_ids) -> str:
    """Digest of a gene space: sha256 over sorted, newline-joined Ensembl ids."""
    return sha256_text("\n".join(sorted(str(x) for x in ensembl_ids)))


# --------------------------------------------------------------------------- #
# Frozen source records (content hashes come from the upstream frozen manifests)
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class SourceRecord:
    role: str
    path: Path
    content_sha256: str | None       # from the upstream frozen manifest
    bytes: int
    note: str = ""
    # "file": content_sha256 is sha256 over the file bytes (re-hashable here)
    # "content": the upstream digest hashes the h5ad arrays, not the file bytes
    # "none": no upstream digest exists
    digest_kind: str = "file"

    def as_dict(self) -> dict:
        return {
            "role": self.role,
            "path": relpath(self.path),
            "content_sha256": self.content_sha256,
            "digest_kind": self.digest_kind,
            "bytes": self.bytes,
            "note": self.note,
        }


def companion_handoff(path: Path, names=("handoff.json", "manifest.json")) -> Path | None:
    """The machine-readable companion document next to a handoff table, if any."""
    for name in names:
        candidate = Path(path).expanduser().resolve().parent / name
        if candidate.exists():
            return candidate
    return None


def read_json(path: Path) -> dict:
    import json as _json
    return _json.loads(Path(path).read_text(encoding="utf-8"))


def _frozen_manifest_row(path: Path, key: str, value: str) -> dict:
    table = pd.read_csv(path, dtype=str, keep_default_na=False)
    row = table.loc[table[key] == value]
    if len(row) != 1:
        raise ContractError(f"{path}: expected exactly one {key}={value!r} row")
    return row.iloc[0].to_dict()


def source_records() -> list[SourceRecord]:
    """The exact frozen bytes this build reads, with adjudicated digests."""
    cm_row = _frozen_manifest_row(CM_BUNDLE_MANIFEST, "view", "temp_env_all-smiles.h5ad")
    th_row = _frozen_manifest_row(TAHOE_VIEW_MANIFEST, "view", "tahoe_cm_schema")
    records = [
        SourceRecord("cm_view", CM_VIEW, str(cm_row["content_sha256"]),
                     int(cm_row["bytes"]),
                     "02 datasets-1 content_sha256 over X + layers",
                     digest_kind="content"),
        SourceRecord("tahoe_view", TAHOE_VIEW, str(th_row["content_sha256"]),
                     int(th_row["bytes"]),
                     "04 pretrain_view content_sha256 over CSR data/indices/indptr",
                     digest_kind="content"),
        SourceRecord("gene_axis", CORPUS_GENES, sha256_file(CORPUS_GENES),
                     CORPUS_GENES.stat().st_size,
                     "corpus v1 genes.csv; symbol order is var everywhere"),
        SourceRecord("source_gene_sets", CORPUS_SOURCE_GENE_SETS,
                     sha256_file(CORPUS_SOURCE_GENE_SETS),
                     CORPUS_SOURCE_GENE_SETS.stat().st_size,
                     "frozen per-source resolved gene sets (CM panel membership)"),
        SourceRecord("gse262419_hash", GSE262419_HASH, sha256_file(GSE262419_HASH),
                     GSE262419_HASH.stat().st_size,
                     "raw temp_env design table (replicate/solvent/chemical index)"),
        SourceRecord("tahoe_corpus", TAHOE_CORPUS, None, TAHOE_CORPUS.stat().st_size,
                     "62,710-gene space for the Tahoe measurement mask; the 04 "
                     "manifest is per-profile, so no file content hash exists — "
                     "the derived gene space is digested instead",
                     digest_kind="none"),
    ]
    for rec in records:
        if not rec.path.exists():                      # pragma: no cover
            raise ContractError(f"missing frozen source {rec.role}: {rec.path}")
    return records


# --------------------------------------------------------------------------- #
# Upstream handoffs (#1 eligibility, #2 per-observation QC)
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Handoff:
    kind: str
    path: Path
    sha256: str
    bytes: int
    rows: int
    note: str = ""
    detail: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {"kind": self.kind, "path": relpath(self.path),
                "sha256": self.sha256, "bytes": self.bytes, "rows": self.rows,
                "note": self.note, "detail": self.detail}


def _current_pointer(directory: Path) -> Path | None:
    """Resolve a versioned artifact set through its `current.json` pointer.

    ADR 0030: a new rule version gets a new artifact identity, the previous
    record is never overwritten, and downstream steps read the pointer instead
    of a version directory.  Ticket #1 writes ``<output>/current.json`` with
    ``table: {path, sha256}``; resolving through it keeps a future rule revision
    consumable and lets this build fail closed if the table drifts from the
    pointer's declared digest.
    """
    pointer = Path(directory) / "current.json"
    if not pointer.exists():
        return None
    block = (read_json(pointer).get("table") or {})
    declared = str(block.get("path") or "")
    if not declared:
        return None
    path = (pointer.parent / declared).expanduser().resolve()
    if not path.exists():
        raise ContractError(f"{pointer} points at a missing table: {path}")
    expected = block.get("sha256")
    if expected and sha256_file(path) != str(expected):
        raise ContractError(
            f"{path} does not match the sha256 recorded in {pointer}; the "
            "versioned handoff and its pointer disagree")
    return path


def current_pointer_record(path: Path | str) -> dict | None:
    """The `current.json` pointer governing `path`, if it lives in a version set."""
    path = Path(path).expanduser().resolve()
    pointer = path.parent.parent / "current.json"
    if not pointer.exists():
        return None
    document = read_json(pointer)
    return {"path": relpath(pointer), "sha256": sha256_file(pointer),
            "artifact_version": document.get("artifact_version"),
            "rule": document.get("rule"), "handoff": document.get("handoff")}


def resolve_handoff(kind: str, override: str | Path | None = None) -> Path:
    if override:
        path = Path(override).expanduser().resolve()
        if not path.exists():
            raise ContractError(f"{kind} handoff not found: {path}")
        return path
    # Versioned artifact sets first: `current.json` is the authoritative pointer.
    searched_dirs: list[Path] = []
    for rel in HANDOFF_CANDIDATES[kind]:
        directory = (STEP_DATA / rel).parent
        if directory in searched_dirs:
            continue
        searched_dirs.append(directory)
        resolved = _current_pointer(directory)
        if resolved is not None:
            return resolved
    for rel in HANDOFF_CANDIDATES[kind]:
        candidate = STEP_DATA / rel
        if candidate.exists():
            return candidate
    searched = "\n  ".join(str(STEP_DATA / rel) for rel in HANDOFF_CANDIDATES[kind])
    raise ContractError(
        f"{kind} handoff not found. Searched:\n  {searched}\n"
        f"Pass --{kind.replace('_', '-')} PATH or have ticket #1/#2 write one.")


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        return pd.read_parquet(path)
    # pandas decompresses .gz transparently; `compression="infer"` keeps both
    # a plain table and #2's `observations.csv.gz` on the same path.
    if suffix in (".csv", ".tsv", ".txt", ".gz"):
        stem = path.name[:-3].lower() if suffix == ".gz" else path.name.lower()
        sep = "\t" if stem.endswith((".tsv", ".txt")) else ","
        return pd.read_csv(path, sep=sep, keep_default_na=False, dtype=str,
                           compression="infer")
    raise ContractError(f"unsupported handoff format: {path}")


def _require_columns(table: pd.DataFrame, required, path: Path, kind: str) -> None:
    missing = [c for c in required if c not in table.columns]
    if missing:
        raise ContractError(
            f"{kind} handoff {path} is missing column(s) {missing}; "
            f"present: {list(table.columns)}")


def names_present(table: pd.DataFrame, names) -> str | None:
    for name in names:
        if name in table.columns:
            return name
    return None


def _first_present(table: pd.DataFrame, names, what: str, path: Path) -> str:
    for name in names:
        if name in table.columns:
            return name
    raise ContractError(f"{path}: no {what} column among {names}")


ELIGIBILITY_ALIASES = {
    "drug_key": ["drug_key", "drugkey", "drug", "parent_inchikey"],
    "eligibility_decision": ["eligibility_decision", "decision", "eligibility"],
    "elig_class": ["elig_class", "eligibility_class", "decision_class", "class",
                   "drug_class"],
}


def load_eligibility(path: Path | str | None = None) -> tuple[pd.DataFrame, Handoff]:
    """#1 handoff -> (drug_key, eligibility_decision, elig_class), validated.

    Semantics frozen here: `eligibility_decision` is exactly one of
    include|exclude|uncertain; it is a *pharmacology* judgement, independent of
    experimental QC and of transcriptional response strength.
    """
    resolved = resolve_handoff("eligibility", path)
    raw = _read_table(resolved)
    names = {}
    for canon, aliases in ELIGIBILITY_ALIASES.items():
        names[canon] = _first_present(raw, aliases, canon, resolved)
    table = pd.DataFrame({
        "drug_key": raw[names["drug_key"]].astype(str).str.strip(),
        "eligibility_decision": raw[names["eligibility_decision"]].astype(str).str.strip().str.lower(),
        "elig_class": raw[names["elig_class"]].astype(str).str.strip(),
    })
    if (table["drug_key"] == "").any():
        raise ContractError(f"{resolved}: empty drug_key row(s)")
    if table["drug_key"].duplicated().any():
        dupes = table.loc[table["drug_key"].duplicated(), "drug_key"].unique()[:5]
        raise ContractError(f"{resolved}: duplicated drug_key(s): {list(dupes)}")
    bad = sorted(set(table["eligibility_decision"]) - set(ELIGIBILITY_DECISIONS))
    if bad:
        raise ContractError(
            f"{resolved}: eligibility_decision values outside "
            f"{ELIGIBILITY_DECISIONS}: {bad}")
    detail: dict = {}
    companion = companion_handoff(resolved, ("handoff.json",))
    if companion is not None:
        doc = read_json(companion)
        detail = {
            "handoff_json": relpath(companion),
            "handoff_json_sha256": sha256_file(companion),
            "schema": doc.get("schema"),
            "rule_version": (doc.get("rule") or {}).get("version"),
            "rule_sha256": (doc.get("rule") or {}).get("sha256"),
            "counts": doc.get("counts"),
            "decision_class_counts": doc.get("decision_class_counts"),
        }
        declared = doc.get("counts") or {}
        observed = table["eligibility_decision"].value_counts().to_dict()
        for decision in ELIGIBILITY_DECISIONS:
            if decision in declared and int(declared[decision]) != int(
                    observed.get(decision, 0)):
                raise ContractError(
                    f"#1 handoff.json declares {declared[decision]} {decision} keys, "
                    f"the table has {observed.get(decision, 0)} (rule "
                    f"{(doc.get('rule') or {}).get('version')})")
        # A newer handoff may add an explicit study/main-set list.  If it does,
        # it must agree with the decision column — the decision is the authority.
        # handoff v2 lists `include`/`uncertain`/`exclude` as records with a
        # `drug_key`; older shapes use bare key lists.
        include_keys = set(table.loc[table["eligibility_decision"] == "include",
                                     "drug_key"])
        for field_name in ("include", "main_set", "study_set", "included_keys",
                           "include_keys"):
            declared_keys = doc.get(field_name)
            if not isinstance(declared_keys, list):
                continue
            if declared_keys and isinstance(declared_keys[0], dict):
                declared_keys = [row.get("drug_key") for row in declared_keys]
            declared_set = {str(k) for k in declared_keys if k is not None}
            if declared_set != include_keys:
                raise ContractError(
                    f"#1 handoff.json {field_name} disagrees with the decision "
                    f"column: only_table={sorted(include_keys - declared_set)[:5]} "
                    f"only_handoff={sorted(declared_set - include_keys)[:5]}")
            detail[f"handoff_{field_name}_reconciled"] = True
    # The versioned-artifact pointer (ADR 0030) is the provenance of *which*
    # handoff was consumed; record it so the manifest names the rule version.
    pointer = current_pointer_record(resolved)
    if pointer is not None:
        detail["current_pointer"] = pointer
    # An explicit boolean column in the table must also agree.
    for column in ("main_set", "study_set", "in_study"):
        if column in raw.columns:
            truthy = raw[column].astype(str).str.strip().str.lower().isin(
                {"true", "1", "yes"})
            expected = table["eligibility_decision"] == "include"
            if (truthy.to_numpy() != expected.to_numpy()).any():
                raise ContractError(
                    f"#1 eligibility table column {column!r} disagrees with the "
                    "decision column")
            detail[f"table_{column}_reconciled"] = True
    handoff = Handoff("eligibility", resolved, sha256_file(resolved),
                      resolved.stat().st_size, len(table),
                      "include/exclude/uncertain per drug_key", detail)
    return table, handoff


def _bool_series(column: pd.Series, path: Path) -> pd.Series:
    if column.dtype == bool:
        return column
    text = column.astype(str).str.strip().str.lower()
    mapping = {"true": True, "1": True, "yes": True, "pass": True, "passed": True,
               "false": False, "0": False, "no": False, "fail": False, "failed": False}
    unknown = sorted(set(text) - set(mapping))
    if unknown:
        raise ContractError(f"{path}: qc_pass values not boolean-like: {unknown[:8]}")
    return text.map(mapping).astype(bool)


QC_ALIASES = {
    "key": ["well", "observation_id", "obs_id", "observation", "row_id"],
    "qc_pass": ["qc_pass", "qc_status", "pass", "qc_ok", "kept"],
    "qc_reason": ["qc_reason", "reason", "qc_note", "note"],
    # #2's handoff also decides matched-control availability; carried when present.
    "matched_control_ok": ["matched_control_ok", "matched_control_available"],
    "matched_control_key": ["matched_control_key", "matched_control_group"],
}


def load_qc(path: Path | str, kind: str = "cm") -> tuple[pd.DataFrame, Handoff]:
    """#2 handoff -> (key, qc_pass, qc_reason), validated.

    `kind` is "cm" (keys = CM well ids) or "tahoe" (keys = Tahoe observation ids).
    A combined table with a `source` column is split by the caller and passed to
    `normalize_qc_table`.
    """
    resolved = Path(path).expanduser().resolve()
    raw = _read_table(resolved)
    table = normalize_qc_table(raw, resolved, kind)
    handoff = Handoff(f"{kind}_qc", resolved, sha256_file(resolved),
                      resolved.stat().st_size, len(table),
                      f"per-observation QC for {kind}")
    return table, handoff


def normalize_qc_table(raw: pd.DataFrame, resolved: Path, kind: str) -> pd.DataFrame:
    """Validate/normalise one per-source #2 QC table (see `load_qc`)."""
    key = _first_present(raw, QC_ALIASES["key"], "observation key", resolved)
    qc_col = _first_present(raw, QC_ALIASES["qc_pass"], "qc_pass", resolved)
    reason_col = _first_present(raw, QC_ALIASES["qc_reason"], "qc_reason", resolved)
    table = pd.DataFrame({
        "well": raw[key].astype(str).str.strip(),
        "qc_pass": _bool_series(raw[qc_col], resolved),
        "qc_reason": raw[reason_col].astype(str).str.strip(),
    })
    column = names_present(raw, QC_ALIASES["matched_control_ok"])
    if column is not None:
        # controls have no matched-control entry; only treated rows are judged.
        # Fail closed on an unknown token instead of silently reading it as NaN.
        values = raw[column].astype(str).str.strip().str.lower()
        known = {"true": True, "1": True, "yes": True,
                 "false": False, "0": False, "no": False, "": None}
        unknown = sorted(set(values) - set(known))
        if unknown:
            raise ContractError(
                f"{resolved}: matched_control_ok has unknown value(s) {unknown}; "
                "expected true/false or empty")
        table["matched_control_ok"] = values.map(known)
    column = names_present(raw, QC_ALIASES["matched_control_key"])
    if column is not None:
        table["matched_control_key"] = raw[column].astype(str).str.strip()
    if (table["well"] == "").any():
        raise ContractError(f"{resolved}: empty observation key row(s)")
    if table["well"].duplicated().any():
        dupes = table.loc[table["well"].duplicated(), "well"].unique()[:5]
        raise ContractError(f"{resolved}: duplicated observation id(s): {list(dupes)}")
    return table


def split_combined_qc(table: pd.DataFrame, path: Path) -> dict[str, pd.DataFrame]:
    """Split a combined per-observation QC table on its `source` column.

    #2's handoff is one `observations.csv.gz` with `source` ∈ {temp_env, tahoe}
    and `qc_status` ∈ {pass, fail}; either column name is accepted.
    """
    source_column = names_present(table, ["source", "dataset_id", "dataset"])
    if source_column is None:
        raise ContractError(f"{path}: combined QC table has no `source` column")
    labels = table[source_column].astype(str).str.strip().str.lower()
    mapping = {"cm": {"cm", "temp_env", "temp-env"}, "tahoe": {"tahoe"}}
    out = {}
    for kind, accepted in mapping.items():
        masks = labels.isin(accepted)
        if masks.any():
            out[kind] = table.loc[masks].reset_index(drop=True)
    return out


# --------------------------------------------------------------------------- #
# Locked config (spec §控制与模型, §切分与评价, §重复、预算与选择)
# --------------------------------------------------------------------------- #

def build_locked_config() -> dict:
    """The locked, hashed experiment config.  No result may redefine it.

    The **top level** is the canonical pipeline schema ticket #4's
    `ptlib.config.normalize_locked_config` reads (`structure`, `recipe`,
    `metric`, `loss`, `optimizer`, `prediction_noise`, `n_folds`, `seeds`,
    `lrs`, `dose_encoding`), so the pipeline binds to these values instead of
    falling back to its own defaults.  The decision record (folds, consensus,
    validation, budgets, derived defaults) is carried beside it.

    Every number is quoted from the approved spec/ADR or marked
    `derived_default` with a one-line rationale; nothing contradicts the spec.
    """
    config = {
        "schema": "06.3_input/locked_config/1",
        "version": CONFIG_VERSION,
        "contract_id": CONTRACT_ID,
        "spec": ".scratch/stage2-prnet-formal-transfer/spec.md (approved 2026-11-02)",
        "frozen_by": "ticket #3 (GitHub issue Moke110/arc#3)",
        "consumer_note": "top-level structure/recipe/metric/loss/optimizer/"
                         "prediction_noise/n_folds/seeds/lrs/dose_encoding are the "
                         "canonical pipeline schema (ticket #4); the blocks below are "
                         "the decision record and must not be re-tuned",
        # ---------------- canonical pipeline schema ----------------
        "structure": {
            "n_genes": N_GENES,
            "hidden": [128],
            "z_dim": 64,
            "c_dim": 64,
            "n_dim": 10,
            "drug_dim": 1024,
            "adaptor": [128],
            "dropout": 0.05,
        },
        "recipe": {
            "batch_size": 512,
            "max_epochs": 100,
            "min_epochs": 20,
            "early_stop_patience": 20,
            "sched_factor": 0.5,
            "sched_patience": 10,
            "min_lr": 1e-08,
            "weight_decay": 1e-08,
        },
        "metric": {
            "name": "direction_auprc",
            "headline": "direction AUPRC@200",
            "column": "AUPRC@200",
            "k": 200,
            "plateau_mode": "max",
            "selection_metric": "validation direction AUPRC@200 (max)",
            "true_set": "top-k up and bottom-k down of the condition's mean true delta; "
                        "stable argsort, ties by corpus gene order; a direction with < k "
                        "non-zero genes falls back to its non-zero support",
            "per_condition_ap": "AP_up = average_precision(y_up, pred); AP_down = "
                                "average_precision(y_down, -pred); condition score = mean "
                                "of the non-NaN directions",
            "condition": "CM: drug_key x dose_uM; Tahoe: drug_key x dose_uM x cell_line_id",
            "replicate_rule": "average predictions within a condition BEFORE scoring",
            "aggregation": "macro mean over conditions",
            "not_a_significance_test": True,
            "auxiliary": ["delta RMSE", "delta correlation",
                          "low-response / toxicity stratification (labels are all "
                          "`unknown` in this round)"],
            "test_usage": "scored once after candidate selection; never used for LR, "
                          "epoch, seed or capacity",
        },
        "loss": "gaussian_nll",
        "optimizer": "adam",
        "prediction_noise": "zero",
        "n_folds": N_FOLDS,
        "seeds": list(MODEL_SEEDS),
        "lrs": {"pretrain": 0.001, "cm": [0.0001, 0.0003, 0.001]},
        "dose_encoding": {
            "fingerprint": "FCFP4",
            "rdkit": "AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, "
                     "useFeatures=True, nBits=1024)",
            "n_bits": 1024,
            "radius": 2,
            "use_features": True,
            "dose_scale": "log10(dose_uM+1)",
            "label": "FCFP4_bits * log10(dose_uM + 1)",
            "binning": "none; the two sources' dose grids are not aligned",
            "documented_scalings": {"0.05": 0.02119, "0.2": 0.07918, "0.5": 0.17609,
                                    "1.0": 0.30103, "5.0": 0.77815, "10.0": 1.04139},
        },
        # ---------------- decision record ----------------
        "gene_axis": {
            "n_genes": N_GENES,
            "order": "stage2/02_cm_data/output/corpus/v1/genes.csv (HGNC symbol, lexicographic)",
            "identity_columns": ["hgnc_id", "entrez_id", "ensembl_id", "match_type"],
            "measurement_mask": "var.gene_measured, per source; never a measured zero",
        },
        "model_detail": {
            "architecture": "PRnet PGM (author structure, re-implemented in-file)",
            "in_dim": N_GENES,
            "comb_num": 1,
            "params_expected": 386 * N_GENES + 173696,
            "params_formula": "386*n_genes + 173696 (measured 7,185,000 at 18,164)",
            "objective": "expression target: predict treated log1p_cp10k from "
                         "(control, structure, dose)",
            "output_nonnegative_mean": True,
            "signed_delta_head": False,
            "finetune": "full-parameter; no frozen layers",
            "plateau_metric": "validation direction AUPRC@200",
            "plateau_mode": "max",
            "cm_lr_candidates": [0.0001, 0.0003, 0.001],
            "pretrain_lr": 0.001,
        },
        "loss_detail": {
            "implementation": "F.gaussian_nll_loss(mean, target, var, reduction='mean')",
            "clamp": False,
            "note": "a log-likelihood may be negative; never clamp or abs it",
        },
        "data_contract": {
            "cm_rows": "QC-pass treated wells with a resolvable drug_key, plus QC-pass DMSO/Media control wells",
            "cm_main_set": "eligibility_decision == include AND qc_pass",
            "cm_expansion_set": "eligibility_decision == exclude AND qc_pass (ticket #11 training only)",
            "cm_target": "layers['log1p_cp10k'] (row-wise log1p(CP10k))",
            "cm_train_control": "mean over QC-pass DMSO rows on the same plate (never the mixed-vehicle control_mean)",
            "cm_validation_test_input": "one fixed per-fold consensus c* (consensus_controls.npz); never a per-plate control",
            "eval_delta": "y - mean(log1p_cp10k of QC-pass matched same-plate DMSO rows), computed only at scoring time",
            "tahoe_target": "layers['log1p_cp10k']",
            "tahoe_control": "mean of the same (plate, cell_line) group's QC-pass control rows; never across cell lines",
            "tahoe_rows_used": "~is_control & has_structure & qc_pass & matched_control_ok (tahoe_usable); time never filters rows",
            "deployment_semantics": "validation/test/inference receive one fixed c*; the experimental truth is computed separately at scoring time",
        },
        "folds": {
            "n_folds": N_FOLDS,
            "split_unit": "drug_key",
            "seed": FOLD_SEED,
            "ratios": "~6:2:2 drugs per fold (test:valid:train)",
            "cluster_rule": "connected components of drug_key <-> inchikey[:14] (salt/connectivity cannot straddle)",
            "assignment": "greedy, capacity = ceil(n_main_drugs / n_folds) drugs per fold",
            "valid_fold": "(test_fold + n_folds - 1) % n_folds",
            "roles": list(FOLD_ROLES),
            "excluded_drugs": "eligibility exclude/uncertain receive NO fold (reserved for the #11 expansion ablation)",
            "plate_confounding": "459/460 temp_env keys sit on one plate; a fold's test drugs bring their plates. A drug split is NOT a plate split.",
        },
        "consensus_control": {
            "name": "CM consensus control c*",
            "per_fold": True,
            "candidate_wells": "CM control rows with vehicle == 'DMSO' AND qc_pass on plates that carry >= 1 TRAINING drug of the fold",
            "held_out_environment_rule": "a plate whose fold drugs are all test/valid contributes no control row",
            "plate_mean": "per-plate mean of log1p_cp10k over candidate wells",
            "cross_plate_weight": "equal weight per plate (1 / n_plates)",
            "well_weight": "plate_weight / n_wells_in_plate",
            "artifacts": ["consensus_controls.npz", "consensus_membership.parquet"],
            "forbidden": "test/valid drug treated expression; held-out plates' control rows; adjusting weights after seeing test data",
        },
        "validation": {
            "cm": "the fold's valid roles only; model input is the fold's fixed c*",
            "tahoe_pretrain": "plate13 + plate14 valid (existing technical default; not a cold-drug claim)",
            "tahoe_conditions": "drug_key x dose_uM x cell_line_id",
            "tahoe_matched_control": "same (plate, cell_line) group only",
        },
        "seeds_detail": {
            "derived_default": "three consecutive seeds in the repo's 20261021 family; "
                               "the spec fixes the count (3), not the values",
            "no_seed_cherry_picking": True,
            "fold_seed": FOLD_SEED,
            "fold_seed_derived_default": "study start date; the spec fixes only that the "
                                         "seed is recorded",
        },
        "lr_selection": {
            "rule": "per initialisation and per fold, choose the LR maximising the mean "
                    "over the three seeds of the best validation AUPRC@200",
            "tie_break": "lower LR first, then candidate order",
            "tie_break_derived_default": "the spec fixes the mean-of-three-seeds rule, "
                                         "not the tie-break",
            "then_evaluate": "the three checkpoints (one per seed) at the selected LR",
        },
        "test_units": {
            "definition": "one (init, fold, seed) selected checkpoint at the fold's chosen LR",
            "count": 30,
            "factors": "2 inits x 5 folds x 3 seeds",
            "policy": "all 30 reported; no seed or fold selection; paired by drug_key",
            "candidates": 90,
            "candidates_formula": "2 inits x 5 folds x 3 seeds x 3 LRs",
        },
        "prediction": {
            "train_noise": "N(0,1) noise vector of width 10 per sample, drawn from a "
                           "seeded generator recorded per run",
            "validation_test_inference_noise": "exactly zero",
            "deterministic": True,
            "repeat_inference_policy": "byte-identical predictions for the same "
                                       "checkpoint and input",
        },
        "environments": {
            "builder": "/home/chang/miniconda3/envs/arc/bin/python (anndata 0.12.16)",
            "reader_trainer": "/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python "
                              "(anndata 0.10.8, torch 2.0.1+cu118)",
            "note": "the builder writes the h5ad; the pinned trainer must read it "
                    "(anndata 0.10.8 cannot write a dataframe/0.2.0 uns frame but "
                    "reads one)",
        },
        "derived_defaults": [
            {"key": "folds.seed", "value": FOLD_SEED,
             "rationale": "the spec requires a recorded seed, not a specific value"},
            {"key": "folds.assignment", "value": "greedy capacity ceil(n/5)",
             "rationale": "realises ~6:2:2 over clusters without a solver"},
            {"key": "seeds", "value": list(MODEL_SEEDS),
             "rationale": "the spec fixes three seeds, not their values"},
            {"key": "lr_selection.tie_break", "value": "lower LR, then candidate order",
             "rationale": "the spec fixes the mean-of-three-seeds rule only"},
            {"key": "consensus_control.well_weight", "value": "plate_weight / n_wells_in_plate",
             "rationale": "equal plate weight is required; spelling it per well makes the pool auditable"},
            {"key": "recipe.min_lr", "value": 1e-08,
             "rationale": "the repo's existing plateau floor (05/03.1); the spec fixes "
                          "factor and patience"},
        ],
    }
    return config


def config_digest(config: dict) -> str:
    return sha256_json(config)
