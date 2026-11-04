"""Executable verification gates G1-G6 for the v2 dataset corpus (ticket #17).

``build_corpus.py --verify`` runs these gates **reading only from disk**: it
never loads a raw source for G1-G4/G6 and never writes into ``output/``.  The
only write outside ``execution/build_record/`` is gate G5's temporary corpus,
which is deleted afterwards.

======================  ===================================================
G1 structure            the five corpora: unique index; non-negative integer
                        ``X`` (dense ``int32`` x4, CSR ``int32`` Tahoe);
                        sampled ``log1p_cp10k`` formula; empty
                        ``uns/obsm/varm/obsp/varp``; no ``pert_idx`` /
                        ``cell_idx`` / ``dose_idx``; the corpus's own obs/var
                        contract; rows / axis / chemicals against the
                        published contract.
G2 counts               the corpus-coverage table in ``output/corpus/v2/
                        README.md`` re-parsed and compared item by item with
                        ``build_record/corpus_coverage.json`` *and* with the
                        counts recomputed from each corpus's own ``obs``.
G3 crosswalk            every ``var`` id of all five corpora resolves in
                        ``tables/gene_crosswalk.csv``; the 18,164
                        ``tables/genes.csv`` symbols are a subset of every CM
                        corpus's ``var``; ``tables/compounds.csv`` covers every
                        non-control/non-blank row's ``drug_key``.
G4 v1 reconstruction    rebuild the five frozen v1 views *from the v2
                        corpora* with the archived writer
                        ``stage2/_frozen/cm-data-v1/execution/views.py`` and
                        require each v1 content digest to equal
                        ``_frozen/cm-data-v1/datasets-1/manifest.csv`` (5/5).
G5 reproducibility      rebuild the smallest corpus twice into a temp dir:
                        identical file sha256 plus the recorded
                        content/file sha256.
G6 path audit           ``script/check_stage2_paths.sh`` reports 0 stale
                        retired step-path literals in tracked ``stage2/``
                        executable code.
======================  ===================================================

The G4 recipe (proven in ``.scratch/stage2-corpus-reorg/impl-notes-g4-proof.md``)
feeds the frozen writer a ``SourceData``-shaped proxy built from a corpus:

* ``counts`` = corpus ``X``;
* ``gene_map`` = an identity map against the frozen 18,164 ``genes.csv`` axis
  (each corpus ``var`` ``hgnc_id`` -> its position in that axis);
* ``meta`` = corpus ``obs`` verbatim, plus the v1-only fields reconstructed:
  ``entity_type`` (``small_molecule`` for non-control/non-blank rows, ``blank``
  for blanks), ``pharm_all_smiles`` joined from ``tables/compounds.csv`` on
  ``drug_key`` (the other two tier booleans are already in ``obs``), and
  ``dataset_id`` remapped back over the rename table
  (``sctox->screen_ctox``, ``tenv->temp_env``; dtoxs ids unchanged).  ``well``
  stays verbatim.  The writer recomputes ``control_group``/``replicate``
  exactly as v1 did.

The archived module keeps its ``from .drugs import ...`` / ``from .genes
import ...`` relative imports (see ``_frozen/cm-data-v1/execution/README.md``),
so it is executed **read-only by path** in a namespace whose ``__package__`` is
the live ``cmdata`` package: the relative imports then resolve to the live
``cmdata.drugs`` / ``cmdata.genes`` (``TIER_LABELS``, ``CompoundTable``,
``GeneUniverse`` are unchanged by the corpus reorganisation).  The frozen file
is never imported as a module, copied or edited.
"""

from __future__ import annotations

import gc
import hashlib
import json
import platform
import re
import resource
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
import types
from dataclasses import dataclass, field
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

from . import paths, tables
from .annotation import load_assets
from .corpus_h5ad import build_corpus, content_digest, write_corpus
from .drugs import build_compound_table
from .genes import GeneUniverse
from .sources import load_source, load_dtoxs_entrez_map

FROZEN_CM_DIR = paths.REPO_ROOT / "stage2" / "_frozen" / "cm-data-v1"
FROZEN_MANIFEST = FROZEN_CM_DIR / "datasets-1" / "manifest.csv"
FROZEN_GENES = FROZEN_CM_DIR / "corpus" / "v1" / "genes.csv"
FROZEN_VIEWS_PY = FROZEN_CM_DIR / "execution" / "views.py"
FROZEN_PATHS_PY = FROZEN_CM_DIR / "execution" / "paths_v1.py"
FROZEN_VIEW_DIR = FROZEN_CM_DIR / "datasets-1"

GATE_ORDER = ("G1", "G2", "G3", "G4", "G5", "G6")
GATE_TITLES = {
    "G1": "structure",
    "G2": "counts",
    "G3": "crosswalk",
    "G4": "v1 reconstruction",
    "G5": "reproducibility",
    "G6": "path audit",
}
GATE_LABELS = {
    "G1": "G1 structure + counts",
    "G2": "G2 counts",
    "G3": "G3 crosswalk coverage",
    "G4": "G4 v1 reconstruction",
    "G5": "G5 reproducibility",
    "G6": "G6 no stale paths",
}
# Spec #12 contract table / ADR 0033 §2.  Tahoe's corpus has no EXPECTED_* dict
# in cmdata.paths because it was moved in by #14 rather than built by #15.
TAHOE_EXPECTED = {
    "rows": 67018, "perturbed": 65618, "control": 1400, "blank": 0,
    "chemicals": 379, "axis": 62710,
}
V1_DATASET_ID = {"sctox": "screen_ctox", "tenv": "temp_env"}
ALL_CORPORA = [*paths.CORPUS_IDS, "tahoe"]
SAMPLE_ROWS = 8
MIN_FREE_GIB = 12.0
COVERAGE_COLUMNS = [
    "corpus", "distinct cell lines", "chemicals (perturbed wells)",
    "distinct doses (µM)", "time window", "(cell line, chemical, dose)",
    "control wells", "blank wells",
]
COVERAGE_KEYS = [
    "rows", "perturbed", "control", "blank", "distinct_cell_lines",
    "chemicals", "distinct_doses", "combinations", "time_window",
]
TEMP_ENV_VIEWS = ("temp_env_atc-moa.h5ad", "temp_env_atc-moa-mesh.h5ad",
                  "temp_env_all-smiles.h5ad")


@dataclass
class GateCheck:
    gate: str
    title: str
    passed: bool
    detail: str
    data: dict = field(default_factory=dict)

    def payload(self) -> dict:
        return {"gate": self.gate, "title": self.title,
                "passed": bool(self.passed), "detail": self.detail,
                "data": self.data}


class VerifyError(RuntimeError):
    """A gate could not even be attempted (missing artifact etc.)."""


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def _iso_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _environment() -> dict:
    return {
        "python": platform.python_version(),
        "anndata": version("anndata"),
        "h5py": version("h5py"),
        "numpy": version("numpy"),
        "pandas": version("pandas"),
        "platform": platform.platform(),
    }


def _available_gib() -> float:
    with open("/proc/meminfo") as handle:
        for line in handle:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 1024 ** 2
    raise VerifyError("cannot read MemAvailable from /proc/meminfo")


def _rel(path: Path) -> str:
    return str(path.relative_to(paths.REPO_ROOT))


def _sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def _corpus_path(corpus_id: str) -> Path:
    if corpus_id == "tahoe":
        return paths.TAHOE_H5AD
    return paths.corpus_h5ad(corpus_id)


def _read_obs(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as handle:
        return read_elem(handle["obs"])


def _read_var(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as handle:
        return read_elem(handle["var"])


def _sample_rows(n_rows: int, n_sample: int = SAMPLE_ROWS) -> list[int]:
    """Deterministic row sample: first, last and evenly spaced rows."""
    if n_rows <= n_sample:
        return list(range(n_rows))
    step = max(1, n_rows // (n_sample - 1))
    rows = sorted({0, n_rows - 1, *range(0, n_rows, step)})
    while len(rows) > n_sample:
        rows.pop(len(rows) // 2)
    return rows


def _frame_counts(obs: pd.DataFrame,
                  chemical_column: str = "drug_key") -> dict:
    control = obs["is_control"].astype(bool).to_numpy()
    blank = (obs["is_blank"].astype(bool).to_numpy()
             if "is_blank" in obs.columns
             else np.zeros(len(obs), dtype=bool))
    perturbed = ~control & ~blank
    return {
        "rows": int(len(obs)),
        "perturbed": int(perturbed.sum()),
        "control": int(control.sum()),
        "blank": int(blank.sum()),
        "chemicals": int(obs.loc[perturbed, chemical_column]
                         .astype(str).nunique()),
    }


def _expect(problems: list[str], condition: bool, message: str) -> None:
    if not condition:
        problems.append(message)


# ---------------------------------------------------------------------------
# G1 — structure
# ---------------------------------------------------------------------------
def _g1_dense(corpus_id: str, path: Path) -> dict:
    problems: list[str] = []
    with h5py.File(path, "r") as handle:
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])
        x = handle["X"][:]
        layer = handle["layers"]["log1p_cp10k"]
        layer_dtype = layer.dtype
        # containers and integer encodings (read from the file, not in memory)
        _expect(problems, not list(handle["uns"].keys()),
                f"{corpus_id}: uns is not empty")
        for container in ("obsm", "varm", "obsp", "varp"):
            _expect(problems, not list(handle[container].keys()),
                    f"{corpus_id}: {container} is not empty")
        # unique row index
        wells = obs.index.astype(str)
        _expect(problems, bool(obs.index.is_unique),
                f"{corpus_id}: row index is not unique")
        _expect(problems, not (wells == "").any(),
                f"{corpus_id}: empty row index")
        _expect(problems, list(wells) == obs["well"].astype(str).tolist(),
                f"{corpus_id}: index is not the well column")
        # X: dense non-negative int32
        _expect(problems, x.dtype == np.int32,
                f"{corpus_id}: X dtype is {x.dtype}, expected int32")
        negative = int((x < 0).sum())
        _expect(problems, negative == 0,
                f"{corpus_id}: {negative} negative X cell(s)")
        lib_size = x.sum(axis=1, dtype=np.int64)
        n_detected = (x > 0).sum(axis=1).astype(np.int32)
        _expect(problems, bool((lib_size > 0).all()),
                f"{corpus_id}: zero library size row")
        _expect(problems,
                (obs["lib_size"].to_numpy() == lib_size).all()
                and (obs["n_genes_detected"].to_numpy() == n_detected).all(),
                f"{corpus_id}: obs lib_size/n_genes_detected disagree with X")
        # layer formula on a row sample (never densify beyond the sample)
        rows = _sample_rows(len(obs))
        stored = layer[rows]
        expected = np.log1p(x[rows] / lib_size[rows, None] * 1e4)
        _expect(problems, layer_dtype == np.float32,
                f"{corpus_id}: layer dtype is {layer_dtype}, expected float32")
        _expect(problems,
                bool(np.isfinite(stored).all())
                and bool(np.allclose(stored, expected, rtol=1e-5, atol=1e-6)),
                f"{corpus_id}: log1p_cp10k formula fails on rows {rows}")
        # the corpus's own obs/var contract
        _expect(problems, list(obs.columns) == paths.CORPUS_OBS,
                f"{corpus_id}: obs columns differ from paths.CORPUS_OBS")
        _expect(problems, obs.index.name == "well",
                f"{corpus_id}: row index name is {obs.index.name!r}, "
                "expected 'well'")
        _expect(problems, list(var.columns) == [
            "hgnc_id", "entrez_id", "ensembl_id", "match_type"],
            f"{corpus_id}: var columns are {list(var.columns)}")
        _expect(problems, var.index.name == "symbol",
                f"{corpus_id}: var index name is {var.index.name!r}")
        _expect(problems, bool(var.index.is_unique),
                f"{corpus_id}: var index is not unique")
        _expect(problems, set(obs["dataset_id"].astype(str))
                == {corpus_id}, f"{corpus_id}: dataset_id is not "
                f"{ {corpus_id} }")
    forbidden = sorted({"pert_idx", "cell_idx", "dose_idx", "pert_dose_idx"}
                       & set(obs.columns))
    _expect(problems, not forbidden,
            f"{corpus_id}: forbidden integer encoding columns {forbidden}")
    counts = _frame_counts(obs)
    return {"problems": problems, "shape": [int(x.shape[0]), int(x.shape[1])],
            "dtype": str(x.dtype), "layer_dtype": str(layer_dtype),
            "index": "well", "sampled_rows": rows, "counts": counts}


def _tahoe_csr_row(handle, group: str, row: int) -> np.ndarray:
    indptr = handle[group]["indptr"]
    start, stop = int(indptr[row]), int(indptr[row + 1])
    data = handle[group]["data"][start:stop]
    indices = handle[group]["indices"][start:stop]
    dense = np.zeros(int(handle[group].attrs["shape"][1]), dtype=data.dtype)
    dense[indices] = data
    return dense


def _g1_tahoe(path: Path) -> dict:
    problems: list[str] = []
    with h5py.File(path, "r") as handle:
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])
        x_group = handle["X"]
        _expect(problems, x_group.attrs.get("encoding-type") == "csr_matrix",
                "tahoe: X is not CSR")
        shape = tuple(int(value) for value in x_group.attrs["shape"])
        data = x_group["data"]
        indptr = x_group["indptr"][:]
        _expect(problems, data.dtype == np.int32,
                f"tahoe: X data dtype is {data.dtype}, expected int32")
        _expect(problems, indptr.dtype == np.int64,
                f"tahoe: X indptr dtype is {indptr.dtype}, expected int64")
        _expect(problems, len(indptr) == shape[0] + 1,
                "tahoe: X indptr length is not n_rows + 1")
        _expect(problems, int(indptr[0]) == 0
                and int(indptr[-1]) == int(data.shape[0])
                and bool((np.diff(indptr) >= 0).all()),
                "tahoe: X indptr is not monotone/complete")
        # Stream the values: no dense materialisation of a 67k x 62k matrix.
        minimum = np.iinfo(np.int32).max
        step = 50_000_000
        for start in range(0, int(data.shape[0]), step):
            block = data[start:start + step]
            if block.size:
                minimum = min(minimum, int(block.min()))
        _expect(problems, minimum >= 0, f"tahoe: negative X value ({minimum})")
        layer_group = handle["layers"]["log1p_cp10k"]
        _expect(problems,
                layer_group.attrs.get("encoding-type") == "csr_matrix"
                and layer_group["data"].dtype == np.float32,
                "tahoe: layer is not a CSR float32 matrix")
        layer_indptr = layer_group["indptr"][:]
        _expect(problems, np.array_equal(layer_indptr, indptr),
                "tahoe: X and layer sparsity patterns differ")
        rows = _sample_rows(len(obs), n_sample=4)
        for row in rows:
            x_row = _tahoe_csr_row(handle, "X", row)
            stored = _tahoe_csr_row(handle, "layers/log1p_cp10k", row)
            lib = int(x_row.sum(dtype=np.int64))
            _expect(problems, lib > 0, f"tahoe: zero library size row {row}")
            if lib > 0:
                expected = np.log1p(x_row / lib * 1e4)
                _expect(problems,
                        bool(np.isfinite(stored).all())
                        and bool(np.allclose(stored, expected, rtol=1e-5,
                                             atol=1e-6)),
                        f"tahoe: log1p_cp10k formula fails on row {row}")
        _expect(problems, not list(handle["uns"].keys()),
                "tahoe: uns is not empty")
        for container in ("obsm", "varm", "obsp", "varp"):
            _expect(problems, not list(handle[container].keys()),
                    f"tahoe: {container} is not empty")
        _expect(problems, bool(obs.index.is_unique),
                "tahoe: row index is not unique")
        _expect(problems, obs.index.name == "profile_id",
                f"tahoe: row index name is {obs.index.name!r}")
        _expect(problems, list(obs.columns) == paths.TAHOE_OBS,
                f"tahoe: obs columns differ from paths.TAHOE_OBS")
        _expect(problems, list(var.columns) == ["ensembl_id", "token_id"],
                f"tahoe: var columns are {list(var.columns)}")
        _expect(problems, bool(var.index.is_unique),
                "tahoe: var index is not unique")
        _expect(problems, shape == (TAHOE_EXPECTED["rows"],
                                    TAHOE_EXPECTED["axis"]),
                f"tahoe: X shape {shape} != contract "
                f"{(TAHOE_EXPECTED['rows'], TAHOE_EXPECTED['axis'])}")
    forbidden = sorted({"pert_idx", "cell_idx", "dose_idx", "pert_dose_idx"}
                       & set(obs.columns))
    _expect(problems, not forbidden,
            f"tahoe: forbidden integer encoding columns {forbidden}")
    counts = _frame_counts(obs, chemical_column="drug")
    return {"problems": problems, "shape": list(shape), "dtype": "int32 (CSR)",
            "layer_dtype": "float32 (CSR)", "index": "profile_id",
            "sampled_rows": rows, "counts": counts}


def check_g1_structure() -> GateCheck:
    problems: list[str] = []
    data: dict = {}
    for corpus_id in ALL_CORPORA:
        path = _corpus_path(corpus_id)
        if not path.exists():
            problems.append(f"{corpus_id}: corpus missing ({_rel(path)})")
            continue
        entry = (_g1_tahoe(path) if corpus_id == "tahoe"
                 else _g1_dense(corpus_id, path))
        problems.extend(entry["problems"])
        published = (TAHOE_EXPECTED if corpus_id == "tahoe"
                     else paths.EXPECTED_CORPUS_COUNTS[corpus_id])
        pairs = {key: (entry["counts"].get(key), published[key])
                 for key in ("rows", "perturbed", "control", "blank",
                             "chemicals")}
        pairs["axis"] = (entry["shape"][1], published["axis"])
        mismatched = {key: pair for key, pair in pairs.items()
                      if pair[0] != pair[1]}
        if mismatched:
            problems.append(f"{corpus_id}: counts (computed, published) "
                            f"{mismatched} differ from the contract")
        entry.pop("problems")
        entry["published"] = published
        data[corpus_id] = entry
    passed = not problems
    detail = (f"{len(ALL_CORPORA)}/{len(ALL_CORPORA)} corpora: unique index, "
              "non-negative int32 X (dense x4, CSR tahoe), sampled "
              "log1p_cp10k formula, empty uns/obsm/varm/obsp/varp, no integer "
              "encodings, obs/var contract and published counts"
              if passed else "; ".join(problems))
    return GateCheck("G1", GATE_TITLES["G1"], passed, detail,
                     {"corpora": data, "problems": problems})


# ---------------------------------------------------------------------------
# G2 — counts
# ---------------------------------------------------------------------------
def _read_readme_table(header: list[str]) -> list[dict[str, str]]:
    """Parse the first Markdown table whose header row matches *header*."""
    lines = paths.CORPUS_V2_README.read_text().splitlines()
    header_line = "| " + " | ".join(header) + " |"
    rows: list[dict[str, str]] = []
    found = False
    for line in lines:
        if not found:
            found = line.strip() == header_line
            continue
        if not line.startswith("|"):
            break
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) != len(header):
            continue
        if all(set(cell) <= {"-", ":"} for cell in cells):
            continue
        rows.append(dict(zip(header, cells)))
    if not found:
        raise VerifyError(f"README coverage table not found: {header_line}")
    return rows


def _readme_number(cell: str) -> int:
    return int(cell.replace(",", ""))


def _coverage_frame(corpus_id: str) -> pd.DataFrame:
    obs = _read_obs(_corpus_path(corpus_id))
    if corpus_id == "tahoe" and "is_blank" not in obs.columns:
        obs = obs.copy()
        obs["is_blank"] = False
    return obs


def check_g2_counts() -> GateCheck:
    problems: list[str] = []
    if not paths.COVERAGE_RECORD.exists():
        raise VerifyError(f"coverage record missing: "
                          f"{_rel(paths.COVERAGE_RECORD)}")
    record = json.loads(paths.COVERAGE_RECORD.read_text())["corpora"]
    readme_rows = {row["corpus"].strip("`"): row
                   for row in _read_readme_table(COVERAGE_COLUMNS)}
    missing = sorted(set(ALL_CORPORA) - set(readme_rows))
    if missing:
        problems.append(f"README coverage table misses {missing}")
    data: dict = {}
    for corpus_id in ALL_CORPORA:
        entry = record.get(corpus_id, {})
        row = readme_rows.get(corpus_id)
        if row is None:
            continue
        # 1. README table vs the coverage record, item by item.
        pairs = {
            "distinct_cell_lines": (row["distinct cell lines"],
                                    entry.get("distinct_cell_lines")),
            "chemicals": (row["chemicals (perturbed wells)"],
                          entry.get("chemicals")),
            "distinct_doses": (row["distinct doses (µM)"],
                               entry.get("distinct_doses")),
            "combinations": (row["(cell line, chemical, dose)"],
                             entry.get("combinations")),
            "control": (row["control wells"], entry.get("control")),
            "blank": (row["blank wells"], entry.get("blank")),
        }
        for key, (readme_value, record_value) in pairs.items():
            if _readme_number(readme_value) != record_value:
                problems.append(
                    f"{corpus_id}: README {key} {readme_value} != "
                    f"corpus_coverage.json {record_value}")
        if row["time window"] != entry.get("time_window"):
            problems.append(f"{corpus_id}: README time window "
                            f"{row['time window']!r} != coverage record "
                            f"{entry.get('time_window')!r}")
        # 2. the corpus's own obs re-computed with the same counting rules.
        computed = tables.five_counts(_coverage_frame(corpus_id), corpus_id)
        for key in COVERAGE_KEYS:
            if computed[key] != entry.get(key):
                problems.append(
                    f"{corpus_id}: corpus obs {key}={computed[key]!r} != "
                    f"coverage record {entry.get(key)!r}")
        data[corpus_id] = {
            "readme": {record_key: row[readme_column] for
                       readme_column, record_key in (
                           ("distinct cell lines", "distinct_cell_lines"),
                           ("chemicals (perturbed wells)", "chemicals"),
                           ("distinct doses (µM)", "distinct_doses"),
                           ("(cell line, chemical, dose)", "combinations"),
                           ("control wells", "control"),
                           ("blank wells", "blank"),
                       )},
            "record": {key: entry.get(key) for key in COVERAGE_KEYS},
            "recomputed": {key: computed[key] for key in COVERAGE_KEYS},
        }
    passed = not problems
    detail = (f"{len(ALL_CORPORA)}/{len(ALL_CORPORA)} corpora: README "
              "coverage table = corpus_coverage.json = the corpora's own obs "
              "for every five-count item"
              if passed else "; ".join(problems[:5]))
    return GateCheck("G2", GATE_TITLES["G2"], passed, detail,
                     {"corpora": data, "problems": problems})


# ---------------------------------------------------------------------------
# G3 — crosswalk coverage
# ---------------------------------------------------------------------------
def check_g3_crosswalk() -> GateCheck:
    problems: list[str] = []
    genes_path = paths.TABLES_DIR / "genes.csv"
    crosswalk_path = paths.TABLES_DIR / "gene_crosswalk.csv"
    compounds_path = paths.TABLES_DIR / "compounds.csv"
    for path in (genes_path, crosswalk_path, compounds_path):
        if not path.exists():
            raise VerifyError(f"crosswalk table missing: {_rel(path)}")
    genes = pd.read_csv(genes_path, dtype=str, keep_default_na=False)
    gene_symbols = set(genes["symbol"].astype(str))
    if len(genes) != 18164 or len(gene_symbols) != 18164:
        problems.append(f"genes.csv is {len(genes)} rows / "
                        f"{len(gene_symbols)} symbols, expected 18,164")
    crosswalk = pd.read_csv(crosswalk_path, dtype=str,
                            keep_default_na=False)
    crosswalk["dataset_id"] = crosswalk["dataset_id"].astype(str)
    compounds = set(pd.read_csv(compounds_path, keep_default_na=False,
                                dtype={"drug_key": str})["drug_key"]
                    .astype(str))
    data: dict = {"genes_symbols": len(gene_symbols),
                  "crosswalk_rows": int(len(crosswalk)),
                  "compounds_drug_keys": len(compounds)}
    for corpus_id in paths.CORPUS_IDS:
        source_id = paths.CORPUS_TO_SOURCE[corpus_id]
        var = _read_var(_corpus_path(corpus_id))
        symbols = set(var.index.astype(str))
        hgnc_ids = set(var["hgnc_id"].astype(str))
        rows = crosswalk[crosswalk["dataset_id"] == source_id]
        rows_symbols = set(rows["symbol"].astype(str))
        rows_hgnc = set(rows["hgnc_id"].astype(str))
        missing_symbols = symbols - rows_symbols
        missing_hgnc = hgnc_ids - rows_hgnc
        if missing_symbols:
            problems.append(f"{corpus_id}: {len(missing_symbols)} var symbols "
                            f"without a crosswalk row, e.g. "
                            f"{sorted(missing_symbols)[:5]}")
        if missing_hgnc:
            problems.append(f"{corpus_id}: {len(missing_hgnc)} var hgnc_ids "
                            f"without a crosswalk row, e.g. "
                            f"{sorted(missing_hgnc)[:5]}")
        if not gene_symbols <= symbols:
            problems.append(f"{corpus_id}: {len(gene_symbols - symbols)} of "
                            "the 18,164 genes.csv symbols are absent from var")
        data[corpus_id] = {
            "var": len(symbols), "crosswalk_rows": int(len(rows)),
            "missing_var_symbols": len(missing_symbols),
            "missing_var_hgnc_ids": len(missing_hgnc),
            "genes_symbols_subset": bool(gene_symbols <= symbols),
        }
    # Tahoe: symbol plus both native id namespaces.
    tahoe_var = _read_var(paths.TAHOE_H5AD)
    tahoe_rows = crosswalk[crosswalk["dataset_id"] == "tahoe"]
    tahoe_symbols = set(tahoe_var.index.astype(str))
    missing_symbols = tahoe_symbols - set(tahoe_rows["symbol"].astype(str))
    if missing_symbols:
        problems.append(f"tahoe: {len(missing_symbols)} var symbols without a "
                        f"crosswalk row, e.g. {sorted(missing_symbols)[:5]}")
    native = {
        "token_id": {str(int(value)) for value in tahoe_var["token_id"]},
        "ensembl_id": {str(value) for value in tahoe_var["ensembl_id"]},
    }
    missing_by_namespace: dict[str, int] = {}
    for column, ids in native.items():
        missing_ids = ids - set(tahoe_rows["source_gene_id"].astype(str))
        missing_by_namespace[column] = len(missing_ids)
        if missing_ids:
            problems.append(f"tahoe: {len(missing_ids)} {column} values "
                            f"without a crosswalk row, e.g. "
                            f"{sorted(missing_ids)[:5]}")
    data["tahoe"] = {
        "var": len(tahoe_symbols),
        "crosswalk_rows": int(len(tahoe_rows)),
        "missing_var_symbols": len(missing_symbols),
        "missing_token_id": missing_by_namespace["token_id"],
        "missing_ensembl_id": missing_by_namespace["ensembl_id"],
    }
    # compounds.csv covers every non-control/non-blank row's drug_key.
    for corpus_id in paths.CORPUS_IDS:
        obs = _read_obs(_corpus_path(corpus_id))
        control = obs["is_control"].astype(bool).to_numpy()
        blank = obs["is_blank"].astype(bool).to_numpy()
        used = set(obs.loc[~control & ~blank, "drug_key"].astype(str))
        missing_keys = used - compounds
        if missing_keys:
            problems.append(f"{corpus_id}: compounds.csv misses "
                            f"{len(missing_keys)} perturbed drug_keys, e.g. "
                            f"{sorted(missing_keys)[:5]}")
        data[corpus_id]["perturbed_drug_keys"] = len(used)
        data[corpus_id]["missing_drug_keys"] = len(missing_keys)
    passed = not problems
    detail = ("every var id of all five corpora resolves in "
              "gene_crosswalk.csv; the 18,164 genes.csv symbols are a subset "
              "of every CM corpus var; compounds.csv covers every perturbed "
              "drug_key" if passed else "; ".join(problems[:5]))
    return GateCheck("G3", GATE_TITLES["G3"], passed, detail,
                     {"totals": data, "problems": problems})


# ---------------------------------------------------------------------------
# G4 — v1 reconstruction from the v2 corpora
# ---------------------------------------------------------------------------
class _GeneMapProxy:
    """Identity map: one corpus var column -> one position on the 18,164 axis."""

    def __init__(self, hgnc_ids):
        self.hgnc_ids = np.asarray(hgnc_ids, dtype=object)

    def targets(self, universe: GeneUniverse) -> np.ndarray:
        return np.array(
            [universe.position.get(str(hgnc_id), -1)
             for hgnc_id in self.hgnc_ids], dtype=np.int32)


class _SourceProxy:
    """The ``SourceData`` surface ``views.build_view`` actually reads."""

    __slots__ = ("dataset_id", "meta", "counts", "gene_map", "mapped_counts")

    def __init__(self, dataset_id, meta, counts, gene_map):
        self.dataset_id = dataset_id
        self.meta = meta
        self.counts = counts
        self.gene_map = gene_map
        self.mapped_counts = None


class _CompoundTableProxy:
    """``CompoundTable.members(tier)`` backed by ``tables/compounds.csv``."""

    def __init__(self, frame: pd.DataFrame):
        self.frame = frame

    def members(self, tier: str) -> set[str]:
        column = {"atc-moa": "pharm_atc_moa",
                  "atc-moa-mesh": "pharm_mesh",
                  "all-smiles": "pharm_all_smiles"}[tier]
        mask = self.frame[column].astype(str).str.lower() == "true"
        return set(self.frame.loc[mask, "drug_key"].astype(str))


def _load_frozen_python(path: Path, module_name: str):
    """Execute a frozen file read-only by path in the live ``cmdata`` package.

    The archived ``views.py`` keeps ``from .drugs import ...`` / ``from .genes
    import ...``; giving the exec namespace ``__package__ = "cmdata"`` resolves
    those to the live (unchanged) symbols.  The file is neither imported nor
    copied nor edited.
    """
    module = types.ModuleType(module_name)
    module.__package__ = "cmdata"
    module.__file__ = str(path)
    sys.modules[module_name] = module
    try:
        exec(compile(path.read_text(), str(path), "exec"), module.__dict__)
    finally:
        sys.modules.pop(module_name, None)
    return module


def _load_g4_universe() -> GeneUniverse:
    if not FROZEN_GENES.exists():
        raise VerifyError(f"frozen v1 gene axis missing: "
                          f"{_rel(FROZEN_GENES)}")
    return GeneUniverse(pd.read_csv(FROZEN_GENES))


def _load_g4_source(corpus_id: str, all_smiles: dict) -> _SourceProxy:
    path = _corpus_path(corpus_id)
    if not path.exists():
        raise VerifyError(f"corpus missing: {_rel(path)}")
    with h5py.File(path, "r") as handle:
        obs = read_elem(handle["obs"])
        var = read_elem(handle["var"])
        counts = np.ascontiguousarray(handle["X"][:])
    obs = obs.copy()
    obs.index = obs.index.astype(str)
    obs["well"] = obs["well"].astype(str)
    obs["dataset_id"] = obs["dataset_id"].astype(str).map(
        lambda value: V1_DATASET_ID.get(value, value))
    control = obs["is_control"].astype(bool).to_numpy()
    blank = obs["is_blank"].astype(bool).to_numpy()
    obs["entity_type"] = np.where(
        blank, "blank", np.where(control, "control", "small_molecule"))
    obs["pharm_all_smiles"] = obs["drug_key"].astype(str).map(
        lambda key: all_smiles.get(key, False)).astype(bool)
    for column in ("compound_name", "drug_key", "inchikey", "smiles",
                   "dose_source", "vehicle", "plate", "cell_type", "donor",
                   "atc", "moa", "mesh_classes", "annotation_source"):
        obs[column] = obs[column].astype(str)
    return _SourceProxy(corpus_id, obs, counts,
                        _GeneMapProxy(var["hgnc_id"].astype(str).values))


def _v1_content_digest(adata) -> str:
    """The frozen v1 recipe, copied verbatim from the v1 evidence record."""
    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(adata.X).tobytes())
    digest.update(np.ascontiguousarray(
        adata.layers["log1p_cp10k"]).tobytes())
    digest.update(adata.obs.to_csv(index=False, lineterminator="\n").encode())
    digest.update(adata.var.to_csv(lineterminator="\n").encode())
    digest.update(adata.uns["control_mean"].to_csv(
        lineterminator="\n").encode())
    return digest.hexdigest()


def _first_differing_cell(got: pd.DataFrame, want: pd.DataFrame) -> str:
    if list(got.columns) != list(want.columns):
        return f"columns {list(got.columns)} vs {list(want.columns)}"
    if len(got) != len(want):
        return f"{len(got)} rows vs {len(want)}"
    for column in got.columns:
        left = got[column].astype(str).to_numpy()
        right = want[column].astype(str).to_numpy()
        different = np.flatnonzero(left != right)
        if different.size:
            row = int(different[0])
            return (f"column {column!r} row {row} ({got.index[row]}): "
                    f"rebuilt {left[row]!r} vs frozen {right[row]!r}")
    if list(got.index.astype(str)) != list(want.index.astype(str)):
        return "row index differs"
    return "values equal; CSV serialisation differs"


def _g4_mismatch_detail(name: str, rebuilt, frozen_path: Path) -> str:
    """Locate the first differing component/row/column against the frozen view."""
    import anndata as ad

    frozen = ad.read_h5ad(frozen_path)
    if rebuilt.X.shape != frozen.X.shape:
        return f"{name}: shape {rebuilt.X.shape} vs frozen {frozen.X.shape}"
    if not np.array_equal(rebuilt.X, frozen.X):
        flat = np.flatnonzero(rebuilt.X.ravel() != frozen.X.ravel())[0]
        row, column = divmod(int(flat), rebuilt.X.shape[1])
        return (f"{name}: component X, first difference at row {row} "
                f"({rebuilt.obs.index[row]}), gene column {column} "
                f"({rebuilt.var.index[column]}): "
                f"{int(rebuilt.X[row, column])} vs "
                f"{int(frozen.X[row, column])}")
    if not np.array_equal(rebuilt.layers["log1p_cp10k"],
                          frozen.layers["log1p_cp10k"]):
        flat = np.flatnonzero(
            rebuilt.layers["log1p_cp10k"].ravel()
            != frozen.layers["log1p_cp10k"].ravel())[0]
        row, column = divmod(int(flat), rebuilt.X.shape[1])
        return (f"{name}: component log1p_cp10k, first difference at row {row}, "
                f"gene column {column}")
    for component, got, want in (
            ("obs", rebuilt.obs, frozen.obs),
            ("var", rebuilt.var, frozen.var),
            ("uns[control_mean]", rebuilt.uns["control_mean"],
             frozen.uns["control_mean"])):
        include_index = component != "obs"     # the v1 recipe drops the well index
        got_csv = got.to_csv(index=include_index, lineterminator="\n")
        want_csv = want.to_csv(index=include_index, lineterminator="\n")
        if got_csv != want_csv:
            return f"{name}: component {component}, " \
                   f"{_first_differing_cell(got, want)}"
    return f"{name}: components equal but the content digest differs"


def check_g4_v1_reconstruction() -> GateCheck:
    problems: list[str] = []
    if not FROZEN_MANIFEST.exists():
        raise VerifyError(f"frozen v1 manifest missing: "
                          f"{_rel(FROZEN_MANIFEST)}")
    if not FROZEN_VIEWS_PY.exists() or not FROZEN_PATHS_PY.exists():
        raise VerifyError("frozen v1 writer is missing from "
                          f"{_rel(FROZEN_CM_DIR / 'execution')}")
    manifest = pd.read_csv(FROZEN_MANIFEST)
    expected = dict(zip(manifest["view"].astype(str),
                        manifest["content_sha256"].astype(str)))
    universe = _load_g4_universe()
    views_module = _load_frozen_python(FROZEN_VIEWS_PY,
                                       "cmdata._frozen_v1_views")
    paths_module = _load_frozen_python(FROZEN_PATHS_PY,
                                       "cmdata._frozen_v1_paths")
    compounds = pd.read_csv(paths.TABLES_DIR / "compounds.csv",
                            keep_default_na=False,
                            dtype={"drug_key": str})
    all_smiles = {str(key): bool(str(flag).lower() == "true") for key, flag in
                  zip(compounds["drug_key"], compounds["pharm_all_smiles"])}
    table = _CompoundTableProxy(compounds)
    sources = {}
    view_data: list[dict] = []
    try:
        for corpus_id in ["dtoxs_primary", "dtoxs_ipsc", "sctox", "tenv"]:
            sources[paths.CORPUS_TO_SOURCE[corpus_id]] = _load_g4_source(
                corpus_id, all_smiles)
        for name, datasets, tier in paths_module.VIEWS:
            view = views_module.build_view(name, datasets, tier, sources,
                                           universe, table)
            adata = view.adata
            rebuilt_hash = _v1_content_digest(adata)
            frozen_hash = expected.get(name, "")
            match = rebuilt_hash == frozen_hash
            manifest_row = manifest[manifest["view"].astype(str) == name].iloc[0]
            count_checks = {
                "n_obs": (int(adata.n_obs), int(manifest_row["n_obs"])),
                "n_perturbed": (view.n_perturbed,
                                int(manifest_row["n_perturbed"])),
                "n_controls": (view.n_controls,
                               int(manifest_row["n_controls"])),
                "n_drugs": (view.n_drugs, int(manifest_row["n_drugs"])),
                "n_genes": (int(adata.n_vars), int(manifest_row["n_genes"])),
            }
            bad_counts = {key: pair for key, pair in count_checks.items()
                          if pair[0] != pair[1]}
            if not match:
                frozen_path = FROZEN_VIEW_DIR / name
                if frozen_path.exists():
                    problems.append(_g4_mismatch_detail(name, adata,
                                                        frozen_path))
                else:
                    problems.append(f"{name}: rebuilt {rebuilt_hash} != frozen "
                                    f"{frozen_hash} (frozen view h5ad absent, "
                                    "cannot localise the difference)")
            if bad_counts:
                problems.append(f"{name}: counts {bad_counts} != manifest row")
            view_data.append({
                "view": name,
                "datasets": list(datasets),
                "tier": tier,
                "rows": int(adata.n_obs),
                "n_perturbed": view.n_perturbed,
                "n_controls": view.n_controls,
                "n_drugs": view.n_drugs,
                "n_genes": int(adata.n_vars),
                "rebuilt_content_sha256": rebuilt_hash,
                "frozen_content_sha256": frozen_hash,
                "match": bool(match),
                "counts_match_manifest": not bad_counts,
            })
            print(f"    {name}: {rebuilt_hash[:16]}... "
                  f"{'MATCH' if match else 'MISMATCH'} "
                  f"({adata.n_obs:,} rows)")
            del view, adata
            gc.collect()
    finally:
        for source in sources.values():
            source.counts = None
        del sources
        gc.collect()
    matched = sum(1 for entry in view_data if entry["match"])
    passed = (matched == len(expected) == len(view_data) and not problems)
    if len(view_data) != 5:
        problems.append(f"rebuilt {len(view_data)} views, expected 5")
    detail = (f"{matched}/5 content_sha256 equal the frozen "
              f"{_rel(FROZEN_MANIFEST)} (rebuilt from the v2 corpora with the "
              "archived v1 writer)" if passed else
              f"{matched}/5 hashes match; " + "; ".join(problems[:4]))
    return GateCheck("G4", GATE_TITLES["G4"], passed, detail,
                     {"frozen_manifest": _rel(FROZEN_MANIFEST),
                      "views": view_data, "matched": matched,
                      "expected": len(expected), "problems": problems})


# ---------------------------------------------------------------------------
# G5 — reproducibility
# ---------------------------------------------------------------------------
def _load_sources(index) -> dict:
    entrez = load_dtoxs_entrez_map()
    return {dataset_id: load_source(dataset_id, index, entrez)
            for dataset_id in paths.SOURCE_DATASETS}


def check_g5_reproducibility() -> GateCheck:
    corpus_id = "dtoxs_ipsc"          # the smallest corpus (312 x 26,328)
    if not paths.CORPUS_BUILD_RECORD.exists():
        raise VerifyError(f"build record missing: "
                          f"{_rel(paths.CORPUS_BUILD_RECORD)}")
    record = json.loads(paths.CORPUS_BUILD_RECORD.read_text())[
        "corpora"][corpus_id]
    index, annotation = load_assets()
    loaded = _load_sources(index)
    table = build_compound_table(
        {dataset_id: loaded[dataset_id].meta
         for dataset_id in paths.SOURCE_DATASETS}, annotation)
    scratch = Path(tempfile.mkdtemp(prefix="g5_", dir=paths.EXECUTION_DIR))
    file_hashes, content_hashes = [], []
    try:
        for run in (1, 2):
            adata = build_corpus(
                corpus_id, loaded[paths.CORPUS_TO_SOURCE[corpus_id]], table,
                index)
            content_hashes.append(content_digest(adata))
            info = write_corpus(adata, scratch / f"run{run}_{corpus_id}.h5ad")
            file_hashes.append(info["file_sha256"])
            del adata
            gc.collect()
    finally:
        for source in loaded.values():
            source.counts = None
        shutil.rmtree(scratch, ignore_errors=True)
        gc.collect()
    identical_file = file_hashes[0] == file_hashes[1]
    content_matches_record = content_hashes[0] == record["content_sha256"]
    file_matches_record = file_hashes[0] == record["file_sha256"]
    checks = {
        "identical_file_sha256": identical_file,
        "content_sha256_matches_record": content_matches_record,
        "file_sha256_matches_record": file_matches_record,
    }
    problems = [name for name, ok in checks.items() if not ok]
    passed = not problems
    detail = (f"`{corpus_id}` rebuilt twice in a temp dir: identical file "
              f"sha256 {file_hashes[0][:12]}...; content_sha256 "
              f"{content_hashes[0][:12]}... = corpus_build.json"
              if passed else f"`{corpus_id}`: " + ", ".join(problems))
    return GateCheck("G5", GATE_TITLES["G5"], passed, detail, {
        "corpus": corpus_id,
        "runs": 2,
        "file_sha256_run1": file_hashes[0],
        "file_sha256_run2": file_hashes[1],
        "content_sha256": content_hashes[0],
        "record_content_sha256": record["content_sha256"],
        "record_file_sha256": record["file_sha256"],
        "checks": checks,
        "problems": problems,
    })


# ---------------------------------------------------------------------------
# G6 — path audit
# ---------------------------------------------------------------------------
def check_g6_paths() -> GateCheck:
    script = paths.REPO_ROOT / "script" / "check_stage2_paths.sh"
    if not script.exists():
        raise VerifyError(f"path audit script missing: {_rel(script)}")
    started = time.time()
    process = subprocess.run(["bash", str(script)],
                             cwd=str(paths.REPO_ROOT), capture_output=True,
                             text=True, timeout=600)
    stdout = process.stdout or ""
    stderr = process.stderr or ""
    stale_match = re.search(r"^stale:\s*(\d+) occurrences,\s*(\d+) distinct",
                            stdout, re.MULTILINE)
    files_match = re.search(r"^tracked stage2 executable files scanned:\s*"
                            r"(\d+)", stdout, re.MULTILINE)
    whitelist_match = re.search(r"^whitelisted:\s*(\d+) occurrences",
                                stdout, re.MULTILINE)
    stale = int(stale_match.group(1)) if stale_match else -1
    distinct = int(stale_match.group(2)) if stale_match else -1
    files_scanned = int(files_match.group(1)) if files_match else -1
    whitelisted = int(whitelist_match.group(1)) if whitelist_match else -1
    passed = process.returncode == 0 and stale == 0
    detail = (f"`script/check_stage2_paths.sh`: 0 stale retired step-path "
              f"literal(s) in {files_scanned} tracked stage2 executable "
              f"file(s); {whitelisted} whitelisted (frozen-record payload)"
              if passed else
              f"`script/check_stage2_paths.sh` exit {process.returncode}: "
              f"{stale} stale / {distinct} distinct literal(s)")
    return GateCheck("G6", GATE_TITLES["G6"], passed, detail, {
        "script": _rel(script),
        "exit_code": process.returncode,
        "stale": stale,
        "distinct": distinct,
        "files_scanned": files_scanned,
        "whitelisted": whitelisted,
        "wall_seconds": round(time.time() - started, 2),
        "stdout_tail": stdout.strip().splitlines()[-12:],
        "stderr": stderr.strip()[-2000:],
    })


GATE_FUNCTIONS = [
    ("G1", check_g1_structure),
    ("G2", check_g2_counts),
    ("G3", check_g3_crosswalk),
    ("G4", check_g4_v1_reconstruction),
    ("G5", check_g5_reproducibility),
    ("G6", check_g6_paths),
]


# ---------------------------------------------------------------------------
# Records, README gate section, entry point
# ---------------------------------------------------------------------------
def render_gate_block(payload: dict) -> list[str]:
    """Render the ``## Gate results`` section from a verify payload.

    Shared by ``--verify`` (which patches the README in place) and
    ``--tables`` (which re-renders the whole README), so the two commands do
    not disagree.  No timestamps are emitted: the block is stable while the
    gate outcomes are.
    """
    gates = payload["gates"]
    lines = ["## Gate results", ""]
    lines.append("| gate | result | evidence |")
    lines.append("|---|---|---|")
    for gate in GATE_ORDER:
        entry = gates[gate]
        result = "**PASS**" if entry["passed"] else "**FAIL**"
        lines.append(f"| {GATE_LABELS[gate]} | {result} (`--verify`) | "
                     f"{entry['detail']} |")
    lines.append("")
    g4 = payload["gates"]["G4"]["data"]
    views = g4.get("views", [])
    if views:
        lines.append("### G4 rebuilt v1 view hashes")
        lines.append("")
        lines.append("| view | rebuilt `content_sha256` | frozen "
                     "`manifest.csv` | n_obs | match |")
        lines.append("|---|---|---|---:|---|")
        for view in views:
            lines.append(
                f"| `{view['view']}` | `{view['rebuilt_content_sha256']}` | "
                f"`{view['frozen_content_sha256']}` | {view['rows']:,} | "
                f"{'yes' if view['match'] else '**NO**'} |")
        lines.append("")
    lines.append("Machine-readable record: "
                 "[`execution/build_record/verify.json`]"
                 "(../../../execution/build_record/verify.json) (ticket #17). "
                 "Regenerate with "
                 "`python stage2/02_dataset_corpus/execution/build_corpus.py"
                 " --verify`.")
    lines.append("")
    return lines


def refresh_readme_gates(payload: dict) -> dict:
    """Replace the README ``## Gate results`` section with the verify block."""
    readme = paths.CORPUS_V2_README
    lines = readme.read_text().splitlines()
    try:
        start = lines.index("## Gate results")
    except ValueError as exc:
        raise VerifyError("README has no '## Gate results' section") from exc
    end = next((index for index in range(start + 1, len(lines))
                if lines[index].startswith("## ")), len(lines))
    block = render_gate_block(payload)
    new_lines = lines[:start] + block + lines[end:]
    return tables.write_text("\n".join(new_lines).rstrip("\n") + "\n", readme)


def refresh_coverage_record(payload: dict, readme_info: dict) -> dict:
    """Refresh the README hash and the G4-G6 gate results in the #16 record."""
    record = json.loads(paths.COVERAGE_RECORD.read_text())
    record["readme"] = readme_info
    record["gate_results"].pop("G4-G6", None)
    record["gate_results"].pop("verified_by", None)
    record["gate_results"]["G4"] = _gate_summary(payload["gates"]["G4"])
    record["gate_results"]["G5"] = _gate_summary(payload["gates"]["G5"])
    record["gate_results"]["G6"] = _gate_summary(payload["gates"]["G6"])
    record["gate_results"]["verified_by"] = _rel(paths.VERIFY_RECORD)
    record["verification"] = {
        "record": _rel(paths.VERIFY_RECORD),
        "checked_at": payload["checked_at"],
        "result": payload["result"],
        "gates": {gate: bool(payload["gates"][gate]["passed"])
                  for gate in GATE_ORDER},
    }
    return tables.write_json(record, paths.COVERAGE_RECORD)


def _gate_summary(entry: dict) -> dict:
    return {"passed": bool(entry["passed"]), "detail": entry["detail"]}


def _write_verify_record(payload: dict) -> dict:
    path = paths.VERIFY_RECORD
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(payload, indent=2) + "\n")
    info = {"file": _rel(path), "sha256": _sha256_file(temp),
            "bytes": temp.stat().st_size}
    temp.replace(path)
    return info


def run_verify(command: str | None = None) -> int:
    """Run G1-G6, print one line per gate, write ``verify.json``; 0 iff all pass."""
    started = time.time()
    print("Stage 2 dataset corpus verification — gates G1-G6 (ticket #17)")
    available = _available_gib()
    print(f"available memory: {available:.1f} GiB "
          f"(need >= {MIN_FREE_GIB:.0f} GiB for G5's rebuild)")
    if available < MIN_FREE_GIB:
        print(f"REFUSING: only {available:.1f} GiB available", file=sys.stderr)
        return 2
    checks: list[GateCheck] = []
    for gate, function in GATE_FUNCTIONS:
        clock = time.time()
        try:
            check = function()
        except Exception as exc:                      # noqa: BLE001 - gate boundary
            check = GateCheck(gate, GATE_TITLES[gate], False,
                              f"gate raised {type(exc).__name__}: {exc}",
                              {"traceback": traceback.format_exc()})
        checks.append(check)
        status = "PASS" if check.passed else "FAIL"
        print(f"{gate} {check.title:<20} {status} "
              f"({time.time() - clock:.1f}s) {check.detail}")
    peak_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    passed = sum(1 for check in checks if check.passed)
    payload = {
        "ticket": "#17",
        "adr": "0033",
        "step": "stage2/02_dataset_corpus",
        "checked_at": _iso_now(),
        "command": command or ("python stage2/02_dataset_corpus/execution/"
                               "build_corpus.py --verify"),
        "gates": {check.gate: check.payload() for check in checks},
        "passed": passed,
        "total": len(checks),
        "result": "PASS" if passed == len(checks) else "FAIL",
        "environment": _environment(),
        "wall_seconds": round(time.time() - started, 1),
        "peak_rss_kib": int(peak_kib),
        "peak_rss_gib": round(peak_kib / 1024 ** 2, 2),
    }
    print("Refreshing the coverage README gate section and the #16 record")
    readme_info = refresh_readme_gates(payload)
    payload["readme"] = readme_info
    payload["coverage_record"] = refresh_coverage_record(payload, readme_info)
    payload["verify_record"] = _write_verify_record(payload)
    print(f"verify record: {_rel(paths.VERIFY_RECORD)}")
    print(f"peak RSS {peak_kib / 1024 ** 2:.2f} GiB, "
          f"wall {time.time() - started:.1f}s")
    if payload["result"] == "PASS":
        print(f"RESULT: PASS — {passed}/{len(checks)} gates")
        return 0
    print(f"RESULT: FAIL — {passed}/{len(checks)} gates passed",
          file=sys.stderr)
    return 1
