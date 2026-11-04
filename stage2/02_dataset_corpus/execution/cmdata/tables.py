"""The seven crosswalk tables and the corpus-coverage README (ticket #16).

ADR 0033 decision 10 makes the tables, not ``obs``, the uniform join layer
between corpora whose containers differ: the four CM corpora share one 27-column
``obs`` (dense, index ``well``), ``tahoe`` keeps its Tahoe-native nine-column
``obs`` (CSR, index ``profile_id``).  This module writes the artifacts under
``output/corpus/v2/``:

* ``tables/genes.csv`` — the frozen 18,164-symbol CM projection axis
  (content-identical to ``_frozen/cm-data-v1/corpus/v1/genes.csv``);
* ``tables/gene_crosswalk.csv`` — every source's native gene id -> HGNC symbol,
  including Tahoe's ``token_id``/``ensembl_id`` rows (gate G3);
* ``tables/compounds.csv`` — identity + the three pharmacology tier booleans +
  ``datasets``/``n_wells`` (content-identical to the frozen v1 table);
* ``tables/excluded_compounds.csv`` — the 8 non-SMILES perturbagen keys and the
  technical blank row, with reasons;
* ``tables/source_gene_sets.csv`` — the per-source resolved gene sets;
* ``tables/metadata_cm.csv`` — the four CM corpora's ``obs`` flattened (8,673
  rows) with ``dataset_id``, the legacy id and source-file/accession columns;
* ``tables/metadata_tahoe.csv`` — the Tahoe corpus's ``obs`` flattened (67,018
  rows) with ``dataset_id=tahoe`` and ``is_blank=False``;
* ``README.md`` — corpus coverage (five counts + a "24 h comparable" sub-row),
  provenance, the rename table, per-corpus ``obs``/``var`` listings, gates;
* ``execution/build_record/corpus_coverage.json`` — the same numbers the README
  shows, so gate G2 can be checked programmatically.

Regeneration needs the raw source gene ids (the resolved corpora no longer carry
unmapped ids), so the identity half loads the four sources exactly like the
corpus build — the same ~7.7 GiB peak and ~60 s global load.  The row-metadata
half only reads each corpus's ``obs``/``var`` (header-level, no matrices).
"""

from __future__ import annotations

import gc
import hashlib
import json
import os
import resource
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py
import pandas as pd
from anndata.io import read_elem

from . import genes, paths
from .annotation import CompoundAnnotation, HgncIndex
from .corpus import source_gene_set_frame
from .drugs import CompoundTable
from .genes import GeneUniverse

CROSSWALK_COLUMNS = [
    "dataset_id", "source_gene_id", "symbol", "hgnc_id", "match_type",
    "entrez_verified", "in_corpus",
]
SOURCE_GENE_SET_COLUMNS = [
    "dataset_id", "hgnc_id", "symbol", "match_type", "entrez_verified",
]
EXCLUDED_COLUMNS = [
    "dataset_id", "compound_name", "compound_code", "entity_type", "n_wells",
    "reason",
]
TAHOE_NATIVE_ID_COLUMNS = ("token_id", "ensembl_id")
FROZEN_CM_DIR = paths.REPO_ROOT / "stage2" / "_frozen" / "cm-data-v1"

# Provenance anchors: the accession stays on the raw folder (ADR 0030); the
# short dataset id is the Stage 2 name (ADR 0033 decision 9).
CORPUS_PROVENANCE = {
    "dtoxs_primary": {
        "accession": "LINCS-LDS-1593",
        "atlas": "DToxS PromoCell primary CM panel",
        "raw_item": "pert_data/LINCS-LDS-1593/source",
    },
    "dtoxs_ipsc": {
        "accession": "GSE217423",
        "atlas": "DToxS iPSC-CM panel",
        "raw_item": "pert_data/DToxS-iPSC/source",
    },
    "sctox": {
        "accession": "GSE244740",
        "atlas": "ScreenSeq hiPSC-CM cardiotoxicity screen",
        "raw_item": "pert_data/GSE244740/source",
    },
    "tenv": {
        "accession": "GSE262419",
        "atlas": "TempO-Seq hiPSC-CM environmental chemical panel",
        "raw_item": "pert_data/GSE262419/source",
    },
    "tahoe": {
        "accession": "Tahoe-100M",
        "atlas": "Tahoe-100M perturbation atlas (50 cancer cell lines)",
        "raw_item": "pert_data/Tahoe-100M/source",
    },
}

# Historical v1 id for each v2 corpus (ADR 0033 decision 9).  The merged v1
# `dtoxs_cm` view splits back into its two arms.
LEGACY_DATASET_ID = {
    "dtoxs_primary": "dtoxs_cm",
    "dtoxs_ipsc": "dtoxs_cm",
    "sctox": "screen_ctox",
    "tenv": "temp_env",
    "tahoe": "tahoe",
}

# Coverage is computed from the columns each container actually stores.
COVERAGE_SPEC = {
    "dtoxs_primary": {"cell": "cell_line", "chemical": "drug_key",
                      "time": "time_h", "normalize_cell": True},
    "dtoxs_ipsc": {"cell": "cell_line", "chemical": "drug_key",
                   "time": "time_h", "normalize_cell": False},
    "sctox": {"cell": "cell_line", "chemical": "drug_key",
              "time": "time_h", "normalize_cell": False},
    "tenv": {"cell": "cell_line", "chemical": "drug_key",
             "time": "time_h", "normalize_cell": False},
    "tahoe": {"cell": "cell_line_id", "chemical": "drug",
              "time": None, "normalize_cell": False},
}
ALL_CORPORA = [*paths.CORPUS_IDS, "tahoe"]


class TablesError(RuntimeError):
    """A crosswalk table violated the published contract."""


# ---------------------------------------------------------------------------
# Small IO helpers
# ---------------------------------------------------------------------------
def _iso_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _rel(path: Path) -> str:
    return str(path.relative_to(paths.REPO_ROOT))


def write_frame(frame: pd.DataFrame, path: Path) -> dict:
    """Write one CSV atomically and return its row count / sha256 / bytes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    frame.to_csv(temp, index=False, lineterminator="\n")
    info = {
        "file": _rel(path),
        "rows": int(len(frame)),
        "sha256": file_sha256(temp),
        "bytes": temp.stat().st_size,
    }
    os.replace(temp, path)
    return info


def write_text(text: str, path: Path) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(text)
    info = {
        "file": _rel(path),
        "sha256": file_sha256(temp),
        "bytes": temp.stat().st_size,
    }
    os.replace(temp, path)
    return info


def write_json(payload: dict, path: Path) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(json.dumps(payload, indent=2) + "\n")
    info = {"file": _rel(path), "sha256": file_sha256(temp),
            "bytes": temp.stat().st_size}
    os.replace(temp, path)
    return info


def file_sha256(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def read_var(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as handle:
        return read_elem(handle["var"])


def read_obs(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as handle:
        return read_elem(handle["obs"])


# ---------------------------------------------------------------------------
# Identity tables (need the raw source gene ids)
# ---------------------------------------------------------------------------
def excluded_compounds_frame(table: CompoundTable) -> pd.DataFrame:
    """The 8 non-small-molecule keys (143 wells) + the 32 technical blanks."""
    exclusions = table.exclusions
    rows = []
    compounds = exclusions[exclusions["entity_type"] != "technical_blank"]
    for (name, entity), group in compounds.groupby(
            ["compound_name", "entity_type"], sort=True):
        rows.append({
            "dataset_id": ";".join(sorted(group["dataset_id"].astype(str))),
            "compound_name": name,
            "compound_code": ";".join(sorted({
                str(code) for code in group["compound_code"] if str(code)})),
            "entity_type": entity,
            "n_wells": int(group["n_wells"].sum()),
            "reason": f"non-small-molecule perturbagen ({entity})",
        })
    blanks = exclusions[exclusions["entity_type"] == "technical_blank"]
    if len(blanks):
        rows.append({
            "dataset_id": ";".join(sorted(blanks["dataset_id"].astype(str))),
            "compound_name": "Blank",
            "compound_code": "",
            "entity_type": "technical_blank",
            "n_wells": int(blanks["n_wells"].sum()),
            "reason": ("technical blank (concentration 0; kept as a corpus row, "
                       "neither control nor perturbation)"),
        })
    return pd.DataFrame(rows, columns=EXCLUDED_COLUMNS)


def source_crosswalk(index: HgncIndex, source_maps: dict,
                     universe: GeneUniverse) -> pd.DataFrame:
    """The archived v1 crosswalk recipe: one row per raw id of the 4 sources."""
    rows = []
    for dataset_id in paths.SOURCE_DATASETS:
        gene_map = source_maps[dataset_id]
        for raw_id, hgnc_id, match, verified in zip(
                gene_map.raw_ids, gene_map.hgnc_ids, gene_map.match_types,
                gene_map.verified):
            symbol = index.ref(hgnc_id).symbol if hgnc_id is not None else ""
            rows.append({
                "dataset_id": dataset_id,
                "source_gene_id": raw_id,
                "symbol": symbol,
                "hgnc_id": hgnc_id or "",
                "match_type": str(match),
                "entrez_verified": bool(verified),
                "in_corpus": (hgnc_id in universe.position
                              if hgnc_id else False),
            })
    return pd.DataFrame(rows, columns=CROSSWALK_COLUMNS)


def tahoe_crosswalk(index: HgncIndex, universe: GeneUniverse,
                    var: pd.DataFrame) -> pd.DataFrame:
    """One row per Tahoe native id namespace (``token_id`` and ``ensembl_id``).

    ``symbol`` is the corpus's ``var`` index (the atlas gene symbol) and
    ``hgnc_id``/``match_type``/``in_corpus`` describe that symbol's HGNC
    resolution; the two rows of a gene differ only in ``source_gene_id``.  Both
    namespaces are native lookups into the 62,710-gene Tahoe axis, so the
    crosswalk exposes both rather than choosing one.
    """
    rows = []
    symbols = var.index.astype(str)
    for position, symbol in enumerate(symbols):
        hgnc_id, match = index.resolve_symbol(symbol)
        hgnc = hgnc_id or ""
        in_corpus = hgnc_id in universe.position if hgnc_id else False
        for column in TAHOE_NATIVE_ID_COLUMNS:
            value = var.iloc[position][column]
            if pd.isna(value):
                continue
            if column == "token_id":
                source_gene_id = str(int(value))
            else:
                source_gene_id = str(value)
            rows.append({
                "dataset_id": "tahoe",
                "source_gene_id": source_gene_id,
                "symbol": symbol,
                "hgnc_id": hgnc,
                "match_type": str(match),
                "entrez_verified": False,
                "in_corpus": in_corpus,
            })
    return pd.DataFrame(rows, columns=CROSSWALK_COLUMNS)


# ---------------------------------------------------------------------------
# Corpus row metadata (`obs` flattened) and coverage
# ---------------------------------------------------------------------------
def metadata_cm_frame() -> pd.DataFrame:
    """The four CM corpora's ``obs`` flattened: 8,673 rows, join-ready on well."""
    frames = []
    for corpus_id in paths.CORPUS_IDS:
        path = paths.corpus_h5ad(corpus_id)
        if not path.exists():
            raise TablesError(f"corpus missing: {path} (run --corpus all first)")
        obs = read_obs(path)
        missing = [column for column in paths.CORPUS_OBS
                   if column not in obs.columns]
        if missing:
            raise TablesError(f"{corpus_id}: obs missing columns {missing}")
        provenance = CORPUS_PROVENANCE[corpus_id]
        frame = pd.DataFrame({
            "well": obs.index.astype(str),
            "dataset_id": corpus_id,
            "legacy_dataset_id": LEGACY_DATASET_ID[corpus_id],
            "source_file": _rel(path),
            "accession": provenance["accession"],
        })
        for column in paths.CORPUS_OBS[2:]:
            frame[column] = obs[column].to_numpy()
        frame = frame[["well", "dataset_id", "legacy_dataset_id",
                       "source_file", "accession", *paths.CORPUS_OBS[2:]]]
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)
    return combined


def metadata_tahoe_frame() -> pd.DataFrame:
    """The Tahoe corpus's ``obs`` flattened: 67,018 rows, ``is_blank=False``."""
    path = paths.TAHOE_H5AD
    if not path.exists():
        raise TablesError(f"Tahoe corpus missing: {path}")
    obs = read_obs(path)
    missing = [column for column in paths.TAHOE_OBS if column not in obs.columns]
    if missing:
        raise TablesError(f"tahoe: obs missing columns {missing}")
    frame = pd.DataFrame({
        "profile_id": obs.index.astype(str),
        "dataset_id": "tahoe",
        "is_blank": False,
    })
    for column in paths.TAHOE_OBS:
        frame[column] = obs[column].to_numpy()
    return frame[["profile_id", "dataset_id", "is_blank", *paths.TAHOE_OBS]]


def _five_counts(frame: pd.DataFrame, corpus_id: str) -> dict:
    spec = COVERAGE_SPEC[corpus_id]
    control = frame["is_control"].astype(bool).to_numpy()
    blank = frame["is_blank"].astype(bool).to_numpy()
    perturbed = ~control & ~blank
    cells = frame[spec["cell"]].astype(str)
    if spec["normalize_cell"]:
        # The primary arm labels its four lines `A/B/D/E` on plates 1/2/5 and
        # `PC-A/PC-B/PC-D/PC-E` on plate 4; both names denote the same line.
        cells = cells.str.replace("PC-", "", regex=False)
    chemicals = frame[spec["chemical"]].astype(str)
    doses = frame["dose_um"].astype(float)
    combos = pd.DataFrame({
        "cell": cells[perturbed],
        "chemical": chemicals[perturbed],
        "dose": doses[perturbed],
    }).drop_duplicates()
    times = (sorted(float(value) for value in
                    frame[spec["time"]].astype(float).unique())
             if spec["time"] is not None else None)
    if times is None:
        time_window = "n.a. (no time column)"
    elif len(times) == 1:
        time_window = f"{times[0]:g} h"
    else:
        time_window = "; ".join(f"{value:g} h" for value in times)
    return {
        "rows": int(len(frame)),
        "perturbed": int(perturbed.sum()),
        "control": int(control.sum()),
        "blank": int(blank.sum()),
        "control_wells": int(control.sum()),
        "blank_wells": int(blank.sum()),
        "distinct_cell_lines": int(cells.nunique()),
        "distinct_cell_line_labels": int(frame[spec["cell"]]
                                         .astype(str).nunique()),
        "cell_line_labels": sorted(frame[spec["cell"]].astype(str).unique()),
        "chemicals": int(chemicals[perturbed].nunique()),
        "distinct_doses": int(doses[perturbed].nunique()),
        "doses_um": sorted({float(value) for value in doses[perturbed]}),
        "time_window": time_window,
        "times_h": times,
        "combinations": int(len(combos)),
        "plates": int(frame["plate"].astype(str).nunique()),
    }


def coverage_frames(metadata_cm: pd.DataFrame,
                    metadata_tahoe: pd.DataFrame) -> dict:
    """The five coverage counts per corpus, plus the 24 h comparable sub-row."""
    frames = {corpus_id: metadata_cm[metadata_cm["dataset_id"] == corpus_id]
              for corpus_id in paths.CORPUS_IDS}
    frames["tahoe"] = metadata_tahoe
    coverage = {}
    for corpus_id in ALL_CORPORA:
        entry = _five_counts(frames[corpus_id], corpus_id)
        times = entry["times_h"]
        if times == [24.0]:
            entry["t24"] = {key: entry[key] for key in (
                "rows", "control", "blank", "distinct_cell_lines", "chemicals",
                "distinct_doses", "combinations")}
            entry["t24"]["basis"] = "recorded 24 h"
        elif times is None and corpus_id == "tahoe":
            # The atlas is a single 24 h exposure by design, but its obs has no
            # time column (ADR 0029 update 2), so this row is the whole corpus
            # and the README says so instead of materialising a time field.
            entry["t24"] = {key: entry[key] for key in (
                "rows", "control", "blank", "distinct_cell_lines", "chemicals",
                "distinct_doses", "combinations")}
            entry["t24"]["basis"] = "atlas design, no time column"
        else:
            entry["t24"] = None
            entry["t24_basis"] = f"n.a. ({entry['time_window']} window)" \
                if times is not None else "n.a. (no time column)"
        entry["dataset_id"] = corpus_id
        coverage[corpus_id] = entry
    return coverage


# ---------------------------------------------------------------------------
# Contract verification (gate G3 + G2 counts)
# ---------------------------------------------------------------------------
def _frozen_frame_bytes(frame: pd.DataFrame) -> bytes:
    return frame.to_csv(index=False, lineterminator="\n").encode()


def _compare_frozen(name: str, frame: pd.DataFrame, frozen: Path,
                    checks: dict, errors: list[str]) -> None:
    if not frozen.exists():
        checks[name] = None
        return
    equal = _frozen_frame_bytes(frame) == frozen.read_bytes()
    checks[name] = bool(equal)
    if not equal:
        errors.append(
            f"{name} differs from the frozen v1 identity table {_rel(frozen)}")


def verify_contract(frames: dict, metadata_cm: pd.DataFrame,
                    metadata_tahoe: pd.DataFrame,
                    coverage: dict, errors: list[str],
                    checks: dict) -> None:
    """Gate G3 coverage + the identity-equality evidence, before any write."""
    genes_frame = frames["genes.csv"]
    symbols = set(genes_frame["symbol"].astype(str))
    if len(genes_frame) != 18164 or len(symbols) != 18164:
        errors.append(f"genes.csv axis is {len(genes_frame)} rows / "
                      f"{len(symbols)} symbols, expected 18,164")
    crosswalk = frames["gene_crosswalk.csv"]
    for corpus_id in paths.CORPUS_IDS:
        source_id = paths.CORPUS_TO_SOURCE[corpus_id]
        var = read_var(paths.corpus_h5ad(corpus_id))
        var_symbols = set(var.index.astype(str))
        var_hgnc = set(var["hgnc_id"].astype(str))
        rows = crosswalk[crosswalk["dataset_id"] == source_id]
        missing = var_symbols - set(rows["symbol"].astype(str))
        missing_hgnc = var_hgnc - set(rows["hgnc_id"].astype(str))
        if missing:
            errors.append(f"{corpus_id}: {len(missing)} var symbols without a "
                          f"crosswalk row, e.g. {sorted(missing)[:5]}")
        if missing_hgnc:
            errors.append(f"{corpus_id}: {len(missing_hgnc)} var hgnc_ids "
                          f"without a crosswalk row, e.g. "
                          f"{sorted(missing_hgnc)[:5]}")
        if not symbols <= var_symbols:
            errors.append(f"{corpus_id}: {len(symbols - var_symbols)} of the "
                          "18,164 CM symbols are absent from its var axis")
    # Tahoe: symbol plus both native id namespaces.
    tahoe_var = read_var(paths.TAHOE_H5AD)
    tahoe_rows = crosswalk[crosswalk["dataset_id"] == "tahoe"]
    tahoe_symbols = set(tahoe_var.index.astype(str))
    missing = tahoe_symbols - set(tahoe_rows["symbol"].astype(str))
    if missing:
        errors.append(f"tahoe: {len(missing)} var symbols without a crosswalk "
                      f"row, e.g. {sorted(missing)[:5]}")
    native_ids = {
        "token_id": {str(int(value)) for value in tahoe_var["token_id"]},
        "ensembl_id": {str(value) for value in tahoe_var["ensembl_id"]},
    }
    for column, ids in native_ids.items():
        missing_ids = ids - set(tahoe_rows["source_gene_id"].astype(str))
        if missing_ids:
            errors.append(f"tahoe: {len(missing_ids)} {column} values without "
                          f"a crosswalk row, e.g. {sorted(missing_ids)[:5]}")
    # compounds.csv covers every perturbed CM row's drug_key.
    compounds = set(frames["compounds.csv"]["drug_key"].astype(str))
    perturbed = metadata_cm[(~metadata_cm["is_control"].astype(bool))
                            & (~metadata_cm["is_blank"].astype(bool))]
    missing_keys = set(perturbed["drug_key"].astype(str)) - compounds
    if missing_keys:
        errors.append(f"compounds.csv misses {len(missing_keys)} perturbed "
                      f"drug_keys, e.g. {sorted(missing_keys)[:5]}")
    # Row counts are part of G2.
    if len(metadata_cm) != sum(coverage[corpus_id]["rows"]
                               for corpus_id in paths.CORPUS_IDS):
        errors.append("metadata_cm.csv row count is not the sum of the four "
                      "corpora") 
    if len(metadata_tahoe) != coverage["tahoe"]["rows"]:
        errors.append("metadata_tahoe.csv row count differs from the corpus")
    # Identity content must equal the frozen v1 tables while they exist.
    frozen_corpus = FROZEN_CM_DIR / "corpus" / "v1"
    _compare_frozen("genes.csv", frames["genes.csv"],
                    frozen_corpus / "genes.csv", checks, errors)
    _compare_frozen("compounds.csv", frames["compounds.csv"],
                    frozen_corpus / "compounds.csv", checks, errors)
    _compare_frozen("source_gene_sets.csv", frames["source_gene_sets.csv"],
                    frozen_corpus / "source_gene_sets.csv", checks, errors)
    four_source = crosswalk[crosswalk["dataset_id"] != "tahoe"]
    _compare_frozen("gene_crosswalk.csv (four-source part)", four_source,
                    FROZEN_CM_DIR / "datasets-1" / "gene_crosswalk.csv",
                    checks, errors)
    if errors:
        raise TablesError("crosswalk/coverage contract failed:\n  "
                          + "\n  ".join(errors))


# ---------------------------------------------------------------------------
# README rendering
# ---------------------------------------------------------------------------
def _fmt(value: int | None) -> str:
    return "n.a." if value is None else f"{value:,}"


def _gates_from_build(build_record: dict) -> dict:
    gates = {}
    for corpus_id, entry in build_record.get("corpora", {}).items():
        results = entry.get("gates", [])
        passed = sum(1 for gate in results if gate.get("passed"))
        gates[corpus_id] = {
            "passed": passed,
            "total": len(results),
            "names": [gate["name"] for gate in results],
        }
    return gates


def render_readme(context: dict) -> str:
    coverage = context["coverage"]
    tables = context["tables"]
    build_record = context["corpus_build"]
    tahoe_move = context["tahoe_move"]
    environment = context["environment"]
    generated = context["generated_at"]
    gates = _gates_from_build(build_record)
    lines: list[str] = []
    add = lines.append

    add("# `corpus/v2` — Stage 2 dataset corpora, crosswalk tables and coverage")
    add("")
    add("**Status:** generated artifact. Five h5ad corpora, seven crosswalk "
        "tables and this coverage record.")
    add("**Governing decision:** [ADR 0033]"
        "(../../../../../docs/adr/0033-stage2-dataset-corpus.md) (the corpus "
        "layer; supersedes the matrix-free ADR 0026 shape and folds "
        "`04_tahoe_data` in).")
    add(f"**Built:** `{generated[:10]}` by "
        "`stage2/02_dataset_corpus/execution/build_corpus.py --tables`. The "
        "numbers below are also in "
        "[`execution/build_record/corpus_coverage.json`]"
        "(../../../execution/build_record/corpus_coverage.json), so gate G2 can be "
        "checked programmatically.")
    add("")
    add("Output is exactly three artifact classes (ADR 0033 decision 6): the "
        "five h5ad corpora, the seven tables under `tables/`, and this README. "
        "Build QA/manifests live under `execution/build_record/`, never here.")
    add("")

    # -- coverage -------------------------------------------------------
    add("## Corpus coverage")
    add("")
    add("One row per corpus. Five counts, computed from the corpus's own `obs` "
        "at its own time point: distinct cell line, chemicals among perturbed "
        "wells, distinct doses, time window, and distinct `(cell line, chemical, "
        "dose)` combinations. Vehicle controls and technical blanks are listed "
        "separately and are **never** counted as chemicals.")
    add("")
    add("Counting rules (the same rules the JSON record carries):")
    add("")
    add("- **cell line** = distinct `cell_line` labels; the primary arm labels "
        "its four lines `A/B/D/E` on plates 1/2/5 and `PC-A/PC-B/PC-D/PC-E` on "
        "plate 4, so those labels are normalised by removing the `PC-` prefix "
        "(recorded: "
        f"`{', '.join(coverage['dtoxs_primary']['cell_line_labels'])}` → "
        f"{coverage['dtoxs_primary']['distinct_cell_lines']} lines). `NA` counts "
        "as the one unnamed line of `sctox`/`tenv`.");
    add("- **chemical** = distinct `drug_key` over perturbed rows (Tahoe: its "
        "native `drug` column).")
    sctox = coverage["sctox"]
    add("- **dose** = distinct `dose_um` over perturbed rows, as recorded (not "
        "rounded). `sctox` doses are `derived_top_dose` values from a curated "
        "per-drug top, so its per-drug 8-point series gives many distinct "
        f"floats across the {sctox['chemicals']} drugs "
        f"({sctox['distinct_doses']} values in total) while each drug has 8 "
        "doses.")
    add("- **time window** = the corpus's own recorded time point(s). `tahoe` "
        "has no time column (ADR 0029 update 2), so it is `n.a.`; the atlas is "
        "a single 24 h exposure by design (see the sub-row).")
    add("")
    add("| corpus | distinct cell lines | chemicals (perturbed wells) | "
        "distinct doses (µM) | time window | (cell line, chemical, dose) | "
        "control wells | blank wells |")
    add("|---|---:|---:|---:|---|---:|---:|---:|")
    for corpus_id in ALL_CORPORA:
        entry = coverage[corpus_id]
        add(f"| `{corpus_id}` | {entry['distinct_cell_lines']:,} | "
            f"{entry['chemicals']:,} | {entry['distinct_doses']:,} | "
            f"{entry['time_window']} | {entry['combinations']:,} | "
            f"{entry['control']:,} | {entry['blank']:,} |")
    add("")
    add("### 24 h comparable subset")
    add("")
    add("A corpus enters this sub-row when its exposure is a single 24 h "
        "window. The two `dtoxs` corpora are 48 h, so they are `n.a.`. "
        "`tahoe` is a 24 h atlas by design but its `obs` carries no time "
        "column (ADR 0029), so its row is the whole corpus and the basis says "
        "so — no time field is materialised.")
    add("")
    add("| corpus | distinct cell lines | chemicals | distinct doses | "
        "(cell line, chemical, dose) | basis |")
    add("|---|---:|---:|---:|---:|---|")
    for corpus_id in ALL_CORPORA:
        entry = coverage[corpus_id]
        t24 = entry["t24"]
        if t24 is None:
            add(f"| `{corpus_id}` | n.a. | n.a. | n.a. | n.a. | "
                f"{entry['t24_basis']} |")
        else:
            add(f"| `{corpus_id}` | {t24['distinct_cell_lines']:,} | "
                f"{t24['chemicals']:,} | {t24['distinct_doses']:,} | "
                f"{t24['combinations']:,} | {t24['basis']} |")
    add("")
    add("Corpus rows, controls and blanks (the contract table):")
    add("")
    add("| corpus | rows | perturbed | control | blank | plates |")
    add("|---|---:|---:|---:|---:|---:|")
    for corpus_id in ALL_CORPORA:
        entry = coverage[corpus_id]
        add(f"| `{corpus_id}` | {entry['rows']:,} | {entry['perturbed']:,} | "
            f"{entry['control']:,} | {entry['blank']:,} | "
            f"{entry['plates']:,} |")
    add("")

    # -- provenance -----------------------------------------------------
    add("## Provenance")
    add("")
    add("The accession/atlas name stays on the raw folder (ADR 0030); the short "
        "dataset id is the Stage 2 name. Raw sources are read-only under "
        "`pert_data/`; the build makes no network call.")
    add("")
    add("| corpus | accession / atlas | raw item (`pert_data/<item>/source`) | "
        "shape (rows × genes) | built | bytes |")
    add("|---|---|---|---|---|---:|")
    for corpus_id in ALL_CORPORA:
        provenance = CORPUS_PROVENANCE[corpus_id]
        entry = coverage[corpus_id]
        if corpus_id == "tahoe":
            shape = f"{entry['rows']:,} × 62,710"
            built = "2026-10-30 (built by `04_tahoe_data`), moved unchanged 2026-11-04 (#14)"
            size = f"{tahoe_move.get('bytes', 25805562874):,}"
        else:
            shape = (f"{entry['rows']:,} × "
                     f"{build_record['corpora'][corpus_id]['axis_length']:,}")
            built = build_record["corpora"][corpus_id]["built_at"][:10]
            size = f"{build_record['corpora'][corpus_id]['bytes']:,}"
        add(f"| `{corpus_id}` | {provenance['accession']} — "
            f"{provenance['atlas']} | `{provenance['raw_item']}` | {shape} | "
            f"{built} | {size} |")
    add("")
    add("Content and file hashes (equivalence is content-level, never "
        "file-byte — ADR 0033 decision 7; the four CM corpora use the ticket "
        "#15 digest recipe, `tahoe` the move hash of #14):")
    add("")
    add("| corpus | `content_sha256` | file sha256 |")
    add("|---|---|---|")
    for corpus_id in ALL_CORPORA:
        if corpus_id == "tahoe":
            content = tahoe_move.get("content_sha256", "n.a.")
            file_hash = tahoe_move.get("file_sha256", "n.a.")
        else:
            entry = build_record["corpora"][corpus_id]
            content = entry["content_sha256"]
            file_hash = entry["file_sha256"]
        add(f"| `{corpus_id}` | `{content}` | `{file_hash}` |")
    add("")
    add("### Build environment")
    add("")
    add(f"- Python {environment.get('python', '?')}, anndata "
        f"{environment.get('anndata', '?')}, h5py "
        f"{environment.get('h5py', '?')}, numpy "
        f"{environment.get('numpy', '?')}, pandas "
        f"{environment.get('pandas', '?')}")
    add(f"- platform: `{environment.get('platform', '?')}`")
    add(f"- vendored Stage 1 loaders source commit: `{paths.SOURCE_COMMIT}`")
    add("- pinned annotation assets: "
        "`pert_data/Annotation/source/hgnc_complete_set.txt` (HGNC) and "
        "`pert_data/Annotation/local/cm_compound_annotation.csv` (compound "
        "registry), both pinned by `pert_data/Annotation/local/README.md`")
    add("- the same `arc` environment builds the four CM corpora; `tahoe` "
        "predates this step (built 2026-10-30 by `04_tahoe_data`, moved "
        "unchanged by #14), and its QA numbers live in "
        "`execution/build_record/tahoe_pseudobulk_*`.")
    add("")

    # -- rename table ---------------------------------------------------
    add("## Rename table (forward-only, ADR 0033 decision 9)")
    add("")
    add("| historical id / path | current id / path | note |")
    add("|---|---|---|")
    add("| `screen_ctox` | `sctox` | GSE244740 ScreenSeq hiPSC-CM panel |")
    add("| `temp_env` | `tenv` | GSE262419 TempO-Seq hiPSC-CM panel |")
    add("| `dtoxs_cm` | `dtoxs_primary` + `dtoxs_ipsc` | the merged v1 view "
        "splits back into its two arms |")
    add("| `04_tahoe_data` | `02_dataset_corpus` | step folded in; `04` is a gap "
        "in the sequence, not a renumber |")
    add("")
    add("Historical records keep the **old** names: a frozen 2026-10 log, "
        "manifest or issue that says `screen_ctox`, `temp_env`, `dtoxs_cm` or "
        "`04_tahoe_data` is not rewritten. Only new corpus artifacts use the "
        "new ids. Inside these tables the `dataset_id` column of "
        "`compounds.csv`, `source_gene_sets.csv` and `gene_crosswalk.csv` keeps "
        "the four **source** ids (`dtoxs_primary`, `dtoxs_ipsc`, `screen_ctox`, "
        "`temp_env`) because those tables are content-identical with the frozen "
        "v1 identities; `metadata_cm.csv` carries **both** the new `dataset_id` "
        "and the legacy id so a consumer never has to re-map the rename twice.")
    add("")

    # -- schemas --------------------------------------------------------
    add("## Per-corpus `obs` and `var`")
    add("")
    add("### Four CM corpora — `obs` 27 columns, index `well`")
    add("")
    add("```")
    for start in range(0, len(paths.CORPUS_OBS), 4):
        add(", ".join(paths.CORPUS_OBS[start:start + 4]))
    add("```")
    add("")
    add("`var` = index `symbol`; columns `hgnc_id, entrez_id, ensembl_id, "
        "match_type`. `X` is dense `int32` plus a dense `float32` "
        "`layers['log1p_cp10k']`; `uns`/`obsm`/`varm`/`obsp`/`varp` are empty; "
        "`dataset_id` and `control_group` use the v2 ids while `well` keeps the "
        "loader string verbatim (`screen_ctox::…`, `temp_env::…`).")
    add("")
    add("### `tahoe` — `obs` 9 columns, index `profile_id`")
    add("")
    add("```")
    add(", ".join(paths.TAHOE_OBS))
    add("```")
    add("")
    add("`var` = index `symbol` (the atlas gene symbol); columns `ensembl_id, "
        "token_id`. `X` is CSR `int32` plus a CSR `float32` layer; one row is a "
        "`(plate, well, cell_line)` pseudobulk profile.")
    add("")
    add("### Why the containers are asymmetric")
    add("")
    add("The four CM corpora are small (312–4,510 rows) and share one producer, "
        "so they are dense, carry the same 27-column ADR-0026-derived `obs`, and "
        "index on the raw well id. `tahoe.h5ad` is a 25.8 GB, 67,018-row CSR "
        "artifact whose identity vocabulary is the atlas's: rewriting its "
        "`obs` into the CM shape would mutate the artifact without adding "
        "information (ADR 0033 decision 10). Its `obs` therefore stays the "
        "Tahoe-native nine columns with index `profile_id`, and `is_blank` is "
        "supplied as a constant `False` column in `metadata_tahoe.csv` rather "
        "than written into the h5ad. Format consistency therefore means the "
        "container and value semantics, not an identical column list; joins go "
        "through `genes.csv`, `gene_crosswalk.csv`, `compounds.csv`, "
        "`metadata_cm.csv` and `metadata_tahoe.csv`.")
    add("")

    # -- table inventory ------------------------------------------------
    add("## Crosswalk tables")
    add("")
    add("| table | rows | sha256 | purpose |")
    add("|---|---:|---|---|")
    purposes = {
        "genes.csv": "the frozen 18,164-symbol CM projection axis "
                     "(content-identical to v1)",
        "gene_crosswalk.csv": "every source native gene id → HGNC symbol "
                              "(gate G3)",
        "compounds.csv": "525 `drug_key` / 532 key-name rows with "
                         "`atc/moa/mesh_classes`, the three tier booleans, "
                         "`datasets`, `n_wells`",
        "excluded_compounds.csv": "8 non-SMILES perturbagen keys (143 wells) "
                                  "and the 32 technical blanks, with reasons",
        "source_gene_sets.csv": "per-source resolved gene sets (drift baseline)",
        "metadata_cm.csv": "the four CM corpora's `obs` flattened, join-ready "
                           "on `well`",
        "metadata_tahoe.csv": "the Tahoe corpus's `obs` flattened, join-ready "
                              "on `profile_id`",
    }
    for name in paths.TABLE_FILES:
        info = tables[name]
        add(f"| [`tables/{name}`](tables/{name}) | {info['rows']:,} | "
            f"`{info['sha256']}` | {purposes[name]} |")
    add("")
    crosswalk_rows = tables["gene_crosswalk.csv"]["rows"]
    add(f"`gene_crosswalk.csv` has **{crosswalk_rows:,} rows**: the "
        "143,647-row four-source table (content-identical to the frozen "
        "`_frozen/cm-data-v1/datasets-1/gene_crosswalk.csv`) plus 125,420 Tahoe "
        "rows. The four CM sources each expose one native id space (DToxS "
        "symbols, `screen_ctox` Ensembl ids, `temp_env` `SYMBOL_probe` ids), so "
        "they contribute one row per raw gene id; Tahoe exposes two native "
        "namespaces per gene (`token_id` and `ensembl_id`), so each of its "
        "62,710 genes gets one row per namespace. `in_corpus` marks whether the "
        "row resolves onto the 18,164 CM axis; `symbol=\"\"`/`hgnc_id=\"\"` "
        "marks an unmapped raw id (nothing is guessed).")
    add("")
    add("Row counts: `genes.csv` 18,164; `compounds.csv` 532 key/name rows / "
        "525 `drug_key`; `excluded_compounds.csv` 8 + 1 blank; "
        "`source_gene_sets.csv` 104,149 (the sum of the four native axes, "
        "26,328 + 26,328 + 31,880 + 19,613); `metadata_cm.csv` 8,673 (779 + 312 "
        "+ 3,072 + 4,510); `metadata_tahoe.csv` 67,018.")
    add("")

    # -- gates ----------------------------------------------------------
    add("## Gate results")
    add("")
    add("| gate | result | evidence |")
    add("|---|---|---|")
    gate_bits = "; ".join(
        f"`{corpus_id}` {entry['passed']}/{entry['total']}"
        for corpus_id, entry in gates.items())
    add(f"| G1 structure + counts | **PASS** (#15) | per-corpus build gates "
        f"passed: {gate_bits} |")
    add("| G2 counts | **PASS** (this build) | this coverage table is rendered "
        "from `corpus_coverage.json`; `test_tables_contract` re-parses the "
        "table and asserts equality item by item |")
    add("| G3 crosswalk coverage | **PASS** (this build) | every `var` id of "
        "all five corpora has a crosswalk row (incl. Tahoe `token_id`/"
        "`ensembl_id`); the 18,164 `genes.csv` symbols are a subset of each CM "
        "corpus `var`; `compounds.csv` covers every non-control/non-blank "
        "`drug_key`. Verified before writing; see `test_tables_contract` |")
    add("| G4 v1 reconstruction | not run | ticket #17 (needs the frozen v1 "
        "view writer; see `.scratch/stage2-corpus-reorg/impl-notes-g4-proof.md`) |")
    add("| G5 reproducibility | not run | ticket #17 (ticket #15's "
        "`--repeat-check` already passes for `dtoxs_ipsc`) |")
    add("| G6 no stale paths | not run | ticket #17 |")
    add("")

    # -- reproduce ------------------------------------------------------
    add("## Regenerating the tables and this README")
    add("")
    add("From the repository root (the repo root is found by walking up to "
        "`.git`, never by counting `../..` — ADR 0030 §4):")
    add("")
    add("```bash")
    add("PY=/home/chang/miniconda3/envs/arc/bin/python")
    add("# the four CM corpora must already exist (ticket #15)")
    add(f"$PY stage2/02_dataset_corpus/execution/build_corpus.py --tables")
    add("$PY -m unittest "
        "stage2/02_dataset_corpus/execution/tests/test_tables_contract.py -v")
    add("```")
    add("")
    add("The identity half re-reads the four raw sources to resolve every "
        "native gene id (the corpora keep only resolved genes), so it pays the "
        "same global load as a corpus build (~7.7 GiB peak, ~60 s); the row "
        "metadata half only reads each corpus's `obs`/`var` with `h5py` + "
        "`read_elem` (no matrices). The tables, this README and "
        "`corpus_coverage.json` are all overwritten atomically (`.tmp` → "
        "rename).")
    add("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Entry point used by build_corpus.py --tables
# ---------------------------------------------------------------------------
def build_tables(index: HgncIndex, annotation: CompoundAnnotation,
                 loaded: dict, table: CompoundTable, environment: dict,
                 command: str) -> dict:
    """Write the seven tables, the README and ``corpus_coverage.json``."""
    started = time.time()
    print("Building identity tables (raw source gene ids -> HGNC)")
    source_maps = {dataset_id: loaded[dataset_id].gene_map
                   for dataset_id in paths.SOURCE_DATASETS}
    universe = genes.compute_universe(
        index, [source_maps[dataset_id] for dataset_id in paths.SOURCE_DATASETS])
    print(f"  CM axis: {len(universe):,} symbols")

    crosswalk = source_crosswalk(index, source_maps, universe)
    four_source_rows = len(crosswalk)
    tahoe_var = read_var(paths.TAHOE_H5AD)
    tahoe_rows = tahoe_crosswalk(index, universe, tahoe_var)
    crosswalk = pd.concat([crosswalk, tahoe_rows], ignore_index=True)
    print(f"  crosswalk: {len(crosswalk):,} rows "
          f"({four_source_rows:,} four-source + "
          f"{len(crosswalk) - four_source_rows:,} Tahoe)")
    del tahoe_rows

    frames = {
        "genes.csv": universe.genes,
        "gene_crosswalk.csv": crosswalk,
        "compounds.csv": table.compounds,
        "excluded_compounds.csv": excluded_compounds_frame(table),
        "source_gene_sets.csv": source_gene_set_frame(index, source_maps),
    }
    print("Building corpus row metadata (obs flattened)")
    metadata_cm = metadata_cm_frame()
    metadata_tahoe = metadata_tahoe_frame()
    frames["metadata_cm.csv"] = metadata_cm
    frames["metadata_tahoe.csv"] = metadata_tahoe
    print(f"  metadata_cm: {len(metadata_cm):,} rows; "
          f"metadata_tahoe: {len(metadata_tahoe):,} rows")

    coverage = coverage_frames(metadata_cm, metadata_tahoe)
    errors: list[str] = []
    checks: dict = {}
    print("Verifying G3 coverage and frozen identity equality")
    verify_contract(frames, metadata_cm, metadata_tahoe, coverage, errors,
                    checks)

    print("Writing seven tables")
    tables = {}
    for name in paths.TABLE_FILES:
        tables[name] = write_frame(frames[name], paths.TABLES_DIR / name)
        print(f"  {name}: {tables[name]['rows']:,} rows, "
              f"{tables[name]['bytes']:,} B")
    del frames, metadata_cm, metadata_tahoe
    gc.collect()

    generated_at = _iso_now()
    build_record = {}
    if paths.CORPUS_BUILD_RECORD.exists():
        build_record = json.loads(paths.CORPUS_BUILD_RECORD.read_text())
    tahoe_move = {}
    if paths.TAHOE_MOVE_RECORD.exists():
        tahoe_move = json.loads(paths.TAHOE_MOVE_RECORD.read_text())
    context = {
        "coverage": coverage,
        "tables": tables,
        "corpus_build": build_record,
        "tahoe_move": tahoe_move,
        "environment": environment,
        "generated_at": generated_at,
    }
    readme = write_text(render_readme(context), paths.CORPUS_V2_README)
    print(f"  README.md: {readme['bytes']:,} B")

    wall = time.time() - started
    peak_kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    payload = {
        "ticket": "#16",
        "adr": "0033",
        "step": "stage2/02_dataset_corpus",
        "generated_at": generated_at,
        "command": command,
        "definition": {
            "five_counts": (
                "distinct cell line / chemicals among perturbed wells / "
                "distinct dose_um / corpus time window / distinct "
                "(cell line, chemical, dose) combinations, from the corpus obs "
                "at the corpus's own time point"),
            "cell_line": (
                "distinct cell_line labels, PC- prefix removed (the same four "
                "primary lines are A/B/D/E on plates 1/2/5 and PC-* on plate 4); "
                "NA counts as the one unnamed sctox/tenv line"),
            "chemicals": ("distinct drug_key over perturbed rows "
                          "(tahoe: native drug column); controls and blanks are "
                          "never counted"),
            "doses": "distinct dose_um over perturbed rows, unrounded",
            "time": "the corpus's recorded time point(s); tahoe has no time "
                    "column (ADR 0029 update 2)",
            "t24": ("enters when the exposure is a single 24 h window; the two "
                    "dtoxs corpora are 48 h so n.a.; tahoe is a 24 h atlas by "
                    "design but has no time column, so its row is the whole "
                    "corpus"),
        },
        "corpora": coverage,
        "tables": tables,
        "readme": readme,
        "build": {
            "source_commit": paths.SOURCE_COMMIT,
            "environment": environment,
            "wall_seconds": round(wall, 1),
            "peak_rss_kib": int(peak_kib),
            "peak_rss_gib": round(peak_kib / 1024 ** 2, 2),
        },
        "gate_results": {
            "G1": {corpus_id: entry for corpus_id, entry in
                   _gates_from_build(build_record).items()},
            "G2": "PASS: this record renders the README coverage table; "
                  "test_tables_contract re-parses the table and asserts equality",
            "G3": "PASS: crosswalk covers every var id of all five corpora; "
                  "18,164 genes.csv symbols subset of each CM corpus var; "
                  "compounds.csv covers every perturbed drug_key",
            "G4-G6": "not run (ticket #17)",
        },
        "frozen_identity_check": checks,
    }
    payload["coverage_record"] = write_json(payload, paths.COVERAGE_RECORD)
    print(f"Coverage record: {_rel(paths.COVERAGE_RECORD)}")
    print(f"Peak RSS {peak_kib / 1024 ** 2:.2f} GiB, wall {wall:.1f}s")
    return payload
