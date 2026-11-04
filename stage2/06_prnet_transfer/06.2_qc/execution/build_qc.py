#!/usr/bin/env python
"""Per-observation experimental-design + QC audit for the two PRnet-transfer sources.

Pipeline (three separate commands, deliberately):

    python build_qc.py --lock     # derive + freeze thresholds into output/qc_rules.json
    python build_qc.py --build    # apply the frozen rules, write the keep/exclude tables
    python build_qc.py --verify   # re-hash every artifact and re-check the headline counts

``--lock`` runs first and its timestamp is what "locked before any model result" means:
the rule payload (thresholds + evidence + instantiated derived values) is hashed into
``qc_rules.json``, and ``--build`` refuses to run if the code's payload hash no longer
matches the locked one. No model, DEG table or response statistic is read anywhere in
this file.

Inputs are read-only: the frozen `02`/`04` outputs and the `pert_data/` source trees.
The full per-observation table is gitignored; its sha256 lives in ``manifest.json``.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import gzip
import hashlib
import io
import json
import platform
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import control_matching as cmatch
import rules as ruleset

# --------------------------------------------------------------------------- #
# Paths / provenance
# --------------------------------------------------------------------------- #

TEMPENV_SRC = "pert_data/GSE262419/source"
TEMPENV_HASH = "pert_data/GSE262419/source/GSE262419_hash.csv.gz"
TEMPENV_PLATE = TEMPENV_SRC + "/GSE262419_Plate{:02d}.csv.gz"
TEMPENV_SERIES = "pert_data/GSE262419/source/GSE262419_series_matrix.txt.gz"
TAHOE_META = "pert_data/Tahoe-100M/source/metadata/sample_metadata.parquet"
TAHOE_GENE_META = "pert_data/Tahoe-100M/source/metadata/gene_metadata.parquet"
TAHOE_VIEW = "stage2/05_prnet_pretrain/output/pretrain_view/tahoe_cm_schema.h5ad"
TAHOE_VIEW_MANIFEST = "stage2/05_prnet_pretrain/output/pretrain_view/manifest.csv"
TAHOE_PB_MANIFEST = ("stage2/02_dataset_corpus/execution/build_record/"
                     "tahoe_pseudobulk_manifest.csv")
CM_VIEW = "stage2/_frozen/cm-data-v1/datasets-1/temp_env_all-smiles.h5ad"
CM_GENES = "stage2/_frozen/cm-data-v1/corpus/v1/genes.csv"

OUT_DIR = "stage2/06_prnet_transfer/06.2_qc/output"
FULL_TABLE = "observations.csv.gz"

N_PLATES_TEMPENV = 16


def repo_root(start: str | Path | None = None) -> Path:
    """Walk up to the directory holding .git (never Path.parents arithmetic)."""
    here = Path(start or __file__).resolve()
    for cand in [here, *here.parents]:
        if (cand / ".git").exists():
            return cand
    raise RuntimeError(f"no .git found above {here}")


def sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def read_h5ad_obs(path: str | Path) -> pd.DataFrame:
    """Decode only an h5ad's ``obs`` frame, with h5py — never ``anndata.read_h5ad``.

    ``anndata.read_h5ad(..., backed="r")`` is **not** a metadata-only read: anndata
    materialises every ``layers`` entry in memory even in backed mode.  On the frozen
    Tahoe view that is 4.3 GB of ``log1p_cp10k`` for a caller that only wants ``obs``
    (the "8.1 GB" in ``docs/agents/execution-environment.md`` is X plus that layer),
    which is the memory trap that killed the 2026-11-02 session.  Decoding just the
    ``obs`` group keeps the peak proportional to the metadata and never touches
    ``X``/``layers``; this is the pattern of
    ``06.3_input/execution/contract.py::read_h5ad_header``, copied rather than
    imported so this substep keeps its independently re-derivable execution tree
    (ADR 0030 §3).
    """
    import h5py

    try:                                            # anndata >= 0.11
        from anndata.io import read_elem
    except ImportError:                             # pragma: no cover - 0.10 reader
        from anndata.experimental import read_elem
    path = Path(path)
    with h5py.File(path, "r") as handle:
        if "obs" not in handle:
            raise RuntimeError(f"{path}: h5ad has no obs group")
        obs = read_elem(handle["obs"])
    if not isinstance(obs, pd.DataFrame):
        raise RuntimeError(f"{path}: obs decoded to {type(obs)!r}, not a DataFrame")
    return obs


def _now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def _git_head(root: Path) -> str:
    import subprocess
    try:
        return subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:  # pragma: no cover - only when git is absent
        return "unknown"


def _env_versions() -> dict:
    from importlib.metadata import version
    return {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "anndata": version("anndata"),
        "pyarrow": version("pyarrow"),
    }


# --------------------------------------------------------------------------- #
# Time provenance (the 24 h claim, per source)
# --------------------------------------------------------------------------- #

TIME_SOURCES = {
    "tempenv_paper_24h": {
        "source": "temp_env",
        "time_h": 24.0,
        "materialised_from": ("frozen view obs.time_h, which is a hardcoded loader constant "
                              "(stage2/_frozen/cm-data-v1/execution/cmdata/sources.py), "
                              "not a GEO field"),
        "basis": {
            "kind": "source_study",
            "citation": ("Tsai et al., Chem Res Toxicol 2024, doi:10.1021/acs.chemrestox.4c00193 "
                         "(PMC11691792)"),
            "url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11691792/",
            "quote": ("Gene Expression Assays: Following chemical exposure for 24 hrs, the wells "
                      "were completely aspirated and 10 uL/well of 2x lysis buffer ... was added"),
            "note": ("the same study's functional/cytotoxicity phenotypes used a separate 90-minute "
                     "exposure, so a 24 h transcriptomic well has no same-well viability readout"),
            "retrieved": "2026-11-02",
        },
    },
    "tahoe_preprint_24h": {
        "source": "tahoe",
        "time_h": None,
        "materialised_from": "NOT materialised: the atlas has no time column, so obs.time_h stays NaN",
        "basis": {
            "kind": "source_study",
            "citation": ("Zhang J, Ubas AA, de Borja R, Svensson V, ... Goodarzi H, Yu J. "
                         "Tahoe-100M ... bioRxiv 2025.02.20.639398v2; formal record Cell 2026, "
                         "doi:10.1016/j.cell.2026.08.035"),
            "url": "https://www.biorxiv.org/content/10.1101/2025.02.20.639398v2.full",
            "quote": ("we cultured spheroids consisting of a mixture of cell lines in suspension "
                      "and subjected each spheroid to individual treatments, including DMSO vehicle "
                      "controls. After 24 hours of drug treatment, spheroids were dissociated, "
                      "fixed, and profiled using the Parse GigaLab kit"),
            "figure_quote": ("Figure 1A legend: Each well receives a drug perturbation, and after "
                             "24h, tumors are dissociated into single cells"),
            "retrieval_attempts": {
                "biorxiv_v2_full_text": "HTTP 200, retrieved 2026-11-02",
                "cell_fulltext_S0092-8674(26)01008-1": "HTTP 403 (Cloudflare) on 2026-11-02",
                "cell_pdf": "HTTP 403 (Cloudflare) on 2026-11-02",
            },
            "note": ("the quoted sentences are the authors' own text in the v2 preprint of the same "
                     "study; the Cell version is the formal record and is cited for the record. "
                     "The atlas itself carries no time field, so no obs column is written."),
            "retrieved": "2026-11-02",
        },
    },
}


# --------------------------------------------------------------------------- #
# temp_env: raw design table + raw per-plate count matrices
# --------------------------------------------------------------------------- #

def tempenv_well_id(plate_id: int, well_raw: str) -> str:
    """The frozen view's own well id format: temp_env::Plate1-B02."""
    m = re.match(r"^([A-Z]+)0*(\d+)$", str(well_raw))
    if not m:
        raise ValueError(f"unparsable well id {well_raw!r}")
    return f"temp_env::Plate{int(plate_id)}-{m.group(1)}{int(m.group(2)):02d}"


def load_tempenv_design(root: Path) -> pd.DataFrame:
    h = pd.read_csv(root / TEMPENV_HASH)
    expect = ["No", "Plate_ID", "Well_ID", "Chemical_Index", "Column order",
              "Chemical_Concentration_uM", "Replicate_Number", "Solvent", "Chemical_name"]
    if list(h.columns) != expect:
        raise RuntimeError(f"temp_env hash table columns changed: {list(h.columns)}")
    h["observation_id"] = [tempenv_well_id(p, w) for p, w in zip(h.Plate_ID, h.Well_ID)]
    if h["observation_id"].duplicated().any():
        raise RuntimeError("duplicate temp_env observation ids")
    return h


def _tempenv_kind(name: str) -> str:
    return {"VEH": "vehicle_control_dmso", "Media": "control_media",
            "Blank": "technical_blank", "ET-1": "positive_control_et1"}.get(
                name, "treated_compound")


def tempenv_plate_counts(root: Path, design: pd.DataFrame):
    """Two streaming passes over the 16 raw plate matrices.

    Returns (per_well_metrics, gene_totals, n_panel_genes, retained_mask, fam153b_per_well).
    Nothing is held dense beyond one plate at a time.
    """
    gene_totals = None
    n_wells = 0
    panel_index = None
    fam_total = 0
    fam_wells_nonzero = 0
    for p in range(1, N_PLATES_TEMPENV + 1):
        df = pd.read_csv(root / TEMPENV_PLATE.format(p), index_col=0)
        if panel_index is None:
            panel_index = list(df.index)
        elif list(df.index) != panel_index:
            raise RuntimeError(f"plate {p}: probe panel differs from plate 1")
        X = df.to_numpy(dtype=np.int64)
        gene_totals = X.sum(1) if gene_totals is None else gene_totals + X.sum(1)
        n_wells += X.shape[1]
        for i in [i for i, g in enumerate(panel_index) if str(g).upper().startswith("FAM153B")]:
            fam_total += int(X[i].sum())
            fam_wells_nonzero += int((X[i] > 0).sum())

    retained = ruleset.author_gene_filter(gene_totals, n_wells,
                                          min_mean=ruleset.RULES["tempenv"]["thresholds"]["gene_filter_min_mean"])

    metrics = []
    for p in range(1, N_PLATES_TEMPENV + 1):
        df = pd.read_csv(root / TEMPENV_PLATE.format(p), index_col=0)
        X = df.to_numpy(dtype=np.int64)
        wells = [re.sub(r"^([A-Z]+)0*(\d+)$", r"\1\2", str(w).split("-")[-1]) for w in df.columns]
        metrics.append(pd.DataFrame({
            "Plate_ID": p,
            "well_norm": wells,
            "plate_well_raw": [str(w).split("-")[-1] for w in df.columns],
            "lib_size_panel": X.sum(0),
            "n_genes_panel": (X > 0).sum(0),
            "n_genes_retained_panel": (X[retained] > 0).sum(0),
            "lib_size_retained": X[retained].sum(0),
        }))
    per_well = pd.concat(metrics, ignore_index=True)
    fam = {"total_counts": fam_total, "n_wells_nonzero": fam_wells_nonzero,
           "panel_probe": [str(g) for g in panel_index if str(g).upper().startswith("FAM153B")],
           "n_wells_total": n_wells}
    return per_well, gene_totals, panel_index, retained, fam


def tempenv_correlations(root: Path, design: pd.DataFrame,
                         retained: np.ndarray) -> pd.DataFrame:
    """Mean Pearson r of every well against (a) all 192 DMSO wells, (b) its plate's DMSO wells.

    The statistic is Pearson correlation of *raw counts* over the author-filtered gene
    set — the only reading of the study's "mean correlation threshold of <0.85" that
    reproduces a >0.85 median among vehicle controls (log-space correlations of
    TempO-Seq counts sit near 0.6 and would discard 80% of the study).
    """
    blocks, meta = [], []
    for p in range(1, N_PLATES_TEMPENV + 1):
        df = pd.read_csv(root / TEMPENV_PLATE.format(p), index_col=0)
        X = df.to_numpy(dtype=np.float64)[retained]
        blocks.append(X)
        wells = [re.sub(r"^([A-Z]+)0*(\d+)$", r"\1\2", str(w).split("-")[-1]) for w in df.columns]
        hh = design[design.Plate_ID == p].set_index("Well_ID")
        meta.append(pd.DataFrame({"Plate_ID": p, "well_norm": wells,
                                  "Chemical_name": [hh.loc[w, "Chemical_name"] for w in wells]}))
    X = np.hstack(blocks)
    m = pd.concat(meta, ignore_index=True)

    def _z(a):
        a = a - a.mean(0, keepdims=True)
        s = a.std(0, keepdims=True)
        s[s == 0] = 1.0
        return a / s

    Z = _z(X)
    dmso = np.where(m.Chemical_name.to_numpy() == "VEH")[0]
    Zd = Z[:, dmso]
    n_genes = X.shape[0]
    mean_r_all = (Z.T @ Zd).mean(1) / n_genes
    # remove self-correlation for the DMSO wells themselves
    self_mask = np.zeros(len(m), dtype=bool)
    self_mask[dmso] = True
    mean_r_all[self_mask] = (mean_r_all[self_mask] * len(dmso) - 1.0) / (len(dmso) - 1)

    plate_r = np.full(len(m), np.nan)
    for p in range(1, N_PLATES_TEMPENV + 1):
        idx = np.where(m.Plate_ID.to_numpy() == p)[0]
        d = np.intersect1d(idx, dmso)
        if len(d) < 2:
            continue
        C = (Z[:, idx].T @ Z[:, d]) / n_genes
        sub = C[:, [list(d).index(x) for x in d]]
        for j, i in enumerate(idx):
            row = C[j].copy()
            if i in d:
                row[list(d).index(i)] = np.nan
            plate_r[i] = np.nanmean(row)
    return pd.DataFrame({"Plate_ID": m.Plate_ID, "well_norm": m.well_norm,
                         "meanr_vs_vehicle_pool": mean_r_all,
                         "meanr_vs_plate_vehicle": plate_r})


def build_tempenv(root: Path) -> tuple[pd.DataFrame, dict]:
    design = load_tempenv_design(root)
    per_well, gene_totals, panel_index, retained, fam = tempenv_plate_counts(root, design)
    corr = tempenv_correlations(root, design, retained)

    df = design.merge(per_well, left_on=["Plate_ID", "Well_ID"],
                      right_on=["Plate_ID", "well_norm"], how="left", validate="1:1")
    df = df.merge(corr, on=["Plate_ID", "well_norm"], how="left", validate="1:1")
    if df["lib_size_panel"].isna().any():
        raise RuntimeError("temp_env: a design well has no count-matrix column")

    df["role"] = [_tempenv_kind(n) for n in df.Chemical_name]
    df["plate"] = ["Plate%d" % p for p in df.Plate_ID]
    df["physical_well"] = df.Well_ID.astype(str)
    df["vehicle_raw"] = df.Solvent.astype(str)
    df["dose_um"] = df.Chemical_Concentration_uM.astype(float)
    df["n_genes_retained_by_filter"] = int(retained.sum())

    refs = {
        "tempenv_gene_filter_min_mean": ruleset.RULES["tempenv"]["thresholds"]["gene_filter_min_mean"],
        "tempenv_n_genes_retained": float(retained.sum()),
        "tempenv_sample_expression_floor": ruleset.RULES["tempenv"]["thresholds"]["sample_expression_floor"],
        "tempenv_vehicle_corr_floor": ruleset.RULES["tempenv"]["thresholds"]["vehicle_corr_floor"],
    }
    blanks = df.loc[df.role == "technical_blank", "n_genes_panel"].to_numpy()
    med, mad, thr = ruleset.blank_separation(blanks, k=ruleset.RULES["tempenv"]["thresholds"]["blank_mad_k"])
    refs.update({
        "tempenv_blank_ngenes_median": med,
        "tempenv_blank_ngenes_mad": mad,
        "tempenv_blank_ngenes_threshold": thr,
    })
    info = {
        "n_panel_genes": len(panel_index),
        "n_genes_retained": int(retained.sum()),
        "blank_ref": {"n": int(len(blanks)), "median": med, "mad": mad, "threshold": thr},
        "fam153b": {
            "panel_probe": fam["panel_probe"],
            "total_counts": float(fam["total_counts"]),
            "n_wells_nonzero": int(fam["n_wells_nonzero"]),
            "n_wells_total": int(fam["n_wells_total"]),
            "measured": bool(fam["total_counts"] > 0),
        },
    }
    return df, {"refs": refs, "info": info}


# --------------------------------------------------------------------------- #
# Tahoe: atlas per-well metadata + the frozen CM-schema view's obs
# --------------------------------------------------------------------------- #

def build_tahoe(root: Path) -> tuple[pd.DataFrame, dict]:
    sm = pd.read_parquet(root / TAHOE_META)
    if len(sm) != 1344:
        raise RuntimeError(f"tahoe sample_metadata rows changed: {len(sm)}")
    sm["plate"] = sm["plate"].astype(str)
    sm["physical_well"] = sm["sample"].astype(str)
    sm["drug_exact"] = sm["drug"].astype(str)

    # per-plate modified z on log10 well depth (H4), reference = the plate's 96 wells
    zs = []
    for plate, g in sm.groupby("plate"):
        x = np.log10(g["mean_tscp_count"].to_numpy(dtype=float))
        zs.append(pd.DataFrame({"physical_well": g.physical_well.to_numpy(),
                                "plate": plate,
                                "well_z_tscp": ruleset.modified_z(x),
                                "well_z_gene_count":
                                    ruleset.modified_z(np.log10(g["mean_gene_count"].to_numpy(dtype=float)))}))
    well_z = pd.concat(zs, ignore_index=True)

    # header-only read: never `anndata.read_h5ad(backed="r")` (see read_h5ad_obs)
    obs = read_h5ad_obs(root / TAHOE_VIEW)
    keep_cols = ["well", "plate", "cell_line_id", "cell_name", "is_control", "dose_um",
                 "drug_key", "compound_name", "n_cells", "lib_size", "n_genes_detected",
                 "replicate", "has_structure", "time_h", "vehicle", "control_group"]
    o = obs[keep_cols].copy()
    o.index = o.index.astype(str)
    # obs.well is "tahoe::<plate>:<physical well>:<cell line id>"
    parsed = o.index.to_series().str.extract(
        r"^tahoe::(?P<plate>plate\d+):(?P<physical_well>smp_\d+):(?P<cell_line>\S+)$")
    if parsed["physical_well"].isna().any():
        raise RuntimeError("tahoe: physical well could not be parsed from obs.well")
    o["physical_well"] = parsed["physical_well"].astype(str).to_numpy()
    if not np.array_equal(parsed["cell_line"].astype(str).to_numpy(),
                          o["cell_line_id"].astype(str).to_numpy()):
        raise RuntimeError("tahoe: parsed cell line disagrees with obs.cell_line_id")
    if not np.array_equal(parsed["plate"].astype(str).to_numpy(),
                          o["plate"].astype(str).to_numpy()):
        raise RuntimeError("tahoe: parsed plate disagrees with obs.plate")
    o["observation_id"] = o.index.astype(str)
    o["plate"] = o["plate"].astype(str)
    o["cell_line_id"] = o["cell_line_id"].astype(str)
    o["cell_name"] = o["cell_name"].astype(str)
    o = o.reset_index(drop=True)

    o = o.merge(sm[["physical_well", "plate", "mean_gene_count", "mean_tscp_count",
                    "mean_mread_count", "mean_pcnt_mito", "drug_exact"]],
                on=["physical_well", "plate"], how="left", validate="m:1")
    o = o.merge(well_z, on=["physical_well", "plate"], how="left", validate="m:1")
    if o["mean_tscp_count"].isna().any():
        raise RuntimeError("tahoe: an observation has no atlas per-well metadata")

    # the atlas's own control definition is the exact string DMSO_TF
    is_ctl_view = o["is_control"].to_numpy().astype(bool)
    is_ctl_atlas = (o["drug_exact"] == "DMSO_TF").to_numpy()
    if not np.array_equal(is_ctl_view, is_ctl_atlas):
        raise RuntimeError("tahoe: view is_control disagrees with the atlas drug == 'DMSO_TF'")
    o["role"] = np.where(is_ctl_view, "vehicle_control_dmso_tf", "treated_compound")

    # FAM153B mask: measured in the atlas or a missing slot?
    gm = pd.read_parquet(root / TAHOE_GENE_META)
    fam_alias = gm[gm["gene_symbol"].astype(str).str.upper() == "FAM153B"]
    cm_fam = "ENSG00000182230"
    info = {
        "n_obs": int(len(o)),
        "n_physical_wells": int(o["physical_well"].nunique()),
        "n_control_groups": int(o["control_group"].nunique()),
        "n_cell_lines": int(o["cell_line_id"].nunique()),
        "fam153b": {
            "cm_ensembl": cm_fam,
            "atlas_rows_for_symbol": fam_alias[["gene_symbol", "ensembl_id"]].to_dict("records"),
            "measured": bool((fam_alias["ensembl_id"].astype(str) == cm_fam).any()),
            "evidence": ("atlas gene_metadata maps the symbol FAM153B to "
                         f"{fam_alias['ensembl_id'].astype(str).tolist()}; {cm_fam} is absent, so "
                         "the CM slot is an unmeasured gene, not a measured zero"),
        },
    }
    return o, info


# --------------------------------------------------------------------------- #
# Assembly
# --------------------------------------------------------------------------- #

TABLE_COLUMNS = [
    "observation_id", "source", "role", "plate", "physical_well",
    "cell_line_id", "cell_name", "cell_type", "donor",
    "compound_name", "drug_key", "dose_um", "is_control", "vehicle",
    "time_h", "time_source_id",
    "raw_replicate", "raw_replicate_source", "frozen_view_replicate_index",
    "lib_size_panel", "n_genes_panel", "n_genes_retained_by_filter",
    "mean_count_per_gene_retained", "meanr_vs_vehicle_pool", "meanr_vs_plate_vehicle",
    "n_cells", "mean_gene_count", "mean_tscp_count", "mean_mread_count", "mean_pcnt_mito",
    "well_z_tscp", "well_z_gene_count",
    "lib_size_view", "n_genes_detected_view", "has_structure",
    "qc_status", "qc_reason", "qc_flags",
    "matched_control_key", "matched_control_n", "matched_control_n_total",
    "matched_control_ok", "matched_control_status",
    "response_label", "toxicity_label", "label_evidence",
]

LABEL_EVIDENCE = {
    "temp_env": ("none:temp_env.study_measured_cytotoxicity_in_a_separate_90min_assay_"
                 "(PMC11691792),not_on_the_24h_transcriptomic_wells;values_not_in_this_repo"),
    "tahoe": ("none:tahoe.atlas_reports_cytotoxicity/proliferation_as_a_quantified_phenotype_"
              "and_ships_n_cells_per_profile,no_per-well_viability_label_in_this_repo"),
}


def assemble_tempenv(root: Path) -> tuple[pd.DataFrame, dict, dict]:
    df, meta = build_tempenv(root)
    refs = meta["refs"]

    # frozen CM view: drug identity + the 18,164-axis QC columns (header-only read)
    vo = read_h5ad_obs(root / CM_VIEW)
    v = pd.DataFrame({
        "observation_id": vo.index.astype(str),
        "drug_key": vo["drug_key"].astype(str).to_numpy(),
        "compound_name_view": vo["compound_name"].astype(str).to_numpy(),
        "inchikey": vo["inchikey"].astype(str).to_numpy(),
        "lib_size_view": vo["lib_size"].to_numpy(),
        "n_genes_detected_view": vo["n_genes_detected"].to_numpy(),
        "frozen_view_replicate_index": vo["replicate"].to_numpy(),
    })
    df = df.merge(v, on="observation_id", how="left", validate="1:1")

    verdict = ruleset.evaluate_tempenv(df, refs)
    df = df.join(verdict)

    # per-plate robust-z diagnostics: reference = the plate's own DMSO wells
    # (informational flags; they never change qc_status)
    for col, flag in (("lib_size_panel", "plate_libsize_robust_z_gt3"),
                      ("n_genes_panel", "plate_ngenes_robust_z_gt3")):
        df["z_" + col] = np.nan
        for _plate, g in df.groupby("plate"):
            ref = np.log10(g.loc[g.role == "vehicle_control_dmso", col].to_numpy(dtype=float))
            mad = np.median(np.abs(ref - np.median(ref))) if len(ref) >= 3 else 0.0
            if not mad:
                continue
            df.loc[g.index, "z_" + col] = (
                (np.log10(g[col].to_numpy(dtype=float)) - np.median(ref)) / (1.4826 * mad))
        hit = df["z_" + col].abs() > 3
        df.loc[hit, "qc_flags"] = df.loc[hit, "qc_flags"].apply(
            lambda s: (s + ";" if s else "") + flag)

    out = pd.DataFrame({
        "observation_id": df.observation_id,
        "source": "temp_env",
        "role": df.role,
        "plate": df.plate,
        "physical_well": df.physical_well,
        "cell_line_id": "",
        "cell_name": "",
        "cell_type": "ipsc_cm",
        "donor": "NA",
        "compound_name": df.Chemical_name.astype(str),
        "drug_key": df.drug_key.fillna("").astype(str),
        "dose_um": df.dose_um,
        "is_control": df.role.isin(["vehicle_control_dmso", "control_media"]),
        "vehicle": df.vehicle_raw,
        "time_h": 24.0,
        "time_source_id": "tempenv_paper_24h",
        "raw_replicate": df.Replicate_Number.astype("Int64"),
        "raw_replicate_source": "GSE262419_hash.csv.gz:Replicate_Number",
        "frozen_view_replicate_index": df.frozen_view_replicate_index,
        "lib_size_panel": df.lib_size_panel,
        "n_genes_panel": df.n_genes_panel,
        "n_genes_retained_by_filter": df.n_genes_retained_by_filter,
        "mean_count_per_gene_retained": df.lib_size_retained / df.n_genes_retained_by_filter,
        "meanr_vs_vehicle_pool": df.meanr_vs_vehicle_pool,
        "meanr_vs_plate_vehicle": df.meanr_vs_plate_vehicle,
        "n_cells": np.nan,
        "mean_gene_count": np.nan,
        "mean_tscp_count": np.nan,
        "mean_mread_count": np.nan,
        "mean_pcnt_mito": np.nan,
        "well_z_tscp": np.nan,
        "well_z_gene_count": np.nan,
        "lib_size_view": df.lib_size_view,
        "n_genes_detected_view": df.n_genes_detected_view,
        "has_structure": df.drug_key.notna() & df.drug_key.fillna("").ne(""),
        "qc_status": df.qc_status,
        "qc_reason": df.qc_reason,
        "qc_flags": df.qc_flags,
    })
    out["response_label"] = "unknown"
    out["toxicity_label"] = "unknown"
    out["label_evidence"] = LABEL_EVIDENCE["temp_env"]

    matched = cmatch.match_controls(out, "temp_env")
    out = out.join(matched)
    return out[TABLE_COLUMNS], meta["info"], refs


def assemble_tahoe(root: Path) -> tuple[pd.DataFrame, dict, dict]:
    o, info = build_tahoe(root)
    refs = {
        "tahoe_n_cells_min": ruleset.RULES["tahoe"]["thresholds"]["n_cells_min"],
        "tahoe_source_n_cells_floor": ruleset.RULES["tahoe"]["thresholds"]["source_n_cells_floor"],
        "tahoe_excluded_cell_lines": ruleset.RULES["tahoe"]["thresholds"]["excluded_cell_lines"],
        "tahoe_mito_max_pcnt": ruleset.RULES["tahoe"]["thresholds"]["mito_max_pcnt"],
        "tahoe_min_mean_gene_count": ruleset.RULES["tahoe"]["thresholds"]["min_mean_gene_count"],
        "tahoe_well_z_abs_max": ruleset.RULES["tahoe"]["thresholds"]["well_z_abs_max"],
    }
    verdict = ruleset.evaluate_tahoe(o, refs)
    o = pd.concat([o, verdict], axis=1)

    out = pd.DataFrame({
        "observation_id": o.observation_id,
        "source": "tahoe",
        "role": o.role,
        "plate": o.plate.astype(str),
        "physical_well": o.physical_well,
        "cell_line_id": o.cell_line_id.astype(str),
        "cell_name": o.cell_name.astype(str),
        "cell_type": np.nan,
        "donor": np.nan,
        "compound_name": o.compound_name.astype(str),
        "drug_key": o.drug_key.astype(str),
        "dose_um": o.dose_um.astype(float),
        "is_control": o.role.eq("vehicle_control_dmso_tf"),
        "vehicle": o.vehicle.astype(str),
        "time_h": np.nan,
        "time_source_id": "tahoe_preprint_24h",
        "raw_replicate": pd.array([pd.NA] * len(o), dtype="Int64"),
        "raw_replicate_source": ("none: the atlas ships no replicate id; the physical well "
                                 "(atlas sample) is the replicate axis"),
        "frozen_view_replicate_index": o.replicate,
        "lib_size_panel": np.nan,
        "n_genes_panel": np.nan,
        "n_genes_retained_by_filter": np.nan,
        "mean_count_per_gene_retained": np.nan,
        "meanr_vs_vehicle_pool": np.nan,
        "meanr_vs_plate_vehicle": np.nan,
        "n_cells": o.n_cells,
        "mean_gene_count": o.mean_gene_count,
        "mean_tscp_count": o.mean_tscp_count,
        "mean_mread_count": o.mean_mread_count,
        "mean_pcnt_mito": o.mean_pcnt_mito,
        "well_z_tscp": o.well_z_tscp,
        "well_z_gene_count": o.well_z_gene_count,
        "lib_size_view": o.lib_size,
        "n_genes_detected_view": o.n_genes_detected,
        "has_structure": o.has_structure.astype(bool),
        "qc_status": o.qc_status,
        "qc_reason": o.qc_reason,
        "qc_flags": o.qc_flags,
    })
    out["response_label"] = "unknown"
    out["toxicity_label"] = "unknown"
    out["label_evidence"] = LABEL_EVIDENCE["tahoe"]
    matched = cmatch.match_controls(out, "tahoe")
    out = out.join(matched)
    return out[TABLE_COLUMNS], info, refs


# --------------------------------------------------------------------------- #
# Counts tables / summary
# --------------------------------------------------------------------------- #

def counts_by_role(table: pd.DataFrame) -> pd.DataFrame:
    g = (table.groupby(["source", "role", "qc_status"], observed=True)
         .agg(observations=("observation_id", "size"),
              physical_wells=("physical_well_uid", "nunique"))
         .reset_index())
    return g.sort_values(["source", "role", "qc_status"]).reset_index(drop=True)


def counts_by_plate_dose(table: pd.DataFrame) -> pd.DataFrame:
    t = table.copy()
    t["dose_um"] = t["dose_um"].fillna(-1.0)
    g = (t.groupby(["source", "plate", "role", "dose_um", "qc_status"], observed=True)
         .agg(observations=("observation_id", "size"),
              physical_wells=("physical_well_uid", "nunique"))
         .reset_index())
    g["dose_um"] = g["dose_um"].replace(-1.0, np.nan)
    return g.sort_values(["source", "plate", "role", "dose_um", "qc_status"]).reset_index(drop=True)


def counts_tahoe_by_line(table: pd.DataFrame) -> pd.DataFrame:
    t = table[table.source == "tahoe"]
    g = (t.groupby(["cell_name", "cell_line_id", "qc_status"], observed=True)
         .agg(observations=("observation_id", "size"),
              physical_wells=("physical_well_uid", "nunique"))
         .reset_index())
    return g.sort_values(["cell_name", "qc_status"]).reset_index(drop=True)


def exclusions_by_rule(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for source, block in ruleset.RULES.items():
        sub = table[table.source == block["dataset_id"]]
        for name, rule in block["rules"].items():
            code = rule.get("reason_code")
            if not code:
                continue
            hit = sub["qc_reason"].fillna("").str.contains(code, regex=False)
            rows.append({
                "source": block["dataset_id"], "rule_id": rule.get("id", name),
                "rule": name, "reason_code": code,
                "label": ruleset.REASON_LABELS.get(code, ""),
                "observations_flagged": int(hit.sum()),
                "physical_wells_flagged": int(sub.loc[hit, "physical_well_uid"].nunique()),
                "observations_total": int(len(sub)),
                "share_of_source": round(float(hit.sum()) / max(len(sub), 1), 6),
            })
    return pd.DataFrame(rows).sort_values(["source", "rule_id"]).reset_index(drop=True)


def primary_reason_counts(table: pd.DataFrame) -> pd.DataFrame:
    t = table[table.qc_status == "fail"].copy()
    t["primary_reason"] = t["qc_reason"].str.split(";").str[0]
    g = (t.groupby(["source", "primary_reason"], observed=True)
         .agg(observations=("observation_id", "size"),
              physical_wells=("physical_well_uid", "nunique"))
         .reset_index())
    return g.sort_values(["source", "primary_reason"]).reset_index(drop=True)


def plate_drug_confounding(root: Path, table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    # temp_env: raw chemical names and collapsed drug keys
    h = load_tempenv_design(root)
    for level, col in (("raw_name", "Chemical_name"), ("drug_key", "drug_key")):
        if level == "raw_name":
            d = h[~h.Chemical_name.isin(["VEH", "Media", "Blank", "ET-1"])]
            per = d.groupby(col)["Plate_ID"].nunique()
            n = int(d[col].nunique())
        else:
            d = table[(table.source == "temp_env") & (table.role == "treated_compound")]
            d = d[d.drug_key.fillna("") != ""]
            per = d.groupby("drug_key")["plate"].nunique()
            n = int(d.drug_key.nunique())
        rows.append({"source": "temp_env", "level": level, "n_entities": n,
                     "on_1_plate": int((per == 1).sum()), "on_2_plates": int((per == 2).sum()),
                     "on_gt2_plates": int((per > 2).sum()), "max_plates": int(per.max())})
    d = table[(table.source == "tahoe") & (table.role == "treated_compound")]
    per = d.groupby("drug_key")["plate"].nunique()
    rows.append({"source": "tahoe", "level": "drug_key", "n_entities": int(d.drug_key.nunique()),
                 "on_1_plate": int((per == 1).sum()), "on_2_plates": int((per == 2).sum()),
                 "on_gt2_plates": int((per > 2).sum()), "max_plates": int(per.max())})
    return pd.DataFrame(rows)


def control_tables(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    te = table[(table.source == "temp_env") & (table.role == "vehicle_control_dmso")]
    te = te[["plate", "physical_well", "observation_id", "raw_replicate", "lib_size_panel",
             "n_genes_panel", "meanr_vs_vehicle_pool", "qc_status", "qc_reason"]]
    te = te.assign(eligible_for_control_pool=te.qc_status.eq("pass"))
    te = te.sort_values(["plate", "physical_well"]).reset_index(drop=True)

    th = table[(table.source == "tahoe") & (table.role == "vehicle_control_dmso_tf")]
    th = th[["plate", "cell_line_id", "cell_name", "physical_well", "observation_id",
             "n_cells", "mean_pcnt_mito", "mean_tscp_count", "qc_status", "qc_reason"]]
    th = th.assign(eligible_for_control_pool=th.qc_status.eq("pass"))
    th = th.sort_values(["plate", "cell_line_id", "physical_well"]).reset_index(drop=True)
    return te, th


def matched_availability(table: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for source in ("temp_env", "tahoe"):
        t = table[table.source == source]
        t = t.assign(cell_line_id=t.cell_line_id.fillna(""))
        frames.append(cmatch.availability_table(t, source).assign(source=source))
    out = pd.concat(frames, ignore_index=True)
    front = ["source", "matched_control_key", "plate"]
    rest = [c for c in out.columns if c not in front]
    return out[front + rest]


# --------------------------------------------------------------------------- #
# Deterministic writers
# --------------------------------------------------------------------------- #

def write_csv_gz(df: pd.DataFrame, path: Path) -> dict:
    buf = io.StringIO()
    df.to_csv(buf, index=False, lineterminator="\n", float_format="%.6g")
    raw = buf.getvalue().encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as fh:
        gz = gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0, compresslevel=6)
        gz.write(raw)
        gz.close()
    return {"bytes": path.stat().st_size, "sha256_gz": sha256_file(path),
            "sha256_uncompressed": hashlib.sha256(raw).hexdigest(),
            "rows": int(len(df)), "columns": int(df.shape[1])}


def write_csv(df: pd.DataFrame, path: Path) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, lineterminator="\n", float_format="%.6g")
    return {"bytes": path.stat().st_size, "sha256": sha256_file(path),
            "rows": int(len(df)), "columns": int(df.shape[1])}


def write_json(obj, path: Path) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(obj, indent=2, sort_keys=False, ensure_ascii=False, default=str) + "\n"
    path.write_text(text, encoding="utf-8")
    return {"bytes": path.stat().st_size, "sha256": sha256_file(path)}


# --------------------------------------------------------------------------- #
# Commands
# --------------------------------------------------------------------------- #

def compute_reference_block(data: Path) -> tuple[dict, dict, dict]:
    """Derive the instantiated thresholds and the per-source provenance block."""
    _, te_meta = build_tempenv(data)
    _, th_info = build_tahoe(data)
    refs = dict(te_meta["refs"])
    refs.update({
        "tahoe_n_cells_min": ruleset.RULES["tahoe"]["thresholds"]["n_cells_min"],
        "tahoe_source_n_cells_floor": ruleset.RULES["tahoe"]["thresholds"]["source_n_cells_floor"],
        "tahoe_excluded_cell_lines": ruleset.RULES["tahoe"]["thresholds"]["excluded_cell_lines"],
        "tahoe_mito_max_pcnt": ruleset.RULES["tahoe"]["thresholds"]["mito_max_pcnt"],
        "tahoe_min_mean_gene_count": ruleset.RULES["tahoe"]["thresholds"]["min_mean_gene_count"],
        "tahoe_well_z_abs_max": ruleset.RULES["tahoe"]["thresholds"]["well_z_abs_max"],
    })
    info = {"temp_env": te_meta["info"], "tahoe": th_info}
    return refs, info, {}


def lock_payload(refs: dict, info: dict) -> dict:
    return {
        "rules_version": ruleset.RULES_VERSION,
        "forbidden_bases": ruleset.FORBIDDEN_BASES,
        "sources": ruleset.RULES,
        "instantiated_thresholds": refs,
        "time_sources": TIME_SOURCES,
        "derived_reference_info": info,
    }


def payload_sha256(payload: dict) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
                      default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def cmd_lock(root: Path, data: Path, out: Path) -> int:
    refs, info, _ = compute_reference_block(data)
    payload = lock_payload(refs, info)
    sha = payload_sha256(payload)
    inputs = {p: sha256_file(data / p) for p in (
        TEMPENV_HASH, TEMPENV_SERIES, TAHOE_META, TAHOE_GENE_META,
        TAHOE_VIEW_MANIFEST, TAHOE_PB_MANIFEST, CM_GENES)}
    doc = {
        "artifact": "stage2/06_prnet_transfer/06.2_qc/output/qc_rules.json",
        "rules_version": ruleset.RULES_VERSION,
        "locked_at": _now(),
        "rules_sha256": sha,
        "hash_scope": ("sha256 over the canonical JSON of {rules_version, forbidden_bases, "
                       "sources, instantiated_thresholds, time_sources, derived_reference_info}"),
        "model_results_seen": False,
        "model_results_note": ("no model has been trained or evaluated for stage2/06 at lock "
                               "time; no DEG table, effect size, response statistic or model "
                               "metric was read by the locking code"),
        "forbidden_bases": ruleset.FORBIDDEN_BASES,
        "data_root": str(data),
        "git_head": _git_head(root),
        "environment": _env_versions(),
        "input_sha256": inputs,
        "time_sources": TIME_SOURCES,
        "instantiated_thresholds": refs,
        "derived_reference_info": info,
        "sources": ruleset.RULES,
    }
    meta = write_json(doc, out / "qc_rules.json")
    print(f"locked rules {ruleset.RULES_VERSION} sha256={sha[:16]}… -> {out/'qc_rules.json'}")
    print(f"  temp_env: {info['temp_env']['n_genes_retained']} genes retained, "
          f"blank floor {refs['tempenv_blank_ngenes_threshold']:.0f} genes")
    print("  " + json.dumps(meta))
    return 0


def load_lock(out: Path) -> dict:
    path = out / "qc_rules.json"
    if not path.exists():
        raise SystemExit(f"no lock at {path}; run --lock first")
    return json.loads(path.read_text(encoding="utf-8"))


def build_all(root: Path, data: Path, out: Path) -> int:
    lock = load_lock(out)
    refs, info, _ = compute_reference_block(data)
    payload = lock_payload(refs, info)
    sha = payload_sha256(payload)
    if sha != lock["rules_sha256"]:
        raise SystemExit("REFUSING TO BUILD: the rule payload changed after the lock "
                         f"(locked {lock['rules_sha256'][:16]}…, now {sha[:16]}…). "
                         "Re-lock deliberately.")
    built_at = _now()
    te, te_info, te_refs = assemble_tempenv(data)
    th, th_info, th_refs = assemble_tahoe(data)
    table = pd.concat([te, th], ignore_index=True)
    table["physical_well_uid"] = table["plate"].astype(str) + ":" + table["physical_well"].astype(str)
    table = table.sort_values(["source", "observation_id"]).reset_index(drop=True)
    if table["observation_id"].duplicated().any():
        raise RuntimeError("duplicate observation_id across sources")

    out.mkdir(parents=True, exist_ok=True)
    files = {}
    files[FULL_TABLE] = write_csv_gz(table, out / FULL_TABLE)
    files["counts_by_source_role.csv"] = write_csv(counts_by_role(table),
                                                   out / "counts_by_source_role.csv")
    files["counts_by_plate_dose.csv"] = write_csv(counts_by_plate_dose(table),
                                                  out / "counts_by_plate_dose.csv")
    files["counts_tahoe_by_cell_line.csv"] = write_csv(counts_tahoe_by_line(table),
                                                       out / "counts_tahoe_by_cell_line.csv")
    files["exclusions_by_rule.csv"] = write_csv(exclusions_by_rule(table),
                                                out / "exclusions_by_rule.csv")
    files["primary_reason_counts.csv"] = write_csv(primary_reason_counts(table),
                                                   out / "primary_reason_counts.csv")
    files["plate_drug_confounding.csv"] = write_csv(plate_drug_confounding(data, table),
                                                    out / "plate_drug_confounding.csv")
    files["matched_control_availability.csv"] = write_csv(matched_availability(table),
                                                          out / "matched_control_availability.csv")
    te_ctrl, th_ctrl = control_tables(table)
    files["temp_env_dmso_controls.csv"] = write_csv(te_ctrl, out / "temp_env_dmso_controls.csv")
    files["tahoe_dmso_controls.csv"] = write_csv(th_ctrl, out / "tahoe_dmso_controls.csv")

    summary = build_summary(table, lock, te_info, th_info, refs, built_at)
    files["qc_summary.json"] = write_json(summary, out / "qc_summary.json")
    handoff = build_handoff(table, lock, files, refs, te_info, th_info)
    files["handoff.json"] = write_json(handoff, out / "handoff.json")

    manifest = {
        "artifact": f"{OUT_DIR}/{FULL_TABLE}",
        "table": files[FULL_TABLE],
        "sha256_uncompressed_csv": files[FULL_TABLE]["sha256_uncompressed"],
        "sha256_of_full_per_observation_table": files[FULL_TABLE]["sha256_uncompressed"],
        "note": ("the full per-observation table is gitignored; this manifest entry is its "
                 "content identity (sha256 over the uncompressed CSV bytes)"),
        "rules_version": lock["rules_version"],
        "rules_sha256": lock["rules_sha256"],
        "locked_at": lock["locked_at"],
        "built_at": built_at,
        "git_head": _git_head(root),
        "environment": _env_versions(),
        "rows": int(len(table)),
        "columns": int(table.shape[1]),
        "column_order": list(table.columns),
        "rows_by_source": {k: int(v) for k, v in table.source.value_counts().items()},
        "input_sha256": lock["input_sha256"],
        "files": files,
    }
    files["manifest.json"] = write_json(manifest, out / "manifest.json")
    print(f"built {len(table):,} observations -> {out}")
    for name, meta in files.items():
        print(f"  {name:34s} {meta.get('rows',''):>8} {meta.get('bytes',0):>12,} B  "
              f"{meta.get('sha256', meta.get('sha256_gz',''))[:12]}")
    return 0


def build_summary(table: pd.DataFrame, lock: dict, te_info: dict, th_info: dict,
                  refs: dict, built_at: str) -> dict:
    te = table[table.source == "temp_env"]
    th = table[table.source == "tahoe"]

    def block(df):
        roles = {}
        for role, g in df.groupby("role", observed=True):
            roles[role] = {
                "observations": int(len(g)),
                "physical_wells": int(g.physical_well_uid.nunique()),
                "qc_pass": int((g.qc_status == "pass").sum()),
                "qc_fail": int((g.qc_status == "fail").sum()),
                "unresolved": int((g.qc_status == "unresolved").sum()),
            }
        return roles

    te_tr = te[te.role == "treated_compound"]
    th_tr = th[th.role == "treated_compound"]
    te_kept = te[te.qc_status == "pass"]
    th_kept = th[th.qc_status == "pass"]

    def unmatched(df):
        arm = df[df.role.isin(cmatch.ARMS[df.source.iloc[0]])]
        kept = arm[arm.qc_status == "pass"]
        return {
            "arm_observations": int(len(arm)),
            "arm_qc_pass": int(len(kept)),
            "qc_pass_without_matched_control": int((kept.matched_control_n.fillna(0) == 0).sum()),
            "qc_pass_with_1_control": int((kept.matched_control_n.fillna(0) == 1).sum()),
            "qc_pass_with_2plus_controls": int((kept.matched_control_n.fillna(0) >= 2).sum()),
        }

    # sensitivity: the source floor of 50 cells instead of the project floor of 100
    th_n = th["n_cells"].astype(float)
    sens50 = int(((th_n >= refs["tahoe_source_n_cells_floor"])
                  & ~th.cell_name.isin(refs["tahoe_excluded_cell_lines"])
                  & (th.mean_pcnt_mito <= refs["tahoe_mito_max_pcnt"])
                  & (th.mean_gene_count >= refs["tahoe_min_mean_gene_count"])
                  & (th.well_z_tscp.abs() <= refs["tahoe_well_z_abs_max"])).sum())

    te_mean_all = te.lib_size_panel / te_info["n_panel_genes"]
    return {
        "artifact": "stage2/06_prnet_transfer/06.2_qc/output/qc_summary.json",
        "rules_version": lock["rules_version"],
        "rules_sha256": lock["rules_sha256"],
        "locked_at": lock["locked_at"],
        "built_at": built_at,
        "model_results_seen": False,
        "observation_unit": {
            "temp_env": "physical well (one TempO-Seq library)",
            "tahoe": "(plate, physical well, cell_line) pseudobulk profile",
        },
        "temp_env": {
            "raw_wells": int(len(te)),
            "physical_wells": int(te.physical_well_uid.nunique()),
            "plates": int(te.plate.nunique()),
            "roles": block(te),
            "treated_observations": int(len(te_tr)),
            "treated_physical_wells": int(te_tr.physical_well_uid.nunique()),
            "treated_qc_pass": int((te_tr.qc_status == "pass").sum()),
            "treated_by_dose": {
                str(d): {
                    "observations": int(len(g)),
                    "qc_pass": int((g.qc_status == "pass").sum()),
                    "qc_fail": int((g.qc_status == "fail").sum()),
                } for d, g in te_tr.groupby("dose_um")},
            "control_arm": unmatched(te),
            "dmso_controls_total": int((te.role == "vehicle_control_dmso").sum()),
            "dmso_controls_eligible": int(((te.role == "vehicle_control_dmso")
                                           & (te.qc_status == "pass")).sum()),
            "plate_vehicle_strata": int(te.groupby(["plate", "vehicle"], observed=True).ngroups),
            "gene_panel": te_info["n_panel_genes"],
            "genes_retained_by_author_filter": te_info["n_genes_retained"],
            "blank_reference": te_info["blank_ref"],
            "sensitivity_T1_all_panel_genes_removed": int((te_mean_all < refs[
                "tempenv_sample_expression_floor"]).sum()),
            "meanr_vs_vehicle_pool": {
                "median": float(te.meanr_vs_vehicle_pool.median()),
                "n_below_source_floor": int((te.meanr_vs_vehicle_pool < refs[
                    "tempenv_vehicle_corr_floor"]).sum()),
            },
            "fam153b": te_info["fam153b"],
            "kept_observations": int(len(te_kept)),
            "sensitivity_kept_without_T3": int(
                ((te.qc_status == "pass")
                 | te.qc_reason.fillna("").str.split(";").apply(
                     lambda p: bool(p) and set(p) == {"derived_blank_separation_floor"})).sum()),
        },
        "tahoe": {
            "raw_wells": int(th.physical_well_uid.nunique()),
            "observations": int(len(th)),
            "plates": int(th.plate.nunique()),
            "control_groups": int(th.matched_control_key.nunique()),
            "cell_lines": int(th.cell_line_id.nunique()),
            "roles": block(th),
            "treated_observations": int(len(th_tr)),
            "treated_qc_pass": int((th_tr.qc_status == "pass").sum()),
            "control_arm": unmatched(th),
            "dmso_control_observations": int((th.role == "vehicle_control_dmso_tf").sum()),
            "dmso_control_observations_eligible": int(((th.role == "vehicle_control_dmso_tf")
                                                       & (th.qc_status == "pass")).sum()),
            "control_groups_with_0_eligible_controls": int(
                (th[th.role == "vehicle_control_dmso_tf"]
                 .groupby("matched_control_key", observed=True)["qc_status"]
                 .apply(lambda s: (s == "pass").sum() == 0)).sum()),
            "n_cells_sensitivity_kept_at_source_floor_50": sens50,
            "well_metrics": {
                "mean_pcnt_mito_max": float(th.mean_pcnt_mito.max()),
                "mean_gene_count_min": float(th.mean_gene_count.min()),
                "n_wells_z_tscp_gt3": int(th.drop_duplicates("physical_well")
                                          .well_z_tscp.abs().gt(3).sum()),
            },
            "fam153b": th_info["fam153b"],
            "kept_observations": int(len(th_kept)),
        },
        "published_aggregate_check": {
            "temp_env_source_study": ("4,379 of 4,558 samples (96%) and 11,039 genes remained "
                                      "after both QC steps (PMC11691792)"),
            "temp_env_this_audit": (f"{int((te.qc_status=='pass').sum())} of {len(te)} wells pass; "
                                    f"{te_info['n_genes_retained']} genes pass the author gene "
                                    "filter"),
            "deviation": ("the study reports only the aggregate; per-sample membership is not "
                          "deposited, so our rule cannot be verified sample-by-sample against "
                          "theirs. Our filter is reproducible from the deposited counts."),
        },
        "plate_drug_confounding_consequence": (
            "459/460 temp_env drug keys sit on exactly one plate, so a drug split is not a plate "
            "split: every fold's training drugs still come from most plates, and the per-plate "
            "DMSO pool of a fold must be taken from the plates its training drugs touch, never "
            "from a held-out plate's controls."),
    }


def build_handoff(table: pd.DataFrame, lock: dict, files: dict, refs: dict,
                  te_info: dict, th_info: dict) -> dict:
    te = table[table.source == "temp_env"]
    th = table[table.source == "tahoe"]
    return {
        "artifact": "stage2/06_prnet_transfer/06.2_qc/output/handoff.json",
        "version": "handoff-v1",
        "consumer": "ticket #3 (frozen two-source input, five folds, consensus controls)",
        "rules_version": lock["rules_version"],
        "rules_sha256": lock["rules_sha256"],
        "locked_at": lock["locked_at"],
        "observation_table": {
            "path": f"{OUT_DIR}/{FULL_TABLE}",
            "gitignored": True,
            "sha256_uncompressed_csv": files[FULL_TABLE]["sha256_uncompressed"],
            "sha256_gz": files[FULL_TABLE]["sha256_gz"],
            "rows": files[FULL_TABLE]["rows"],
            "key": "observation_id",
            "read": "pandas.read_csv(path) — no index column",
        },
        "kept_observation_ids": {
            "definition": "qc_status == 'pass' (QC only; pharmacology eligibility is ticket #1)",
            "temp_env": {
                "n": int((te.qc_status == "pass").sum()),
                "by_role": {r: int(((te.role == r) & (te.qc_status == "pass")).sum())
                            for r in sorted(te.role.unique())},
            },
            "tahoe": {
                "n": int((th.qc_status == "pass").sum()),
                "by_role": {r: int(((th.role == r) & (th.qc_status == "pass")).sum())
                            for r in sorted(th.role.unique())},
            },
        },
        "matched_control_join_keys": {
            "temp_env": {
                "key": ["plate", "vehicle"],
                "vehicle_value": "DMSO",
                "control_rows": "role == 'vehicle_control_dmso' and qc_status == 'pass'",
                "table": f"{OUT_DIR}/temp_env_dmso_controls.csv",
                "rule": ("same physical plate and same vehicle; Media is a different arm and "
                         "Blank is not a control"),
            },
            "tahoe": {
                "key": ["plate", "cell_line_id"],
                "control_rows": "role == 'vehicle_control_dmso_tf' and qc_status == 'pass'",
                "table": f"{OUT_DIR}/tahoe_dmso_controls.csv",
                "rule": ("same plate and same cell line only; never average across the ~50 lines "
                         "of a plate; the two control profiles of a group are two wells, not "
                         "two donors"),
            },
        },
        "eligible_dmso_rows_per_plate": {
            "temp_env": {
                "columns": ["plate", "physical_well", "observation_id", "raw_replicate",
                            "qc_status", "eligible_for_control_pool"],
                "note": "12 wells per plate before QC; ticketed counts per plate in the table",
            },
            "tahoe": {
                "columns": ["plate", "cell_line_id", "physical_well", "observation_id", "n_cells",
                            "qc_status", "eligible_for_control_pool"],
                "note": "2 wells per (plate, cell line) group before QC; 700 groups",
            },
        },
        "control_availability": {
            "temp_env": {
                "arm_rows": int(te.role.isin(cmatch.ARMS["temp_env"]).sum()),
                "arm_qc_pass_without_matched_control": int(
                    ((te.role.isin(cmatch.ARMS["temp_env"])) & (te.qc_status == "pass")
                     & (te.matched_control_n.fillna(0) == 0)).sum()),
            },
            "tahoe": {
                "arm_rows": int((th.role == "treated_compound").sum()),
                "arm_qc_pass_without_matched_control": int(
                    ((th.role == "treated_compound") & (th.qc_status == "pass")
                     & (th.matched_control_n.fillna(0) == 0)).sum()),
            },
            "table": f"{OUT_DIR}/matched_control_availability.csv",
        },
        "qc_flags": {
            "columns": ["lib_size_panel", "n_genes_panel", "mean_count_per_gene_retained",
                        "meanr_vs_vehicle_pool", "meanr_vs_plate_vehicle", "n_cells",
                        "mean_gene_count", "mean_tscp_count", "mean_mread_count",
                        "mean_pcnt_mito", "well_z_tscp", "well_z_gene_count",
                        "lib_size_view", "n_genes_detected_view", "qc_status", "qc_reason",
                        "qc_flags"],
            "notes": [
                "lib_size_panel / n_genes_panel are computed here from the 22,537-probe raw "
                "plate matrices; lib_size_view / n_genes_detected_view are the frozen 18,164-gene "
                "view columns — different gene axes, never interchangeable",
                "qc_reason is ';'-joined reason codes; the first element is the primary reason",
                "qc_flags are informational and never change qc_status",
                "raw_replicate is the source's own replicate field (temp_env); Tahoe has none — "
                "frozen_view_replicate_index is a generated index and must never be cited as "
                "biology",
            ],
        },
        "gene_measurement_mask": {
            "FAM153B": {
                "temp_env": {
                    "status": "measured",
                    "evidence": te_info["fam153b"],
                    "meaning": "a measured, non-zero gene: a zero in a well is a real zero",
                },
                "tahoe": {
                    "status": "missing" if not th_info["fam153b"]["measured"] else "measured",
                    "evidence": th_info["fam153b"],
                    "meaning": ("the CM axis slot is unmeasured in the atlas; it must be masked "
                                "and never read as a measured zero"),
                },
            },
            "consumer_action": ("carry a boolean gene_measured vector for the 18,164 axis; "
                                "FAM153B is False for Tahoe rows and True for temp_env rows"),
        },
        "routing": {
            "temp_env": {
                "treated_compound": "training/scoring arm",
                "vehicle_control_dmso": "the only vehicle control; source of the consensus pool",
                "control_media": "Media is NOT a vehicle control — never merge it into a DMSO pool",
                "technical_blank": "no cells; never a control",
                "positive_control_et1": ("ET-1 at 0.01 uM is a peptide positive control, not a "
                                         "small-molecule treatment; excluded from the training arm"),
            },
            "tahoe": {
                "treated_compound": "training/scoring arm",
                "vehicle_control_dmso_tf": "vehicle control (exact drug == 'DMSO_TF')",
                "no_structure": "2 drugs have no resolvable structure: dropped by the dose "
                                "encoder, counted as an input-layer exclusion, not a QC failure",
            },
        },
        "labels": {
            "response_label": "unknown for every row (no evidence source)",
            "toxicity_label": "unknown for every row (no per-well viability readout in this repo)",
            "rule": ("never infer toxicity or low response from transcriptomic magnitude; the "
                     "study-level readouts exist upstream but are not deposited here, see "
                     "README"),
        },
        "unresolved": {
            "count": int((table.qc_status == "unresolved").sum()),
            "reason_codes": ["unresolved_missing_metric"],
        },
        "not_covered_by_this_ticket": [
            "pharmacological eligibility (ticket #1) — join on drug_key",
            "the fold assignment (ticket #3)",
            "the consensus c* value itself (ticket #3 builds it from the eligible DMSO rows)",
            "any dose-response or toxicity label (tickets #3/#9, must be pre-registered)",
        ],
    }


def cmd_verify(root: Path, out: Path) -> int:
    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    ok = True
    for name, meta in manifest["files"].items():
        p = out / name
        sha = sha256_file(p)
        want = meta.get("sha256") or meta.get("sha256_gz")
        good = sha == want and p.stat().st_size == meta.get("bytes")
        ok &= good
        print(f"[{'PASS' if good else 'FAIL'}] {name} sha256={sha[:12]} bytes={p.stat().st_size}")
    with gzip.open(out / FULL_TABLE, "rb") as fh:
        raw = fh.read()
    content = hashlib.sha256(raw).hexdigest()
    good = content == manifest["sha256_uncompressed_csv"]
    ok &= good
    print(f"[{'PASS' if good else 'FAIL'}] table content sha256={content[:12]}")
    table = pd.read_csv(out / FULL_TABLE, low_memory=False)
    print(f"       parsed rows={len(table):,} columns={table.shape[1]}")
    invariants = [
        ("temp_env physical wells == 4558", int((table.source == "temp_env").sum()) == 4558),
        ("tahoe observations == 67018", int((table.source == "tahoe").sum()) == 67018),
        ("tahoe physical wells == 1344",
         table[table.source == "tahoe"].physical_well_uid.nunique() == 1344),
        ("tahoe control groups == 700",
         table[(table.source == "tahoe")
               & (table.role == "vehicle_control_dmso_tf")].matched_control_key.nunique() == 700),
        ("tahoe DMSO_TF observations == 1400",
         int(((table.source == "tahoe") & (table.role == "vehicle_control_dmso_tf")).sum()) == 1400),
        ("no observation is both pass and unexplained",
         not table[table.qc_status == "pass"].qc_reason.fillna("").ne("").any()),
        ("every failure has a reason code",
         table[table.qc_status == "fail"].qc_reason.fillna("").ne("").all()),
    ]
    for label, good in invariants:
        ok &= bool(good)
        print(f"[{'PASS' if good else 'FAIL'}] {label}")
    print("verify:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--lock", action="store_true", help="freeze the rule payload")
    ap.add_argument("--build", action="store_true", help="apply the locked rules")
    ap.add_argument("--verify", action="store_true", help="re-hash and re-check every artifact")
    ap.add_argument("--root", default=None,
                    help="code root, used for .git provenance (default: walk up to .git)")
    ap.add_argument("--data-root", default=None,
                    help="root that holds pert_data/ and stage2/*/output (default: --root; "
                         "the worktree has no data, so pass the main checkout)")
    ap.add_argument("--out-root", default=None,
                    help="root the output/ directory hangs from (default: the data root)")
    args = ap.parse_args(argv)
    root = Path(args.root).resolve() if args.root else repo_root()
    data = Path(args.data_root).resolve() if args.data_root else root
    if not (data / TEMPENV_HASH).exists():
        raise SystemExit(f"no inputs under --data-root {data}: expected {TEMPENV_HASH}. "
                         "Pass --data-root /path/to/main/checkout.")
    out_root = Path(args.out_root).resolve() if args.out_root else data
    out = out_root / OUT_DIR
    n = sum([args.lock, args.build, args.verify])
    if n != 1:
        ap.error("pass exactly one of --lock / --build / --verify")
    if args.lock:
        return cmd_lock(root, data, out)
    if args.build:
        return build_all(root, data, out)
    return cmd_verify(root, out)


if __name__ == "__main__":
    sys.exit(main())
