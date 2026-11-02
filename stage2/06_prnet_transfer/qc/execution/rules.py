"""Locked QC rule set for the two PRnet-formal-transfer sources.

Everything here is *data* plus two pure decision functions. Nothing in this module
reads a count matrix, a DEG table, a model metric or a response statistic: the
rules are depth / detection / reproducibility / cell-count rules only.

A threshold may exist only under one of two evidence kinds:

* ``source_study`` — quoted from the deposited study record (GEO series matrix) or
  the study's publication/supplement, with URL, quote and retrieval date;
* ``source_convention`` / ``project_derived`` — an independently documented
  derivation (author's own outlier convention applied to a new statistic, a
  technical-blank control distribution, or a density measurement recorded in an
  ADR), with the derivation written out and the alternative reported.

Forbidden bases (never a threshold, never a filter): DEG counts, effect size,
response strength, model performance. They are listed in ``FORBIDDEN_BASES`` and
asserted absent from every rule's quoted evidence by the test suite.

``--lock`` (see ``build_qc.py``) instantiates the derived thresholds from the QC
inputs once, writes them to ``output/qc_rules.json`` with a ``locked_at``
timestamp and the canonical sha256 of this module's rule payload, and every later
build asserts the payload hash is unchanged.
"""

from __future__ import annotations

import hashlib
import json

import numpy as np
import pandas as pd

RULES_VERSION = "qc-rules-v1.0.0"

FORBIDDEN_BASES = [
    "DEG counts / number of differentially expressed genes",
    "response strength / effect size / E-distance",
    "model performance (loss, AUPRC, correlation with truth)",
]

# --------------------------------------------------------------------------- #
# Evidence records (quoted, with retrieval date: 2026-11-02 Asia/Shanghai)
# --------------------------------------------------------------------------- #

EV_TEMPENV_GEO = {
    "kind": "source_study",
    "citation": ("Tsai H-H, Ford LC, Burnett SD, Dickey AN, Wright FA, Chiu WA, Rusyn I. "
                 "Informing Hazard Identification and Risk Characterization of Environmental "
                 "Chemicals ... Chem Res Toxicol 2024. doi:10.1021/acs.chemrestox.4c00193 "
                 "(PMID 39046974; GSE262419)"),
    "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262419 "
           "(series matrix field !Sample_data_processing)",
    "retrieved": "2026-11-02",
    "source_file": "pert_data/GSE262419/source/GSE262419_series_matrix.txt.gz",
}

EV_TEMPENV_PAPER = {
    "kind": "source_study",
    "citation": ("Tsai et al., Chem Res Toxicol 2024, doi:10.1021/acs.chemrestox.4c00193, "
                 "PMC11691792 (full text via NCBI E-utilities efetch, db=pmc)"),
    "url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11691792/",
    "retrieved": "2026-11-02",
}

EV_TAHOE_PREPRINT = {
    "kind": "source_study",
    "citation": ("Zhang J, Ubas AA, de Borja R, Svensson V, ... Goodarzi H, Yu J. "
                 "Tahoe-100M: A Giga-Scale Single-Cell Perturbation Atlas for Context-Dependent "
                 "Gene Function and Cellular Modeling. bioRxiv 2025.02.20.639398v2 "
                 "(author version of the formal record: Cell 2026, doi:10.1016/j.cell.2026.08.035)"),
    "url": "https://www.biorxiv.org/content/10.1101/2025.02.20.639398v2.full",
    "retrieved": "2026-11-02",
    "note": ("The Cell version returned HTTP 403 to every fetch attempt from this host on "
             "2026-11-02; the quoted sentences are verified in the authors' own v2 text of the "
             "same study, and the formal publication is cited for the record."),
}

EV_ADR0029 = {
    "kind": "project_derived",
    "citation": ("ADR 0029 §Why pseudobulk, measured: the 100-cell view threshold is set from "
                 "measured pseudobulk density against CM wells"),
    "url": "docs/adr/0029-tahoe-pseudobulk-pretrain-corpus.md",
    "retrieved": "2026-11-02",
    "quote": ("Tahoe pseudobulk, < 100 cells (9.5 %): ~50 cells, 11,485 / 17,935 genes detected, "
              "0.36 zero fraction, 0.04 M counts — 10–20× shallower than a CM well "
              "(temp_env 0.35–0.61 M counts)"),
}

EV_IGLEWICZ = {
    "kind": "source_convention",
    "citation": ("Iglewicz B, Hoaglin DC. How to Detect and Handle Outliers. ASQC Quality Press, "
                 "1993: modified z-score M = 0.6745·(x − median)/MAD with |M| > 3.5 as the flag "
                 "(equivalently |x − median| > 3·1.4826·MAD)"),
    "url": "https://doi.org/10.1080/00401706.1994.10485840 (review of the ASQC volume)",
    "retrieved": "2026-11-02",
}

# --------------------------------------------------------------------------- #
# The rule set
# --------------------------------------------------------------------------- #

RULES = {
    "tempenv": {
        "dataset_id": "temp_env",
        "observation_unit": "physical well (one TempO-Seq library); 4,558 raw wells",
        "rules": {
            "gene_filter": {
                "scope": "preprocessing step, not an exclusion rule",
                "statistic": "mean raw count per panel gene across all 4,558 wells",
                "threshold": 10.0,
                "op": "keep gene if mean >= 10",
                "evidence": {
                    **EV_TEMPENV_GEO,
                    "quote": ("Genes with fewer than 10 mean counts across all samples "
                              "(n = 4558) were excluded"),
                    "reported_by_authors": "11,039 genes remained after both QC steps",
                },
                "reason_code": None,
            },
            "sample_expression_floor": {
                "id": "T1",
                "statistic": "mean raw count per gene, computed after the gene filter "
                             "(lib_size_panel / n_genes_retained)",
                "threshold": 15.0,
                "op": "exclude well if mean < 15",
                "evidence": {
                    **EV_TEMPENV_GEO,
                    "quote": ("Samples with an average expression across all genes lower than 15 "
                              "were also removed."),
                    "reading": ("the two steps are listed in order, so the sample mean is taken "
                                "over the genes retained by the author gene filter"),
                    "alternative_reading": {
                        "statistic": "mean raw count over all 22,537 panel probes",
                        "threshold": 15.0,
                        "consequence": "stricter; the concrete well count is reported in "
                                       "qc_summary.json as sensitivity_T1_all_panel_genes",
                    },
                },
                "reason_code": "src_sample_expression_floor",
            },
            "vehicle_reproducibility": {
                "id": "T2",
                "statistic": ("mean Pearson r of the well's raw-count vector (author-filtered "
                              "genes) against every other DMSO vehicle control well, all 16 "
                              "plates pooled"),
                "threshold": 0.85,
                "op": "exclude DMSO vehicle-control well if mean r < 0.85",
                "evidence": {
                    **EV_TEMPENV_GEO,
                    "quote": ("Outliers in vehicle control samples were identified and excluded "
                              "using a mean correlation threshold of <0.85"),
                    "paper_quote": ("chemical effects were evaluated by comparing them to the "
                                    "plate-specific DMSO controls ... Outliers in vehicle control "
                                    "samples were identified and excluded using a mean correlation "
                                    "threshold of <0.85 ... a total of 4,379 samples (96% of "
                                    "4,558 samples) and 11,039 genes remained"),
                    "extension_policy": ("the study applies this floor to vehicle controls only; "
                                         "applying it to treated wells would be a "
                                         "response-correlated filter (a strong responder genuinely "
                                         "sits further from vehicle), so treated wells are "
                                         "flagged, never excluded, and the count is reported"),
                },
                "reason_code": "src_vehicle_reproducibility",
            },
            "blank_separation_floor": {
                "id": "T3",
                "statistic": "n_genes_detected over the 22,537-probe panel",
                "derivation": ("threshold = median + 3 · 1.4826 · MAD over the 32 technical Blank "
                               "wells (no cells, same plates and same sequencing run); a real cell "
                               "library must detect materially more genes than a no-cell well"),
                "threshold": None,
                "op": "exclude well if n_genes_detected <= threshold",
                "evidence": {
                    "kind": "project_derived",
                    "citation": ("technical-blank control distribution of GSE262419, derived in "
                                 "this audit; 3-MAD band from " +
                                 EV_IGLEWICZ["citation"]),
                    "url": "pert_data/GSE262419/source/GSE262419_Plate*.csv.gz "
                           "(Blank wells, Chemical_name == 'Blank')",
                    "retrieved": "2026-11-02",
                    "why_not_the_paper": ("the study names no detection floor per well; the "
                                          "blank wells are the only same-run no-cell reference"),
                },
                "reason_code": "derived_blank_separation_floor",
            },
        },
        "thresholds": {
            "gene_filter_min_mean": 10.0,
            "sample_expression_floor": 15.0,
            "vehicle_corr_floor": 0.85,
            "blank_ngenes_threshold": None,   # instantiated by --lock
            "blank_ngenes_median": None,
            "blank_ngenes_mad": None,
            "blank_mad_k": 3.0,
            "n_genes_retained": None,         # instantiated by --lock
        },
        "roles": {
            "treated_compound": "Solvent == DMSO and Chemical_name not in {VEH, ET-1}",
            "vehicle_control_dmso": "Chemical_name == VEH (Solvent DMSO)",
            "control_media": "Chemical_name == Media (is not a vehicle control)",
            "technical_blank": "Chemical_name == Blank (no cells; never a control)",
            "positive_control_et1": "Chemical_name == ET-1 at 0.01 uM (peptide positive control, "
                                    "not a small-molecule treatment)",
        },
        "control_rule": {
            "matched_control_key": "(plate, vehicle)",
            "text": ("a treated temp_env well is matched only to DMSO vehicle wells on the same "
                     "physical plate; Media is not a vehicle control, Blank is not a control"),
            "evidence": {
                **EV_TEMPENV_PAPER,
                "quote": ("the chemical effects were evaluated by comparing them to the "
                          "plate-specific DMSO controls"),
            },
        },
    },
    "tahoe": {
        "dataset_id": "tahoe",
        "observation_unit": ("(plate, physical well, cell_line) pseudobulk profile; 1,344 "
                             "physical wells -> 67,018 observations"),
        "rules": {
            "min_cells_per_profile": {
                "id": "H1",
                "statistic": "n_cells of the (plate, well, cell_line) profile",
                "threshold": 100,
                "op": "exclude observation if n_cells < 100",
                "evidence": {
                    **EV_ADR0029,
                    "source_floor": {
                        **EV_TAHOE_PREPRINT,
                        "quote": ("each cell line-drug treatment condition was required to have a "
                                  "minimum of 50 cells for inclusion in downstream analyses"),
                    },
                },
                "reason_code": "project_min_cells_per_profile",
            },
            "cell_line_excluded_by_source": {
                "id": "H2",
                "statistic": "cell line of the observation",
                "threshold": ["NCI-H2122", "NCI-H596", "NCI-H661"],
                "op": "exclude every observation of a cell line the source study removed",
                "evidence": {
                    **EV_TAHOE_PREPRINT,
                    "quote": ("three cell lines (NCI-H661, NCI-H596, NCI-H2122) that consistently "
                              "exhibited low counts across conditions were excluded from "
                              "subsequent analyses"),
                },
                "reason_code": "src_cell_line_excluded",
            },
            "mito_fraction_ceiling": {
                "id": "H3a",
                "statistic": "per-well mean mitochondrial read fraction",
                "threshold": 20.0,
                "op": "exclude well if mean_pcnt_mito > 20 (%)",
                "evidence": {
                    **EV_TAHOE_PREPRINT,
                    "quote": ("The full filters required cells to have at least 700 unique "
                              "molecular identifiers (UMIs), less than 20% mitochondrial reads, "
                              "a UMI z-score within ±3 ..., a mitochondrial percentage z-score "
                              "within ±3 ..., and at least 250 genes detected."),
                    "role": "verification gate: the atlas's own per-cell filters are already "
                            "applied upstream, so the per-well mean must sit inside them",
                },
                "reason_code": "src_mito_fraction_ceiling",
            },
            "min_genes_per_cell": {
                "id": "H3b",
                "statistic": "per-well mean detected genes per cell",
                "threshold": 250.0,
                "op": "exclude well if mean_gene_count < 250",
                "evidence": {
                    **EV_TAHOE_PREPRINT,
                    "quote": "and at least 250 genes detected.",
                    "role": "verification gate, as H3a",
                },
                "reason_code": "src_min_genes_per_cell",
            },
            "plate_depth_outlier": {
                "id": "H4",
                "statistic": ("per-plate modified z-score of log10(mean_tscp_count) over the "
                              "plate's 96 wells"),
                "threshold": 3.0,
                "op": "exclude well if |z| > 3",
                "evidence": {
                    **EV_IGLEWICZ,
                    "author_convention": ("the source study removes cells with a UMI z-score "
                                          "outside ±3 (quoted under H3a); this rule applies the "
                                          "same ±3 convention at well level with a MAD-based "
                                          "scale so the outlier cannot inflate its own reference"),
                },
                "reason_code": "derived_plate_depth_outlier",
            },
        },
        "thresholds": {
            "n_cells_min": 100,
            "source_n_cells_floor": 50,
            "excluded_cell_lines": ["NCI-H2122", "NCI-H596", "NCI-H661"],
            "mito_max_pcnt": 20.0,
            "min_mean_gene_count": 250.0,
            "well_z_abs_max": 3.0,
        },
        "roles": {
            "treated_compound": "drug != 'DMSO_TF' (exact equality; 'Trametinib (DMSO_TF solvate)' "
                                "is a treated drug)",
            "vehicle_control_dmso_tf": "drug == 'DMSO_TF'",
        },
        "control_rule": {
            "matched_control_key": "(plate, cell_line_id)",
            "text": ("a treated Tahoe observation is matched only to DMSO_TF profiles of the same "
                     "plate and the same cell line; a plate pools ~50 lines, so a cross-line mean "
                     "is the Simpson's-paradox mixture ADR 0029 rejects; the two control profiles "
                     "of one (plate, line) group are two wells, not two donors"),
            "evidence": {
                "kind": "project_derived",
                "citation": ("ADR 0031 D1 (control keying) and ADR 0029 decision 4 "
                             "(Delta stays inside one cell line)"),
                "url": "docs/adr/0031-tahoe-cm-schema-view-and-prnet-pretrain.md",
                "retrieved": "2026-11-02",
            },
        },
    },
}

# Reason codes are stable identifiers consumed by the handoff; this map is the
# human label written into the counts tables.
REASON_LABELS = {
    "src_sample_expression_floor": "T1 source sample expression floor (mean count/gene < 15)",
    "src_vehicle_reproducibility": "T2 source vehicle-control mean correlation (r < 0.85)",
    "derived_blank_separation_floor": "T3 technical-blank detection floor",
    "project_min_cells_per_profile": "H1 fewer cells than the 100-cell profile floor",
    "src_cell_line_excluded": "H2 cell line excluded by the source study",
    "src_mito_fraction_ceiling": "H3a per-well mitochondrial fraction above 20%",
    "src_min_genes_per_cell": "H3b per-well mean detected genes below 250",
    "derived_plate_depth_outlier": "H4 per-plate well depth outlier (|z| > 3)",
    "unresolved_missing_metric": "unresolved: a required QC metric is missing",
}

FLAG_LABELS = {
    "vehicle_corr_below_source_floor": "mean r vs vehicle pool below the source 0.85 floor "
                                       "(informational; not applied to non-vehicle wells)",
    "plate_libsize_robust_z_gt3": "per-plate log10(lib_size) modified z outside ±3",
    "plate_ngenes_robust_z_gt3": "per-plate log10(n_genes) modified z outside ±3",
    "n_cells_between_source_and_project_floor": "n_cells >= 50 (source floor) but < 100 "
                                                "(project floor): kept, sensitivity strata",
    "no_structure": "no resolvable structure: the dose encoder cannot consume this row "
                    "(input-layer exclusion, not a QC failure)",
    "well_mito_above_15pct": "per-well mean mitochondrial fraction above 15% "
                             "(inside the 20% source ceiling, reported as a flag)",
}


# --------------------------------------------------------------------------- #
# Pure statistics
# --------------------------------------------------------------------------- #

def modified_z(values) -> np.ndarray:
    """median/MAD modified z-score; a zero MAD yields NaN (never inf)."""
    x = np.asarray(values, dtype=float)
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad == 0:
        return np.full(x.shape, np.nan)
    return (x - med) / (1.4826 * mad)


def author_gene_filter(total_counts_per_gene, n_wells: int, min_mean: float = 10.0):
    """The study's gene filter: keep genes whose mean count across all wells is >= min_mean."""
    totals = np.asarray(total_counts_per_gene, dtype=float)
    return totals / float(n_wells) >= min_mean


def blank_separation(blank_n_genes, k: float = 3.0):
    """median + k·1.4826·MAD of the technical-blank detection counts."""
    x = np.asarray(blank_n_genes, dtype=float)
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med)))
    thr = med + k * 1.4826 * mad
    return med, mad, thr


# --------------------------------------------------------------------------- #
# Decision functions
# --------------------------------------------------------------------------- #

def _resolve(df: pd.DataFrame, required, reason_code="unresolved_missing_metric"):
    """Rows with a missing required metric are unresolved, not pass and not fail."""
    missing = pd.Series(False, index=df.index)
    for col in required:
        if col not in df:
            missing |= True
            continue
        s = df[col]
        missing |= s.isna() if s.dtype.kind in ("f", "c") else s.isin([None, ""])
    return missing.fillna(True)


def _finish(df, status, reasons, flags):
    out = pd.DataFrame({"qc_status": status,
                        "qc_reason": [";".join(r) for r in reasons],
                        "qc_flags": [";".join(f) for f in flags]}, index=df.index)
    return out


def _hit(mask, code, reasons, flags, into="reasons"):
    """Append a reason/flag code to every row where mask is True."""
    target = reasons if into == "reasons" else flags
    for i in np.flatnonzero(np.asarray(mask, dtype=bool)):
        target[i].append(code)


def evaluate_tempenv(df: pd.DataFrame, refs: dict) -> pd.DataFrame:
    """Apply T1/T2/T3 to a per-well temp_env frame.

    Required columns: role, lib_size_panel, n_genes_panel, meanr_vs_vehicle_pool.
    ``refs`` carries the instantiated thresholds (RULES['tempenv']['thresholds']).
    """
    n = len(df)
    reasons, flags = [[] for _ in range(n)], [[] for _ in range(n)]
    status = np.array(["pass"] * n, dtype=object)

    unresolved = _resolve(df, ["lib_size_panel", "n_genes_panel", "meanr_vs_vehicle_pool"])
    unresolved = unresolved.to_numpy()
    status[unresolved] = "unresolved"
    _hit(unresolved, "unresolved_missing_metric", reasons, flags)
    live = ~unresolved

    n_retained = float(refs["tempenv_n_genes_retained"])
    mean_per_gene = df["lib_size_panel"].to_numpy(dtype=float) / n_retained
    _hit(live & (mean_per_gene < refs["tempenv_sample_expression_floor"]),
         "src_sample_expression_floor", reasons, flags)

    is_vehicle = df["role"].to_numpy() == "vehicle_control_dmso"
    r = df["meanr_vs_vehicle_pool"].to_numpy(dtype=float)
    corr_floor = refs["tempenv_vehicle_corr_floor"]
    _hit(live & is_vehicle & (r < corr_floor), "src_vehicle_reproducibility", reasons, flags)

    # T3 asks "does this cell library look like a no-cell well?", so technical blanks
    # are its reference population and are never judged by it.
    not_blank = df["role"].to_numpy() != "technical_blank"
    _hit(live & not_blank & (df["n_genes_panel"].to_numpy(dtype=float)
                             <= refs["tempenv_blank_ngenes_threshold"]),
         "derived_blank_separation_floor", reasons, flags)

    _hit(live & (r < corr_floor), "vehicle_corr_below_source_floor", reasons, flags,
         into="flags")

    status[live & np.array([bool(x) for x in reasons])] = "fail"
    return _finish(df, status, reasons, flags)


def evaluate_tahoe(df: pd.DataFrame, refs: dict) -> pd.DataFrame:
    """Apply H1/H2/H3a/H3b/H4 to a per-observation Tahoe frame."""
    n = len(df)
    reasons, flags = [[] for _ in range(n)], [[] for _ in range(n)]
    status = np.array(["pass"] * n, dtype=object)

    unresolved = _resolve(df, ["n_cells", "mean_gene_count", "mean_pcnt_mito",
                               "mean_tscp_count", "cell_name"]).to_numpy()
    status[unresolved] = "unresolved"
    _hit(unresolved, "unresolved_missing_metric", reasons, flags)
    live = ~unresolved

    n_cells = df["n_cells"].to_numpy(dtype=float)
    n_min, src_floor = refs["tahoe_n_cells_min"], refs["tahoe_source_n_cells_floor"]
    _hit(live & (n_cells < n_min), "project_min_cells_per_profile", reasons, flags)
    _hit(live & (src_floor <= n_cells) & (n_cells < n_min),
         "n_cells_between_source_and_project_floor", reasons, flags, into="flags")

    excluded = set(refs["tahoe_excluded_cell_lines"])
    _hit(live & df["cell_name"].astype(str).isin(excluded).to_numpy(),
         "src_cell_line_excluded", reasons, flags)

    mito = df["mean_pcnt_mito"].to_numpy(dtype=float)
    mito_max = refs["tahoe_mito_max_pcnt"]
    _hit(live & (mito > mito_max), "src_mito_fraction_ceiling", reasons, flags)
    _hit(live & (mito <= mito_max) & (mito > 15.0), "well_mito_above_15pct", reasons, flags,
         into="flags")

    _hit(live & (df["mean_gene_count"].to_numpy(dtype=float)
                 < refs["tahoe_min_mean_gene_count"]),
         "src_min_genes_per_cell", reasons, flags)

    _hit(live & (np.abs(df["well_z_tscp"].to_numpy(dtype=float))
                 > refs["tahoe_well_z_abs_max"]),
         "derived_plate_depth_outlier", reasons, flags)

    if "has_structure" in df:
        _hit(live & ~df["has_structure"].to_numpy(dtype=bool), "no_structure", reasons, flags,
             into="flags")

    status[live & np.array([bool(x) for x in reasons])] = "fail"
    return _finish(df, status, reasons, flags)


# --------------------------------------------------------------------------- #
# Rule identity
# --------------------------------------------------------------------------- #

def rules_payload(rules: dict | None = None) -> dict:
    return {
        "rules_version": RULES_VERSION,
        "forbidden_bases": FORBIDDEN_BASES,
        "sources": rules if rules is not None else RULES,
    }


def canonical_rules_sha256(rules: dict | None = None) -> str:
    blob = json.dumps(rules_payload(rules), sort_keys=True, separators=(",", ":"),
                      ensure_ascii=False, default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()
