#!/usr/bin/env python
"""Build, freeze and verify the single two-source input formal training reads.

Stage 2 step 06.3 (`feat/pt/t03-input`, ticket #3).  See `../README.md` for the
contract, the counts and the reproduce commands; `contract.py` is the source of
truth for every frozen decision.

    # validate the frozen sources and (if present) the #1/#2 handoffs, no writes
    execution/build_input.py --check

    # the real build: writes <out-root>/06.3_input/output/
    execution/build_input.py --out-root /home/chang/arc/stage2/06_prnet_transfer

    # offline re-verification of a frozen build (hashes + every invariant)
    execution/build_input.py --verify --out-root /home/chang/arc/stage2/06_prnet_transfer

The build refuses to overwrite an existing frozen record without `--force` and
never writes to a frozen source.  `--stub-handoffs` exercises the whole pipeline
with synthetic include-all/pass-all handoffs (tests and smoke only); a stub
manifest is labelled `status: stub` and `--verify` rejects it unless
`--allow-stub`.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import consensus as consensus_mod
import folds as folds_mod
import writers
from contract import (
    Handoff as HandoffRecord, _frozen_manifest_row as _frozen_row,
    companion_handoff, read_json,
    ARTIFACTS, CM_BUNDLE_MANIFEST, CM_OBS, CM_VIEW, CONFIG_VERSION, CONTRACT_ID,
    ContractError, DATA_ROOT, ELIGIBILITY_DECISIONS, FOLD_SEED, MANIFEST_SCHEMA,
    N_FOLDS, N_GENES, STEP_DATA, TAHOE_OBS, TAHOE_VIEW, TAHOE_VIEW_MANIFEST,
    TIME_SOURCE_NOTE, TIME_SOURCE_TAHOE, TOOL_VERSION, TRACKED_RECORD,
    build_locked_config, canonical_json, config_digest, gene_measured_cm,
    gene_measured_tahoe, git_state, load_eligibility, load_gene_axis, load_qc,
    output_dir, relpath, resolve_handoff, runtime_env, sha256_file,
    source_records, split_combined_qc, tahoe_gene_space_ensembl,
    read_h5ad_header,
)

# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #

def log(message: str) -> None:
    print(f"[06.3_input] {message}", flush=True)


def _require_writer_env() -> dict:
    """The h5ad writer must be anndata >= 0.11 (0.12.16 is the pinned writer)."""
    env = runtime_env()
    major, minor = (int(x) for x in env["anndata"].split(".")[:2])
    if (major, minor) < (0, 11):
        raise ContractError(
            f"the builder writes dataframe/0.2.0 and needs anndata >= 0.11; "
            f"found {env['anndata']} in {env['executable']}. "
            "Run with /home/chang/miniconda3/envs/arc/bin/python.")
    return env


def _flatten_counts(prefix: str, value, out: dict[str, str]) -> None:
    if isinstance(value, dict):
        for key, sub in value.items():
            _flatten_counts(f"{prefix}.{key}" if prefix else str(key), sub, out)
    elif isinstance(value, (list, tuple)):
        out[prefix] = ";".join(str(v) for v in value)
    elif isinstance(value, bool):
        out[prefix] = "true" if value else "false"
    elif value is None:
        out[prefix] = ""
    else:
        out[prefix] = str(value)


def write_counts_csv(path: Path, counts: dict) -> None:
    flat: dict[str, str] = {}
    _flatten_counts("", counts, flat)
    table = pd.DataFrame({"key": list(flat.keys()), "value": list(flat.values())})
    table.to_csv(path, index=False)


def _utc_now() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


# --------------------------------------------------------------------------- #
# Handoffs (real or stub)
# --------------------------------------------------------------------------- #

def resolve_qc_path(kind: str) -> Path:
    """Per-source QC table if present, else #2's single combined table."""
    from contract import resolve_handoff
    try:
        return resolve_handoff(kind)
    except ContractError:
        return resolve_handoff("qc_combined")


def read_qc_handoff(path: Path | str, kind: str):
    """Normalise one per-source QC table; `path` may be #2's combined table."""
    from contract import Handoff, _read_table, normalize_qc_table, sha256_file, split_combined_qc
    from contract import relpath
    resolved = Path(path).expanduser().resolve()
    raw = _read_table(resolved)
    combined_column = next((c for c in ("source", "dataset_id", "dataset")
                            if c in raw.columns), None)
    note = f"per-observation QC for {kind}"
    if combined_column is not None:
        parts = split_combined_qc(raw, resolved)
        key = "cm" if kind == "cm" else "tahoe"
        part = parts.get(key)
        if part is None:
            raise ContractError(f"{resolved}: no {kind} rows in the combined QC table")
        raw = part.reset_index(drop=True)
        note = f"per-observation QC for {kind}, split from the combined table"
    table = normalize_qc_table(raw, resolved, kind)
    handoff = Handoff(f"{kind}_qc", resolved, sha256_file(resolved),
                      resolved.stat().st_size, len(table), note)
    return table, handoff


def qc_companion_checks(path: Path, kind: str, payload: dict) -> dict:
    """Extract and validate the #2 handoff companion facts for one source."""
    companion = companion_handoff(Path(path), ("handoff.json",))
    if companion is None:
        return {}
    doc = read_json(companion)
    out = {"handoff_json": str(companion), "sha256": sha256_file(companion)}
    key = "temp_env" if kind == "cm" else "tahoe"
    kept = doc.get("kept_observation_ids", {}).get(key)
    if kept is not None:
        out["kept_n_handoff"] = int(kept["n"])
    routing = doc.get("routing", {}).get(key, {})
    if routing:
        out["routing"] = routing
    mask = doc.get("gene_measurement_mask", {}).get("FAM153B", {}).get(key)
    if mask is not None:
        out["fam153b_measured_handoff"] = bool(mask.get("evidence", {}).get("measured"))
    return out


def load_handoffs_for_cm(args, cm_obs: pd.DataFrame):
    """Return (eligibility, cm_qc, records, stub) for the CM side."""
    if args.stub_handoffs:
        keys = sorted(set(cm_obs.loc[~cm_obs["is_control"].astype(bool),
                                     "drug_key"].astype(str)) - {""})
        eligibility = pd.DataFrame({
            "drug_key": keys,
            "eligibility_decision": "include",
            "elig_class": "stub",
        })
        cm_qc = pd.DataFrame({
            "well": cm_obs.index.astype(str),
            "qc_pass": True,
            "qc_reason": "stub: all observations pass",
        })
        return eligibility, cm_qc, [], True

    eligibility, eligibility_record = load_eligibility(args.eligibility)
    records = [eligibility_record]

    if args.qc:
        cm_qc, record = read_qc_handoff(args.qc, "cm")
    elif args.cm_qc:
        cm_qc, record = read_qc_handoff(args.cm_qc, "cm")
    else:
        cm_qc, record = read_qc_handoff(resolve_qc_path("cm_qc"), "cm")
    records.append(record)
    # #1's handoff.json carries the adjudicated class counts; reconcile them so a
    # silently regenerated table cannot pass as the frozen audit.
    companion = companion_handoff(record.path, ("handoff.json",))
    if companion is not None:
        doc = read_json(companion)
        counts_doc = doc.get("counts", {})
        observed = eligibility["eligibility_decision"].value_counts().to_dict()
        for decision in ELIGIBILITY_DECISIONS:
            expected = counts_doc.get(decision)
            if expected is not None and int(expected) != int(observed.get(decision, 0)):
                raise ContractError(
                    f"#1 handoff.json says {expected} {decision} keys, "
                    f"the table has {observed.get(decision, 0)}")
        records.append(HandoffRecord("eligibility_handoff", companion,
                                     sha256_file(companion),
                                     companion.stat().st_size,
                                     int(len(eligibility)),
                                     "eligibility handoff companion (counts reconciled)"))
    return eligibility, cm_qc, records, False


def load_tahoe_qc(args, tahoe_index, *, stub: bool):
    """Return (qc table, handoff record or None) for the Tahoe side."""
    if stub:
        return pd.DataFrame({
            "well": np.asarray(tahoe_index, dtype=str),
            "qc_pass": True,
            "qc_reason": "stub: all observations pass",
        }), None
    if args.qc:
        return read_qc_handoff(args.qc, "tahoe")
    if args.tahoe_qc:
        return read_qc_handoff(args.tahoe_qc, "tahoe")
    return read_qc_handoff(resolve_qc_path("tahoe_qc"), "tahoe")


# --------------------------------------------------------------------------- #
# CM side
# --------------------------------------------------------------------------- #

def load_cm_view(need_matrix: bool = True) -> dict:
    """The frozen CM view: obs/var/uns always, X + log layer only on demand.

    `need_matrix=False` is the planning/inspection path (`--check`, the build
    preamble); it decodes the metadata with h5py so a fail-closed handoff check
    does not pay for the 289 MB view, let alone its 650 MB of dense matrices.
    """
    if not CM_VIEW.exists():
        raise ContractError(f"missing frozen CM view: {CM_VIEW}")
    if not need_matrix:
        header = read_h5ad_header(CM_VIEW)
        obs, var, uns = header.obs, header.var, header.uns
        out: dict = {}
    else:
        import anndata as ad
        adata = ad.read_h5ad(CM_VIEW)
        try:
            obs, var, uns = adata.obs.copy(), adata.var.copy(), dict(adata.uns)
            out = {"X": np.asarray(adata.X),
                   "log1p": np.asarray(adata.layers["log1p_cp10k"])}
        finally:
            if adata.file is not None:
                adata.file.close()
        if out["X"].dtype != np.int32 or out["log1p"].dtype != np.float32:
            raise ContractError(
                f"{CM_VIEW}: X {out['X'].dtype} layer {out['log1p'].dtype}, "
                "expected int32/float32")
    if list(var.index.astype(str)) != list(load_gene_axis().index):
        raise ContractError(f"{CM_VIEW}: var is not the frozen gene axis")
    out.update({"obs": obs, "var": var, "uns": uns})
    return out


def build_cm(handoffs, gene_axis, *, write: bool, out: Path) -> dict:
    """Read the frozen CM view, attach the flags and build folds + consensus."""
    from contract import GSE262419_HASH
    view = load_cm_view(need_matrix=True)
    hash_table = pd.read_csv(GSE262419_HASH, dtype={"Well_ID": str})

    obs = writers.assemble_cm_obs(view["obs"], hash_table,
                                  handoffs["eligibility"], handoffs["cm_qc"])
    keep = writers.cm_keep_mask(obs)
    counts = {
        "view_rows": int(len(obs)),
        "view_treated": int((~obs["is_control"].astype(bool)).sum()),
        "view_controls": int(obs["is_control"].astype(bool).sum()),
        "qc_pass_rows": int(obs["qc_pass"].astype(bool).sum()),
        "qc_fail_rows": int((~obs["qc_pass"].astype(bool)).sum()),
        "qc_fail_treated_view": int(((~obs["is_control"].astype(bool))
                                     & (~obs["qc_pass"].astype(bool))).sum()),
        "qc_fail_controls_view": int((obs["is_control"].astype(bool)
                                      & (~obs["qc_pass"].astype(bool))).sum()),
        "dropped_qc_fail_treated": int(((~obs["is_control"].astype(bool))
                                        & (~obs["qc_pass"].astype(bool))).sum()),
        "dropped_qc_fail_controls": int((obs["is_control"].astype(bool)
                                         & (~obs["qc_pass"].astype(bool))).sum()),
        "kept_rows": int(keep.sum()),
        "kept_treated": int((keep & ~obs["is_control"].astype(bool)).sum()),
        "kept_controls": int((keep & obs["is_control"].astype(bool)).sum()),
        "kept_dmso_controls": int((keep & obs["is_control"].astype(bool)
                                   & (obs["vehicle"].astype(str) == "DMSO")).sum()),
        "kept_media_controls": int((keep & obs["is_control"].astype(bool)
                                    & (obs["vehicle"].astype(str) == "Media")).sum()),
        "drug_keys_in_view": int(obs.loc[~obs["is_control"].astype(bool),
                                         "drug_key"].astype(str).nunique()),
    }
    decision = obs["eligibility_decision"].astype(str)
    counts["eligibility_include_drugs"] = int(
        obs.loc[~obs["is_control"].astype(bool) & (decision == "include"),
                "drug_key"].astype(str).nunique())
    counts["eligibility_exclude_drugs"] = int(
        obs.loc[~obs["is_control"].astype(bool) & (decision == "exclude"),
                "drug_key"].astype(str).nunique())
    counts["eligibility_uncertain_drugs"] = int(
        obs.loc[~obs["is_control"].astype(bool) & (decision == "uncertain"),
                "drug_key"].astype(str).nunique())

    kept = obs.loc[keep].copy()
    X = view["X"][keep.to_numpy()]
    log1p = view["log1p"][keep.to_numpy()]

    main = kept.loc[kept["main_set"].astype(bool)]
    main_drugs = (main.assign(drug_key=main["drug_key"].astype(str),
                              plate=main["plate"].astype(str),
                              dose=main["dose_um"].astype(float))
                  .groupby(["drug_key", "plate"], observed=True)
                  .agg(n_wells=("well", "size"), n_doses=("dose_um", "nunique"))
                  .reset_index())
    include_keys = set(kept.loc[kept["eligibility_decision"].astype(str) == "include",
                               "drug_key"].astype(str))
    counts["main_set_drugs"] = int(main_drugs["drug_key"].nunique())
    counts["main_set_wells"] = int(len(main))
    counts["main_set_plates"] = int(main_drugs["plate"].nunique())
    counts["main_set_without_matched_control"] = int(
        (~main["matched_control_ok"].astype(bool)).sum())
    counts["matched_control_source"] = (
        "handoff" if "matched_control_ok" in handoffs["cm_qc"].columns else "derived")
    counts["qc_pass_vehicle"] = {
        "dmso": int((kept["is_control"].astype(bool)
                     & (kept["vehicle"].astype(str) == "DMSO")).sum()),
        "media": int((kept["is_control"].astype(bool)
                      & (kept["vehicle"].astype(str) == "Media")).sum()),
    }
    counts["eligibility_include_without_usable_well"] = int(
        len(include_keys - set(main_drugs["drug_key"])))
    expansion = kept.loc[kept["expansion_set"].astype(bool)]
    counts["expansion_set_wells"] = int(len(expansion))
    counts["expansion_set_drugs"] = int(expansion["drug_key"].astype(str).nunique())
    counts["expanded_train_drugs"] = int(
        len(set(main_drugs["drug_key"].astype(str))
            | set(expansion["drug_key"].astype(str))))
    counts["main_set_doses"] = sorted(float(d) for d in main["dose_um"].unique())

    # Representative structure per drug (all rows of a drug share the parent key).
    rep = (kept.loc[kept["main_set"].astype(bool),
                    ["drug_key", "inchikey", "compound_name", "smiles"]]
           .assign(drug_key=lambda d: d["drug_key"].astype(str))
           .sort_values("drug_key", kind="mergesort")
           .drop_duplicates("drug_key"))
    drugs = main_drugs.groupby("drug_key", observed=True)["n_wells"].sum().reset_index()
    drugs = drugs.merge(rep[["drug_key", "inchikey"]], on="drug_key", how="left",
                        validate="one_to_one")
    drugs["inchikey"] = drugs["inchikey"].fillna("").astype(str)

    assign = folds_mod.cluster_drugs(drugs, n_folds=N_FOLDS, seed=FOLD_SEED)
    fold_long = folds_mod.expand_fold_table(assign, main_drugs, n_folds=N_FOLDS)
    fold_table = folds_mod.collapse_fold_table(fold_long)
    fold_facts = folds_mod.assert_fold_invariants(
        fold_table, assign, expected_drugs=set(main_drugs["drug_key"].astype(str)),
        n_folds=N_FOLDS)
    sizes = folds_mod.fold_sizes(fold_long)
    counts["folds"] = {
        "n_folds": N_FOLDS,
        "seed": FOLD_SEED,
        "test_drugs": [int(v) for v in sizes["test_drugs"]],
        "valid_drugs": [int(v) for v in sizes["valid_drugs"]],
        "train_drugs": [int(v) for v in sizes["train_drugs"]],
        "test_wells": [int(v) for v in sizes["test_wells"]],
        "valid_wells": [int(v) for v in sizes["valid_wells"]],
        "train_wells": [int(v) for v in sizes["train_wells"]],
        "max_drugs_per_fold": int(sizes[["test_drugs", "valid_drugs"]].to_numpy().max()),
        "n_fold_rows": int(len(fold_table)),
    }

    c_star, membership, summary, plate_names, plate_weights = (
        consensus_mod.build_consensus(log1p, kept, fold_long,
                                      n_folds=N_FOLDS, n_genes=N_GENES))
    provenance = consensus_mod.assert_consensus_provenance(
        membership, kept, fold_long, n_folds=N_FOLDS)
    counts["consensus"] = {
        "fold_n_plates": [int(v) for v in summary["n_plates"]],
        "fold_n_wells": [int(v) for v in summary["n_wells"]],
        "members_total": int(len(membership)),
        "per_fold_plates": provenance["per_fold"],
    }

    result = {
        "obs": kept, "X": X, "log1p": log1p, "folds": fold_table,
        "folds_long": fold_long,
        "assign": assign, "sizes": sizes, "c_star": c_star,
        "membership": membership, "consensus_summary": summary,
        "plate_names": plate_names, "plate_weights": plate_weights,
        "counts": counts, "fold_facts": fold_facts,
    }
    return result


# --------------------------------------------------------------------------- #
# Tahoe side
# --------------------------------------------------------------------------- #

def tahoe_gene_nnz(path: Path, n_genes: int = N_GENES,
                   chunk: int = 20_000_000) -> np.ndarray:
    """Per-gene non-zero counts from the CSR indices (X is never densified)."""
    import h5py
    counts = np.zeros(n_genes, dtype=np.int64)
    with h5py.File(path, "r") as fh:
        indptr = np.asarray(fh["X/indptr"][:], dtype=np.int64)
        n_rows = len(indptr) - 1
        if n_rows <= 0:
            raise ContractError(f"{path}: empty CSR matrix")
        if indptr[-1] != fh["X/indices"].shape[0]:
            raise ContractError(f"{path}: CSR indptr/nnz disagree")
        indices = fh["X/indices"]
        for start in range(0, indptr[-1], chunk):
            stop = min(start + chunk, indptr[-1])
            block = np.asarray(indices[start:stop], dtype=np.int64)
            counts += np.bincount(block, minlength=n_genes)
    return counts


def build_tahoe(qc: pd.DataFrame, gene_axis, *, tahoe_rows: int | None,
                tahoe_measured: pd.Series) -> dict:
    header = read_h5ad_header(TAHOE_VIEW)
    obs_all = header.obs
    var = header.var
    if list(var.index.astype(str)) != list(gene_axis.index):
        raise ContractError(f"{TAHOE_VIEW}: var is not the frozen gene axis")

    if tahoe_rows is not None:
        obs_all = obs_all.iloc[: int(tahoe_rows)].copy()
    keep_ids = set(obs_all.index.astype(str))
    qc = qc.loc[qc["well"].astype(str).isin(keep_ids)].reset_index(drop=True)
    missing = keep_ids - set(qc["well"].astype(str))
    if missing:
        raise ContractError(
            f"QC handoff does not cover {len(missing)} Tahoe rows, "
            f"e.g. {sorted(missing)[:5]}")
    obs = writers.assemble_tahoe_obs(obs_all, qc)
    control_groups = writers.tahoe_control_groups(obs)
    pairs = writers.tahoe_pairs(obs, control_groups)

    counts = {
        "rows": int(len(obs)),
        "treated": int((~obs["is_control"].astype(bool)).sum()),
        "controls": int(obs["is_control"].astype(bool).sum()),
        "qc_fail_rows": int((~obs["qc_pass"].astype(bool)).sum()),
        "qc_fail_treated": int(((~obs["is_control"].astype(bool))
                                & (~obs["qc_pass"].astype(bool))).sum()),
        "qc_fail_controls": int((obs["is_control"].astype(bool)
                                 & (~obs["qc_pass"].astype(bool))).sum()),
        "has_structure_false": int((~obs["has_structure"].astype(bool)).sum()),
        "usable_treated": int(obs["tahoe_usable"].astype(bool).sum()),
        "usable_train": int((obs["tahoe_usable"].astype(bool)
                             & (obs["tahoe_split"].astype(str) == "train")).sum()),
        "usable_valid": int((obs["tahoe_usable"].astype(bool)
                             & (obs["tahoe_split"].astype(str) == "valid")).sum()),
        "control_groups": int(control_groups["group_key"].nunique()),
        "groups_with_zero_qc_pass": int(
            obs.loc[obs["is_control"].astype(bool), "control_group"].nunique()
            - control_groups["group_key"].nunique()),
        "control_members": int(len(control_groups)),
        "pairs": int(len(pairs["row"])),
        "pairs_skipped_no_control": int(pairs["skipped_no_control"]),
        "usable_without_pair": int(
            int(obs["tahoe_usable"].astype(bool).sum())
            - int(np.asarray(pairs["usable"], dtype=bool).sum())),
        "split_train_rows": int((obs["tahoe_split"].astype(str) == "train").sum()),
        "split_valid_rows": int((obs["tahoe_split"].astype(str) == "valid").sum()),
        "cell_lines": int(obs["cell_line_id"].astype(str).nunique()),
        "plates": int(obs["plate"].astype(str).nunique()),
        "drug_keys": int(len(set(obs.loc[~obs["is_control"].astype(bool),
                                         "drug_key"].astype(str)) - {""})),
        "valid_conditions": int(obs.loc[(obs["tahoe_split"].astype(str) == "valid")
                                        & (~obs["is_control"].astype(bool)),
                                        "tahoe_valid_group"].astype(str).nunique()),
        "time_source": TIME_SOURCE_TAHOE,
        "matched_control_source": (
            "handoff" if "matched_control_ok" in qc.columns else "derived"),
    }
    nnz = tahoe_gene_nnz(TAHOE_VIEW)
    counts["all_zero_measured_genes"] = int(
        ((nnz == 0) & tahoe_measured.to_numpy()).sum())
    counts["unmeasured_genes"] = int((~tahoe_measured).sum())
    return {"obs": obs, "var": var, "control_groups": control_groups,
            "pairs": pairs, "nnz": nnz, "counts": counts}


# --------------------------------------------------------------------------- #
# Assemble + write
# --------------------------------------------------------------------------- #

def _cm_uns(args, records, handoff_records) -> dict:
    from contract import (CM_REPLICATE_PROVENANCE, CM_TIME_PROVENANCE)
    return {
        "schema": "06.3_input/1",
        "contract_id": CONTRACT_ID,
        "tool_version": TOOL_VERSION,
        "gene_axis": "stage2/02_cm_data/output/corpus/v1/genes.csv (order preserved)",
        "column_provenance": {
            "replicate": CM_REPLICATE_PROVENANCE,
            "replicate_raw": "pert_data/GSE262419/source/GSE262419_hash.csv.gz "
                             "Chemical_Index/Replicate_Number/Solvent joined by "
                             "(Plate_ID, normalised Well_ID)",
            "time_h": CM_TIME_PROVENANCE,
            "vehicle": "inherited loader value; equals the raw Solvent on every retained row",
            "control_group": "inherited temp_env::Plate<N>; MIXES DMSO and Media rows — "
                             "never use it as a matched baseline",
            "main_set": "eligibility_decision == include AND qc_pass (treated rows only)",
            "matched_control_ok": "carried from ticket #2's per-observation QC (same-plate "
                                  "DMSO control availability); a control row is its own "
                                  "matched control.  Derived when the handoff omits it.",
            "matched_control_key": "ticket #2's matched-control group key (plate x vehicle)",
        },
        "control_cache": "NO uns['control_mean'] is carried: the inherited 02 cache is a "
                         "DMSO+Media raw-count mean and is superseded. The correctly "
                         "scaled per-fold consensus c* is in consensus_controls.npz with "
                         "row membership in consensus_membership.parquet. X keeps raw "
                         "int32 counts, so any raw-count mean is reproducible.",
        # uns must stay h5ad-writable: flat strings, never lists of dicts
        "handoffs": [f"{r.kind}|{relpath(r.path)}|{r.sha256}" for r in handoff_records],
        "sources": [f"{r.role}|{relpath(r.path)}|{r.content_sha256 or 'no-content-hash'}"
                    for r in records],
    }


def _tahoe_uns(args, records, handoff_records) -> dict:
    return {
        "schema": "06.3_input/1",
        "contract_id": CONTRACT_ID,
        "tool_version": TOOL_VERSION,
        "gene_axis": "stage2/02_cm_data/output/corpus/v1/genes.csv (order preserved)",
        "column_provenance": {
            "time_h": "NaN on every row: no atlas table carries a time column "
                      "(ADR 0029 update 2026-11-01, ADR 0031 D3). " + TIME_SOURCE_NOTE,
            "time_source": TIME_SOURCE_TAHOE,
            "physical_well": "the smp_* token of the composite obs.well "
                             "(tahoe::<plate>:<smp>:<cell_line>); the physical well is "
                             "not unique on its own",
            "observation_id": "obs index = Tahoe observation = plate x physical well x "
                              "cell line; never treat it as an independent culture well",
            "replicate": "generated within (control_group, parent_key, dose) index, not an "
                         "atlas replicate id",
            "vehicle": "'' on every treated row; 'DMSO_TF' on controls; no solvent column exists",
            "has_structure": "False on 2 drugs / 300 treated rows: kept for auditability, "
                             "excluded from training and evaluation via tahoe_usable",
            "matched_control_ok": "carried from ticket #2's per-observation QC (same "
                                  "(plate, cell_line) control availability); a control row "
                                  "is its own matched control.",
            "matched_control_key": "ticket #2's matched-control group key (plate x cell line)",
            "tahoe_usable": "~is_control AND has_structure AND qc_pass AND "
                            "matched_control_ok",
            "tahoe_split": "train | valid by plate13/plate14 (existing technical default)",
            "tahoe_valid_group": "drug_key x dose_uM x cell_line_id (CM-style grouping for "
                                 "direction AUPRC@200; controls keyed by plate x cell line)",
        },
        "control_cache": "NO uns['control_mean'] is carried: the source cache is a raw-count "
                         "mean while training targets are log1p_cp10k (the ADR 0031 recorded "
                         "deviation). Matched controls are computed at load time as the mean "
                         "of the control group's log1p_cp10k rows, using "
                         "tahoe_control_groups.parquet (700 groups x 2 rows) and "
                         "tahoe_pairs.npz; 67k x 18k control vectors are never materialised.",
        "handoffs": [f"{r.kind}|{relpath(r.path)}|{r.sha256}" for r in handoff_records],
        "sources": [f"{r.role}|{relpath(r.path)}|{r.content_sha256 or 'no-content-hash'}"
                    for r in records],
    }


def build(args) -> dict:
    env = _require_writer_env()
    out = output_dir(args.out_root)
    if (out / ARTIFACTS["manifest"]).exists() and not args.force:
        raise ContractError(
            f"{out} already holds a frozen manifest; refusing to overwrite it "
            "without --force")

    gene_axis = load_gene_axis()
    records = source_records()
    if [r.role for r in records][:3] != ["cm_view", "tahoe_view", "gene_axis"]:
        raise ContractError("source record order changed unexpectedly")
    log(f"gene axis: {N_GENES} symbols, sha256 {sha256_file(records[2].path)[:12]}…")
    for record in records:
        log(f"source {record.role}: {relpath(record.path)} ({record.bytes:,} B)")

    # ---- CM side -----------------------------------------------------------
    log("reading the frozen CM view (4,478 x 18,164 dense)…")
    cm_preview = load_cm_view(need_matrix=False)
    eligibility, cm_qc, handoff_records, stub = load_handoffs_for_cm(
        args, cm_preview["obs"])
    log(f"handoffs: {len(eligibility)} eligibility rows, {len(cm_qc)} CM QC rows"
        + (", STUB mode (not a frozen build)" if stub else ""))
    del cm_preview
    # Only create the artifact directory once every handoff has resolved: a
    # fail-closed build must leave no half-initialised output tree behind.
    out.mkdir(parents=True, exist_ok=True)
    handoff_checks: dict = {"stub": bool(stub)}
    handoff_checks["fold_seed"] = FOLD_SEED
    if not stub:
        qc_record = next(r for r in handoff_records if r.kind == "cm_qc")
        handoff_checks["cm_qc"] = qc_companion_checks(qc_record.path, "cm", {})

    cm = build_cm({"eligibility": eligibility, "cm_qc": cm_qc},
                  gene_axis, write=True, out=out)
    log(f"CM: kept {cm['counts']['kept_rows']} rows "
        f"({cm['counts']['kept_treated']} treated / {cm['counts']['kept_controls']} control), "
        f"main_set {cm['counts']['main_set_drugs']} drugs / "
        f"{cm['counts']['main_set_wells']} wells")

    cm_measured = gene_measured_cm(gene_axis)
    cm_nnz = (cm["X"] != 0).sum(axis=0)
    cm_var = writers.var_frame(gene_axis, cm_measured, cm_nnz)

    # Reconcile the artifact population against #2's own counts: the CM
    # artifact must be exactly the QC-pass rows of the frozen view.
    cm_kept_pass = int(cm["obs"]["qc_pass"].astype(bool).sum())
    if cm_kept_pass != len(cm["obs"]):
        raise ContractError("a kept CM row is not QC-pass")
    if cm_kept_pass != int(cm["counts"]["qc_pass_rows"]):
        raise ContractError(
            f"kept CM rows {cm_kept_pass} != the view's QC-pass rows "
            f"{cm['counts']['qc_pass_rows']}")
    handoff_checks["cm_kept_rows"] = cm_kept_pass
    handoff_checks["cm_qc_pass_rows_in_view"] = int(cm["counts"]["qc_pass_rows"])
    handoff_checks["cm_qc_fail_rows_in_view"] = int(cm["counts"]["qc_fail_rows"])
    handoff_checks["cm_handoff_pass_rows_all_wells"] = int(cm_qc["qc_pass"].sum())
    fam = cm_measured.loc[contract_fam153b()]
    if "cm_qc" in handoff_checks:
        expected = handoff_checks["cm_qc"].get("fam153b_measured_handoff")
        if expected is not None and bool(expected) != bool(fam):
            raise ContractError(
                "#2's handoff.json says FAM153B is measured in temp_env but the "
                "gene-panel derivation disagrees")

    # ---- Tahoe side --------------------------------------------------------
    tahoe = None
    tahoe_measured = None
    tahoe_var = None
    if args.skip_tahoe:
        log("Tahoe side skipped (--skip-tahoe)")
    else:
        log("reading the frozen Tahoe view (h5py, obs only — layers never loaded)…")
        tahoe_index = read_h5ad_header(TAHOE_VIEW, with_uns=False).obs.index.astype(str).to_numpy()
        if args.tahoe_rows is not None:
            tahoe_index = tahoe_index[: int(args.tahoe_rows)]
        tahoe_qc, tahoe_qc_record = load_tahoe_qc(args, tahoe_index, stub=stub)
        if tahoe_qc_record is not None:
            handoff_records = handoff_records + [tahoe_qc_record]
        ens_space = tahoe_gene_space_ensembl()
        tahoe_measured = gene_measured_tahoe(gene_axis, ens_space)
        tahoe = build_tahoe(tahoe_qc, gene_axis, tahoe_rows=args.tahoe_rows,
                            tahoe_measured=tahoe_measured)
        if not stub:
            tahoe_qc_record = next(
                (r for r in handoff_records if r.kind == "tahoe_qc"), None)
            if tahoe_qc_record is not None:
                handoff_checks["tahoe_qc"] = qc_companion_checks(
                    tahoe_qc_record.path, "tahoe", {})
                kept_n = handoff_checks["tahoe_qc"].get("kept_n_handoff")
                observed = int(tahoe["obs"]["qc_pass"].astype(bool).sum())
                if kept_n is not None and int(kept_n) != observed:
                    raise ContractError(
                        f"#2 handoff.json says {kept_n} kept Tahoe observations, "
                        f"the view yields {observed}")
            fam = tahoe_measured.loc[contract_fam153b()]
            expected = handoff_checks.get("tahoe_qc", {}).get(
                "fam153b_measured_handoff")
            if expected is not None and bool(expected) != bool(fam):
                raise ContractError(
                    "#2's handoff.json says FAM153B is missing from Tahoe but the "
                    "gene-space derivation disagrees")
        log(f"Tahoe: {tahoe['counts']['rows']} rows, "
            f"{tahoe['counts']['usable_treated']} usable treated "
            f"({tahoe['counts']['usable_train']} train / "
            f"{tahoe['counts']['usable_valid']} valid)")
        tahoe_var = writers.var_frame(gene_axis, tahoe_measured, tahoe["nnz"])

    # ---- Write the small artifacts ----------------------------------------
    folds_path = out / ARTIFACTS["folds"]
    cm["folds"].to_parquet(folds_path, index=False)
    cm["sizes"].to_csv(out / ARTIFACTS["fold_sizes"], index=False)
    cm["consensus_summary"].to_csv(out / ARTIFACTS["consensus_summary"], index=False)
    cm["membership"].to_parquet(out / ARTIFACTS["consensus_membership"], index=False)
    np.savez_compressed(
        out / ARTIFACTS["consensus_controls"],
        schema=np.array("06.3_input/consensus/1"),
        c_star=cm["c_star"],
        folds=np.arange(N_FOLDS, dtype=np.int16),
        fold_ids=np.arange(N_FOLDS, dtype=np.int16),
        plate_names=np.asarray(cm["plate_names"], dtype=str),
        plate_weights=cm["plate_weights"],
        n_plates=cm["consensus_summary"]["n_plates"].to_numpy(dtype=np.int16),
        n_wells=cm["consensus_summary"]["n_wells"].to_numpy(dtype=np.int32),
        n_genes=np.array(N_GENES, dtype=np.int64),
    )
    if tahoe is not None:
        tahoe["control_groups"].to_parquet(out / ARTIFACTS["tahoe_control_groups"],
                                          index=False)
        np.savez_compressed(
            out / ARTIFACTS["tahoe_pairs"],
            schema=np.array("06.3_input/tahoe_pairs/1"),
            row=tahoe["pairs"]["row"],
            group_index=tahoe["pairs"]["group_index"],
            usable=tahoe["pairs"]["usable"],
            split=tahoe["pairs"]["split"].astype(str),
        )

    var_checks = pd.DataFrame({
        "symbol": gene_axis.index.astype(str),
        "ensembl_id": gene_axis["ensembl_id"].astype(str).to_numpy(),
        "match_type": gene_axis["match_type"].astype(str).to_numpy(),
        "cm_measured": cm_measured.to_numpy(dtype=bool),
        "tahoe_measured": tahoe_measured.to_numpy(dtype=bool) if tahoe else False,
        "cm_n_nonzero": cm_nnz.astype(np.int64),
        "tahoe_n_nonzero": tahoe["nnz"].astype(np.int64) if tahoe else -1,
    })
    var_checks.to_csv(out / ARTIFACTS["var_checks"], index=False)

    locked = build_locked_config()
    (out / ARTIFACTS["locked_config"]).write_text(canonical_json(locked), encoding="utf-8")

    # ---- h5ad artifacts ----------------------------------------------------
    log("writing cm_transfer.h5ad…")
    writers.write_cm_h5ad(out / ARTIFACTS["cm"], cm["obs"], cm["X"], cm["log1p"],
                          cm_var, _cm_uns(args, records, handoff_records))
    if tahoe is not None:
        log("writing tahoe_transfer.h5ad (streaming CSR X + layer)…")
        writers.write_tahoe_h5ad(TAHOE_VIEW, out / ARTIFACTS["tahoe"],
                                 tahoe["obs"], tahoe_var,
                                 _tahoe_uns(args, records, handoff_records),
                                 row_limit=None)
        facts = writers.verify_tahoe_h5ad(out / ARTIFACTS["tahoe"],
                                          expected_rows=None)
        log(f"tahoe_transfer.h5ad: {facts['rows']} x {facts['cols']}, "
            f"nnz {facts['nnz']:,}")

    # ---- counts, lock and manifest ----------------------------------------
    counts = {"genes": N_GENES,
              "tahoe_unmeasured_genes": int((~tahoe_measured).sum())
              if tahoe_measured is not None else None,
              "cm": cm["counts"],
              "folds": cm["counts"]["folds"],
              "consensus": cm["counts"]["consensus"],
              "tahoe": tahoe["counts"] if tahoe else {"skipped": True}}
    write_counts_csv(out / ARTIFACTS["counts"], counts)
    log(f"counts written ({sum(1 for _ in (out / ARTIFACTS['counts']).read_text().splitlines())} rows)")

    write_contract_json(out, cm, tahoe, cm_var, tahoe_var, counts, handoff_checks)
    config_path = out / ARTIFACTS["locked_config"]
    artifacts = []
    for name in [ARTIFACTS["cm"], ARTIFACTS["tahoe"], ARTIFACTS["folds"],
                 ARTIFACTS["consensus_controls"], ARTIFACTS["consensus_membership"],
                 ARTIFACTS["tahoe_control_groups"], ARTIFACTS["tahoe_pairs"],
                 ARTIFACTS["locked_config"], ARTIFACTS["contract"],
                 ARTIFACTS["counts"], ARTIFACTS["fold_sizes"],
                 ARTIFACTS["consensus_summary"], ARTIFACTS["var_checks"]]:
        path = out / name
        if not path.exists():
            continue
        artifacts.append({
            "path": name,
            "bytes": int(path.stat().st_size),
            "sha256": sha256_file(path),
            "tracked": name in TRACKED_RECORD,
            "role": _artifact_role(name),
        })
        log(f"hashed {name} ({path.stat().st_size:,} B)")

    manifest = {
        "schema": MANIFEST_SCHEMA,
        "tool_version": TOOL_VERSION,
        "contract_id": CONTRACT_ID,
        "config_version": CONFIG_VERSION,
        "status": "stub" if (stub or args.tahoe_rows) else "frozen",
        "stub_handoffs": bool(stub),
        "tahoe_rows_limited": args.tahoe_rows,
        "built_at": _utc_now(),
        "git": git_state(),
        "env": env,
        "out_dir": relpath(out),
        "data_root": str(DATA_ROOT),
        "build": {
            "fold_seed": FOLD_SEED,
            "n_folds": N_FOLDS,
            "tahoe_valid_plates": ["plate13", "plate14"],
            "argv": sys.argv[1:],
        },
        "config": {"path": ARTIFACTS["locked_config"],
                   "sha256": sha256_file(config_path),
                   "digest": config_digest(build_locked_config())},
        "sources": [record.as_dict() for record in records],
        "source_manifests": [
            {"path": relpath(CM_BUNDLE_MANIFEST),
             "sha256": sha256_file(CM_BUNDLE_MANIFEST)},
            {"path": relpath(TAHOE_VIEW_MANIFEST),
             "sha256": sha256_file(TAHOE_VIEW_MANIFEST)},
        ],
        "handoffs": [record.as_dict() for record in handoff_records],
        "handoff_checks": handoff_checks,
        "artifacts": artifacts,
        "counts": counts,
        "gates": {
            "gene_axis_equal": True,
            "fold_invariants": True,
            "consensus_provenance": True,
            "handoff_coverage": True,
            "frozen_sources_unmodified": True,
        },
    }
    (out / ARTIFACTS["manifest"]).write_text(canonical_json(manifest), encoding="utf-8")
    log(f"manifest written: {relpath(out)}  status={manifest['status']}")
    return manifest


def contract_fam153b() -> str:
    """The one CM axis slot with no Tahoe measurement (ADR 0031 D2 / #2 mask)."""
    from contract import FAM153B_SYMBOL
    return FAM153B_SYMBOL


def _artifact_role(name: str) -> str:
    return {
        ARTIFACTS["cm"]: "cm frozen input",
        ARTIFACTS["tahoe"]: "tahoe frozen input",
        ARTIFACTS["folds"]: "five drug-key folds",
        ARTIFACTS["consensus_controls"]: "per-fold consensus c*",
        ARTIFACTS["consensus_membership"]: "consensus control row membership",
        ARTIFACTS["tahoe_control_groups"]: "Tahoe matched-control membership",
        ARTIFACTS["tahoe_pairs"]: "Tahoe treated row -> control group",
        ARTIFACTS["locked_config"]: "locked experiment config",
        ARTIFACTS["counts"]: "frozen counts",
        ARTIFACTS["fold_sizes"]: "per fold/role drug and well counts",
        ARTIFACTS["consensus_summary"]: "per-fold consensus plate pool",
        ARTIFACTS["var_checks"]: "per-gene measurement audit",
        ARTIFACTS["contract"]: "machine-readable interface contract for downstream tickets",
    }[name]


# --------------------------------------------------------------------------- #
# Verify (offline)
# --------------------------------------------------------------------------- #

def verify(args) -> dict:
    import anndata as ad

    out = output_dir(args.out_root)
    manifest_path = out / ARTIFACTS["manifest"]
    if not manifest_path.exists():
        raise ContractError(f"no manifest at {manifest_path}; nothing frozen to verify")
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ContractError(f"unknown manifest schema {manifest.get('schema')!r}")
    if manifest.get("status") != "frozen" and not args.allow_stub:
        raise ContractError(
            f"manifest status is {manifest.get('status')!r}; a stub build is not a "
            "frozen artifact (use --allow-stub to verify it as a smoke test)")

    checks: list[tuple[str, bool, str]] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        checks.append((name, bool(ok), detail))

    # 1. Every artifact hash and size.  Stop before parsing anything: a
    #    tampered artifact must not be read as if it were evidence.
    tampered = []
    for artifact in manifest["artifacts"]:
        path = out / artifact["path"]
        if not path.exists():
            tampered.append(f"{artifact['path']} (missing)")
            continue
        size_ok = path.stat().st_size == artifact["bytes"]
        hash_ok = sha256_file(path) == artifact["sha256"]
        check(f"artifact {artifact['path']} sha256", size_ok and hash_ok)
        if not (size_ok and hash_ok):
            tampered.append(artifact["path"])
    if tampered:
        raise ContractError(
            "artifact content does not match the frozen manifest: "
            f"{tampered}; refusing to verify a tampered build")

    # 2. Locked config re-hash + digest.
    config_path = out / ARTIFACTS["locked_config"]
    locked = json.loads(config_path.read_text())
    expected_digest = config_digest(build_locked_config())
    check("locked_config sha256", sha256_file(config_path) == manifest["config"]["sha256"])
    check("locked_config decisions unchanged", expected_digest == manifest["config"]["digest"],
          "the config is a fixed decision document, not a per-build output")
    for key in ("structure", "recipe", "metric", "loss", "optimizer",
                "prediction_noise", "n_folds", "seeds", "lrs", "dose_encoding"):
        check(f"locked_config has canonical {key}", key in locked)
    # The decision record the pipeline must not re-tune travels beside the schema.
    for key in ("folds", "consensus_control", "validation", "lr_selection",
                "test_units", "derived_defaults", "data_contract"):
        check(f"locked_config carries the decision record {key}", key in locked)

    # 3. Gene axis equality across corpus and both artifacts.
    gene_axis = load_gene_axis()
    cm = ad.read_h5ad(out / ARTIFACTS["cm"])
    axis_ok = list(cm.var.index.astype(str)) == list(gene_axis.index)
    check("cm var order == corpus/v1 genes.csv", axis_ok)
    check("cm obs columns == frozen contract", list(cm.obs.columns) == CM_OBS)
    tahoe_path = out / ARTIFACTS["tahoe"]
    header = None
    if tahoe_path.exists():
        header = read_h5ad_header(tahoe_path)
        check("tahoe var order == corpus/v1 genes.csv",
              list(header.var.index.astype(str)) == list(gene_axis.index))
        check("tahoe obs columns == frozen contract",
              list(header.obs.columns) == TAHOE_OBS)
        check("tahoe carries no inherited uns['control_mean']",
              "control_mean" not in header.uns)
        check("tahoe X is int32 CSR", header.x_dtype == "int32"
              and int(header.x_nnz) >= 0)
        check("tahoe layers == ['log1p_cp10k']", header.layer_names == ["log1p_cp10k"])
    check("cm carries no inherited uns['control_mean']", "control_mean" not in cm.uns)

    # 4. CM value contract.
    X = np.asarray(cm.X)
    layer = np.asarray(cm.layers["log1p_cp10k"])
    check("cm X int32 non-negative", X.dtype == np.int32 and bool((X >= 0).all()))
    check("cm log layer finite float32",
          layer.dtype == np.float32 and bool(np.isfinite(layer).all()))
    check("cm lib_size == X.sum(1)",
          bool(np.array_equal(cm.obs["lib_size"].to_numpy(dtype=np.int64),
                              X.sum(axis=1, dtype=np.int64))))

    # 5. Fold invariants + main-set population.
    fold_table = pd.read_parquet(out / ARTIFACTS["folds"])
    main = cm.obs.loc[cm.obs["main_set"].astype(bool)]
    main_drugs = set(main["drug_key"].astype(str))
    try:
        facts = folds_mod.assert_fold_invariants(
            fold_table, fold_table[["drug_key", "cluster_id", "test_fold",
                                    "valid_fold"]].drop_duplicates(),
            expected_drugs=main_drugs, n_folds=N_FOLDS)
        check("fold invariants (test once / valid once / no straddle)", True,
              f"{facts['n_drugs']} drugs")
    except ContractError as exc:
        check("fold invariants (test once / valid once / no straddle)", False, str(exc))
    decisions = cm.obs["eligibility_decision"].astype(str)
    excluded = set(cm.obs.loc[~cm.obs["is_control"].astype(bool)
                             & decisions.isin(["exclude", "uncertain"]),
                             "drug_key"].astype(str))
    check("exclude/uncertain drugs have no fold",
          not (excluded & set(fold_table["drug_key"].astype(str))))
    check("fold drugs == main_set drugs", set(fold_table["drug_key"].astype(str)) == main_drugs)
    counts = manifest["counts"]
    check("main_set drug count reconciles",
          counts["cm"]["main_set_drugs"] == len(main_drugs),
          f"{len(main_drugs)}")
    check("main_set well count reconciles", counts["cm"]["main_set_wells"] == len(main))

    # 6. Consensus: provenance, weights, recomputation.
    membership = pd.read_parquet(out / ARTIFACTS["consensus_membership"])
    npz = np.load(out / ARTIFACTS["consensus_controls"], allow_pickle=False)
    try:
        prov = consensus_mod.assert_consensus_provenance(
            membership, cm.obs, fold_table, n_folds=N_FOLDS)
        check("consensus provenance (control DMSO qc_pass rows on training plates, "
              "equal plate weights)", True,
              f"{prov['total_members']} member rows")
    except ContractError as exc:
        check("consensus provenance", False, str(exc))
    recomputed = consensus_mod.recompute_c_star(layer, membership, n_folds=N_FOLDS)
    check("consensus c* recomputes from membership",
          bool(np.allclose(recomputed, npz["c_star"], atol=1e-6, equal_nan=False)),
          f"max|delta|={float(np.max(np.abs(recomputed - npz['c_star']))):.2e}")
    check("consensus c* is finite float32",
          npz["c_star"].dtype == np.float32 and bool(np.isfinite(npz["c_star"]).all()))
    check("consensus carries no gene axis mismatch", npz["c_star"].shape == (N_FOLDS, N_GENES))

    # 7. Handoff reconciliation and frozen-source immutability.
    if manifest["handoffs"]:
        for record in manifest["handoffs"]:
            path = DATA_ROOT / record["path"]
            ok = path.exists() and sha256_file(path) == record["sha256"]
            check(f"handoff {record['kind']} unchanged", ok)
        eligibility_record = next((r for r in manifest["handoffs"]
                                   if r["kind"] == "eligibility"), None)
        if eligibility_record is None:
            check("eligibility handoff recorded", False, "missing from manifest")
        else:
            eligibility, _ = load_eligibility(DATA_ROOT / eligibility_record["path"])
            check("eligibility decision values valid",
                  set(eligibility["eligibility_decision"]) <= set(ELIGIBILITY_DECISIONS))
            resolvable = (set(cm.obs.loc[~cm.obs["is_control"].astype(bool),
                                          "drug_key"].astype(str)) - {""})
            uncovered = resolvable - set(eligibility["drug_key"].astype(str))
            check("eligibility covers every resolvable CM drug_key", not uncovered,
                  f"missing={sorted(uncovered)[:5]}")

        # Per-observation QC reconciliation: the CM artifact holds exactly the
        # QC-pass rows, and every Tahoe row carries the handoff's decision.
        cm_ids = pd.Index(cm.obs.index.astype(str))
        tahoe_obs = header.obs if header is not None else None
        cm_qc_table = None
        for record in manifest["handoffs"]:
            if record["kind"] not in ("cm_qc", "tahoe_qc"):
                continue
            kind = "cm" if record["kind"] == "cm_qc" else "tahoe"
            table, _ = read_qc_handoff(DATA_ROOT / record["path"], kind)
            qc_map = table.set_index("well")
            if kind == "cm":
                cm_qc_table = qc_map
                covered = cm_ids.isin(set(qc_map.index))
                check("QC handoff covers every CM artifact row", bool(covered.all()),
                      f"{int((~covered).sum())} uncovered")
                passed = qc_map["qc_pass"].reindex(cm_ids)
                check("every CM artifact row is QC-pass", bool(passed.fillna(False).all()))
                check("CM kept-row count reconciles with the handoff",
                      int(len(cm.obs)) == int(passed.fillna(False).sum()),
                      f"{len(cm.obs)} rows")
                if "matched_control_ok" in qc_map.columns:
                    expected = qc_map["matched_control_ok"].reindex(cm_ids)
                    treated = ~cm.obs["is_control"].astype(bool).to_numpy()
                    got = cm.obs["matched_control_ok"].astype(bool).to_numpy()
                    agree = np.where(expected.isna().to_numpy(), got, 
                                     expected.to_numpy(dtype=object) == got)
                    check("CM matched_control_ok matches the handoff",
                          bool(np.asarray(agree, dtype=bool).all()))
            else:
                if tahoe_obs is None:
                    continue
                obs_ids = pd.Index(tahoe_obs.index.astype(str))
                covered = obs_ids.isin(set(qc_map.index))
                check("QC handoff covers every Tahoe artifact row", bool(covered.all()),
                      f"{int((~covered).sum())} uncovered")
                expected = qc_map["qc_pass"].reindex(obs_ids)
                got = tahoe_obs["qc_pass"].astype(bool).to_numpy()
                check("Tahoe qc_pass matches the handoff row by row",
                      bool((expected.to_numpy(dtype=object) == got).all()))
                check("Tahoe qc_pass count reconciles",
                      int(got.sum()) == int(len(tahoe_obs)
                                            - counts["tahoe"]["qc_fail_rows"]),
                      f"{int(got.sum())} pass")
                if "matched_control_ok" in qc_map.columns:
                    expected = qc_map["matched_control_ok"].reindex(obs_ids)
                    got = tahoe_obs["matched_control_ok"].astype(bool).to_numpy()
                    agree = np.where(expected.isna().to_numpy(), got,
                                     expected.to_numpy(dtype=object) == got)
                    check("Tahoe matched_control_ok matches the handoff",
                          bool(np.asarray(agree, dtype=bool).all()))

        # The raw temp_env replicate id must agree between the two independent
        # routes: our re-join from the design table and #2's per-observation QC.
        if cm_qc_table is not None:
            from contract import _read_table
            source_record = next((r for r in manifest["handoffs"]
                                  if r["kind"] == "cm_qc"), None)
            raw = _read_table(DATA_ROOT / source_record["path"])
            if "raw_replicate" in raw.columns:
                if "source" in raw.columns:
                    raw = raw.loc[raw["source"].astype(str) == "temp_env"]
                elif "dataset_id" in raw.columns:
                    raw = raw.loc[raw["dataset_id"].astype(str) == "temp_env"]
                key_column = "observation_id" if "observation_id" in raw.columns else "well"
                raw = raw.set_index(raw[key_column].astype(str))
                expected = pd.to_numeric(raw["raw_replicate"].reindex(cm_ids),
                                         errors="coerce").to_numpy()
                got = cm.obs["replicate_raw"].to_numpy(dtype=np.float64)
                same = np.isclose(expected, got, equal_nan=False)
                check("cm replicate_raw agrees with #2's raw_replicate",
                      bool(same.all()),
                      f"{int((~same).sum())} mismatches")
            check("matched-control provenance recorded",
                  counts["cm"]["matched_control_source"] in ("handoff", "derived"))
    upstream = {
        "cm_view": (CM_BUNDLE_MANIFEST, "view", "temp_env_all-smiles.h5ad"),
        "tahoe_view": (TAHOE_VIEW_MANIFEST, "view", "tahoe_cm_schema"),
    }
    for record in manifest["sources"]:
        path = DATA_ROOT / record["path"]
        if not path.exists():
            check(f"source present {record['role']}", False, "missing")
            continue
        ok = path.stat().st_size == record["bytes"]
        kind = record.get("digest_kind", "file")
        if kind == "file" and record["content_sha256"]:
            ok = ok and sha256_file(path) == record["content_sha256"]
        elif kind == "content":
            # the upstream digest hashes the h5ad arrays; re-deriving it here
            # would duplicate 02/04 QA.  Verify the transcription against the
            # frozen upstream manifest *and* the recorded byte count instead.
            manifest_path, key, value = upstream[record["role"]]
            row = _frozen_row(manifest_path, key, value)
            ok = (ok and str(row["content_sha256"]) == record["content_sha256"]
                  and int(row["bytes"]) == record["bytes"])
        check(f"source unchanged {record['role']}", ok)
    for record in manifest["source_manifests"]:
        check(f"source manifest unchanged {Path(record['path']).name}",
              sha256_file(DATA_ROOT / record["path"]) == record["sha256"])

    # 8. Tahoe counts + control membership consistency.
    if header is not None:
        obs = header.obs
        usable = obs["tahoe_usable"].astype(bool)
        check("tahoe usable treated reconciles",
              counts["tahoe"]["usable_treated"] == int(usable.sum()),
              f"{int(usable.sum())}")
        check("tahoe split covers only train|valid",
              set(obs["tahoe_split"].astype(str)) == {"train", "valid"})
        check("tahoe time_h is NaN everywhere",
              bool(obs["time_h"].isna().all()))
        check("tahoe time_source is the cited literature provenance",
              set(obs["time_source"].astype(str)) == {TIME_SOURCE_TAHOE})
        check("tahoe structure-less rows retained",
              int((~obs["has_structure"].astype(bool)).sum()) ==
              counts["tahoe"]["has_structure_false"])
        groups = pd.read_parquet(out / ARTIFACTS["tahoe_control_groups"])
        pairs = np.load(out / ARTIFACTS["tahoe_pairs"], allow_pickle=False)
        plate = obs["plate"].astype(str).to_numpy()
        line = obs["cell_line_id"].astype(str).to_numpy()
        is_control = obs["is_control"].astype(bool).to_numpy()
        qc_pass = obs["qc_pass"].astype(bool).to_numpy()
        member_rows = groups["row"].to_numpy(dtype=np.int64)
        check("tahoe control membership rows are controls",
              bool(is_control[member_rows].all()))
        check("tahoe control membership rows are QC-pass",
              bool(qc_pass[member_rows].all()))
        check("tahoe control membership rows are DMSO vehicle",
              bool((obs["vehicle"].astype(str).to_numpy()[member_rows]
                    == "DMSO_TF").all()))
        expected_pass = np.flatnonzero(is_control & qc_pass)
        check("tahoe control membership == the QC-pass control rows",
              sorted(int(r) for r in member_rows) == [int(r) for r in expected_pass],
              f"{len(member_rows)} members vs {len(expected_pass)} QC-pass controls")
        group_key = (groups["plate"].astype(str) + "|"
                     + groups["cell_line_id"].astype(str))
        check("tahoe group_index maps one group_key",
              groups.groupby("group_index")["group_key"].nunique().eq(1).all())
        check("tahoe group_key matches its rows",
              bool((group_key.to_numpy() == groups["group_key"].astype(str).to_numpy()).all()))
        usable_mask = obs["tahoe_usable"].astype(bool).to_numpy()
        key_of_row = np.array([f"{p}|{l}" for p, l in zip(plate, line)], dtype=object)
        group_keys = set(groups["group_key"].astype(str))
        usable_without = [i for i in np.flatnonzero(usable_mask)
                          if key_of_row[i] not in group_keys]
        check("every usable Tahoe row has a QC-pass matched-control group",
              not usable_without, f"{len(usable_without)} row(s)")
        pair_rows = pairs["row"].astype(np.int64)
        pair_groups = pairs["group_index"].astype(np.int64)
        check("tahoe pairs reference treated rows",
              bool((~is_control[pair_rows]).all()))
        index_of = dict(zip(groups["group_index"].astype(np.int64),
                            groups["group_key"].astype(str)))
        mismatched = sum(1 for row, gid in zip(pair_rows, pair_groups)
                         if index_of.get(int(gid)) != key_of_row[int(row)])
        check("tahoe pairs agree with the control-group keys", mismatched == 0,
              f"{mismatched} mismatch(es)")
        check("tahoe pairs cover every treated row with a matched control",
              len(pair_rows) == counts["tahoe"]["pairs"],
              f"{len(pair_rows)} pairs")
        check("tahoe usable rows are all paired",
              int(np.asarray(pairs["usable"], dtype=bool).sum()) == int(usable_mask.sum()),
              f"{int(np.asarray(pairs['usable'], dtype=bool).sum())} paired usable")
        # measurement mask is membership, not counts
        ens = tahoe_gene_space_ensembl()
        measured = gene_measured_tahoe(gene_axis, ens)
        check("tahoe gene_measured == Ensembl membership",
              bool((header.var["gene_measured"].to_numpy() == measured.to_numpy()).all()))
        check("only FAM153B is unmeasured in Tahoe",
              int((~measured).sum()) == 1 and not bool(measured.loc["FAM153B"])
              if "FAM153B" in measured.index else False)

    # 8b. contract.json describes every field downstream code may read.
    contract_path = out / ARTIFACTS["contract"]
    if not contract_path.exists():
        check("contract.json present", False, "missing")
    else:
        document = json.loads(contract_path.read_text())
        described = document.get("artifacts", {})
        required_artifacts = [ARTIFACTS["cm"], ARTIFACTS["folds"],
                              ARTIFACTS["consensus_controls"],
                              ARTIFACTS["consensus_membership"],
                              ARTIFACTS["locked_config"]]
        if header is not None:
            required_artifacts += [ARTIFACTS["tahoe"],
                                   ARTIFACTS["tahoe_control_groups"],
                                   ARTIFACTS["tahoe_pairs"]]
        check("contract.json describes every frozen artifact",
              all(name in described for name in required_artifacts),
              f"missing={[n for n in required_artifacts if n not in described]}")
        cm_fields = {entry["name"] for entry in described[ARTIFACTS["cm"]]["obs_columns"]}
        for field in ("main_set", "expansion_set", "eligibility_decision", "qc_pass",
                      "observation_id", "physical_well", "replicate_raw", "solvent",
                      "chemical_index", "matched_control_ok", "drug_key", "dose_um",
                      "is_control", "vehicle", "plate", "replicate"):
            check(f"contract.json carries cm field {field}", field in cm_fields)
        if header is not None:
            ta_fields = {entry["name"]
                         for entry in described[ARTIFACTS["tahoe"]]["obs_columns"]}
            for field in ("tahoe_split", "tahoe_usable", "physical_well",
                          "cell_line_id", "observation_id", "qc_pass", "time_source",
                          "matched_control_ok"):
                check(f"contract.json carries tahoe field {field}", field in ta_fields)
            check("contract.json exposes c* as an array",
                  "c_star" in described[ARTIFACTS["consensus_controls"]]["arrays"])
            check("contract.json exposes tahoe_pairs arrays",
                  {"row", "group_index", "usable", "split"}
                  == set(described[ARTIFACTS["tahoe_pairs"]]["arrays"]))
        fold_columns = {entry["name"]
                        for entry in described[ARTIFACTS["folds"]]["columns"]}
        check("contract.json exposes folds.parquet keys",
              {"drug_key", "fold", "role", "plate"} <= fold_columns)
        unmeaning = [entry["name"]
                     for name in required_artifacts
                     for entry in described[name].get("obs_columns",
                                                      described[name].get("columns", []))
                     if not entry["meaning"]]
        check("contract.json gives every field a meaning", not unmeaning,
              f"{unmeaning[:5]}")

    # 9. counts.csv re-flatten agrees with the manifest.
    table = pd.read_csv(out / ARTIFACTS["counts"], keep_default_na=False)
    flat = dict(zip(table["key"].astype(str), table["value"].astype(str)))
    check("counts.csv rows == manifest counts keys",
          len(flat) == len(_flatten_dict(manifest["counts"])),
          f"{len(flat)} rows")
    check("manifest counts match counts.csv",
          all(flat.get(k) == v for k, v in _flatten_dict(manifest["counts"]).items()))

    failed = [c for c in checks if not c[1]]
    for name, ok, detail in checks:
        log(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail and not ok else ""))
    log(f"verify: {len(checks) - len(failed)}/{len(checks)} checks passed")
    if failed:
        raise ContractError(f"{len(failed)} verification check(s) failed")
    return {"checks": len(checks), "failed": 0}


def write_contract_json(out: Path, cm: dict, tahoe: dict | None,
                        cm_var: pd.DataFrame, tahoe_var: pd.DataFrame | None,
                        counts: dict, handoff_checks: dict) -> None:
    """The machine-readable interface ticket #4/#5/#6/#7/#8/#9/#11 read.

    It lists every frozen file with its keys, columns, shapes and dtypes and the
    meaning of every field, plus the read rules that must not be re-derived.
    """
    from contract import FIELD_MEANINGS, N_GENES, N_FOLDS

    def columns(frame: pd.DataFrame) -> list[dict]:
        return [{"name": str(column), "dtype": str(dtype),
                 "meaning": FIELD_MEANINGS.get(str(column), "")}
                for column, dtype in frame.dtypes.items()]

    membership = cm["membership"]
    summary = cm["consensus_summary"]
    fold_table = cm["folds"]
    document = {
        "schema": "06.3_input/contract/1",
        "version": TOOL_VERSION,
        "contract_id": CONTRACT_ID,
        "purpose": "the single frozen two-source 18,164-gene input formal training reads; "
                   "downstream tickets read these files and never re-derive folds, c* or config",
        "gene_axis": {
            "n_genes": N_GENES,
            "order": "stage2/02_cm_data/output/corpus/v1/genes.csv (HGNC symbol, frozen order)",
            "identity_columns": ["hgnc_id", "entrez_id", "ensembl_id", "match_type"],
            "measurement_mask": "var.gene_measured is PER SOURCE: all True in cm_transfer.h5ad; "
                                "False only for FAM153B in tahoe_transfer.h5ad.  It is gene-space "
                                "membership, never a non-zero-count test (18 all-zero Tahoe "
                                "columns are measured).",
        },
        "artifacts": {},
        "fold_semantics": {
            "n_folds": N_FOLDS,
            "seed": FOLD_SEED,
            "split_unit": "drug_key",
            "roles": ["train", "valid", "test"],
            "primary_key": ["drug_key", "fold"],
            "population": "main_set drugs only; every main-set drug is test in exactly one "
                          "fold and valid in exactly one fold",
            "excluded": "exclude/uncertain drugs have no fold row and are never test/valid",
            "plate_confounding": "a drug split is not a plate split; a fold's test drugs bring their plates",
        },
        "control_semantics": {
            "cm_train": "mean of QC-pass DMSO rows on the same plate (obs.control_group mixes "
                        "DMSO+Media and is NOT a matched baseline)",
            "cm_validation_test_inference": "one fixed per-fold c* from consensus_controls.npz; "
                                            "never a per-plate control",
            "cm_truth_delta": "treated log1p_cp10k minus the same-plate QC-pass DMSO log mean, "
                              "computed only at scoring time",
            "tahoe_train_validation": "mean of the same (plate, cell_line) group's QC-pass control "
                                      "rows on the log layer (tahoe_control_groups.parquet + "
                                      "tahoe_pairs.npz); never across cell lines",
            "consensus_provenance": "only control rows of plates carrying a TRAINING drug; "
                                    "per-plate log mean then equal plate weight",
        },
        "role_selection": {
            "main_set": "cm obs.main_set == True (include AND qc_pass) is the fold/score population",
            "expansion_set": "cm obs.expansion_set == True (exclude AND qc_pass) is ticket #11's "
                             "expansion training arm only; never scored, no fold row",
            "uncertain": "never trained, never scored, no fold row",
            "tahoe": "tahoe obs.tahoe_usable == True is trainable/evaluable; tahoe_split selects "
                     "train vs valid; time_h never filters",
        },
        "locked_config": {"file": ARTIFACTS["locked_config"],
                          "meaning": "the only recipe: model, dose encoding, loss, optimiser, "
                                     "seeds, LR rule, metric, prediction policy, 30-unit policy"},
        "renames": {},
        "artifacts_tracked_in_git": sorted(TRACKED_RECORD),
        "artifacts_local_only": [ARTIFACTS["cm"], ARTIFACTS["tahoe"]],
        "handoff_checks": handoff_checks,
        "frozen_counts": counts,
    }

    document["artifacts"][ARTIFACTS["cm"]] = {
        "format": "anndata h5ad (dense)",
        "local_only": True,
        "shape": [int(len(cm["obs"])), N_GENES],
        "X": {"encoding": "dense ndarray", "dtype": "int32", "meaning": "raw counts, non-negative"},
        "layers": {"log1p_cp10k": {"encoding": "dense ndarray", "dtype": "float32",
                                   "meaning": "row-wise log1p(counts / lib_size * 1e4)"}},
        "obs_index": {"name": "well", "meaning": FIELD_MEANINGS["well"]},
        "obs_columns": columns(cm["obs"]),
        "var_index": {"name": "symbol", "meaning": FIELD_MEANINGS["symbol"]},
        "var_columns": columns(cm_var),
        "uns": {"schema": "input contract version", "column_provenance": "per-column notes",
                "control_cache": "why no uns['control_mean'] is carried"},
    }
    if tahoe is not None:
        document["artifacts"][ARTIFACTS["tahoe"]] = {
            "format": "anndata h5ad (CSR; X and layers streamed from the 04 view)",
            "local_only": True,
            "shape": [int(len(tahoe["obs"])), N_GENES],
            "X": {"encoding": "csr_matrix", "dtype": "int32", "meaning": "raw counts, non-negative"},
            "layers": {"log1p_cp10k": {"encoding": "csr_matrix", "dtype": "float32",
                                       "meaning": "row-wise log1p(counts / lib_size * 1e4)"}},
            "obs_index": {"name": "well", "meaning": FIELD_MEANINGS["well"]},
            "obs_columns": columns(tahoe["obs"]),
            "var_index": {"name": "symbol", "meaning": FIELD_MEANINGS["symbol"]},
            "var_columns": columns(tahoe_var),
            "uns": {"schema": "input contract version", "column_provenance": "per-column notes",
                    "control_cache": "why no uns['control_mean'] is carried"},
        }
    document["artifacts"][ARTIFACTS["folds"]] = {
        "format": "parquet", "tracked": True, "rows": int(len(fold_table)),
        "primary_key": ["drug_key", "fold"],
        "meaning": "one row per (drug_key, fold): role is train|valid|test and every "
                   "main-set drug appears in all five folds exactly once as test and once "
                   "as valid. `plate` holds the drug's plate(s) ';'-joined (n_plates counts "
                   "them); the per-plate long form is not exported because a consumer keys "
                   "the split by drug_key alone",
        "columns": columns(fold_table),
    }
    document["artifacts"][ARTIFACTS["consensus_controls"]] = {
        "format": "numpy npz (uncompressed keys, zlib compressed file)", "tracked": False,
        "arrays": {
            "c_star": {"shape": [N_FOLDS, N_GENES], "dtype": "float32",
                       "meaning": "per-fold consensus control on the log1p_cp10k scale"},
            "fold_ids": {"shape": [N_FOLDS], "dtype": "int16", "meaning": "0..4"},
            "plate_names": {"shape": [int(len(cm["plate_names"]))], "dtype": "str",
                            "meaning": "union of contributing plates across folds"},
            "plate_weights": {"shape": [N_FOLDS, int(len(cm["plate_names"]))], "dtype": "float32",
                              "meaning": "1/n_plates in that fold, 0 where the plate is not used"},
            "n_plates": {"shape": [N_FOLDS], "dtype": "int16", "meaning": "contributing plates"},
            "n_wells": {"shape": [N_FOLDS], "dtype": "int32", "meaning": "contributing control rows"},
            "n_genes": {"shape": [], "dtype": "int64", "meaning": "18,164"},
        },
    }
    document["artifacts"][ARTIFACTS["consensus_membership"]] = {
        "format": "parquet", "tracked": True, "rows": int(len(membership)),
        "primary_key": ["fold", "observation_id"],
        "columns": columns(membership),
    }
    if tahoe is not None:
        document["artifacts"][ARTIFACTS["tahoe_control_groups"]] = {
            "format": "parquet (long form: one row per (group, member observation))",
            "tracked": True,
            "rows": int(len(tahoe["control_groups"])),
            "primary_key": ["group_index", "row"],
            "meaning": "matched-control membership: only QC-pass DMSO_TF control rows; "
                       "a group's log-layer mean is the matched control",
            "columns": columns(tahoe["control_groups"]),
        }
        document["artifacts"][ARTIFACTS["tahoe_pairs"]] = {
            "format": "numpy npz", "tracked": True,
            "meaning": "treated row -> control group; rows without a QC-pass matched "
                       "control are absent (never re-paired)",
            "arrays": {
                "row": {"shape": [int(len(tahoe["pairs"]["row"]))],
                        "dtype": "int32", "meaning": FIELD_MEANINGS["row"]},
                "group_index": {"shape": [int(len(tahoe["pairs"]["group_index"]))],
                                "dtype": "int32",
                                "meaning": FIELD_MEANINGS["group_index"]},
                "usable": {"shape": [int(len(tahoe["pairs"]["usable"]))], "dtype": "bool",
                           "meaning": FIELD_MEANINGS["usable"]},
                "split": {"shape": [int(len(tahoe["pairs"]["split"]))], "dtype": "str",
                          "meaning": FIELD_MEANINGS["split"]},
            },
        }
    for name, source in (("counts", None), ("fold_sizes", None),
                         ("consensus_summary", None), ("var_checks", None)):
        frame = {
            "counts": None, "fold_sizes": cm["sizes"],
            "consensus_summary": summary, "var_checks": None,
        }[name]
        entry = {"format": "csv", "tracked": True}
        if frame is not None:
            entry["rows"] = int(len(frame))
            entry["columns"] = columns(frame)
        document["artifacts"][ARTIFACTS[name]] = entry
    document["artifacts"][ARTIFACTS["locked_config"]] = {
        "format": "json", "tracked": True,
        "meaning": "locked experiment config; hashed in manifest.json",
    }
    document["artifacts"][ARTIFACTS["manifest"]] = {
        "format": "json", "tracked": True,
        "meaning": "sha256/bytes of every artifact, source and handoff hashes, counts, git SHA",
    }
    document["artifacts"][ARTIFACTS["contract"]] = {
        "format": "json", "tracked": True, "meaning": "this document",
    }
    (out / ARTIFACTS["contract"]).write_text(canonical_json(document), encoding="utf-8")
    log(f"contract.json written ({len(document['artifacts'])} artifacts described)")


def _flatten_dict(counts: dict) -> dict[str, str]:
    flat: dict[str, str] = {}
    _flatten_counts("", counts, flat)
    return flat


# --------------------------------------------------------------------------- #
# --check
# --------------------------------------------------------------------------- #

def _discover_qc(args, kind: str):
    """Resolve one per-source QC table, falling back to #2's combined table."""
    override = {
        "cm_qc": args.qc or args.cm_qc,
        "tahoe_qc": args.qc or args.tahoe_qc,
    }[kind]
    if override:
        return Path(override).expanduser().resolve(), "override"
    try:
        return resolve_handoff(kind), "per-source"
    except ContractError:
        pass
    return resolve_handoff("qc_combined"), "combined"


def check(args) -> dict:
    """Validate the frozen sources and the resolved handoffs.  Writes nothing.

    This is the pre-freeze gate: it reads the 289 MB CM view's metadata and both
    handoffs (never the 16 GB Tahoe matrix), validates their schema and reports
    what the build would do.  A handoff that is absent is reported as blocking
    ticket #1/#2; a handoff that is present but malformed raises.
    """
    gene_axis = load_gene_axis()
    records = source_records()
    log(f"gene axis: {N_GENES} symbols, {gene_axis['ensembl_id'].nunique()} Ensembl ids")
    for record in records:
        log(f"source {record.role}: {relpath(record.path)} ({record.bytes:,} B)")

    cm_preview = load_cm_view(need_matrix=False)
    cm_obs = cm_preview["obs"]
    log(f"CM view: {len(cm_obs)} rows x {N_GENES:,} genes "
        f"({int((~cm_obs['is_control'].astype(bool)).sum())} treated)")
    report: dict = {"handoffs": {}, "blocked": []}

    # --- #1 eligibility -----------------------------------------------------
    try:
        path = resolve_handoff("eligibility", args.eligibility)
        table, record = load_eligibility(path)
        counts = table["eligibility_decision"].value_counts().to_dict()
        log(f"handoff eligibility: {relpath(record.path)} "
            f"({len(table)} keys; " + ", ".join(
                f"{k}={counts.get(k, 0)}" for k in ELIGIBILITY_DECISIONS) + ")")
        report["handoffs"]["eligibility"] = record.as_dict()
    except ContractError as exc:
        log(f"handoff eligibility: BLOCKED — {exc}")
        report["blocked"].append("eligibility")

    # --- #2 per-observation QC ---------------------------------------------
    for kind, label in (("cm_qc", "CM"), ("tahoe_qc", "Tahoe")):
        try:
            path, how = _discover_qc(args, kind)
            table, record = read_qc_handoff(path, "cm" if kind == "cm_qc" else "tahoe")
            passed = int(table["qc_pass"].astype(bool).sum())
            log(f"handoff {kind}: {relpath(record.path)} [{how}] "
                f"({len(table)} rows, {passed} QC-pass)")
            report["handoffs"][kind] = {**record.as_dict(), "resolved_as": how}
        except ContractError as exc:
            log(f"handoff {kind}: BLOCKED — {exc}")
            report["blocked"].append(kind)

    log("tahoe valid plates: plate13, plate14 (existing technical default)")
    log("folds: 5 by drug_key, seed 20261102, capacity ceil(n_main_drugs/5)")
    if report["blocked"]:
        log("check: PASS for the frozen sources — still blocked on "
            f"{', '.join(report['blocked'])}")
    else:
        log("check: PASS — frozen sources readable and both handoffs valid; "
            "ready to freeze")
    return report


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Build/verify the frozen two-source 18,164-gene 06.3 input.")
    parser.add_argument("--out-root", default=None,
                        help="directory whose 06.3_input/output/ receives the artifacts "
                             f"(default: {DATA_ROOT / 'stage2/06_prnet_transfer'})")
    parser.add_argument("--check", action="store_true",
                        help="validate frozen sources and handoffs, write nothing")
    parser.add_argument("--verify", action="store_true",
                        help="offline verification of an existing frozen build")
    parser.add_argument("--eligibility", default=None, help="#1 eligibility table")
    parser.add_argument("--cm-qc", default=None, help="#2 CM per-observation QC table")
    parser.add_argument("--tahoe-qc", default=None, help="#2 Tahoe per-observation QC table")
    parser.add_argument("--qc", default=None,
                        help="#2 combined per-observation QC table with a `source` column")
    parser.add_argument("--stub-handoffs", action="store_true",
                        help="synthetic include-all/pass-all handoffs (tests/smoke only)")
    parser.add_argument("--allow-stub", action="store_true",
                        help="let --verify accept a manifest whose status is `stub`")
    parser.add_argument("--tahoe-rows", type=int, default=None,
                        help="smoke only: write the first N Tahoe rows")
    parser.add_argument("--skip-tahoe", action="store_true",
                        help="build the CM side only (development)")
    parser.add_argument("--force", action="store_true",
                        help="overwrite an existing output directory")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        if args.check:
            check(args)
        elif args.verify:
            verify(args)
        else:
            build(args)
    except ContractError as exc:
        print(f"[06.3_input] FAIL: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
