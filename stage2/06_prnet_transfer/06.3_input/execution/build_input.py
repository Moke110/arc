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
    _frozen_manifest_row as _frozen_row,
    ARTIFACTS, CM_BUNDLE_MANIFEST, CM_OBS, CM_VIEW, CONFIG_VERSION, CONTRACT_ID,
    ContractError, DATA_ROOT, ELIGIBILITY_DECISIONS, FOLD_SEED, MANIFEST_SCHEMA,
    N_FOLDS, N_GENES, STEP_DATA, TAHOE_OBS, TAHOE_VIEW, TAHOE_VIEW_MANIFEST,
    TIME_SOURCE_NOTE, TIME_SOURCE_TAHOE, TOOL_VERSION, TRACKED_RECORD,
    build_locked_config, canonical_json, config_digest, gene_measured_cm,
    gene_measured_tahoe, git_state, load_eligibility, load_gene_axis, load_qc,
    output_dir, relpath, runtime_env, sha256_file, source_records,
    split_combined_qc, tahoe_gene_space_ensembl,
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
        combined = pd.read_csv(args.qc, keep_default_na=False, dtype=str)
        parts = split_combined_qc(combined, Path(args.qc))
        cm_part = parts.get("cm")
        if cm_part is None:
            raise ContractError(f"{args.qc}: combined QC table has no CM rows")
        scratch = Path(args.qc).with_suffix(".cm_scratch.csv")
        cm_part.to_csv(scratch, index=False)
        try:
            cm_qc, record = load_qc(scratch, "cm")
        finally:
            scratch.unlink()
        records.append(record)
    elif args.cm_qc:
        cm_qc, record = load_qc(args.cm_qc, "cm")
        records.append(record)
    else:
        from contract import resolve_handoff
        cm_qc, record = load_qc(resolve_handoff("cm_qc"), "cm")
        records.append(record)
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
        combined = pd.read_csv(args.qc, keep_default_na=False, dtype=str)
        parts = split_combined_qc(combined, Path(args.qc))
        part = parts.get("tahoe")
        if part is None:
            raise ContractError(f"{args.qc}: combined QC table has no Tahoe rows")
        scratch = Path(args.qc).with_suffix(".tahoe_scratch.csv")
        part.to_csv(scratch, index=False)
        try:
            return load_qc(scratch, "tahoe")
        finally:
            scratch.unlink()
    if args.tahoe_qc:
        return load_qc(args.tahoe_qc, "tahoe")
    from contract import resolve_handoff
    return load_qc(resolve_handoff("tahoe_qc"), "tahoe")


# --------------------------------------------------------------------------- #
# CM side
# --------------------------------------------------------------------------- #

def load_cm_view(need_matrix: bool = True) -> dict:
    import anndata as ad
    if not CM_VIEW.exists():
        raise ContractError(f"missing frozen CM view: {CM_VIEW}")
    adata = ad.read_h5ad(CM_VIEW)
    var = adata.var
    if list(var.index.astype(str)) != list(load_gene_axis().index):
        raise ContractError(f"{CM_VIEW}: var is not the frozen gene axis")
    out = {"obs": adata.obs.copy(), "var": var.copy(), "uns": dict(adata.uns)}
    if need_matrix:
        X = np.asarray(adata.X)
        layer = np.asarray(adata.layers["log1p_cp10k"])
        if X.dtype != np.int32 or layer.dtype != np.float32:
            raise ContractError(
                f"{CM_VIEW}: X {X.dtype} layer {layer.dtype}, expected int32/float32")
        out["X"] = X
        out["log1p"] = layer
    adata.file.close() if adata.file is not None else None
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
                                         "drug_key"].astype(str).nunique() - 1),
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
    counts["eligibility_include_without_usable_well"] = int(
        len(include_keys - set(main_drugs["drug_key"])))
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
    fold_table = folds_mod.expand_fold_table(assign, main_drugs, n_folds=N_FOLDS)
    fold_facts = folds_mod.assert_fold_invariants(
        fold_table, assign, expected_drugs=set(main_drugs["drug_key"].astype(str)),
        n_folds=N_FOLDS)
    sizes = folds_mod.fold_sizes(fold_table)
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
        consensus_mod.build_consensus(log1p, kept, fold_table,
                                      n_folds=N_FOLDS, n_genes=N_GENES))
    provenance = consensus_mod.assert_consensus_provenance(
        membership, kept, fold_table, n_folds=N_FOLDS)
    counts["consensus"] = {
        "fold_n_plates": [int(v) for v in summary["n_plates"]],
        "fold_n_wells": [int(v) for v in summary["n_wells"]],
        "members_total": int(len(membership)),
        "per_fold_plates": provenance["per_fold"],
    }

    result = {
        "obs": kept, "X": X, "log1p": log1p, "folds": fold_table,
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
    import anndata as ad
    view = ad.read_h5ad(TAHOE_VIEW, backed="r")
    try:
        obs_all = view.obs.copy()
        var = view.var.copy()
    finally:
        view.file.close()
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
    pairs = writers.tahoe_pairs(obs)

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
        "control_groups": int(len(control_groups)),
        "groups_with_zero_qc_pass": int((control_groups["n_qc_pass"] == 0).sum()),
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
            "tahoe_usable": "~is_control AND has_structure AND qc_pass AND the control group "
                            "keeps at least one QC-pass control",
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
    out.mkdir(parents=True, exist_ok=True)
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

    cm = build_cm({"eligibility": eligibility, "cm_qc": cm_qc},
                  gene_axis, write=True, out=out)
    log(f"CM: kept {cm['counts']['kept_rows']} rows "
        f"({cm['counts']['kept_treated']} treated / {cm['counts']['kept_controls']} control), "
        f"main_set {cm['counts']['main_set_drugs']} drugs / "
        f"{cm['counts']['main_set_wells']} wells")

    cm_measured = gene_measured_cm(gene_axis)
    cm_nnz = (cm["X"] != 0).sum(axis=0)
    cm_var = writers.var_frame(gene_axis, cm_measured, cm_nnz)

    # ---- Tahoe side --------------------------------------------------------
    tahoe = None
    tahoe_measured = None
    if args.skip_tahoe:
        log("Tahoe side skipped (--skip-tahoe)")
    else:
        log("reading the frozen Tahoe view (backed, obs only)…")
        import anndata as ad
        preview = ad.read_h5ad(TAHOE_VIEW, backed="r")
        try:
            tahoe_index = preview.obs.index.astype(str).to_numpy()
        finally:
            preview.file.close()
        if args.tahoe_rows is not None:
            tahoe_index = tahoe_index[: int(args.tahoe_rows)]
        tahoe_qc, tahoe_qc_record = load_tahoe_qc(args, tahoe_index, stub=stub)
        if tahoe_qc_record is not None:
            handoff_records = handoff_records + [tahoe_qc_record]
        ens_space = tahoe_gene_space_ensembl()
        tahoe_measured = gene_measured_tahoe(gene_axis, ens_space)
        tahoe = build_tahoe(tahoe_qc, gene_axis, tahoe_rows=args.tahoe_rows,
                            tahoe_measured=tahoe_measured)
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
            treated_row=tahoe["pairs"]["treated_row"],
            control_group=tahoe["pairs"]["control_group"].astype(str),
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

    config_path = out / ARTIFACTS["locked_config"]
    artifacts = []
    for name in [ARTIFACTS["cm"], ARTIFACTS["tahoe"], ARTIFACTS["folds"],
                 ARTIFACTS["consensus_controls"], ARTIFACTS["consensus_membership"],
                 ARTIFACTS["tahoe_control_groups"], ARTIFACTS["tahoe_pairs"],
                 ARTIFACTS["locked_config"], ARTIFACTS["counts"],
                 ARTIFACTS["fold_sizes"], ARTIFACTS["consensus_summary"],
                 ARTIFACTS["var_checks"]]:
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
    for key in ("model", "dose_encoding", "loss", "optimisation", "metric",
                "test_units", "prediction"):
        check(f"locked_config has {key}", key in locked)

    # 3. Gene axis equality across corpus and both artifacts.
    gene_axis = load_gene_axis()
    cm = ad.read_h5ad(out / ARTIFACTS["cm"])
    axis_ok = list(cm.var.index.astype(str)) == list(gene_axis.index)
    check("cm var order == corpus/v1 genes.csv", axis_ok)
    check("cm obs columns == frozen contract", list(cm.obs.columns) == CM_OBS)
    tahoe_path = out / ARTIFACTS["tahoe"]
    tahoe = None
    if tahoe_path.exists():
        tahoe = ad.read_h5ad(tahoe_path, backed="r")
        check("tahoe var order == corpus/v1 genes.csv",
              list(tahoe.var.index.astype(str)) == list(gene_axis.index))
        check("tahoe obs columns == frozen contract",
              list(tahoe.obs.columns) == TAHOE_OBS)
        check("tahoe carries no inherited uns['control_mean']",
              "control_mean" not in tahoe.uns)
        X = tahoe.X
        check("tahoe X is int32 CSR", getattr(X, "dtype", None) == np.int32)
        check("tahoe layers == ['log1p_cp10k']",
              sorted(map(str, tahoe.layers.keys())) == ["log1p_cp10k"])
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
        eligibility = pd.read_csv(DATA_ROOT / manifest["handoffs"][0]["path"],
                                  keep_default_na=False, dtype=str)
        decisions_csv = eligibility.iloc[:, 1].astype(str).str.strip().str.lower()
        check("eligibility decision values valid",
              set(decisions_csv) <= set(ELIGIBILITY_DECISIONS))
        resolvable = (set(cm.obs.loc[~cm.obs["is_control"].astype(bool),
                                      "drug_key"].astype(str)) - {""})
        uncovered = resolvable - set(eligibility.iloc[:, 0].astype(str))
        check("eligibility covers every resolvable CM drug_key", not uncovered,
              f"missing={sorted(uncovered)[:5]}")
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
    if tahoe is not None:
        obs = tahoe.obs
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
        check("tahoe control groups have >= 1 QC-pass control",
              bool((groups["n_qc_pass"] >= 1).all()))
        expected_rows = sorted(int(r) for rows in groups["control_rows"] for r in rows)
        actual_controls = np.flatnonzero(obs["is_control"].astype(bool).to_numpy())
        check("tahoe control membership == the control rows",
              expected_rows == [int(r) for r in actual_controls])
        check("tahoe pairs cover every treated row",
              len(pairs["treated_row"]) == int((~obs["is_control"].astype(bool)).sum()))
        check("tahoe pairs groups exist",
              set(pairs["control_group"].astype(str)) <=
              set(groups["control_group"].astype(str)))
        # measurement mask is membership, not counts
        ens = tahoe_gene_space_ensembl()
        measured = gene_measured_tahoe(gene_axis, ens)
        check("tahoe gene_measured == Ensembl membership",
              bool((tahoe.var["gene_measured"].to_numpy() == measured.to_numpy()).all()))
        check("only FAM153B is unmeasured in Tahoe",
              int((~measured).sum()) == 1 and not bool(measured.loc["FAM153B"])
              if "FAM153B" in measured.index else False)
        tahoe.file.close()

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


def _flatten_dict(counts: dict) -> dict[str, str]:
    flat: dict[str, str] = {}
    _flatten_counts("", counts, flat)
    return flat


# --------------------------------------------------------------------------- #
# --check
# --------------------------------------------------------------------------- #

def check(args) -> dict:
    """Validate the frozen sources and report the build plan.  Writes nothing."""
    gene_axis = load_gene_axis()
    records = source_records()
    log(f"gene axis: {N_GENES} symbols, {gene_axis['ensembl_id'].nunique()} Ensembl ids")
    for record in records:
        log(f"source {record.role}: {relpath(record.path)} ({record.bytes:,} B)")

    cm_preview = load_cm_view(need_matrix=False)
    log(f"CM view: {len(cm_preview['obs'])} rows x 18,164 genes")
    qc_available = {}
    from contract import HANDOFF_CANDIDATES
    for kind in ("eligibility", "cm_qc", "tahoe_qc"):
        found = None
        for rel in HANDOFF_CANDIDATES[kind]:
            if (STEP_DATA / rel).exists():
                found = STEP_DATA / rel
                break
        qc_available[kind] = relpath(found) if found else None
        log(f"handoff {kind}: {qc_available[kind] or 'NOT PRESENT (blocked on #1/#2)'}")
    log("tahoe valid plates: plate13, plate14 (existing technical default)")
    log("folds: 5 by drug_key, seed 20261102, capacity ceil(n_main_drugs/5)")
    log("check: PASS — frozen sources readable and on the corpus gene axis")
    return {"handoffs": qc_available}


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
