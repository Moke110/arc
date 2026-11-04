#!/usr/bin/env python
"""Build the Tahoe-100M pseudobulk corpus (ADR 0029).

Unit = (plate, well, cell_line); well is upstream `sample`.
Values = summed raw UMI counts, plus a derived log1p_cp10k layer.
Gene space = Tahoe native 62,710 (token_id 3..62712 -> column 0..62709).

Stages (each resumable; `execution/_parts/` and its `_state/` are the resumable
part-cache and bookkeeping state, not results; QA records land in
`execution/build_record/`):

    python stage2/02_dataset_corpus/execution/tahoe/build_pseudobulk.py scan
    python stage2/02_dataset_corpus/execution/tahoe/build_pseudobulk.py aggregate
    python stage2/02_dataset_corpus/execution/tahoe/build_pseudobulk.py assemble
    python stage2/02_dataset_corpus/execution/tahoe/build_pseudobulk.py qa

`--check` prints the repository root, the Tahoe source it reads and the corpus
output directory, then exits. `--plates plate1,plate2` limits any stage;
`--shards N` limits the scan/aggregate to the first N shards (smoke tests on a
partial atlas). No network: the atlas is already downloaded under
`pert_data/Tahoe-100M/source/`; the importer lives in
`pert_data/Tahoe-100M/script/`.
"""
from __future__ import annotations

import argparse
import ast
import json
import sys
import time
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq


def find_repo_root(start: Path) -> Path:
    """Walk up from *start* to the repository root (the directory holding `.git`).

    Never counts directory levels: this step's layout has already changed once
    (ADR 0030) and a hard-coded depth fails silently when it changes again.
    """
    for candidate in (start, *start.parents):
        if (candidate / ".git").exists():
            return candidate
    raise SystemExit(f"cannot find the repository root: no .git at or above {start}")


ROOT = find_repo_root(Path(__file__).resolve().parent)
assert (ROOT / "stage2").is_dir(), f"resolved repo root has no stage2/: {ROOT}"
DATA_HOME = ROOT / "pert_data" / "Tahoe-100M" / "source"
# Ticket #14 moved this builder out of the retired step `04` into the
# `02_dataset_corpus` step: the corpus is `output/corpus/v2/tahoe.h5ad`, the
# resumable per-plate cache is `execution/_parts/`, and the tracked QA/manifests
# are `execution/build_record/tahoe_pseudobulk_*`.
STEP = ROOT / "stage2" / "02_dataset_corpus"
OUT = STEP / "output" / "corpus" / "v2"
RECORD = STEP / "execution" / "build_record"
PARTS = STEP / "execution" / "_parts"
STATE = PARTS / "_state"

N_GENES = 62_710
TOKEN_OFFSET = 3          # token_id 3 is the first gene
NG = N_GENES              # column index = token_id - TOKEN_OFFSET, max 62_709
N_SHARDS = 3_388
N_PLATES = 14
WELLS_PER_PLATE = 96
CONTROL_DRUG = "DMSO_TF"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def shard_path(idx: int) -> Path:
    return DATA_HOME / "data" / f"train-{idx:05d}-of-{N_SHARDS:05d}.parquet"


def all_shards() -> list[int]:
    return [i for i in range(N_SHARDS) if shard_path(i).exists()]


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def load_sample_meta() -> dict:
    """well (sample id) -> plate, drug, dose, label."""
    t = pq.read_table(DATA_HOME / "metadata" / "sample_metadata.parquet")
    out = {}
    for r in t.to_pylist():
        pert = ast.literal_eval(r["drugname_drugconc"])
        if len(pert) != 1:
            raise SystemExit(f"{r['sample']}: {len(pert)} perturbagens, expected 1")
        name, conc, unit = pert[0]
        out[r["sample"]] = {
            "plate": r["plate"],
            "drug": name,
            "dose": float(conc),
            "unit": unit,
            "label": r["drugname_drugconc"],
            "is_control": name == CONTROL_DRUG,
        }
    return out


def load_gene_meta() -> list[dict]:
    t = pq.read_table(DATA_HOME / "metadata" / "gene_metadata.parquet")
    rows = t.to_pylist()
    assert len(rows) == N_GENES, len(rows)
    assert rows[0]["token_id"] == TOKEN_OFFSET and rows[-1]["token_id"] == N_GENES + TOKEN_OFFSET - 1
    return rows


def load_cell_lines() -> dict:
    """CVCL id -> cell_name (Vevo)."""
    t = pq.read_table(DATA_HOME / "metadata" / "cell_line_metadata.parquet",
                      columns=["Cell_ID_Cellosaur", "cell_name"])
    out = {}
    for r in t.to_pylist():
        out.setdefault(r["Cell_ID_Cellosaur"], r["cell_name"])
    return out


# ---------------------------------------------------------------------------
# scan: which shard belongs to which plate
# ---------------------------------------------------------------------------
def stage_scan(args) -> int:
    """Map each expression shard to the plate(s) its cells belong to.

    A shard is *mostly* one plate, but not always: cell-level `plate` must be
    honoured (observed: shard 236 mixes plate4 and plate5). Aggregate therefore
    reads a shard once per plate it touches and filters cells by `plate`.
    """
    STATE.mkdir(parents=True, exist_ok=True)
    cache = STATE / "shard_plate.json"
    raw = json.loads(cache.read_text()) if cache.exists() else {}
    mapping = {k: (v if isinstance(v, list) else [v]) for k, v in raw.items()}
    todo = [i for i in all_shards() if str(i) not in mapping]
    if args.shards:
        todo = todo[: args.shards]
    t0 = time.time()
    for k, idx in enumerate(todo):
        col = pq.read_table(shard_path(idx), columns=["plate"])["plate"]
        mapping[str(idx)] = sorted(set(col.unique().to_pylist()))
        if (k + 1) % 200 == 0:
            print(f"  scanned {k+1}/{len(todo)}  {time.time()-t0:.0f}s", flush=True)
            cache.write_text(json.dumps(mapping, sort_keys=True))
    cache.write_text(json.dumps(mapping, sort_keys=True))
    per_plate: dict[str, int] = {}
    mixed = {i: v for i, v in mapping.items() if len(v) > 1}
    for idx, plates in mapping.items():
        for p in plates:
            per_plate[p] = per_plate.get(p, 0) + 1
    print(f"scan: {len(mapping)} shards on disk -> {len(per_plate)} plates")
    for p in sorted(per_plate, key=lambda s: int(s.replace("plate", ""))):
        print(f"  {p:8s} {per_plate[p]:5d} shard-touches")
    print(f"mixed-plate shards: {len(mixed)}"
          + (f"  e.g. {dict(list(mixed.items())[:5])}" if mixed else ""))
    return 0


# ---------------------------------------------------------------------------
# aggregate: one dense int32 accumulator per plate, streamed to a CSR part
# ---------------------------------------------------------------------------
def _flatten(table):
    """Return (gene_cols, counts, lens) for the non-marker entries of every cell."""
    genes = table["genes"].combine_chunks()
    expr = table["expressions"].combine_chunks()
    off = genes.offsets.to_numpy().astype(np.int64)
    lens = np.diff(off) - 1                      # marker token dropped
    keep = np.ones(len(genes.values), dtype=bool)
    keep[off[:-1]] = False
    gcols = genes.values.to_numpy()[keep].astype(np.int64) - TOKEN_OFFSET
    if gcols.min() < 0 or gcols.max() >= N_GENES:
        raise SystemExit("gene token out of range")
    counts = expr.values.to_numpy()[keep].astype(np.int32)
    return gcols, counts, lens, off


def stage_aggregate(args) -> int:
    PARTS.mkdir(parents=True, exist_ok=True)
    STATE.mkdir(parents=True, exist_ok=True)
    raw = json.loads((STATE / "shard_plate.json").read_text())
    mapping = {int(k): (v if isinstance(v, list) else [v]) for k, v in raw.items()}
    sample_meta = load_sample_meta()
    by_plate: dict[str, list[int]] = {}
    for idx in sorted(mapping):
        for plate in mapping[idx]:
            by_plate.setdefault(plate, []).append(idx)
    if args.plates:
        keep = {p.strip() for p in args.plates.split(",")}
        by_plate = {p: v for p, v in by_plate.items() if p in keep}
    for plate, shards in sorted(by_plate.items(), key=lambda kv: int(kv[0][5:])):
        part = PARTS / f"{plate}.npz"
        if part.exists() and not args.force:
            print(f"{plate}: part exists, skip")
            continue
        wells = sorted(s for s, m in sample_meta.items() if m["plate"] == plate)
        if len(wells) != WELLS_PER_PLATE:
            raise SystemExit(f"{plate}: {len(wells)} wells != {WELLS_PER_PLATE}")
        well_ix = {w: i for i, w in enumerate(wells)}
        n_wells = len(wells)
        acc = np.zeros((n_wells * 50, N_GENES), dtype=np.int32)   # row = well*50 + line
        ncell = np.zeros(n_wells * 50, dtype=np.int32)
        line_ix: dict[str, int] = {}

        t0 = time.time()
        for k, idx in enumerate(shards):
            mixed = len(mapping[idx]) > 1
            cols = ["genes", "expressions", "sample", "cell_line_id"]
            if mixed:
                cols.append("plate")
            tb = pq.read_table(shard_path(idx), columns=cols)
            gcols, counts, lens, off = _flatten(tb)
            samples = tb["sample"].to_pylist()
            lines = tb["cell_line_id"].to_pylist()
            if mixed:
                here = tb["plate"].to_pylist()
                cell_ix = [i for i, p in enumerate(here) if p == plate]
            else:
                cell_ix = range(len(samples))
            kept_lines = [lines[i] for i in cell_ix]
            for ln in set(kept_lines):
                if ln not in line_ix:
                    if len(line_ix) >= 50:
                        raise SystemExit(f"{plate}: more than 50 cell lines")
                    line_ix[ln] = len(line_ix)
            wrow = np.fromiter((well_ix[samples[i]] * 50 for i in cell_ix),
                               dtype=np.int64, count=len(kept_lines))
            lrow = np.fromiter((line_ix[ln] for ln in kept_lines),
                               dtype=np.int64, count=len(kept_lines))
            row = wrow + lrow
            np.add.at(ncell, row, 1)
            src = off[:-1] + 1
            for j, i in enumerate(cell_ix):
                a = src[i]
                n = lens[i]
                if n:
                    np.add.at(acc[row[j]], gcols[a:a + n], counts[a:a + n])
            if (k + 1) % 40 == 0:
                el = time.time() - t0
                rate = (k + 1) / el
                print(f"{plate}: {k+1}/{len(shards)} shards  {el:5.0f}s  "
                      f"eta {((len(shards)-k-1)/rate):4.0f}s", flush=True)
        import scipy.sparse as sp
        rows = np.nonzero(ncell)[0]
        M = sp.csr_matrix(acc[rows])
        line_ids = [k for k, _ in sorted(line_ix.items(), key=lambda kv: kv[1])]
        np.savez_compressed(
            part,
            data=M.data.astype(np.int32),
            indices=M.indices.astype(np.int32),
            indptr=M.indptr.astype(np.int64),
            row_well_ix=(rows // 50).astype(np.int16),
            row_line_ix=(rows % 50).astype(np.int16),
            row_ncell=ncell[rows].astype(np.int32),
            wells=np.array(wells, dtype="U16"),
            line_ids=np.array(line_ids, dtype="U16"),
            plate=np.array(plate),
        )
        print(f"{plate}: {len(shards)} shards -> {len(rows)} profiles, {M.nnz:,} nonzeros "
              f"({time.time()-t0:.0f}s) -> {part.name}", flush=True)
    return 0


# ---------------------------------------------------------------------------
# assemble: parts -> sparse h5ad (streamed, never holds the corpus in RAM)
# ---------------------------------------------------------------------------
def _csr_group(f, name, n_rows, nnz, dtype, n_genes_out=N_GENES):
    """Create an anndata-compatible CSR group and return its datasets."""
    import h5py
    g = f.create_group(name)
    g.attrs["encoding-type"] = "csr_matrix"
    g.attrs["encoding-version"] = "0.1.0"
    g.attrs["shape"] = np.array([n_rows, n_genes_out], dtype=np.int64)
    data = g.create_dataset("data", shape=(nnz,), dtype=dtype, chunks=(2_000_000,))
    indices = g.create_dataset("indices", shape=(nnz,), dtype=np.int32, chunks=(2_000_000,))
    indptr = g.create_dataset("indptr", shape=(n_rows + 1,), dtype=np.int64)
    return g, data, indices, indptr


def stage_assemble(args) -> int:
    import h5py
    import pandas as pd
    import scipy.sparse as sp
    import anndata as ad

    parts = sorted(PARTS.glob("plate*.npz"), key=lambda p: int(p.stem[5:]))
    if args.plates:
        keep = {p.strip() for p in args.plates.split(",")}
        parts = [p for p in parts if p.stem in keep]
    if not parts:
        raise SystemExit("no parts to assemble")
    sample_meta = load_sample_meta()
    gene_meta = load_gene_meta()
    cell_names = load_cell_lines()
    t0 = time.time()

    obs_rows = []
    n_rows = nnz = 0
    for p in parts:
        d = np.load(p, allow_pickle=False)
        plate = str(d["plate"])
        wells = d["wells"]
        line_ids = d["line_ids"]
        obs_rows.append(pd.DataFrame({
            "plate": plate,
            "well": wells[d["row_well_ix"]],
            "cell_line_id": line_ids[d["row_line_ix"]],
            "n_cells": d["row_ncell"],
        }))
        n_rows += len(d["row_ncell"])
        nnz += len(d["data"])
    obs = pd.concat(obs_rows, ignore_index=True)
    obs["cell_name"] = obs["cell_line_id"].map(cell_names).fillna("")
    obs["drug"] = obs["well"].map(lambda w: sample_meta[w]["drug"])
    obs["dose_um"] = obs["well"].map(lambda w: sample_meta[w]["dose"])
    obs["is_control"] = obs["drug"] == CONTROL_DRUG
    obs["control_group"] = "tahoe::" + obs["plate"]
    obs["profile_id"] = (obs["plate"] + ":" + obs["well"] + ":" + obs["cell_line_id"])
    obs = obs.set_index("profile_id")
    obs["n_cells"] = obs["n_cells"].astype(np.int32)

    var = pd.DataFrame({
        "ensembl_id": [g["ensembl_id"] for g in gene_meta],
        "token_id": [g["token_id"] for g in gene_meta],
    }, index=[g["gene_symbol"] for g in gene_meta])

    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "tahoe.h5ad"
    skeleton = ad.AnnData(X=sp.csr_matrix((n_rows, N_GENES), dtype=np.int32),
                          obs=obs, var=var)
    skeleton.write_h5ad(dest)
    del skeleton
    print(f"assemble: {len(parts)} parts, {n_rows:,} profiles, {nnz:,} nonzeros; "
          f"skeleton written ({time.time()-t0:.0f}s)", flush=True)

    with h5py.File(dest, "a") as f:
        for stale in ("X", "layers"):
            if stale in f:
                del f[stale]
        f.attrs["pseudobulk_unit"] = "(plate, well, cell_line)"
        f.attrs["pseudobulk_counts"] = "summed raw UMI counts"
        _, xd, xi, xp = _csr_group(f, "X", n_rows, nnz, np.int32)
        ly = f.create_group("layers")
        _, ld, li, lp = _csr_group(ly, "log1p_cp10k", n_rows, nnz, np.float32)
        ip = np.zeros(n_rows + 1, dtype=np.int64)
        row0 = 0
        nz0 = 0
        for k, p in enumerate(parts):
            d = np.load(p, allow_pickle=False)
            indptr = d["indptr"]
            n_part = len(indptr) - 1
            M = sp.csr_matrix((d["data"].astype(np.float64), d["indices"], indptr),
                              shape=(n_part, N_GENES))
            rowsum = np.asarray(M.sum(axis=1)).ravel()
            if (rowsum <= 0).any():
                raise SystemExit(f"{p.stem}: {int((rowsum<=0).sum())} profiles with zero counts")
            scale = np.repeat(1e4 / rowsum, np.diff(M.indptr))
            xd[nz0:nz0 + len(d["data"])] = d["data"]
            xi[nz0:nz0 + len(d["indices"])] = d["indices"]
            ip[row0:row0 + n_part] = indptr[:-1] + nz0
            ld[nz0:nz0 + len(d["data"])] = np.log1p(M.data * scale).astype(np.float32)
            li[nz0:nz0 + len(d["indices"])] = d["indices"]
            lp[row0:row0 + n_part] = indptr[:-1] + nz0
            row0 += n_part
            nz0 += len(d["data"])
            print(f"  {p.stem}: {n_part:,} profiles ({time.time()-t0:.0f}s)", flush=True)
        ip[-1] = nz0
        lp[-1] = nz0
        xp[:] = ip
        lp[:] = ip

    manifest = pd.DataFrame({
        "plate": obs["plate"], "well": obs["well"],
        "cell_line_id": obs["cell_line_id"], "drug": obs["drug"],
        "dose_um": obs["dose_um"], "n_cells": obs["n_cells"],
        "is_control": obs["is_control"],
    })
    RECORD.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(RECORD / "tahoe_pseudobulk_manifest.csv", index=False)
    print(f"assemble: wrote {dest} ({dest.stat().st_size/1e9:.1f} GB) "
          f"in {time.time()-t0:.0f}s", flush=True)
    # round-trip check
    back = ad.read_h5ad(dest, backed="r")
    print(f"round-trip: shape={back.shape} obs={list(back.obs.columns)[:6]}...")
    return 0


# ---------------------------------------------------------------------------
# qa: does the aggregation reproduce the authors' DESeq2 pseudobulk contrasts?
# ---------------------------------------------------------------------------
def _de_runs(table):
    """Group a DE shard's rows into contiguous per-contrast runs."""
    import pyarrow.compute as pc
    cols = ["plate", "Cell_Name_Vevo", "drug", "concentration"]
    codes = np.stack([pc.dictionary_encode(table[c]).combine_chunks().indices.to_numpy()
                      for c in cols], axis=1)
    if len(codes) == 0:
        return []
    new = np.any(codes[1:] != codes[:-1], axis=1)
    starts = np.concatenate([[0], np.nonzero(new)[0] + 1])
    ends = np.concatenate([starts[1:], [len(codes)]])
    uniq = [tuple(int(x) for x in codes[s]) for s in starts]
    return list(zip(starts, ends, uniq))


def stage_qa(args) -> int:
    import glob
    import scipy.sparse as sp
    from scipy.stats import spearmanr, pearsonr

    de_home = Path(args.de_home) if args.de_home else DATA_HOME
    de_files = sorted(glob.glob(str(de_home / "metadata" / "pseudobulk_differential_expression" / "*.parquet")))
    if not de_files:
        raise SystemExit(f"no pseudobulk DE shards under {de_home}")
    parts = sorted(PARTS.glob("plate*.npz"), key=lambda p: int(p.stem[5:]))
    if args.plates:
        keep = {p.strip() for p in args.plates.split(",")}
        parts = [p for p in parts if p.stem in keep]
    parts = {p.stem: p for p in parts}

    sample_meta = load_sample_meta()
    gene_meta = load_gene_meta()
    # The authors' DE axis is NOT in metadata order and mixes symbols with
    # Ensembl IDs (observed in shard 240), so alignment must be by name.
    name2col: dict[str, int] = {}
    for i, g in enumerate(gene_meta):
        name2col[g["gene_symbol"]] = i
        if g.get("ensembl_id"):
            name2col.setdefault(g["ensembl_id"], i)
    cvcl2name = load_cell_lines()
    name2cvcl = {}
    for cvcl, name in cvcl2name.items():
        name2cvcl.setdefault(name, cvcl)

    # --- collect contrasts from a sample of DE shards, grouped by plate ---
    rng = np.random.default_rng(0)
    picks = rng.choice(len(de_files), size=min(args.de_shards, len(de_files)), replace=False)
    per_plate: dict[str, list[dict]] = {}
    t0 = time.time()
    for f in sorted(picks):
        t = pq.read_table(de_files[f],
                          columns=["gene_name", "baseMean", "log2FoldChange",
                                   "plate", "Cell_Name_Vevo", "drug", "concentration"])
        names = t["plate"].to_pylist()
        lfc = t["log2FoldChange"].to_numpy(zero_copy_only=False)
        base = t["baseMean"].to_numpy(zero_copy_only=False)
        runs = _de_runs(t)
        for i, (a, b, key) in enumerate(runs):
            if b - a < 60_000:          # contrast straddles the shard edge
                continue
            plate = f"plate{names[a]}"
            if plate not in parts:
                continue
            cell, drug, conc = t["Cell_Name_Vevo"][a].as_py(), t["drug"][a].as_py(), float(t["concentration"][a].as_py())
            cvcl = name2cvcl.get(cell)
            if cvcl is None:
                continue
            gname = t["gene_name"][a:b].to_pylist()
            idx = np.array([name2col.get(n, -1) for n in gname], dtype=np.int64)
            valid = idx >= 0
            per_plate.setdefault(plate, []).append({
                "cvcl": cvcl, "drug": drug, "dose": round(conc, 6),
                "idx": idx[valid],
                "lfc": lfc[a:b][valid].astype(np.float64),
                "base": base[a:b][valid].astype(np.float64),
                "plate": plate,
            })
        del t
    got = sum(len(v) for v in per_plate.values())
    print(f"qa: {len(per_plate)} plates / {got} usable contrasts from {len(picks)} DE shards "
          f"({time.time()-t0:.0f}s)", flush=True)

    # --- per-plate: rebuild cp10k, controls, and compare ---
    results = []
    detect = []
    ctrl_checked = set()
    ctrl_missing = []
    for plate, path in sorted(parts.items(), key=lambda kv: int(kv[0][5:])):
        d = np.load(path, allow_pickle=False)
        wells = d["wells"]
        lines = d["line_ids"]
        obs_well = wells[d["row_well_ix"]]
        obs_line = lines[d["row_line_ix"]]
        plate_w = np.array([sample_meta[w]["drug"] for w in obs_well], dtype=object)
        plate_dose = np.array([sample_meta[w]["dose"] for w in obs_well], dtype=np.float64)
        n_rows = len(obs_well)
        nnz_row = np.diff(d["indptr"])
        detect.append(nnz_row)
        M = sp.csr_matrix((d["data"].astype(np.float64), d["indices"], d["indptr"]),
                          shape=(n_rows, N_GENES))

        is_ctl = plate_w == CONTROL_DRUG
        cmean = {}
        for ln in np.unique(obs_line):
            sel = is_ctl & (obs_line == ln)
            if sel.any():
                cmean[ln] = np.asarray(M[sel].mean(axis=0)).ravel()
        treated = {}
        for i in np.nonzero(~is_ctl)[0]:
            treated.setdefault((obs_line[i], plate_w[i], round(float(plate_dose[i]), 6)), i)
        for ln in np.unique(obs_line):
            if (not is_ctl[obs_line == ln].any()) and ((~is_ctl & (obs_line == ln)).any()):
                ctrl_missing.append(f"{plate}/{ln}")
            ctrl_checked.add((plate, ln))

        block = per_plate.get(plate, [])
        for c in block:
            row = treated.get((c["cvcl"], c["drug"], c["dose"]))
            ctrl = cmean.get(c["cvcl"])
            if row is None or ctrl is None:
                continue
            ours = M.getrow(row).toarray().ravel()
            base, lfc, idx = c["base"], c["lfc"], c["idx"]
            ours_lfc = np.log2((ours[idx] + 1.0) / (ctrl[idx] + 1.0))
            rec = {"plate": plate, "cell_line": c["cvcl"], "drug": c["drug"], "dose": c["dose"],
                   "n_cells_treated": int(d["row_ncell"][row]),
                   "n_cells_control": int(d["row_ncell"][sel].sum())}
            finite = np.isfinite(lfc)
            for thr in (300, 500, 800):
                m = finite & (base > thr)
                if m.sum() >= 50:
                    rec[f"pearson_{thr}"] = float(pearsonr(ours_lfc[m], lfc[m]).statistic)
                    rec[f"spearman_{thr}"] = float(spearmanr(ours_lfc[m], lfc[m]).statistic)
                    rec[f"n_{thr}"] = int(m.sum())
            m_all = finite & (base > 0)
            if m_all.sum() >= 50:
                rec["spearman_all"] = float(spearmanr(ours_lfc[m_all], lfc[m_all]).statistic)
            results.append(rec)
        for ln in np.unique(obs_line):
            if (is_ctl & (obs_line == ln)).any():
                ctrl_checked.add((plate, ln))
        del M
    rhos = np.array([r["spearman_all"] for r in results if "spearman_all" in r])
    p500 = np.array([r["pearson_500"] for r in results if "pearson_500" in r])
    p800 = np.array([r["pearson_800"] for r in results if "pearson_800" in r])
    n800 = np.array([r["n_800"] for r in results if "n_800" in r])
    s500 = np.array([r["spearman_500"] for r in results if "spearman_500" in r])
    p300 = np.array([r["pearson_300"] for r in results if "pearson_300" in r])
    n500 = np.array([r["n_500"] for r in results if "n_500" in r])
    det = np.concatenate(detect) if detect else np.array([0])
    all_ncells = np.concatenate([np.load(p, allow_pickle=False)["row_ncell"] for p in parts.values()])
    summary = {
        "checked_at": now(),
        "contrasts_compared": len(results),
        "pearson_median_base800": float(np.median(p800)) if len(p800) else None,
        "pearson_p10_base800": float(np.percentile(p800, 10)) if len(p800) else None,
        "genes_per_contrast_median_base800": int(np.median(n800)) if len(n800) else 0,
        "pearson_median_base500": float(np.median(p500)) if len(p500) else None,
        "pearson_p10_base500": float(np.percentile(p500, 10)) if len(p500) else None,
        "spearman_median_base500": float(np.median(s500)) if len(s500) else None,
        "pearson_median_base300": float(np.median(p300)) if len(p300) else None,
        "genes_per_contrast_median_base500": int(np.median(n500)) if len(n500) else 0,
        "spearman_median_all_genes": float(np.median(rhos)) if len(rhos) else None,
        "gate_pearson_ge_0.60_base500": bool(len(p500) >= 50 and np.median(p500) >= 0.60),
        "gate_pearson_ge_0.70_base500": bool(len(p500) >= 50 and np.median(p500) >= 0.70),
        "detected_genes_median": float(np.median(det)),
        "gate_detected_ge_12000": bool(np.median(det) >= 12_000),
        "profiles": int(len(all_ncells)),
        "profiles_ge_100_cells": int((all_ncells >= 100).sum()),
        "control_missing": sorted(set(ctrl_missing)),
        "gate_controls_present": not ctrl_missing,
    }
    RECORD.mkdir(parents=True, exist_ok=True)
    (RECORD / "tahoe_pseudobulk_qa.json").write_text(json.dumps(summary, indent=2))
    import csv as _csv
    keys = sorted({k for r in results for k in r})
    with open(RECORD / "tahoe_pseudobulk_qa_contrasts.csv", "w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in results:
            w.writerow(r)
    lines_md = ["# Tahoe pseudobulk corpus — QA", "",
                f"- profiles: **{summary['profiles']:,}** "
                f"(>=100 cells: **{summary['profiles_ge_100_cells']:,}**, "
                f"{100*summary['profiles_ge_100_cells']/max(summary['profiles'],1):.1f}%)",
                f"- median detected genes per profile: **{summary['detected_genes_median']:.0f}** "
                f"(gate >=12,000: {'PASS' if summary['gate_detected_ge_12000'] else 'FAIL'})",
                f"- LFC vs authors' DESeq2 DE over **{summary['contrasts_compared']}** contrasts; "
                f"genes with baseMean>800 (median {summary['genes_per_contrast_median_base800']}/contrast): "
                f"Pearson {summary['pearson_median_base800']:.3f} "
                f"(p10 {summary['pearson_p10_base800']:.3f}) - reported, not gated",
                f"- baseMean>500 ({summary['genes_per_contrast_median_base500']}/contrast): "
                f"Pearson {summary['pearson_median_base500']:.3f}, "
                f"Spearman {summary['spearman_median_base500']:.3f} - gate >=0.60: "
                f"{'PASS' if summary['gate_pearson_ge_0.60_base500'] else 'FAIL'}",
                f"- baseMean>300: Pearson {summary['pearson_median_base300']:.3f}; "
                f"Spearman over all genes: {summary['spearman_median_all_genes']:.3f} "
                f"(reported, not gated - near-null genes dominate ranks)",
                f"- (plate, cell_line) pairs with a vehicle control: "
                f"**{len(ctrl_checked)}** checked, "
                f"{'all present' if not ctrl_missing else f'{len(ctrl_missing)} missing'}",
                "",
                "Method: our log2((treated counts + 1) / (same-plate same-cell-line DMSO mean + 1)) "
                "vs DESeq2 `log2FoldChange`, per-contrast Pearson/Spearman restricted to the "
                "authors' well-measured genes (baseMean threshold). The full-gene Spearman is "
                "reported for context only: on genes with no real effect both estimates are noise, "
                "so rank agreement there is uninformative (measured ~0.2).",
                ""]
    (RECORD / "tahoe_pseudobulk_qa_report.md").write_text("\n".join(lines_md))
    print(json.dumps({k: v for k, v in summary.items() if k != "control_missing"}, indent=2))
    return 0 if (summary["gate_pearson_ge_0.60_base500"]
                 and summary["gate_detected_ge_12000"]
                 and summary["gate_controls_present"]) else 1


def check_layout() -> dict:
    """Report what the resolver found: repo root, Tahoe input, corpus output."""
    return {
        "repo_root": str(ROOT),
        "tahoe_source": str(DATA_HOME),
        "tahoe_source_present": DATA_HOME.is_dir(),
        "expression_shards_on_disk": len(all_shards()),
        "corpus_output": str(OUT),
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=f"resolved from {Path(__file__).resolve()}:\n"
               f"  repo root   {ROOT}\n"
               f"  tahoe input {DATA_HOME}\n"
               f"  output      {OUT}")
    p.add_argument("stage", nargs="?",
                   choices=["scan", "aggregate", "assemble", "qa"],
                   help="scan | aggregate | assemble | qa")
    p.add_argument("--check", action="store_true",
                   help="print the resolved repository root, Tahoe source and output dir, then exit")
    p.add_argument("--plates", default=None)
    p.add_argument("--shards", type=int, default=None)
    p.add_argument("--force", action="store_true")
    p.add_argument("--de-home", default=None, help="data home holding metadata/pseudobulk_differential_expression")
    p.add_argument("--de-shards", type=int, default=30)
    args = p.parse_args(argv)
    if args.check:
        print(json.dumps(check_layout(), indent=2))
        return 0 if DATA_HOME.is_dir() else 1
    if args.stage is None:
        p.print_help(sys.stderr)
        return 2
    if not DATA_HOME.is_dir():
        raise SystemExit(f"Tahoe source missing: {DATA_HOME} "
                         f"(download/verify it with pert_data/Tahoe-100M/script/)")
    return {"scan": stage_scan, "aggregate": stage_aggregate,
            "assemble": stage_assemble, "qa": stage_qa}[args.stage](args)


if __name__ == "__main__":
    sys.exit(main())
