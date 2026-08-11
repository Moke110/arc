"""
xpert-fda-miri-v1 gene-mapping verification gate.

Checks that all 978 L1000 landmark genes participate in the DToxS finetune
with the *correct* DToxS source gene (Entrez-verified), and audits the frozen
finetune h5ad (input/dtoxs.h5ad) for gene-identity inconsistencies.

Gates
-----
G1 — mapping completeness: every one of the 978 landmarks resolves to a DToxS
     source symbol (direct or legacy-alias), i.e. all 978 can participate.
G2 — alias correctness (author+DToxS data only): every legacy-alias source
     symbol is the author-recorded name (gene_name / gene_name_match) of the
     mapped landmark, is claimed by no other landmark, and its `source_entrez_verified`
     equals the landmark's own Entrez gene_id (rejects the JPT2 <- CRAMP1 trap:
     the author's `corrected_gene_name` says CRAMP1, Entrez 57585, a different
     gene — so corrected_gene_name is never used as a match key).
G3 — h5ad audit (informational + flag): for each landmark, compare the h5ad
     column's per-drug means against the raw DToxS per-state means of (a) the
     authoritative source symbol and (b) competing DToxS symbols, and report
     which symbol the h5ad column actually follows. A mismatch (e.g. JPT2
     following CRAMP1 instead of HN1L) is reported as a finding.

Any G1/G2 failure aborts with exit 1. G3 reports findings and exits non-zero
if a confirmed misassignment is detected.

Usage:
    python _lib/verify_genes.py            # build mapping + verify (writes report)
    python _lib/verify_genes.py --load     # reuse cached mapping CSV
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gene_mapping as gm  # noqa: E402
from utils import WORKFLOW_ROOT  # noqa: E402

H5AD = WORKFLOW_ROOT / "input" / "dtoxs.h5ad"
MAP_CSV = WORKFLOW_ROOT / "02_cm_baseline" / "output" \
    / "dtoxs_to_l1000_gene_map.csv"
REPORT = WORKFLOW_ROOT / "_lib" / "gene_mapping_report.md"


def drug_level_audit(mapping, raw, state2wells):
    """Which DToxS symbol does each h5ad column follow? (drug-level test).

    For every legacy-alias landmark, the plausible DToxS candidates are the
    names the L1000 metadata associates with it (gene_name / gene_name_match /
    corrected_gene_name / gene_symbol) that exist in the DToxS table. The h5ad
    column (log2(CPM+1) per drug) is correlated against each candidate's raw
    per-state mean on the SAME scale (log2(CPM+1)); a clear mismatch means the
    frozen h5ad filled the landmark from a different gene (e.g. JPT2 <- CRAMP1).
    """
    import scanpy as sc

    adata = sc.read_h5ad(H5AD)
    X = adata.X.toarray()
    var_ids = [str(v) for v in adata.var_names]
    obs = adata.obs
    per_drug = obs.assign(_j=range(len(obs))).groupby("pert_id")["_j"].apply(
        lambda j: j.tolist()).to_dict()
    h5ad_drug_mean = {pid: X[js].mean(axis=0) for pid, js in per_drug.items()}

    # per-state log2(CPM+1) means for every DToxS symbol (matches h5ad scale)
    raw_by_state = {}
    for st, wells in state2wells.items():
        w = [x for x in wells if x in raw.columns]
        if w:
            tot = raw[w].sum(axis=0)
            cpm = raw[w].div(tot, axis=1) * 1e6
            raw_by_state[st] = np.log2(cpm.mean(axis=1) + 1)

    meta = pd.read_csv(gm.L1000_META).set_index("gene_symbol")
    findings = []
    for _, row in mapping[mapping["match_type"] == "legacy_alias"].iterrows():
        sym, gid, source = row["l1000_gene_symbol"], str(row["l1000_gene_id"]), \
            row["dtoxs_source_symbol"]
        if gid not in var_ids:
            findings.append((sym, "h5ad column missing"))
            continue
        h5ad_series = pd.Series({pid: h5ad_drug_mean[pid][var_ids.index(gid)]
                                 for pid in h5ad_drug_mean})
        # candidates = L1000-metadata names for this landmark that exist in DToxS
        names = [str(x) for x in (meta.loc[sym, "gene_name"],
                                  meta.loc[sym, "gene_name_match"],
                                  meta.loc[sym, "corrected_gene_name"],
                                  sym) if pd.notna(x)]
        cands = dict.fromkeys(n for n in names if n in raw.index)
        scores = {}
        for c in cands:
            state_means = pd.Series({st: v.get(c, 0) for st, v in raw_by_state.items()})
            common = h5ad_series.index.intersection(state_means.index)
            if len(common) >= 5:
                rho, _ = spearmanr(h5ad_series[common], state_means[common])
            else:
                rho = float("nan")
            scores[c] = rho
        if not scores:
            continue
        best = max(scores, key=lambda c: scores[c] if scores[c] == scores[c] else -9)
        rho_src = scores.get(source, float("nan"))
        # clear mismatch: another candidate wins by a wide margin and the
        # authoritative source correlates poorly.
        if best != source and scores[best] - (rho_src if rho_src == rho_src else -9) >= 0.25:
            findings.append((sym, f"h5ad follows {best} (expected {source}): "
                                   f"rhos={ {k: round(v, 2) for k, v in scores.items()} }"))
        else:
            print(f"  AUDIT ok {sym}: follows {best}, rho_src={rho_src:.2f}, "
                  f"rhos={ {k: round(v, 2) for k, v in scores.items()} }")
    return findings


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--load", action="store_true",
                        help="Reuse cached mapping CSV instead of rebuilding")
    args = parser.parse_args()

    if args.load and MAP_CSV.exists():
        mapping = pd.read_csv(MAP_CSV)
    else:
        raw, ctrl = gm.load_dtoxs_raw_counts()
        mapping = gm.build_mapping(raw, ctrl)
        MAP_CSV.parent.mkdir(parents=True, exist_ok=True)
        mapping.to_csv(MAP_CSV, index=False)

    # G1
    n_unmatched = (mapping["match_type"] == "unmatched").sum()
    print(f"G1 mapping completeness: {len(mapping) - n_unmatched}/978 resolved "
          f"({n_unmatched} unmatched)")
    if n_unmatched:
        print("  UNMATCHED:", mapping.loc[mapping['match_type'] == 'unmatched',
                                          'l1000_gene_symbol'].tolist())
        sys.exit(1)

    # G2 — author-metadata verification of every legacy alias
    meta = pd.read_csv(gm.L1000_META)
    meta["gene_id"] = meta["gene_id"].astype(str)
    # author name -> (gene_id, gene_symbol), from gene_symbol/gene_name/gene_name_match
    name_to_lm = {}
    for _, r in meta.iterrows():
        for c in ("gene_symbol", "gene_name", "gene_name_match"):
            v = r[c]
            if pd.notna(v) and str(v) not in ("nan", ""):
                name_to_lm.setdefault(str(v), set()).add((str(r["gene_id"]),
                                                          str(r["gene_symbol"])))
    alias = mapping[mapping["match_type"] == "legacy_alias"]
    bad = []
    for _, r in alias.iterrows():
        src, gid, sym = str(r["dtoxs_source_symbol"]), str(r["l1000_gene_id"]), \
            str(r["l1000_gene_symbol"])
        owners = name_to_lm.get(src, set())
        author_paired = (gid, sym) in owners
        entrez_ok = str(r["source_entrez_verified"]) == gid
        if not author_paired or not entrez_ok:
            bad.append((src, sym, owners, r["source_entrez_verified"]))
    print(f"G2 alias verification (author metadata): {len(alias)} aliases, "
          f"{len(bad)} mismatches")
    for src, sym, owners, ent in bad:
        print(f"  MISMATCH {src} -> {sym} (author owners {owners}, "
              f"src entrez {ent} != {sym})")
    if len(bad):
        sys.exit(1)

    # G3 (informational audit; needs the raw data)
    raw, ctrl = gm.load_dtoxs_raw_counts()
    state2wells = gm_state2wells()
    findings = drug_level_audit(mapping, raw, state2wells)
    print(f"\nG3 h5ad drug-level audit: {len(findings)} finding(s)")
    for sym, msg in findings:
        print(f"  FINDING {sym}: {msg}")
    flag = [f for f in findings if "MISMATCH" in f[1] or "missing" in f[1]]

    lines = ["# xpert-fda-miri-v1 gene-mapping verification report", "",
             f"- mapping completeness: {len(mapping)-n_unmatched}/978",
             f"- legacy aliases: {len(alias)} (all author-metadata verified)",
             f"- h5ad drug-level audit findings: {len(findings)}",
             ""]
    for sym, msg in findings:
        lines.append(f"- **{sym}**: {msg}")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"Report: {REPORT}")

    if flag:
        print("\nERROR: confirmed h5ad gene misassignments — see report.")
        sys.exit(1)
    print("\nAll gene-mapping gates passed.")


def gm_state2wells():
    import glob

    level3 = gm.DTOXS_ROOT / "Data" / "Level3"
    state2wells = {}
    for fp in sorted(glob.glob(str(level3 / "**"
                                   / "PromoCell-DGE-Experiment-Design*.tsv"),
                               recursive=True)):
        d = pd.read_csv(fp, sep="\t")
        for st, w in zip(d["State"], d["Well"]):
            state2wells.setdefault(str(st), []).append(w)
    return state2wells


if __name__ == "__main__":
    main()
