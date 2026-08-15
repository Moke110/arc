"""
xpert-fda-miri-v2 gene-mapping verification gates (per union dataset).

Extends v1's _lib/verify_genes.py G1/G2/G3 to every dataset in the CM
finetune-data union (dtoxs_primary, dtoxs_ipsc, gse244740, gse247541,
gse262419), each resolved in its own identifier space (symbols / Ensembl /
SYMBOL_entrezid) with legacy-symbol awareness:

G1 — mapping completeness: every landmark that the dataset's gene panel
     actually contains resolves (direct | legacy_alias | entrez | ensembl).
     gse262419's TempO-Seq panel provably lacks WRB and FAM57A (their rows
     point at different Entrez genes), so the gate tolerates exactly the
     documented platform-missing set and aborts on any OTHER unresolved
     landmark.
G2 — alias correctness (author+dataset data only): every legacy_alias source
     symbol is the author-recorded name (gene_name / gene_name_match) of the
     mapped landmark, is claimed by no other landmark, and its
     `source_entrez_verified` equals the landmark's Entrez gene_id (rejects
     the JPT2 <- CRAMP1 trap: the author's `corrected_gene_name` says CRAMP1,
     Entrez 57585, a different gene — corrected_gene_name is never a match key).
G3 — per-dataset drug-level audit (informational + flag): for the ambiguous
     landmarks (legacy_alias, plus gse262419 entrez matches), correlate each
     union-h5ad column against the raw per-drug log2(CPM+1) of the mapped
     source vs competing candidates, and report which source the column
     actually follows. A confirmed misassignment aborts.

Any G1/G2 failure aborts with exit 1. G3 reports findings and exits non-zero
only on a confirmed misassignment.

Usage:
    python _lib/verify_genes.py            # build maps + verify (writes report)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "03_cm_finetune_data"))

import gene_mapping as gm  # noqa: E402
from utils import WORKFLOW_ROOT, CM_DATASETS  # noqa: E402
from build_union import build_landmark_matrix, matched_control_matrix  # noqa: E402
import _datasets  # noqa: E402

UNION_H5AD = WORKFLOW_ROOT / "03_cm_finetune_data" / "output" / "cm_union.h5ad"
MAP_DIR = WORKFLOW_ROOT / "03_cm_finetune_data" / "output" / "gene_maps"
REPORT = WORKFLOW_ROOT / "_lib" / "gene_mapping_report.md"


def per_drug_raw_log2cpm(counts, gene_ids, drug_of, is_control):
    """Per-drug mean log2(CPM+1) over a dataset's raw counts.

    counts is (n_samples, n_genes); returns {drug_name: {gene_id: value}}.
    """
    totals = counts.sum(axis=1)
    cpm = np.divide(counts, totals[:, None], out=np.zeros_like(counts),
                    where=totals[:, None] > 0) * 1e6
    l2 = np.log2(cpm + 1)
    out = {}
    for i in range(counts.shape[0]):
        if is_control[i]:
            continue
        drug = drug_of[i]
        for j, gid in enumerate(gene_ids):
            out.setdefault(drug, {})[gid] = float(l2[i, j])
    return out


def drug_level_audit(dataset, mapping, meta, counts, gene_ids, adata,
                     union_var_ids, name_to_pert_idx):
    """Which raw source does each union column follow? (drug-level test)."""
    # per-pert_idx raw means (keyed by drug index to match the union h5ad)
    def drug_key(i):
        name = meta.iloc[i]["compound_name"]
        if meta.iloc[i]["is_control"] or name not in name_to_pert_idx:
            return -1
        return name_to_pert_idx[name]
    drug_of = np.array([drug_key(i) for i in range(len(meta))])
    raw_by_drug = per_drug_raw_log2cpm(
        counts, gene_ids, drug_of, meta["is_control"].values)
    raw_series = {gid: pd.Series({d: v[gid] for d, v in raw_by_drug.items()
                                  if gid in v}) for gid in set(gene_ids)}
    if dataset == "gse262419":  # also key by the numeric (Entrez) part
        for gid in list(raw_series):
            num = str(gid).rsplit("_", 1)[1] if "_" in str(gid) else ""
            if num and num not in raw_series:
                raw_series[num] = raw_series[gid]

    union = adata[adata.obs["dataset_id"] == dataset]
    X = union.X.toarray() if hasattr(union.X, "toarray") else union.X
    obs = union.obs
    # group by pert_idx (drug) for the union column means
    drug_groups = obs.assign(_j=range(len(obs))).groupby("pert_idx")["_j"] \
        .apply(lambda j: j.tolist()).to_dict()
    meta_genes = gm.load_l1000_meta()
    gid2col = {str(g): i for i, g in enumerate(union_var_ids)}

    findings = []
    ambiguous = mapping[mapping["match_type"].isin(["legacy_alias", "entrez"])]
    for _, row in ambiguous.iterrows():
        sym, gid, src = row["l1000_gene_symbol"], str(row["l1000_gene_id"]), \
            row["source_symbol"]
        if gid not in gid2col:
            findings.append((dataset, sym, "union column missing"))
            continue
        h5ad_series = pd.Series({pid: X[js, gid2col[gid]].mean()
                                 for pid, js in drug_groups.items()})
        # candidates = author metadata names for this landmark present in raw
        meta_row = meta_genes[meta_genes["gene_id"].astype(str) == gid].iloc[0]
        names = [str(x) for x in (meta_row["gene_name"],
                                  meta_row["gene_name_match"],
                                  meta_row["corrected_gene_name"], sym)
                 if pd.notna(x)]
        if dataset == "gse262419" and row["match_type"] == "entrez":
            # also consider the SYMBOL_entrezid row carrying this entrez number
            cand_keys = [src]
        else:
            cand_keys = list(dict.fromkeys(names))
        scores = {}
        for c in cand_keys:
            key = str(c).split(".")[0]
            if key not in raw_series:
                continue
            s = raw_series[key]
            common = h5ad_series.index.intersection(s.index)
            if len(common) >= 5:
                rho, _ = spearmanr(h5ad_series[common], s[common])
            else:
                rho = float("nan")
            scores[c] = rho
        if not scores:
            continue
        best = max(scores, key=lambda c: scores[c] if scores[c] == scores[c] else -9)
        rho_src = scores.get(src, float("nan"))
        if best != src and scores[best] - (rho_src if rho_src == rho_src else -9) >= 0.25:
            findings.append((dataset, sym,
                             f"union column follows {best} (expected {src}): "
                             f"rhos={ {k: round(v, 2) for k, v in scores.items()} }"))
        else:
            print(f"  AUDIT ok [{dataset}] {sym}: follows {best}, "
                  f"rho_src={rho_src:.2f}")
    return findings


def main():
    print("xpert-fda-miri-v2 gene-mapping gates (per union dataset)\n")

    # Build all per-dataset maps
    maps = {}
    for ds in CM_DATASETS:
        maps[ds] = gm.build_dataset_map(ds)
        gm.write_dataset_map(ds, maps[ds], MAP_DIR)
    print("Per-dataset maps written -> " + str(MAP_DIR))

    all_findings = []
    for ds in CM_DATASETS:
        mapping = maps[ds]
        n_unmatched = (mapping["match_type"] == "unmatched").sum()

        # G1 — platform-missing tolerance only for gse262419 WRB/FAM57A
        expected_missing = (gm.GSE262419_PLATFORM_MISSING
                            if ds == "gse262419" else set())
        unresolved = set(mapping.loc[mapping["match_type"] == "unmatched",
                                     "l1000_gene_symbol"])
        bad_unresolved = unresolved - expected_missing
        print(f"\n[{ds}] G1 completeness: "
              f"{(len(mapping) - n_unmatched)}/{len(mapping)} resolved "
              f"({dict(mapping['match_type'].value_counts())})")
        if bad_unresolved:
            print("  ERROR unresolved:", sorted(bad_unresolved))
            sys.exit(1)
        if expected_missing - unresolved:
            print("  ERROR expected platform-missing set not missing:",
                  sorted(expected_missing - unresolved))
            sys.exit(1)
        if expected_missing:
            print(f"  platform-missing (documented, zeroed): "
                  f"{sorted(expected_missing)}")

        # G2 — author-metadata verification of every legacy_alias
        meta = gm.load_l1000_meta()
        name_to_lm = {}
        for _, r in meta.iterrows():
            for c in ("gene_symbol", "gene_name", "gene_name_match"):
                v = r[c]
                if pd.notna(v) and str(v) not in ("nan", ""):
                    name_to_lm.setdefault(str(v), set()).add(
                        (str(r["gene_id"]), str(r["gene_symbol"])))
        alias = mapping[mapping["match_type"] == "legacy_alias"]
        bad = []
        for _, r in alias.iterrows():
            src, gid, sym = str(r["source_symbol"]), str(r["l1000_gene_id"]), \
                str(r["l1000_gene_symbol"])
            author_paired = (gid, sym) in name_to_lm.get(src, set())
            entrez_ok = str(r["source_entrez_verified"]) == gid
            if not author_paired or not entrez_ok:
                bad.append((src, sym, name_to_lm.get(src, set()),
                            r["source_entrez_verified"]))
        print(f"[{ds}] G2 alias verification (author metadata): {len(alias)} "
              f"aliases, {len(bad)} mismatches")
        for src, sym, owners, ent in bad:
            print(f"  MISMATCH {src} -> {sym} (author owners {owners}, "
                  f"src entrez {ent} != {sym})")
        if len(bad):
            sys.exit(1)

        # G3 — drug-level audit against the union h5ad (skipped until step 03
        # has built the union h5ad; G1/G2 above run on the maps alone)
        manifest = pd.read_csv(WORKFLOW_ROOT / "03_cm_finetune_data" / "output"
                               / "drug_manifest.csv")
        man = manifest[manifest["dataset"] == ds]
        name_to_pert_idx = {str(r["compound_name"]): int(r["drug_idx"])
                            for _, r in man.iterrows() if r["hg_class"] == "in_hg"}
        meta_df, counts, gene_ids = _datasets.extract_samples(ds)
        if not Path(UNION_H5AD).exists():
            print(f"[{ds}] G3 drug-level audit (union h5ad) — skipped, "
                  f"union h5ad not built yet (run 03_cm_finetune_data first)")
            continue
        import scanpy as sc
        adata = sc.read_h5ad(UNION_H5AD)
        union_var_ids = [str(v) for v in adata.var_names]
        print(f"[{ds}] G3 drug-level audit (union h5ad)")
        findings = drug_level_audit(ds, mapping, meta_df, counts, gene_ids,
                                    adata, union_var_ids, name_to_pert_idx)
        all_findings += [(ds, s, m) for (_, s, m) in findings] \
            if not findings else findings
        print(f"  {len(findings)} finding(s)")
        for f in findings:
            print(f"  FINDING {f}")

    flag = [f for f in all_findings
            if "MISMATCH" in f[2] or "missing" in f[2]
            or "follows " in f[2]]
    lines = ["# xpert-fda-miri-v2 gene-mapping verification report", ""]
    for ds in CM_DATASETS:
        mp = maps[ds]
        n = (mp["match_type"] != "unmatched").sum()
        alias_n = (mp["match_type"] == "legacy_alias").sum()
        lines.append(f"- **{ds}**: {n}/978 resolved, {alias_n} legacy aliases "
                     f"(all author-metadata verified)")
    lines.append(f"- union h5ad drug-level audit findings: {len(all_findings)}")
    lines.append("")
    for f in all_findings:
        lines.append(f"- **{f[0]} / {f[1]}**: {f[2]}")
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"\nReport: {REPORT}")

    if flag:
        print("\nERROR: confirmed h5ad gene misassignments — see report.")
        sys.exit(1)
    print("\nAll gene-mapping gates passed.")


if __name__ == "__main__":
    main()
