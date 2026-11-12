# Gene-axis 08 — review-fix acceptance record (#183)

Scope: the blocking quality gate for ticket **#183** (parent **#177**; the
documentation/acceptance surface is #184). An adversarial code review of the
Stage 2 gene-axis evidence found six real numerical bugs — all index-space errors
in which consumers read gene-indexed arrays with the row positions of the tracked
`fused_rank.csv`, which is sorted by `fused_rank` — plus three provenance /
documentation findings. This record is the gate evidence for the repair round. It
re-runs and re-derives the evidence and records the documentation, path-audit
and inventory state at acceptance.

This is **decision preparation, not the axis decision** (ADR 0065). It chooses no
`N`, per-source quota, direction rule or animal weight (#173), builds no corpus,
trains no model, publishes no `genes.csv` (#174), and opens no truth or Heliste
label file. The six evidence answers below come from
`stage2/01_bases/output/gene-axis/`; the tracked resolution asset is
`docs/research/stage2-stage3-gene-axis-rank-aggregation.md`.

## What this round changed

The six numerical findings (details and old → new headline values in
`stage2/01_bases/log/gene-axis/record.md → "Review fixes at the root"`):

1. `tenv_non_hvg_curve.csv` was the non-HVG share of the first N genes in
   `gene_index` order, not the fused top-N.
2. `overlap.csv` `fused_rank_min/median/max` were `gene_index + 1`, not fused
   ranks.
3. `direction_summary.csv` was the first N genes in `gene_index` order, and the
   `|log2FC|`-weighted row duplicated the symmetric row by construction; it now
   reports the true fused top-N and a genuinely different ranking
   (`max_abs_log2fc`: descending maximum cross-source `|log2FC|`, fused-rank
   tie-break).
4. The #175-Q3 DEG-count tie-break is by **ascending fused rank** (the tracked
   values already came from that ranking); the prose wrongly said `gene_index`
   and now states the fused-rank rule.
5. `hvg_share` was identically `2000 / n_axis_expression`; it is now
   `n_hvg_on_axis / n_hvg_total` over the source's whole measured transcriptome,
   with `n_hvg`, `n_hvg_total`, `n_hvg_on_axis` and `hvg_density` reported.
6. `overlap.csv` now adds tenv-HVG × source-HVG at K = 1000/2000/5000 and the
   multi-way intersection rows report the intersection size and true fused ranks
   with a `note` instead of a fake `jaccard = 1.0` (`n_a`/`n_b`/`n_union`/
   `jaccard` blank).

Also fixed: the DEG-CSV sha256 is recorded in
`disease_data/GSE141910/local/deg/GSE141910_deg.json` and
`disease_data/GSE116250/local/deg/GSE116250_deg.json` (the three other items
already carried it; the CSVs are unchanged), the `GLOSSARY.md` "Source rank"
definition now matches the product, and
`stage2/01_bases/execution/build_mouse_orthologs.R` writes a candidate under
`execution/gene_axis/` and replaces the accepted table only after a sha256 check
(`--args --accept` for a deliberate re-pin).

`fused_rank.csv` and `inflection.json` are byte-identical to the pre-fix product:
the fusion and the sorted-curve inflection were never wrong. Everything
downstream of the four fixed consumers was regenerated.

## Source → item mapping

The ten fusion sources of `output/gene-axis/manifest.json → constants.fused_sources`,
each projected onto the unchanged frozen 18,164-gene axis
(`stage2/01_bases/output/gene-identity/genes.csv`, sha256
`5877a3cd46a137f739ca985c5357578b1e77d11c78ce5acf0af93249592dbe57`).

| source id | item / accession | role | interface read by the build |
|---|---|---|---|
| `tenv` | published tenv descriptive HVG ranking | variability reference; **not recomputed** | `stage2/01_bases/output/tenv-hvg/descriptive.csv` |
| `gse141910` | `disease_data/GSE141910/` | human LV HF 200/166; HF primary anchor | `local/deg/GSE141910_deg.csv` + `local/expr/GSE141910_log_expr.tsv.gz` |
| `gse116250` | `disease_data/GSE116250/` | human LV HF 50/14 incl. ICM; cheapest HF anchor | `local/deg/GSE116250_deg.csv` + `local/expr/…` |
| `gse243674` | `disease_data/GSE243674/` | human hiPSC-CM anthracyclines + vehicle, 72 wells | `local/deg/GSE243674_deg.csv` + `local/expr/…` |
| `gse128188` | `disease_data/GSE128188/` | human AF 5/5, paired LA+RA; >5% AF category | `local/deg/GSE128188_deg.csv` + `local/expr/…` |
| `gse282618` | `disease_data/GSE282618/` | human HHD 3 vs 2 subjects, exploratory | `local/deg/GSE282618_deg.csv` + `local/expr/…` |
| `gse160516` | `disease_data/GSE160516/` | mouse MIRI, Affymetrix Clariom S; DEG-only | `source/GSE160516_deg_results.csv` |
| `gse108940` | `disease_data/GSE108940/` | mouse MIRI, GPL7202, corrected #176; DEG-only | `source/GSE108940_WT_deg_results.csv` |
| `dtoxs_ipsc` | `pert_data/DToxS-iPSC/` (GSE217421 tables) | local anthracycline evidence; DEG-only | `source/GSE217421_Conv-DEG-Calc-FDR-0.1.tar.gz` |
| `gse244740` | `pert_data/GSE244740/` | local ScreenSeq hiPSC-CM cardiotox panel; DEG-only | `source/GSE244740_processed_data_counts.txt.gz` + series matrix |

Listed only, never fused: `disease_data/GSE137293/` (pig scar vs border — an
injury-severity gradient) and `pert_data/LINCS-LDS-1593/` (L1000-978 coverage
note only). Both appear in `source_summary.csv` with `supplementary_only = True`.

The animal mapping chain is `output/gene-axis/mouse_human_orthologs.csv`
(homologene `mouse2human`, `inTax=10090`), generated by
`stage2/01_bases/execution/build_mouse_orthologs.R` on the pinned `arc` R. In
this round the generator's candidate was byte-identical to the accepted table
(sha256 `cc22358e…`), so the accepted path was not replaced.

## Six §C answers (real numbers and pointers)

1. **Per-source HVG share, DEG up/down counts and shares, axis coverage** —
   `output/gene-axis/source_summary.csv` (12 rows, columns
   `n_axis_measured, axis_coverage_share, n_hvg, n_hvg_total, n_hvg_on_axis,
   hvg_share, hvg_density, deg_up, deg_down, deg_up_share, deg_down_share,
   deg_share`). At the declared reporting size `K=2000`: axis coverage 0.6803
   (GSE108940) – 1.0000 (tenv, GSE244740); `hvg_share` (share of the source's own
   top-2,000 HVG set on the axis) 0.1465 (GSE116250) – 0.9375 (GSE128188) for the
   expression sources, tenv 1.0000 by construction; `hvg_density` 0.0162
   (GSE116250) – 0.1389 (GSE128188), tenv 0.1101. DEG up/down counts range from
   0/1 (GSE108940) to 8,376/7,076 (DToxS-iPSC). GSE141910 4,538/5,543;
   GSE116250 4,936/2,806; GSE243674 6,248/6,250; GSE128188 234/382; GSE282618
   35/53; GSE160516 437/375; GSE244740 1,137/7,246. Full table in the resolution
   asset and the record.
2. **Per-source and fused rank distributions** — `output/gene-axis/fused_rank.csv`
   (18,164 rows; byte-identical to the pre-fix product). Every source column is a
   complete 1…18,164 rank (unmeasured = `n_measured + 1`). `rank_sum` spans
   **0.606225418233 … 9.85504691298**; `rank_sum_animal_half`
   **0.565792831323 … 8.85504691298**; `n_sources_measured` is **2–10** for every
   gene; `n_sources_deg_significant` is **0–7** (16,746 genes ≥ 1, 15,231 ≥ 2).
3. **Source overlap** — `output/gene-axis/overlap.csv` (81 rows: 27 tenv-HVG ×
   DEG at K = 1000/2000/5000, 15 tenv-HVG × source-HVG at K = 1000/2000/5000,
   36 DEG pairs, 3 multi-way intersections). tenv-HVG(2000) ∩ source DEG:
   972 / 735 / 1,192 / 59 / 10 / 82 / 0 / 1,522 / 666 for GSE141910 / GSE116250 /
   GSE243674 / GSE128188 / GSE282618 / GSE160516 / GSE108940 / DToxS-iPSC /
   GSE244740. tenv-HVG(2000) ∩ source HVG(2000): **155 / 29 / 254 / 234 / 180**
   (Jaccard 0.0467 / 0.0128 / 0.0721 / 0.0643 / 0.0558). HF pair
   GSE141910 ∩ GSE116250 = **5,320** genes, Jaccard **0.4255**, median true fused
   rank **4,954.5**; AF/HHD pair GSE128188 ∩ GSE282618 = **3** genes, Jaccard
   **0.0043**; the four HF×remodelling cross pairs average Jaccard **0.0233**.
   All-nine-DEG intersection empty; the eight-source intersection without
   GSE108940 is also empty; the core four sources intersect in **3,966** genes
   (median fused rank 4,178). The three intersection rows carry no Jaccard (the
   operands are identical) and have a `note`; `n_a`/`n_b`/`n_union`/`jaccard` are
   blank. Rank positions are true fused ranks, e.g. tenv-HVG(1000) ∩
   GSE243674_DEG spans fused ranks 2–17,156 (median 4,987).
4. **Inflection — explicit "no inflection"** — `output/gene-axis/inflection.json`
   (byte-identical to the pre-fix product): `"inflection": false`,
   `"n": null`. Kneedle online `S=1.0`: candidate `n = 1314`,
   `arc = 0.2533404287737608`; second-difference of the same normalized curve
   (`min_n = 100`, `minimum_curvature = 1e-6`): candidate `n = 18156`,
   `curvature = 0.003506475921318475`. The two candidates disagree by far more
   than the declared `agreement_tolerance = 0.25`, so no `N` is accepted. All ten
   leave-one-source-out replays also report `"inflection": false`, `n: null`.
5. **Tenv non-HVG share vs N** — `output/gene-axis/tenv_non_hvg_curve.csv`
   (18,164 rows; `n_non_hvg_k{1000,2000,5000}` and
   `non_hvg_share_k{1000,2000,5000}`; the top-N are the true fused top-N). At
   `K=2000` the share is 0.782 at N=500, 0.802 at N=1,000, 0.807 at N=1,314,
   **0.823 at N=2,000**, 0.839 at N=3,000, 0.856 at N=5,000, 0.880 at N=10,000
   and 0.890 at N=18,164. The curve is not flat: the fused top is enriched for
   tenv HVGs relative to the N=18,164 background (K=1000 0.884→0.945, K=5000
   0.540→0.725), but 82.3% of the fused top-2,000 is still outside tenv's
   top-2,000.
6. **Top-N direction composition** — `output/gene-axis/direction_summary.csv`
   (10 rows = 5 sizes × 2 rules). Symmetric up share 0.5260 / 0.5250 / 0.5380 /
   0.5343 / 0.5202 at N = 500 / 1,000 / 2,000 / 3,000 / 5,000; the
   `max_abs_log2fc` ranking's up share is 0.8680 / 0.8840 / 0.8625 / 0.7900 /
   0.7082 (the two rules select different genes). Exact direction ties are
   0 / 0 / 0 / 8 / 16; `n_unsupported` is 0 for the symmetric rule and 0/0/0/1/2
   for the max-|log2FC| rule. Up/down symmetry: 8,889 up-significant vs 7,857
   down-significant genes; across all significant genes the **down** set's median
   fused rank is 8,031 vs 8,779 for up.

The three #175 questions are answered with the same numbers and pinned in
`manifest.json → answers`: (1) HF sources agree (Jaccard 0.4255) while the
chronic-remodelling AF/HHD pair is nearly disjoint (0.0043) and HF×remodelling
cross pairs average 0.0233; (2) the fused top-N up share is 52–54% (symmetric) and
71–88% (max-|log2FC| ranking), mildly up-biased at the top and not symmetric —
across all significant genes the down set ranks better; (3) the
`n_sources_deg_significant` count ranking (descending count, **ties by ascending
fused rank**) agrees with the fused ranking at Spearman 0.8095 and top-N overlap
0.616 / 0.555 / 0.4775 / 0.4477 / 0.7604 at N = 500 / 1,000 / 2,000 / 3,000 /
5,000. A frozen `gene_index` tie-break instead gives Spearman 0.6791 and top-500
overlap 0.204; the product's declared tie-break is the fused-rank one.

## Independent re-derivation (every published number)

Every number quoted in `record.md`, this file, `stage2/01_bases/README.md`,
`GLOSSARY.md` and `docs/research/stage2-stage3-gene-axis-rank-aggregation.md`
was recomputed from the tracked evidence by a snippet that imports **no build
code**. It reads the four tracked CSVs and the tracked `manifest.json` /
`inflection.json`; the per-gene direction and significance arrays (needed for the
Q2 medians, the direction composition and the DEG-set overlaps) are re-derived
from the regenerated untracked per-source `*_deg.csv` intermediates with an
independent wiring. Result: **85/85 checks PASS**. The snippet:

```python
import json, numpy as np, pandas as pd
OUT = 'stage2/01_bases/output/gene-axis'
f = pd.read_csv(f'{OUT}/fused_rank.csv'); s = pd.read_csv(f'{OUT}/source_summary.csv')
o = pd.read_csv(f'{OUT}/overlap.csv'); d = pd.read_csv(f'{OUT}/direction_summary.csv')
c = pd.read_csv(f'{OUT}/tenv_non_hvg_curve.csv')
m = json.load(open(f'{OUT}/manifest.json')); inf = json.load(open(f'{OUT}/inflection.json'))
N = len(f); g = f.gene_index.to_numpy()
fr = np.empty(N, int); fr[g] = f.fused_rank.to_numpy()
tenv = np.empty(N, int); tenv[g] = f.tenv_hvg_rank.to_numpy()
cnt = np.empty(N, int); cnt[g] = f.n_sources_deg_significant.to_numpy()
order = np.argsort(fr); checks = []
def eq(label, a, b, tol=1e-12):
    assert abs(float(a) - float(b)) <= tol * max(1.0, abs(float(b))), (label, a, b)
    checks.append(label)
for k in (1000, 2000, 5000):                       # corrected non-HVG curve
    for n in (500, 1000, 1314, 2000, 3000, 5000, 10000, 18164):
        eq(f'curve k{k} n{n}', np.mean(tenv[order[:n]] > k),
           c.loc[c.n == n, f'non_hvg_share_k{k}'].iloc[0])
co = np.lexsort((fr, -cnt)); cr = np.empty(N, int); cr[co] = np.arange(1, N + 1)
q3 = m['answers']['q3_effect_replication_vs_deg_count']
eq('q3 rho', pd.Series(fr).rank().corr(pd.Series(cr).rank()),
   q3['spearman_fused_vs_count_rank'])
for n in (500, 1000, 2000, 3000, 5000):            # Q3 fused-rank tie-break
    eq(f'q3 top{n}', len(set(order[:n]) & set(co[:n])) / n,
       q3['top_n_agreement'][str(n)]['share'])
ids = ['gse141910', 'gse116250', 'gse243674', 'gse128188', 'gse282618',
       'gse160516', 'gse108940', 'dtoxs_ipsc', 'gse244740']
eff = np.zeros((N, len(ids))); l2 = np.zeros((N, len(ids)))
sig = np.zeros((N, len(ids)), bool)
for i, name in enumerate(ids):
    t = pd.read_csv(f'stage2/01_bases/execution/gene_axis/sources/{name}_deg.csv')
    ix = t.gene_index.to_numpy(); eff[ix, i] = t.effect; l2[ix, i] = t.log2fc
    sig[ix, i] = t.padj < 0.05
be = eff.max(1); b = eff.argmax(1)
dir_ = np.where(be > 0, np.sign(l2[np.arange(N), b]).astype(int), 0)
mab = np.nanmax(np.abs(l2), 1); mab = np.where(np.isfinite(mab), mab, 0.0)
any_sig = sig.any(1); up = any_sig & (dir_ > 0); down = any_sig & (dir_ < 0)
q2 = m['answers']['q2_up_down_rank_symmetry']
eq('q2 up', up.sum(), q2['n_up_significant']); eq('q2 down', down.sum(), q2['n_down_significant'])
eq('q2 med up', np.median(fr[up]), q2['median_fused_rank_up'])
eq('q2 med down', np.median(fr[down]), q2['median_fused_rank_down'])
eo = np.lexsort((fr, -mab))
ties = (np.isclose(eff, be[:, None], rtol=0, atol=1e-12)).sum(1) - 1
for rule, seq in (('symmetric', order), ('max_abs_log2fc', eo)):   # both direction rules
    for n in (500, 1000, 2000, 3000, 5000):
        row = d[(d.n == n) & (d.rule == rule)].iloc[0]; ch = seq[:n]
        eq(f'{rule} n{n} up', (dir_[ch] > 0).sum(), row.n_up)
        eq(f'{rule} n{n} ties', ties[ch].sum(), row.n_direction_ties)
disease = {name: set(np.flatnonzero(sig[:, i])) for i, name in enumerate(ids)}
for name in ids:                                   # overlap sizes + true fused ranks
    inter = set(np.flatnonzero(tenv <= 2000)) & disease[name]
    row = o[(o.comparison == 'tenv_hvg_vs_source_deg') & (o.set_a == 'tenv_hvg_2000')
            & (o.set_b == f'{name}_deg')].iloc[0]
    eq(f'overlap {name} size', len(inter), row.n_intersection)
    if inter:
        eq(f'overlap {name} median', float(np.median(fr[sorted(inter)])), row.fused_rank_median)
ordf = f.sort_values('fused_rank'); x = ordf.fused_rank.to_numpy(float); y = ordf.rank_sum.to_numpy()
xn = (x - x[0]) / (x[-1] - x[0]); yn = (y - y.min()) / (y.max() - y.min()); arc = yn - xn
eq('kneedle n', x[arc.argmax()], inf['kneedle']['n'])
curv = yn[2:] - 2 * yn[1:-1] + yn[:-2]
curv = np.where(np.arange(2, len(yn)) >= 100, curv, np.inf)
eq('second-difference n', np.argmin(curv) + 2, inf['second_difference']['n'])
for _, row in s.iterrows():                        # hvg share/density arithmetic
    if pd.isna(row.n_hvg_total):
        continue
    eq(f"{row.source} hvg_share", row.n_hvg_on_axis / row.n_hvg_total, row.hvg_share)
    eq(f"{row.source} hvg_density", row.n_hvg_on_axis / row.n_axis_expression, row.hvg_density)
print(f'{len(checks)} independent checks PASS')
```

## Tests (fixture regressions that would have caught each bug)

`execution/tests/test_gene_axis.py` — **31 passed**. `IndexSpaceTests` pins the
contract with fixtures where fused-rank order differs from `gene_index` order:

- `test_align_to_gene_index_scatters_sorted_rows` — scattering a fused-sorted
  column back to `gene_index` order.
- `test_tenv_curve_uses_fused_order_not_gene_index` — the curve fixture gives a
  different, asserted answer for fused order vs `gene_index` order.
- `test_overlap_rank_stats_use_gene_indexed_fused_rank` — rank stats are looked up
  by `gene_index`, not by row position.
- `test_direction_composition_rules_select_different_genes` and
  `test_direction_composition_uses_fused_rank_values_not_row_order` — the two
  rules differ on a fixture and the symmetric rule follows fused-rank values.
- `test_deg_count_tie_break_is_fused_rank` — the fused-rank tie-break and the
  `gene_index` tie-break give different overlap/spearman.
- `test_transcriptome_hvg_sets_project_top_k_onto_axis` and
  `test_source_summary_hvg_share_and_density` — the new HVG-set/share definition.

## Commands and gate results (run 2026-11-13, alone, on this branch)

Pre-flight: `free -g` → 47 GiB total / 42 GiB available (≥ 12 GiB required); one
heavy job at a time; `OMP_NUM_THREADS=2`, `ulimit -v 8388608` for the build.

```bash
script/check_paths.sh
# PASS: 1,034 tracked text files scanned; path literals 6,521 occurrences /
# 289 distinct; 6,491 / 269 resolved; 30 / 20 whitelisted (20 whitelist
# entries); 0 unresolved. stage2 retired-step audit PASS, stage3 source-first
# audit PASS.

script/check_tracked.sh
# PASS: 322 code files scanned, 351 tracked, 13 documented expected-ignored
# (stage1/xpert-fda-miri-v0), 0 holes.

MPLCONFIGDIR=/tmp/mpl python3 -m pytest -q \
  stage2/01_bases/execution/tests/test_gene_axis.py
# 31 passed in ~1.3 s (one harmless RequestsDependencyWarning on import).

ulimit -v 8388608 && OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 /usr/bin/time -v \
  python3 stage2/01_bases/execution/build_gene_axis.py
# exit 0; wall 0:11.0-0:11.1; Maximum resident set size 507,640-528,712 KiB.

ulimit -v 8388608 && OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  /usr/bin/time -v python3 stage2/01_bases/execution/build_gene_axis.py --verify
# {"result": "PASS", "inputs": 20, "outputs": 7, "n_axis": 18164,
#  "inflection": false, "inflection_n": null, "schema": "stage2-gene-axis-evidence-v1",
#  "manifest_sha256": "cf9dc2dec7cd5bb72a7b6d503a2bb246675e540403c8af9ba4e012b0fb56f89b"}
# exit 0; wall 0:09.8-0:10.9; Maximum resident set size 496,988-505,896 KiB.

# Figure CSV-only replay in a directory containing only the copied CSV and renderer:
python3 render_gene_axis_source_summary.py --csv plotted.csv --png out.png
# rendered out.png from plotted.csv; sha256
#   916738888f1f25e3954e48489661edb2b6bd7a018f19f1ea6391debb4e37c255
#   == the committed figure/gene_axis_source_summary.png.

/usr/bin/time -v ./script/verify_all.sh
# RESULT: PASS — verified 39/39 item(s); expected items: 39 (includes all eight
# disease_data items). exit 0; wall 4:23.61; Maximum resident set size 79,500 KiB.
# (Pre-flight `free -g`: 47 GiB total / 42 GiB available.)

# No accepted product / corpus / truth / Stage 1 byte changed (712b144 base):
git diff 712b144 -- stage2/01_bases/output/gene-identity \
  stage2/01_bases/output/tenv-hvg 'stage2/01_bases/output/tenv-truth*' \
  stage2/01_bases/output/training-corpora stage1
# (no output; exit 0)
```

## Tracked digests (sha256)

| file | sha256 |
|---|---|
| `stage2/01_bases/output/gene-axis/source_summary.csv` | `673f58d66a6b35ce93d721e654175e48ddf1d175bb92f31bbf60cf49a9bc9dc5` |
| `stage2/01_bases/output/gene-axis/fused_rank.csv` | `c1c632f3f166b93f206d441cbd96a7e031fcbd226d3779dc5a610a4c90be1a47` |
| `stage2/01_bases/output/gene-axis/overlap.csv` | `45b8eafbd27c69e39d5a385354c375aa569c0d6f1ba8174d58fab59fba1a74f6` |
| `stage2/01_bases/output/gene-axis/inflection.json` | `0f2c955f89b489aaf781e5724cdf2ed44487f8051b0ebcc8f15fe8a1c32201f4` |
| `stage2/01_bases/output/gene-axis/tenv_non_hvg_curve.csv` | `b421a4b64e7648cda48d5766f014db127ae419b0f2503b99855906899e7981ae` |
| `stage2/01_bases/output/gene-axis/direction_summary.csv` | `23a6b8c6bb3c5f0534f44f3b77d4eed2a95f669a264bd5ea0282f8518eb108f5` |
| `stage2/01_bases/output/gene-axis/mouse_human_orthologs.csv` | `cc22358e0ef0d886bb304d940fd8e6dac70c7c5008811f9dfb1f524fc31cbb6a` |
| `stage2/01_bases/output/gene-axis/manifest.json` | `cf9dc2dec7cd5bb72a7b6d503a2bb246675e540403c8af9ba4e012b0fb56f89b` |
| `stage2/01_bases/figure/data/gene_axis_source_summary.csv` | `7361eae4170b8601c5144e214446c15277ee087fe0221396a2d0daa5f1e05160` |
| `stage2/01_bases/figure/gene_axis_source_summary.png` | `916738888f1f25e3954e48489661edb2b6bd7a018f19f1ea6391debb4e37c255` |

`--verify` additionally checks all 20 pinned input sha256 values
(`manifest.json → inputs`); none of them changed in this round (the two DEG-CSV
JSON summaries are not build inputs, and the accepted ortholog table is
byte-identical).

## "No accepted product / corpus / truth / Stage 1 byte changed" proof

The command in the gate block above is empty (exit 0): `output/gene-identity/`,
`output/tenv-hvg/`, `output/tenv-truth/`, `output/tenv-truth-successor/`,
`output/training-corpora/` and `stage1/` are byte-identical to `712b144`.
`output/shared-expression/`, `pert_data/` and `disease_data/*/source/` are
untouched as well (no build wrote to them; every `script/verify.sh` passes over
the unchanged manifests). The evidence product changes are confined to
`output/gene-axis/` (regenerated CSVs/JSON + manifest), the figure CSV/PNG, the
tests, the build script and the R generator; the two `disease_data/*/local/deg/
*_deg.json` files only gain the already-true DEG-CSV digest, and the CSVs are
unchanged.

## Deviations and reading notes

- **Index-space contract.** `fused_rank.csv` is sorted by `fused_rank`; all
  consumers now scatter it back with `align_to_gene_index` before indexing by
  `gene_index`, and the fixture tests pin that. `fused_order(fused_rank)` on the
  sorted frame is the identity permutation and must never be treated as a gene
  ordering — that was the bug class.
- **`overlap.csv` has 81 data rows, not 71.** 27 tenv-HVG × DEG at
  K = 1000/2000/5000 + 15 tenv-HVG × source-HVG at K = 1000/2000/5000 + 36 DEG
  pairs + 3 intersections. The multi-way rows have blank
  `n_a`/`n_b`/`n_union`/`jaccard` and a `note`.
- **`hvg_share` is per source, not per axis.** A source's HVG ranking is over its
  whole measured transcriptome; `hvg_share` is the share of its declared top-2,000
  set that lands on the axis (`hvg_density` is the share of its on-axis measured
  genes in that set). tenv's share is 1.0 by construction.
- **Direction is not significance-filtered.** The composition rule takes the sign
  of the maximum-effect per-gene DEG row; `n_unsupported` counts genes with no
  reported effect. The #175 up/down symmetry counts separately require
  `padj < 0.05` in at least one source.
- **#175 question 3 tie-break.** The DEG-count ranking breaks ties by **ascending
  fused rank** (`np.lexsort((fused_rank, -count))` on gene-indexed arrays). This
  is declared in `manifest.json → constants.deg_count_tie_break` and in the
  record; a frozen `gene_index` tie-break gives Spearman 0.6791 and top-500
  overlap 0.204 instead of the tracked 0.8095 / 0.616.
- **R ortholog candidate policy.** `build_mouse_orthologs.R` writes
  `execution/gene_axis/mouse_human_orthologs.candidate.csv` (candidates stay
  outside formal `output/`, ADR 0062 / stage2 AGENTS "build candidates outside
  formal output, verify, then replace"), compares the candidate sha256 with the
  accepted table and replaces only when byte-identical or when `--args --accept`
  is given. This round's run was byte-identical; the accepted table and its
  digest are unchanged.
- **DEG JSON digests.** `GSE141910_deg.json` and `GSE116250_deg.json` now record
  `digests.local/deg/<GSE>_deg.csv` (sha256
  `3beed5df35bde8e748d03a0c76e2001ffa8f9dff963d60f5472363d88869e2da` and
  `23d89527ca1f220cad172d2bbf47c74c4509c181b8c0d6515039be3e1c04bb6d`), matching
  the unchanged CSVs and #177 acceptance 2. The three other items already
  carried the digest in the same shape (GSE243674 under
  `outputs.deg_csv_sha256`).
- **Pre-existing `GSE108940/source/` anomaly, preserved not introduced.** That old
  item's `source/` holds the locally generated limma DEG table instead of the raw
  download; ADR 0030 purity was already violated before this round. The
  provenance correction keeps `source/` byte-identical (its manifest still
  passes) and records the anomaly in
  `docs/research/stage2-stage3-gene-axis-gse108940-provenance.md`.

## Open decisions left to #173 / #174

- `N` (axis size) — there is no inflection; choose from
  `tenv_non_hvg_curve.csv` plus the Stage 3 sample arithmetic.
- Per-source quota and the tenv non-HVG cap — `source_summary.csv` coverage and
  `tenv_non_hvg_curve.csv` give the trade-off; the corrected curve is increasing
  in N, so the fused top is enriched for tenv HVGs relative to the background,
  but a capped non-HVG share still rejects most of the disease evidence.
- Direction rule — symmetric vs `max_abs_log2fc` composition in
  `direction_summary.csv`; the symmetric ranking is mildly up-biased (52–54% up)
  while the effect-magnitude ranking is strongly up-biased at the top.
- Animal weighting and the >5% category rule — `rank_sum_animal_half` is a
  reported sensitivity, not a choice; the AF/HHD overlap evidence informs
  whether chronic-remodelling categories belong.

#174 owns any release (`genes.csv` or otherwise); the old 18,164-gene axis stays
read-only and no new corpus, model or truth is produced.
