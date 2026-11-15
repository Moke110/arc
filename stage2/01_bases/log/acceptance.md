# Acceptance record — Stage 2 `01_bases` rebuild (#192, ADR 0066 / ADR 0068 / ADR 0069 / ADR 0062)

> **v2 republish (#208, 2026-11-15).** The product digests in §1, §5b and
> §6 are the #192/#194 record and are superseded by the ADR 0071 v2
> republish: see **§10** for the six-category disease-DEG, six-disease
> fusion, `f2000` axis and re-pinned views, and `inputs.json` for the
> machine-readable old-to-new digest mapping.
>
> **Figure re-pin (#210, 2026-11-15).** The two gene-axis PNG digests and the
> panel labels in §4 and in the 2026-11-14 figure paragraphs are superseded:
> both PNGs are re-rendered and re-pinned from their committed v2 plotted
> values, the coverage panel-b and rank-sum x labels are corrected, and both
> views are re-published on the new gene-axis manifest — see **§10.10** and
> `inputs.json.figure_repin`.

Accepted 2026-11-13 on branch `ticket/192-acceptance` (worktree `.scratch/worktrees/t192`,
base `dev=cb72b92`). Re-accepted after the post-review fixes of #194 on branch
`fix/194-products` (worktree `.scratch/worktrees/fix194a`, base `dev=d5512dd`); the
product table below records that post-#194 replay, and §5 records the re-sealed
raw-count provenance. Environment: `/home/chang/miniconda3/envs/arc/bin/python`
(3.13.13; anndata 0.12.16, h5py 3.16.0, numpy 2.4.6, pandas 2.3.3), R worker
`/home/chang/miniconda3/envs/biomni_e1/bin/Rscript` for the cardiotox DESeq2/edgeR
products, `OMP_NUM_THREADS=2`, `ulimit -v 10485760` fallback, one heavy job at a
time. Every replay below rebuilt from the pinned inputs into
`execution/checkpoints/` and compared the manifest and every output digest; the
accepted release was never written. `free -g` pre-flight on this host: 42 GiB
available, 12 GiB floor respected.

## 1. Product acceptance (replay gate)

| product | driver / mode | manifest sha256 | result | wall | peak RSS |
|---|---|---|---:|---:|---:|
| `output/raw_counts/` | `build_raw_counts.py verify` (`source-accepted-replay`) | `e02e5a1b…f6f406` | PASS | 15.9 s | 367 MiB |
| `output/corpus/` | `build_corpus.py verify` | `f8c35ccc…9d5d9f6d` | PASS | 1 m 10.5 s | 462 MiB |
| `output/disease-deg/{icm,nhf,genetic_cardiomyopathy}/` | `build_disease_deg.py verify --all` | ICM `1f453553…`, nHF `afb50564…`, genetic `304d8182…` | PASS (3/3) | 6.9 s | 282 MiB |
| `output/disease-deg/{cardiotox,cardiotox_TOP2,cardiotox_other}/` | `build_disease_deg_cardiotox.py verify` (R: DESeq2/edgeR) | `06e4e08b…`, `5ef7fa40…`, `12a3e1ba…` | PASS (3/3) | 1 m 40.1 s | 1.47 GiB |
| `output/gene-axis/` | `build_gene_axis.py --verify` | `641057dd…aa3ac2` | PASS | 0.43 s | 140 MiB |
| `output/view-2k/` | `build_view.py verify` (payload published first) | `72be67e1…f73d72` | PASS | 3 m 34.6 s (both views) | 449 MiB |
| `output/view-5k/` | `build_view.py verify` | `70298831…2d5c60` | PASS | (both views) | (both views) |

The two view payloads are gitignored and were absent from the fresh worktree, so
`build_view.py publish` materialised them first (3 m 32.2 s, 427 MiB); the
published manifests were **byte-identical to the tracked ones** (`git diff` empty),
and the subsequent `verify` replayed both views from the pinned corpus + axis.
`output/tenv-hvg/` is a retained pre-0066 product, not rebuilt by this ticket
(manifest `953103a1…`; see README *Retained tenv-HVG ranking*).

**2026-11-14 gene-axis update.** The release gained
`output/gene-axis/source_topn_coverage.csv` (each source's top-`N` coverage of
both axes, percent columns) and the two-panel evidence figure (gray candidate
histogram plus a colored selected-gene overlay, replacing the
interval/median/inflection panel), and gained the human-facing
`gene_axis_fill_rule.docx` (a reading aid, **not** a manifest output). The replay
gate was re-run: `build_gene_axis.py --verify` PASS, 32 inputs / **6** outputs,
manifest sha256 `df1179fd…4350fa`, 0.45 s / 142 MiB; the recorded
`gene_axis_verify.json` is that run. The CSV-only renderer reproduces the new PNG
byte-identically (`5c54608f…`), and `execution/make_gene_axis_docx.py`
reproduces the DOCX byte-identically (`37882ca4…`) from the coverage CSV + PNG
(both are test cases in `test_gene_axis`). Because both views pin the gene-axis
release manifest as an input, they were re-published and re-verified against the
new manifest: the view **payload digests are unchanged** (byte-identical h5ad and
`library_totals.csv`; only the pinned input differs), the view manifests are now
`1ea15c9c…` / `1747adfa…`, and both replays passed (`publish` 3 m 23.2 s /
449 MiB, `verify` 3 m 32.5 s / 430 MiB; `log/view_publish_verify.json` and
`log/view_verify.json` are that run).

**2026-11-14 AF removal (ADR 0069) and final re-seal.** The `AF` category was
removed from 01: `output/disease-deg/af/` is deleted, the DEG roster is
`nHF`/`ICM`/`genetic_cardiomyopathy`, and the gene-axis roster is five categories
plus tenv-HVG (candidate table 12 columns, `axis_composition.csv` 40 rows,
`source_topn_coverage.csv` 7 columns, `fusion.excluded_categories: []`). The axes
and rank sums are unchanged (`G_2000` `cb07535f…`, `G_5000` `15e5c49e…`,
intersection 1,972, `rank_sum` 0.0237–5.8103). Final replays: gene axis
`--verify` PASS, **30 inputs / 6 outputs**, manifest `641057dd…aa3ac2`, 0.43 s /
140 MiB; disease DEG `verify --all` PASS (3/3), 6.9 s / 282 MiB; views
`publish` 3 m 28.8 s / 429 MiB then `verify` 3 m 34.6 s / 449 MiB PASS with
**unchanged payload digests** and manifests `72be67e1…` / `70298831…`. The rule
document (`6bb82423…`) and the PNG (`5c54608f…`) replay byte-identically. The
manifest's `fusion.rank_sum` wording is now
`"smaller is better (a descending value is better); ties on gene_index"` — the
same rule as before, but no longer readable as "higher is better".

Record files: `raw_counts_verify.json`, `corpus_verify.json`,
`disease_deg_verify.json`, `cardiotox_verify.json`, `gene_axis_verify.json`,
`view_publish_verify.json`, `view_verify.json`; the machine-readable digest
inventory is `inputs.json`.

**2026-11-14 source-coverage figure.** The release gained
`figure/data/gene_axis_source_topn_coverage.csv` and
`figure/gene_axis_source_topn_coverage.png`, a two-panel grouped bar figure: one
x group per fused source (`tenv-HVG` + the five categories; `af` is not a
source, ADR 0069) and one bar per released axis (`G_2000` blue, `G_5000`
orange). Panel a ("Covered Gene Count") draws the absolute covered gene count
and panel b ("Top-N Coverage Rate") the same cells as percentages, both from
`bases.gene_axis.source_topn_covered_counts` / `source_topn_coverage`; the
integer count is the new public companion of the published percentage. The
plotted values are a manifest output (14 rows). The replay gate was re-run:
`build_gene_axis.py --verify` PASS, **30 inputs / 7 outputs**, manifest
`cfd3ba11…1a6317`, 0.48 s / 139 MiB; the CSV-only renderer reproduces the PNG
byte-identically (`51ea4bc8…`). Both views pin that manifest, so they were
re-published and re-verified: payload digests **unchanged**, view manifests now
`44f0c0df…` / `d7ac0227…`, both replays PASS (`publish` 3 m 27.5 s / 430 MiB,
`verify` 3 m 30.5 s / 435 MiB; `log/view_publish_verify.json` and
`log/view_verify.json` are those runs). `test_gene_axis` is now 49 tests.

### Structural and semantic checks exercised by the replays

- `raw_counts`: integer-only dense tenv (4,510 × 19,613) and CSR Tahoe
  (67,018 × 18,164); no normalized layers; probe bridge `exact_hgnc` over the
  frozen axis in order; token mapping covers 18,164 frozen genes.
- `corpus`: tenv 4,478 rows = 4,206 treated + 272 vehicle; `tahoe_bulk` 66,718 =
  65,318 treated + 1,400 `DMSO_TF`; `role` booleans; matched-control groups
  (16 tenv groups / 700 Tahoe groups; 0 treated without control; 332 exclusions
  with reasons); no normalization/Δ columns.
- `disease-deg`: per-dataset tables/ranks rebuilt; rank rule (descending
  `|log2FC|·−log10(padj)`, max collapse, `NA` padj non-significant, unmeasured =
  `n_measured + 1`, `rank_norm = rank/n_measured` so an unmeasured gene is
  `(n_measured+1)/n_measured`, `gene_index` tie-break) re-verified against the
  pinned item DEGs; catalogs' self-consistent `catalog_sha256`; every record
  carries explicit `informative`/`informativeness` evidence (ADR 0068) and only
  informative datasets enter level 1; level-1 inputs are the released rank frames.
- `gene-axis`: level-1 dataset fusion inside each category, level-2 six-source
  `rank_sum` (five categories + tenv-HVG, `fusion.excluded_categories: []`),
  `pN ∈ 0..6`; axes exactly 2,000/5,000 ascending **original** `gene_index`
  (never renumbered), fill `0.5·N` tenv-HVG + `0.5·N` pN/rank_sum remainder;
  non-nested axes with reported intersection **1,972**; the source top-`N`
  coverage table (two percent rows, 7 columns); the two-panel figure values
  (`sum(selected) = N` per panel, `sum(count) = 18,164`); manifest digests and
  the 18,164-row candidate table.
- `view-*`: normalization on the **full-axis** library total then subset (never
  the view-local denominator); `control` = mean of the named matched-control
  rows; `Δ = X − control` exactly; `X/control/Δ` float32; obs carries row
  identity + `library_total`; no `{N}k_eval_truth/` anywhere.

## 2. Downstream reader compatibility

Every retained reader module was exercised through its own test file under the
`arc` env (all green): `test_raw_counts` (8), `test_corpus` (10), `test_view`
(17), `test_gene_axis` (45), `test_deg_rank` (15), `test_disease_deg_human_sc`
(18), `test_disease_deg_isogenic` (35), `test_drug_classes` (17),
`test_disease_deg_cardiotox` (11), `test_identity`, `test_fold_seed`,
`test_palette`, `test_hvg` (2, after removing the retired metrics-QA assertion), `test_native`,
`test_tahoe`, `test_raw_provenance`.
`02_models_tenv` is deliberately **not** re-pointed (ADR 0066 *Deferred*); its
target contract (`view-{N}k` `X`/`control`/`delta`, treated-row selection) is
documented in the README, and its stale imports of the deleted evaluation
modules are the deferred work item, not a live consumer of any accepted product.

## 3. Cross-category catalog invariant (new #192 gate)

`execution/tests/test_acceptance_invariants.py` (9 tests) passes over the six
released catalogs. It pins the canonical `stage2-disease-deg-catalog-v1` shape
(`datasets` list of full records, per-record `method`, self-consistent 64-hex
`catalog_sha256`, no legacy `test_method`/`dataset_records`), positive
`n_case`/`n_control` **or** an explicit `*_basis`, 64-hex `derived_inputs` where
declared, an explicit boolean `informative` plus its `informativeness` reason
(the degenerate-source field of ADR 0068), and — explicitly, not silently — the documented cross-ticket variance:
#187/#188 declare the published `gene-identity/genes.csv` axis, a
`rank_contract`/`rejected` block, bare-digest `inputs` and gitignored
`derived_inputs`; #189 declares the `tables/genes.csv` axis, `{bytes,sha256}`
`inputs` and no `derived_inputs`. A drift to a third shape fails the gate.

## 4. CSV-only figure reproduction (products unreadable)

With every `output/` product directory made unreadable (`chmod a-r`; only
`output/api` readable), the renderer
`figure/script/render_gene_axis_rank_sum_histogram.py --csv
figure/data/gene_axis_rank_sum_histogram.csv` re-rendered the committed PNG
byte-identically: `3184b37f…3fcbe134` = `3184b37f…3fcbe134` (exit 0; pre-redesign
digest). The 2026-11-14 two-panel figure still replays from its CSV alone with
`output/` unreadable (test case; digest `5c54608f…`). Separately, the rule
document replays byte-identically from the coverage CSV + PNG —
`execution/make_gene_axis_docx.py` → `gene_axis_fill_rule.docx` (`6bb82423…`) —
but it is a document, not the evidence figure, and it does read the published
coverage CSV. The retired diagnostic PR figure and its renderer are deleted.

## 5. "No retired path read" evidence

Input sets were read from the accepted manifests (`inputs.json`), not assumed.
`output/raw_counts/` is built only from `gene-identity/{genes,crosswalk}`,
`pert_data/GSE262419` + `pert_data/Tahoe-100M` source manifests and the source
`gene_metadata.parquet`, and the tracked
`log/raw_counts_provenance.json` (the two pre-0066 source-replay acceptance
records plus the native tenv build manifest, restored verbatim from git history
and digest-pinned). `output/corpus/` reads only
`raw_counts/` payloads + `gene-identity` + Tahoe drug metadata.
`disease-deg` reads item DEG JSONs/CSVs under `disease_data/`/`pert_data/` and
`gene_identity`. `gene-axis` reads `disease-deg/`, `tenv-hvg/descriptive.csv` and
`gene-identity`; `view-*` reads `corpus/` + `gene-axis/` only. No retained
builder or reader opens `shared-expression/`, `training-corpora/`, the old
raw-count directories, `tenv-truth*` or `historical-training-release/`
(verified by reviewing every manifest `inputs` key and by
`git grep` over retained code after deletion).

The one accepted, one-way transition is documented instead of hidden:
`raw_counts`' default **rehome** mode pins the tracked
`log/raw_counts_provenance.json` record, whose `acceptance.tenv` /
`acceptance.tahoe` sections are the two pre-0066 source-replay acceptance
records and whose `payload_sha256` map pins the rehomed bytes (tenv
`efd07632…`, Tahoe `a830a34a…`; the bytes in `raw_counts/` are byte-identical,
`source-accepted-replay`). Every other pinned input resolves in a merged
checkout: the two source item manifests, the source `gene_metadata.parquet`,
the frozen axis/crosswalk and the provenance record. The retired
`log/records/` and `output/{tenv,tahoe}-raw-counts/` paths are **not** pinned by
any live manifest; their digests remain in `inputs.json.retired_evidence_digests`
as history. A from-scratch replay still uses the documented stronger mode
`build_raw_counts.py --full-rebuild publish|verify` (Tahoe build 1 h 35 m,
independent source replay 2 h 35 m). The retained `tenv-hvg` ranking has the
analogous deferred re-point documented in the README.

## 5b. Post-review re-acceptance (#194)

The #194 fixes changed three products: every disease-DEG rank table (`rank_norm`
for unmeasured genes), the gene-axis fusion (`af` excluded, `pN ∈ 0..6`, both axes
and composition/histogram rebuilt; superseded by the 2026-11-14 AF removal, ADR
0069) and both views (rebuilt from the new axes); the
raw-count release was re-sealed against the tracked provenance record without
changing a payload byte. All six replay gates were re-run green at the times and
peaks in §1, and the per-slice `*_verify.json` records were regenerated. The full
`stage2/01_bases/execution/tests` suite passed (17 files, each run directly from
`stage2/01_bases`; the two files that import `disease_deg_fixtures` resolve it from
`execution/tests/`).

`script/check_tracked.sh` is green. `script/check_paths.sh` is **red only in a
payload-free merged checkout**, where it reports the same 68 pre-existing,
payload-dependent unresolved literals on pristine `dev=d5512dd` and on
`fix/194-products` (identical distinct-literal sets: the branch adds none); with
the external `source/`/`local/` trees materialised it reports
`RESULT: PASS — every extracted path literal resolves`, the stage2 retired-path
audit `PASS` and the stage3 source-first audit `PASS`.

## 6. Old → new mapping (deletions)

| removed (commit group) | accepted substitute / rationale |
|---|---|
| `output/shared-expression/` (products) | `output/corpus/` integer counts + view normalization |
| `output/training-corpora/` (products) | `output/corpus/` + `output/view-2k/`/`view-5k/` |
| `output/tenv-raw-counts/`, `output/tahoe-raw-counts/` (products) | `output/raw_counts/` — same payload bytes, one manifest |
| `output/tenv-truth/`, `output/tenv-truth-successor/`, `output/historical-training-release/` (products) | none; evaluation redesign deferred; each view will own `{N}k_eval_truth/` |
| old `output/gene-axis/` file set: `source_summary.csv`, `fused_rank.csv`, `overlap.csv`, `inflection.json`, `tenv_non_hvg_curve.csv`, `direction_summary.csv`, `mouse_human_orthologs.csv` (product) | `candidate_rank_sum.csv`, `gene_axis_{2000,5000}.csv`, `axis_composition.csv` (retired fusion semantics not resurrected); added later: `source_topn_coverage.csv` and the human-facing `gene_axis_fill_rule.docx` (2026-11-14) |
| truth/evaluation chain: `build_significance.py`, `bases/{significance,truth_build,truth_successor,evaluation,metrics,plots}.py`, `deseq2_contrasts.R` + sealed copy, `evaluate_predictions.py`, `replay_pr.py`, `verify_hvg_metrics.py`, `export_diagnostic_predictions.py`, `publish_diagnostic_figure.py`, `hvg_metrics.md`, `truth_successor.md`, `figure/diagnostic_tenv_pr.*`, their tests/fixtures (code) | none; reworked with the deferred evaluation design |
| one-off migration/audit: `migrate_products.py`, `audit_metadata_bindings.py`, `verify_payload_digests.py` (code) | none; one-shot, records in git |
| superseded build chain: `build_shared.py`, `build_corpora.py`, `verify_training_source.py`, `bases/corpora.py`, `build_native.py`, `prepare_tahoe.py`, `gene_axis_sources.py`, `figure/gene_axis_source_summary.*`, `build_mouse_orthologs.R` (code) | `build_raw_counts.py` (`publish`/`--full-rebuild`), `build_corpus.py`, `bases.raw_counts`, `bases.corpus`; mouse ortholog table has no surviving eligible source |
| pre-0066 `log/` tree (`log/{155,156,157,161,162}-*`, `log/{records,health,diagnostic,historical,gene-axis}/`, old `log/README.md`) | this record + `inputs.json` + per-slice `*_verify.json` (ADR 0066: `log/` deleted then regenerated, never archived) |

Deletion commits are separate by class (products, truth/eval code, one-off
scripts + superseded chain, pre-0066 log) so the diff is reviewable. Git retains
all deleted code and compact files; the multi-GB ignored payloads are re-creatable
from the immutable `source/` trees (see `docs/data-inventory.md`).

## 7. Recorded errata (evidence; no ADR text edited)

### GSE270788 arm roster

ADR 0066's roster says "`GSE270788` NICM arm (7 donors vs 6 Donor)" and "All three
`GSE270788` arms use the same six Donor controls". The deposit contradicts both.
Measured directly from the item's `source/GSE270788_metadata.csv.gz`
grouped by the metadata `sample` column (46,384 cells): **Donor 7 / NICM 6 / ICM 5
/ AMI 5**. The series-matrix `!Sample_characteristics_ch1` disease state agrees.
Therefore the NICM arm is 6 case vs **7** shared Donor controls and the ICM/AMI
arms are 5 vs 7. The `orig.ident` column is a single value `HF` across all arms and
must not be used for grouping. No donor is relabelled: all 23 donors are used,
and the released `nHF` table already records 6/7 while `ICM`/`AMI` record 5/7
(`n_case`/`n_control` in the catalogs). ADR 0066's counts are superseded by the
deposit (design unchanged). **The ADR erratum is to be landed as a separate
docs-only ADR (0067) by the Lead**; this ticket deliberately edits no ADR text.

### View basis figures

ADR 0066 *Corpus and view semantics* quotes "12.9% / 36.3%" of full-axis CP10K
from an interim pre-#190 axis. On the **released** axes after the #194 fusion
amendment and both view rebuilds, measured from the view manifests (tenv, all
4,478 profiles): **G_2000 = 15.0548953%** and **G_5000 = 37.4273748%**
(treated-only 15.0892505% / 37.3846128%; `tahoe_bulk` 10.8317768% /
24.7366390%). The 2.5× baseline-change argument and the full-axis decision are
unaffected (measured denominator ratio ≈2.49×); the released values are
authoritative and recorded in both view manifests.

### Cosmetic, not fixed (out of scope)

`output/api/bases/disease_deg.py` emits a pre-existing `SyntaxWarning: invalid
escape sequence '\'` from a #187 docstring. It changes no behaviour and is owned
by a follow-up; recorded here so it is not mistaken for a new regression.

## 8. Decision closure for #173 (answers from the released products)

1. **Axis size — `G_2000` and `G_5000`, both released, because the size was never
   identifiable from the evidence and the model comparison is cheap.** The
   inflection diagnostic on the rank-sum curve is honest (`"inflection": false`
   with full Kneedle + second-difference parameters in the gene-axis manifest), so
   `N` is not data-selected; two sizes bracket the plausible range and an axis
   comparison can be run under one fixed hyperparameter protocol or a per-axis
   tuning budget. The released axes are 2,000/5,000 rows with 1,919 shared genes,
   ascending original `gene_index`, fill 0.5·N tenv-HVG + 0.5·N pN/rank_sum.
2. **Source quota — 50/50 tenv-HVG / candidate (`0.5·N` each).** The tenv-HVG half
   buys general chemical-perturbation variability (no direction); the candidate
   half buys disease-mechanism evidence ranked by `pN` descending then `rank_sum`
   ascending then `gene_index`. Both axes use the same quota, so the size
   comparison is not confounded by composition.
3. **Direction rule — selection on `|effect|`; no direction composition.** The
   per-dataset effect is `|log2FC|·−log10(padj)` with a per-gene maximum collapse
   across contrasts; the axis is a magnitude selection and `rank_sum` deliberately
   erases sign. Direction composition (the retired pre-0066 semantics) is not
   implemented anywhere.
4. **Animal and `>5%` inclusion — none; no animal source survives eligibility.**
   Whole-tissue bulk is ineligible even when it contains cardiomyocytes (cell
   composition/fibrosis dominates). The rejected tier is documented, not silently
   dropped: the six whole-tissue series (`GSE141910`, `GSE116250`, `GSE128188`,
   `GSE282618`, `GSE160516`, `GSE108940`) are marked retained-but-not-axis-sources
   in `disease_data/README.md`/`AGENTS.md`; `mouse_human_orthologs.csv` and
   `build_mouse_orthologs.R` were deleted with no surviving animal source. The
   `>5%`-prevalence idea is superseded by ADR 0066's category rules, not by a
   prevalence threshold.

## 9. Health (final tree)

Run from the worktree root after the deletions and after `git merge dev`
(`dev=cb72b92`; "Already up to date" — no drift).

- `script/check_paths.sh`: **PASS** — 987 tracked text files scanned, 6,425 path
  literals (0 unresolved; the two stage2 retired-root literals are whitelisted
  frozen evidence), stage2 retired-step-path audit **PASS**, stage3 source-first
  audit **PASS**, 10.0 s. The check resolves literals against the filesystem, so
  the worktree used read-only hardlinks to the gitignored external `source/`/
  `local/` trees (a fresh clone is not a merged checkout); nothing in the tracked
  tree was changed by that.
- `script/check_tracked.sh`: **PASS** — all scanned executable files are tracked
  (or in a documented expected-ignored tree).
- `script/verify_all.sh` (run alone, the documented ~9.7 GiB page-cache job):
  **PASS — 48/48 items**, 4 m 48.7 s, peak RSS 79,124 KiB. `free -g` available
  memory: 42 GiB before and 42 GiB after (the kernel reclaimed page cache; no
  pressure on the 12 GiB floor).

**2026-11-15 rule-document source table, excluded-AF block and code-review
follow-up (#198/#199).** `output/gene-axis/gene_axis_fill_rule.docx` replayed
byte-identically from its inputs — the coverage CSV, the figure PNG, the five
`output/disease-deg/<category>/catalog.json` files and the tracked
`execution/gene_axis_source_metadata.py` — with the documented no-arg command,
twice: 0.41 s wall, peak RSS 125,736 KiB, final sha256
`2bbe1f0e1fc8e328c90288f27c1e338d9c50e564eefe30c2b012fece0091ed88`. Earlier
docx digests (`37882ca4…`, `6bb82423…`, `044f7ece…`) are superseded and match no
committed file. 表 1 has 19 rows (18 released datasets + tenv-HVG; columns
`category / accession / sample type / groups / sample sizes / sequencing
technology / CM DEG count`), and the six cardiotox `groups` cells are now derived
from each catalog's own `shared_control_group` rather than a hard-coded token.
Its footnote defines the CM DEG count (`n_significant`, best contrast, BH
padj<0.05), the contrast semantics and the `GSE255612` 19/19-vs-18/16 mismatch,
and states the five non-obvious sample-size bases (GSE270788 46 GSMs / DToxS
1,171 GSMs / GSE243674 72 GSMs / GSE244740 3,072 wells / Reichart-CM 79 donors);
表 2 lists the five audited-and-excluded AF sources and the closing paragraph
declares `AF` a gap. The generator's five `--catalog` inputs and its new
`--metadata` module input are proven read by wiring tests (a changed catalog
`n_significant` and a changed `technology` each move the rendered cell). Test
inventory: 266 step tests / 55 `test_gene_axis.py` / 22
`test_gene_axis_source_metadata.py`; `test_gene_axis.py` alone 1.845 s, peak RSS
194,324 KiB; the full directory 8.384 s, peak RSS 323,476 KiB. The DOCX is not a
manifest output, so `build_gene_axis.py --verify` and `build_view.py verify` are
unaffected; the integrator re-ran both PASS in the main checkout, which holds the
gitignored payloads this worktree lacks. **Erratum:** the 2026-11-14 entries
describe the docx replay as reading only the coverage CSV + PNG; since
#197/#198 it also reads the five catalogs and the tracked metadata module, and
their docx digests (`37882ca4…`, `6bb82423…`) are historical. No release
product, manifest, axis or view changed. `docs/data-inventory.md`'s
`output/gene-axis/` row was refreshed (`du -sb` 1,461,704; docx 92,412 B;
`du -sb`, 2026-11-15).

## 10. Gene-axis v2 republish gate (#208, ADR 0071)

Accepted 2026-11-15 on branch `ticket/208-republish-acceptance` (base `dev=ed31945`).
Contract: ADR 0071 (DEG-only tiered sources, `DEG2000`/`f2000`) supersedes the
ADR 0066 measured-gene rank and `pN`; ADR 0068's non-informative-source rule is
kept; publication is ADR 0062 candidate-then-replace. The product-schema label in
`bases.deg_rank` moves `stage2-disease-deg-v1 -> stage2-disease-deg-v2`: the
published compact table now carries the authoritative pass-any-contrast `deg`
flag and `n_contrasts_deg`, the DEG-only `rank`/`rank_norm = rank/n_DEG`, and no
measured-gene denominator. The gene-axis release keeps its schema label and adds
`release.algorithm = f2000-deg2000-v2`; a pre-v2 manifest is refused by
`build_gene_axis.py --verify` with a self-describing error. Every driver was run
under `execution/capped.sh` with `/usr/bin/time -v`, one heavy job at a time;
the R products used `/home/chang/miniconda3/envs/biomni_e1/bin/Rscript`
(R 4.5.3, edgeR 4.8.2, DESeq2 1.50.2).

### 10.1 Replay gate

| product | driver / mode | result | publish wall / peak RSS | verify wall / peak RSS |
|---|---|---:|---:|---:|
| disease-deg registry (nHF, ICM, AF, genetic_cardiomyopathy) | `build_disease_deg.py publish --all / verify --all` | PASS | 0:13.24 / 345 MiB | 0:14.54 / 338 MiB |
| disease-deg cardiotox (cardiotox, cardiotox_TOP2, cardiotox_other) | `build_disease_deg_cardiotox.py publish / verify (R DESeq2/edgeR)` | PASS | 1:39.77 / 1501 MiB | 1:35.79 / 1502 MiB |
| disease-fusion (six diseases) | `build_disease_fusion.py publish / verify` | PASS | 0:01.13 / 138 MiB | 0:01.13 / 138 MiB |
| gene-axis (G_2000, G_5000, candidate, composition, coverage) | `build_gene_axis.py (publish) / --verify` | PASS | 0:00.35 / 129 MiB | 0:00.33 / 131 MiB |
| view-2k / view-5k | `build_view.py publish / verify` | PASS | 3:28.81 / 422 MiB | 3:28.69 / 432 MiB |

Six of the ten steps were `publish` (the accepted products were replaced) and
then replayed with `verify` byte-identically; every recorded
`*_verify.json` is the verify run. The cardiotox step samples the whole process
group (Python + R): peak group RSS 1,502,868 KiB during publish and 1,427,692 KiB
during verify; the views stream the corpus in bounded batches (peak group RSS
451,708 / 452,700 KiB). `raw_counts/` and `corpus/` were **not** rebuilt: their
manifests still hash `e02e5a1b…f6f406` and `f8c35ccc…9d5d9f6d`, the recorded
#192 digests, so no raw count or corpus byte changed.

### 10.2 Six-category disease-DEG products

| category | datasets | DEG union | DEG2000 | tiers entered |
|---|---:|---:|---:|---|
| `nhf` | 6 | 1661 | 1661 | [1, 2, 3] |
| `icm` | 3 | 1647 | 1647 | [1] |
| `af` | 4 | 1 | 1 | [1, 3] |
| `genetic_cardiomyopathy` | 8 | 1482 | 1482 | [1, 2] |
| `cardiotox_TOP2` | 3 | 4514 | 2000 | [2] |
| `cardiotox_other` | 3 | 297 | 297 | [2] |

The unions reproduce the ADR 0071 ledger exactly for `nHF` (1,056 -> 1,553 ->
1,661), `ICM` (1,647), `AF` (1, `SLC9A7`, tiers 1 -> 3) and
`genetic_cardiomyopathy` (59 -> 1,482). The two cardiotox classes deviate and the
deviation is the flag source, not a data change (§10.3): 4,363 -> **4,514** and
279 -> **297**, each reproduced exactly by re-deriving the retired
winning-contrast convention from the new tables. `DEG2000` is therefore 7,088
rows over 5,282 distinct genes (1,806 duplicate rows), and only
`cardiotox_TOP2` is cut (4,514 -> 2,000). All nine previously unpublished
sources are now published (`gse35108__dcm_ipsc_cm`, `gse165303__dcm_lv_bulk`,
`gse116250__hf_lv_bulk`, `gse141910__hf_lv_bulk`, `kuppe-mi__iz`, and the four
AF sources), and the cardiotox classes publish their own `tiering` record.

### 10.3 Recorded deviations from the ADR 0071 ledger

| product | ADR 0071 ledger | measured | delta | explanation |
|---|---:|---:|---:|---|
| cardiotox_TOP2 (DEG union) | 4363 | 4514 | +151 | the ADR ledger was measured on the pre-#208 published tables, which carried only the collapsed winning-contrast row per gene; the v2 product publishes the authoritative pass-any-contrast deg flag, and a gene whose non-winning contrast passed the threshold is now a DEG (flag_source published_deg_flag). Re-deriving the old convention from the new tables reproduces exactly 4,363, so the delta is the flag source, not a data change. |
| cardiotox_other (DEG union) | 279 | 297 | +18 | same flag-source correction as cardiotox_TOP2. |
| DEG2000 rows (six diseases) | 7070 | 7088 | +18 | the same two cardiotox flag corrections; every other category matches the ADR ledger exactly (nHF 1,661 / ICM 1,647 / AF 1 / genetic 1,482, cardiotox_TOP2 still cut to 2,000). |

### 10.4 Released axes `G_2000` / `G_5000` (f2000 fill)

| axis | tenv-HVG first half | second half from DEG2000 | tenv-HVG fallback | axis `f2000` distribution |
|---|---:|---:|---:|---|
| `G_2000` | 1,000 | 1,000 | 0 | 0:624, 1:258, 2:817, 3:248, 4:48, 5:5, 6:0 |
| `G_5000` | 2,500 | 2,500 | 0 | 0:1548, 1:2005, 2:1146, 3:248, 4:48, 5:5, 6:0 |

Both axes fill their second half entirely from `DEG2000` (4,906 / 4,330
candidates available for 1,000 / 2,500 slots), so the tenv-HVG fallback is
implemented, exercised by tests, and 0 in both axes. `f2000` is 0 for 12,882 of
18,164 candidates, 1 for 3,835, 2 for 1,146, 3 for 248, 4 for 48 and 5 for 5
(no candidate is a DEG in all six diseases). Because `f2000` and `rank_sum` do
not depend on `N`, the axes are nested: `G_2000` is a subset of `G_5000`
(|intersection| = 2,000, Jaccard 0.4); the retired v1 `pN` fill made them
non-nested. `AF` is present in the release (candidate column, catalog, coverage
row); its single `DEG2000` gene is in `G_5000` but not `G_2000` because the
1,000 `G_2000` disease slots are all `f2000 >= 2` genes.

### 10.5 DEG2000 coverage and composition

`axis_composition.csv` carries the declared halves, the fallback and the `f2000`
distribution with explicit denominators; `source_topn_coverage.csv` is long
format with numerator / denominator / rate per fusion group:

| axis | group | numerator | denominator | rate (%) |
|---|---|---:|---:|---:|
| `G_2000` | `tenv_hvg` | 1079 | 2000 | 53.95 |
| `G_2000` | `nhf` | 638 | 1661 | 38.41 |
| `G_2000` | `icm` | 550 | 1647 | 33.39 |
| `G_2000` | `af` | 0 | 1 | 0.00 |
| `G_2000` | `genetic_cardiomyopathy` | 635 | 1482 | 42.85 |
| `G_2000` | `cardiotox_TOP2` | 831 | 2000 | 41.55 |
| `G_2000` | `cardiotox_other` | 199 | 297 | 67.00 |
| `G_5000` | `tenv_hvg` | 3068 | 5000 | 61.36 |
| `G_5000` | `nhf` | 1151 | 1661 | 69.30 |
| `G_5000` | `icm` | 1042 | 1647 | 63.27 |
| `G_5000` | `af` | 1 | 1 | 100.00 |
| `G_5000` | `genetic_cardiomyopathy` | 1068 | 1482 | 72.06 |
| `G_5000` | `cardiotox_TOP2` | 1740 | 2000 | 87.00 |
| `G_5000` | `cardiotox_other` | 256 | 297 | 86.20 |

### 10.6 No retired path is read or published

* `candidate_rank_sum.csv` columns: symbol, gene_index, tenv_hvg_rank, nhf_rank, icm_rank, af_rank, genetic_cardiomyopathy_rank, cardiotox_TOP2_rank, cardiotox_other_rank, rank_sum, f2000 —
  no `p2000`/`p5000`, no `n_sources_measured`.
* the gene-axis manifest carries no `pN` field (checked over the serialised
  manifest) and the driver reads only `output/disease-fusion/` through
  `DiseaseFusion`; it never opens a within-dataset rank frame.
* every published `<dataset>.deg.csv` carries `deg` + `n_contrasts_deg`, and every
  `.rank.csv` ranks only the DEGs (`n_DEG + 1` for every other gene), checked for
  all 23 published datasets.
* the two cardiotox class catalogs now carry `n_deg`/`deg2000`/`tiering` from
  `bases.disease_deg.tier_ledger`; the registry catalogs carry the same record
  plus `rank_contract = stage2-disease-deg-v2`.
* views are re-pinned on the new axis manifest; the corpus and raw-count
  manifests are byte-unchanged.

### 10.7 Superseded digests and files

`inputs.json.superseded_digests` holds every old digest of a file whose bytes
changed (88 tracked paths: 82 product payloads, 2 figure plotted-value
CSVs and the four public API modules; plus 21 newly published product files), the
view payload digests and the figure/document entries. `log/README.md` lists the record files, including
the new `disease_fusion_verify.json`.

### 10.8 Handover

* **#209**: rebuild `gene_axis_fill_rule.docx` (v2 rule, per-disease table, six
  categories, AF included, new source rows) and extend
  `execution/gene_axis_source_metadata.py` with the five new source rows plus the
  AF sources; re-pin the six `DocxTests` and the two
  `test_gene_axis_source_metadata` coverage tests.
* **#210**: re-render and re-pin the source-coverage PNG from the regenerated
  DEG2000 plotted values. The rank-sum histogram PNG has **no owning ticket** and
  needs the same re-render/re-pin or retirement.
* **#211**: docs sweep (README/AGENTS/GLOSSARY still describe `pN`, five
  categories and the measured-gene rank), the ADR 0071 ledger correction
  (cardiotox 4,514 / 297), and the spec #200 / #207 "seven sources" wording.

### 10.9 Health

On the merged tree, `python -m unittest discover -s stage2/01_bases/execution/tests`
gives **398 tests: 5 failures + 5 errors**, all in the #209
docx/source-metadata and #210 figure byte re-pin suites (exactly the ten tests
listed in `inputs.json.health.red_tests`); every other suite is green, including
the two suites that were red on `dev` before this ticket
(`test_disease_deg_cardiotox`, `test_disease_deg_human_sc`,
`test_disease_deg_isogenic` — now v2 assertions). `script/check_paths.sh` PASS
(all three audits), `script/check_tracked.sh` PASS, `script/verify_all.sh` PASS
54/54. The tree is clean after the merge; no item `source/` or `local/` payload
was written by the republish.

### 10.10 Figure re-pin and label correction (#210, 2026-11-15)

Both committed gene-axis PNGs were stale against the v2 plotted values that
§10.1 regenerated. They are re-rendered from their own tracked CSVs with the
committed renderers and re-pinned; the CSV-only replay assertion now passes for
both, and no plotted value changed.

| figure | plotted values | PNG (old → new) | replay |
|---|---|---|---|
| `gene_axis_source_topn_coverage.png` | `figure/data/gene_axis_source_topn_coverage.csv` `103073c1…` → `00a04341…` (panel-b label only) | `51ea4bc8…` → `37c24d17…` | byte-identical (`SourceCoverageFigureTests`) |
| `gene_axis_rank_sum_histogram.png` | unchanged `c402554a…` | `5c54608f…` → `874f6020…` | byte-identical (`FigureTests`) |

**Labels corrected.** Two labels that the v1 wording made untrue:

- Panel b of the coverage figure was `Top-N Coverage Rate` (y = `coverage rate
  (%)`). The v1 denominator was `N` for every source; the v2 denominator is the
  group's own set — `|DEG2000_d|` for a disease (`af` 1, `cardiotox_other` 297,
  `genetic_cardiomyopathy` 1,482, `icm` 1,647, `nhf` 1,661, `cardiotox_TOP2`
  2,000) and `N` only for `tenv-HVG` (ADR 0071 decision 8). Under the old title
  `AF`'s 1-of-1 `G_5000` cell read as a 100% "Top-N" coverage. The panel is now
  `Coverage of Each Group` (y = `share of group covered (%)`). The label lives
  in `SOURCE_COVERAGE_FIGURE_PANELS`, so the plotted-values CSV and
  `manifest.coverage.figure.panels` carry it too: gene-axis manifest `078f604a…`
  → `c247dcc0…`, `--verify` PASS, 10 inputs / 7 outputs, 0.33 s / 131 MiB. Every
  axis and coverage payload digest is byte-unchanged (`source_topn_coverage.csv`
  stays `8db5fcbf…`).
- The rank-sum figure's x annotation was `rank_sum (descending is better)`,
  contradicting the released `RANK_SUM_RULE` ("Ascending is better": rank 1 is
  the strongest DEG and a gene outside `DEG2000_d` pays `(n_d + 1) / n_d`); the
  plotted selected genes sit at low `rank_sum`. It now reads `rank_sum (lower is
  better)`. Only the renderer constant changed.

Because both views pin the gene-axis release manifest, they were re-published
and re-verified against the new one: the view **payload digests are
byte-unchanged** (`tahoe_bulk.h5ad`, `tenv.h5ad`, `library_totals.csv`; only the
pinned input differs). Manifests: `view-2k` `30d2447f…` → `bb313a64…`,
`view-5k` `4da06c5d…` → `88f1774f…`; `publish` 3 m 29.1 s / 430 MiB, `verify`
3 m 31.9 s / 427 MiB, both PASS (`log/view_publish_verify.json` and
`log/view_verify.json` are those runs). The gitignored corpus payloads were
materialised into the worktree by symlink to the main checkout — the documented
#195 trap — and died with the worktree; no item `source/` or `local/` byte was
written.

**Observation (not a rendering change).** The six-term `rank_sum` is nearly
constant for the 12,882 `f2000 = 0` candidates: a gene outside a disease's
`DEG2000` pays `(n_d + 1) / n_d` in every term, and `AF`'s `n_d = 1` makes its
term exactly 2.0 for every non-AF gene, so those candidates pile into one bin at
≈7.0 — the figure's single spike. The figure plots the committed values
honestly; the thinness is the rule's arithmetic (ADR 0071), not the renderer.

**Health.** On the merged tree, `python -m unittest discover -s
stage2/01_bases/execution/tests` gives **398 tests: 3 failures + 5 errors**,
exactly the eight #209-owned tests (six `DocxTests`, two
`test_gene_axis_source_metadata.CoverageContractTest`, listed in
`inputs.json.health.red_tests`); the two #210 figure replay tests are green and
every other suite is green. `script/check_paths.sh`, `script/check_tracked.sh`
and `script/verify_all.sh` PASS.

## 11. Rule document v2 rebuild (#209, ADR 0071)

Accepted 2026-11-15 on branch `ticket/209-docx-v2` (base `dev=e1a08da`, the tree
that carries #210's re-pinned figure bytes). Contract: ADR 0071. The document is
still a reading aid, not a manifest output, so no release product, manifest, axis
or view changed.

**What changed.** `execution/make_gene_axis_docx.py` is rebuilt around the v2
rule: §1 取用规则 opens the document (DEG = `padj < 0.05` and `|log2FC| ≥ 2`,
pass-any-contrast; DEG-only `1..n_DEG` source rank with `rank_norm = rank/n_DEG`;
level-1 per-disease fusion; tier 1→2→3 collection with the axis-resolved 2,000
stop; `f2000`/`rank_sum` level-2 fill with the tenv-HVG fallback). The retired
`pN`, the measured-gene rank and the five-category phrasing are gone; `pN` appears
only in the "已退役" explanation. 表 1 is the source dataset table (27 datasets +
tenv-HVG, 28 rows), 表 2 is the new per-disease table (source-dataset count, tier
path, DEG union, `DEG2000`, and a note column that states the cut or the
shortfall), §4 the four misreading clarifications, §5 the coverage table in the
published long format. The audited-and-excluded AF block and its closing
paragraph are removed because `AF` returns; 表 2 records AF's single-gene
`DEG2000` (`SLC9A7`, named from `output/disease-fusion/af.deg2000.csv` through
the frozen axis identity table) and the three AF sources with `n_DEG = 0`.

**Inputs.** The generator now reads 12 tracked inputs: the coverage CSV, the
figure PNG, the six `output/disease-deg/<category>/catalog.json` files,
`output/disease-fusion/catalog.json` (表 2's source, cross-checked against each
category catalog's `n_deg`/`deg2000` and `len(datasets)`),
`output/gene-identity/genes.csv` (the frozen identity table the catalogs declare
as `axis.path`) and `output/gene-axis/manifest.json`. The manifest guard refuses
any manifest whose `release.algorithm` is not `f2000-deg2000-v2`, so a pre-ADR
0071 release cannot be rendered under the v2 rule. `make_gene_axis_docx.py`
gained `--fusion-catalog`, `--axis` and `--manifest`; the six `--catalog` and the
`--metadata` overrides are unchanged.

**Coverage denominators.** `make_gene_axis_docx.COVERAGE_COLUMNS` is now
imported from `bases.gene_axis.COVERAGE_COLUMNS` (the published long format), and
the build refuses a coverage row whose denominator is not the v2 one: `N` for
`tenv_hvg`, the disease's own `|DEG2000_d|` otherwise (nHF 1,661 / ICM 1,647 /
AF 1 / genetic_cardiomyopathy 1,482 / cardiotox_TOP2 2,000 / cardiotox_other 297).

**Source metadata.** `execution/gene_axis_source_metadata.py` grows from 19 to 28
records (`gse35108__dcm_ipsc_cm`, `gse165303__dcm_lv_bulk`,
`gse116250__hf_lv_bulk`, `gse141910__hf_lv_bulk`, `kuppe-mi__iz`,
`gse255612__af`, `gse138252__af`, `gse238242__af`, `gse269882__af`). Accession,
`n_samples` and the `GPL` record were re-read from NCBI E-utilities `esummary`
(`db=gds`) and the sample-size reconciliation from each item's tracked
`local/deg/*_deg.json` / `author_log.json`; 18 of the 28 records now carry a
`sample_size_note`. Two item labels disagree with the GPL record they cite and
the `platform_note` keeps both readings: `GSE165303`'s DEG summary calls
`GPL20795` "Illumina NovaSeq 6000" while the record is **HiSeq X Ten**, and
`GSE116250`'s calls `GPL16791` "Illumina HiSeq 2000" while the record is
**Illumina HiSeq 2500**. `platform` carries the GPL record; the item label stays
visible in the note. `Kuppe-MI` is a cellxgene object with no GEO platform.

**Digest and replay.** The committed
`output/gene-axis/gene_axis_fill_rule.docx` replayed byte-identically from its
inputs: `2bbe1f0e1fc8e328c90288f27c1e338d9c50e564eefe30c2b012fece0091ed88`
(92,412 B, #198/#199) → `e7fbc26a38ad005297989096ea2e4a5b087a69e7daf11d95757edbb7b6d37c61`
(94,590 B). The build is byte-deterministic: three runs (two 2 s apart, one from
`/tmp` with all absolute paths) produced the same digest, the DOCX zip is
re-timed by `normalize_docx`, and the `python-docx` template's fixed
`docProps/core.xml` timestamps (2013-12-23) leave no clock in the bytes. The
replay is pinned by `test_gene_axis.DocxTests.test_committed_docx_replays_from_its_inputs`,
which rebuilds from the 12 inputs into a temp dir and compares sha256.

**Tests.** `test_gene_axis.DocxTests` and
`test_gene_axis_source_metadata.CoverageContractTest` are re-pinned to the v2
release: 28 source rows, 27 published dataset ids + tenv, 18 reconciled
sample-size notes, the long-format coverage contract with the v2 denominators, the
per-disease 表 2 contract, and a guard that the v2 rule is the document's opening
section. The `test_af_excluded_table_contract` case is retired with the block and
replaced by `test_disease_table_matches_the_published_catalogs`. Test inventory:
**402 tests, 0 failures, 0 errors** (398 at base, −1 retired AF test, +2 表 2 /
rule-section tests, +3 metadata spot checks); `test_gene_axis.py` and
`test_gene_axis_source_metadata.py` together 11.0 s. The DOCX is not a manifest
output, so `build_gene_axis.py --verify`, `build_view.py verify` and the two
figure replays are untouched by this ticket. `script/check_paths.sh`,
`script/check_tracked.sh` and `script/verify_all.sh` PASS on the merged tree.

**Flags (not decided here).** (1) The #208 handover and spec #200 say "23
published datasets"; the released catalogs publish **27** dataset rows (18 before
#202–#206, +5 from #202/#203/#206, +4 AF from #204), which is what
`metadata.validate` is driven by. The number in the handover is a pre-AF count.
(2) `test_build_reads_the_five_catalog_paths_it_is_given` keeps its name (the red
list in `inputs.json.health.red_tests` refers to it) although it now reads six.
(3) The two item-labelled platform disagreements above are corrected in the
document metadata only; the `disease_data/` item summaries still carry them and
were not edited (out of this ticket's scope).
