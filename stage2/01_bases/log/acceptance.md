# Acceptance record — Stage 2 `01_bases` rebuild (#192, ADR 0066 / ADR 0068 / ADR 0069 / ADR 0062)

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

**2026-11-15 rule-document source table (#198).** The rule document
`output/gene-axis/gene_axis_fill_rule.docx` gained a source dataset information
table as its first section: 18 rows (one per released dataset, in the five
catalogs' order) plus a tenv-HVG fusion row, columns `category / accession /
sample type / groups / sample sizes / sequencing technology / CM DEG count`, and
one footnote defining the CM DEG count (`catalog.json` `n_significant`, best
contrast, BH padj<0.05), the cardiotox vehicle-vs-compound and isogenic contrast
semantics, and the `GSE255612` GEO `overall_design` 19/19-vs-18/16 mismatch (the
release's numbers are unchanged). The existing five sections are renumbered
2..6. The document was regenerated twice byte-identically from its inputs; new
sha256 `044f7ece07d57a8130b5ea3d5b5ad0941a18df35e16fc1c3fdaf5ba08d8d0402`. The
generator now reads the five `output/disease-deg/<category>/catalog.json` inputs
plus `execution/gene_axis_source_metadata.py` (#197): ADR 0070's Consequences
line says "two catalog inputs", which was an early estimate — the implementation
uses five (recorded here; the accepted ADR is not edited). No release product,
manifest or view changed: the docx is a reading aid, not a manifest output.
