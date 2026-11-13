# Acceptance record — Stage 2 `01_bases` rebuild (#192, ADR 0066 / ADR 0062)

Accepted 2026-11-13 on branch `ticket/192-acceptance` (worktree `.scratch/worktrees/t192`,
base `dev=cb72b92`). Environment: `/home/chang/miniconda3/envs/arc/bin/python`
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
| `output/raw_counts/` | `build_raw_counts.py verify` (`source-accepted-replay`) | `6d45b171…cbff4f` | PASS | 14.4 s | 379 MiB |
| `output/corpus/` | `build_corpus.py verify` | `f8c35ccc…9d5d9f6d` | PASS | 1 m 10.5 s | 462 MiB |
| `output/disease-deg/{af,icm,nhf,genetic_cardiomyopathy}/` | `build_disease_deg.py verify --all` | AF `397c37c1…`, ICM `f05936a8…`, nHF `9cc58c08…`, genetic `18590976…` | PASS (4/4) | 8.5 s | 286 MiB |
| `output/disease-deg/{cardiotox,cardiotox_TOP2,cardiotox_other}/` | `build_disease_deg_cardiotox.py verify` (R: DESeq2/edgeR) | `06e4e08b…`, `496b4977…`, `57529b8c…` | PASS (3/3) | 1 m 39.9 s | 1.47 GiB |
| `output/gene-axis/` | `build_gene_axis.py --verify` | `42baa74f…aeb6fcb` | PASS | 0.5 s | 138 MiB |
| `output/view-2k/` | `build_view.py verify` (payload published first) | `b48eeddb…040c18db` | PASS | 3 m 38.4 s (both views) | 435 MiB |
| `output/view-5k/` | `build_view.py verify` | `37af263b…f2f6ae1e` | PASS | (both views) | (both views) |

The two view payloads are gitignored and were absent from the fresh worktree, so
`build_view.py publish` materialised them first (3 m 32.5 s, 428 MiB); the
published manifests were **byte-identical to the tracked ones** (`git diff` empty),
and the subsequent `verify` replayed both views from the pinned corpus + axis.
`output/tenv-hvg/` is a retained pre-0066 product, not rebuilt by this ticket
(manifest `953103a1…`; see README *Retained tenv-HVG ranking*).

Record files: `raw_counts_verify.json`, `corpus_verify.json`,
`disease_deg_verify.json`, `cardiotox_verify.json`, `gene_axis_verify.json`,
`view_publish_verify.json`, `view_verify.json`; the machine-readable digest
inventory is `inputs.json`.

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
  `n_measured + 1`, `rank_norm = rank/n_measured`, `gene_index` tie-break)
  re-verified against the pinned item DEGs; catalogs' self-consistent
  `catalog_sha256`; level-1 inputs are the released rank frames.
- `gene-axis`: level-1 dataset fusion inside each category, level-2 seven-source
  `rank_sum`, `pN ∈ 0..7`; axes exactly 2,000/5,000 ascending **original**
  `gene_index` (never renumbered), fill `0.5·N` tenv-HVG + `0.5·N` pN/rank_sum
  remainder; non-nested axes with reported intersection **1,919**; manifest
  digests and 18,164-row candidate table.
- `view-*`: normalization on the **full-axis** library total then subset (never
  the view-local denominator); `control` = mean of the named matched-control
  rows; `Δ = X − control` exactly; `X/control/Δ` float32; obs carries row
  identity + `library_total`; no `{N}k_eval_truth/` anywhere.

## 2. Downstream reader compatibility

Every retained reader module was exercised through its own test file under the
`arc` env (all green): `test_raw_counts` (8), `test_corpus` (10), `test_view`
(17), `test_gene_axis` (38), `test_deg_rank` (12), `test_disease_deg_human_sc`
(18), `test_disease_deg_isogenic` (35), `test_drug_classes` (17),
`test_disease_deg_cardiotox` (11), `test_identity`, `test_fold_seed`,
`test_palette`, `test_hvg` (2, after removing the retired metrics-QA assertion), `test_native`,
`test_tahoe`, `test_raw_provenance`.
`02_models_tenv` is deliberately **not** re-pointed (ADR 0066 *Deferred*); its
target contract (`view-{N}k` `X`/`control`/`delta`, treated-row selection) is
documented in the README, and its stale imports of the deleted evaluation
modules are the deferred work item, not a live consumer of any accepted product.

## 3. Cross-category catalog invariant (new #192 gate)

`execution/tests/test_acceptance_invariants.py` (8 tests) passes over the six
released catalogs. It pins the canonical `stage2-disease-deg-catalog-v1` shape
(`datasets` list of full records, per-record `method`, self-consistent 64-hex
`catalog_sha256`, no legacy `test_method`/`dataset_records`), positive
`n_case`/`n_control` **or** an explicit `*_basis`, 64-hex `derived_inputs` where
declared, and — explicitly, not silently — the documented cross-ticket variance:
#187/#188 declare the published `gene-identity/genes.csv` axis, a
`rank_contract`/`rejected` block, bare-digest `inputs` and gitignored
`derived_inputs`; #189 declares the `tables/genes.csv` axis, `{bytes,sha256}`
`inputs` and no `derived_inputs`. A drift to a third shape fails the gate.

## 4. CSV-only figure reproduction (products unreadable)

With every `output/` product directory made unreadable (`chmod a-r`; only
`output/api` readable), the renderer
`figure/script/render_gene_axis_rank_sum_histogram.py --csv
figure/data/gene_axis_rank_sum_histogram.csv` re-rendered the committed PNG
byte-identically: `49d1dc5d…c1f2581` = `49d1dc5d…c1f2581` (exit 0). The retired
diagnostic PR figure and its renderer are deleted.

## 5. "No retired path read" evidence

Input sets were read from the accepted manifests (`inputs.json`), not assumed.
`output/raw_counts/` is built only from `gene-identity/{genes,crosswalk}`,
`pert_data/GSE262419` + `pert_data/Tahoe-100M` source metadata, and the two
pre-0066 source-replay acceptance records. `output/corpus/` reads only
`raw_counts/` payloads + `gene-identity` + Tahoe drug metadata.
`disease-deg` reads item DEG JSONs/CSVs under `disease_data/`/`pert_data/` and
`gene_identity`. `gene-axis` reads `disease-deg/`, `tenv-hvg/descriptive.csv` and
`gene-identity`; `view-*` reads `corpus/` + `gene-axis/` only. No retained
builder or reader opens `shared-expression/`, `training-corpora/`, the old
raw-count directories, `tenv-truth*` or `historical-training-release/`
(verified by reviewing every manifest `inputs` key and by
`git grep` over retained code after deletion).

The one accepted, one-way exception is documented instead of hidden:
`raw_counts`' default **rehome** mode pins the pre-0066 raw-count payload bytes
and the two pre-0066 source-replay records as its inputs. Those bytes are
preserved verbatim in `raw_counts/` (identical sha256: tenv
`efd07632…`, Tahoe `a830a34a…`) and inside git history, so the retired paths are
not read by any live product; a from-scratch replay uses the documented stronger
mode `build_raw_counts.py --full-rebuild publish|verify` (Tahoe build 1 h 35 m,
independent source replay 2 h 35 m). The retired evidence digests are preserved
in `inputs.json.retired_evidence_digests`. The retained `tenv-hvg` ranking has
the analogous deferred re-point documented in the README.

## 6. Old → new mapping (deletions)

| removed (commit group) | accepted substitute / rationale |
|---|---|
| `output/shared-expression/` (products) | `output/corpus/` integer counts + view normalization |
| `output/training-corpora/` (products) | `output/corpus/` + `output/view-2k/`/`view-5k/` |
| `output/tenv-raw-counts/`, `output/tahoe-raw-counts/` (products) | `output/raw_counts/` — same payload bytes, one manifest |
| `output/tenv-truth/`, `output/tenv-truth-successor/`, `output/historical-training-release/` (products) | none; evaluation redesign deferred; each view will own `{N}k_eval_truth/` |
| old `output/gene-axis/` file set: `source_summary.csv`, `fused_rank.csv`, `overlap.csv`, `inflection.json`, `tenv_non_hvg_curve.csv`, `direction_summary.csv`, `mouse_human_orthologs.csv` (product) | `candidate_rank_sum.csv`, `gene_axis_{2000,5000}.csv`, `axis_composition.csv` (retired fusion semantics not resurrected) |
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
from an interim pre-#190 axis. On the **released** axes, measured from the view
manifests (tenv, all 4,478 profiles): **G_2000 = 12.9299987%** and **G_5000 =
37.2063189%** (treated-only 12.9566820% / 37.1590893%; `tahoe_bulk` 10.4600758% /
24.7819987%). The 2.5× baseline-change argument and the full-axis decision are
unaffected (measured denominator ratio ≈2.88×); the released values are
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
