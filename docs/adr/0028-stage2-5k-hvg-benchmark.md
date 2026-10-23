# Stage 2 5,000-gene HVG sensitivity benchmark

**Status:** accepted (2026-10-22)

ADR 0027 froze the cold-drug benchmark at all 18,164 corpus genes and explicitly
excluded HVG selection. That contract answered one question: how do the seven
cold-drug-capable models rank on the whole transcriptome? It left a second
question open: **is that ranking an artefact of the gene space?** Most author
pipelines select 5,000 HVGs, and Squidiff did not converge at 18,164 genes. This
ADR defines a deliberately different, single-axis sensitivity run. It does not
amend, supersede or reinterpret ADR 0027.

## Decisions

### 1. A separate archive, not a change to the 18k contract

- `stage2/03_bench_models/temp_env_all_5k/` is a fork of the executed
  `temp_env_all_18k/` archive.
- `temp_env_all_18k/` stays read-only. ADR 0027 is not edited. The 18k
  leaderboard remains the primary result; the 5k table is a sensitivity view.
- Everything except the gene axis is held fixed: wells, Δ, controls, groups,
  five cold-drug folds, the 12 entrants, seeds, budgets, adapters and metrics.

### 2. Gene-axis selection

The axis is the 5,000 genes with the smallest sum of per-source HVG ranks over
the three curated CM views in `stage2/02_cm_data/datasets-1/`:
`dtoxs_cm`, `screen_ctox`, `temp_env_all-smiles`.

Per source:

1. `scanpy.pp.highly_variable_genes(flavor="seurat", n_top_genes=None)` on
   `layers["log1p_cp10k"]`.
2. Genes with finite mean in `[0.0125, 3]` are ranked `1..M` by descending
   `dispersions_norm`, ties by corpus order.
3. Every other gene, including zero-mean genes whose dispersion is NaN, gets the
   worst rank `18,165`.

The three ranks are summed; the 5,000 smallest sums are kept, with ties broken
by corpus order, and written in corpus order. The rule deliberately prevents
low-expression/NaN genes from masquerading as highly variable. The complete
per-gene table and the input hashes are recorded by
`scripts/build_gene_set.py` in `data/gene_set_5k_manifest.json` (see also
`data/gene_set_5k_ranks.csv`).

The selection is unsupervised and uses **all wells**, including wells of drugs
that are test in some fold. This matches the authors' pipelines and keeps one
fixed gene axis across folds; per-fold gene axes would make prevalence and
cross-fold summaries incomparable. No expression labels or Δ values enter the
selection. This is declared, not hidden.

### 3. Data derivation and parity

`scripts/build_data.py` rebuilds from the same source h5ad, keeps only the
selected columns, and recomputes `true_sets_k{100,200,500}` inside the 5k axis.
Before writing the manifest it verifies that every rebuilt array equals the 18k
array at the same gene columns, and that `wells.parquet`, `groups.parquet` and
`splits.parquet` are identical to the 18k archive. The five folds are therefore
paired with the 18k run, per fold and per drug.

### 4. Metrics and prevalence

The metric definition is unchanged: direction-macro, group-macro AUPRC at
k = 100 / 200 / 500 absolute genes, with k = 200 as the headline. Because the
denominator shrinks by 3.6×, random prevalence rises from 1.1% to **4.0%** at
k = 200. Absolute k is kept so the two runs ask the same "top-200 genes"
question; the higher chance floor is reported explicitly and is visible in the
`no_change` row.

### 5. Entrants and execution

The same 12 registry entries (5 CPU baselines + 7 real models) run the same five
folds under the same GPU lock, 8 GB RAM scope and per-unit budget. No model is
tuned, dropped or added for the 5k axis, and no seed changes. A model that
converges at 5k but not at 18k is evidence about the gene space, not a new
entrant.

## Consequences

- A paired 5k-vs-18k comparison is possible for every model and fold.
- The 5k leaderboard answers only whether the ranking is robust to a
  five-fold smaller, HVG-selected gene space. It does not replace the 18k
  conclusion and must not be cited as the benchmark's primary result.
- Any future change to the selection rule requires a new ADR; the two archives
  and their manifests are immutable evidence.
- The 5k run is not a test of the authors' exact preprocessing: the HVG
  statistic is Seurat on cp10k log1p, chosen to avoid a new dependency, and it
  is shared across all three views rather than fitted per training fold.
