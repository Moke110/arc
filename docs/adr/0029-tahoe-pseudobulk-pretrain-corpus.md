# Tahoe-100M enters Stage 2 as a pseudobulk pretrain corpus

**Status:** accepted (2026-10-30)

**Update (2026-10-30, ADR 0030):** paths in this record were rewritten by the external-data / step-layout restructure; the decision itself is unchanged.

Stage 2's CM corpus is bulk — one physical library per well — frozen at 18,164
genes (ADR 0026). The candidate pretrain corpus is Tahoe-100M: 97,840,590
QC-passed cells from 50 cancer cell lines, in which **every well pools most of
the 50 lines**. Feeding single cells to a model that must later emit bulk CM
profiles would hand it inputs roughly an order of magnitude sparser than its
finetune data. Tahoe therefore enters the repository as **pseudobulk profiles**,
and this ADR freezes how they are built. It does not choose a model: the corpus
is model-agnostic by construction.

## Decision

One corpus, built by `stage2/04_tahoe_data/`:

1. **Unit = `(plate, well, cell_line)`.** Upstream `sample` *is* the well
   (14 plates × 96 wells = 1,344). `drug` and `dose` are well-level attributes,
   not key components: no well carries more than one perturbagen, and the
   compound set is 379 drugs × 3 doses (0.05 / 0.5 / 5 µM) plus 28 `DMSO_TF`
   wells (2 per plate). The unit is one notch finer than the authors' shipped
   DE, which pools the 22 `(plate, drug, dose)` pairs that occupy two wells —
   our QA must fold those 22 pairs back to compare.
2. **Gene space = Tahoe native 62,710, not the 18,164 CM intersection.** A
   corpus freezes identity, a bundle applies filter rules (ADR 0026); 62,710 is
   Tahoe's identity, 18,164 is CM's. Freezing at 18,164 would bake an
   as-yet-unmade model decision into 428.8 GB of input data. Projections
   (18,164 / 978 landmarks / HVG) are views, and must go through **Ensembl IDs**
   (18,163 of 18,164 CM genes map) rather than symbols (17,935) — the 229-gene
   gap is alias/previous-symbol drift, not absence.
3. **Values = summed raw counts**, plus a derived `log1p_cp10k` layer. No DESeq2
   size factors, no quantile alignment, no HVG selection — those are
   model-specific (ADR 0026 §2).
4. **Controls are rows.** `DMSO_TF` wells are profiles too, with
   `control_group = tahoe::plate`. **Δ is not materialized**: it is
   `log1p_cp10k(treated) − mean(same-plate, same-cell-line DMSO_TF)`, one
   groupby away. Materializing it would freeze a pairing rule that belongs to
   the eventual model.
5. **No silent filtering.** Every unit with ≥ 1 QC-passed cell is a row and
   carries `n_cells`; the ≥ 100-cell cut belongs to a view, not the corpus.

**Rejected.** (a) Summing a whole well into one bulk profile — a 50-cell-line
mixture, whose DE is a Simpson's-paradox artefact and whose control is a mixture
too. (b) Using only the shipped `pseudobulk_differential_expression` — it
carries LFC/Wald statistics, not levels, so it cannot be re-normalized or
re-filtered; it is retained as a **QA yardstick** instead. (c) Freezing at
18,164 genes. (d) Materializing Δ.

## Why pseudobulk, measured

Against the four CM matrices restricted to the same 18,164-gene space
(2026-10-30, `perf_data` raw counts vs Tahoe pseudobulk, bytes on disk):

| profile | cells | genes detected | zero fraction | counts |
|---|---|---|---|---|
| CM `dtoxs_primary` well | 1 library | 10,604–10,774 / 17,112 | 0.38 | 0.73–0.86 M |
| CM `screen_ctox` well | 1 library | 11,066–11,320 / 18,164 | 0.377 | 0.35–0.43 M |
| CM `temp_env` well | 1 library | 6,920–7,569 / 20,301 | 0.63 | 0.35–0.61 M |
| Tahoe pseudobulk, median | ~1,111 | 14,450 / 17,935 | 0.19 | 0.37 M |
| Tahoe pseudobulk, > 1,500 cells | > 1,500 | 15,602 | 0.13 | 1.3 M |
| Tahoe pseudobulk, < 100 cells (9.5 %) | ~50 | 11,485 | 0.36 | 0.04 M |

Pseudobulk is **denser** than a CM well at every size above 100 cells and
depth-matched to `screen_ctox`; the 9.5 % tail matches CM sparsity but is
10–20× shallower — which is exactly why the view threshold is 100 cells. The
authors' own headline analysis is at pseudobulk level ("we assessed the
transcriptional impact of treatment at the level of pseudobulk profiles",
*Cell* 2026, doi 10.1016/j.cell.2026.08.035), so this is the sanctioned
aggregation, not an invention.

## Consequences

- **Pseudobulk fixes sparsity, not biology.** The corpus is cancer cell lines at
  a single 24 h timepoint with 24/525 drug-key overlap against the CM corpus.
  Those remain pretrain→finetune domain gaps and must be stated wherever the
  corpus is used, not engineered away.
- **Provenance is an import, not a fetch.** The atlas arrives from a Windows-side
  HF snapshot (`E:\Desktop\tahoe\data\hf`) and is verified against HF's
  `lfs.sha256` before the source copy is deleted. `fetch_tahoe.py` stays as the
  verifier and full-fetch path; the 36-shard vendored subset it originally
  produced is a strict subset of the imported tree and is superseded.
- **Imported tree contradicts two vendored notes.** `Tahoe_summary.md` and
  `DATASET.json` record `pseudobulk_included: false`; the imported tree contains
  all 1,026 DE shards. Both files are rewritten as part of the import.
- **QA gates.** Median detected genes ≥ 12,000 per profile; every
  `(plate, cell_line)` pair carries `DMSO_TF` cells; and — on genes with the
  authors' `baseMean > 500` — the median per-contrast Pearson between our
  `log2((treated+1)/(same-plate same-cell-line vehicle mean+1))` and their
  DESeq2 `log2FoldChange` is ≥ 0.60. A corpus that fails any gate is not
  promoted.

  This gate replaces the originally proposed "LFC Spearman ≥ 0.8 against the
  shipped DE", which measurement showed to be the wrong statistic twice over.
  First, Spearman over all 62,710 genes sits at ~0.16 because most genes have no
  real effect in a single well, so both estimates are noise and their rank
  agreement is uninformative; the correlation rises monotonically with
  expression and only becomes meaningful above `baseMean` ≈ 500. Second, even
  there agreement tops out at Pearson ≈ 0.7–0.88 depending on cell line
  (A549 0.88, A-427 0.77, C32 0.68 on plate 4), because the shipped DE is not a
  naive count-sum pseudobulk: its per-gene `baseMean` sits 1.5–2.8× below our
  summed counts by a factor that varies gene to gene, so it was produced under a
  normalization/aggregation this repository cannot reconstruct. The gate is
  therefore a **direction-and-magnitude concordance floor**, reported per cell
  line, not a claim of byte-identical reproduction; the full threshold sweep
  lives in `qa_report.md`.

  Measured on the full corpus (2026-10-30, 3,755 contrasts sampled from 60 of the
  1,026 DE shards): Pearson **0.779** at `baseMean > 500` (644 genes/contrast),
  **0.812** at `baseMean > 800`, Spearman 0.713 and 0.173 over all genes; median
  25,290 detected genes per profile; all 700 `(plate, cell_line)` pairs carry a
  vehicle. All gates pass. The plate-4 pilot alone read 0.706–0.715, which is why
  the gate floor was not set from a single plate.
