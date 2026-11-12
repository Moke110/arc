# Prepare the Stage 2 gene-axis decision without deciding the axis

**Status:** accepted (2026-11-13).

The Stage 2 exit is now fixed as **the predicted Δ on a selected gene axis, as accurate as
possible** (`docs/roadmap.md`). The frozen 18,164-gene axis (ADR 0061) is not that selected
axis. Narrowing it is an open decision (#173), tracked on the decision map #171; this record
fixes only how the evidence for that decision is produced.

The user re-scoped the work on 2026-11-13: this round is **decision preparation, not the
decision**. The public myocardial-disease transcriptomes the decision needs are ingested into
`disease_data/` under the ADR 0030 external-item contract, the cardiotoxicity evidence the
repository already owns is integrated read-only from `pert_data/`, and the HVG computation,
rank-sum fusion and per-source summary tables are built in `stage2/01_bases/` under the ADR
0062 step layout. No new corpus, no `genes.csv` release, no model, no truth, no Stage 2
re-run.

## Boundary: preparation, not release

The outcome of this decision is an **evidence product**, not an axis. It reports per-source
variability and differential-expression evidence on the unchanged frozen gene space, the
source overlap, a declared inflection answer (including an honest "no inflection"), the tenv
non-HVG share and the top-`N` direction composition.

Which is deliberately left open, in #173:

- the axis size `N`;
- a per-source quota;
- the direction rule;
- the weighting of animal sources relative to human ones.

#174 owns any gene-axis release (`genes.csv` or otherwise), and a release is out of scope
here. The old axis, every accepted product, the corpora, the measured truth and Stage 1 stay
byte-identical (ADR 0028/0061 are not amended, only used).

## Sources

Five open GEO series are ingested, each as its own external 项 under `disease_data/`:

| accession | disease / role | species, model | why it is in the set |
|---|---|---|---|
| GSE141910 | heart failure, primary anchor | human, LV bulk RNA-seq | largest verified human HF LV transcriptome (200 HF / 166 NF) |
| GSE116250 | heart failure, cheapest anchor + human post-MI proxy | human, LV bulk RNA-seq | 50 HF (37 DCM, 13 ICM) / 14 NF; its ICM arm is the smallest human post-MI-proxy LV download |
| GSE243674 | chemotherapy cardiotoxicity, independent | human, hiPSC-CM | 72 wells, anthracyclines + vehicle; the only new series whose cell model matches what Stage 2 predicts a Δ in |
| GSE128188 | atrial fibrillation (>5% prevalence category) | human, LA + RA appendage | 5 AF / 5 SR; the AF category has no other processed series locally |
| GSE282618 | hypertensive heart disease (exploratory) | human, LV + RV | 3 vs 2 subjects; the only HHD series found, and explicitly flagged exploratory |

Local animal MIRI joins as gap-fill only, with its ortholog chain recorded: the mouse
`disease_data/GSE160516/` item and the corrected `disease_data/GSE108940/` item (mouse MuRF1
I/R; the #176 provenance correction). Both are **DEG-only** sources — no expression matrix is
shipped for them, so they contribute a `deg_rank` only.

The cardiotoxicity evidence is integrated from data the repository already owns, read-only:
`pert_data/DToxS-iPSC/` (GSE217424 full-transcriptome counts and per-drug DEG calculation
tables; anthracyclines present) and `pert_data/GSE244740/` (ScreenSeq hiPSC-CM counts with
`DMSO` vehicle and doxorubicin / idarubicin / mitoxantrone). `pert_data/` bytes are never
modified and nothing is re-downloaded. `pert_data/LINCS-LDS-1593/` is a **coverage note
only**: it is L1000-978, and that ceiling is exactly what this work replaces, so it is never
a fusion source.

`disease_data/GSE137293/` (pig scar-versus-border zone) is listed as **supplementary evidence
only** and never enters the fusion: it contrasts two regions of the same injured heart — a
severity gradient, not a disease-versus-healthy contrast.

The **optional replication tier is deliberately excluded**: GSE57338 (1.3 GB CEL),
GSE5406 (series matrix only), GSE36961, GSE262690 (915 MB RDS), GSE79768 (125 MB CEL) and
GSE1145 (665 MB CEL, no series matrix). The dataset survey
(`docs/research/stage2-stage3-gene-axis-datasets.md`) says not to fetch them before the
axis-selection rule is fixed, they add no disease category the five series do not already
cover, and each adds a different platform and preprocessing cost to a decision-preparation
round whose job is evidence, not exhaustive coverage. They stay a documented option for a
later decision round.

## Where the artefacts live

`disease_data/` holds **measurements**; `stage2/01_bases/` holds the **evidence computed on
the frozen axis**. Neither owns the other's content.

Each new series is a full ADR 0030 项 — `script/`, `source/` (immutable download),
`local/` (declared additions), tracked `author_log.json` and `manifest.sha256` over
`source/` only — with these item-relative artefacts:

- `local/deg/<GSE>_deg.csv` — the full-transcriptome DEG table (gitignored, regenerable);
- `local/deg/<GSE>_deg.json` — the tracked compact summary (counts, method, thresholds,
  input digests);
- `local/expr/<GSE>_log_expr.tsv.gz` + `local/expr/<GSE>_samples.csv` — the normalized,
  gene-level log expression matrix and its case/control sample sheet. The DEG table is the
  source's interface for its DEG rank and the expression pair its interface for its HVG
  rank; the axis step reads nothing else from the item.

The ADR 0030 rules apply unchanged: `source/` is never edited, every derived path is declared
in `local_additions[]`, `manifest.sha256` covers `source/` only, and `script/verify.sh`
verifies the download offline.

The evidence product is `stage2/01_bases/output/gene-axis/`, with reusable code, driver,
tests, figure and records at the ADR 0062 paths: `output/api/bases/gene_axis.py` (pure
functions), `execution/build_gene_axis.py` (build and `--verify` replay),
`execution/gene_axis_sources.py` (the `pert_data/` adapters), `execution/tests/`,
`figure/gene_axis_source_summary.png` with its tracked `figure/data/` plotted values and
CSV-only renderer, `log/gene-axis/` records, and `manifest.json` pinning every input path and
digest. Compact evidence files are force-added like the other accepted evidence; the rest of
`output/gene-axis/` stays untracked payload.

## Full-transcriptome differential expression, never L1000-filtered

Every new DEG table is **full-transcriptome and gene-level**, carries a signed effect
(`log2FC`), test statistic, p-value, BH `padj`, expression level and case/control `n`, and
declares its contrast (disease versus control; positive `log2FC` = up in disease). A
multi-contrast source keeps one long table with a `contrast` column.

The existing `disease_data/*/local/deg/*.json` lists are L1000-978-filtered and remain
Stage-1-only. They are never reused as axis sources, and no new artefact is filtered to the
978 landmarks: an L1000-filtered list cannot rank 18,164 genes, so reusing one would make the
whole preparation circular.

## The frozen gene space and the tenv input

The evidence projects onto the unchanged frozen 18,164 `gene_index` rows
(`stage2/01_bases/output/gene-identity/genes.csv`, ADR 0061). No new identity table, no
reordering and no new axis are created; source symbols and aliases are resolved against the
identity tables with the resolver and its one-to-many/lost counts recorded.

The only tenv input is `stage2/01_bases/output/tenv-hvg/descriptive.csv`, the published
descriptive HVG ranking (ADR 0061) — it is read, **not recomputed**. Because this product is
evidence about the axis rather than a model-facing selection, it does not consume or replace
the train-only rankings that ADR 0061 reserves for model gene selection.

## Source ranks and the rank-sum fusion

Per source, and **within that source only**:

- `hvg_rank` — where an expression matrix exists: Seurat-style normalized dispersion on the
  item's normalized log expression (sample variance/mean, log dispersion z-binned in 20
  equal-width mean bins, unaffected/singleton bins scored as in the tenv recipe,
  `gene_index` tie-break). Same family as the published tenv recipe; deliberate differences
  are documented beside the implementation.
- `deg_rank` — from the continuous effect score `|log2FC| × -log10(padj)` over the source's
  measured genes; a multi-contrast source collapses per gene by the **maximum over contrasts**
  and records `n_contrasts_supporting`.

A gene the source does not measure gets the worst rank `n_measured + 1`; it is never
silently dropped. Ranks are never compared across sources: a microarray effect size and an
RNA-seq effect size are not the same quantity.

The fusion is **unweighted and quota-free**, as the 2026-11-13 instruction requires:

- `rank_norm = rank / n_measured` puts sources of different coverage on one scale, and an
  unmeasured source contributes `rank_norm = 1.0` (the worst rank);
- `rank_sum = Σ_sources rank_norm` over all sources;
- the primary fused ordering is ascending `rank_sum`, ties broken by the frozen
  `gene_index`.

No source weight is applied and no per-source quota is imposed. The animal-source
half-weighting variant is written as a **sensitivity column only** — it is a decision input
for #173, not a choice made here. The build is deterministic: the same inputs produce
byte-identical CSVs, and `manifest.json` records the recipe string and every input digest.

## Reported evidence

The product answers six quantities; each is one compact tracked table (or the figure that
renders from it). The inflection answer may honestly be "no inflection": the preparation
reports what the criterion finds and never invents an `N`.

1. per-source HVG share, DEG up/down counts and shares, and axis coverage
   (`output/gene-axis/source_summary.csv`);
2. per-source rank distributions and the fused-rank distribution
   (`output/gene-axis/fused_rank.csv`);
3. source overlap — tenv HVG versus each disease source, disease sources pairwise, and the
   all-source intersection — with sizes, Jaccard and rank positions
   (`output/gene-axis/overlap.csv`);
4. inflection detection on the fused-score curve with a declared criterion (Kneedle
   `S=1.0`, online variant, plus a second-difference peak, with parameters and
   leave-one-source-out stability) and an explicit `"inflection": false` when none passes
   (`output/gene-axis/inflection.json`);
5. the tenv non-HVG share versus `N` curve for declared tenv HVG set sizes
   `K ∈ {1000, 2000, 5000}` (`output/gene-axis/tenv_non_hvg_curve.csv`);
6. top-`N` direction composition under symmetric inclusion versus `|log2FC|` weighting
   (`output/gene-axis/direction_summary.csv`).

## Leakage red line

The preparation is **label-free**: it never opens a tenv truth file (`output/tenv-truth*`)
and never reads Heliste phenotype labels. The only tenv input is the published descriptive
HVG ranking, which carries no perturbation direction (ADR 0061); the expression-only use of
tenv rows that ADR 0028 permits is not needed here. A test poisons the truth path and asserts
that opening it fails, so the red line is enforced rather than promised.

## Non-goals

- No axis decision: `N`, the source quota, the direction rule and the animal weighting stay
  open in #173.
- No corpus, no `genes.csv` release and no selected-axis product (#174).
- No model, no training, no truth rebuild and no Stage 2 comparison re-run.
- No Stage 1 change; no accepted product, corpus or truth byte changes.
- No optional GEO replication tier and no controlled-access data; no new L1000-978 lists.

## Consequences and limits

- A decision can cite this evidence, but this record chooses nothing. Every number is
  descriptive of the sources, not of biology.
- Ranks erase effect magnitude, so cross-source comparability is bought at the price of
  losing how large an effect is. That is the intended trade for an unweighted fusion.
- Animal sources are DEG-only and ortholog-mapped; the mapping's losses are recorded, not
  corrected. Half-weighting exists as a sensitivity column because the bias is real.
- The disease sources differ in what they measure: bulk LV tissue carries cell-composition
  shifts that a hiPSC-CM monolayer does not; GSE282618 is exploratory at n=5 subjects;
  GSE243674 has 72 wells. The per-source summary records these, and the ADR does not
  pretend they are interchangeable.
- The tenv ranking is descriptive over all wells, so the non-HVG share bounds coverage of
  tenv variability evidence, not model trainability under the train-only rankings.
- `check_paths.sh` audits `disease_data/` path literals. The five new item directories are
  created by the sibling gene-axis tickets, so this record cites the `disease_data/`
  category and item-relative members (`local/deg/…`, `local/expr/…`) rather than literal
  per-item paths; that keeps the path audit green with no whitelist exception before those
  items exist.

## See

- [External data provenance and the Stage 2 step layout](0030-external-data-provenance-and-step-layout.md)
- [Stage 2 accepted products and independently reproducible figures](0062-stage2-accepted-products-and-reproducible-figures.md)
- [Shared Stage 2 bases and indexed corpora](0061-stage2-shared-bases-and-indexed-corpora.md)
- [Stage 2 5,000-gene HVG sensitivity benchmark](0028-stage2-5k-hvg-benchmark.md)
- `docs/roadmap.md` — the fixed Stage 2 exit
- `docs/research/stage2-stage3-gene-axis-datasets.md` — the dataset survey behind the source set
- `stage2/AGENTS.md`, `stage2/01_bases/AGENTS.md`
- `GLOSSARY.md` — the gene-axis decision-preparation vocabulary
- Spec working copy: `.scratch/stage2-gene-axis/spec.md` (local, not tracked)
