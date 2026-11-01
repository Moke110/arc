# A CM-schema view of the Tahoe corpus, and a PRnet pretrain bring-up

**Status:** accepted (2026-11-01)

ADR 0029 froze the Tahoe corpus at 62,710 native genes and deliberately refused to
freeze pairing or projection rules: "projections (18,164 / 978 landmarks / HVG) are
views". Stage 2 also had no model that had ever been trained on that corpus. This
record captures the two decisions that close that gap — what the projection view is,
and what the first training run on it is allowed to claim. Sibling records: ADR 0026
(CM corpus and views), ADR 0027 (benchmark contract), ADR 0029 (the corpus).

## Decision

### 1. The pretrain view is schema-isomorphic to the CM views, not layout-isomorphic

`stage2/04_tahoe_data/execution/build_pretrain_view.py` emits
`output/pretrain_view/tahoe_cm_schema.h5ad`: **67,018 × 18,164** — **1,400** `DMSO_TF`
controls + **65,618** treated rows, **379** drugs, **700** `(plate, cell_line)` control
groups. Built 2026-11-01 in **175 s** cold / **80 s** warm-cache and accepted by
**nine consumer-side gates, 9/9 PASS** (`qa.json` `all_pass: true`):

- `X` = the corpus's summed raw counts restricted to the 18,164 retained columns, as
  **CSR `int32`**, nnz **1,010,892,356**; `layers['log1p_cp10k']` as CSR `float32` on
  the same sparsity pattern. The artifact is **16,319,543,984 bytes (16.32 GB)** on
  disk; its content sha256 `476812c0…bf407` (full value in `manifest.csv`) is computed
  over the CSR `data`/`indices`/`indptr` arrays, never over a densified matrix.
- `var` = the CM corpus axis verbatim (`corpus/v1/genes.csv`, lexicographic HGNC
  symbols, `entrez_id`/`ensembl_id`/`match_type`), joined on **Ensembl only**:
  **18,163 of 18,164** map; the symbol route matches **17,935** and loses **229**
  alias-drifted genes. `FAM153B` / `ENSG00000182230` keeps an all-zero slot.
- `obs` = CM's 24 columns, plus Tahoe-native extras, plus three author shims
  (`SMILES`, `dose`, `paired_control_index`).
- `uns['control_mean']` = **700 × 18,164 `float64`**, anndata `dataframe/0.2.0`,
  indexed by `control_group`.

**D0 — layout. Rejected: re-encoding the view in the CM views' dense `int32`
layout.** The plan's stated reason — "a dense 67,018 × 18,164 `int32` matrix is
4.87 TB" — was a **units error**. Dense is **4.87 GB** (`int32` `X`) and **9.74 GB**
with the `float32` layer. At **83 %** density (1,010,892,356 nonzeros of
1,217,314,952 cells) the CSR artifact is **16.32 GB — larger than either dense
form** — so the layout is *not* justified by size, and the size cost is recorded as
part of this deviation. CSR is kept because the builder streams the corpus's 14
per-plate CSR parts into preallocated CSR arrays with slice writes (a dense buffer is
never materialised) and because it matches the corpus container. Everything a
consumer asserts — shape, gene axis, value dtype, the derived layer, the control
frame — is parity; only the storage layout differs, and the content digest hashes the
CSR structural arrays rather than a densified copy.

**D1 — control keying. Rejected: the corpus's plate-level `tahoe::plate` (14
groups).** A Tahoe plate pools ~50 cell lines, so a plate mean is exactly the
Simpson's-paradox mixture ADR 0029's rejected option (a) refuses, and Δ must stay inside
one cell line. The view's key is `tahoe::<plate>::<cell_line>` — **700 groups**, every one
carrying a vehicle, **0** treated groups without a control. The key keeps the CM key's
*name* for loader compatibility; only its granularity differs. That granularity is the
view's: the corpus's own `control_group` stays plate-level.

**D2 — gene axis. Rejected: an 18,163-gene axis.** The one CM gene with no Tahoe
measurement keeps an all-zero slot, so the axis is the frozen 18,164-gene CM axis
verbatim and a later 18,164-gene adapter loads it without a remap. Dropping the slot
would shift one position for every consumer.

**D3 — provenance gaps. Rejected: inventing constants for the unsourceable
columns.** `time_h`, `cell_type`, `donor` are NaN; `atc`, `moa`, `mesh_classes`,
`annotation_source` are `""`; `pharm_atc_moa`/`pharm_mesh` are False. None of these has
a source anywhere in the atlas — notably there is no time column at all, so "24 h"
stays ADR 0029 prose and is **not** materialised as `24.0`. The columns exist so the
shape is uniform; the deviation table in `qa_report.md` names each empty field and why.

**Rejected: silently dropping the two compounds with no resolvable structure.**
`Sacubitril/Valsartan` and `Verteporfin` have no `canonical_smiles`; their rows stay in
the view, flagged `has_structure=False` (the **300** non-control rows of 2 of 379
drugs), consistent with ADR 0029 decision 5 ("no silent filtering"). The trainer drops
those rows and **counts** them.

**What the view is not.** It is **not the corpus**: the corpus is 62,710 native genes
with 14 plate-level control groups and no projection, and it is frozen; the view is
derived from the corpus's 14 per-plate parts, reproducible and disposable — the corpus
is none of those. It is **not a `datasets-1` bundle**: bundles are built from the CM
corpus and carry per-view membership booleans, which is why the view ships its own
`compounds.csv` (379 rows, `overlap_cm`) rather than extending the frozen bundle's
compound table.

### 2. The first run on the view is a bring-up, and claims nothing about transfer

`stage2/05_prnet_pretrain/` trains PRnet's `PGM` — **7,185,000 parameters** — on the
view from random initialisation: Gaussian NLL on the treated `log1p_cp10k` expression
from (control mean, rFCFP4 × log10-dose), with the **same recipe** as the 03.1 adapter
(Adam 1e-3 / wd 1e-8, batch 512, `ReduceLROnPlateau(0.5, patience 10)`, 100-epoch cap,
early stop after 20 stale epochs with a 20-epoch floor, seed 20261021), split by plate
(`plate13` + `plate14` valid).

Measured 2026-11-01 (RTX 4070 Ti SUPER, torch 2.0.1+cu118): **65,318 samples = 55,922
train / 9,396 valid**, **300** rows dropped for having no resolvable structure;
**68 epochs**, early stop fired at epoch 68 after 21 stale epochs with the 20-epoch
floor met; best valid **Gaussian NLL −2.628080 at epoch 47**; lr halved to **5e-4 at
epoch 58**; **883.2 s** wall end to end; **15,150 MiB** peak host RSS.

Its only claims are: **the pipeline runs on GPU and the loss-epoch curve exists.** A
negative Gaussian NLL is the correct value of a log-likelihood and is never clamped; the run
makes no transfer, finetune, warm-start, usefulness or benchmark-comparison claim.

**Rejected for this record:** a warm-start / finetune claim (it needs a random-init vs
Tahoe-init comparison on the cold-drug benchmark — a separate effort, out of scope);
starting from PRnet's released checkpoints (locked to 978 / 5,000 genes, and 03.1
already rejected public-weight initialisation for the benchmark); and any comparison of
this NLL with 03.1's best valid −1.174120 (its PRnet fold-0 run; different corpus,
rows and split — not commensurable).

## Consequences

- The view is **1,010,892,356 nonzeros / 16.32 GB** and local-only; the tracked record
  is `manifest.csv` + `compounds.csv` + `qa.json` + `qa_report.md`. Rebuilding it
  re-derives everything from the corpus's 14 per-plate parts — never from the raw
  atlas shards (428,811,080,261 bytes imported; `pert_data/Tahoe-100M/author_log.json`)
  — and a warm gather cache rebuilds it in 80 s.
- The corpus's identity rules now have a second consumer. Measured drug overlap against
  the CM corpus's **525** `drug_key`s: **31** by the corpus's desalted-parent InChIKey
  rule, **24** by raw (non-desalted) InChIKey — the figure ADR 0029 originally recorded
  — **36** by PubChem CID, **34** by alphanumerically normalised name. The union of the
  three routes is **42**, a derived number, not a contract figure. ADR 0029 is amended
  to state which route it quotes.
- The obs columns the atlas cannot supply stay present and empty rather than being
  dropped or given invented values, so schema parity never becomes a fabricated fact
  (D3); the two structure-less drugs keep their rows and are flagged, so the corpus's
  no-silent-filtering rule survives the projection.
- Pretrain→finetune domain gaps stay stated, per ADR 0029: 50 cancer cell lines, one
  timepoint, 31/525 drugs. A loss curve does not narrow them, and this run is not
  evidence that pretraining helps.
- Readers are version-sensitive in one direction only: anndata 0.12.16 writes the
  18,164-column `uns` DataFrame (`dataframe/0.2.0`); the pinned benchmark reader
  (anndata 0.10.8) cannot write it but reads it back as the declared frame. The builder
  is pinned to `arc`, the trainer to `arc-stage2-bench-prnet`.

## See

- `docs/adr/0026-stage2-cm-data-contract.md` — the CM corpus and views whose schema this view follows.
- `docs/adr/0029-tahoe-pseudobulk-pretrain-corpus.md` — the corpus this view projects.
- `docs/adr/0027-stage2-benchmark-contract.md` — the recipe mirrored, and the losses this run is not comparable to.
- `stage2/04_tahoe_data/execution/build_pretrain_view.py` — the view builder.
- `stage2/04_tahoe_data/output/pretrain_view/{manifest.csv,compounds.csv,qa.json,qa_report.md}` — the view's tracked record (the h5ad is local-only and gitignored).
- `stage2/04_tahoe_data/README.md` — the step's human record: contract, gates, D0–D3.
- `stage2/05_prnet_pretrain/execution/train_prnet.py` — the trainer.
- `stage2/05_prnet_pretrain/output/{config.json,fit_history.json,summary_table.csv,summary.png,model_card.md}` — the run's tracked record (weights and raw log are local-only).
- `stage2/05_prnet_pretrain/README.md` — the step's human record: recipe, status, reproduce commands.
- `.scratch/stage2-05-prnet-pretrain/spec.md` — the contract and out-of-scope as planned.
- `.scratch/stage2-05-prnet-pretrain/issues/` — the tickets that produced the measurements.
