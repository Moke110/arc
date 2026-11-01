# 04_tahoe_data — Tahoe-100M pseudobulk pretrain corpus

**Status:** built and verified (2026-11-01). The corpus below is frozen at its
2026-10-30 build; the pretrain view was built on 2026-11-01 and passed all nine
consumer-side gates (9/9), and PRnet was pretrained on it the same day (68
epochs, best valid Gaussian NLL −2.628080; ADR 0031).

The **whole-transcriptome pretraining input** for Stage 2: a bulk, well-level
corpus derived from Tahoe-100M and built once, then read by the pretrain /
adaptation work. CM data is bulk (one physical library per well), so the
pretraining substrate is pseudobulk profiles rather than single cells;
[ADR 0029](../../docs/adr/0029-tahoe-pseudobulk-pretrain-corpus.md) freezes how
they are built. The corpus is model-agnostic: identity (genes, units) and raw
counts are frozen; thresholds, Δ pairing and gene-space projections belong to a
view or to the model.

**This step performs no download.** The atlas is an external item
([`pert_data/Tahoe-100M/`](../../pert_data/Tahoe-100M/) and ADR 0030); its
importer and verifiers live in `pert_data/Tahoe-100M/script/` (`fetch_tahoe.py`,
`tahoe-dl.sh`, `run_fetch.sh`, `verify_import.py`, offline `verify.sh`). This
step only aggregates bytes that are already on disk.

## Layout

```
04_tahoe_data/
├── AGENTS.md  README.md
├── execution/build_pseudobulk.py      the four resumable stages (corpus, below)
├── execution/build_pretrain_view.py   the model-facing view (built, below)
├── output/pseudobulk/
│   ├── tahoe_pseudobulk.h5ad          25.8 GB corpus (local-only)
│   ├── manifest.csv                   67,018 profiles x plate/well/line/drug/dose/n_cells
│   ├── qa.json  qa_report.md          the gates and their measured values
│   ├── qa_contrasts.csv               3,755 per-contrast agreement rows
│   ├── _parts/                        14 per-plate CSR parts (resumable intermediate)
│   └── _state/shard_plate.json        scan bookkeeping (machine state)
└── output/pretrain_view/              built + verified 2026-11-01 (9/9 gates)
    ├── tahoe_cm_schema.h5ad           67,018 x 18,164 CM-schema view (local-only)
    ├── manifest.csv  compounds.csv  qa.json  qa_report.md   (tracked build record)
    └── _view_parts/                   per-plate gather cache (local-only)
```

## Corpus

| | |
|---|---|
| unit | `(plate, well, cell_line)` — well is the upstream `sample`; **67,018 profiles** |
| gene space | Tahoe native **62,710** (`token_id - 3`), not the 18,164 CM intersection |
| `X` | summed raw UMI counts (`int32`) |
| `layers['log1p_cp10k']` | `log1p(1e4 * x / row_sum)` |
| controls | 1,400 `DMSO_TF` profiles, `control_group = tahoe::plate` |
| cells consumed | 95,624,334 — every cell exactly once (matches the atlas card) |
| file | `tahoe_pseudobulk.h5ad`, 25.8 GB |

Δ is deliberately **not** materialized: it is
`log1p_cp10k(treated) − mean(same-plate, same-cell-line DMSO_TF)`, one groupby
away, and freezing it would freeze a pairing rule that belongs to the model.

## Pretrain view — built and verified (2026-11-01)

`output/pretrain_view/tahoe_cm_schema.h5ad` is the **model-facing** projection
ADR 0029 left to a view: the frozen corpus restricted to the CM corpus's
18,164-gene axis and shaped so that a CM-shaped loader reads it with no special
cases. Contract:
[`spec.md`](../../.scratch/stage2-05-prnet-pretrain/spec.md) §Implementation
Decisions ("the view's contract"); the accepted decision record is
[ADR 0031](../../docs/adr/0031-tahoe-cm-schema-view-and-prnet-pretrain.md)
(accepted 2026-11-01).

**Built and verified 2026-11-01** — 16,319,543,984 bytes (**16.32 GB**), nine
consumer-side gates **9/9 PASS** (`qa.json` `all_pass: true`, QA run at
2026-11-01T16:59:42+0800). Full build + QA took **175 s** (gather 90 s, assemble
34 s, QA ~50 s); a warm-cache rebuild is **80 s** (pass 1 skipped, `_view_parts/`
is the resumable gather cache). The build streams the corpus's 14 per-plate
parts — gather the retained columns, then fill the preallocated CSR arrays with
slice writes — so it never materialises a dense buffer (peak RSS ~2 GiB).

The layout is **not** a size win, and the plan's figure was a units error: a
dense `int32` `X` is **4.87 GB** (**9.74 GB** with the `float32` layer) — not
terabytes — and at **83% density** (1,010,892,356 of 1,217,314,952 cells) this
16.32 GB CSR artifact is the *larger* container. D0 keeps CSR because it mirrors
the corpus's per-plate CSR parts and lets the builder assemble with slice
writes — a recorded container deviation, not a size claim.

PRnet was pretrained on the view on 2026-11-01 — **68 epochs** (early stop fired
at 68), best valid Gaussian NLL **−2.628080 at epoch 47**, **883.2 s** on an
RTX 4070 Ti SUPER; artifacts in `stage2/05_prnet_pretrain/output/` (ADR 0031).
That run claims only that the pipeline runs.

### Contract

| element | this view |
|---|---|
| shape | **67,018 × 18,164** (`well` × CM `symbol`) |
| rows | **1,400** `DMSO_TF` controls + **65,618** treated; **379** drugs |
| `X` | CSR **`int32`** counts, nnz **1,010,892,356**, stored max 1,658,817; the corpus's summed raw counts on the retained columns, no renormalisation |
| `layers['log1p_cp10k']` | CSR **`float32`** on the same sparsity pattern: `log1p(1e4 * x / row_sum)` |
| `var` | the frozen **18,164**-symbol CM axis verbatim, order included; index `symbol`; columns `entrez_id`, `ensembl_id`, `match_type` |
| `obs` | **35 columns**: CM's 24 in the CM order, then 8 Tahoe-native (`cell_line_id`, `cell_name`, `n_cells`, `moa_broad`, `moa_fine`, `pubchem_cid`, `has_structure`, `has_matched_control`), then 3 author shims |
| `obs` index | `well` — unique, equal to the `well` column |
| `obs` per-row | `lib_size` = `X.sum(1)` and `n_genes_detected` = `(X>0).sum(1)`, both exact |
| `uns` | one key, `control_mean`: **700 × 18,164** `float64` DataFrame, index `control_group`, encoding `dataframe/0.2.0` |
| forbidden | no `pert_idx`/`cell_idx`/`dose_idx`/`pert_dose_idx`; `obsm`/`varm`/`obsp`/`varp` empty |
| file | `tahoe_cm_schema.h5ad`; `content_sha256` `476812c0…` (full value in `manifest.csv`), hashed over the CSR `data`/`indices`/`indptr` arrays, chunked |

The three shims exist for the PRnet author code, not as contract columns:
`SMILES` = `smiles` (its `Dataset.py` reads `SMILES`); `dose` = `dose_um`, so
`0.0` marks a control (`pert_model/PRnet/source/data/Dataset.py:39`);
`paired_control_index` = the obs label of the first control row of the row's
`control_group`, and the row's own label on control rows.

`obs` column order as written: `well`, `dataset_id`, `compound_name`, `drug_key`,
`inchikey`, `smiles`, `dose_um`, `dose_source`, `time_h`, `is_control`,
`vehicle`, `plate`, `control_group`, `cell_type`, `donor`, `replicate`, `atc`,
`moa`, `mesh_classes`, `pharm_atc_moa`, `pharm_mesh`, `annotation_source`,
`lib_size`, `n_genes_detected`, `cell_line_id`, `cell_name`, `n_cells`,
`moa_broad`, `moa_fine`, `pubchem_cid`, `has_structure`, `has_matched_control`,
`SMILES`, `dose`, `paired_control_index` — 24 CM + 8 Tahoe-native + 3 shims.
`var` column order `entrez_id`, `ensembl_id`, `match_type`; both frames carry
`dataframe/0.2.0` encoding.

`manifest.csv` is the authoritative record: `n_obs` 67,018, `n_perturbed` 65,618,
`n_controls` 1,400, `n_drugs` 379, `n_genes` 18,164, `n_control_groups` 700,
`corpus_version` `v1`, `source_commit` `c72518e`, `content_sha256`, `bytes`.

### Deviations from the CM contract (D0–D3)

| id | deviation | why |
|---|---|---|
| D0 | layout is **CSR**, not CM's dense `int32` — CSR `int32` counts + a CSR `float32` layer | mirrors the corpus's per-plate CSR parts and the slice-write assembly. **Not** a size claim: dense `int32` `X` is 4.87 GB (9.74 GB with the layer), and at 83% density this 16.32 GB CSR artifact is larger than dense; the plan's figure was a units error (terabytes for gigabytes) |
| D1 | `control_group` = `tahoe::<plate>::<cell_line_id>` (**700** groups), not the corpus's `tahoe::<plate>` (14) | a plate pools ~50 cell lines, so a plate-level mean is the Simpson's-paradox mixture ADR 0029 rejected; Δ stays inside one cell line, and all 700 groups carry a vehicle |
| D2 | 18,164 CM slots with `FAM153B` all-zero | the gene has no Tahoe measurement; dropping it would break axis identity |
| D3 | `time_h` NaN; `cell_type`/`donor` NaN; `atc`/`moa`/`mesh_classes`/`annotation_source` `""`; `pharm_atc_moa`/`pharm_mesh` False | the atlas has no time column and no CM annotation vocabulary — recorded here rather than invented |

### Gene mapping (the join is on `ensembl_id` only)

- CM axis: **18,164** slots, values and order copied from `corpus/v1/genes.csv`.
- Measured by Ensembl: **18,163/18,164**. By symbol: **17,935/18,164** — the
  symbol route loses **229** of the 18,164 CM genes (228 relative to the
  18,163 that map by Ensembl), so the join uses Ensembl only.
- Unmeasured: `ENSG00000182230` (`FAM153B`) — all-zero row (D2).

### Nine gates (9/9 PASS, 2026-11-01)

| # | gate | result | key measured |
|---|---|---|---|
| 1 | shape, controls, perturbed, control groups, row identity | PASS | 67,018 × 18,164; 1,400 controls; 65,618 treated; 379 drugs; 700 groups; `well` index unique and equal to the `well` column |
| 2 | `var` identical to `corpus/v1/genes.csv` (values, order, dtypes) | PASS | 18,164 rows, index `symbol`; columns `entrez_id`/`ensembl_id`/`match_type`; identical |
| 3 | `lib_size == X.sum(1)`; `n_genes_detected == (X>0).sum(1)` | PASS | max abs diff **0** / **0** |
| 4 | `log1p_cp10k` recomputes from `X` (atol 1e-5, 200 rows) | PASS | max abs diff **2.37 × 10⁻⁷** |
| 5 | `X` non-negative `int32` counts | PASS | dtype `int32`, min 1, max 1,658,817, nnz 1,010,892,356 |
| 5.1 | `X` is a pure subset of the corpus counts (3 random rows) | PASS | 3 rows compared against the corpus `_parts`, **0 mismatches** |
| 6 | control mean covers 700 groups, matches recomputed means, round-trips | PASS | 700 × 18,164 `float64`, index `control_group`; 0 treated groups without a control; sampled group `tahoe::plate12::CVCL_1055` diff **0.0**; plain `anndata.read_h5ad` round-trip returns the declared frame |
| 7 | no forbidden obs columns; structure columns empty iff control or no-structure | PASS | 0 mismatches on `drug_key`/`inchikey`/`smiles`/`SMILES`/`pubchem_cid`; 1,700 rows expected empty; no forbidden column |
| 8 | `paired_control_index` resolves; shim `dose` is 0.0 exactly on controls | PASS | 0 unresolved; all 1,400 control rows point to themselves; `dose` 0.0 on all controls and > 0 on all treated |

`qa.json` records `gate_count: 9` / `gates_passed: 9`; `qa_report.md` carries the
same record for humans. `--qa` re-runs these gates over the produced artifact
and exits non-zero if any fails, so it is the view's acceptance test.

### Drug overlap against the CM corpus's 525 `drug_key`s

| route | distinct CM keys | definition |
|---|---|---|
| desalted parent InChIKey | **31** | RDKit `LargestFragmentChooser` + `MolToInchiKey` on `canonical_smiles`; the corpus's own identity rule |
| raw, non-desalted InChIKey | **24** | `MolToInchiKey` without desalting — the "24/525" quoted in ADR 0029 |
| PubChem CID | **36** | joined on `pubchem_cid` |
| alphanumerically normalised name | **34** | lowercase, every non-`[a-z0-9]` character removed; joined on `compound_name` |
| union of the three identity routes | **42** | derived: parent ∪ CID ∪ name — **not** a spec figure |

Counting is on **distinct CM drug keys**, not on Tahoe rows: **33** Tahoe rows
carry `overlap_cm = True` on the desalted route, but 2 of them share their
parent key with another row (`Trametinib` / `Trametinib (DMSO_TF solvate)`),
which is why the route yields 31 distinct keys. `compounds.csv` stores
`raw_inchikey`, `parent_inchikey`, `pubchem_cid`, `has_structure` and
`overlap_cm` per drug; the four recorded routes are asserted against the
builder's expectations (31/24/36/34) so a drift is loud rather than silent.

**379 drugs, 2 with no resolvable structure:** `Sacubitril/Valsartan` and
`Verteporfin` keep their rows with `has_structure=False` and empty
`drug_key`/`inchikey`/`smiles`/`SMILES`/`pubchem_cid` — no silent filtering, the
corpus's rule. The other 377 resolve to a structure.

### What it is not

- **Not the corpus.** The corpus is 62,710 native genes with 14 plate-level
  control groups and no projection; the view is a derived projection of it onto
  18,164 genes with 700 cell-line-scoped groups, and it is regenerable from the
  corpus's `_parts/`.
- **Not a `datasets-1` bundle.** Bundles are built from the CM corpus and carry
  per-view membership booleans; this view ships its own `compounds.csv` with
  `overlap_cm` instead of extending the frozen bundle.
- **Not evidence of transfer.** The view has been consumed once — the PRnet
  pretrain bring-up (ADR 0031) — which claims only that the pipeline runs: no
  finetune, no warm start, no cold-drug comparison.

### Interpreter pinning (a fact of the artifact)

`uns['control_mean']` is an 18,164-column `dataframe/0.2.0` frame that only the
newer anndata can **write**. The builder writes with
`/home/chang/miniconda3/envs/arc/bin/python` (anndata 0.12.16); the trainer and
every other consumer read with the pinned benchmark environment
`/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python` (anndata 0.10.8),
which cannot write that encoding but reads it correctly. That read is exercised
by `--qa` as a subprocess against the pinned interpreter and recorded in
`qa.json.interop_anndata_0_10_8` (recorded, not one of the nine gates).

## Reproduce

Step-local program: `python`, with `numpy`, `pyarrow`, `scipy`, `h5py`,
`anndata`, `pandas` available (the build ran on the shared `arc` interpreter).
Run from the repository root; the script resolves the root by walking up to
`.git`, so any working directory works.

```bash
# report the resolved root / Tahoe input / output dir, then exit
python stage2/04_tahoe_data/execution/build_pseudobulk.py --check

# the four stages (resumable, in order)
python stage2/04_tahoe_data/execution/build_pseudobulk.py scan        # shard -> plate map  (~20 s)
python stage2/04_tahoe_data/execution/build_pseudobulk.py aggregate   # 14 per-plate parts  (~50 min)
python stage2/04_tahoe_data/execution/build_pseudobulk.py assemble    # sparse h5ad + manifest  (~4 min)
python stage2/04_tahoe_data/execution/build_pseudobulk.py qa          # vs the authors' DESeq2 DE
```

`--plates plate4,plate5` limits any stage; an existing plate part is skipped
(`--force` rebuilds it); `--shards N` limits scan/aggregate to the first N
shards for a smoke run. `qa` exits non-zero if any gate fails, so it doubles as
the corpus acceptance test.

### Pretrain view — produce and validate are separate commands

```bash
ARC=/home/chang/miniconda3/envs/arc/bin/python              # writer, anndata 0.12.16
CONSUMER=/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python   # reader, anndata 0.10.8

# report the resolved corpus / corpus parts / CM gene axis / output, then exit
$ARC stage2/04_tahoe_data/execution/build_pretrain_view.py --check

# produce: gather _view_parts/, assemble the h5ad + side files, then run the gates
$ARC stage2/04_tahoe_data/execution/build_pretrain_view.py

# validate only: re-run the nine gates over the existing artifact and rewrite
# qa.json + qa_report.md; exits non-zero on any failed gate
$ARC stage2/04_tahoe_data/execution/build_pretrain_view.py --qa

# re-gather _view_parts/ from the corpus parts (a missing part is rebuilt anyway)
$ARC stage2/04_tahoe_data/execution/build_pretrain_view.py --force

# consumer-side read with the pinned trainer env: prints (700, 18164) control_group
$CONSUMER -c "import anndata as ad; a = ad.read_h5ad('stage2/04_tahoe_data/output/pretrain_view/tahoe_cm_schema.h5ad', backed='r'); m = a.uns['control_mean']; print(m.shape, m.index.name)"
```

The writer/reader split is not interchangeable: the h5ad's `uns['control_mean']`
frame can only be written by anndata 0.12.16, and the trainer/consumer must read
it with the pinned 0.10.8 environment above. `--qa` is the view's acceptance
test; it does not build anything and its exit code is the gate result.

## QA (2026-10-30 build)

From [`output/pseudobulk/qa.json`](output/pseudobulk/qa.json) and
[`qa_report.md`](output/pseudobulk/qa_report.md):

- **67,018 profiles**, of which **61,767 (92.2%) have ≥ 100 cells**.
- median detected genes per profile **25,290** — gate ≥ 12,000: **PASS**.
- **3,755 contrasts** compared against the authors' shipped DESeq2
  `log2FoldChange`:
  - `baseMean > 800` (median 328 genes/contrast): Pearson **0.812** (p10 0.669) — reported, not gated;
  - `baseMean > 500` (644 genes/contrast): Pearson **0.779**, Spearman **0.713** — gate ≥ 0.60: **PASS**;
  - `baseMean > 300`: Pearson **0.741**; all-gene Spearman **0.173** — reported only.
- **700 `(plate, cell_line)` pairs** carry a vehicle control; none missing — **PASS**.

The all-gene Spearman is context, not a gate: on genes with no real effect both
estimates are noise, so rank agreement there is uninformative. The originally
proposed "Spearman ≥ 0.8 over all genes" gate was retired after measurement
(full reasoning in ADR 0029 §QA): even on well-measured genes the achievable
agreement is capped because the shipped DE was produced under a
normalization/aggregation that is not a naive count sum.

## What is tracked, what is local

| path | git |
|---|---|
| `execution/build_pseudobulk.py`, `AGENTS.md`, `README.md` | tracked |
| `execution/build_pretrain_view.py` | tracked |
| `output/pseudobulk/{qa.json,qa_report.md,qa_contrasts.csv,manifest.csv}` | tracked build record (force-added past the `output/` and `*.csv` rules) |
| `output/pseudobulk/tahoe_pseudobulk.h5ad` | local-only, 25.8 GB, ignored |
| `output/pseudobulk/_parts/` | local-only resumable cache, ignored |
| `output/pseudobulk/_state/` | local-only machine state, ignored |
| `output/pretrain_view/{manifest.csv,compounds.csv,qa.json,qa_report.md}` | tracked view record (force-added past the same rules) |
| `output/pretrain_view/tahoe_cm_schema.h5ad` | local-only, 16.32 GB, ignored |
| `output/pretrain_view/_view_parts/` | local-only gather cache, ignored |

`_parts/`, `_state/` and `_view_parts/` are **intermediates, not results** —
they exist only so that a failed or partial rebuild can resume. `_state/` is the
scan's map from expression shard to plate(s); deleting it is safe, `scan`
recreates it. `_view_parts/` is the view's per-plate gather cache: a missing part
is rebuilt automatically, and `--force` re-gathers all of it. All three trees are
hidden from git by the repository-wide `output/` rule, which is why they need
declaring here:

```console
$ git check-ignore -v stage2/04_tahoe_data/output/pseudobulk/_state/shard_plate.json
.gitignore:62:output/	stage2/04_tahoe_data/output/pseudobulk/_state/shard_plate.json
```

For the view the tracked record is exactly `manifest.csv` + `compounds.csv` +
`qa.json` + `qa_report.md`; the h5ad and the gather cache are local-only, so a
clone rebuilds the artifact from the corpus parts rather than checking it out.

## Two traps this build encodes

1. **A well is not a bulk profile.** Every well pools most of the 50 cell lines,
   so summing a well gives a mixture whose DE is a Simpson's-paradox artefact.
   The corpus is cell-line-stratified; `(plate, well, cell_line)` is the unit.
2. **A shard is not always one plate.** 13 of the 3,388 shards straddle a plate
   boundary (e.g. shard 236 holds plate4 + plate5 cells), so aggregation routes
   each *cell* by its own `plate` column and reads the shard once per plate it
   touches.

## Scope

Pseudobulk fixes sparsity, not biology: this is a cancer-cell-line corpus at a
single 24 h timepoint with a 24/525 drug-key overlap against the CM corpus
(ADR 0029 §Consequences). It is a pretraining substrate, not a Stage 2
evaluation dataset. Atlas provenance, the HF revision and the verification
record are in
[`pert_data/Tahoe-100M/local/Tahoe_summary.md`](../../pert_data/Tahoe-100M/local/Tahoe_summary.md).
