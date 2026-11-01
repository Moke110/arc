# Pretrain view — QA report

QA run at 2026-11-01T17:40:14+0800 by `stage2/04_tahoe_data/execution/build_pretrain_view.py`.
Contract: `.scratch/stage2-05-prnet-pretrain/spec.md` §Implementation Decisions
("the view's contract"). Source corpus:
`stage2/04_tahoe_data/output/pseudobulk/tahoe_pseudobulk.h5ad` (never rebuilt; the view only projects the frozen rows onto the frozen
18,164-gene CM axis).

## Gates (nine)

| # | gate | result |
|---|---|---|
| 1 | shape, controls, perturbed, control groups, row identity | PASS |
| 2 | var identical to corpus/v1/genes.csv (values, order, dtypes) | PASS |
| 3 | lib_size == X.sum(1); n_genes_detected == (X>0).sum(1) | PASS |
| 4 | log1p_cp10k recomputes from X (atol=1e-5, 200 rows) | PASS |
| 5 | X non-negative int32 counts (subset spot-checked separately) | PASS |
| 5.1 | X is a pure subset of the corpus counts (3 random rows) | PASS |
| 6 | control_mean covers 700 groups, matches recomputed means, round-trips | PASS |
| 7 | no forbidden obs columns; structure + D3 emptiness policy exact | PASS |
| 8 | paired_control_index resolves; shim dose is 0.0 exactly on controls | PASS |

Nine gates: **9/9 PASS** (`all_pass=True`).
`--qa` exits non-zero when any gate fails.
The machine-readable measurements are in [`qa.json`](qa.json).

## Deviations from the CM contract (D0–D3)

| id | topic | this view | why |
|---|---|---|---|
| D0 | matrix layout | CM views are dense int32; this view is CSR int32 (+ CSR float32 layer) | mirrors the corpus's per-plate CSR parts and the slice-write assembly; the plan's 'dense = 4.87 TB' was a units error — dense is 4.87 GB (9.74 GB with the layer) and the matrix is 83% dense, so this 16.32 GB CSR artifact is larger than a dense one |
| D1 | control keying | `control_group` = `tahoe::<plate>::<cell_line_id>` (700 groups), `uns['control_mean']` has 700 rows | a plate pools ~50 cell lines; a plate-level mean is the Simpson's-paradox mixture ADR 0029's "Rejected (a)" rejects, so Δ stays inside one cell line |
| D2 | gene axis | 18,164 CM slots; `FAM153B` / ENSG00000182230 is an all-zero row | the gene has no Tahoe measurement; dropping it would break axis identity |
| D3 | provenance gaps | `time_h` NaN; `cell_type`/`donor` NaN; `atc`/`moa`/`mesh_classes`/`annotation_source` ""; `pharm_atc_moa`/`pharm_mesh` False; and the derived or constant columns `dataset_id` "tahoe", `dose_source` "nominal_atlas", `vehicle` `DMSO_TF` on controls / "" on treated, `inchikey` derived from the SMILES rather than looked up | the atlas carries no time column and no CM annotation vocabulary, and it has no InChIKey; recorded rather than invented |

## Gene mapping (verbatim CM axis, joined on `ensembl_id`)

- CM axis: **18164** slots, order and values copied from `corpus/v1/genes.csv`.
- Measured by Ensembl: **18163/18164**; by symbol:
  **17935/18164** — the symbol route loses
  229 of the 18164 CM genes
  (228 relative to the 18163 that map
  by Ensembl), so the join is on Ensembl only.
- Unmeasured: `ENSG00000182230` (`FAM153B`) — kept as an
  all-zero row so the axis stays byte-identical to every CM view (D2).

## Null-structure drugs

2 of 379 Tahoe drugs have no `canonical_smiles` (the other
377 resolve to a structure):
`Sacubitril/Valsartan`, `Verteporfin`. Their `has_structure` is False and
`drug_key`/`inchikey`/`smiles`/`SMILES`/`pubchem_cid` are `""`; the matching says
nothing about biology.

## Drug overlap with the CM corpus

`compounds.csv` has 379 rows and `overlap_cm` = *the desalted parent
InChIKey is one of the 525 CM `drug_key`s*.
Raw and desalted keys are both stored (`raw_inchikey`, `parent_inchikey`).

| route | distinct CM drug keys (of 525) | definition |
|---|---|---|
| desalted parent InChIKey | **31** | `LargestFragmentChooser` + `MolToInchiKey` on `canonical_smiles`; the corpus's own identity rule (ticket/spec: 31) |
| raw, non-desalted InChIKey | **24** | `MolToInchiKey` without desalting (ADR 0029's "24/525") |
| PubChem CID | **36** | joined on `pubchem_cid` (ticket/spec: 36) |
| alphanumerically normalised name | **34** | lowercase, every non-`[a-z0-9]` character removed; joined on `compound_name` (ticket/spec: 34) |
| **union of the three identity routes** | **42** | derived: parent ∪ CID ∪ name — **not** stated in the spec or ticket |

Counting is on **distinct CM drug keys**, not on Tahoe rows: `33`
Tahoe drug rows match by desalted parent key but
2 of them share a key with another row
(e.g. `Trametinib` / `Trametinib (DMSO_TF solvate)`), which is why `overlap_cm`
has 33 True rows and 31
distinct keys. The recorded routes are reproduced exactly (expected/measured:
desalted 31/31,
raw 24/24,
CID 36/36,
name 34/34).
The union 42 is a *derived* number recorded from those routes; it is not a spec figure.

## Author-compatibility shims (not contract columns)

- `SMILES` = `smiles` (the author `Dataset.py` reads `SMILES`).
- `dose` = `dose_um`, so `0.0` marks a control, matching
  `pert_model/PRnet/source/data/Dataset.py:39`.
- `paired_control_index` = obs label of the first control row of the row's
  `control_group`; on control rows it is the row's own label.

## `uns['control_mean']` encoding (consumer compatibility)

`uns['control_mean']` is an anndata **dataframe/0.2.0**
DataFrame (column-wise storage, `_index = control_group`):
700 control groups × 18164 genes, float64, index named
`control_group`, columns = the 18,164 var symbols in var order. Gate 6 re-reads it through a
plain `anndata.read_h5ad` consumer read, requires the declared 700 x 18,164 shape and
`control_group` index name, and compares it with a group mean recomputed from `X`; the writer
here is anndata 0.12.16. The PRnet benchmark environment pins anndata 0.10.8,
which cannot *write* an 18k-column DataFrame but reads this one: verified — anndata 0.10.8 read it back as a `DataFrame` [700, 18164] with index `control_group` (interpreter `/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python`) (recorded,
not one of the nine gates).

## Shape

`67,018 x 18,164` CSR int32 + `layers['log1p_cp10k']` CSR float32,
`uns['control_mean']` float64 700 x 18164 (D1),
16.32 GB on disk. Content SHA-256 in `manifest.csv` covers the
CSR `data`/`indices`/`indptr` arrays separately (chunked; never densified).
