# 02_cm_data — Stage 2 CM perturbation corpus and dataset bundles

**Status:** built (2026-10-21). `corpus/v1/` is frozen and `datasets-1/` holds
the five h5ad views plus bundle metadata. The build is reproducible: a rerun
produces an identical `manifest.csv` and identical content hashes. The design
contract is recorded in ADR 0026; this file documents what was actually built.

## Purpose

Give Stage 2 a **model-agnostic** cardiomyocyte perturbation corpus: frozen gene
and drug universes (`corpus/`) plus load-ready h5ad views derived from them
(`datasets-<n>/`).

It is **not** Stage 1's `Union finetune data`: no 978-landmark subset, no
per-plate quantile alignment to the L1000 pretrain, no `drug_idx < 8981` / HG
filter, no per-sample `pert_idx`/`cell_idx` encoding.

## Inputs (read-only, under `pert_data/`)

| source | perturbed wells | drugs | controls | exposure |
|---|---|---|---|---|
| `dtoxs_primary` | 677 | 47 | 102 | 48 h, literature dose |
| `dtoxs_ipsc` | 296 | 37 | 16 | 48 h, series-matrix dose |
| `screen_ctox` | 2,520 | 42 | 552 | 24 h, 8-point half-log series from a curated top dose |
| `temp_env` | 4,206 | 460 | 272 | 24 h, measured 3-point series |

`gse247541` is out of scope (4 antivirals, 144 h, outside the author time map).

### Dataset ids

The two GEO sources carry a short Stage 2 id; the accession remains the
provenance anchor and stays on the raw folder:

| dataset_id | accession | assay | panel | raw folder |
|---|---|---|---|---|
| `screen_ctox` | GSE244740 | ScreenSeq | 33 cardiotoxicants + 9 non-cardiotoxicants | `pert_data/gse244740` |
| `temp_env` | GSE262419 | TempO-Seq | 464 environmental chemicals | `pert_data/gse262419` |

The mapping is defined once in `cmdata/paths.py` (`SOURCE_ACCESSIONS`) and is
repeated in `corpus/v1/README.md` and `qa_report.md`. The build writes the
dataset id into every view file name, `obs` value and QA gate; only the raw
inputs and the read-only Stage 1 baseline keep the accession (QA translates the
Stage 1 well ids before comparing). Historical Stage 1 artifacts and the
`.scratch/cm-data-v1/` work orders still use the GSE accessions.

## Three layers

| layer | content | where |
|---|---|---|
| **source** | the four raw datasets | `pert_data/` (read-only) |
| **CM perturbation corpus** | frozen gene set + drug set, **no expression matrices** | `corpus/v1/` |
| **dataset bundle** | one build output: h5ad views + metadata | `datasets-1/` |

A view = corpus + one filter rule. `datasets-2` and later bundles reuse the same
corpus; the build hard-fails if a source introduces a gene or drug outside it
(see *Drift*).

## Corpus (`corpus/v1/`)

```
corpus/v1/
├── genes.csv               # fixed-order HGNC symbols + entrez/ensembl + match provenance
├── compounds.csv           # every SMILES-resolvable small molecule, all pharmacology tiers as boolean columns
├── excluded_compounds.csv  # perturbagens/technical blanks that never enter the corpus, with reasons
├── source_gene_sets.csv    # frozen per-source resolved genes (drift baseline)
└── README.md               # version, build date, source commit, intersection rule
```

- **Gene set** = intersection of all four sources, named by the **latest HGNC
  symbol**, matched alias/previous-symbol-aware and Entrez-verified where the
  source carries an Entrez ID. Unresolvable or ambiguous genes are dropped and
  counted, never guessed. The frozen set is **18,164 genes** (the TempO-Seq
  panel gates the universe); the authoritative file is `genes.csv`.
- **Drug set** = all small molecules with a resolvable structure (**525**
  `drug_key`s across the four sources; 532 key/name rows), **unfiltered by
  pharmacology**; the tiers are boolean columns, not pre-applied filters:
  `atc-moa` 177, `atc-moa-mesh` 316, `all-smiles` 525.

## Bundle `datasets-1/`

```
datasets-1/
├── dtoxs_cm.h5ad                  # primary + iPSC merged; cell_type kept as a covariate
├── screen_ctox.h5ad
├── temp_env_atc-moa.h5ad         # pharmacology tier: ATC ∪ ChEMBL MoA
├── temp_env_atc-moa-mesh.h5ad    # ATC ∪ MoA ∪ PubChem MeSH
├── temp_env_all-smiles.h5ad      # every SMILES-resolvable small molecule
├── compounds.csv                  # bundle-local copy of the corpus drug table + view membership
├── gene_crosswalk.csv             # source gene id -> HGNC symbol for this bundle
├── manifest.csv                   # per-view counts, corpus version, source commit
└── qa_report.md                   # exclusions, gene losses, tier counts, gate results
```

| view | rule | rows (pert + control) | drugs |
|---|---|---|---|
| `dtoxs_cm.h5ad` | `atc-moa` + cited manual supplement | 957 + 118 = 1,075 | 46 |
| `screen_ctox.h5ad` | `atc-moa` | 2,520 + 552 = 3,072 | 42 |
| `temp_env_atc-moa.h5ad` | `atc-moa` | 1,047 + 272 = 1,319 | 113 |
| `temp_env_atc-moa-mesh.h5ad` | `atc-moa-mesh` | 2,310 + 272 = 2,582 | 252 |
| `temp_env_all-smiles.h5ad` | `all-smiles` | 4,206 + 272 = 4,478 | 460 |

The `dtoxs` counts depend on three code-only DToxS states that resolve through the
code table (`DIC` → diclofenac, `SUN+` → sunitinib, same key as `Sunitinib`)
while the `AMI` state carries a **pentyl-benzofuran homolog structure**
(CID 16218975), not amiodarone. Structure-first identity keeps that arm in
`all-smiles` but not in `atc-moa`. `manifest.csv` is authoritative.

`datasets-1/` is a local data-lake output: the five `h5ad` files and
`gene_crosswalk.csv` are large and gitignored. `manifest.csv`, `compounds.csv`
and `qa_report.md` are committed as the build record; rerunning the build
recreates every file.

The three `temp_env` views are filtered projections of **one** corpus: the
retained rows are value-identical across files, and they share the same 272
control rows, so any downstream difference comes only from the compound filter.

## h5ad contract (identical across all five views)

- `X` = **raw counts** (`int32`). Duplicate probes collapse to the symbol by
  **summing counts** before any normalization.
- `layers["log1p_cp10k"]` = `log1p(counts / lib_size * 1e4)` per sample.
  No per-plate alignment, no batch correction, no HVG selection.
- `var` = HGNC symbols in the corpus order, identical in every file; carries
  `entrez_id`, `ensembl_id`, `match_type`.
- `obs`:

| column | meaning |
|---|---|
| `well` | globally unique well id (`dataset::raw well`), for leakage gating |
| `dataset_id` | `dtoxs_primary` / `dtoxs_ipsc` / `screen_ctox` / `temp_env` |
| `compound_name` | normalized display name |
| `drug_key` | desalted parent InChIKey, stereo retained |
| `inchikey` | full InChIKey |
| `smiles` | canonical SMILES |
| `dose_um` | real µM |
| `dose_source` | `measured` / `series_matrix` / `literature` / `derived_top_dose` |
| `time_h` | 24 / 48 |
| `is_control` | control well |
| `vehicle` | DMSO / EtOH / VEH / Media / NoTreatment |
| `plate`, `control_group` | physical plate id; `dataset::plate`. `screen_ctox` qualifies the plate with the assay (`optimization_plate1` … `screening_plate3`) because the GEO series numbers the two plate types independently |
| `cell_type` | `primary_cm` / `ipsc_cm` |
| `donor` | recorded donor/line, else NA |
| `replicate` | within-group replicate index |
| `atc`, `moa`, `mesh_classes` | pharmacology annotations |
| `pharm_atc_moa`, `pharm_mesh` | auditable tier membership |
| `annotation_source` | `chembl` / `pubchem` / `manual:<citation>` |
| `lib_size`, `n_genes_detected` | count-matrix quality |

- `uns["control_mean"]`: per-`control_group` mean raw-count profile
  (`n_groups × n_genes`), so mean-based models need not recompute it.
- No `pert_idx`, `cell_idx`, `dose_idx`, or any integer encoding: those belong
  to model adapters.

**Physical plates.** `screen_ctox` (GSE244740) numbers Optimization and
Screening plates independently, and the series matrix uses the same number for
both. The loader therefore puts the assay in the plate id
(`optimization_plate1` … `optimization_plate5`, `screening_plate1` …
`screening_plate3`); `control_group = dataset::plate` is one physical plate, not
a merge of two. This is why `screen_ctox.h5ad` has 8 control groups, not 5.

## Identity rules

- **Gene**: latest HGNC symbol; alias- and previous-symbol-aware; Entrez-verified
  where available; ambiguous matches dropped and reported.
- **Drug**: `drug_key` = desalted parent InChIKey **with stereochemistry retained**.
  Connectivity-only keys must not be used: doxorubicin and epirubicin are
  4′-epimers and share the same InChIKey connectivity block. Salt/free-base
  pairs stay separate rows (different wells and doses) but share one identity,
  so `Tamoxifen` / `Tamoxifen citrate`, `Diltiazem` / `Diltiazem hydrochloride`,
  `Estradiol` / `17beta-Estradiol` cannot leak across a cold-drug split.
- **Control**: `control_group = dataset::plate` (the ADR 0016 per-plate scope;
  one physical plate — see *Physical plates*), with `vehicle` exposed so an
  adapter may use same-vehicle controls instead.

## Pharmacology tiers

| tier | rule |
|---|---|
| `atc-moa` | compound has an ATC code **or** a ChEMBL mechanism |
| `atc-moa-mesh` | `atc-moa` **or** a PubChem MeSH class |
| `all-smiles` | any small molecule with a resolvable SMILES |

`atc-moa` is the drug-like tier. MeSH admits herbicides, insecticides,
carcinogens and disinfectants, which is why it is a separate, explicitly looser
tier rather than the default.

**Manual supplement** applies only to the two curated drug panels
(`dtoxs`, `screen_ctox`): a compound that demonstrably has a mechanism but is
missing from ChEMBL's mechanism table may be annotated with a citation
(`annotation_source=manual:<ref>`). Current entry: **SB-431542** (TGF-βR1
inhibitor, 16 wells; Inman et al. 2002, *Mol Pharmacol* 62:65-74). The three
`temp_env` tiers stay strictly automatic — that boundary *is* the comparison.

## Pinned annotation assets

The build makes **no network call**. Two read-only assets are pinned under
`pert_data/annotation/` by `fetch_annotation.py`:

| asset | content | provenance |
|---|---|---|
| `hgnc_complete_set.txt` | 45,277 HGNC rows: symbol, aliases, previous symbols, Entrez, Ensembl | HGNC complete set, URL/date/checksum in `annotation/README.md` |
| `cm_compound_annotation.csv` | 547 `(name, source-structure)` rows: parent/exact InChIKey, SMILES, PubChem CID, ChEMBL id, ATC, MoA, MeSH | one-time public-registry snapshot derived from the Stage 1 external audit (PubChem/UniChem/ChEMBL/RxNav) |

`fetch_annotation.py --hgnc` downloads and pins HGNC; `--compounds
--from-audit <csv>` snapshots the registry table; `--refresh` re-queries the
registries when the upstream audit is unavailable. Gene and drug identity is
resolved from these assets; the build never calls an API.

## Build

```bash
python stage2/02_cm_data/build_cm_datasets.py --bundle datasets-1
```

Single entry point: sources → in-memory corpus → filtered views. Deterministic
(fixed var order, stable row order, frozen corpus). A rerun writes an identical
`manifest.csv` (same content hashes and file sizes). The Stage 2 package:

```
02_cm_data/
├── build_cm_datasets.py     # CLI: --check, --corpus-only, --rebuild-corpus
├── fetch_annotation.py      # pins HGNC + compound registry snapshot
├── cmdata/
│   ├── sources.py           # vendored Stage 1 loaders + Stage 2 fixes
│   ├── genes.py             # HGNC identity, crosswalk, intersection, drift
│   ├── drugs.py             # InChIKey identity, tiers, exclusions
│   ├── views.py             # one writer + schema gate for all five views
│   ├── corpus.py            # freeze/verify corpus/v1
│   ├── qa.py                # cross-view gates, Stage 1 regression, report
│   └── primary_doses.py     # DToxS literature concentrations (copied)
├── tests/test_cmdata.py     # unit + built-bundle contract tests
├── corpus/v1/               # frozen genes + compounds (no matrices)
└── datasets-1/              # the five h5ad views + bundle metadata
```

Useful modes:

```bash
python stage2/02_cm_data/build_cm_datasets.py --check        # source counts
python stage2/02_cm_data/build_cm_datasets.py --corpus-only  # freeze/verify
python stage2/02_cm_data/build_cm_datasets.py --rebuild-corpus  # deliberate v2 edit
python -m unittest discover -s stage2/02_cm_data/tests -v
```

The loaders are a **copy** of
`stage1/xpert-fda-miri-v3/03_cm_finetune_data/_datasets.py` (commit `c921f7d`)
with a provenance header and the name/dose/code fixes; Stage 1 is read-only
reference, never imported.

### Stage 2 loader fixes carried by the copy

- Underscore and trailing-space names normalize before any lookup
  (`Acetylsalicylic_acid`, `5_Fluorouracil`, `Mitomycin_C`, `Amphotericin_B`,
  `3_Azido_3_Deoxythymidine`, `Acetaminophen `, `Citric acid `).
- `screen_ctox` top doses are looked up **after** normalization, so the four
  audited underscore drugs get their curated Table-1 top instead of silently
  falling back to the sunitinib top. The Stage 1 regression reports these 48
  corrected wells.
- `DIC` → diclofenac, 5 µM (`dose_source=measured`, DToxS DCIC workbook);
  `SUN+` shares `SUN`'s structure; `AMI` keeps its own ZINC structure.
- `temp_env` probes are keyed by **symbol**, not by the numeric suffix (a probe
  id, not an Entrez id). Stage 1 preferred the suffix and zeroed or hijacked
  probes; Stage 2 recovers those genes. The regression compares on the
  unambiguous subset and reports the 290 entrez-mapped and 35 suffix-collision
  genes it cannot compare landmark-by-landmark.

### Drift

If a rerun finds a gene or drug outside `corpus/v1/`, the build **fails** and
lists the diff. Promoting a new corpus version (`v2`) is a deliberate edit, not
a silent rebuild.

## QA gates

1. `var` identical across all five views; no duplicate or empty symbols.
2. `well` unique; perturbed rows have `drug_key`/`smiles`/`dose_um > 0`; control
   rows have `dose_um == 0`; every `control_group` has ≥1 control in-file.
3. Counts non-negative, no NaN; `lib_size > 0`, `n_genes_detected > 0`.
4. Tier nesting `atc-moa ⊆ atc-moa-mesh ⊆ all-smiles`; retained rows are
   value-identical; the 272 control rows are identical across the three views.
5. Counts in `manifest.csv` match the files.
6. **Stage 1 regression check**: for each source, the well set and sampled
   counts are compared with
   `stage1/xpert-fda-miri-v3/03_cm_finetune_data/output/*.h5ad` (read-only; no
   Stage 1 code imported). This catches loader drift after the copy. The gate
   also reports the deliberate Stage 2 divergences: 48 `screen_ctox` dose
   corrections for underscore names, and the `temp_env` landmarks whose suffix
   is a probe id (290 Stage 1 entrez-mapped + 35 suffix-collision genes), where
   Stage 2's symbol-first mapping recovers the counts.

`qa_report.md` records every corpus-excluded compound with its reason, every
view-level tier exclusion, gene-side losses (unmapped/ambiguous/alias/entrez),
the intersection size, salt collapses, and the final per-view counts.

### Build results (2026-10-21)

- All gates pass: five-view `var` identity, tier nesting and value identity,
  272 shared `temp_env` controls, per-source well/dose/count regression.
- `corpus/v1`: 18,164 genes; 525 drug keys; tiers 177 / 316 / 525.
- `datasets-1`: 1,075 / 3,072 / 1,319 / 2,582 / 4,478 rows and 46 / 42 / 113 /
  252 / 460 drug keys. `manifest.csv` is authoritative for bytes and content
  hashes; a rerun reproduces it exactly.

## Known exclusions and caveats

- **Non-small-molecule perturbagens**: bevacizumab, cetuximab, rituximab,
  trastuzumab (antibodies); endothelin-1/ET-1, IGF-1, TNF-alpha (peptide/protein).
- **Technical blanks**: 32 `temp_env` `Blank` wells (2 per plate), concentration
  0, distinct from `Media`/`VEH`; neither perturbation nor control.
- **Dose provenance is not uniform**: `screen_ctox` doses are derived from a curated
  top concentration (paper Table 1), `dtoxs_primary` doses are literature/clinical
  Cmax values except `DIC` (5 µM, recorded in the DCIC design workbook),
  `dtoxs_ipsc` comes from the series matrix, and `temp_env` is measured. Always
  read `dose_source` before fitting a dose-response model.
- **Cell type is a covariate**, not a collapsed constant: `dtoxs_cm.h5ad` mixes
  `primary_cm` and `ipsc_cm`.
- **The DToxS `AMI` arm is not amiodarone.** Its metadata row (`ZINC150345845`)
  carries a pentyl-benzofuran homolog structure (CID 16218975); structure-first
  identity keeps it in `all-smiles` only. Do not relabel it by name.
- **`drug_key` must keep stereochemistry.** Connectivity-only keys merge
  doxorubicin and epirubicin (4′-epimers); the build uses full parent InChIKeys.
- **Dose ranges are wide** (nM to mM, e.g. digoxin 0.1 µM top vs acetaminophen
  4,000 µM top in `screen_ctox`); no clipping is applied.
