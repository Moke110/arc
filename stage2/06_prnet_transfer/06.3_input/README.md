# 06.3_input — the frozen two-source 18,164-gene input

**Status:** planned — builder built and tested, **artifacts not frozen yet**
(2026-11-02). The freeze is blocked on the #1 eligibility handoff and the #2
per-observation QC handoff, and the counts below cannot be final until both
land. `output/` holds no frozen evidence. Flip this line to *built and verified
(YYYY-MM-DD)* only when `build_input.py --verify` passes on the real handoffs.

This substep is the only place the formal-transfer study's training input is
derived. It freezes **one** two-source input on the frozen 18,164-gene axis, the
five drug-key folds, the per-fold consensus control `c*` and the locked
experiment config. Nothing downstream may re-derive any of them (ADR 0032).

It does **not** train, select or score anything, and it does not rewrite the
frozen inputs it reads: `02_cm_data/output/` (corpus + `datasets-1` view),
`04_tahoe_data/output/` (pseudobulk corpus + pretrain view) and every `source/`
tree stay byte-identical. It performs no network access and no model run.

## Layout

```
06.3_input/
├── AGENTS.md  README.md
├── execution/
│   ├── build_input.py     # --check | build | --verify | --stub-handoffs
│   ├── contract.py        # every frozen path, constant, schema and the locked config
│   ├── folds.py           # five drug-key folds (pure functions)
│   ├── consensus.py       # per-fold consensus c* and its membership
│   ├── writers.py         # observation assembly + the two h5ad writers
│   ├── Makefile  run.sh   # `make test`, `make check`, `make build`, `make verify`
│   └── tests/             # stdlib unittest: contract, folds, consensus, h5ad, manifest, config
└── output/                # frozen artifacts (h5ad local-only, tracked record force-added)
```

## The frozen contract

### 1. Gene axis

18,164 HGNC symbols in `stage2/02_cm_data/output/corpus/v1/genes.csv` order. Both
sources already live on this axis; the build asserts symbol-order equality
against the corpus for each artifact and records a per-source check
(`var_checks.csv`). `var` carries `hgnc_id`, `entrez_id`, `ensembl_id`,
`match_type`, `gene_measured` and `n_nonzero` (index = `symbol`).

`gene_measured` is **per source** and derived from membership, never from
non-zero counts:

| artifact | derivation | result |
|---|---|---|
| `cm_transfer.h5ad` | symbol ∈ the frozen `temp_env` resolved gene set (`corpus/v1/source_gene_sets.csv`) | all 18,164 measured |
| `tahoe_transfer.h5ad` | `ensembl_id` ∈ the Tahoe corpus's 62,710-gene space (`pseudobulk/tahoe_pseudobulk.h5ad` var) | 18,163 measured, **`FAM153B` / `ENSG00000182230` unmeasured** |

Measured fact worth keeping: **18 Tahoe columns are all-zero although measured**
(`NKX2-1`, `KRTAP10-1`, …) and two CM columns are all-zero although measured
(`CDCA2`, `TEX38`), so an all-zero slot is not evidence of a missing gene. Only
`FAM153B` is absent from the atlas, and it must never be read as a measured zero
(ADR 0031 D2).

### 2. CM artifact — `cm_transfer.h5ad`

Built from `stage2/02_cm_data/output/datasets-1/temp_env_all-smiles.h5ad`
(4,478 × 18,164, dense `int32` `X` + dense `float32` `layers['log1p_cp10k']`,
24 obs columns). Row rule: keep the QC-pass treated wells whose `drug_key` is
resolvable **plus** the QC-pass control wells needed for matched controls and the
consensus; rows failing #2's QC are dropped and counted. `X` stays raw `int32`,
`layers['log1p_cp10k']` stays the row-wise `log1p(CP10k)` `float32` — neither
values nor gene axis change.

The 24 inherited columns are unchanged; these are added:

| column | meaning |
|---|---|
| `observation_id` | the CM well id (also the obs index `well`) |
| `physical_well` | the raw well token (`B2`), the join key of the design table |
| `plate_no` | 1–16 |
| `replicate_raw` | the raw design `Replicate_Number`, re-joined by `(plate, normalised well)` |
| `solvent` | raw design `Solvent` (`DMSO` / `Media` / `Blank`) |
| `chemical_index` | raw design `Chemical_Index` (the source's own chemical id) |
| `eligibility_decision` | `include` / `exclude` / `uncertain` from #1 (NaN on controls) |
| `elig_class` | #1's eligibility class |
| `qc_pass`, `qc_reason` | #2's per-observation QC |
| `main_set` | treated **and** `eligibility_decision == include` **and** `qc_pass` |

`obs.replicate` is a **generated** within-group index (`cmdata/views.py:50`), not
the raw design; the raw value is `replicate_raw`, and both the `uns`
`column_provenance` block and this table say so. `obs.time_h = 24.0` on every row
is a loader constant (`cmdata/sources.py:539`) whose 24 h design is recorded in
GSE262419, not a raw metadata field.

### 3. Tahoe artifact — `tahoe_transfer.h5ad`

Built from `stage2/04_tahoe_data/output/pretrain_view/tahoe_cm_schema.h5ad`
(67,018 × 18,164 CSR `int32` + CSR `float32` log layer). `X` and
`layers['log1p_cp10k']` are copied through h5py without ever densifying the
16 GB matrix; the 35 inherited obs columns are kept. Added:

| column | meaning |
|---|---|
| `observation_id` | the Tahoe observation id (`tahoe::plate1:smp_1495:CVCL_0399`) |
| `physical_well` | the `smp_*` token — `obs.well` is composite, the physical well is not unique alone |
| `tahoe_split` | `train` \| `valid` by `plate13`/`plate14` (the existing technical default) |
| `tahoe_valid_group` | validation grouping key: `drug_key × dose_µM × cell_line_id` (controls: `control × plate × cell_line`) |
| `qc_pass`, `qc_reason` | #2's per-observation QC |
| `tahoe_usable` | `~is_control AND has_structure AND qc_pass AND` its control group keeps ≥ 1 QC-pass control |
| `time_source` | `literature:10.1016/j.cell.2026.08.035` on every row |

`time_h` stays **NaN**: no atlas table carries a time column (ADR 0029 update
2026-11-01, ADR 0031 D3), so the 24 h design is cited as literature provenance
and never written as a constant and never used to filter rows. No row is dropped
from Tahoe: the 1,400 controls, the QC-failed observations and the 300 treated
rows of the two structure-less drugs stay in the artifact for auditability and
are excluded from use through `tahoe_usable`.

Matched controls are **not** materialised as 67k × 18k vectors. Two side files
carry the membership: `tahoe_control_groups.parquet` (700 `(plate, cell_line)`
groups → their control row positions, `n_controls`, `n_qc_pass`) and
`tahoe_pairs.npz` (treated row → `control_group`, `usable`, `split`). The
pipeline computes the matched control vector as the mean of that group's
`log1p_cp10k` rows — never across cell lines, never a 50-line average.

### 4. Consensus control `c*` — `consensus_controls.npz` + `consensus_membership.parquet`

Per fold, `c*` is one 18,164 `float32` vector built only from **control** rows:

1. candidate = CM control wells with `vehicle == 'DMSO'` **and** `qc_pass` **on
   plates that carry at least one training drug of that fold** — a plate whose
   fold drugs are all test/valid is a held-out environment and contributes
   nothing (`Media` is not DMSO, `Blank` is not a control, and neither is used);
2. per plate, the arithmetic mean of `log1p_cp10k` over that plate's candidates;
3. equal weight across plates: `plate_weight = 1 / n_plates`, and
   `well_weight = plate_weight / n_wells_in_plate`, so
   `c* = Σ well_weight · log1p(row)` and the weights of a fold sum to exactly 1.

`consensus_membership.parquet` lists every contributing row
(`fold, plate, observation_id, row_ix, weight, plate_weight, n_wells_plate,
qc_pass, vehicle`), so the pool is auditable and reproducible;
`consensus_controls.npz` holds `c_star` (5 × 18,164), `plate_names`,
`plate_weights`, `n_plates`, `n_wells`.

**Leak assertions (fail closed, re-run by `--verify`):** every member row is a
control row; every member row is DMSO and QC-pass; every member row's plate
carries a training drug of that fold. Treated expression of a test/valid drug
therefore cannot be in `c*`, and a held-out plate's controls cannot be either.

**Deployment semantics.** Training sees matched same-plate DMSO controls.
Validation, test and inference receive this single fixed `c*` for their fold.
The experimental truth Δ = `y − mean(log1p_cp10k of matched same-plate DMSO
rows)` is computed **only at scoring time** and never enters a model input.

### 5. Folds — `folds.parquet`

One row per `(drug_key, fold, role, plate)` for the `main_set` drugs only
(`role ∈ train|valid|test`). Five folds by `drug_key`, seed **20261102**:
clusters are connected components of `drug_key ↔ inchikey[:14]` (so salt/alias
pairs and shared skeletons cannot straddle a boundary), placed greedily with a
capacity of `ceil(n_main_drugs / 5)` drugs per fold, and
`valid_fold = (test_fold + 4) % 5`. Every main-set drug is test in exactly one
fold and valid in exactly one; every dose and replicate of a drug moves with it.

`exclude` and `uncertain` drugs get **no fold row** — they are reserved for the
#11 expansion ablation and are never test or valid. `fold_sizes.csv` reports the
per-`(fold, role)` drug/well counts.

**Drug vs plate.** 459 of the 460 temp_env keys sit on a single plate (one key
spans two), so a fold's test drugs bring their plates with them: a **drug split
is not a plate split**, and a validation/test plate may share the culture
background of training plates. That is why the consensus is scoped by
"plates that carry a training drug", not by "plates without a test drug".

### 6. Locked config — `locked_config.json`

Versioned and hashed (`1.0.0`, digest recorded in `manifest.json`). It fixes the
model structure (hidden 128, drug latent 64, state latent 64, noise 10, in_dim
18,164, adaptor 128, FCFP4 1024, dropout 0.05, `comb_num` 1), the dose encoding
(FCFP4 radius 2 / `useFeatures=True` / 1024 bits × `log10(dose_µM + 1)`, no
binning), the Gaussian NLL loss (never clamped), Adam / batch 512 / weight decay
1e-8 / 100-epoch cap / min 20 epochs / early-stop patience 20, the plateau
scheduler (factor 0.5, patience 10, **max** semantics on validation direction
AUPRC@200), pretrain LR 1e-3, CM LR candidates {1e-4, 3e-4, 1e-3}, the three
seeds, the per-init per-fold common-LR rule (mean of the three seeds' validation
AUPRC@200), the direction-AUPRC@200 definition (k = 200; per-condition AP up and
down, macro over conditions; within-condition replicate aggregation first;
Tahoe conditions = drug × dose × cell line), the deterministic prediction policy
(noise = 0 at validation/test/inference) and the 30-test-unit policy. Entries the
approved spec does not state are marked `derived_default` with a one-line
rationale; no value contradicts the spec.

### 7. Verification — `manifest.json` + `build_input.py --verify`

`manifest.json` records sha256/bytes of every artifact, the frozen source
records (content hashes from the 02/04 manifests, `genes.csv` file hash, the
GSE262419 design-table hash, the derived Tahoe gene-space digest), the #1/#2
handoff hashes and row counts, all counts, the locked-config hash, the git SHA
and the runtime environment. `--verify` is offline and re-computes every hash,
then asserts the gene axis, the fold invariants, the consensus provenance and
recomputation, the count reconciliation against the handoffs, the absence of
forbidden NaN (Tahoe `time_h` is the one declared NaN), the control membership
and that no frozen source changed. A stub build is labelled `status: stub` and
rejected unless `--allow-stub`.

## Artifacts and who reads them

| file | git | #4 reads |
|---|---|---|
| `cm_transfer.h5ad` (dense, ~242 MB) | local-only | training/validation/test rows, flags, `main_set` |
| `tahoe_transfer.h5ad` (CSR, ~16.3 GB) | local-only | pretrain rows, `tahoe_usable`, `tahoe_split`, identity columns |
| `folds.parquet`, `fold_sizes.csv` | tracked (force-added) | fold/role assignment |
| `consensus_controls.npz`, `consensus_membership.parquet`, `consensus_summary.csv` | tracked | fixed `c*` per fold + audit |
| `tahoe_control_groups.parquet`, `tahoe_pairs.npz` | tracked | matched-control mean at load time |
| `locked_config.json`, `manifest.json`, `counts.csv`, `var_checks.csv` | tracked | recipe, identity, counts, measurement mask |

Neither h5ad carries `uns['control_mean']` (see D0). Raw counts are in `X`, so
any raw-count mean stays reproducible; the correctly scaled caches are the
explicit ones above.

## Counts

Structural counts measured from the frozen sources (independent of #1/#2):

| fact | value | source |
|---|---|---|
| gene axis | 18,164 symbols | `corpus/v1/genes.csv` |
| CM view rows | 4,478 (4,206 treated / 272 control) | `datasets-1/manifest.csv` |
| CM drug keys | 460 | `datasets-1/manifest.csv` |
| CM controls per plate | 12 DMSO + 5 Media | measured |
| Tahoe rows | 67,018 (65,618 treated / 1,400 control) | `pretrain_view/manifest.csv` |
| Tahoe control groups | 700 × 2 control rows | measured |
| Tahoe drugs | 379 names / 355 `drug_key`, 2 structure-less | `pretrain_view/manifest.csv` |
| Tahoe cell lines / plates | 50 / 14 | measured |
| unmeasured genes | 1 (`FAM153B`, Tahoe only) | derived, this build |

The tier-1 counts that depend on the handoffs — main-set drugs/wells, folds per
role, `c*` plate/well pools, Tahoe usable/train/valid — are **pending #1/#2**
and are written to `counts.csv` by the build. No number here is a placeholder
presented as a result.

## Reproduce

Run from the repository root. The builder must use the `arc` environment
(anndata 0.12.16); the pinned trainer reads the result with anndata 0.10.8.

```bash
ARC=/home/chang/miniconda3/envs/arc/bin/python
B=stage2/06_prnet_transfer/06.3_input/execution/build_input.py

# 1. validate frozen sources + handoffs (writes nothing)
$ARC $B --check

# 2. build and freeze into the main checkout's artifact tree
$ARC $B --out-root /home/chang/arc/stage2/06_prnet_transfer

# 3. offline re-verification (hashes + every invariant)
$ARC $B --verify --out-root /home/chang/arc/stage2/06_prnet_transfer

# 4. tests (offline; stub build + verify end-to-end)
make -C stage2/06_prnet_transfer/06.3_input/execution test
```

Inside an implementation worktree (tracked files only) prefix the commands with
`ARC_DATA_ROOT=/home/chang/arc` so the frozen 16 GB sources are read from the
main checkout. `run.sh` resolves the repository root by walking up to `.git` and
execs the same CLI.

## Upstream handoffs (#1 / #2)

The build consumes two handoffs and **fails closed** when either is missing,
partial or malformed. Paths are searched under
`stage2/06_prnet_transfer/{06.1_eligibility,06.2_qc}/output/` (overridable with
`--eligibility`, `--cm-qc`, `--tahoe-qc`, or a single combined `--qc` table with
a `source` column):

| handoff | required columns (aliases accepted) | coverage required |
|---|---|---|
| #1 eligibility | `drug_key`, `eligibility_decision` ∈ {include, exclude, uncertain}, `elig_class` | every resolvable CM `drug_key` |
| #2 CM QC | `well`, `qc_pass` (boolean-like), `qc_reason` | every CM view row |
| #2 Tahoe QC | `well`, `qc_pass`, `qc_reason` | every Tahoe observation |

## Environments

| role | interpreter |
|---|---|
| build / write / tests | `/home/chang/miniconda3/envs/arc/bin/python` (Python 3.13, anndata 0.12.16) |
| read / train | `/home/chang/miniconda3/envs/arc-stage2-bench-prnet/bin/python` (Python 3.9, anndata 0.10.8, torch 2.0.1+cu118) |

The writer/reader split is one-directional: anndata 0.12.16 writes the
`dataframe/0.2.0` obs/var encoding and the pinned 0.10.8 reads it; 0.10.8 cannot
write that encoding. The test suite includes a subprocess check that the pinned
reader opens the CM artifact.

## Deviations (D0–D3)

| id | deviation | why |
|---|---|---|
| D0 | Neither artifact carries `uns['control_mean']`; the frozen sources keep their caches untouched. | The 02 cache is a **DMSO+Media raw-count** mean inside one group and the 04 cache is a raw-count mean against log-scale targets (the ADR 0031 recorded deviation). Carrying either would silently invite a raw/log mismatch. Raw counts stay in `X`, and the correctly scaled controls are explicit: per-fold `c*` with row membership, and Tahoe control-group membership computed on the log layer. |
| D1 | Asymmetric row retention: CM keeps QC-pass rows only; Tahoe keeps **all** 67,018 rows flagged. | The CM artifact is the training input and #2's QC is per observation, so dropped rows are counted and listed by #2; Tahoe keeps the structure-less and QC-failed rows so every exclusion stays auditable in one frozen file (`tahoe_usable`). |
| D2 | `gene_measured` is per source and derived from gene-space membership, not from non-zero counts. | 18 Tahoe columns and 2 CM columns are all-zero yet measured; only `FAM153B` is absent from the Tahoe gene space. Deriving measurement from counts would mislabel 20 measured genes. |
| D3 | Tahoe `time_h` stays NaN with a `time_source` literature citation; CM `time_h = 24.0` is recorded as a loader constant. | No atlas table carries a timepoint (ADR 0029/0031 D3), so writing `24.0` would fabricate metadata; the CM constant is not raw GSE metadata either, so its provenance is recorded rather than implied. |

## What #4 may rely on

- One gene axis, one fold table, one `c*` per fold, one config — read them, do
  not recompute. `folds.parquet` is the only fold source; `consensus_controls.npz`
  is the only inference-time control; `locked_config.json` is the only recipe.
- `main_set` selects CM training/validation/test rows as treated wells.
- Controls: use `qc_pass & is_control & vehicle == 'DMSO'` on the same plate for
  CM training, the fold's `c*` for CM validation/test, and
  `tahoe_control_groups.parquet` + `tahoe_pairs.npz` (mean on the log layer) for
  Tahoe.
- Do not read `obs.control_group` as a matched baseline (it mixes vehicles), do
  not use `replicate` as a raw replicate (use `replicate_raw`), do not treat
  `tahoe_split` as a cold-drug split, and do not use `time_h` to filter Tahoe.

## Out of scope

No training, no model selection, no evaluation, no eligibility or QC decisions
(those are #1/#2), no regeneration of 02/04/05 artifacts, no ScreenSeq or other
CM source, no dose extrapolation claim, and no change to the frozen 18,164-gene
axis.
