# 06.3_input — the frozen two-source 18,164-gene input

**Status: built and verified, frozen (2026-11-02).** Ticket #3 of the PRnet
formal-transfer effort. `build_input.py --verify` passes **125/125** checks on
the real handoffs, so the record below is frozen evidence, not a plan.

This substep is the only place the formal-transfer study's training input is
derived. It freezes **one** two-source input on the frozen 18,164-gene axis, the
five drug-key folds, the per-fold consensus control `c*` and the locked
experiment config. Nothing downstream may re-derive any of them (ADR 0032).

It consumes two upstream handoffs — ticket #1's **rule v1.3 / handoff v2**
(resolved through `06.1_eligibility/output/current.json`, never a hard-coded
version directory) and ticket #2's per-observation QC — and it does not rewrite
the frozen inputs it reads: `02_cm_data/output/` (corpus + `datasets-1` view),
`04_tahoe_data/output/` (pseudobulk corpus + pretrain view) and every `source/`
tree stay byte-identical. It performs no network access and no model run.

## Frozen counts

Every number below is written by the build into `output/counts.csv`, recorded in
`output/manifest.json`, and re-checked by `--verify`.

| fact | value |
|---|---|
| gene axis | 18,164 HGNC symbols, corpus order |
| CM artifact | 3,947 rows (3,745 treated / 202 control) × 18,164 |
| CM main set (`main_set`) | **118 drugs / 951 wells** / 16 plates, doses 0.2, 1, 10 µM |
| CM expansion set (`expansion_set`, #11 only) | 314 drugs / 2,553 wells — never scored |
| folds (test drugs) | 23 / 23 / 24 / 24 / 24 |
| folds (train drugs) | 72 / 71 / 70 / 70 / 71 |
| folds (wells) | test 188–192, valid 188–192, train 569–572 per fold |
| consensus `c*` | 5 × 18,164 float32; **15–16 plates, 127–138 wells** per fold; 667 member rows |
| Tahoe artifact | 67,018 rows (65,618 treated / 1,400 control) × 18,164 |
| Tahoe usable | 59,044 treated (50,520 train / 8,524 valid) |
| Tahoe matched controls | 649 control groups, 60,861 pairs, 1,294 QC-pass controls; 4,757 treated rows without a control are **reported, never imputed** |
| unmeasured genes | 1 (`FAM153B`, Tahoe only) |
| artifact sizes | `cm_transfer.h5ad` 218,559,498 B; `tahoe_transfer.h5ad` 16,207,298,308 B |

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

Measured fact worth keeping: **17 Tahoe columns are all-zero although measured**
(`NKX2-1`, `KRTAP10-1`, `SFTA3`, …) and two CM columns are all-zero although
measured (`CDCA2`, `TEX38`), so an all-zero slot is not evidence of a missing
gene. Only `FAM153B` is absent from the atlas, and it must never be read as a
measured zero (ADR 0031 D2).

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
| `eligibility_decision` | `include` / `exclude` / `uncertain` from #1 rule v1.3 (NaN on controls) |
| `elig_class` | #1's decision class (`approved`, `clinical`, `mechanism_tool`, `other_use`, `unknown`) |
| `qc_pass`, `qc_reason` | #2's per-observation QC |
| `matched_control_ok` | #2's matched-control availability; a control row is its own control |
| `matched_control_key` | #2's matched-control group key (`temp_env::Plate<N>::vehicle=<solvent>`) |
| `main_set` | treated **and** `eligibility_decision == include` **and** `qc_pass` |
| `expansion_set` | treated **and** `eligibility_decision == exclude` **and** `qc_pass` — #11 training only, never scored |

`obs.replicate` is a **generated** within-group index (`cmdata/views.py:50`), not
the raw design; the raw value is `replicate_raw`, and both the `uns`
`column_provenance` block and this table say so. `obs.time_h = 24.0` on every row
is a loader constant (`cmdata/sources.py:539`) whose 24 h design is recorded in
GSE262419, not a raw metadata field.

`matched_control_ok` / `matched_control_key` are carried from #2's
`observations.csv.gz` (`matched_control_status ∈ matched|is_control|not_an_arm`).
When a handoff lacks the column the build derives the rule documented here and
records `matched_control_source: derived` in the manifest; a handoff claiming a
matched control for a group that has none fails the build.

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
| `matched_control_ok`, `matched_control_key` | #2's matched-control decision (`tahoe::<plate>::cell_line_id=<line>`) |
| `tahoe_usable` | `~is_control AND has_structure AND qc_pass AND matched_control_ok` |
| `time_source` | `literature:10.1016/j.cell.2026.08.035` on every row |

`time_h` stays **NaN**: no atlas table carries a time column (ADR 0029 update
2026-11-01, ADR 0031 D3), so the 24 h design is cited as literature provenance
and never written as a constant and never used to filter rows. No row is dropped
from Tahoe: the 1,400 controls, the QC-failed observations and the treated rows
of the two structure-less drugs stay in the artifact for auditability and are
excluded from use through `tahoe_usable`.

Matched controls are **not** materialised as 67k × 18k vectors. Two side files
carry the membership: `tahoe_control_groups.parquet` (the contracted **long**
form — one row per `(group, QC-pass member observation)`: `group_index`,
`group_key` = `plate|cell_line_id`, `plate`, `cell_line_id`, `row`, `n_controls`,
`n_qc_pass`; 649 groups over 1,294 QC-pass control rows) and `tahoe_pairs.npz`
(treated `row` → `group_index`, plus `usable`/`split`; the 4,757 treated rows
with no QC-pass control are absent and counted in the manifest, never
re-paired). The pipeline computes the matched control vector as the mean of that
group's `log1p_cp10k` rows — never across cell lines, never a 50-line average.

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

The held-out rule is not decorative: in the frozen build folds 1–3 each drop one
plate (a plate whose drugs are all test/valid in that fold), so `c*` differs by
fold and uses 15–16 of the 16 plates.

`consensus_membership.parquet` lists every contributing row
(`fold, plate, observation_id, row, weight, plate_weight, n_wells_plate,
qc_pass, vehicle`), where `row` is the 0-based position in `cm_transfer.h5ad`
obs/`X`/layers — the exact name ticket #4 reads. The pool is therefore auditable
and reproducible; `consensus_controls.npz` holds `c_star` (5 × 18,164),
`plate_names`, `plate_weights`, `n_plates`, `n_wells`.

**Leak assertions (fail closed, re-run by `--verify`):** every member row is a
control row; every member row is DMSO and QC-pass; every member row's plate
carries a training drug of that fold. Treated expression of a test/valid drug
therefore cannot be in `c*`, and a held-out plate's controls cannot be either.

**Deployment semantics.** Training sees matched same-plate DMSO controls.
Validation, test and inference receive this single fixed `c*` for their fold.
The experimental truth Δ = `y − mean(log1p_cp10k of matched same-plate DMSO
rows)` is computed **only at scoring time** and never enters a model input.

### 5. Folds — `folds.parquet`

**One row per `(drug_key, fold)`** — the primary key ticket #4 validates
(`drug_key`, `fold`, `role`) — for the `main_set` drugs only
(`role ∈ train|valid|test`). Five folds by `drug_key`, seed **20261102**:
clusters are connected components of `drug_key ↔ inchikey[:14]` (so salt/alias
pairs and shared skeletons cannot straddle a boundary), placed greedily with a
capacity of `ceil(n_main_drugs / 5)` drugs per fold, and
`valid_fold = (test_fold + 4) % 5`. Every main-set drug is test in exactly one
fold and valid in exactly one; every dose and replicate of a drug moves with it.

A drug can sit on more than one plate (459 of the 460 temp_env keys are on one
plate; one spans two), so `plate` holds the drug's plate(s) `;`-joined and
`n_plates` counts them. The per-plate long form is an internal build object —
`fold_sizes.csv` is computed from it, so its plate counts are real — and is not
exported, because a consumer keys the split by `drug_key` alone.

`exclude` and `uncertain` drugs get **no fold row** — they are reserved for the
#11 expansion ablation and are never test or valid. `fold_sizes.csv` reports the
per-`(fold, role)` drug/well/plate counts.

**Drug vs plate.** 459 of the 460 temp_env keys sit on a single plate (one key
spans two), so a fold's test drugs bring their plates with them: a **drug split
is not a plate split**, and a validation/test plate may share the culture
background of training plates. That is why the consensus is scoped by
"plates that carry a training drug", not by "plates without a test drug".

### 6. Locked config — `locked_config.json`

Versioned and hashed (`1.0.0`, digest recorded in `manifest.json`). The **top
level is the canonical pipeline schema ticket #4's
`ptlib.config.normalize_locked_config` reads** (`structure`, `recipe`, `metric`,
`loss`, `optimizer`, `prediction_noise`, `n_folds`, `seeds`, `lrs`,
`dose_encoding`), so the pipeline binds to these values instead of falling back
to its own defaults; the decision record (folds, consensus, validation, budgets,
derived defaults, data contract) is carried beside them.

It fixes the model structure (hidden 128, drug latent 64, state latent 64, noise
10, in_dim 18,164, adaptor 128, FCFP4 1024, dropout 0.05, `comb_num` 1), the
dose encoding (FCFP4 radius 2 / `useFeatures=True` / 1024 bits ×
`log10(dose_µM + 1)`, no binning), the Gaussian NLL loss (never clamped), Adam /
batch 512 / weight decay 1e-8 / 100-epoch cap / min 20 epochs / early-stop
patience 20, the plateau scheduler (factor 0.5, patience 10, **max** semantics
on validation direction AUPRC@200), pretrain LR 1e-3, CM LR candidates
{1e-4, 3e-4, 1e-3}, the three seeds, the per-init per-fold common-LR rule (mean
of the three seeds' validation AUPRC@200), the direction-AUPRC@200 definition
(k = 200; per-condition AP up and down, macro over conditions; within-condition
replicate aggregation first; Tahoe conditions = drug × dose × cell line), the
deterministic prediction policy (noise = 0 at validation/test/inference) and the
30-test-unit policy. Entries the approved spec does not state are marked
`derived_default` with a one-line rationale; no value contradicts the spec.

### 7. Verification — `manifest.json` + `contract.json` + `build_input.py --verify`

`manifest.json` records sha256/bytes of **every** artifact (including
`contract.json`), the frozen source records (content hashes from the 02/04
manifests, `genes.csv` file hash, the GSE262419 design-table hash, the derived
Tahoe gene-space digest), the #1/#2 handoff hashes and row counts (with the
`current.json` pointer and the rule version/sha256), all counts, the
locked-config hash, the git SHA and the runtime environment.

`contract.json` is the machine-readable interface for tickets #4–#11: every
artifact with its shape, keys, columns, dtypes and the meaning of every field,
plus the read rules that must not be re-derived.

`--verify` is offline and re-computes every hash, then asserts the gene axis, the
fold invariants, the consensus provenance and recomputation, the count
reconciliation against the handoffs, the absence of forbidden NaN (Tahoe `time_h`
is the one declared NaN), the control membership/data-plane agreement, and that
no frozen source changed. A stub build is labelled `status: stub` and rejected
unless `--allow-stub`. **The real frozen build passes 125/125 checks.**

## Artifacts and who reads them

| file | git | #4 reads |
|---|---|---|
| `cm_transfer.h5ad` (dense, ~218 MB) | local-only | training/validation/test rows, flags, `main_set`, `expansion_set` |
| `tahoe_transfer.h5ad` (CSR, ~16.2 GB) | local-only | pretrain rows, `tahoe_usable`, `tahoe_split`, identity columns |
| `folds.parquet`, `fold_sizes.csv` | tracked (force-added) | fold/role assignment keyed by `(drug_key, fold)` |
| `consensus_controls.npz`, `consensus_membership.parquet`, `consensus_summary.csv` | tracked | fixed `c*` per fold + its row membership |
| `tahoe_control_groups.parquet`, `tahoe_pairs.npz` | tracked | matched-control mean at load time |
| `locked_config.json`, `contract.json`, `manifest.json`, `counts.csv`, `var_checks.csv` | tracked | recipe, interface, identity, counts, measurement mask |

Neither h5ad carries `uns['control_mean']` (see D0). Raw counts are in `X`, so
any raw-count mean stays reproducible; the correctly scaled caches are the
explicit ones above.

## Reading the 16 GB artifacts (memory)

`anndata.read_h5ad(..., backed="r")` is **not safe** on these files: backed mode
still materialises `layers`, i.e. 8.1 GB for the Tahoe view (1.01e9 nnz) and
12.9 GB for the 25.8 GB corpus. That alone exceeds the WSL2 envelope and is what
killed the 2026-11-02 session. The build therefore reads metadata with
`contract.read_h5ad_header()` (h5py + the anndata element reader; obs/var/uns,
shape, nnz and layer names only), streams the CSR groups with h5py, and never
opens a matrix it does not write. `--tahoe-rows N` reads only N CSR rows.

Measured on this WSL2 VM (31.3 GiB, 8 vCPU; `systemd-run --scope` caps in
parentheses):

| run | wall | peak RSS |
|---|---|---|
| real freeze (`--force`, cap 12 GiB) | 2:15 | **1.73 GiB** |
| `--verify` on the frozen build (cap 8 GiB) | 0:12 | 0.90 GiB |
| `--check` (cap 8 GiB) | 0:04 | 0.55 GiB |
| full unittest suite (cap 8 GiB) | 0:53 | 1.82 GiB |

The CM artifact is dense by source contract (`X` `int32` + log layer `float32`,
3,947 × 18,164 ≈ 574 MB), which is the largest single array the build holds; the
Tahoe side is never densified.

## Reproduce

Run from the repository root. The builder must use the `arc` environment
(anndata 0.12.16); the pinned trainer reads the result with anndata 0.10.8.

```bash
ARC=/home/chang/miniconda3/envs/arc/bin/python
B=stage2/06_prnet_transfer/06.3_input/execution/build_input.py

# 1. validate frozen sources + both handoffs (writes nothing)
$ARC $B --check

# 2. build and freeze into the main checkout's artifact tree
$ARC $B --out-root /home/chang/arc/stage2/06_prnet_transfer

# 3. offline re-verification (hashes + every invariant; 125 checks)
$ARC $B --verify --out-root /home/chang/arc/stage2/06_prnet_transfer

# 4. tests (offline; stub build + verify end-to-end)
make -C stage2/06_prnet_transfer/06.3_input/execution test

# smoke only: synthetic include-all/pass-all handoffs, never a frozen build
$ARC $B --stub-handoffs --out-root /tmp/06_3_stub
```

Inside an implementation worktree (tracked files only) prefix the commands with
`ARC_DATA_ROOT=/home/chang/arc` so the frozen 16 GB sources are read from the
main checkout. `run.sh` resolves the repository root by walking up to `.git` and
execs the same CLI.

## Upstream handoffs (#1 / #2)

The build consumes two handoffs and **fails closed** when either is missing,
partial or malformed.

**#1 eligibility — rule v1.3 / handoff v2.** Resolved through
`06.1_eligibility/output/current.json` (`schema prnet-eligibility-current-v1`):
the pointer declares the artifact version, the rule version/sha256
(`eligibility-rule-v1.3`) and the table path + sha256. The build verifies the
table against the pointer's digest, reconciles the decision counts and the
`include` record list against the table, and records the pointer and rule
identity in the manifest. A future rule revision is consumable by updating the
pointer — no code change, no hard-coded version directory. Record: include 118 /
uncertain 28 / exclude 314 of 460 keys. `v1.2/` is preserved byte-identical and
never rewritten.

**#2 per-observation QC** — `06.2_qc/output/observations.csv.gz`, one 71,576-row
table covering both sources (`source ∈ {temp_env, tahoe}`, `qc_status ∈
{pass, fail}`, plus `matched_control_*` and `raw_replicate`), split by the build.
Per-source tables and `--cm-qc` / `--tahoe-qc` overrides are still accepted.

| handoff | required columns (aliases accepted) | coverage required |
|---|---|---|
| #1 eligibility | `drug_key`, `eligibility_decision` (alias `decision`) ∈ {include, exclude, uncertain}, `elig_class` (alias `decision_class`) | every resolvable CM `drug_key` |
| #2 CM QC | `well` (alias `observation_id`), `qc_pass` (alias `qc_status`), `qc_reason` | every CM view row |
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

## Deviations (D0–D5)

| id | deviation | why |
|---|---|---|
| D0 | Neither artifact carries `uns['control_mean']`; the frozen sources keep their caches untouched. | The 02 cache is a **DMSO+Media raw-count** mean inside one group and the 04 cache is a raw-count mean against log-scale targets (the ADR 0031 recorded deviation). Carrying either would silently invite a raw/log mismatch. Raw counts stay in `X`, and the correctly scaled controls are explicit: per-fold `c*` with row membership, and Tahoe control-group membership computed on the log layer. |
| D1 | Asymmetric row retention: CM keeps QC-pass rows only; Tahoe keeps **all** 67,018 rows flagged. | The CM artifact is the training input and #2's QC is per observation, so dropped rows are counted and listed by #2; Tahoe keeps the structure-less and QC-failed rows so every exclusion stays auditable in one frozen file (`tahoe_usable`). |
| D2 | `gene_measured` is per source and derived from gene-space membership, not from non-zero counts. | 17 Tahoe columns and 2 CM columns are all-zero yet measured; only `FAM153B` is absent from the Tahoe gene space. Deriving measurement from counts would mislabel 19 measured genes. |
| D3 | Tahoe `time_h` stays NaN with a `time_source` literature citation; CM `time_h = 24.0` is recorded as a loader constant. | No atlas table carries a timepoint (ADR 0029/0031 D3), so writing `24.0` would fabricate metadata; the CM constant is not raw GSE metadata either, so its provenance is recorded rather than implied. |
| D4 | `folds.parquet` exports one row per `(drug_key, fold)`, with a multi-plate drug's plates `;`-joined; the per-plate long form is not exported. | Ticket #4 validates `drug_key` as the split unit and rejects a fold where a drug appears twice; a long export would look like two roles. `fold_sizes.csv` still reports per-plate counts, computed from the internal long form. |
| D5 | Metadata is read with `h5py`/`read_elem`, never `anndata.read_h5ad(backed="r")`. | Backed mode materialises the 8.1 GB Tahoe `layers` (12.9 GB for the corpus), which exceeds the WSL2 memory envelope (see *Reading the 16 GB artifacts*). |

## What #4 may rely on

- One gene axis, one fold table, one `c*` per fold, one config — read them, do
  not recompute. `folds.parquet` is the only fold source; `consensus_controls.npz`
  is the only inference-time control; `locked_config.json` is the only recipe;
  `contract.json` is the interface of record.
- `main_set` selects CM training/validation/test rows as treated wells;
  `expansion_set` is #11's training-only arm and has no fold row.
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
