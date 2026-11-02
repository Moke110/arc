# 06.2 · qc — per-observation experimental design + QC audit (ticket #2)

Status: **built and verified 2026-11-02** — two sources audited from raw plate/well
provenance to a keep/exclude/unresolved observation list, with the rule set and its
thresholds **locked before any model result** (`model_results_seen: false`; no model has
been trained or evaluated anywhere in `stage2/06`). This ticket produces counts and a QC
verdict; it does **not** decide pharmacology eligibility (ticket #1) and does not build
folds or the consensus control (ticket #3).

Headline: **Tahoe 1,344 physical wells → 67,018 `(plate, well, cell_line)` observations**
(60,815 kept) and **temp_env 4,558 physical wells = 4,558 observations** (4,021 kept;
3,745 of 4,206 treated wells). Counts, every exclusion and every unmatched control are
in `output/qc_summary.json` + the small tables; the full per-observation table is
`output/observations.csv.gz` (gitignored, sha256 in `output/manifest.json`).

## 1. Physical wells are not observations

| source | physical wells | observation | why the difference |
|---|---:|---|---|
| temp_env | 4,558 (16 plates × 286 wells; plates 4/14 have 277) | 4,558 | one TempO-Seq library per well; well ≡ observation |
| Tahoe | 1,344 (14 plates × 96 wells) | 67,018 | each well pools ~50 cancer lines and is pseudobulked **per line**: 1,344 × 50 = 67,200 minus unassigned combinations; 28 DMSO_TF wells → 1,400 control observations |

Consequences that this audit enforces and that ticket #3 must not undo:

- The **50 lines of one Tahoe well are not 50 biological replicates**. A well is one
  culture; its 50 profiles share that well's cells, seeding and sequencing depth.
- **Two control observations of one `(plate, line)` group are two wells, not two
  donors.** The atlas carries no donor field at all (`donor` is NaN).
- Tahoe `obs.well` is a composite (`tahoe::plate1:smp_1495:CVCL_0399`); the physical well
  is `smp_1495`, carried separately in the table.
- temp_env's frozen-view `obs.replicate` is a **generated index** (equal to
  `Replicate_Number − 1` on 99.665 % of rows). The table carries the source's own
  `Replicate_Number` as `raw_replicate` and the frozen index only as
  `frozen_view_replicate_index`, flagged as generated.

## 2. The locked rule set

Full text, per-rule evidence (citation + URL + quote + retrieval date 2026-11-02) and the
instantiated derived values are in `output/qc_rules.json` (`rules_version`
`qc-rules-v1.0.0`, `rules_sha256`; the file's own sha256 is in `manifest.json`).
`--build` recomputes the payload and refuses to run if any threshold moved after locking.

| id | statistic | threshold | evidence | removes |
|---|---|---|---:|---:|
| gene filter (pre-step) | mean raw count per probe over all 4,558 wells | keep gene if ≥ 10 | source: GSE262419 series matrix `!Sample_data_processing` | 22,537 → 11,497 probes |
| **T1** sample expression floor | mean count per gene after the gene filter | exclude if < 15 | source: same field ("Samples with an average expression across all genes lower than 15 were also removed") | 17 wells |
| **T2** vehicle reproducibility | mean Pearson r of a DMSO well vs the other 191 DMSO wells, raw counts | exclude control if < 0.85 | source: same field + paper ("mean correlation threshold of <0.85") | 41 wells |
| **T3** blank separation | `n_genes_detected` over the 22,537-probe panel | exclude if ≤ 9,242 = blank median 7,067.5 + 3·1.4826·MAD 489 | derived: the 32 technical Blank wells are the only same-run no-cell reference; 3-MAD band from Iglewicz & Hoaglin (1993) | 499 wells |
| **H1** profile size | `n_cells` | exclude if < 100 | derived: ADR 0029's measured density argument; source floor 50 (paper) kept as the sensitivity bound | 5,251 obs |
| **H2** source line exclusion | cell line | exclude NCI-H661, NCI-H596, NCI-H2122 | source: paper ("three cell lines ... consistently exhibited low counts ... were excluded") | 3,995 obs |
| **H3a** mito ceiling | per-well mean mitochondrial fraction | exclude if > 20 % | source: paper's full cell filter (< 20 % mito) — verification gate | 0 (max 14.1 %) |
| **H3b** genes per cell | per-well mean detected genes | exclude if < 250 | source: paper ("at least 250 genes detected") — verification gate | 0 (min 732) |
| **H4** plate depth outlier | per-plate modified z of log10(mean transcripts/cell) | exclude well if \|z\| > 3 | derived: the paper's own ±3 UMI z convention applied at well level with a MAD scale | 490 obs (10 wells) |

H3a/H3b removing nothing is the expected result: the atlas applies those filters per
cell upstream, so the per-well means must sit inside them. They are recorded as
verification gates rather than dropped from the rule set.

Thresholds *never* used (asserted absent from every rule's evidence by the test suite):
DEG counts, effect size / response strength, model performance. QC-passing **low-response
and toxicity/stress wells are kept** — nothing in the rule set reads a response statistic,
so a weak or strong responder is treated identically to any other well of the same
technical quality.

Two readings were possible and are both reported: T1 over the retained genes (locked,
17 wells) vs over all 22,537 probes (290 wells, `sensitivity_T1_all_panel_genes_removed`);
and dropping T3 (`sensitivity_kept_without_T3` = 4,502 wells). Neither is used for the
locked verdict.

**The correlation statistic had to be fixed to reproduce the study.** TempO-Seq
correlations in log space sit at 0.60–0.75 among vehicle controls (a 0.85 floor would
discard ~80 % of the study), while raw-count correlations sit at 0.90 median — so the
locked rule uses raw counts over the author gene filter. Even so, the study applies it to
vehicle controls only; applying it to treated wells would make the filter
response-correlated (a strong responder genuinely sits further from vehicle), so treated
wells are **flagged, not excluded** (1,277 flags, `vehicle_corr_below_source_floor`).

## 3. Counts

### temp_env (4,558 physical wells, 16 plates)

| role | wells | pass | fail | what it is |
|---|---:|---:|---:|---|
| treated compound | 4,206 | 3,745 | 461 | 3 doses × 1,402 wells |
| DMSO vehicle control | 192 | 138 | 54 | the only vehicle control (12/plate before QC) |
| Media | 80 | 64 | 16 | **not** a vehicle control: different arm, never merged into a DMSO pool |
| Blank | 32 | 29 | 3 | technical blank, no cells; **never** a control |
| ET-1 (0.01 µM) | 48 | 45 | 3 | peptide **positive control**, not a small-molecule treatment |

Treated kept by dose: 0.2 µM 1,249 / 1 µM 1,226 / 10 µM 1,270. Per-plate eligible DMSO
after QC: 3–11 (median 10); every plate keeps at least one, so **zero** kept treated wells
lack a matched control.

### Tahoe (1,344 wells, 14 plates)

| role | observations | pass | fail |
|---|---:|---:|---:|
| treated | 65,618 | 59,521 | 6,097 |
| DMSO_TF control | 1,400 | 1,294 | 106 |

Kept by plate 4,030–4,477. Sensitivity at the source's 50-cell floor: 61,679 kept
(`n_cells_sensitivity_kept_at_source_floor_50`). 1,498 kept-because-large-enough profiles
sit between 50 and 100 cells and carry the flag
`n_cells_between_source_and_project_floor`.

## 4. Matched-control availability (the handoff keys)

- **temp_env**: key `(plate, vehicle=DMSO)` — a treated well is matched only to DMSO wells
  of the same physical plate. Media and Blank never enter the pool. 3,790 kept arm rows
  (treated + ET-1) all resolve to ≥ 1 eligible control; per-row counts in
  `matched_control_n`, per-plate totals in `output/temp_env_dmso_controls.csv`.
- **Tahoe**: key `(plate, cell_line_id)` — never across lines. 59,521 kept treated
  observations: 59,148 have 2 eligible controls, 169 have 1, and **204 have none**
  because 51 of the 700 `(plate, line)` groups lost both control profiles to QC. Those 204
  are reported as unmatched, **not** re-paired to another group and **not** imputed.

## 5. Plate ↔ drug confounding

`459 of 460` temp_env `drug_key`s appear on exactly **one** plate (462/462 raw names on
one plate; one key spans two). Tahoe is the opposite — 345 of 355 drug keys appear on more
than two plates (`output/plate_drug_confounding.csv`). Consequence, stated for ticket #3:
a drug-key split is **not** a plate split. Every fold's training drugs still come from most
physical plates, so the per-fold DMSO pool is "plates touched by the training drugs" (all
sixteen here), and a held-out plate's controls must never be described as a held-out
environment. Any effect attributed to a fold is an effect of the drug set, not of a plate.

## 6. Labels: low response and toxicity stay `unknown`

`response_label` and `toxicity_label` are `unknown` on **every** row, with the reason in
`label_evidence`. No per-well viability or response readout exists in this repository:

- temp_env's study did measure cytotoxicity, but in a **separate 90-minute** phenotypic
  assay, not on the 24 h transcriptomic wells (PMC11691792); those values are not deposited
  here either.
- Tahoe's atlas reports proliferation/cytotoxicity as a quantified phenotype and ships
  `n_cells` per profile, but no per-well viability label is deposited.

Neither label may be inferred from transcriptomic magnitude, and none is. Anything that
wants a toxicity stratum must pre-register its criterion (tickets #3/#9) and say so.

## 7. `FAM153B`: measured zero vs missing

| source | status | evidence |
|---|---|---|
| temp_env | **measured** (a zero is a real zero) | probe `FAM153B_22242` in the 22,537-probe panel; 54,881 counts, non-zero in 2,755 / 4,558 wells |
| Tahoe | **missing** (never a measured zero) | the CM slot maps to `ENSG00000182230`; the atlas `gene_metadata` maps the symbol `FAM153B` to `ENSG00000289731` and does not contain `ENSG00000182230` |

Ticket #3 must carry a boolean `gene_measured` vector on the 18,164-gene axis (False for
Tahoe rows at this slot, True for temp_env rows) instead of reading the slot as zero.

## 8. Time provenance (24 h, both sources)

- temp_env: `time_h = 24.0` comes from the frozen view, where it is a **hardcoded loader
  constant**, not a GEO field. Basis: "Following chemical exposure for 24 hrs, the wells
  were completely aspirated ..." (Chem Res Toxicol 2024, PMC11691792). Recorded as
  `time_source_id = tempenv_paper_24h`.
- Tahoe: the atlas has **no time column**, so `time_h` stays **NaN** — no constant is
  written. Basis: "After 24 hours of drug treatment, spheroids were dissociated, fixed, and
  profiled using the Parse GigaLab kit" (bioRxiv 2025.02.20.639398v2, the authors' text of
  the same study as Cell 2026 doi:10.1016/j.cell.2026.08.035). Recorded as
  `time_source_id = tahoe_preprint_24h`.
- Both records, with URL, quote and retrieval date, are in `output/qc_rules.json`
  (`time_sources`). The Cell full text returned HTTP 403 to every fetch from this host on
  2026-11-02; the quote is verified in the authors' v2 preprint, and the DOI is cited for
  the formal record. Tahoe is **not** filtered on `time_h`.

## 9. Handoff to ticket #3

`output/handoff.json` is the machine-readable interface: kept ids per source, the
matched-control join keys and eligible control rows, the QC flag columns, the gene mask,
the routing of Media/Blank/ET-1, the label policy, and what this ticket deliberately does
not cover. Companion tables:

| file | rows | content |
|---|---:|---|
| `observations.csv.gz` | 71,576 | the full per-observation table (gitignored; sha256 in `manifest.json`) |
| `temp_env_dmso_controls.csv` | 192 | every DMSO well per plate + `eligible_for_control_pool` |
| `tahoe_dmso_controls.csv` | 1,400 | every DMSO_TF profile per `(plate, line)` + eligibility |
| `matched_control_availability.csv` | 748 | per stratum: controls total/pass, treated kept, kept-but-unmatched |
| `counts_by_source_role.csv`, `counts_by_plate_dose.csv`, `counts_tahoe_by_cell_line.csv` | 14 / 250 / 97 | the per-source / per-plate / per-dose / per-line counts |
| `exclusions_by_rule.csv`, `primary_reason_counts.csv` | 8 / 5 | every exclusion, by rule and by primary reason |
| `plate_drug_confounding.csv` | 3 | the confounding statement as counts |
| `qc_rules.json`, `qc_summary.json`, `manifest.json` | — | locked rules, headline counts, sha256 of everything |

Reading note: `qc_status` is `pass` / `fail` / `unresolved` (`unresolved` = a required QC
metric was missing; 0 rows today). `qc_reason` is a `;`-joined list of reason codes whose
first element is the primary one; `qc_flags` never changes a verdict. Empty-string fields
(`vehicle` on treated Tahoe rows, `drug_key` on controls, `donor = 'NA'` on temp_env)
become NaN under a plain `pandas.read_csv` — same convention as the frozen views.

## 10. Reproduce

```bash
cd stage2/06_prnet_transfer/06.2_qc/execution
PYTHONDONTWRITEBYTECODE=1 /home/chang/miniconda3/envs/arc/bin/python -m unittest discover -s tests -v   # 47 tests
./run.sh --lock      # derive + freeze thresholds      (~20 s)
./run.sh --build     # apply them, write every table    (~40 s)
./run.sh --verify    # re-hash + count invariants       (~5 s)
make record          # copy the tracked record into 06.2_qc/output/ (then git add -f)
```

`QC_DATA_ROOT` (default `/home/chang/arc`, the main checkout that holds the gitignored
`pert_data/` and `stage2/02|04` outputs) and `QC_OUT_ROOT` (default = data root) are the
two overrides; `run.sh` is a thin wrapper over the three `build_qc.py` modes. Peak RSS is
under 1.5 GB: the 16 GB Tahoe view is opened with `backed='r'` and only its `obs` is read —
no matrix is ever materialised.

## 11. Deviations

| id | deviation | why / effect |
|---|---|---|
| D0 | temp_env kept **4,021 of 4,558**; the study reports 4,379 (96 %) | T3 has no published counterpart and the study publishes only the aggregate, not per-sample membership, so the two cannot be matched sample-by-sample. Effect: stricter QC; dropping T3 keeps 4,502 |
| D1 | T1 reads "average expression across all genes" as the genes retained by the author gene filter | the two steps are listed in that order; the all-panel reading removes 290 wells and is reported as sensitivity, not used |
| D2 | the 0.85 correlation floor is applied to DMSO controls only, in raw-count space | log-space correlations cannot reproduce the study's >0.85 vehicle median; extending the floor to treated wells would filter on response strength, so those wells are flagged only |
| D3 | Tahoe profile floor is 100 cells, not the source's 50 | ADR 0029's measured density argument (below 100 cells the pseudobulk is 10–20× shallower than a CM well). Count at 50 reported (61,679) |
| D4 | Tahoe `time_h` remains NaN | the atlas has no time field; the 24 h claim is recorded as a sourced `time_source`, never as an obs constant |
| D5 | all response/toxicity labels are `unknown` | no per-well readout in the repo; the upstream readouts exist but are not deposited (see §6) |
| D6 | 204 kept Tahoe treated observations have **no** matched control | 51 of 700 `(plate, line)` groups lost both control profiles to QC; reported, not imputed |

## 12. Open gaps

- Per-sample membership of the study's own 179 removed temp_env samples is not deposited,
  so our stricter filter cannot be reconciled row by row (D0).
- `temp_env` per-plate physical plate format is not verified beyond the 286-well layout
  (384-well coordinates, rows B–N); nothing in this audit depends on it.
- The Cell version of the Tahoe paper could not be machine-read from this host; the 24 h
  and QC quotes are from the authors' v2 preprint of the same study (D4).
- This ticket does not decide pharmacology eligibility, folds, dose inclusion or the
  consensus control — those are tickets #1 and #3.
