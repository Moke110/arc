# dtoxs_primary: delete raw GEO deposit, release LINCS LDS-1593 flat, freeze v1 data

On 2026-08-18 we finished collapsing `pert_data/dtoxs_primary/` to a single
deposit:

1. **Deleted the GSE146095/GSE146096 raw GEO deposit** — every
   `GSE146095*`/`GSE146096*` file (read-counts, experiment designs, DEG-calc,
   configs, params, series matrices, 2× ~1 GB `Reference-Library.tar.gz`).
   The LINCS LDS-1593 well-level package is a strict superset of that content
   (the 23-KI DGE arm + the CER/IBR/LEN/NIN Conv arm), so nothing unique to
   the pipeline was lost (verified well-level: only non-pipeline content like
   DEG-calc artifacts and reference libraries was dropped).
2. **Released the LINCS LDS-1593 package flat into `dtoxs_primary/`** — the
   `lincs-lds-1593/` wrapper was removed; `Data/`, `Metadata/`,
   `ProcessingPipeline/`, `ReadMe.txt`, `LDS-1593.tar.gz`, and
   `DToxS_summary.md` now live directly under `pert_data/dtoxs_primary/`.
   Workflow paths changed accordingly
   (`DTOXS_LINCS_ROOT = pert_data/dtoxs_primary`; the finetune-data readers
   point at `dtoxs_primary/Data/Level3`).
3. **Froze the v1 data at `pert_data/dtoxs_primary/v1-frozen/`** — a
   byte-identical copy of the well-level package (`Data/` `Metadata/`
   `ProcessingPipeline/` `ReadMe.txt`) plus a copy of v1's frozen finetune
   h5ad (`dtoxs.h5ad`, 521 samples / 36 drugs). This is the data
   `workflow/xpert-fda-miri-v1` consumed; v1's step-02 baseline now reads
   `v1-frozen/Data/Level3` and reproduces its cached baseline byte-identically
   (`58c8610…`). v1's code/benchmarks stay untouched.

Consequences for the siblings:

- v1.1/v2 step-02 baselines now read the released LINCS package with the same
  per-experiment well-position semantics v1 used; they reproduce the v1
  baseline **byte-identically** (previously they differed by float rounding,
  being built from the GSE deposit's CTRL-only wells).
- v1.1 step 03 still rebuilds the DToxS-LINCS finetune data (662 samples / 45
  drugs) and v2 step 03 the 7,165-sample union, unchanged, from
  `dtoxs_primary/Data/Level3`.
- `load_dtoxs_primary_raw` (GSE reader) was removed; the baseline and gene
  index now read the LINCS Level3 tables.

This is the natural end state of ADR 0011 (re-home) — one deposit per study,
with the v1-era data preserved for provenance and v1 reproducibility.
