# CLAUDE.md — disease_data/miri/

## Purpose

Shared MIRI disease-specific inputs, consumed by the repositioning workflow (`miri_reposition/xpert-fda-v1/`) to build the integrated MIRI reference DEG list. Raw GEO acquisitions and per-GSE DEG lists are reusable across future workflows rather than buried inside one pipeline.

## Structure

```
miri/
├── scripts/
│   ├── map_to_human.R               # Shared ortholog mapping library (mouse/rat/pig -> human L1000)
│   └── map_orthologs.R              # CLI: map rat/pig DEG results to human L1000 lists
├── GSE160516/                       # Mouse MIRI (Affymetrix Clariom S), primary reference
│   ├── data/GSE160516_deg_results.csv   # limma DEG results (GEO acquisition)
│   └── deg/
│       ├── GSE160516_up_l1000.csv       # per-GSE DEG list -> human L1000 (UP)
│       ├── GSE160516_down_l1000.csv     # per-GSE DEG list -> human L1000 (DOWN)
│       ├── GSE160516_MIRI.json          # per-GSE reference gene set
│       └── build_deg.R                  # per-GSE build script
├── GSE108940/                       # Rat MIRI (Agilent GE2-v5)
│   ├── data/GSE108940_WT_deg_results.csv
│   └── deg/{GSE108940_human_l1000.csv, GSE108940_MIRI.json, build_deg.R}
└── GSE137293/                       # Pig scar vs border (RNA-seq; no sham control)
    ├── data/{GSE137293_COUNTS.txt.gz, GSE137293_deg_results.csv}
    └── deg/{GSE137293_human_l1000.csv, GSE137293_scar_vs_border.json, build_deg.R}
```

## Key Files

- scripts/map_to_human.R — function library: `load_landmark_genes()`, `mouse2human()`, `rat2human()`, `species2human()`, `map_deg_to_l1000()`. Paths resolve relative to the checkout via `ARC_L1000_CSV` (defaults to `pert_model/xpert/metadata/l1000_gene_info_978.csv`).
- scripts/map_orthologs.R — CLI wrapper mapping a per-GSE DEG CSV to a human L1000 list (rat/pig).
- GSE160516/deg/build_deg.R — mouse limma DEG analysis + ortholog mapping.
- GSE108940/deg/build_deg.R — rat Agilent DEG analysis + ortholog mapping.
- GSE137293/deg/build_deg.R — pig Welch t-test DEG analysis + ortholog mapping (reproduces the committed per-GSE reference: 2 UP + 1 DOWN).

## Dependencies

- R packages: homologene, limma, org.Hs.eg.db, org.Mm.eg.db, org.Rn.eg.db, annotate, jsonlite (conda env `arc`).
- L1000 landmark gene list: `pert_model/xpert/metadata/l1000_gene_info_978.csv`.

## Entry Points

```bash
# Rebuild a per-GSE DEG list (from the GSE data/):
cd disease_data/miri/GSE160516
Rscript deg/build_deg.R

# Map a rat/pig DEG CSV to human L1000:
cd disease_data/miri
Rscript scripts/map_orthologs.R --deg ... --out ... --species rat|pig
```

## Notes

- Per-GSE DEG lists are the inputs to workflow stage 07 (`build_miri_deg_reference.py`).
- The single 3-GSE R script from the old pipeline was split into per-GSE `build_deg.R` scripts.
- Pig caveat: GSE137293 is scar vs border (no sham control) — use with caution.
- The workflow regenerates the integrated `MIRI_up-down.json`/`MIRI_rank.json` from these per-GSE lists.
