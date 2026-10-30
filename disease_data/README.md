# disease_data — shared MIRI disease inputs

The disease-side gene sets the Stage 1 MIRI screen scores against, split into one
directory per **GEO series** (项) under the external-data provenance contract
([ADR 0030](../docs/adr/0030-external-data-provenance-and-step-layout.md)). Each item is:

```
<GSE>/
├── script/          per-series build_deg.R + shared map_to_human.R / map_orthologs.R + verify.sh  [tracked]
├── source/          upstream GEO-derived acquisition files (bytes unchanged)                        [gitignored]
├── local/deg/       locally derived DEG lists (JSON tracked; regenerable CSV ignored)              [mixed]
└── author_log.json  GEO accession, derivation chain, local additions, verification                 [tracked]
```

`manifest.sha256` (per-file hashes over `source/` only) exists for every item — 453.0 GB hashed
across the four categories, 27/27 verified. `script/verify.sh` checks it offline, and a
missing manifest is a failure. Run `script/verify_all.sh` from the repository root to check
every item.

## The three 项

| item | organism | comparison | source/ | local/deg/ |
|---|---|---|---|---|
| `GSE108940/` | rat (*Rattus norvegicus*), Agilent GE2-v5 | WT MIRI (I/R) vs WT Sham, limma | `GSE108940_WT_deg_results.csv` (2,766,695 B) | `GSE108940_MIRI.json` (tracked), `GSE108940_human_l1000.csv` |
| `GSE137293/` | pig (*Sus scrofa*), Illumina NovaSeq | scar vs border zone, Welch t-test | `GSE137293_COUNTS.txt.gz` + `GSE137293_deg_results.csv` (2,582,175 B) | `GSE137293_scar_vs_border.json` (tracked), `GSE137293_human_l1000.csv` |
| `GSE160516/` | mouse (*Mus musculus*), Affymetrix Clariom S | pooled MIRI (I/R) vs Sham, limma | `GSE160516_deg_results.csv` (2,004,508 B) | `GSE160516_MIRI.json` (tracked), `GSE160516_{up,down}_l1000.csv` |

## Rat vs human (read before using GSE108940)

MIRI (myocardial ischaemia-reperfusion) is modelled in rodents and pigs, not humans.
`GSE108940` is a **rat** series: its human gene set is a derived artefact, produced by
rat→human ortholog mapping (homologene via `org.Rn.eg.db` Entrez IDs, see
`script/map_to_human.R:rat2human`) and then filtered to the L1000 978 landmark genes.
`GSE160516` is mouse (the workflow's primary rank reference) and `GSE137293` is pig
(scar vs border — an injury-severity gradient, **not** disease vs healthy). All three
`local/deg/*.json` lists are human-symbol L1000 sets; the species is recorded in the JSON
metadata and in each `author_log.json`.

## How a new GEO series would be added

1. Create `disease_data/<GSE>/` with `script/`, `source/`, and `local/deg/`.
2. Put the upstream-derived acquisition files (GEO supplement / series matrix output) in
   `source/`; do not edit downloaded bytes.
3. Copy both shared helpers (`script/map_to_human.R`, `script/map_orthologs.R`) into
   `script/` and add a per-series `script/build_deg.R` that reads `source/` and writes
   `local/deg/` (the helpers resolve the repo root by walking up to `.git`).
4. Write `author_log.json` (accession, organism, derivation chain, `local_additions[]`) and
   `script/verify.sh`; declare every generated `local/deg/` artefact in `local_additions[]`.
5. Add one line to `AGENTS.md` and one row to the table above; hash the new `source/` with
   `python3 .scratch/data-provenance-restructure/ops/hash_items.py --only disease_data/<GSE>`.

The downstream consumer is the Stage 1 MIRI step `08_miri_phenotype/miri_phenotype.py`,
which reads `disease_data/<GSE>/local/deg/*.json` and `source/*_deg_results.csv`.
