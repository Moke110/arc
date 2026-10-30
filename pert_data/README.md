# pert_data — CM drug-perturbation data lake

The cardiomyocyte perturbation datasets the pipeline reads, split into one directory
per **download source** (项) under the external-data provenance contract
([ADR 0030](../docs/adr/0030-external-data-provenance-and-step-layout.md)). Each item is:

```
<Item>/
├── script/          downloader (download.sh|.py|.md) + offline verifier (verify.sh)   [tracked]
├── source/          the download: bytes unchanged, upstream layout                    [gitignored]
├── local/           declared local additions (generated artifacts, frozen copies)     [gitignored]
└── author_log.json  source URL/version/retrieval time, local additions, verification  [tracked]
```

`manifest.sha256` (per-file hashes over `source/` only) arrives in ticket 09; `script/verify.sh`
already expects it and is a declared no-op until then. Run `script/verify_all.sh` from the
repository root to check every item offline.

## The eight 项

| item | what it is |
|---|---|
| `Tahoe-100M/` | HuggingFace `tahoebio/Tahoe-100M` full atlas (3,388 expression shards, 95.6M cells, 428.8 GB, rev `349e2a55…`, CC0-1.0); the Stage 2 whole-transcriptome pretraining source. |
| `LINCS-LDS-1593/` | The DToxS center's LINCS **LDS-1593** well-level Level-3 package + metadata suite for primary human CMs (manual clue.io download). |
| `DToxS-iPSC/` | The DToxS iPSC-CM GEO deposit, SuperSeries **GSE217424** (GSE217421 Conv + GSE217423 DGE). |
| `GSE244740/` | Rosell-Hidalgo ScreenSeq hiPSC-CM drug-cardiotoxicity panel (2 GEO files). |
| `GSE247541/` | Xu COVID-19 antiviral hiPSC-CM panel; provenance-only, outside the Stage 2 union scope (2 GEO files). |
| `GSE262419/` | Tsai TempO-Seq hiPSC-CM environmental-chemical panel (16 plate CSVs + hash table + series matrix). |
| `XPert-L1000/` | The three L1000 landmark-gene pretrain h5ads used by the workflow pretrain recipes; upstream provenance still unknown. |
| `Annotation/` | Pinned HGNC complete set (download) + the public-registry compound snapshot and CM drug-audit report (local additions). |

Each item's `author_log.json` carries its provenance; each `local/` file is declared in that
log's `local_additions[]`. The downstream corpus built from these sources lives in
`stage2/02_cm_data/` (ADR 0026), and the Tahoe pseudobulk corpus built from `Tahoe-100M/` lives
in `stage2/04_tahoe_pretrain/` (ADR 0029).
