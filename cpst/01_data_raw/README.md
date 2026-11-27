# CarPS-Transcriptome raw data

Step 01 publishes nonnegative integer counts in each dataset's native gene
order. It preserves source observation identities, vehicles, technical blanks
and unmatched treatments. The release records measured coverage of the frozen
18,164-gene common axis and requires complete coverage for tenv and Tahoe bulk.
Stored counts remain on the native axis. [Data ownership decision](../../docs/adr/0085-stage2-raw-and-training-data.md)

## Accepted corpus

| Dataset key | Observation unit | Observations | Native genes | Common genes measured |
| --- | --- | ---: | ---: | ---: |
| `tenv` | GSE262419 plate well | 4,510 | 19,613 | 18,164 |
| `tahoe_bulk` | Tahoe-100M plate, well and cell line | 67,018 | 62,710 | 18,164 |
| `dtoxs_ipsc_conv` | Published DToxS-iPSC conventional sample | 1,171 | 38,478 | 17,112 |
| `dtoxs_ipsc_dge` | Published DToxS-iPSC DGE sample | 344 | 38,478 | 17,112 |
| `dtoxs_primary` | LINCS LDS-1593 observation across eight publisher tables | 842 | 38,478 | 17,112 |
| `sctox` | GSE244740 ScreenSeq well | 3,072 | 44,154 | 18,164 |
| `gse247541` | Published sample replicate | 15 | 60,656 | 18,160 |

tenv retains 32 technical blanks and keeps the same-plate DMSO and Media
groups distinct. Tahoe bulk uses same-plate, same-cell-line DMSO_TF controls.
DToxS-primary matches treatments within their original table and compatible
solvent context. Its `source_observations.csv` traces 842 accepted
observations to all 956 publisher rows, including the original table and
control context. Shared observations passed count equality after alignment
by native gene ID. The accepted primary corpus contains 41 matched treatments,
699 unmatched treatments and 102 vehicles. [Migration receipt](log/migration-303.json)

The installed `output/manifest.json` pins the seven dataset manifests and
their H5AD payloads. The tenv inventory also pins the probe map and source-build
manifest. XPert-L1000 lacks verified native integer counts. Annotation supplies
reference metadata only. [Perturbation data sources](../../pert_data/README.md) ·
[Release inventory and digests](log/migration-303.json)

The original source replay verified tenv against its pinned counts and
Tahoe bulk across all 14 plates and all native genes. The v2 migration
compared counts and checked observation and gene metadata for the six
datasets carried forward from v1, and verified the new DToxS-primary
membership ledger. The relocation receipt records the original audit paths.
[Original acceptance](log/acceptance.json) ·
[Migration acceptance](log/migration-303.json) ·
[Relocation record](log/migration-306.md)

## Downstream use

Step 02 projects native counts onto ordered prediction axes. Step 03 uses
the accepted corpus for eight plate-pooled tenv training releases and
the Tahoe well-level HVG-DEG 5k release. Use the source-native reader for
bounded count windows, observation metadata and matched vehicle rows.
The step workflow gives the reader setup and the accepted-release and
candidate verification routes. [Prediction axes](../02_gene_axis/README.md) ·
[Training data](../03_data_train/README.md) ·
[Reader](output/api/native_corpus.py) · [Step workflow](AGENTS.md)
