# Step 03 plate-matched training data

Step 03 publishes expression with measured matched vehicle controls from
[source-native integer counts](../01_data_raw/README.md) on the ordered
[prediction axes](../02_gene_axis/README.md). The tenv releases pool wells
by condition. The Tahoe release retains individual well and cell-line profiles.

## Source-specific releases

The first tenv release used HVG-DEG 5k. It pools DMSO wells by plate and
treated wells by plate, canonical drug parent and measured dose. Among 4,510
source observations, 112 were excluded, leaving 16 pooled controls and 1,383
pooled treatment conditions. The original build and its independent replay
matched byte for byte. [03.1 tenv report](03.1_tenv_hvg-5k/README.md) ·
[Source-specific acceptance](log/312-source-specific-releases.json)

The Tahoe release preserves 66,718 eligible well and cell-line profiles on
HVG-DEG 5k, including vehicle profiles. It pools distinct matched DMSO_TF
observations from the same plate and cell line in integer counts, then applies
CP10K and natural `log1p` separately to the well and its matched pool on the
selected axis. Δ subtracts the normalized control from the normalized well
profile. The source-specific acceptance records byte-exact replay of the
original combined view and digest equality for the three Tahoe assets.
[03.2 Tahoe report](03.2_tahoe_hvg-5k/README.md) ·
[Tahoe manifest](03.2_tahoe_hvg-5k/output/release/manifest.json) ·
[Source-specific acceptance](log/312-source-specific-releases.json)

The native corpus publisher uses Tahoe to check the gene order and endpoint
profiles before replacing the accepted corpus. The release also passed the
archive-isolated active-input smoke test. The formal Step 04 TranSiGen and
Step 09 CarPS-Transcriptome axis comparisons used the tenv releases below.
[Corpus compatibility check](../01_data_raw/execution/build_corpus.py) ·
[Archive-isolated smoke receipt](../04_tsg_axis/log/323-archive-isolated-smoke.json)

## Tenv training axes

Each release selects its genes before pooling counts and uses each pooled
profile's selected-axis total for CP10K and natural `log1p`. The released Δ
subtracts the normalized same-plate DMSO control from the normalized treatment.
Changing the axis can therefore change the normalized values of shared genes.
Every release binds the same native counts and the same three cold-drug and
three cold-dose-drug folds. The cold-dose-drug test sets use the intersection
of held parents and held doses. Rows matching exactly one held criterion stay
unassigned, with the reason recorded.
[Tenv release contract](03.1_tenv_hvg-5k/README.md) ·
[Frozen tenv folds](output/frozen-folds/tenv/manifest.json)

Five releases expanded the original HVG-DEG 5k baseline to six axes at 3k,
4k and 5k. Acceptance recorded byte-exact replay and confirmed that the ordered
condition keys, matched controls, raw membership, exclusions, fold assignments
and strip reasons matched the baseline. Step 04 selected Denoise-DEG 5k from
these six axes for the subsequent TranSiGen loss and architecture studies.
The CarPS-Transcriptome alpha and chemical representation study then used
the same release.
[Six-axis acceptance](log/313-six-tenv-axes.json) ·
[Step 04 axis selection](../04_tsg_axis/README.md#gene-axis-selection) ·
[Downstream use of Denoise-DEG 5k](03.7_tenv_denoise-5k/README.md#downstream-use)

Two 2k releases then passed the same replay and population checks. An
independent audit recomputed every pooled row from native counts and checked
gene order, library totals, normalized expression, matched Δ and the
dataset-wide control vector. Step 09 evaluated all eight axes with UniMol
and selected Denoise-DEG 5k for CarPS-Transcriptome. Step 10 held that axis
fixed for its dose comparison. [2k acceptance](log/366-two-tenv-axes.json) ·
[Step 09 formal results](../09_cpst_axis/README.md#formal-results) ·
[Step 10 formal results](../10_cpst_dose/README.md#formal-results)

| Tenv axis key | Release |
| --- | --- |
| `hvg-5k` | [03.1 HVG-DEG 5k](03.1_tenv_hvg-5k/README.md) |
| `hvg-3k` | [03.3 HVG-DEG 3k](03.3_tenv_hvg-3k/README.md) |
| `hvg-4k` | [03.4 HVG-DEG 4k](03.4_tenv_hvg-4k/README.md) |
| `denoise-3k` | [03.5 Denoise-DEG 3k](03.5_tenv_denoise-3k/README.md) |
| `denoise-4k` | [03.6 Denoise-DEG 4k](03.6_tenv_denoise-4k/README.md) |
| `denoise-5k` | [03.7 Denoise-DEG 5k](03.7_tenv_denoise-5k/README.md) |
| `hvg-2k` | [03.8 HVG-DEG 2k](03.8_tenv_hvg-2k/README.md) |
| `denoise-2k` | [03.9 Denoise-DEG 2k](03.9_tenv_denoise-2k/README.md) |

The named tenv reader resolves all eight axes through the shared pooled
training reader. The original six-axis lookup retains the roster and runtime
pins used by the sealed Step 04 export and the Step 07 and 08 releases.
[Named reader](output/api/tenv_axes.py) ·
[Build and replay procedure](AGENTS.md#build-and-replay-a-tenv-release)

## Historical evidence

The paired dose study predates these source-specific releases and used the
axis-local views specified in ADR 0077. ADR 0085 records the original
publication boundary for plate-matched data. The CarPS-Transcriptome index
links the tracked historical records and current active-input checks.
[Dual-source dose study decision](../../docs/adr/0078-dual-source-axis-local-dose-round1.md) ·
[Original publication decision](../../docs/adr/0085-stage2-raw-and-training-data.md) ·
[Historical evidence and active-input checks](../README.md#historical-evidence-and-active-input-checks)
