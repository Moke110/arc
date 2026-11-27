# 02.1 HVG-DEG prediction axes

This step combines variable genes from tenv with disease differential-expression evidence to publish four nested prediction axes for CarPS-Transcriptome. Each axis takes half its genes from the frozen tenv-HVG ranking and fills the remaining slots from eligible tier-first disease DEG1000 sets. The 2k, 3k, 4k and 5k releases retain the original order and symbols on the 18,164-gene common axis. All four filled their disease slots without fallback to additional tenv-HVG genes. [Selection decision](../../../docs/adr/0073-full-source-tier-first-deg1000-and-four-axes.md) · [Release manifest](output/release/manifest.json) · [Axis composition](output/release/evidence/gene-axis/axis_composition.csv)

The 5k axis covered 86.9% of the nHF DEG1000 list, 79.9% of ICM and all 1,000 genes each for genetic cardiomyopathy and TOX. AF evidence remains in the release and did not meet the minimum full-union size for inclusion in fusion. These coverage values describe overlap with the source lists. [Coverage figure](figure/disease-DEG1000-coverage.png) · [Plotted values](figure/data/disease-DEG1000-coverage_bars.csv) · [Render receipt](../log/342-figure-render.json)

[Step 03](../../03_data_train/README.md) published tenv training views for all four widths and the Tahoe well-level view for 5k. [Step 04](../../04_tsg_axis/README.md) compared three HVG-DEG widths with TranSiGen. [Step 09](../../09_cpst_axis/README.md) included all four with CarPS-Transcriptome and selected Denoise-DEG 5k.

The [prediction-axis reader](output/api/prediction_axis.py) projects bounded integer count reads from a [native corpus](../../01_data_raw/README.md) with complete common-gene coverage. See the [operating guide](AGENTS.md) for publication, replay and figure verification.
