# Step 02 prediction axes

Step 02 publishes two families of ordered gene subsets of the fixed 18,164-gene
cardiomyocyte common axis for perturbation prediction. The shared
`prediction_axis.open_axis` reader verifies and returns the ordered genes by
axis ID. Projection preserves source-native integer counts in prediction-axis
order and requires complete common-gene coverage. Step 03 owns the derived
training views. [Prediction-axis API](02.1_hvg-deg/output/api/prediction_axis.py) ·
[Step 01 native counts](../01_data_raw/README.md) · [Step 02 guide](AGENTS.md)

## Axis selection evidence

Step 02.1 **HVG-DEG** published four nested axes, `hvg-2k` through `hvg-5k`.
Each draws half its genes from the frozen tenv-HVG ranking and fills the
remaining slots from eligible tier-first disease DEG1000 sets. All four
filled the disease slots without an HVG fallback. The 5k axis included all
1,000 genes from each of the genetic cardiomyopathy and TOX lists.
[HVG-DEG report](02.1_hvg-deg/README.md) ·
[Release manifest](02.1_hvg-deg/output/release/manifest.json)

Step 02.2 **Denoise-DEG** screened for expression detection and DMSO
variability, retaining 15,071 genes. It ranked those genes by tenv top-DEG
presence weighted by membership in the eligible Step 02.1 disease lists.
The accepted release published five nested axes, `denoise-1k` through
`denoise-5k`. Acceptance verified nesting, all 20 disease-coverage values
and bounded tenv count projection for each axis.
[Denoise-DEG report](02.2_denoise-deg/README.md) ·
[Release manifest](02.2_denoise-deg/output/release/manifest.json) ·
[Acceptance receipt](02.2_denoise-deg/log/acceptance.json)

## Figures

The two disease-coverage plots show overlap with the full eligible DEG1000
lists on a 0–100% scale, with consistent hue by axis size and color depth by
disease. The Denoise-DEG diagnostics show the screening thresholds, score
distributions and axis membership by score. The rendering receipt records
two independent renders and binds all five PNGs to their plotted data and
renderers. [HVG-DEG coverage](02.1_hvg-deg/figure/disease-DEG1000-coverage.png) ·
[Denoise-DEG figures and data](02.2_denoise-deg/README.md#from-screening-to-ordered-axes) ·
[Denoise-DEG coverage](02.2_denoise-deg/figure/disease-DEG1000-coverage.png) ·
[Rendering receipt](log/342-figure-render.json) ·
[Figure verification guide](AGENTS.md#verify-the-current-figures)

## Downstream use

Step 03 published plate-pooled tenv training views for both families at 2k
through 5k and a Tahoe well-level view on `hvg-5k`. Step 04 compared six
axes at 3k through 5k and selected `denoise-5k` for the subsequent TranSiGen
studies. Step 09 compared eight axes at 2k through 5k and selected the same
ordered axis for CarPS-Transcriptome with UniMol.
[Step 03 training views](../03_data_train/README.md#tenv-training-axes) ·
[Step 04 selection](../04_tsg_axis/README.md#gene-axis-selection) ·
[Step 09 selection](../09_cpst_axis/README.md#formal-results)
