# Step 08 CarPS-Transcriptome chemical representations

Step 08 tested whether changing the chemical representation improved held-out perturbation predictions from the fixed [CarPS-Transcriptome alpha](../07_cpst_alpha/README.md). The seven arms were ECFP4C, ECFP6C, FCFP4C, FCFP6C, E3FP, KPGT and UniMol. Each arm used the same 1,383 plate-pooled tenv conditions from 459 parent compounds, the Denoise-DEG 5k axis, measured same-plate DMSO controls, and the alpha dose transform and objective. Frozen vectors entered at their native widths (512 to 2,304) and each first projection produced 400 features. Downstream parameter initialization and the CPU and CUDA random streams were paired within each cold-drug fold and training seed. [Study design in ADR 0099](../../docs/adr/0099-carps-vae-alpha-and-chemical-representations.md)

## Bounded CUDA gate

All seven arms passed the two-epoch gate on cold-drug fold 0 and training seed 17, using 64 training, 24 validation and 24 held-out conditions from the inherited protocol. Each arm passed independent checkpoint replay and the bounded report passed CSV-only figure reproduction. The KPGT arm reproduced the Step 06.4 bounded validation and held-out predictions, condition and gene ordering, signed Δ arrays and traces exactly. The preflight listed all 63 formal coordinates. [Bounded gate receipt](log/363-gate.json)

## Formal results

The accepted comparison covers all 63 coordinates across folds 0–2 and training seeds 17, 42 and 2026. It includes 54 new fits with exact independent checkpoint replay and nine verified Step 06.4 KPGT reuses. The audit independently verified all 567 scores across nine metrics. Each table entry gives the equal-weight mean ± sample SD across nine fold and seed coordinates. Paired differences subtract the KPGT score for the same fold and seed. Higher Δ PCC and lower Δ RMSE indicate improvement. [Formal acceptance receipt](log/364-formal.json)

| Representation | Δ PCC | Δ RMSE | Paired Δ PCC vs KPGT | Paired Δ RMSE vs KPGT |
| --- | ---: | ---: | ---: | ---: |
| ECFP4C | 0.424806 ± 0.014891 | 0.242734 ± 0.001844 | −0.004024 ± 0.019202 | +0.001982 ± 0.003108 |
| ECFP6C | 0.421715 ± 0.016410 | 0.242380 ± 0.001609 | −0.007115 ± 0.015305 | +0.001628 ± 0.001952 |
| FCFP4C | 0.428691 ± 0.024270 | 0.240448 ± 0.001044 | −0.000139 ± 0.019621 | −0.000304 ± 0.002016 |
| FCFP6C | 0.423095 ± 0.013477 | 0.241313 ± 0.001584 | −0.005735 ± 0.016716 | +0.000560 ± 0.002741 |
| E3FP | 0.435100 ± 0.021071 | 0.238626 ± 0.002190 | +0.006270 ± 0.014760 | −0.002126 ± 0.002559 |
| KPGT | 0.428830 ± 0.018574 | 0.240752 ± 0.002598 | +0.000000 ± 0.000000 | +0.000000 ± 0.000000 |
| UniMol | **0.441262 ± 0.013620** | **0.238599 ± 0.003729** | +0.012432 ± 0.013493 | −0.002153 ± 0.002016 |

**Step 08 retained KPGT as its formal working default.** Changing the default required one representation to lead both primary metrics and each mean paired improvement to exceed its sample SD across the nine pairs. UniMol led both metrics and its mean paired Δ PCC gain (0.012432) fell below its sample SD (0.013493). [Formal configuration decision](log/364-formal.json)

**The successor configuration uses UniMol from Step 09 onward.** ADR 0100 adopted the frozen UniMol representation by project decision under the separate adoption route in ADR 0099. The Step 09 gene-axis comparison verified and reused all nine Step 08 UniMol Denoise-DEG 5k fits. [Successor decision in ADR 0100](../../docs/adr/0100-carps-vae-unimol.md) · [Step 09 results](../09_cpst_axis/README.md#formal-results)

[Fold and seed scores](figure/data/fold-scores.csv) · [Paired summaries by fold](figure/data/paired-by-fold.csv) · [Overall paired summaries](figure/data/paired-summary.csv) · [Fit and checkpoint ledger](log/runs.csv) · [Epoch objectives](figure/data/epoch-objectives.csv) · [MMD bandwidths](log/bandwidths.csv)

## Comparison figures

The learning curves show the accepted training and validation traces for all 63 coordinates, including the nine KPGT reuses. The two metric figures show eight held-out metrics and their paired differences from KPGT. The exporter checks JAC-top200 against the sealed condition and fold scores before deriving wJAC-2000 from the fixed held-out predictions and matched controls. The accepted receipt records two independent CSV-only renders. [Figure protocol in ADR 0102](../../docs/adr/0102-stage2-axis-selection-and-comparison-figures.md) · [Figure receipt](log/figure-render.json)

| Figure | Plot data |
| --- | --- |
| [Learning curves](figure/learning-curve.png) | [Training and validation traces](figure/data/learning-curve.csv) |
| [Held-out metrics](figure/metrics.png) | [Scores by representation, fold and seed](figure/data/metrics.csv) |
| [Paired differences from KPGT](figure/metrics-paired.png) | [Differences by representation, fold and seed](figure/data/metrics-paired.csv) |

The installed `output/formal` bundle contains the sealed matrix, report and KPGT reuse plan. Its manifests are bound by the [formal acceptance receipt](log/364-formal.json). The study guide covers the bounded gate, installed report verification and figure regeneration. [Execution and reproduction guide](AGENTS.md)
