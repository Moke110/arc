# FDA divergence report — v2 (HG-aligned) vs v1

Captured after re-running the v2 FDA screening tail (steps 06-11) on the author-pinned FDA-extension HG rows. Divergence is the point of the HG-aligned feature, not a failure.

## Predicted DEG divergence (step 06)

- Per-drug Pearson r between v1 and v2 FDA DEG: **0.982** (mean)
- Max |diff|: 0.541; mean |diff|: 0.0054

## Screening funnel

| Stage | v1 | v2 |
|---|---|---:|
| FDA drugs screened | 3127 | 3127 |
| WCS < 0 | 347 | 410 |
| WCS < 0 + >= 6/8 critical pathways | 40 | 62 |

## Integrated ranking divergence (step 11)

- Ranked candidates: v1 **40**, v2 **62**; overlap **5 / 97**.
- Integrated-score correlation across all drugs: **0.64**.
- Mean |delta integrated|: 0.35; max: 2.18.

Top 10 candidates:

| # | v1 | v2 |
|---:|---|---|
| 1 | PROCHLORPERAZINE (+1.220) | PYRVINIUM PAMOATE (+1.279) |
| 2 | DISOPYRAMIDE (+1.126) | QUINIDINE SULFATE (+1.167) |
| 3 | DIPHEMANIL (+1.064) | DIPHEMANIL (+1.147) |
| 4 | BISOPROLOL (+0.974) | DOXEPIN (+0.933) |
| 5 | EFINACONAZOLE (+0.945) | MAPROTILINE (+0.906) |
| 6 | ESCITALOPRAM (+0.889) | BROMPHENIRAMINE (+0.828) |
| 7 | METHAPYRILENE FUMARATE (+0.859) | GALANTAMINE (+0.777) |
| 8 | METHOXAMINE (+0.827) | METHAPYRILENE FUMARATE (+0.755) |
| 9 | FROVATRIPTAN SUCCINATE (+0.795) | NALOXONE (+0.715) |
| 10 | MIGLUSTAT (+0.766) | BENFLUOREX (+0.706) |
