# Reviewed viability figure sources

## Linn 2024 (#34)

DOI10.3390/ijms25031427, PMC10855687.1, original article/metadata/Fig3
archived in `pheno_data/Linn-2024/source/`. CC BY4.0 evidence is retained in
the original article. Figure3b is48h cytotoxicity, separate from infected-cell
antiviral experiments. MTT adds3h color development after drug treatment.

All48 nonzero black data points (three drugs × two models × eight points)
were digitized on the original797×636 JPEG. The red fitted CC50 lines are not
measurements. Values are approximate author-normalized technical-well means,
rounded to1percent with±3percent extraction precision. No technical wells were
invented from these means. The reported pointwise n3 is technical; separate
CC50 n3 independent experiments does not establish a pointwise independent n.
Exact dose grids are absent from text: every point keeps log-axis coordinates,
calibration and approximate dose±15%, while canonical dose remains null.
The six curves, controls or absent zero-control plot, failed exact-dose mapping,
source image hash and independent review are retained in versioned curation.

Independent reviewer `/root` inspected original images and each extracted value.
Review corrected a missing leftmost hiPSC remdesivir point and its43pixel/decade
calibration, yielding eight points in every curve. This verifies approximate
summary digitization, not exact source values, raw replicates or numerical p.

Dose uncertainty prevents all48 comparisons from admission; matched vehicle
composition is also unavailable, so a non-treated100percent reference is not
silently treated as proven matched vehicle. Exact administered salt is unknown.
No pointwise numerical p/BH evidence was published. Public build/verify produces
two typed empty-after-QC cell_viability partitions with all comparisons retained.

Reproduce (offline after acquisition):

```bash
python pheno_data/Linn-2024/script/download.py
python stage3/01_data_corpus/execution/viability_adapters.py --source linn_2024 --prepare --output /tmp/linn-bundle.json
python stage3/01_data_corpus/execution/corpus.py build --input /tmp/linn-bundle.json --output /tmp/linn-release
python stage3/01_data_corpus/execution/corpus.py verify --release /tmp/linn-release
```

`--prepare` explicitly materializes the tracked extraction as a hash-declared
local addition. Ordinary adapter generation/build/verify read only frozen files.
