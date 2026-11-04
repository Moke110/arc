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

## Liu 2021 icariin (#35)

DOI10.1002/2211-5463.13276, PMC8564343.1, article XML and actual Fig2
`FEB4-11-2966-g004.jpg` were acquired and archived byte-for-byte. The deceptively
numbered g002 is Fig5 (Znpp mechanism combinations), excluded from this adapter.
CCBY4.0 is recorded in original XML. Archive timestamps refer to the actual prior
primary downloads; a subsequent stalled live refetch did not replace these bytes.

Fig2A four nonzero doses2.5/5/10/20µM,24h normoxia, H9C2 rat cell line:
approximate percent values99/96/100/98 from pixel calibration (±3percent).
Control is actual0µM panel-A control100, and Methods explicitly says absorbance
is expressed as control percentage. Its precise vehicle composition remains
unknown; no PBS assumption is made for normoxia. Source small effects admit four
operational class0 rows without p. These are **magnitude_failed_without_p**,
not demonstrated nonsignificant or equivalent effects.

Fig2B includes all four doses, including20µM (early research only confirmed the
lower range):24h pretreatment→4h glucose-free95%N2/5%CO2→24hreoxygenation
95%air/5%CO2. Drug persistence remains unknown. Model comparator is the explicit
H/R+PBS bar~40percent, not normal reference100. Drug means~52/63/77/80percent
cross the magnitude threshold against this comparator. No numerical p or valid
independent n is established, so all four are excluded rather than labeled0.
Published SD and three repetitions do not establish independent unit structure;
ANOVA/Tukey stars do not become exact p/BH. A conflicting Results sentence says
normal viability was unaffected but gives P<0.01; this is retained as an
uninterpretable source statement and supplies no numerical test evidence.

All eight summaries were independently checked by `/root` against the original
figure and coordinate table; normalized summaries remain one measurement each.
CCK8 maps to cell_viability while preserving method and units. Unknown vehicle
or injury persistence scopes condition identity to this source. Exact decimal
micromolar conversion avoids binary-arithmetic source key drift.

Use the same commands above with `--source liu_2021`. Public build/verify and
real-source regression produce one4-row normoxic partition and an explicit empty
HR partition, with complete comparison/exclusion evidence.
