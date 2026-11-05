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
unknown; no PBS assumption is made for normoxia. The apparent small effects do
not establish a legal FC: matching administered vehicle is unconfirmed. Review
therefore excludes all four normoxic comparisons with matched_vehicle_unconfirmed;
none receives class0. Same-panel0µM establishes a plotted baseline, not verified
vehicle matching. Original v1 extraction bytes remain unchanged; the adapter
interpretation is corrected without rewriting its measurement evidence.

Fig2B includes all four doses, including20µM (early research only confirmed the
lower range):24h pretreatment→4h glucose-free95%N2/5%CO2→24hreoxygenation
95%air/5%CO2. Drug persistence remains unknown. Model comparator is the explicit
H/R+PBS bar~40percent, not normal reference100. Drug means~52/63/77/80percent
would cross the magnitude threshold against this plotted comparator. Matching
administered vehicle is also unconfirmed in the HR panel, so FC is unevaluable.
No numerical p or valid independent n is established either; all four are excluded.
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
real-source regression produce two explicit empty-after-QC normoxic/HR partitions
with all eight comparison/exclusion records. Earlier four class0 admissions were
withdrawn after review identified the unsupported vehicle-matching assumption.

## Duan 2013 primary cardiomyocytes (#36)

Correct publication year is2013 (DOI10.1038/aps.2012.185; PMID23396376),
replacing the earlier candidate inventory spelling Duan-2014. Original PMC HTML
and NCBI-hosted Figure2 are archived unchanged with actual earlier retrieval
provenance. Article identity markers and JPEG format/hash are verified offline.
Failed EuropePMC XML/publisher PDF attempts were error documents and were not
archived as originals. Copyright2013CPS/SIMM is preserved; no CC license invented.

Figure2B contains four target comparisons: NAC200µmol/L and SYA40/60/80nmol/L.
All four approximate author mean values0.82/0.41/0.55/0.81 (±0.02 extraction
precision) are independently checked against the original pixels. The model
bar0.27 is the A/R drug comparator, while untreated normal1.0 is only injury
reference. Mean±SD of three experiments does not establish independent units;
later figures explicitly say independent, but that statement is not imported
into this MTT panel. No summary is expanded into invented replicates and no star
becomes numerical p. All four treatment/model ratios cross the effect threshold,
but necessary adjusted p remains unavailable. Source drug names and exact
administered dose units are preserved, with decimal canonical unit conversion.

The source protocol is24h pretreatment→3h anoxia→6hreoxygenation in neonatal
1–3day Sprague-Dawley rat primary cardiomyocytes. Anoxia uses95%N2/5%CO2 and
specified hypoxia solution, then maintenance medium replacement; actual drug
persistence is not stated. Unknown context is source-scoped. Panel2A is only4h
normoxic toxicity, explicitly excluded from the24h target. SYA chemical-parent
identity (including distinction from hydroxysafflor yellowA) and NAC precise
form/stereochemistry remain unresolved rather than supplied from an assumed
alias. This is visible as an additional identity exclusion.

Offline public build/verify therefore yields four retained comparisons with
explicit exclusions and one typed empty primary/A/R cell_viability partition.
Use `--source duan_2013` with the same preparation/build/verify workflow. The
source archive and real-source public-workflow test both pass. No numerical
corpus sample or independently replicated experiment is claimed for this source.
