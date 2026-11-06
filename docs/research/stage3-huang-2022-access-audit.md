# Huang 2022 full-screen access audit

Research date: 2026-11-06. Scope: obtain the complete baseline chemical–cardiomyocyte imaging-viability screen after imaging viability was included in `cell_viability`. Read-only web research and inspection of the already downloaded combined article/supplement; no email or author contact.

## Result

**The experimental design is a strong baseline corpus candidate, but the full numerical screen has not been obtained.** The available combined original publication contains methods, plots, Supplementary Figures S1–S6 and Tables S1–S7. It is not a compound-indexed numerical release of the 2,374/2,375-compound screen. The source explicitly makes the underlying data available by request to its lead contact; bounded searches for the exact title/DOI with dataset, repository, Mendeley and Figshare did not establish a public full dataset. This is an access audit, not proof that a public copy cannot exist.

## Obtained and inspected evidence

- [Combined published article plus supplementary information, hosted by Axion](https://www.axionbiosystems.com/sites/default/files/resources/mmc2%20%281%29.pdf): already downloaded as `/tmp/huang2022.pdf`; text inspected at `/tmp/huang2022.txt`. Includes the original article and supplement as one PDF.
- [Original publication record](https://pubmed.ncbi.nlm.nih.gov/35385754/), DOI [10.1016/j.celrep.2022.110643](https://doi.org/10.1016/j.celrep.2022.110643), original PII `S2211124722003953`.
- [Publisher article](https://www.sciencedirect.com/science/article/pii/S2211124722003953): indexed data-availability text agrees with the PDF. Direct web open returned HTTP403; it must not be described as successfully fetched in this run. The data statement promises sharing through the lead contact on request and reports no original code or repository accession.
- Supplement S1/S2 concerns donor HLA information. S3 contains potential cardiac hits as a compound-by-donor boolean hit matrix; S4 classifies validated cardiac hits by clinical indication. S5/S6 concerns neuronal hits; S7 lists primers. These do not supply continuous viability values for the full chemical collection. No hit flag or blank is converted to an exact viability, no-change class or missing physical replicate.

## Fit and missing inputs

The baseline primary screen treats day44-derived hiPSC-CMs, replated for 3 days, for 24 h at 100 µM with individual compounds, then stains total nuclei with Hoechst33342 and dead cells with propidium iodide. The system images the whole well to derive live/dead percentages. All 13 donor-derived cell lines were screened. The original methods give 2,374 compounds while results give 2,375; neither number is yet a count of unique resolved chemical structures. The methods literally state DMSO 10%; preserve this as an unresolved source value rather than silently correcting it. Check the actual source plate protocol with the authors.

A full corpus needs the complete compound library manifest and code-to-name/structure mapping; donor/plate/well/replicate identifiers; per-well total and PI-positive counts or supplied continuous live/dead percentages; vehicle/reference controls and normalization rules; missing/failure flags; and the actual dose/solvent protocol. Raw image pixels are optional if a traceable per-well analysis export with the counts is available. Public figures and the hit list cannot recover this matrix without fabrication. The15-hit CTG dose-response validation is separate and must not be expanded to the complete collection; the 2,080-compound neuronal screen is outside the requested cell type.

## Erratum check

[Published erratum, PMID 39235943](https://pubmed.ncbi.nlm.nih.gov/39235943/), DOI [10.1016/j.celrep.2024.114748](https://doi.org/10.1016/j.celrep.2024.114748), corrects a duplicated calcium-imaging panel in Figure S4. The [author university publication record](https://hub.tmu.edu.tw/en/publications/erratum-population-based-high-throughput-toxicity-screen-of-human/) reproduces the correction text. Figure S4 in the original supplement is neuronal differentiation/activity, so this correction is not evidence of a corrected primary CM viability matrix or a DMSO correction. The corrected figure should accompany the historical source interpretation, but does not resolve full-screen access.

A search-result PII `S2211124722018113` is unrelated: it belongs to *Axonal plasticity in response to active forces generated through magnetic nano-pulling*, DOI 10.1016/j.celrep.2022.111912, as confirmed by the [author institutional repository](https://ricerca.sns.it/handle/11384/128972). It is not a second Huang2022 data release.

## Next acquisition route

Request the full primary cardiomyocyte screen, including non-hits, from the lead contact Patrick C.H. Hsieh, `phsieh@ibms.sinica.edu.tw` (address provided in the publication and PubMed). No contact has been made. Include the above missing inputs and the compound-count/DMSO clarification in one request. Until receipt, classify full-screen acquisition as unavailable-by-public-download/request-needed, not acquired or verified empty.

## Local acquisition and build verdict

The combined article/supplement is now archived byte-for-byte under `pheno_data/Huang-2022/source/`, with a source-only SHA256 manifest, retrieval provenance and an offline verifier. Its active registry status describes this public evidence archive only. The declared `local/source_annotation.json` records experimental meaning, missing inputs and `corpus_build_ready=false`; `local/data-request.md` is an unsent author request. A visual check of original PDF page 21 confirmed 100µM and the literal 10% DMSO, so these are not merely PDF text-extraction substitutions.

**Verdict: suitable experimental design, insufficient acquired data to build the full baseline corpus today.** No measurement rows or corpus release were fabricated from hit flags. The missing all-compound numerical export and chemical mapping are the decisive access gaps; significance tests and endpoint labels are not required. Raw image files are optional if the authors supply traceable per-well counts/percentages.
