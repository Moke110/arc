# Public baseline imaging viability search — Burnett 2021

Research date: 2026-11-06. Bounded primary-source check for a baseline cardiomyocyte `cell_viability` source with more than 689 chemicals and publicly downloadable full numerical measurements. Imaging viability is now eligible; 24 h is a preference, not an absolute filter. No full assay dataset was acquired in this audit.

## Burnett 2021: strong size/model match, full measurement download not established

*Cardiotoxicity Hazard and Risk Characterization of ToxCast Chemicals Using Human Induced Pluripotent Stem Cell-Derived Cardiomyocytes from Multiple Donors* reports 1,029 chemicals, five human iPSC-CM donor lines, four final doses (0.2, 1, 10 and 100 µM), 0.5% DMSO, no oxygen-injury intervention. Chemicals include pharmaceuticals and environmental/industrial compounds, not 1,029 therapeutic drugs. [Original article](https://pmc.ncbi.nlm.nih.gov/articles/PMC8762671/), [publisher](https://pubs.acs.org/doi/abs/10.1021/acs.chemrestox.1c00203).

Original methods were retrieved through the [NCBI BioC full-text endpoint](https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC8762671/unicode), saved temporarily at `/tmp/burnett2021.xml`. Ca flux is measured after **90 min** of chemical exposure. Immediately afterwards, medium/dye is removed; Hoechst33342 and MitoTracker Orange stain is applied for 15 min, washed out and imaged using DAPI/TRITC/FITC filters at 10×. This is an acute imaging dataset, not 24 h viability. The study's cytotoxicity phenotype is a decrease in total cell number relative to vehicle, with a modeled EC10 as the derived POD. The cell-count measurements would be eligible under imaging viability but are distinct from a direct PI-based live/dead fraction.

The [publisher Supporting Information listing](https://pubs.acs.org/doi/abs/10.1021/acs.chemrestox.1c00203) advertises a PDF containing plate maps, a test-chemical identifier/class/source list, donor characteristics, efficacy values for active compounds and Hill-model definitions. This listing does **not** establish a full all-chemical, all-dose raw or normalized viability table. The full-text BioC rendering does not contain an explicit raw-data download/repository link. Publisher Figshare API metadata lookup returned403 in this run. Thus full measurement access remains **unverified**, rather than proven unavailable.

## Do not confuse related accessible files with this full screen

- [Dryad 10.5061/dryad.612jm641k](https://datadryad.org/dataset/doi%3A10.5061/dryad.612jm641k) offers a17.32MB ZIP of raw screening data and R scripts, but belongs to Blanchette2020 and covers136 compounds across43 donors, not this1,029-compound experiment. Its accessibility cannot establish this source's full data availability.
- A2026 deep-learning reanalysis of the1,029-compound screen has a [public author code repository](https://github.com/Danny-Vu/Syracuse-TAMU), listing Autoencoder.py, Full_Pipeline_Code.ipynb, Preprocessing.py and README.md. The reported task analyzes calcium transients and reconstruction-error scores, not imaging viability. Code availability is not full viability-data availability. [Peer-reviewed reanalysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC13547733/).

## Verdict

Burnett2021 meets the chemical-count and cardiomyocyte-model criteria and could support a separate acute-baseline imaging corpus if continuous per-well/chemical data can be downloaded. This audit has not established that required full numerical download, so it must not yet replace Heliste as an acquired corpus source. Never substitute active/inactive calls, modeled PODs or calcium scores for the original viability readout values.

## Supplement fetch audit

The indexed original PMC record lists exactly two attachments: `NIHMS1767198-supplement-Supplemental_Table_1.pdf` (182.5KB, chemical manifest) and `NIHMS1767198-supplement-All_supplemental_data.docx` (523.5KB). Its advertised content matches the publisher's plate-map/donor/Emax/Hill-model description above. Direct publisher SI PDF fetch returned403. The observed PMC DOCX URL downloaded1,817bytes of HTML, not a ZIP/DOCX, and was rejected as non-data; an official S3 asset attempt returned404. These are access failures, not evidence that the SI is empty or a successful inspection of its binary contents. No full continuous viability matrix was verified. [Original indexed PMC attachment inventory](https://pmc.ncbi.nlm.nih.gov/articles/PMC8762671/).
