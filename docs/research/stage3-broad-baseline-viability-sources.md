# Broad baseline cardiomyocyte viability sources

Research date: 2026-11-06. Target: more than Heliste's 689 chemicals, baseline cardiomyocytes, MTT/CCK8/CTG. This bounded search does not prove that no larger source exists. No assay datasets were loaded or modified.

## Important large-screen near miss: Huang 2022

Huang et al., *Population-based high-throughput toxicity screen of human iPSC-derived cardiomyocytes and neurons*, Cell Reports 39:110643, doi:10.1016/j.celrep.2022.110643, is a strong broad cardiomyocyte toxicity lead but does **not** establish >689 CTG chemicals. [PubMed](https://pubmed.ncbi.nlm.nih.gov/35385754/), [full primary article and supplement hosted by Axion](https://www.axionbiosystems.com/sites/default/files/resources/mmc2%20%281%29.pdf).

The complete PDF was read on 2026-11-06 (`/tmp/huang2022.pdf`, extracted text `/tmp/huang2022.txt`). Results p.4 says 2,375 approved compounds; STAR Methods p.e5 says a proprietary collection of 2,374 approved/pharmacopoeia drugs. This discrepancy must remain explicit. All 13 HLA-homozygous donor hiPSC-CM lines entered the primary screen, after 24 h treatment at 100 µM. The primary assay uses Hoechst 33342 total nuclei plus propidium iodide dead-cell staining and whole-well high-content imaging. It is **not CTG**.

Fifteen hit compounds entered CTG dose-response validation, at eight concentrations: 0.6, 1, 3, 6, 10, 30, 60 and 100 µM, in triplicate, after 24 h. Nine showed decreased viability across the representative lines. Do not multiply the 15-hit CTG assay across the entire 2,374/2,375-compound screen. The 2,080-compound CTG library described elsewhere in the paper is **neuronal**, not cardiomyocyte data. Likewise, 13 donors increase donor breadth, not chemical count.

Data availability states all data will be shared by the lead contact upon request. A complete publicly downloadable, chemical-indexed primary-screen numerical matrix was not verified. Published supplementary Tables S3/S4 are cardiac hits/validation classifications, not evidence that all screened compounds have public CTG data. No author contact was made.

## Other broader leads and exclusions

- **Liu 2025 Darutigenol screen:** 1,271 natural compounds, one screening dose (10 µM), 48 h, HL-1, CCK8. Methods state 0.2% FBS and a standard 5% CO2 incubator; they do not specify hypoxia for the broad screen. Authors call low-serum treatment “ischemic.” Thus this is an assay-eligible, larger **low-serum normoxic candidate**, but whether it belongs to desired baseline rather than a starvation condition needs an explicit scope decision. Complete compound-indexed numerical results were not verified; Fig.1B is a screen plot. [Primary publication and methods](https://link.springer.com/article/10.1186/s13020-025-01141-x).
- **Zhang 2024 / Nature Protocols 2026 antifibrotic screen:** ~5,000 compounds are screened in iPSC-derived cardiac fibroblasts, then promising agents are counter-screened in cardiomyocytes. The primary paper explicitly says only **20 ADMET-filtered hits** were counter-screened in iPSC-CMs and iPSC-ECs, with a separate valinomycin follow-up. Therefore the ~5,000 library size cannot be assigned to cardiomyocyte viability. [Original Cell study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11645214/), [author protocol](https://www.nature.com/articles/s41596-026-01445-8).
- Broad Tox21 hERG screens use U2OS and channel-function readouts; cardiotoxicity predictions are not cardiomyocyte viability measurements. They do not satisfy the requested cell/phenotype scope. [NCATS primary publication registry](https://tripod.nih.gov/tox/publications).

## Recommendation

Keep **Heliste 2021** as the largest currently established standard-normoxia, assay-eligible source with numerical data already archived. Audit the 1,271-compound low-serum CCK8 screen next if low-serum baseline is acceptable; do not declare its library recoverable until chemical-to-value mapping is demonstrated. Huang 2022 is the strongest new human iPSC-CM breadth lead if the user later includes imaging live/dead viability, but should not be misrepresented as a broad CTG dataset.
