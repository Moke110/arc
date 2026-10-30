# Download — LINCS LDS-1593 (manual, browser-only)

**method: manual** (no programmatic endpoint works from this host; ADR 0011).
**source:** <https://clue.io/releases/data-dashboard> — LINCS Data Portal releases.
**release:** `LDS-1593` (ReadMe.txt: LINCS Dataset ID `LDS-1593`; DToxS center,
`LDG-1452`). **license: unknown.**

The package is the DToxS center's full well-level Level-3 transcriptome deposit
plus its LINCS metadata suite.

## Manual steps

1. Open <https://clue.io/releases/data-dashboard> in a browser and select the
   **LDS-1593** release.
2. Use the portal's **bulk-download** action to download the dataset package.
   Do **not** use the per-level download endpoint: as recorded in ADR 0011 it
   returns empty bodies (Content-Length 0) from this WSL host; the
   `bulk-download` endpoint is the one that works.
3. Save the archive as `LDS-1593.tar.gz` and extract it so the release lands
   **flat** at `source/` (the `LDS-1593/` wrapper removed, ADR 0013):

   ```
   pert_data/LINCS-LDS-1593/source/
     Data/Level3/                                  per-experiment Level3 tables
     Data/Level3.tar.gz                            the Level3 matrix archive
     Data/PromoCell-DGE-RNAseq-Experiment-Design.DCIC.xlsx
     Metadata/                                     SMILES, batch codes, cell/antibody/protein metadata
     ProcessingPipeline/ProcessingPipelineSpecifications.txt
     ReadMe.txt
     LDS-1593.tar.gz                               the downloaded package itself
   ```

4. `script/verify.sh` checks the resulting bytes against `../manifest.sha256`
   once ticket 09 produces it (offline, `sha256sum -c` from `source/`).

`local/` holds the hand-written `DToxS_summary.md` and the byte-identical
`v1-frozen/` copy of the v1-era package (ADR 0013); neither is part of the
download and neither is hashed by the manifest.
