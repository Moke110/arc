# Data-lake restructure: flatten pert_data, re-home LINCS LDS-1593, delete dtoxs-xpert

The CM perturbation data lake (`pert_data/`) grouped the five CM datasets under
a `cm/` subdirectory and kept the DToxS-primary study as two overlapping
deposits: the raw GEO deposit (`dtoxs_primary/`, 23 KIs + CTRL) and the
"XPert-curated" well-level deposit (`dtoxs-xpert/`). On 2026-08-17 we removed
the `cm/` grouping (the five datasets live flat under `pert_data/`), downloaded
and verified the LINCS **LDS-1593** well-level package from the LINCS portal
(byte-identical md5/diff match to the local `dtoxs-xpert/`), re-homed it at
`pert_data/dtoxs_primary/lincs-lds-1593/`, and deleted `pert_data/dtoxs-xpert/`.

The LINCS LDS-1593 package is the DToxS center's full well-level Level-3
deposit (956 wells: 23 KIs + non-KI drugs + anthracyclines + antibody/protein
wells) plus the LINCS metadata suite (SMILES, batch codes) — a strict superset
of what the raw GEO deposit holds (GSE146096/146095 cover only the KIs). The
raw GEO deposit stays as-is for the step-02 baseline and provenance; the
finetune-data prep (v1.1 step 03, v2 union step 03) now reads the LINCS wells.
The direct LINCS portal per-level download endpoint returns empty bodies
(Content-Length 0); the `bulk-download` endpoint works and was used instead.
