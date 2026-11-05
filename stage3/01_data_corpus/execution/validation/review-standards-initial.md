## Standards

**One documented-standard finding (medium).** `script/phenotype_archive.py:120–125` creates only per-asset provenance rows, and all five new `pheno_data/*/author_log.json` files omit a source-tree file count and tree digest. ADR 0030 §1 requires the author log to record “per-source file count / bytes / tree digest”; clarification 5 defines `tree_sha256` as SHA-256 of the source manifest bytes. Individual file hashes and byte counts provide useful integrity evidence but do not supply the required aggregate provenance record. Add a source-tree summary with file count, aggregate source bytes and `tree_sha256`, and have offline verification check it against `manifest.sha256`. Preserve original source bytes and existing retrieval records.

No other material documented-standard breaches or baseline smells found in the reviewed diff. Repository-root discovery walks to `.git`; original assets stay separate from declared local additions; public build/verify workflows remain offline; new routing files stay short; source/condition statistics retain their separate provenance. Small explicit source adapters are justified by different source formats and scientific evidence, so their similar record assembly is not flagged as a material duplication smell.

Reviewed fixed diff: `git diff f9209e1859a4153ef4a8b0434e7f83abd81b21cc...8ae83cb01f60c928dc7aa6664d18cd56b8fbf6d7`. Pending final handoff/coverage records were excluded as instructed. Read-only review; no tests or heavy jobs run.

Standards: 1 finding; highest severity medium (missing required source-tree provenance).
