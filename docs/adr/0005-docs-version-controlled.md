# docs/ and domain docs are version-controlled

`.gitignore` previously ignored `docs/`, which would have kept the domain documentation (CONTEXT.md, ADRs) and the agent docs out of version control. `docs/` is now un-ignored so these files are committed. `reference/` (papers and transcripts) is tracked as dev-tool documentation. The `xpert-fda-v1_history/` snapshot remains ignored as it contains large regenerable artifacts.
