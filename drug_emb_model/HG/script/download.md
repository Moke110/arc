# Download — HG

**method: author-provided (manual)** — there is no public download: the HG
library is author-given ([ADR 0006](../../../docs/adr/0006-hg-aligned-fda-extension.md))
and its upstream version is unknown (`author_log.json` records `version: "unknown"`).

## What is here

`source/` holds the HeteroGNN drug biological-embedding library as provided:

- `source/src/` — model + training/inference code; `source/AGENTS.md` is the
  library's own file-content index
- `source/run.sh` — training launcher
- `source/run.ipynb` — interactive notebook

No download command exists; the library was authored/assembled locally and has no
recorded upstream commit. It produces the Stage 1 HG drug embeddings that the
XPert workflows consume (see the category `README.md`).
