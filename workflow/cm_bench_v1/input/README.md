# input/ — frozen benchmark assets

Generated once by the `_lib/fetch_*.py` scripts and read by the benchmark (the
files themselves are large data and are gitignored, per repository convention).

- `ensembl_to_symbol.tsv` — Ensembl gene ID → HGNC symbol for `gse244740`'s
  genes, derived from the Ensembl GRCh38 release-104 GTF.
  Regenerate: `python _lib/fetch_gene_symbol_map.py`
- `gse244740_smiles.tsv` — SMILES for `gse244740`'s compounds (recovered from
  the DToxS small-molecule metadata / PubChem; the series matrix carries no
  SMILES).
  Regenerate: `python _lib/fetch_gse244740_smiles.py`
