# corpus/v1 — frozen Stage 2 CM gene and drug universes

**Frozen:** 2026-10-21

**Source commit (vendored Stage 1 loaders):** `c921f7d`

**Genes:** 18,164

**Compounds:** 525 drug keys (532 key/name rows)

## Dataset ids

The GSE sources use a short Stage 2 id; the raw folder keeps the GEO
accession:

| dataset_id | accession | raw source |
|---|---|---|
| `dtoxs_primary` | — | `pert_data/dtoxs_primary` |
| `dtoxs_ipsc` | — | `pert_data/dtoxs_ipsc` |
| `screen_ctox` | GSE244740 | `pert_data/gse244740` |
| `temp_env` | GSE262419 | `pert_data/gse262419` |

## Gene axis

Four-source intersection (`dtoxs_primary`, `dtoxs_ipsc`, `screen_ctox`,
`temp_env`). Each source identifier resolves to exactly one current
HGNC symbol:

1. DToxS symbols are Entrez-verified through the package's own
   `UCSC-hg38-Gene-HGNC-Entrez-Map.tsv` (an Entrez ID wins when HGNC
   has reassigned the symbol).
2. `screen_ctox` Ensembl IDs (version stripped) map through
   `ensembl_gene_id`.
3. `temp_env` `SYMBOL_probeid` rows are keyed by symbol; duplicate
   probes collapse to one gene and the numeric suffix is only an
   Entrez fallback (it is a probe id, not an Entrez id).

Unresolvable or ambiguous identifiers are dropped and counted in
`../README.md` and the bundle `qa_report.md`; nothing is guessed.
The order is the HGNC symbol, lexicographic, and is authoritative for
every bundle `var`.

Match types: {"alias": 1086, "direct": 17067, "entrez": 11}

## Drug set

`drug_key` is the desalted parent InChIKey **with stereochemistry**.
`inchikey` is the exact registered structure key, so salt and free-base
rows share a `drug_key` but keep distinct `inchikey` values. The three
tiers are boolean columns:

- `pharm_atc_moa`: ATC code or ChEMBL mechanism;
- `pharm_mesh`: the above or a PubChem MeSH class;
- `pharm_all_smiles`: a resolvable small-molecule structure.

Tier sizes: atc-moa 177, atc-moa-mesh 316, all-smiles 525 drug keys.

## Inputs

- HGNC complete set: `pert_data/annotation/hgnc_complete_set.txt`
- compound registry snapshot:
  `pert_data/annotation/cm_compound_annotation.csv`
- source data: the raw folders in the table above

Both annotation assets are pinned with URL, retrieval date and
checksum in `pert_data/annotation/README.md`; the build makes no
network call.

`source_gene_sets.csv` is the frozen per-source resolved gene set: a
rerun that maps a new gene in any source fails, even if the gene is
outside the four-way intersection.

## Per-source gene resolution

| source | direct | alias | previous | entrez | ensembl | unmapped | ambiguous | resolved |
|---|---|---|---|---|---|---|---|---|
| dtoxs_primary | 23,410 | 70 | 877 | 2,047 | 0 | 12,073 | 1 | 26,328 |
| dtoxs_ipsc | 23,410 | 70 | 877 | 2,047 | 0 | 12,073 | 1 | 26,328 |
| screen_ctox | 0 | 0 | 0 | 0 | 31,880 | 12,274 | 0 | 31,880 |
| temp_env | 21,733 | 3 | 705 | 17 | 0 | 79 | 0 | 19,613 |

See the bundle `gene_crosswalk.csv` for the raw identifier of every
gene and the bundle `qa_report.md` for per-view gates.
