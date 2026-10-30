# drug_emb_model — drug biological-embedding model libraries

The drug-side representation models: two external 项 under the external-data provenance
contract ([ADR 0030](../docs/adr/0030-external-data-provenance-and-step-layout.md)). Both
are **code libraries**, so their `source/` trees are tracked; only the large downloaded
assets are gitignored. Each item is:

```
<Item>/
├── script/          acquisition note (download.md) + offline verifier (verify.sh)        [tracked]
├── source/          the author-provided/vendored library as given                          [tracked]
│   └── models|feature_extract/   downloaded assets only in KPGT                           [gitignored]
├── local/           declared local additions (neither item has one today)                 [gitignored]
└── author_log.json  origin/version/licence, asset origins, verification                    [tracked]
```

`manifest.sha256` (per-file hashes over `source/`) arrives in ticket 09; `script/verify.sh`
already expects it and is a declared no-op until then. Run `script/verify_all.sh` from the
repository root to check every item offline.

## The two 项

| item | what it is | version | what it produces |
|---|---|---|---|
| `HG/` | HeteroGNN drug biological-embedding library: contrastive pretraining over a drug–gene graph (PPI, DTI, DDS edges; InfoNCE; 3×SAGEConv, 512→256). **Author-provided** (ADR 0006): no public download, no recorded upstream commit. | `unknown` | The 256-d HG drug embeddings. Pretraining writes `all_drugs_idx2HG.npy` (11,104×256); the finalized Stage 1 workflow instead reads the author-produced `HG_drug_embeddings.npy` in its own `input/` and gates it byte-for-byte (ADR 0008), so the library is dormant in the current pipeline. |
| `KPGT/` | KPGT (knowledge-guided pre-training of graph transformer) drug representation, vendored Apache-2.0. | `unknown` (no reproducible upstream URL/commit) | The 512-d KPGT base drug embeddings used when HG builds FDA-extension node features (`HG/source/src/infer_fda_hg.py`, Ridge KPGT→512). `source/feature_extract/` caches the extracted LINCS2020 / GSE92742 / FDA_drugs tensors (34 MB); `source/models/pretrained/base/` holds the 894 MB pretrained weights. |

## How they feed Stage 1 HG embeddings

Stage 1 `01_drug_embed/` consumes the **author HG embedding table** (`HG_drug_embeddings.npy`,
drug × 256) as the drug representation, copying it into the workflow with a byte-parity gate;
it does not retrain. The HG library in `HG/source/` is the code that produced that table
(`src/pretrain_hg.py`, launcher `run.sh`). KPGT is upstream of HG in the FDA-extension path:
`HG/source/src/infer_fda_hg.py` uses KPGT base embeddings (Ridge, 512-d) as the initial node
features for drugs outside the L1000 graph. Neither item is a Stage 2 model — see
`../pert_model/` for the perturbation-prediction models and `../pert_data/` for the corpora.

## Ignored assets and verification

Only these trees are outside git (ADR 0030): `KPGT/source/models/` (2 weight files,
894,215,942 B) and `KPGT/source/feature_extract/` (6 files, 34,997,345 B). All source code —
including `KPGT/source/src/data/`, which the old global `data/` rule silently ignored — is
tracked. `script/verify.sh` checks each `source/` tree against `manifest.sha256` once ticket 09
produces it.
