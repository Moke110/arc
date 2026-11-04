# `_frozen/cm-data-v1/` — frozen v1 CM evidence (TEMPORARY BRIDGE / 临时设施)

> **Temporary bridge. 临时设施。**
> This directory exists only because the frozen `03_models_cmbench` benchmark and
> the concluded `06_prnet_transfer` study still name the retired v1 CM bundle as
> their input evidence, while `02_cm_data` no longer produces it (ADR 0033
> decision 8). It is **deleted once `03` / `05` / `06` are redone** — it is a
> bridge, not an archive.
>
> **Read-only. Nothing new may be built here.** No build writes into this tree,
> no new artifact may be derived from it, and it must not grow. The data and
> evidence files below (`datasets-1/`, `corpus/v1/`) are byte-identical moved
> originals. The only additions are this README and the archival code copy in
> [`execution/`](execution/README.md), whose files carry an added three-line
> provenance header (the copied bodies are unchanged).

## What is here

| path | contents | tracked |
|---|---|---|
| [`datasets-1/`](datasets-1/) | the 5 retired v1 view h5ads (18,164 genes each), plus `compounds.csv`, `gene_crosswalk.csv`, `manifest.csv`, `qa_report.md` | h5ads + `gene_crosswalk.csv` are gitignored (`*.h5ad` / `*.csv`) and untracked by design; the rest are tracked |
| [`corpus/v1/`](corpus/v1/) | the frozen v1 identity tables: `genes.csv` (18,164 symbols), `compounds.csv`, `excluded_compounds.csv`, `source_gene_sets.csv`, `README.md` | tracked |
| [`execution/`](execution/README.md) | archival copy of the retired view-construction code (`views.py`, the view/tier constants from `paths.py`, the tier-nesting gates from `qa.py`) — evidence only, not an entry point | tracked |
| `README.md` | this file | tracked |

Moved from `stage2/02_cm_data/output/{datasets-1,corpus/v1}/` on 2026-11-04;
tracked files with `git mv`, the untracked h5ads and `gene_crosswalk.csv` with
`mv` (same filesystem, a rename — content untouched). The bundle move is `dev`
commit `58b4d84`; the retired step's Tahoe outputs moved to
`02_dataset_corpus/` and `05_prnet_pretrain/` in `27c1180`/`d6c275a`. `7e7d689`
is the **pre-move fixed point** — the state the archival code in
[`execution/`](execution/README.md) was copied from — not the commit that
performed the moves.

## The five views

`manifest.csv` is authoritative and was written by the v1 build. Sizes and
`content_sha256` values below were **recomputed after the move** with the exact
v1 recipe and equal the manifest row-for-row (verified 2026-11-04; peak RSS
1.5 GiB).

| view | shape | bytes | `content_sha256` | matches `manifest.csv` |
|---|---|---|---:|---|
| `dtoxs_cm.h5ad` | 1,075 × 18,164 | 101,966,077 | `17198963e925c74231eb0ffaf5d876b691f15fd458359cc543091ec367743ac1` | yes |
| `screen_ctox.h5ad` | 3,072 × 18,164 | 180,624,560 | `c0fffab49d0a2a91acabb307b51a4899d2f1a7c804d09ab7d9e34436ae88f749` | yes |
| `temp_env_atc-moa.h5ad` | 1,319 × 18,164 | 115,488,071 | `91591d58dd63aed0866508ca098c2ffb4ff8423208f3cb63f2d0e3f5dccc0af5` | yes |
| `temp_env_atc-moa-mesh.h5ad` | 2,582 × 18,164 | 178,860,257 | `4633eaef47ea3860ff2b79b180e7968418fa1f0ce9d62dfe589fa329964b2c05` | yes |
| `temp_env_all-smiles.h5ad` | 4,478 × 18,164 | 289,629,094 | `f84066c4dfd2c017785d85a3f6bd54cb2a8eb40bfba28e9e99d650e4b8a613fb` | yes |

Total 866,568,059 bytes (826.4 MiB). The content digest is the v1 recipe from
`02_cm_data/execution/build_cm_datasets.py:_content_digest`, unchanged by the
move:

```text
sha256( X.tobytes()
      + layers["log1p_cp10k"].tobytes()
      + obs.to_csv(index=False, lineterminator="\n")
      + var.to_csv(lineterminator="\n")
      + uns["control_mean"].to_csv(lineterminator="\n") )
```

File bytes are not reproducible (h5py `track_times`); only the content digest is
a contract. The views are dense-int32 `X` + dense-float32
`layers["log1p_cp10k"]`; the schema is the 24-column ADR-0026 `obs` contract.

## Who reads this directory

| consumer | what it reads | why |
|---|---|---|
| `03.1_cold_drug_18k/execution/script/build_data.py` and its tests | `datasets-1/temp_env_all-smiles.h5ad`, `datasets-1/manifest.csv`, `corpus/v1/genes.csv` | the frozen benchmark's input; its `output/` stays untouched |
| `03.2_hvg_sensitivity_5k/execution/script/build_gene_set.py`, `build_data.py` and their tests | `datasets-1/{dtoxs_cm,screen_ctox,temp_env_all-smiles}.h5ad`, `corpus/v1/genes.csv` | the frozen 5k-gene HVG axis was ranked on these views |
| `06.1_eligibility/execution/{evidence.py,build_eligibility.py}` and its tests | `datasets-1/compounds.csv` | the 460-key eligibility source; the view-membership column is `temp_env_all-smiles` |

`03`/`06` read this tree **read-only**. Their frozen `output/` records and
hashes are not rewritten, and no result is recomputed. Later tickets
(`#14`/`#18`/`#19`) repoint the remaining `06.2`/`06.3`/`04` producers; the old
`stage2/02_cm_data/output/...` strings inside frozen records stay historical.

## Provenance

- v1 corpus frozen 2026-10-21, source commit `c921f7d` (vendored Stage 1
  loaders); five views built 2026-10-30 from the four CM raw sources under
  `pert_data/`.
- Moved read-only here on 2026-11-04 by ticket #13 (parent spec #12, ADR 0033
  decision 8). No data or evidence file was regenerated, edited or re-hashed
  except as recorded above; the archival `execution/` copies add only a
  three-line provenance header, leaving the copied bodies unchanged.
- Deletion condition: `03`/`05`/`06` are redone and no record points here any
  more. Do not delete earlier: the frozen benchmark would lose its evidence.
