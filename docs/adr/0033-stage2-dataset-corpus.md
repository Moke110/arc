# Stage 2 dataset corpus — one corpus layer, views downstream

**Status:** accepted (2026-11-04)

Stage 2 grew two data steps that both called their product a "corpus" but meant
different things: `stage2/02_cm_data/` owned a matrix-free gene/drug universe
(ADR 0026) plus five filtered h5ad views, and `stage2/04_tahoe_data/` owned a
25.8 GB pseudobulk corpus (ADR 0029) plus a 16.3 GB CM-schema view (ADR 0031).
Neither step's product was the thing the other step's consumers actually wanted:
the substrate between a raw downloaded dataset and a model's training input. The
CM step froze identity but shipped matrices only as pre-filtered views; the Tahoe
step shipped the substrate but only for one dataset. Every adaptation decision
(pharmacology tier, gene-space projection, author shims) therefore lived in a
data step instead of with the model that needed it.

**Decisions.**

1. **One corpus step.** `stage2/02_dataset_corpus/` owns every source dataset as
   an h5ad corpus. `stage2/04_tahoe_data/` is folded in and removed;
   `05_prnet_pretrain/` and `06_prnet_transfer/` keep their numbers, so 04 is a
   gap in the sequence, not a renumber. A missing number is a cheaper
   misdirection than renaming a concluded study's nine substeps.
2. **The corpus is the information-maximal substrate.** One row per raw
   well/profile; `X` = summed **raw counts** (`int32`); the source's **native**
   gene axis; vehicle controls as rows and technical blanks kept; no
   pharmacology filter, no Δ, no integer encoding (`pert_idx`/`cell_idx`/
   `dose_idx`), no 978 subset, no author shims. Five corpora:

   | corpus | rows | of which control / blank | chemicals | gene axis |
   |---|---|---|---|---|
   | `dtoxs_primary` | 779 | 102 / 0 | 47 | 26,328 |
   | `dtoxs_ipsc` | 312 | 16 / 0 | 37 | 26,328 |
   | `sctox` (was `screen_ctox`) | 3,072 | 552 / 0 | 42 | 31,880 |
   | `tenv` (was `temp_env`) | 4,510 | 272 / 32 | 460 | 19,613 |
   | `tahoe` | 67,018 | 1,400 / 0 | 379 | 62,710 |

3. **Native gene axes; the 18,164 intersection becomes a projection target.**
   Each corpus keeps the axis its source can actually measure (the four CM
   sources resolve to 26,328 / 26,328 / 31,880 / 19,613 HGNC symbols). The
   frozen four-source intersection (`genes.csv`, 18,164 symbols) remains the
   **CM gene axis** that a CM benchmark or CM-model view projects onto, not a
   corpus universe. This follows ADR 0029, which already refused to project the
   Tahoe corpus onto the CM axis for exactly this reason, and it lets a future
   model consume genes the 18,164-gene benchmark never had. Rejected: keeping
   the four CM corpora on 18,164 (the old view axis) — it saves ~1 GB and caps
   the substrate at the benchmark axis, so the information loss is permanent
   while the saving is not.
4. **Row set: controls + blanks + every SMILES-resolvable compound.** No
   pharmacology tier is applied; `atc-moa`/`atc-moa-mesh`/`all-smiles` stay
   boolean columns in `compounds.csv` for a consumer to filter on. The 8
   non-SMILES perturbagen keys (four dtoxs biologics, four iPSC antibodies, the
   `ET-1` peptide; 143 wells) remain excluded because they have no structure
   identity, and stay recorded with reasons in `excluded_compounds.csv`.
5. **Views are consumer-owned.** The five CM views of the retired bundle belong
   to `03_models_cmbench/` (currently frozen; they are rebuilt when the CM
   benchmark is redone), and the CM-schema pretrain view belongs to
   `05_prnet_pretrain/` (moved there unchanged, not rebuilt — its
   `content_sha256` is unchanged and `06.3`'s contract stays valid). A consumer
   builds its own frozen input, records the corpus it derives from, and picks
   the tier/gene-space rule itself.
6. **Output is exactly three artifact classes.** `output/corpus/v2/` holds the
   five h5ad corpora, seven crosswalk tables, and one `README.md` carrying
   corpus coverage. Build QA and machine manifests live in
   `execution/build_record/`, not in `output/`. The v1 identity tables keep
   their content (`genes.csv`, `compounds.csv`, `excluded_compounds.csv`,
   `source_gene_sets.csv`); the version bump marks the artifact-class change,
   not an identity change.
7. **Equivalence is content-level, never file-byte.** `.h5ad` bytes are not
   reproducible (h5py defaults `track_times=True`; two identical writes 1.1 s
   apart differ), so a gate compares `content_sha256` — `X` bytes, layer bytes,
   `obs.to_csv`, `var.to_csv` — and records the file sha256 separately. The
   corpus writer passes `track_times=False` so a corpus rebuild is
   byte-reproducible; consumers' view writers need not.
8. **Frozen CM evidence is a temporary facility.** The five v1 views, their
   bundle tables and the retired view builder move to
   `stage2/_frozen/cm-data-v1/` read-only, because the frozen benchmark and the
   concluded transfer study name them as input evidence and their inputs stop
   being produced here. `03`/`06` point at it, nothing new is built from it, and
   **the directory is deleted once 03/05/06 are redone** — it is a bridge, not
   an archive.
9. **Renames are forward-only.** From this decision the dataset ids are `sctox`
   (was `screen_ctox`), `tenv` (was `temp_env`), and `dtoxs_primary`/`dtoxs_ipsc`
   (was the merged `dtoxs_cm`). Frozen records keep their old strings; the
   corpus README carries the rename table, so a reader can translate a 2026-10
   record instead of silently rewriting it.
10. **The tables, not `obs`, are the uniform join layer.** The four CM corpora
    share one `obs` contract (the 24 ADR-0026 CM columns plus `cell_line`,
    `n_raw_rows`, `is_blank`); the `tahoe` corpus keeps its Tahoe-native
    nine-column `obs` (`plate`, `well`, `cell_line_id`, `n_cells`, `cell_name`,
    `drug`, `dose_um`, `is_control`, `control_group`) because it is a 25.8 GB
    artifact whose identity vocabulary is the atlas's, and rewriting it into a
    CM-shaped `obs` would mutate the artifact without adding information.
    Identity joins go through `genes.csv`, `gene_crosswalk.csv`,
    `compounds.csv`, `metadata_cm.csv` and `metadata_tahoe.csv`. "Format
    consistent" therefore means container and value semantics, not an identical
    column list; each corpus's columns are listed in the README.

**Consequences.**

- Until `03` is redone the four CM corpora have **no consumer**. Fitness is
  proved by a reconstruction gate: rebuild the five frozen v1 views from the v2
  corpora (project onto the CM axis, apply the old tier/row rules, map
  `sctox`→`screen_ctox` and `tenv`→`temp_env`), and require each
  `content_sha256` to equal the archived v1 manifest. Any mismatch stops the
  migration instead of being absorbed by a re-run.
- `02` becomes the only Stage 2 step that reads the four CM raw sources; `05`
  gains a data role (it owns the CM-schema view builder and artifact); `06`
  stays path-only — no result, hash or digest key changes.
- The 18,164 gene axis is still frozen, but it is owned by `genes.csv` (a table),
  so a corpus rebuild cannot silently move it.
- Per-plate pseudobulk parts move to `execution/_parts/` (ADR 0030 intermediate
  location) and are kept because a rebuild is measured in hours; the view
  gather cache is deleted because a rebuild is 175 s.

**Considered options.** Keeping two data steps (two "corpora", two owners of the
same conceptual layer — the duplication this record removes); keeping the corpus
matrix-free with views inside the data step (ADR 0026's shape — every model's
adaptation decision then lives in a data step, and a second corpus cannot join);
corpus on 18,164 (decision 3); rebuilding the moved Tahoe view instead of moving
it (identical content, hours of risk for no information).

**See** [`stage2/02_dataset_corpus/README.md`](../../stage2/02_dataset_corpus/README.md),
[`stage2/_frozen/cm-data-v1/README.md`](../../stage2/_frozen/cm-data-v1/README.md),
ADR [0026](0026-stage2-cm-data-contract.md) (amended by this record),
[0029](0029-tahoe-pseudobulk-pretrain-corpus.md),
[0030](0030-external-data-provenance-and-step-layout.md),
[0031](0031-tahoe-cm-schema-view-and-prnet-pretrain.md),
[0027](0027-stage2-benchmark-contract.md),
[0028](0028-stage2-5k-hvg-benchmark.md),
[0032](0032-prnet-formal-transfer.md).
