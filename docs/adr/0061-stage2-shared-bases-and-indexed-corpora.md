# Consolidate Stage 2 foundations and indexed training corpora

**Status:** accepted (2026-11-11). Amended by 0066: the single normalized `X`/`control`/`delta` training corpus is replaced by a raw `corpus/` plus per-axis views; gene identity is unchanged. The user confirmed the complete design and authorized issue publication followed by implementation and deployment in a new chat.

## Ownership and migration

Merge `01_setup` and `02_dataset_corpus` into `01_bases`, which owns shared gene identity, training corpora, HVG ranking, metric calculation and report rendering. Move the current shared evaluator and its significance-building dependencies into this boundary. Central ownership replaces consumer-specific rebuilding of the same foundation, amending ADR 0033's consumer-owned-view boundary and ADR 0059's evaluator placement.

Remove Step 07 and its exclusive code and local artifacts after relocating shared dependencies. Step 08 has no implemented directory; retire its reserved design. Git retains tracked historical evidence. Preserve historical ADRs with supersession pointers when implemented. Remaining old steps 03, 04, 05, 06 and 09 become 02, 03, 04, 05 and 06 respectively. Preserve existing native-axis source corpora and model setup/environment records within separate areas of bases. External source trees remain unchanged.

## Gene identity

Retain the existing frozen CM axis and order, assigning `gene_index` 0 through 18163. This axis originated as the four-CM-source intersection; Tahoe coverage was subsequently verified. The shared identity boundary maps symbols, HGNC, Ensembl and other available source identifiers to these indices, retaining source-specific identifiers, aliases and mapping evidence. One HGNC identity may have different source Ensembl identifiers (for example FAM153B). Ambiguous mappings fail explicitly.

Formal Stage 2 gene references and metric internals use the shared indices. Raw-source identities, native model interfaces and human-facing display labels cross this boundary through the identity table. Historical sealed evidence retains its original identity representation.

## Two separate training corpora

Publish separate `tahoe_bulk.h5ad` and `tenv.h5ad` files. Retain one treated well per sample; a Tahoe sample is one cell-line pseudobulk within one well. Each sample includes treated expression, matched-control expression, delta, cell-line identity, desalted drug-parent identity and dose, together with source provenance.

Tenv contains one cell identity, `hiPSC-CM`, as confirmed by the user. Assign each retained sample a stable unique `sample_id`, preserving source identity so separate wells remain distinguishable. The tenv corpus contains expression and sample metadata. Statistical significance labels live in a separate file keyed by `(desalted drug parent, dose)`, not by sample ID.

Unify drug-parent identity while preserving individual wells, source salt forms and experimental conditions. Normalize reliably convertible doses to micromolar, retaining the original values and conversion evidence; exclude unresolved doses with reasons. Parent identity alone does not make different times, vehicles or experimental designs compatible for a joint statistical comparison.

Project raw counts to the full shared axis before per-profile library-size normalization to 10000 and natural `log1p`. Average the independently normalized matched controls: same-plate DMSO for tenv, same plate and cell line DMSO_TF for Tahoe. Delta is treated expression minus this control mean. Exclude samples without a compatible matched control and record the reason. Existing `shared-cm-v1` supplies an already compatible expression foundation; some legacy Tahoe experiments used a different normalization and remain historical results rather than evidence for the new inputs. This migration does not automatically retrain all historical experiments.

## Significance and evaluation

Only tenv has ranked significant upregulated and downregulated gene-index lists, stored in a separate label file keyed by `(desalted drug parent, dose)`; the h5ad itself carries no significance labels. `tahoe_bulk.h5ad` has no significance-label build. A real significance list describes the compatible replicate-well contrast for one drug-parent–dose group. Use a plate effect for multi-plate contrasts. Use raw-count DESeq2 two-sided Wald tests, within-contrast BH adjustment and `padj < 0.05`. Direction is the DESeq2 log2 fold-change sign. Within each direction sort by ascending adjusted p-value, descending absolute log2 fold change and ascending gene index. Insufficient replication or failed testing is untestable, distinct from a tested contrast with no significant genes.

The test and BH family is exactly the 18164-gene axis, replacing ADR 0059's CM-plus-HDGFL3 family. Rebuild affected tenv labels and rescore affected retrieval reports. Retrieval remains tenv-only and retains ADR 0059's AUPRC, AUPRC_positive and AUPRC_negative calculation semantics. Aggregate predicted well deltas by mean within `(desalted drug parent, dose)` groups and rank the group predictions: absolute delta for combined retrieval, positive delta for positive retrieval and negative-delta magnitude for negative retrieval. Prediction-side significance testing and significance filtering are not performed. Only the measured truth lists require statistical significance. Keep direction matching, tied-magnitude handling, zero abstention, full-axis recall denominators, exact group AP, macro averaging and fixed-population pooled PR curves. Indexed arrays support shared metric computation without repeated symbol joins.

Ordinary training loaders do not load significance-label files. Keep all replicate wells of a statistical contrast in one train/validation/test partition; contrasts do not span partitions. Evaluation accesses the relevant groups' labels. Future significance-supervised training would require a separately generated training-partition release. Model-facing HVG rankings also use training partitions only.

## HVG

HVG rankings carry no perturbation direction. Supply a full-tenv descriptive ranking and training-partition-only rankings for model gene selection. Use Seurat normalized dispersion, breaking ties by gene index. The shared ranking function owns the calculation; held-out expression does not participate in model-facing rankings.
