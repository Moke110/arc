# Stage 1 v2 cutover verification

The comparison baseline was the materialized v2 screen at Git commit
f0d5bb673435b1a0341301c802ceac5094708a54. Verification ran on 2026-11-27
in a separate worktree whose Stage 1 root had no versioned sibling directories.
Numbered steps used the shared arc conda environment and existing local data
provenance. Outputs regenerated during testing were isolated from the original
formal products. The validation environment used Python 3.13.13, NumPy 2.4.6,
pandas 2.3.3, AnnData 0.12.16, h5py 3.16.0, Matplotlib 3.10.9,
PyTorch 2.12.0, Scanpy 1.12.1 and umap-learn 0.5.12.

## Frozen inputs and paths

- The 286 original artifact hashes are recorded in
  [pre_cutover.sha256](pre_cutover.sha256). The
  [checkpoint inventory](pre_cutover_weights.tsv) records 52 working links to
  the external XPert weight archive.
- The two predecessor H5ADs match [the local reference manifest](input/references/manifest.sha256).
  The schema input has 521 rows and the fold source has 662 rows.
- The local Annotation seed matches [its offline checksum](../pert_data/Annotation/script/audit_seed.sha256).
  `python pert_data/Annotation/script/fetch_annotation.py --compounds` regenerated
  547 rows across 535 names. All keys and fields matched the installed
  Annotation table.
- The historical GSE108940 input is Git blob
  3ebed362193f9eed49616b7e648ad775e61112bf from commit
  712b144cf1161e07ef9d2256c126b7453eba865e. The corrected shared source changes
  Stage 1's union from 78 to 79 upregulated landmarks and changes WCS results.
  The Stage 1 pin reproduces the accepted analysis and its original labels.
- All 36 workflow configs resolve their inputs through the XPert config loader.
  Their scientific values match the v2 configs.

## Data, official folds, and final ranking

The [reproduction guide](README.md) lists the commands for each step.

| Check | Observed result |
|---|---|
| Drug list and HG embeddings | 1,004 FDA compounds. The output HG table matched the original author reference. |
| Gene mapping gates | Passed for all five CM datasets, including DToxS primary and iPSC drug-level audits. |
| Rebuilt baseline | The Step 02 control CSV matched the original byte for byte. |
| Rebuilt quantile reference | All 4,890,000 sampled quantile values matched the original exactly. |
| Rebuilt union and schema gate | 958 rows, 978 landmarks, 45 drugs. Sample/gene order, the expression array, the matched control array, and shared observation fields matched exactly. The complete H5AD matched byte for byte after the official folds were regenerated. |
| Official fold replay | All 21 assignment columns matched for every sample, including ten warm-start folds, ten cold-drug folds and the all-train assignment. |
| Stage 07 z-scoring and Steps 08–11 | Using the frozen raw prediction, the z-score table, both phenotype JSONs, all score/filter/ranking CSVs, and the seven rerendered PNGs matched exactly (tables by bytes, PNGs by pixels). Funnel counts were 1,004 screened, 130 with WCS below zero, and 64 finalists. The Step 11 report differed only in its generation date. |
| Step 12 survey HTML | All 20 reviewed dossiers assembled with zero flags. The 14-dossier subset rebuilt from those frozen inputs. |
| Composition report | The cached normalization matched the prior calculation bit for bit on 130 pinned L1000 and CM profiles. The full report matched byte for byte after processing 336,852 L1000 profiles and 958 finetune profiles, including its two embedded figures. |
| Data ablation | Recomputed the masks and aggregated all 33 frozen evaluations. Both CSVs, the PNG and the report matched byte for byte. |
| Warm-start PR | Recomputed the tables, figure and report from ten frozen fold profiles. All outputs matched byte for byte. Mean direction-aware AUPRC was 0.2234, with sample SD 0.0048. |
| Step 13 limitations | Both CSVs and all three PNGs matched byte for byte. The report differed only in its generation date. |
| Step 07 CUDA checkpoint replay | Re-predicted all 1,004 drugs across 978 landmarks using the archived final XPert checkpoint. The full raw prediction CSV matched byte for byte. |
| Step 05 CUDA benchmark replay | Re-evaluated all 20 finetuned folds and 20 zero-shot fold anchors from archived weights. All 40 fold CSVs, the combined 40-row benchmark CSV and the benchmark PNG matched byte for byte. |
| Step 06 one-epoch CUDA training | Reproduced the original opening objective receipt at the logged three-decimal precision. Total 1.059, terms 0.677, -999.000 (disabled), 0.638 and 0.886. The original four-epoch checkpoint was restored after the test. |

The CUDA replay and short fit ran serially in a scope capped at 10 GiB of WSL
RAM. Its peak process RSS was 2,194,260 KiB. The CPU verification scope was
capped at 8 GiB and peaked at 6,238,124 KiB RSS. The initial prediction
reached inference but failed when the vendor code loaded gene labels from its
old metadata location. The Stage 1 prediction adapter supplies the declared
metadata path without changing the vendor code. The XPert and Annotation
offline verifiers and the repository path and tracked-code checks passed.

The regenerated finetune delta report changed only the displayed schema
reference path. The original report was restored with the other sealed
products. All 286 entries in the pre-cutover checksum inventory then passed
in the promoted worktree.
