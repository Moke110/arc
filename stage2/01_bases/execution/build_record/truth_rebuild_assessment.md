# Bases-only tenv truth input assessment

The user deleted other Stage 2 trees during #148 and then limited work to
`01_bases`. This branch preserves new bases code only. It does not restore any
deleted step, historical output, native corpus or independent experiment.
No R process or actual significance refit was started. Three public contract
tests passed; real release acceptance remains outstanding.

## Required data and existing sources

DESeq2 needs raw integer counts on exactly the frozen 18164-gene axis, not
`tenv.h5ad`'s CP10K/log1p expression, controls or delta. The deleted ignored
native file was 4510 × 19613 (4206 treated,272 controls,32 blanks), approximately
255MB compressed. It cannot be recovered from Git or inferred from normalized
expression.

The immutable external source **still exists** at
`pert_data/GSE262419/source/`:16 compressed plate CSVs,`GSE262419_hash.csv.gz`,
and the series matrix; total approximately80MiB. HGNC annotation
`pert_data/Annotation/source/hgnc_complete_set.txt` and compound evidence
`pert_data/Annotation/local/cm_compound_annotation.csv` also remain present.
The frozen indexed axis remains in bases. No external source was modified.

A bases-native raw reader/rebuilder can stream the plate files, map probes to
HGNC with the frozen annotation, retain raw well counts, and join structures
from preserved annotation. The original tracked source loader is readable in
Git at commit`aca92c2`, at
`stage2/02_dataset_corpus/execution/cmdata/sources.py:load_temp_env`; reuse must
be relocated into bases rather than restoring the deleted step. The old
loader assembles multiple dense matrices, so it should not be invoked unchanged
under the bounded-access rule. A plate/chunk writer avoids a full
observation-by-axis matrix.

No rebuild was attempted. Data size suggests a few hundred MB of new compressed
raw corpus and modest sequential CPU work rather than a Tahoe-atlas scan;
actual runtime/RSS must be measured. A whole native int32 matrix would be about
338MiB before intermediate allocations, but this is a size estimate, not
permission to materialize it. Keep the10GiB cap and process one plate/chunk at
a time. The later1377 contrast DESeq2 fits are the dominant cost:an old50-contrast
batch took about652s, so complete refitting can take hours. Checkpointed output
plus final statistic files needs approximately2GiB.

## Git recoverability and cohort boundary

Tracked code, the original axis/crosswalk, Step03 well/fold CSVs, the R worker
and evaluator inventory remain recoverable from Git history **if relevant to
a later authorized task**. The worker has already been copied into this new
bases branch. Ignored native h5ads, expression memmaps, rerun case caches,
checkpoints and prediction arrays are not recoverable merely by Git checkout.
Their deletion must not be hidden by fabricating replacement historical inputs.

Historical partitions are **not a prerequisite for global parent-dose truth**.
A bases evaluator accepts an explicitly supplied cohort/role table and fails on
incomplete groups or cross-role repeats. Future train-only HVG and new training
need a defined bases population/partition manifest; a new deterministic split
can be designed in the new scope, but must not be represented as an old
experiment or used to trigger retraining. There is no need to restore deleted
consumer trees to provide this boundary.

## Preserved interface

`bases.significance.PartitionTruth(release, assignments, role)` opens only that
role's ranked group lists. Assignments may use stable sample IDs or original
source/case IDs. `truth(parent,dose,sample_ids)` returns the full indexed direction
array or`None` for untestable contrasts; tested-empty returns an all-zero array.
`ranked(parent,dose)` returns measured up/down lists ordered by adjusted P,
absolute log2FC and index. Training readers never import this module.

`build_significance.py` builds checkpointed axis-only raw-count DESeq2 truth;
`--verify` replays every raw count digest, group/control membership, BH value
and ranked list, plus representative independent fits. The CLI requires an explicit`--native` path to a **new bases-native raw corpus**;
it does not default to a deleted step. This preserved branch
is implementation work, not a completed true-label release.
