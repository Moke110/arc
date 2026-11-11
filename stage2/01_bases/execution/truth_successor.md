# Exact truth successor binding

`bind_truth_successor.py` creates a versioned envelope linking the accepted
`output/tenv-truth` statistics to `output/training-corpora`. It never changes
accepted truth bytes, reruns training, estimates prediction significance or
relaxes a numerical tolerance. The production population is exactly 18164
indexed genes, 4206 treated wells and 1377 normalized parent–dose contrasts.
Micro fixtures are explicitly marked and cannot be used by the production CLI or
reader.

Because the envelope records **absolute input paths** that
`bases.truth_successor.check_successor()` resolves and hashes at runtime, a path
move requires a re-issue, never an edit. `execution/reissue_truth_successor.py`
does that: it moves the superseded bytes to
`log/historical/seals/truth-successor-envelope-v1/`, creates and fully verifies a
candidate envelope, and only then swaps it into `output/tenv-truth-successor/`.
The new digest is recorded in `log/156-successor-reissue.json`.

The source truth must first have a complete
`stage2-parent-dose-truth-verification-v1` report. The report binds its truth
manifest and native/training/worker seals, covers all1377 contrasts, and records
three actual independent R refits. Representatives cover ordinary nonempty,
multi-plate and tested-empty contrasts whenever those classes exist.

## Proof and immutable inputs

The creator independently replays every original native count digest,
treated/control membership, BH value, test status and ranked list. It then
compares both native inputs as exact integers in the identical18164 order,
joining samples and controls by stable profile IDs rather than assuming raw
row positions. Every plate/design covariate and parent–dose membership must
agree. Different physical row orders and native file bytes are allowed only
when the actual mathematical data are equal.

The original and successor **tenv.h5ad** files must contain the same sample identities and
matched-control IDs. All X/control/delta values are compared with
`numpy.array_equal`, not `allclose`, and are also independently recomputed from
raw counts, full-axis CP10K/log1p and independently normalized matched-control
means. Row batches and per-contrast reads keep memory bounded. A one-ULP
float32 difference fails; do not widen a tolerance to make an envelope pass.
Even proportional integer-count changes that leave normalized expression
identical require a new fit and cannot be hidden by expression equivalence.

The envelope seals both training manifests/files, both native files and their
producer manifests, the accepted verification report, all original truth
payload hashes and the statistical recipe. It contains its own immutable copy
of the exact R worker. R/DESeq2 versions are retained from the original sealed
fit and compared to the current executable; its binary hash is recorded at
successor creation. This is not a retrospective claim that an R binary hash
was recorded before the original fit. Original statistics and directional
lists remain byte-identical references to v1, without copying or relabeling
old training provenance.

`SuccessorTruth(envelope, training, assignments, role)` validates these seals
and delegates to the partition-authorized reader. Only that role's directional
lists are opened. Full-release verification explicitly checks all statistics
and group payloads; the reader does not open unrelated partitions' truth
lists or the full statistic matrices. Ordinary training readers never import
this evaluator module. Consumers should record
`successor_manifest_sha256` and the source truth manifest seal alongside
prediction and cohort identities.

## Actual acceptance order

Complete 1377-contrast verification with three independent R refits and the
accepted exact successor envelope have passed; their actual records are
`log/records/truth_complete_verify.json` and
`log/records/truth_successor_acceptance.json`. The commands below document
that acceptance sequence. Reuse the sealed accepted evidence when inputs are
unchanged; they are not instructions to repeat an expensive fit or refit.
Obtain a resource slot for any necessary verification, keep the `arc` Python,
threads2, MemoryMax10G/MemorySwapMax8G and `/usr/bin/time -v`.

The original fit uses only accepted v1 **tenv**. V1 Tahoe failed independent
source replay and is rejected/superseded; it is never read by this tenv proof.
Both historical files and their original manifest remain unchanged. Current
training, HVG and diagnostics use `output/training-corpora`, and labels use the
re-issued `output/tenv-truth-successor` envelope.

The original complete verifier selects three distinct representative classes:

```bash
systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 /usr/bin/time -v /home/chang/miniconda3/envs/arc/bin/python /home/chang/arc/stage2/01_bases/execution/build_significance.py --native /home/chang/arc/stage2/01_bases/output/tenv-raw-counts/tenv.h5ad --training /home/chang/arc/stage2/01_bases/output/historical-training-release --output /home/chang/arc/stage2/01_bases/output/tenv-truth --verify --replay-groups 3 --record /home/chang/arc/stage2/01_bases/log/records/truth_complete_verify.json
```

Only after actual shared-expression/training-corpora build and source acceptance:

```bash
systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 /usr/bin/time -v /home/chang/miniconda3/envs/arc/bin/python /home/chang/arc/stage2/01_bases/execution/bind_truth_successor.py --truth /home/chang/arc/stage2/01_bases/output/tenv-truth --original-training /home/chang/arc/stage2/01_bases/output/historical-training-release --training /home/chang/arc/stage2/01_bases/output/training-corpora --original-native /home/chang/arc/stage2/01_bases/output/tenv-raw-counts/tenv.h5ad --native /home/chang/arc/stage2/01_bases/output/tenv-raw-counts/tenv.h5ad --accepted-verification /home/chang/arc/stage2/01_bases/log/records/truth_complete_verify.json --output /home/chang/arc/stage2/01_bases/output/tenv-truth-successor --record /home/chang/arc/stage2/01_bases/log/records/truth_successor_build.json
```

Independently replay the complete equivalence envelope:

```bash
systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 /usr/bin/time -v /home/chang/miniconda3/envs/arc/bin/python /home/chang/arc/stage2/01_bases/execution/bind_truth_successor.py --verify --training /home/chang/arc/stage2/01_bases/output/training-corpora --output /home/chang/arc/stage2/01_bases/output/tenv-truth-successor --record /home/chang/arc/stage2/01_bases/log/records/truth_successor_verify.json
```

A count, axis, control, contrast, recipe/runtime seal, identity or exact layer
mismatch fails closed and leaves no accepted output. Report the mismatch and
fix/rebuild the affected release or refit where necessary; never label a
partial checkpoint or successful fixture as complete production truth.
