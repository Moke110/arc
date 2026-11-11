# Migration inventory — 2026-11-11

Baseline: dev `00346d3`, ADR 0061, issue #146. This ticket only adds bases;
legacy paths, external source trees and historical evidence remain unchanged.
Rollback is removal/revert of the additive identity release; it requires no
restoration or retraining of prior experiments.

## Independent Step09 work

The coordinator checked #144/#140 and obtained the originating chat's migration
idle-window confirmation at 17:34 Asia/Shanghai: three completed runs, five
cancelled runs, no training process. These issues remain independent and open.
Preserve `cancelled-stop-verification.json` and history commit `25e6a41`.
The Step09 `execution/tmp/rerun/formal/`, `execution/rerun_record/`,
`output/corpus/rerun/` G512/G18164, `output/round1/`, `output/retrieval/`, and
`execution/tmp/rerun/superseded-six-model/` are preservation boundaries.
Any later move requires old/new size/hash binding verification, checkpoint
readability and preserved cancellation markers before publication.

## Shared assets and evidence

- `02_dataset_corpus/output/corpus/v2/` contains the native raw h5ads and source
  tables. Native Tahoe is 25,805,562,874 bytes. Preserve the native axes and seals.
- `02_dataset_corpus/output/corpus/shared-cm-v1/` is already compatible projected
  CP10K/log1p expression. Bounded builder/replayer: `cmdata/shared_cm.py`.
  Tahoe has 67,018 profiles (65,618 treated/1,400 control); tenv has 4,478
  nonblank profiles (4,206 treated/272 control). Existing manifests bind the
  native counts and mapping/float32 expression payloads.
- Identity pins and real source-specific FAM153B checks are recorded in
  `output/identity-v1/manifest.json`. The existing crosswalk and two native
  mapping tables preserve identifiers/evidence rather than flattening ENSG IDs.
- Existing evaluator depends on Step09's raw-count/significance worker and
  Step03's sealed wells/groups. Relocate dependencies before retiring code.
  Its strict inventory path checker rejects symlinks. Historical Git source
  paths and JSON seal bytes must survive in migration records unchanged.
- Step03 frozen `output/input/`, Step09 corpus/round1/rerun seals, evaluator
  inventory and selected checkpoint bindings describe original training inputs;
  new normalization/significance provenance does not relabel those inputs.
- Step07 exclusively owned assets may be removed only after consumers of shared
  expression/evaluator code are identified and relocated; this ticket removes none.

## Resource preflight

Read-only exploration measured 60 GiB disk free (95% used), replacing the old
207 GB snapshot for scheduling. Approximate tree sizes: setup 6.4 GB; corpus
34 GB; Steps03/04/05/06 75/59/26/15 GB; Step07 554 MB; Step09 35 GB;
evaluation 1.4 GB. Worktrees carry tracked files only; large ignored assets are
read through explicit paths. The coordinator owns 10-second GPU/RAM/WSL
monitoring and heavy-job admission. This ticket's standard-library identity
build reads small metadata only and does not load expression or train models.
