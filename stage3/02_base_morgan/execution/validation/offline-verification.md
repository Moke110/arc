# Repository offline checks — #64

Both release-time routing/tracking checks passed after the tracked release
metadata was added:

```bash
bash script/check_paths.sh
bash script/check_tracked.sh
```

The full `script/verify_all.sh` provenance sweep was attempted from the root in
the `arc` environment, two library threads and hidden CUDA devices. It ran inside
`arc-morgan-offline-verification.scope`, with `MemoryMax=10G` and
`MemorySwapMax=8G`, under `/usr/bin/time -v` and the production resource guardian.
The guardian sampled host RAM, GPU and processes every five seconds and could
terminate only its owned process group.

**The full sweep did not finish and is not a passing check.** Six verifiers had
passed when an external formal Stage 2 trainer appeared:
`stage2/03_modes_tenv/execution/script/perturbdiff_v2_adapter.py`, PID 1651473,
for `joint_v2_budget/perturbdiff_fold0_seed17/request-0.json`. The guardian stopped
only this task's checksum sweep with `external_long_training` (child exit -15).
The external trainer was left running. At the stopping snapshot, available RAM
was 44,196,880,384 bytes and GPU occupancy 1,433 MiB; sampled owned process-tree
RSS peaked at 12,529,664 bytes. That RSS excludes file cache and is not a kernel
peak. GNU time was interrupted with its process group and left no final peak
record. The systemd journal recorded 2 min 23.936 sec CPU time under the cap.

Completed items in the interrupted sweep: Annotation, DToxS-iPSC, GSE244740,
GSE247541, GSE262419 and LINCS-LDS-1593. Remaining items were not certified by
this attempt. Raw evidence remains in `/tmp/arc-morgan-verify-all-driver.log` and
`/tmp/arc-morgan-offline-verification/{resources.jsonl,resource_summary.json}`.

The repository execution rules require the broad sweep to run without another
heavy job. The other session entered formal training, so this run does not keep
retrying or interfere with it. Once that training is idle, the remaining gate is
a complete guarded rerun, for example under the same limits:

```bash
systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- /usr/bin/time -v bash script/verify_all.sh
```

Check host RAM/GPU/training visibility before that command and retain the live
guardian used by this run. No automatic rerun is queued. This unfinished gate is
recorded on #64 before closure under the repository's end-of-run tracker rule.

The model release itself has separate passing evidence: authenticated corpus
preparation, all 60 real fits, independent numerical/provenance audit, complete
public evaluation/verification with unchanged artifacts, and final Standards/Spec
reviews. See [formal-execution.md](formal-execution.md) and
[review-final.md](review-final.md).
