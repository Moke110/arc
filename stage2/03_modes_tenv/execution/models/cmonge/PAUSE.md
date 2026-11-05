# CMonge #50 pause handoff

The user requested a pause after all already-started programs finish. The existing
256-sample probe finished naturally, session 82155 exit 0; no new test, qualifier,
formal training, scoring or review was launched after the pause instruction.
Ticket #50 stays open.

Current committed/source behavior streams one validation group at a time with
fixed per-group seeds. This fixes the material memory defect that cached 138
full-expression baseline/target batches. The profile evaluated all 138 validation
groups after one training step. This expanded validation is **not** the final
accepted author recipe. The parent's accepted decision is the native author's
one randomly sampled validation condition per 100 iterations, preserving the
exact native Sinkhorn + Monge-gap objective and seeded per-run sampling. That
source change and its rerun remain TODO; the source was not edited while the
probe ran or after the user requested the pause.

The initial batch 512 exceeded the enforced 10 GiB scope memory limit. The user
journal identified `oom-kill`; the old scheduler had misclassified wrapper exit
−15. The parent fixed scope-journal classification. The profile's old state was
explicitly corrected with journal evidence, retaining its 20.708 s allowance,
setting `oom_retries=1` and `batch_scale=.5`. Its single permitted batch-256 retry
completed. No batch below 256 was attempted. See `formal_profile_state.json` and
`formal_profile_oom_journal.log`.

The successful retry produced finite 1401 × 18164 predictions with exact reload
agreement and nonzero adapter/baseline/dose checks. It measured 5,874,356 KiB
`time -v` peak RSS; the more conservative sampled process-tree peak was
5.704 GiB. CPU-only JAX used no model GPU memory. Its native training-plus-full-
validation cost was 620.330 s; inference 0.541 s. Cumulative scheduler training
allowance including the failed attempt was 641.038 s. The completed batch-4
qualification peak must not be used as a formal-batch resource estimate.
The full 138-group validation cost must not be extrapolated to the accepted
native sampled-condition recipe. No test score was computed.

On explicit resume, first implement native sampled validation in
`execution/script/cmonge_adapter.py`: keep an independent RNG seeded per run,
sample one held-out condition at each validation event, compute the same native
total objective and select the minimum author validation objective. Record the
sampled condition in history and update the fidelity/checkpoint policy text.
Then run the previously accepted artifact test under the repository memory gate:

```bash
systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=8G -- /usr/bin/time -v \
  python -m unittest discover -s stage2/03_modes_tenv/execution/tests -p test_cmonge.py -v
```

Source changes invalidate old completed identities. Generate a fresh qualification
output rather than trying to reuse the existing output directories:

```bash
python - <<'PY'
import json
from pathlib import Path
root = Path.cwd()
model = root / 'stage2/03_modes_tenv/execution/models/cmonge'
queue = json.loads((model / 'qualification_queue_v2.json').read_text())
queue['input_root'] = str(root / 'stage2/03_modes_tenv/output/input')
queue['output_root'] = str(model / 'qualification_native_output')
queue['jobs'][0]['command'] = [str(root / 'stage2/_envs/cmonge/bin/python'),
    str(root / 'stage2/03_modes_tenv/execution/script/cmonge_adapter.py'),
    '--config', str(model / 'qualification.json')]
(model / 'qualification_native_queue.json').write_text(json.dumps(queue, indent=2) + '\n')
PY
python stage2/03_modes_tenv/execution/script/schedule.py run \
  stage2/03_modes_tenv/execution/models/cmonge/qualification_native_queue.json
```

After the corrected native-path qualifier succeeds, measure 100 training steps
at the already justified 256-sample batch, with native validation timing. Preserve
the original 512-OOM evidence and the consumed smaller-batch retry; do not permit
an additional reduction below 256. This profiling step has not started. Its exact
new queue/identity must be generated after the source correction and retry guard
are finalized. Formal config remains 10000 iterations, initial author batch 512,
evaluation every 100 iterations and one evidenced smaller-batch retry. Final
admission, independent review and the nine-run formal matrix remain pending.

The replayable synthetic pipeline #43 is independently complete at commits
58ebdec / b156bd3. The CMonge adapter and first qualification are 2986ab7 /
a231599. `pause_record.json` and the final pause commit preserve the current
source, completed probe, source identities, timing, terminal statuses and pending
native-validation correction. Generated binary outputs stay ignored in the
three local `*_output` directories.
