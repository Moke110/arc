# PerturbDiff qualification

Ticket #51, part of #40. Source, native objective and declared adaptations are in
`FIDELITY.md`. The environment lock captures the isolated local overlay over arc.
Native setup is separately owned by `stage2/01_setup/execution/perturbdiff/`.

`config.json` retains the documented native scratch200,000-step ceiling; `pilot.config.json` is
explicitly a12-step resource/whole-axis qualification, never a formal fit. The
source's1000 diffusion training steps and default100 DDIM reverse steps remain
unchanged. A full-gene artifact pilot requires scheduler ownership and an unknown
profile to run alone:

```bash
/home/chang/miniconda3/envs/arc/bin/python \
  stage2/03_modes_tenv/execution/script/schedule.py run \
  stage2/03_modes_tenv/execution/models/perturbdiff/pilot.queue.json
```

Admission is pending actual native setup, compact artifact tests and whole-axis
runtime evidence. Throughput must permit the original training recipe within the
four-hour budget:200,000 steps require at least13.89 steps/s before validation
and checkpoint overhead. Insufficient throughput is a measured resource gap; a
shortened pilot cannot establish formal feasibility. No test scoring is authorized
before the experiment's final persisted-prediction barrier.
