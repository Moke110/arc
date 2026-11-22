# tenv KPGT input and bounded fit acceptance

[ADR 0080](../../../../docs/adr/0080-tenv-kpgt-dose-round1.md) and [#274](https://github.com/Moke110/arc/issues/274) define this input generation. Verified on 2026-11-21 (Asia/Shanghai).

The [input manifest](input-manifest.json) binds the unchanged Step 04 tenv input manifest, the accepted Step 03 representation provenance, the parent roster and the exact KPGT vector bytes. Preparation checked 459 unique parent InChIKeys and canonical SMILES, 2304 finite float32 values per parent and exact axis-local split and control identity. The accepted split ledger includes the excluded strips. Preparation opened the tenv freeze and the three KPGT source artifacts independently of Tahoe. Peak preparation RSS in its isolated check was 190,844 KiB.

The [CUDA smoke receipt](acceptance.json) covers fold 0, seed 42, two epochs and up to eight rows per partition in each scenario. Both fits initialized the response model from scratch, updated weights and reproduced validation and test predictions exactly in fresh processes.

The full acceptance sequence took 38.27 seconds with peak RSS 1,828,604 KiB under a 10 GiB memory cap. The Step 05 fit regression exercised independent replay of a fresh native checkpoint. The unchanged Step 03 API reproduced the accepted formal KPGT fold 0 predictions and scores in a temporary copy (see the [test-replay receipt](legacy-step03-replay.json)). The parent alignment test checked ordering, rejection of duplicate or missing parent identities, canonical SMILES mismatches and non-finite vectors.

Run the preparation and smoke commands singly under the [scope wrapper](../../execution/capped.sh) and sample RAM and GPU every 10 seconds.

```bash
D=stage2/04_tsg_dose
P=stage2/03_tsg_drug_emb
"$D/execution/capped.sh" /usr/bin/time -v python "$D/execution/kpgt/prepare.py" \
  --tenv "$D/execution/checkpoints/round1-axis-local-inputs/tenv" \
  --representations "$P/output/representations" --candidate "$INPUTS"
"$D/execution/capped.sh" /usr/bin/time -v python "$D/execution/kpgt/smoke.py" \
  --inputs "$INPUTS" --source pert_model/TranSiGen/source \
  --destination "$SMOKE" --device cuda
```
