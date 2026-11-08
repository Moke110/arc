# Step 05 — PerturbDiff transfer from Tahoe to tenv

> **Cross-dataset dose encoding limitation.** Tahoe uses 0.05, 0.5 and 5 µM.
> Tenv uses 0.2, 1 and 10 µM. The current bridge interpolates the learned
> categorical embeddings on the log-dose axis and clamps outside Tahoe's
> range. Zero-shot Tahoe receives its 5 µM embedding for tenv 10 µM.
> Transfer starts from the same mapping and can update it with tenv data.
> Improving cross-dataset dose encoding is a priority for follow-up.

[ADR 0047](../../docs/adr/0047-perturbdiff-tahoe-transfer-to-tenv.md) defines the
warm-start, cold-cm-drug and cold-drug comparisons. Each tests tenv responses
and compares CM, Tahoe, Transfer and Mix with two deterministic references.

## Shared cold tenv comparison

Both cold regimes use the frozen Step 04 tenv fold 0. The comparison CLI
verifies all three membership lists, their parent-connectivity blocks and
case order, the matched controls and the ordered 18,164-gene axis. It verifies
the three indexed CM checkpoints and their existing held-out predictions
against the Step 04 recipe, validation replay and prediction seals. Scoring
reads held-out expression only after these checks. The two seedless
references and the three CM predictions are scored once. The replay checks
all eight metrics and both directional PR curves against the Step 04 record.

Run from the repository root using the `arc` interpreter. Submit the full
replay to the dedicated scheduling subagent as a CPU job under the repository
resource limits.

```bash
step04=stage2/04_pdiff_hyperparameter
step05=stage2/05_pdiff_transfer
python "$step05/execution/script/cold_tenv.py" \
  "$step04/execution/output/tenv-formal-ticket72" \
  "$step04/output/checkpoints/tenv" \
  "$step05/execution/output/cold-tenv-replay"
```

`comparison.json` binds the two regime names to the shared artifacts, with
CM evidence indexed by seed and a SHA-256 inventory of the output.
`split.json` gives the frozen input identity, predictor paths and hashes,
case and group counts, and each role's parent-connectivity blocks.
`predictor_request.json` is the verified Step 04 request and lists only
predictor inputs. `per_run.csv`, `per_group.csv` and `tenv/pr_curves.json`
use the Step 04 evaluator format, where the model name `perturbdiff` means
the CM arm for this comparison.

The [scheduled replay record](execution/cold_record/replay.json) confirms
exact agreement with the Step 04 fold 0 scores and directional PR curves.
The frozen train, validation and test roles contain 2,388, 417 and 1,401
cases, respectively. The CPU replay finished in 14.50 seconds and peaked
at 459,572 KiB RSS. The replay record binds the local comparison and raw
scheduler evidence by SHA-256.
