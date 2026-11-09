# Fold 0 Tahoe exposure and shared-parent Δ audit

Run from the repository root after admission by the Step 06 scheduler. Use a
fresh output directory for an independent replay.

```bash
step=stage2/06_pdiff_transfer_optimize
python "$step/execution/script/source_audit.py" "$step/output/source-audit"
python "$step/execution/script/source_audit.py" "$step/output/source-audit-replay"
diff -r "$step/output/source-audit" "$step/output/source-audit-replay"
```

The command uses the frozen Step 05 tenv cold fold 0 split, the historical
Step 04 Tahoe fold 0 request, and the pinned Step 05 globally cold Tahoe
release. It verifies the file hashes of the consumed data, the common
18,164-gene order and the ordered Tahoe fitting membership after global
exclusion. It joins cell lines to both Tahoe releases by case ID from the
pinned historical `wells.csv`. It reads tenv and Tahoe Δ from the frozen
float32 arrays, one selected row at a time.

| File | Contents |
| --- | --- |
| `exposure.csv` | Each Tahoe fitting role against each tenv role in both releases. Includes connectivity block identities and counts, and profile or well and group support for their intersections. |
| `delta_panel.csv` | Globally cold Tahoe fitting profiles paired with tenv training wells by shared parent-connectivity block. One row per tenv drug, Tahoe drug, Tahoe cell line and fitting role, tenv dose and bracket endpoint. |
| `exact_smiles_panel.csv` | The same panel restricted to identical standardized parent SMILES. |
| `lineage.json` | Input paths and SHA-256 digests, predecessor identities, gene axis, parent inventories and table checksums. |

The nine shared training-side parents define the descriptive panel. The
sensitivity panel covers seven parent-connectivity blocks and retains only
pairs with an exact standardized parent SMILES match. At 10 µM, inclusion in
the sensitivity panel means that the tenv structure matched a Tahoe fitting
structure at some observed Tahoe dose.

Each profile is the mean of the available Δ observations for the given drug,
dose and cell context in its recorded role. Amplitude is the root mean square
of the mean Δ across the frozen gene axis. Signed direction concordance is
the fraction of genes with matching positive or negative signs among genes
with nonzero mean Δ in both profiles. The table reports that denominator and
the numbers of contributing tenv wells and Tahoe profiles. A missing dose
leaves that side's amplitude blank. Concordance is blank when either dose or
joint nonzero gene support is missing.

The lower and upper observations are separate rows. For tenv 0.2 µM the
Tahoe endpoints are 0.05 and 0.5 µM. For tenv 1 µM they are 0.5 and 5 µM.
The tenv 10 µM rows record only tenv support and amplitude. The assays have
no dose in common, and the comparisons use the observed dose values.

Tenv measures bulk RNA sequencing of hiPSC-derived cardiomyocytes at 24 h.
Its frozen input projects raw counts to the CM gene axis before CP10K and
log1p normalization and subtracts the same-plate DMSO mean. Tahoe profiles
aggregate single-cell counts within each well and cell line. Its frozen
input projects the native-axis log1p CP10K layer to the CM axis and subtracts
the same-plate, same-cell-line DMSO_TF mean. Tahoe lacks per-sample exposure
time. Assay, normalization, dose, exposure time, cell background and
coverage vary together in this panel. Amplitude and sign agreement are
descriptive and cannot isolate their contributions.

Historical Tahoe training overlaps one tenv validation block. Globally cold
Tahoe fitting overlaps zero. The audit reports coverage and support for
those sparse strata. Model selection uses the separately frozen paired
validation comparisons.
