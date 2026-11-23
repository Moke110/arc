# tenv representation freeze acceptance

The 2026-11-21 CPU freeze covers all 459 accepted tenv parents and all 4,206
treated profiles across the existing three cold-drug folds. Each fold contains
275 train parents, 31 validation parents and 153 held-out parents. A fresh
process verified parent identities and structures, the exact historical
ECFP4 bit pattern, original-scale float32 values, test rows, feature widths
and all digests by reading through `RepresentationInputs`.

| Representation | Parents | Width |
|---|---:|---:|
| ECFP4C | 459 | 2048 |
| ECFP6C | 459 | 2048 |
| FCFP4C | 459 | 2048 |
| FCFP6C | 459 | 2048 |
| E3FP | 459 | 1024 |
| KPGT | 459 | 2304 |
| UniMol | 459 | 512 |

The [manifest](representation-freeze.json) records feature and structure
digests, provenance, runtime versions and per-fold partition coverage. Its
SHA-256 is `bc10c6af897b1b0d939a685aa198e03b22a3392e3891b721f70cc49f1f028678`.
The bounded KPGT smoke passed first for three actual parents and produced
three vectors of 2304 finite values. A separate check of four parents
reproduced the frozen KPGT settings bitwise across two fresh processes.
An independent full UniMol extraction produced an identical matrix digest.

E3FP used the standard adaptive conformer pool, RMSD cutoff 0.5 Å, UFF and
seed 42. Its 1024-bit fingerprints used level 5, radius multiplier 1.718,
stereochemistry and the first three available conformers. The per-parent
ledger records 78 parents with one conformer, 49 with two and 332 with
three. UniMol's official molecular CLS readout used the v1 all-hydrogen
weights and generated 3D conformers for the complete roster.

The capped extraction took 10 minutes 12.16 seconds. GNU time recorded
1,415,020 KiB maximum resident set size and exit status zero. See the
[execution guide](../AGENTS.md) for the reconstruction command.
