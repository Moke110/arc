# Stage 2 ECFP4 contract migration

[Specification #252](https://github.com/Moke110/arc/issues/252) and
[ticket #253](https://github.com/Moke110/arc/issues/253) establish ECFP4 as the
name of the existing Stage 2 representation. The glossary defines its
parameters. The mapping applies only to Stage 2.

| Sealed field or identity | Canonical successor representation |
|---|---|
| `morgan.npy` (02 and ADR 0078 03) | `ecfp4.npy` with the exact same file bytes and `parents.csv` order |
| `config.json:morgan` | `manifest.json:ecfp4` with radius, bits, binary, atom invariants, chirality, and RDKit version |
| `includeChirality` (02 input configuration) | `include_chirality` |
| `MLP_Morgan` (02 label) | `ECFP4_MLP` in current display and new comparisons |
| `mlp-morgan-project-f32-lr1e4-v1` (02 identity) | preserved in 02 replay and its receipts |
| RDKit `GetMorganGenerator` | unchanged API name |

From the repository root, create a canonical fingerprint candidate for each 02
or 03 source input directory with an existing `manifest.json`:

```bash
python stage2/01_bases/output/api/ecfp4.py --inputs <source-freeze> --candidate <new-directory>
```

The Stage 2 ECFP4 manifest binds the complete source manifest digest, original
schema, sorted JSON representation metadata, and the two file digests.
`open_ecfp4(directory, source_manifest=...)` verifies the binding to the
original freeze when its manifest path is supplied. The reader verifies the
version-specific paths and fields and
returns a memory-mapped binary array in that directory's parent order. Existing
02 and 03 replay entry points still consume their original freezes. The
fingerprint candidate is a companion contract for new comparisons and it
carries no expression, split or checkpoint data. The older full-axis 03
historical generation still replays through its original runtime.

## Verification

On 2026-11-20, the companion contracts were created from the real 02 tenv and
ADR 0078 tenv and Tahoe freezes. RDKit 2026.03.2 independently regenerated every
bit from each parent SMILES. All bits and parent rows matched, including the
shared set of 459 tenv parents and the distinct 352 Tahoe parents.

| Freeze | Source manifest SHA-256 | ECFP4 file SHA-256 |
|---|---|---|
| 02 tenv | `d0f0101ca112e5d3b6325db6e80e1a3c90c42332ac46c8211e10e1cc34ed5ade` | `cc60144180a0e0a13486ce162ab0a5c65fbdb83ea75a88ba19ce187aa516ded7` |
| 03 tenv | `e08ea61360eab0b3b0d1a4762652d9070acda3cbe00347285c1a62e6266a427e` | `cc60144180a0e0a13486ce162ab0a5c65fbdb83ea75a88ba19ce187aa516ded7` |
| 03 Tahoe | `aace570ca8d066cebc7152ba989b7653f1aa794d1b39bbe96014834524082b1d` | `f61b3979ed6cb234b7af6a19bf48008cb873029c63dcdba2359e36278be4bf3e` |

The [02 display figures](../../02_models_tenv-5k/figure/ecfp4/) are a
CSV-only successor that displays ECFP4_MLP. The existing 02 plotting CSVs, figure
PNGs and renderer remain at their recorded paths. A second run with only the
three sealed plotting CSVs in its input directory reproduced each new PNG byte
for byte. Visual inspection checked the canonical method label in the metrics,
paired comparisons and learning curves.

| 02 successor figure | SHA-256 |
|---|---|
| metrics.png | `b05015b20d064be3e45507674b30dda048741a472fd38e2550d1fc01b4b59a2d` |
| paired_differences.png | `524587a5fb9a95e00042338085687c9900ae1f65a8fc28546ccaca499d702a69` |
| learning_curve.png | `26887c513bcfab9da25d14dbce5dc9a1d8fb23db1865cb5fe85bbe26df9bbd3e` |
