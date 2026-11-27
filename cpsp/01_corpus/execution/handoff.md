# Heliste desalted-parent corpus — #54 / #55

The user requested desalted forms for every chemical, including the six previously
pending labels. Both conditions now have 2,744 admitted values and 12 pending values.
All 689 source labels and 5,512 numerical cells remain accounted for; 686 labels
resolve to 663 unique desalted parent structures. 264 previously admitted reference structures changed representation after
desalting; 23 excess label identities now share parent structures without row
merging. Source labels/cells, doses,
measurements, null replicates and condition mapping are unchanged, with no merging
of rows that now share a parent. Prior reference-form outputs remain in Git.

The immutable chemical_resolution_v1.json and original downloaded bytes remain
unchanged. The declared chemical_resolution_v2.json freezes their hash, original
reference records, the parent policy, extra raw PubChem evidence and deterministic
RDKit 2026.03.2 derivation. Ordinary drugs use the largest organic fragment and
removable-proton uncharging; stereochemistry, isotopes and permanent charges remain.
Eight active metal compounds preserve their full coordination components and
charges, removing only water and alkali counterions. This avoids replacing B12 or
platinum drugs with their isolated organic ligand. Bacitracin zinc is desalted only
after an independent exact Bacitracin parent query matches the complete original
ligand stereochemistry. No tautomer or stereoisomer selection is performed.

Three newly resolved source labels are Quinine hydrochloride dihydrate (candidate
hydrates converge), Methscopolamine (Pamine) (cation/bromide converge, retaining
permanent positive charge), and Vitamin D2 (independent Ergocalciferol evidence
identifies the ordinary parent rather than a mislabeled deuterated reference).
Pending labels are Verteporfin (Visudyne), Formoterol hemifumarate and Sevelamer HCl.
The former two retain competing stereochemical specifications; Sevelamer is a
polymer without a supported single finite structure. Salt removal cannot justify
choosing an arbitrary isomer or inventing a polymer structure.

Public build/verify and all three acceptance tests pass. The desalting test first
failed against reference-form outputs. Tests retain source fidelity, complete
accounting, CSV/Parquet equivalence, idempotence, rehashed tampering rejection and
now assert desalted Amiloride, permanent Methscopolamine charge, metal retention,
new admissions and continued Formoterol ambiguity. Heliste source integrity, path
and tracking gates are recorded in validation/desalt-*.log. Full verify_all.sh is
not run alongside the active Stage2 training worker, per repository resource rules.

Reproduce with the arc Python environment:

```bash
python cpsp/01_corpus/execution/corpus.py build
python cpsp/01_corpus/execution/corpus.py verify
python -m unittest discover -s cpsp/01_corpus/execution/tests -p test_corpus.py -v
bash pheno_data/Heliste-2021/script/verify.sh
bash script/check_paths.sh
bash script/check_tracked.sh
```

Canonical acquisition is separate: `python pheno_data/Heliste-2021/script/desalt_chemicals.py --acquire`.
It uses the immutable v1 snapshot and only targeted pending/parent reference queries.
Source CTG normalization denominator, drug persistence and physical-replicate
semantics remain unresolved. The corpus contains no inferred controls or statistics.
Measured peak RSS: build 223,580 KiB, offline verify 240,020 KiB, final tests
240,360 KiB. Independent integration review and tracker closure are recorded after
implementation.

Final independent review of `7d2ae7c...e288e49` completed with zero remaining
Standards or Spec findings; see `validation/desalt-review-final.md`. The reviewer
independently checked all 5,512 cells and 689 identity records. All implementation
commits are on dev; no push or implementation worktrees. Spec #54 and ticket #55
resolved with final coverage, validation and the full-sweep resource limitation.
