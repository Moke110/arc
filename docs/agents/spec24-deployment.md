# Spec 24 installed release acceptance

The materialized checkout passed acceptance on 2026-11-27 (Asia/Shanghai).
Local `main` fast-forwarded from
`9e7c85e38d2ab09c9b509ff8b2f30c67ea07ea85` to the independently reviewed
`dev` revision `8f4d890af988f44b7539c89cb92adb7fa16409ed`. The public stage
roots at that revision are `cpst/` and `cpsp/`.
[Spec #384](https://github.com/Moke110/carps/issues/384) ·
[Acceptance ticket #390](https://github.com/Moke110/carps/issues/390) ·
[ADR 0106](../adr/0106-cpst-and-cpsp-public-release-boundaries.md)

## Cutover and public checks

The cutover preserved both original stage directories by same-filesystem rename.
The local Stage 1 directory returned to its original inode after cutover and
after the public checks with Stage 1 temporarily parked. Its 535-entry inventory,
two local checksum manifests, reference manifest, checkpoint links and
Annotation audit seed passed. The only pre-existing working-tree change was the
unstaged deletion of `tableA_cardiomyocyte_hypoxia_compound_datasets.csv`.
Its binary diff SHA-256 was
`39319b06914d61468614655882f30fe5699f7d33f9736acafee8c32177d1c6d8`
before and after acceptance.

Both layouts passed the repository path check, CPST and CPSP path audits,
tracked-executable check, and offline verification of all 66 external items.
The path, tracking and provenance outputs matched byte for byte across layouts.
There were zero unresolved path literals or tracking holes. The scoped
predecessor check skipped Stage 1 when absent and passed after restoration.
Canonical API imports, six CLI help checks and the root-only ignore fixture
passed in both layouts. The document audit verified 72 documents, 72 distinct
authors and 1,368 local links in each layout.
[Document author record](spec24-documentation-authors.md)

## CarPS-Transcriptome

The installed migration verified 1,669 metadata files against their pinned
originals. The sealed archive SHA-256 is
`edf4a127f16d55e216efacb6b0c77e2bbb06ab9b365eeb27254a5326772b5875`
and the correspondence record SHA-256 is
`87b32cf9fe3ba723d9a0f16acb8247c8a75696765738d872f41a0df101c50d76`.
The 8,996-path inventory comprised 8,053 files retaining their original inodes
and sizes, 92 preserved symlink targets and 851 reissued ignored metadata files
with verified digests. It matched after installation and before and after
scientific verification.
[Metadata correspondence](../../cpst/util/log/385-release-migration.json)

Both receipt ledgers passed, covering 11 TranSiGen and 25 CPST receipt entries.
The active tenv and Tahoe input checks passed. The alpha verifier replayed
all nine coordinates exactly. The shared verifier reproduced all 45 registered
PNGs in two independent CSV-only renders. Three separate CUDA replays used
the sealed Step 08 UniMol, Step 09 HVG-DEG 2k and Step 10 Zmol-FiLM fits at
fold 0 and seed 17. All three reproduced validation predictions, the validation
objective and held-out predictions exactly.
[CPST verification inventory](../../cpst/util/training/cpst-verification.md) ·
[Figure verification](../../cpst/AGENTS.md#verify-the-current-figures)

| Installed report | Coordinates | Scores | Condition rows checked | Receipt SHA-256 |
| --- | ---: | ---: | ---: | --- |
| Step 08 | 63 | 567 | | `feb3d9a0a4bddf3f66258ec8bc99c91b9f68de14671dd073325a0a09c0dc462c` |
| Step 09 | 72 | 720 | 33,192 | `0e2694ecfa30b8f2743d4055f219f031f9ca15b2832852fcc7ec1fcbd13c304b` |
| Step 10 | 90 | 900 | 13,830 | `3549fafc89188d7d626d8890543e4433481e3464203ac04249ca83a2c0d35eee` |

## CarPS-Phenotype

The installed baseline and H-R corpora verified 2,744 admitted measurements per
arm. The original and migrated inventories each verified 541 entries,
including 60 symlink targets. The migrated inventory and release proof passed
before and after the benchmark, and the original inventory still matched after
acceptance. The sequential `evaluate --require-complete` and
`verify --require-complete` commands both verified all 60 formal runs. They
returned identical results from saved predictions and artifacts. The final
files matched the release record by SHA-256.
[Corpus manifest](../../cpsp/01_corpus/output/manifest.json) ·
[Morgan release record](../../cpsp/02_base_morgan/execution/release_record.json)

| Formal artifact | SHA-256 |
| --- | --- |
| `benchmark.json` | `ce057dccac6313117c9c345eb380854928b4e6155b98a15f3eb1a94918dab4b6` |
| `report.md` | `c01ba8b98f65bdfc65d6794ee7b9cdd3dc50578edb28bbd861e55f429aca2efe` |
| `summary.json` | `517fdc7034725f9a61da96cc8048eb09af6233cec002ca4e6de319d700b14362` |

## Execution measurements

Heavy checks ran sequentially with a 10 GiB memory cap and 8 GiB swap cap.
Scientific commands used the ARC Python environment with numerical thread
caps of two. Both public passes used system GNU sed 4.9 because ARC sed
exited 139 even for `--version`. The provenance scopes reached observed memory
peaks of 6,895,906,816 bytes with Stage 1 absent and 7,654,047,744 bytes after
restoration. Across 159 complete resource samples, the lowest available memory
was 20,580,980 KiB in WSL, 16,020,740 KiB on Windows and 14,233 MiB on the GPU.
The raw logs, inventories and recovery copies were retained locally.
[WSL resource limits](execution-environment.md)

| Verification | Elapsed | Peak process RSS (KiB) |
| --- | ---: | ---: |
| Offline provenance, Stage 1 absent | 5:15.40 | 79,116 |
| Offline provenance, Stage 1 restored | 5:10.03 | 79,024 |
| CPST alpha | 0:12.49 | 1,559,356 |
| CPST Step 08 report | 2:10.94 | 1,640,312 |
| CPST Step 09 report | 2:37.84 | 878,376 |
| CPST Step 10 report | 1:19.61 | 839,588 |
| CPST figure verification | 1:02.84 | 251,280 |
| CPSP full evaluate | 3:36.85 | 854,316 |
| CPSP full verify | 3:34.83 | 854,928 |
