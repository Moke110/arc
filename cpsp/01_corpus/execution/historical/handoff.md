# Stage 3 first-wave implementation handoff

Validated 2026-11-05 on `dev`. Scope: specification #29 and implementation
tickets #30–#39. Model views, Stage 2 predictions, training and author contact
were outside this run.

## Delivered result

The current immutable local release is `stage3/01_data_corpus/output/corpus/v2/`.
Its public offline verifier passes: five archived source products, nine typed
condition/phenotype Parquet partitions, 6,096 original comparison units and 6,610
finite readout values across 6,620 measurement entries. The ten extra entries
preserve Sirenko raw blanks, including two supplied normalized zeros.

**All 6,096 comparisons are excluded; there are zero admitted samples and zero
validated fold changes.** This is a structurally valid `empty_after_QC` release,
not a training dataset or evidence that the compounds have no biological effects.
Class counts are −1: 0, 0: 0, +1: 0. Independent experimental-unit counts remain
unknown, not zero. No source/endpoint-expanded count establishes the later 500
model-view goal, which remains `not_run`; no model was trained.

| Source | Comparisons | Finite readouts | Main unresolved evidence |
| --- | ---: | ---: | --- |
| Sirenko 2013 | 524 | 1,038 | Parent/use identities, parent dose basis, exact normalization/control moments and numerical statistics |
| Heliste 2021 | 5,512 | 5,512 | CTG denominator, canonical identities, control/independent statistics; 96 h normoxia is isolated |
| Linn 2024 | 48 | 48 | Exact dose grid, matched vehicle and independent pointwise statistics |
| Liu 2021 | 8 | 8 | Administered vehicle matching; required H/R significance |
| Duan 2013 | 4 | 4 | SYA identity, vehicle matching and independent numerical statistics |

Source coverage, identities, overlapping exclusion reasons and original
raw/normalized distinctions remain available. No published summary was expanded
into invented wells. The three figure pilots retain coordinate precision and
independent visual review. Unresolved parent identities and uses remain explicit;
this work did not complete those missing scientific facts.

An earlier implementation admitted four Liu zero-class samples; final Spec review
revoked them because a zero-dose baseline did not confirm matched solvent.
`v1` was then retained unchanged and superseded by `v2` after the same vehicle
check corrected Duan's FC-computability report. The latest release changes the
curation interpretation without overwriting previous published bytes.

## Reproduction and evidence

Follow the step [README](../README.md) for acquisition and explicit curation setup.
After that, `release.py build --release-id v2` and `release.py verify --release-id
v2` consume frozen local evidence without network access. Identical rebuilds are
idempotent; conflicting payloads require a new release identity.

The complete source author logs, source-manifest digests, adapter hashes and
software versions are frozen in the release input and source registry. Rule
`fc-bh-v1` uses FC ≤0.8 / ≥1.2 and BH-adjusted p <0.05. Production endpoints are
nonnull int8 and readouts are native float64 lists with aligned evidence.

Tracked [release_record.json](release_record.json) records exact input/output
checksums, software and corpus-code identity. Human-facing source and partition
coverage are at the step root. Original downloads and full release payloads stay
local under the provenance-aware ignore rules; scripts, logs and small evidence
records are tracked. Liu and Duan archives record earlier acquired primary bytes
truthfully; Heliste's ancillary shRNA S2 timeout is recorded separately.

## Validation

- 21 public workflow/source/provenance/integration tests passed. Peak RSS:
  382,712 KiB; elapsed 26.03 s. Statistical fixtures cover all three endpoint
  labels even though production admits none.
- Root inventory and source corruption/HTML regressions passed. Missing assets,
  unknown items, missing manifests/verifiers and hash mismatches fail closed.
- `script/check_paths.sh` passed, including the Stage 2 retired-path audit;
  `script/check_tracked.sh` passed.
- Full `script/verify_all.sh` passed **32/32** items, preserving all 27 legacy
  sources. It ran under `MemoryMax=10G`, `MemorySwapMax=8G`; elapsed 4:50.77,
  measured peak process RSS 21,616 KiB. RSS does not measure the scope's page cache.
- The release build and independent offline verification passed. Source bytes
  and source manifests remained unchanged throughout curation and review fixes.

Exact command output is retained in `validation/`. Regression fixes caught a
mixed raw-p/author-BH subset adjustment, stale local provenance, duplicated source
cells, decimal dose-key jitter and unconfirmed vehicle matches. Temporary fixture
paths are constructed under their temporary repository rather than presented as
permanent inventory paths.

## Standards

The independent Standards review identified a missing aggregate source-tree
record. All five logs and acquisition/verification scripts now enforce file
count, bytes and the SHA256 of manifest bytes. Final recheck: **zero unresolved
findings**; see `validation/review-standards-final.md`.

## Spec

The independent Spec review identified unconfirmed vehicle matching and stale
taxonomy documentation. The affected comparisons are excluded and the glossary
now records five source categories and 32 items. Final cross-adapter recheck:
**zero unresolved findings**; see `validation/review-spec-final.md`.

All implementation commits are on `dev`; no implementation worktree was created,
and no branch was pushed. GitHub carries the final per-ticket evidence and the
explicit unresolved scientific limitations.
