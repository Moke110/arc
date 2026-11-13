# Errata to ADR 0066: GSE270788 arm sizes and released-axis CP10K shares

**Status:** accepted (2026-11-13). Corrects two factual figures in ADR 0066 without changing
any design decision. ADR 0066 remains the contract and is read together with this record;
its text is not rewritten, because accepted records are historical.

## Why this record exists

Two numbers in ADR 0066 came from interim evidence. Both were corrected while rebuilding
`stage2/01_bases` (spec #185, tickets #186–#192). The decisions they supported are unchanged;
only the evidence-level figures move. The corrections were first recorded in
`stage2/01_bases/log/acceptance.md` §7 and are restated here so the ADR set does not carry
superseded numbers.

## Erratum 1 — `GSE270788` arm roster and shared controls

ADR 0066 describes the `GSE270788` NICM arm as "7 donors vs 6 Donor" and states that all
three `GSE270788` arms use "the same six Donor controls". The deposit contradicts both.

Measured directly from the item's `source/GSE270788_metadata.csv.gz`, grouped by the metadata
`sample` column (46,384 cells), the post-CM-gate arms are **Donor 7 / NICM 6 / ICM 5 / AMI
5**. The series-matrix `!Sample_characteristics_ch1` disease state agrees. The `orig.ident`
column is a single value `HF` across all arms and must not be used for grouping.

Consequences: the NICM arm is 6 case against **7** shared Donor controls, and the ICM and AMI
arms are 5 against 7. No donor is relabelled — all 23 donors are used. The released tables
record `n_case`/`n_control` = 6/7 (`nHF`), 5/7 (`ICM`), 5/7 (`AMI`) in the disease-deg
catalogs. The decision itself is unchanged: the six-category roster, the CM gates, the
pseudobulk/control design and the tiering all stand.

## Erratum 2 — released-axis CP10K shares

ADR 0066's *Corpus and view semantics* quotes "12.9% / 36.3%" of full-axis CP10K. The 12.9%
figure reproduces on the released `G_2000` axis; **36.3% does not** — it was measured on an
interim pre-rebuild axis.

Measured on the **released** axes (view manifests, tenv, all 4,478 profiles): `G_2000` =
**12.9299987%** and `G_5000` = **37.2063189%** all-rows pooled (treated-only 12.9566820% /
37.1590893%; `tahoe_bulk` 10.4600758% / 24.7819987%). The released values are authoritative
and are recorded in each view manifest's `basis` block; ADR 0066's figures are retained there
only as `basis.design_record`.

The decision is unchanged: views are normalized on the full-axis library total and then
subset. The measured denominator ratio is ≈2.88× (the ADR's interim-axis ratio was ≈2.81×),
which does not alter the one-normalization choice.

## What this record does not change

No design decision, contract, layout rule, category roster, axis rule, view semantics,
removal list or deferred scope in ADR 0066 changes. This is a factual errata companion; where
the two records disagree on a number, this record is authoritative and ADR 0066 is read with
the correction applied.

## See

- `docs/adr/0066-stage2-corpus-view-and-gene-axis.md` — the amended contract.
- `stage2/01_bases/log/acceptance.md` §7 — recorded errata with the measurement method.
- `disease_data/GSE270788/author_log.json` and the item's `source/` metadata.
- `stage2/01_bases/output/view-2k/manifest.json`, `stage2/01_bases/output/view-5k/manifest.json`
  — `basis` blocks with the measured shares and the `design_record`.
- Issues #186–#192 and spec #185 — the rebuild that produced the evidence.
