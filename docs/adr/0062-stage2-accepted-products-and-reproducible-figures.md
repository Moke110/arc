# Stage 2 accepted products and independently reproducible figures

**Status:** accepted (2026-11-12).

Stage 2 steps separate reusable accepted products, report figures, execution code and
audit evidence. This replaces ADR 0030 contract two and its Stage 2 figure/log/checkpoint
placement, while preserving the external-item provenance contract. The authoritative
step layout and acceptance rules are in `stage2/AGENTS.md`.

Stable semantic paths replace accumulated build-version directories. Public readers and
metric implementations belong to the output contract; build and audit commands consume
that single implementation. Reusable checkpoints are outputs; transient checkpoints are
working intermediates. Figure panels retain complete plotted CSV values and render from
those values alone. This trades extra small plot tables for report reproduction without
large data access or repeating scientific calculations.

Formal publication requires structural/data checks, scientific acceptance, downstream
read compatibility and CSV-only figure reproduction. Record tolerances for stochastic
results separately from byte integrity. Construct and accept candidates before replacing
stable products; an unsuccessful build leaves the accepted release usable. Delete rejected
or superseded payloads only after resolving all accepted dependencies. Preserve original
sealed records, necessary source provenance, compact failure evidence and an explicit
old-to-new mapping rather than rewriting historical evidence. If dependency substitution
cannot be proved, rebuild and accept the affected result before deleting its old inputs.

The initial migration and comprehensive audit cover only `stage2/01_bases`. Stage 1 is
untouched. First inspect code, metadata, non-expression products and existing acceptance
records without a complete expression scan. Path/document/figure organization defects
can be repaired directly. Suspected identity, axis, matching, normalization, truth or
metric defects trigger bounded checks of affected expression inputs; confirmed errors
require rebuilding affected products and downstream results. Unaffected products remain
eligible for reuse. Scientific definitions in ADR 0061 and applicable evaluation ADRs
remain in force; new scientific trade-offs require a separate decision.

## See

- `stage2/AGENTS.md`
- `stage2/01_bases/AGENTS.md`
- `docs/adr/0030-external-data-provenance-and-step-layout.md`
- `docs/adr/0061-stage2-shared-bases-and-indexed-corpora.md`
