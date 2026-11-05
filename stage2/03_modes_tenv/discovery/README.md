# Tenv candidate audit — issue #41

Frozen source-level audit for [spec #40](https://github.com/Moke110/arc/issues/40),
2026-11-05 (Asia/Shanghai). All 14 downloaded model items and six additional
candidates are accounted for in `candidate_registry.json`. The old entrant
registry is historical evidence only. At this source-audit snapshot, no candidate had passed the study's
18,164-gene runtime gate; its `runtime_qualified_models` stays empty as frozen
research evidence. Subsequent runtime admission is recorded separately in
[`../execution/formal/registry.json`](../execution/formal/registry.json).

## Effective discovery and limitations

Biomni Layer 2 used its required interpreter
`/home/chang/miniconda3/envs/biomni_e1/bin/python`, source commit
`400c1f366b96a35ca253e13c9b06c5076af41d65`, and unmodified
`biomni.tool.literature.query_pubmed`. Each attempt was bounded by `timeout 45`.
The six raw logs preserve three sandbox DNS failures and their successful
network-enabled retries. `retrieval_record.json` binds query, result count,
completion time and log hash. No credentials were loaded or recorded.

The broad query returned five records, mostly tangential. The model-name query
returned eight, including the relevant PertDiff paper and a new scDPR paper;
PRnet acronym collisions also returned renal-tumour and unrelated computing
papers. The precise-title query still returned unrelated results. This is
successful remote retrieval with limited precision, not evidence that every
result is relevant or that the search is exhaustive. Biomni's formatter omits
PMIDs/DOIs, so relevant hits were separately checked against primary publisher
pages. The prior failed HTTP 500/429 and cancelled attempts remain recorded in
[ADR 0042](../../../docs/adr/0042-tenv-model-benchmark-redesign.md).

Primary-source checks recovered formal publication for
[PRnet](https://doi.org/10.1038/s41467-024-53457-1),
[PertDiff](https://doi.org/10.1093/bioinformatics/btag482),
[PerturbDiff](https://proceedings.mlr.press/v306/yuan26d.html),
[scDPR](https://doi.org/10.1016/j.neunet.2026.109005), and
[ExPO](https://doi.org/10.1186/s13321-026-01226-1).
The official scDFM repository cites ICLR 2026, but direct OpenReview access met a
browser challenge. The CellFlow README still links a preprint; a formal venue
has not been verified. Publisher direct fetches sometimes failed (OUP/Nature
redirects, Elsevier 403), while publisher-indexed search records were available;
these limitations do not establish model incompatibility.

## Decisions and per-model work list

`admitted-to-runtime-validation` means that the inspected source has a plausible
existing baseline, chemical and dose path. It does not mean admitted to formal
training. Instantiate one qualification ticket per entry, attach all as native
blockers of #46, and settle their outcomes before formal dispatch:

| model | qualification work |
|---|---|
| PRnet — preferred pilot | Native PGM, Gaussian NLL, Morgan 2048→1024 linear adapter and native log10(dose+1) scaling; prove baseline/adapter participation, full-axis checkpoint reload, host/GPU peak and author-validation checkpoint selection. |
| CPA | Native RDKit embedding/projection and basal encoder; select/disclose a native dose setting that does not require fitting held-out drug-specific parameters; preserve author objectives and validate full-axis behavior. |
| chemCPA | Existing continuous embedding input and amortized dose option; jointly train the shared-Morgan linear adapter with native predictor; retain adversarial/reconstruction training and validation fidelity. |
| CMonge | Native conditional control transport with continuous chemical context and log dose; place the linear adapter inside joint JAX training, avoid integer dose truncation, and test whole-axis transport/Sinkhorn feasibility without a prelearned reduced gene substitute. |
| PerturbDiff | Acquire pinned source with ADR 0030 provenance, isolated native setup, then replace the frozen chemical-vector input with the agreed joint linear adapter; retain native dose/control branches and diffusion objective; prove full-axis feasibility. |

`unresolved` entries remain visible and are not formal entrants. ExPO has native
basal expression and configurable gene dimension, so the old 978-only exclusion
is insufficient under the current contract. Its Morgan boundary, objective
fidelity and native setup need verification. scDFM's continuous perturbation
vector signature does not itself prove native chemical/dose conditioning;
the inspected run uses a categorical lookup. scDPR was found through Biomni,
but a public model/train/test repository could not be verified in bounded
model-name and exact-title searches. CellFlow lacks verified formal publication
and interface gates. These are named gaps, not architectural exclusions.

The remaining downloaded items are explicitly excluded with source pointers:
biolord relies on sample-index latent codes rather than encoding a new baseline;
DECANT's released mechanism forward requires treated bags; scVIDR lacks a
chemical cold-drug interface; WAVE and TranSiGen lack native dose paths; the
pinned State transition uses perturbation identity. Squidiff, PerturbNet and
XPert are user exclusions. PertDiff lacks an explicit dose input; its formal
publication does not override that failure. The newly noticed
[MAP](https://doi.org/10.1038/s42256-026-01286-w) pipeline requires independently
pretrained cell and molecular representations and fails the agreed initialization
policy. No resource failure is claimed for an unrun model.

## Reproduction and evidence boundary

The exact queries and bounded Python commands are in `retrieval_record.json`.
Replaying needs network access and can return different literature; the saved
logs are the frozen evidence. `source_evidence_hashes.json` binds locally
inspected files and read-only remote code observations. Local source versions
come from each item's `author_log.json`; remote commit IDs accompany the new
candidate rows. No upstream `source/` bytes, historical outputs, environment,
or source provenance were changed. Acquisition of a new runtime candidate must
create its own ADR 0030 item; temporary read-only audit snippets under `/tmp`
are not a substitute for vendoring or runtime reproducibility.

This ticket changes research records only; it loads no datasets and imports no
model packages. Validation checks registry coverage, unique models, accepted
statuses, local evidence hashes and retrieval content. Software TDD belongs to
the downstream artifact CLI/runtime qualification, at the agreed seam.
