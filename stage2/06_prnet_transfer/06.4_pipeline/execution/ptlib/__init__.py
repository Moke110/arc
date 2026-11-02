"""``ptlib`` — the PRnet formal-transfer pipeline library (ticket #4).

``stage2/06_prnet_transfer/06.4_pipeline/execution/`` is the **only** place the
study's training and evaluation code lives (ADR 0032; step ``AGENTS.md`` for
``06_prnet_transfer``).  Tickets #5–#11 invoke this package through the thin
CLIs beside it and record results in their own substeps; they never edit it.

Provenance (same-source copies, not imports — the frozen upstream tree and the
frozen 05/03.1 steps are never imported):

* ``ptlib/model.py`` — the PRnet ``PGM`` / ``PAdaptor`` / ``PEncoder`` /
  ``PDecoder`` re-implementation from
  ``stage2/05_prnet_pretrain/execution/train_prnet.py`` (which mirrors
  ``pert_model/PRnet/source/models/PRnet.py``).  **Same-source copy — do not
  diverge.**  Layer shapes are parameterised here so #10 can change capacity
  without editing the pipeline.
* ``ptlib/metrics.py`` — the direction-AUPRC primitives from
  ``stage2/03_models_cmbench/03.1_cold_drug_18k/execution/script/metrics.py``
  (``true_set_indices`` / ``direction_labels`` / ``group_metrics``), extended
  with a vectorised, tie-aware AP matrix so per-epoch validation is cheap.
* ``ptlib/data.py`` — ``load_view`` / ``resolve_control_mean`` /
  ``morgan_fingerprints`` patterns from the same 05 file, re-targeted at the
  frozen 06.3 input contract.

Design rules enforced here: model selection and early stopping monitor
validation **direction AUPRC@200 (max mode)** in both stages; validation and
inference scores use zero noise and are deterministic; prediction is
target-free, so only the evaluator reads test treated expression; every run
writes a manifest with data/checkpoint hashes and a stop reason.
"""

__all__ = [
    "checkpoint",
    "config",
    "data",
    "evaluate",
    "h5io",
    "manifest",
    "metrics",
    "model",
    "modelcard",
    "paths",
    "predict",
    "queue",
    "train",
]

PIPELINE_VERSION = "06.4.0"
