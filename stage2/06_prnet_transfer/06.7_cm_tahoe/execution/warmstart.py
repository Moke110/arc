"""Strict Tahoe warm-start loader and provenance for ticket #7 (06.7_cm_tahoe).

The study's only trainer is ``06.4_pipeline/execution/train.py``; this module does
**not** re-implement training or the state-dict load.  It wraps the pipeline's own
:mod:`ptlib.checkpoint` primitives:

* :func:`ptlib.checkpoint.validate_checkpoint` — payload contract, ``n_genes``
  gene-axis equality and structure compatibility, fail closed;
* :func:`ptlib.checkpoint.load_provenance` — build the run's model, prove every
  parameter/buffer key is present with an identical shape/dtype, strict-load, and
  prove the loaded tensors are bit-identical to the checkpoint.  A partial or
  silently-randomised load raises before training.

On top of that, ticket #7 adds the provenance **record** the run needs: the
checkpoint sha256, the recorded config digest, the layer/key list, the empty
missing/unexpected/ mismatch lists, the parameter count, and two ticket-#7 guards:

1. the checkpoint must be the matching seed's **stage ``pretrain``** checkpoint
   (path ``tahoe_seed<seed>/checkpoint/model.pt``, payload seed == matrix seed), so
   a fold can never warm-start the wrong seed's weights;
2. an anti-random-init canary compares the loaded weights with a fresh
   ``random_init(model, seed)``; equal weights mean the "pretrain" checkpoint was
   never trained and the candidate would be a silent random-init run, so it fails
   closed instead.

The formal matrix is gated on :func:`verify_seeds` accepting all three seed
checkpoints; failures are reported as :class:`WarmStartError`, never as a partial
"some seeds loaded" success.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Iterable, Mapping

import common  # noqa: E402  (module-local helpers)

sys.path.insert(0, str(common.PIPELINE_EXEC))

from ptlib import checkpoint as pt_checkpoint  # noqa: E402
from ptlib import config as pt_config  # noqa: E402
from ptlib import model as pt_model  # noqa: E402

import torch  # noqa: E402

__all__ = [
    "WarmStartError",
    "CHECKPOINT_RELPATH",
    "checkpoint_path",
    "state_dict_digest",
    "verify_checkpoint",
    "verify_seeds",
    "verify_manifest_init",
]


class WarmStartError(RuntimeError):
    """A warm-start checkpoint is missing, mismatched, or not a usable pretrain."""


#: canonical run-relative checkpoint location (``ptlib.checkpoint.CHECKPOINT_RELPATH``)
CHECKPOINT_RELPATH = pt_checkpoint.CHECKPOINT_RELPATH


def checkpoint_path(checkpoint_root: Path | str, seed: int) -> Path:
    """Ticket #5's formal checkpoint for one seed: ``tahoe_seed<seed>/checkpoint/model.pt``."""
    return Path(checkpoint_root) / f"tahoe_seed{int(seed)}" / CHECKPOINT_RELPATH


def state_dict_digest(state: Mapping[str, "torch.Tensor"]) -> str:
    """Format-independent digest of a state dict (keys, dtypes, shapes, raw bytes)."""
    import hashlib

    digest = hashlib.sha256()
    for key in sorted(state):
        tensor = state[key]
        digest.update(key.encode("utf-8"))
        digest.update(b"\0")
        digest.update(str(tensor.dtype).encode("utf-8"))
        digest.update(b"\0")
        digest.update(repr(tuple(tensor.shape)).encode("utf-8"))
        digest.update(b"\0")
        digest.update(tensor.detach().cpu().contiguous().numpy().tobytes())
        digest.update(b"\0")
    return digest.hexdigest()


def _structure_problems(
    checkpoint_struct: Mapping[str, Any], expected_struct: Mapping[str, Any]
) -> list[str]:
    ok, mismatches = pt_config.structure_compatible(checkpoint_struct, expected_struct)
    return [] if ok else list(mismatches)


def verify_checkpoint(
    seed: int,
    path: Path | str,
    structure: Mapping[str, Any],
    *,
    require_pretrain_stage: bool = True,
) -> dict[str, Any]:
    """Validate one seed's pretrain checkpoint and return its provenance record.

    Raises :class:`WarmStartError` on the first failed guard (missing file,
    foreign payload, gene-axis mismatch, structure mismatch, partial/unexpected
    load, non-identical weights, wrong stage/seed, or weights identical to a fresh
    random init).  Never returns a partial record.
    """
    seed = int(seed)
    path = Path(path)
    try:
        if not path.is_file():
            raise WarmStartError(f"warm-start checkpoint is missing: {path}")
        size_bytes = path.stat().st_size
        sha256 = common.sha256_file(path)
        payload = pt_checkpoint.read_checkpoint(path)
        checkpoint_struct = pt_checkpoint.validate_checkpoint(
            payload, int(structure["n_genes"]), structure
        )
    except pt_checkpoint.CheckpointError as exc:
        raise WarmStartError(f"seed {seed}: {exc}") from exc

    config = dict(payload.get("config") or {})
    stage = config.get("stage")
    if require_pretrain_stage and stage != "pretrain":
        raise WarmStartError(
            f"seed {seed}: checkpoint config.stage is {stage!r}; ticket #7 warm-starts only "
            "from a stage 'pretrain' checkpoint"
        )
    payload_seed = payload.get("seed")
    if payload_seed is None or int(payload_seed) != seed:
        raise WarmStartError(
            f"seed {seed}: checkpoint payload seed is {payload_seed!r}, expected {seed} "
            f"(path {path})"
        )

    # The decisive strict-load proof; the pipeline's own loader is the only loader.
    try:
        model = pt_model.build_model(structure)
        provenance = pt_checkpoint.load_provenance(model, payload)
    except pt_checkpoint.CheckpointError as exc:
        raise WarmStartError(f"seed {seed}: {exc}") from exc

    # Anti-random-init canary: a checkpoint that was never trained would be
    # bit-identical to a fresh random_init(model, seed) and must not be called
    # a warm start.  Record how many tensors differ, not just a boolean.
    canary = pt_model.build_model(structure)
    pt_model.random_init(canary, seed)
    loaded_state = {key: value.detach().cpu() for key, value in model.state_dict().items()}
    canary_state = canary.state_dict()
    identical = [key for key in loaded_state if torch.equal(loaded_state[key], canary_state[key])]
    n_tensors = len(loaded_state)
    if len(identical) == n_tensors:
        raise WarmStartError(
            f"seed {seed}: checkpoint weights are bit-identical to a fresh random init for "
            f"seed {seed}; this is not a trained pretrain checkpoint"
        )

    record = {
        "seed": seed,
        "path": str(path),
        "checkpoint_relpath": str(
            path.relative_to(path.parents[1]) if len(path.parents) > 1 else path.name
        ),
        "sha256": sha256,
        "size_bytes": int(size_bytes),
        "payload_model": payload.get("model"),
        "payload_seed": int(payload_seed),
        "seed_matches": True,
        "checkpoint_epoch": int(payload.get("epoch", -1)),
        "checkpoint_best_valid_auprc": float(payload.get("best_valid_auprc", 0.0)),
        "checkpoint_stage": stage,
        "checkpoint_config": config,
        "checkpoint_config_sha256": pt_checkpoint.config_digest(payload),
        "checkpoint_structure": dict(checkpoint_struct),
        "structure_problems": _structure_problems(checkpoint_struct, structure),
        "structure_matches_run": True,
        "gene_axis": {
            "checkpoint": int(payload["n_genes"]),
            "expected": int(structure["n_genes"]),
            "matches": True,
        },
        "parameter_count": int(pt_model.parameter_count(model)),
        "state_dict_sha256": state_dict_digest(loaded_state),
        "load_provenance": provenance,
        "anti_random_init": {
            "canary": "fresh random_init(model, seed) with the same structure and seed",
            "n_tensors": int(n_tensors),
            "n_tensors_identical_to_random_init": len(identical),
            "matches_fresh_random_init": len(identical) == n_tensors,
        },
    }
    if record["structure_problems"]:  # defensive: validate_checkpoint already raised
        raise WarmStartError(
            f"seed {seed}: structure mismatch after validation: {record['structure_problems']}"
        )
    return record


def _provenance_identity(record: Mapping[str, Any]) -> dict[str, Any]:
    """The stable part of a record (no timestamp), used as the queue's digest."""
    return {key: value for key, value in record.items() if key != "verified_utc"}


def verify_seeds(
    seeds: Iterable[int],
    checkpoint_root: Path | str,
    structure: Mapping[str, Any],
) -> dict[str, Any]:
    """Verify every requested seed checkpoint; fail closed if any seed is unusable."""
    seeds = [int(seed) for seed in seeds]
    if not seeds:
        raise WarmStartError("no seeds requested for the warm start")
    if len(set(seeds)) != len(seeds):
        raise WarmStartError(f"duplicate seeds requested: {seeds}")
    root = Path(checkpoint_root)
    records: dict[str, dict[str, Any]] = {}
    for seed in seeds:
        record = verify_checkpoint(seed, checkpoint_path(root, seed), structure)
        records[str(seed)] = record
    shas = [record["sha256"] for record in records.values()]
    if len(set(shas)) != len(shas):
        raise WarmStartError(
            "two seeds resolve to the same checkpoint bytes: " + str({s: r["sha256"] for s, r in records.items()})
        )
    body = {
        "checkpoint_root": str(root),
        "seeds": seeds,
        "checkpoints": records,
        "all_verified": True,
    }
    body["provenance_sha256"] = common.content_sha256(
        {
            "checkpoint_root": str(root),
            "seeds": seeds,
            "checkpoints": {
                seed: _provenance_identity(record) for seed, record in records.items()
            },
        }
    )
    return body


# --------------------------------------------------------------------------- #
# per-run gate: does a run manifest prove the warm start we required?
# --------------------------------------------------------------------------- #
def verify_manifest_init(
    init_record: Mapping[str, Any] | None,
    expected: Mapping[str, Any] | None,
) -> list[str]:
    """Every way a run's ``manifest.init`` disagrees with the verified checkpoint."""
    problems: list[str] = []
    if not isinstance(init_record, Mapping):
        return ["run manifest has no init record"]
    if not isinstance(expected, Mapping):
        return ["no verified checkpoint record for this run's seed"]
    if init_record.get("mode") != "pretrained":
        problems.append(
            f"init.mode is {init_record.get('mode')!r}; a ticket-#7 run must be 'pretrained'"
        )
    source = init_record.get("source_path")
    if source is None or Path(str(source)).resolve() != Path(str(expected["path"])).resolve():
        problems.append(
            f"init.source_path is {source!r}, expected the seed checkpoint {expected['path']!r}"
        )
    if init_record.get("checkpoint_sha256") != expected["sha256"]:
        problems.append("init.checkpoint_sha256 differs from the verified checkpoint sha256")
    if init_record.get("checkpoint_config_sha256") != expected["checkpoint_config_sha256"]:
        problems.append("init.checkpoint_config_sha256 differs from the verified config digest")
    if int(init_record.get("checkpoint_seed", -1)) != int(expected["seed"]):
        problems.append(
            f"init.checkpoint_seed is {init_record.get('checkpoint_seed')!r}, "
            f"expected {expected['seed']}"
        )
    gene_axis = init_record.get("gene_axis") or {}
    if (
        gene_axis.get("matches") is not True
        or int(gene_axis.get("checkpoint", -1)) != int(expected["gene_axis"]["checkpoint"])
        or int(gene_axis.get("expected", -1)) != int(expected["gene_axis"]["expected"])
    ):
        problems.append(f"init.gene_axis does not prove the shared gene axis: {gene_axis}")
    if init_record.get("structure_matches_run") is not True:
        problems.append("init.structure_matches_run is not true")
    provenance = init_record.get("load_provenance")
    if not isinstance(provenance, Mapping):
        problems.append("init.load_provenance is missing (the strict-load proof)")
    else:
        if provenance.get("strict_load") is not True:
            problems.append("init.load_provenance.strict_load is not true")
        if provenance.get("weights_identical") is not True:
            problems.append("init.load_provenance.weights_identical is not true")
        for key in (
            "missing_keys",
            "missing_parameter_keys",
            "missing_buffer_keys",
            "unexpected_keys",
            "shape_mismatches",
            "dtype_mismatches",
            "nonfinite_keys",
        ):
            value = provenance.get(key)
            if value != []:
                problems.append(f"init.load_provenance.{key} is {value!r}, expected []")
        if int(provenance.get("n_loaded_tensors", -1)) != int(expected["load_provenance"]["n_loaded_tensors"]):
            problems.append("init.load_provenance.n_loaded_tensors differs from the verified load")
        if provenance.get("loaded_layers") != expected["load_provenance"]["loaded_layers"]:
            problems.append("init.load_provenance.loaded_layers differ from the verified load")
    return problems
