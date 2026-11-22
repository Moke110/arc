#!/usr/bin/env python3
"""Frozen, view-independent evaluator truth and the noise diagnostics (plan C3).

The truth product is what every view is scored against, so it may depend only on the
accepted full ``view-5k`` release, the integer corpus and the sealed Step 04 cold-drug-dose
folds -- never on a fitted model or on a derived view. It publishes:

* ``primary/fold-<f>/`` -- exactly the cold-drug-dose test profiles of the sealed splits,
  each treated profile's ``delta = full-view X - matched control`` and the matched control
  used for ``no-dose`` inference, plus byte-identical ``parents.csv``/``morgan.npy`` copied
  from the canonical full freeze so evaluation needs only this product and a model.
* ``condition_mean/fold-<f>/`` -- one row per **replicated** test condition: the pooled
  delta computed by the ADR 0082 ``merge`` recipe from corpus integers (sum the members'
  retained counts, normalize the sum on its own total; pool the union of their matched
  vehicle wells the same way), with the member profile ids and per-fold counts.
* ``diagnostics.json`` + ``noise_floor.csv`` -- per replicated condition of the whole
  source the delta noise floor and the within-condition repeat correlation, and per source
  the heuristic single-sample delta PCC ceiling.

Nothing in the product depends on a clock, a fitted model or a random draw, so building it
twice yields byte-identical output.
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from views import (BATCH, CONTROL_RULE, RECIPE, SIZE, SOURCES, HarnessError,  # noqa: E402
                   bind, condition_table, environment, merge_profile_id,
                   pooled_condition_blocks, pooled_condition_counts, pooled_control,
                   relative, repo_root, sha256, vehicle_group_sums, write_json)
from inputs import (read_sealed_splits, sealed_binding, sealed_reference)  # noqa: E402

SCHEMA = "tsg-view-truth-v1"
FOLD_SCHEMA = "tsg-view-truth-fold-v1"
DIAGNOSTICS_SCHEMA = "tsg-view-truth-diagnostics-v1"
SCENARIO = "cold-drug-dose"
HEURISTIC = "heuristic diagnostic"

PRIMARY_COLUMNS = ("profile_id", "cell_line", "drug_parent", "dose_um", "partition",
                   "fold", "scenario")
CONDITION_COLUMNS = ("profile_id", "cell_line", "drug_parent", "dose_um", "members",
                     "count", "partition")
DIAGNOSTIC_COLUMNS = ("source", "condition_key", "condition_index", "cell_line",
                      "drug_parent", "dose_um", "member_count", "n_pairs",
                      "noise_floor", "noise_floor_population", "rho_mean", "sqrt_rho",
                      "weight_profile", "fold0_partition", "fold1_partition",
                      "fold2_partition")

# Frozen released-population expectations: a released-input change must fail closed here
# rather than publish a silently different truth population (plan C3).
EXPECTED_TEST_COUNTS = {"tenv": (469, 468, 465), "tahoe_bulk": (6643, 6650, 7950)}
EXPECTED_CONDITIONS = {"tenv": 1377, "tahoe_bulk": 52777}
EXPECTED_REPLICATED = {"tenv": 1377, "tahoe_bulk": 9741}

# ADR 0082 Context reference values (the study's pre-registered sanity check).
ADR_REFERENCE = {
    "tenv": {"noise_floor_median": 0.343, "repeat_correlation_median": 0.121,
             "pcc_ceiling": "~0.38"},
    "tahoe_bulk": {"noise_floor_median": 0.070, "repeat_correlation_median": 0.207,
                   "pcc_ceiling": "0.43-0.52"},
}


# --------------------------------------------------------------------------- #
# diagnostic quantities (pure, gene-axis sized)
# --------------------------------------------------------------------------- #
def _members(deltas):
    deltas = np.asarray(deltas, dtype=np.float64)
    if deltas.ndim != 2:
        raise HarnessError("within-condition diagnostics need a members x genes block")
    if len(deltas) < 2:
        raise HarnessError("a replicated condition needs at least two members")
    return deltas


def noise_floor(deltas):
    """sqrt(mean over genes of the within-condition per-gene variance, ddof=1)."""
    return float(np.sqrt(np.mean(np.var(_members(deltas), axis=0, ddof=1))))


def noise_floor_population(deltas):
    """Population-variance variant (ddof=0); reproduces the ADR 0082 0.343/0.070 refs."""
    return float(np.sqrt(np.mean(np.var(_members(deltas), axis=0, ddof=0))))


def pair_correlations(deltas):
    """Every member pair's Pearson correlation over genes (NaN for a degenerate pair)."""
    deltas = _members(deltas)
    centered = deltas - deltas.mean(axis=1, keepdims=True)
    norms = np.sqrt(np.sum(centered * centered, axis=1))
    values = []
    for left in range(len(deltas)):
        for right in range(left + 1, len(deltas)):
            scale = float(norms[left] * norms[right])
            values.append(float(np.dot(centered[left], centered[right]) / scale)
                          if scale > 0 else float("nan"))
    return np.asarray(values, dtype=np.float64)


def repeat_correlation(deltas):
    """Mean pairwise Pearson correlation of the member delta vectors (genes = samples).

    NaN when every pair is degenerate (a member delta constant across genes), the same
    convention the Step 04 vector metrics use for a constant truth row.
    """
    values = pair_correlations(deltas)
    return float(np.mean(values)) if np.isfinite(values).any() else float("nan")


def pcc_ceiling(rhos, weights):
    """Profile-weighted mean of sqrt(rho_c) over the replicated conditions.

    A heuristic diagnostic, not a measurement-error model: ``rho_c`` may fall below zero,
    so the square root is taken on the clipped correlation; non-finite ``rho_c`` rows drop
    out of both the numerator and the weight.
    """
    rhos = np.asarray(rhos, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    if rhos.shape != weights.shape:
        raise HarnessError("the ceiling needs one weight per condition")
    keep = np.isfinite(rhos) & (weights > 0)
    if not keep.any():
        return float("nan")
    root = np.sqrt(np.clip(rhos[keep], 0.0, 1.0))
    return float(np.sum(weights[keep] * root) / np.sum(weights[keep]))


def condition_diagnostics(deltas):
    """Every published within-condition diagnostic for one replicated condition."""
    deltas = _members(deltas)
    return {
        "member_count": int(len(deltas)),
        "n_pairs": int(len(deltas) * (len(deltas) - 1) // 2),
        "noise_floor": noise_floor(deltas),
        "noise_floor_population": noise_floor_population(deltas),
        "rho_mean": repeat_correlation(deltas),
    }


def partition_key(condition):
    return f"{condition.cell_line}|{condition.drug_parent_inchikey}|{condition.dose_um:g}"


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _record(path):
    path = Path(path)
    return {"file": path.name, "sha256": sha256(path), "bytes": path.stat().st_size}


def _read_population(path):
    frame = pd.read_csv(path, dtype={"profile_id": str})
    frame["profile_id"] = frame["profile_id"].astype(str)
    return frame


def _fold_mapping(sealed_splits):
    """``fold -> {(parent, dose): partition}`` from the sealed cold-drug-dose splits."""
    block = sealed_splits[sealed_splits["scenario"] == SCENARIO]
    folds = sorted({int(fold) for fold in block["fold"]})
    mapping = {}
    for fold in folds:
        rows = block[block["fold"] == fold]
        table = {}
        for key, partition in zip(
                zip(rows["drug_parent"].astype(str), rows["dose_um"].astype(float)),
                rows["partition"].astype(str)):
            previous = table.get(key)
            if previous is not None and previous != partition:
                raise HarnessError(
                    f"{SCENARIO} fold {fold}: condition {key} sits in two partitions "
                    f"({previous}, {partition})")
            table[key] = partition
        mapping[int(fold)] = table
    return folds, mapping


def _condition_partitions(conditions, mapping, folds):
    """``fold -> partition per condition index`` for every condition of the source."""
    lookup = {}
    parents = conditions["drug_parent_inchikey"].astype(str).to_numpy()
    doses = conditions["dose_um"].astype(float).to_numpy()
    for fold in folds:
        table = mapping[fold]
        partitions = []
        for parent, dose in zip(parents, doses):
            partition = table.get((parent, float(dose)))
            if partition is None:
                raise HarnessError(
                    f"condition ({parent}, {dose}) is absent from the sealed "
                    f"{SCENARIO} fold {fold} mapping")
            partitions.append(partition)
        lookup[fold] = np.asarray(partitions, dtype=object)
    return lookup


# --------------------------------------------------------------------------- #
# builder
# --------------------------------------------------------------------------- #
class _Truth:
    def __init__(self, *, source, output_dir, repo_root_path=None, size=SIZE,
                 fixture=False, corpus_dir=None, axis_dir=None, accepted_dir=None,
                 sealed_dir=None, views_root=None, full_freeze_dir=None):
        from views import Context

        if source not in SOURCES:
            raise HarnessError(f"unknown source: {source}")
        self.source = source
        self.size = int(size)
        self.fixture = bool(fixture)
        self.root = repo_root(repo_root_path)
        self.modules = bind(self.root)
        self.output = Path(output_dir) / source
        self.views_root = Path(views_root) if views_root else (
            self.root / "stage2/tsg_view/output/views")
        self.context = Context(
            view="full", source=source, output_dir=self.output,
            staging_dir=self.output / ".staging", repo_root_path=self.root,
            size=self.size, fixture=self.fixture, corpus_dir=corpus_dir,
            axis_dir=axis_dir, accepted_dir=accepted_dir)
        self.genes = self.context.genes
        self.sealed = Path(sealed_dir) if sealed_dir else sealed_reference(
            self.root, source)
        self.sealed_binding = sealed_binding(self.sealed, self.root)
        self.full_freeze = self._full_freeze(full_freeze_dir)
        self.folds = []
        self.receipt = {}

    def _full_freeze(self, explicit):
        if explicit is not None:
            path = Path(explicit)
            if not (path / "parents.csv").is_file() or not (path / "morgan.npy").is_file():
                raise HarnessError(
                    f"--full-freeze does not carry parents.csv/morgan.npy: {path}")
            return path
        candidate = self.root / "stage2/tsg_view/execution/freezes/full" / self.source
        if (candidate / "parents.csv").is_file() and (candidate / "morgan.npy").is_file():
            return candidate
        return None

    # -- provenance ------------------------------------------------------- #
    def _fingerprint_binding(self):
        provider = self.full_freeze if self.full_freeze is not None else self.sealed
        parents, morgan = provider / "parents.csv", provider / "morgan.npy"
        if not parents.is_file() or not morgan.is_file():
            raise HarnessError(
                f"no parents.csv/morgan.npy under {provider}; point --full-freeze at the "
                "canonical full freeze")
        identical = None
        sealed_parents = self.sealed / "parents.csv"
        sealed_morgan = self.sealed / "morgan.npy"
        if sealed_parents.is_file() and sealed_morgan.is_file():
            identical = (sha256(parents) == sha256(sealed_parents)
                         and sha256(morgan) == sha256(sealed_morgan))
            if not identical:
                raise HarnessError(
                    "the full-freeze parents/morgan differ from the sealed Step 04 freeze")
        return {"dir": relative(provider, self.root),
                "parents_sha256": sha256(parents),
                "morgan_sha256": sha256(morgan),
                "byte_identical_to_sealed_freeze": identical}

    def _fingerprint_files(self):
        provider = self.full_freeze if self.full_freeze is not None else self.sealed
        return provider / "parents.csv", provider / "morgan.npy"

    # -- split-derived populations ---------------------------------------- #
    def _populations(self, sealed_population, sealed_splits, conditions, mapping):
        view_rows = np.asarray(sealed_population["view_row"], dtype=np.int64)
        if len(np.unique(view_rows)) != len(view_rows):
            raise HarnessError("the sealed population view_row mapping is not injective")
        partitions = _condition_partitions(conditions, mapping, self.folds)
        test_rows, condition_mean = {}, {}
        for fold in self.folds:
            block = sealed_splits[(sealed_splits["scenario"] == SCENARIO)
                                  & (sealed_splits["fold"] == fold)
                                  & (sealed_splits["partition"] == "test")]
            rows = np.asarray(block["row"], dtype=np.int64)
            profiles = block["profile_id"].astype(str).to_numpy()
            population_profiles = sealed_population.iloc[rows]["profile_id"].astype(
                str).to_numpy()
            if not np.array_equal(profiles, population_profiles):
                raise HarnessError(f"{self.source}: sealed splits/population rows disagree")
            test_rows[fold] = rows
            replicated = ((partitions[fold] == "test")
                          & (conditions["member_count"].to_numpy(dtype=np.int64) > 1))
            condition_mean[fold] = np.flatnonzero(replicated)
        return test_rows, condition_mean, partitions

    def _partition_purity(self, conditions, selected, sealed_splits, sealed_population,
                          folds):
        """Assert no condition-mean row mixes members from different partitions."""
        view_rows = np.asarray(sealed_population["view_row"], dtype=np.int64)
        split_of_view = np.full(self._view_row_count, -1, dtype=np.int64)
        split_of_view[view_rows] = np.arange(len(view_rows), dtype=np.int64)
        lookup = {}
        for fold in folds:
            block = sealed_splits[(sealed_splits["scenario"] == SCENARIO)
                                  & (sealed_splits["fold"] == fold)]
            array = np.empty(len(sealed_population), dtype=object)
            array[np.asarray(block["row"], dtype=np.int64)] = block["partition"].astype(
                str).to_numpy()
            lookup[fold] = array
        for fold in folds:
            for index in selected[fold]:
                condition = conditions.iloc[int(index)]
                members = np.asarray(condition["members"], dtype=np.int64)
                sealed_rows = split_of_view[members]
                if (sealed_rows < 0).any():
                    raise HarnessError(
                        f"condition {partition_key(condition)} has a member outside the "
                        "sealed population")
                partitions = set(lookup[fold][sealed_rows].tolist())
                if partitions != {"test"}:
                    raise HarnessError(
                        f"condition {partition_key(condition)} mixes members from "
                        f"partitions {sorted(partitions)} in fold {fold}")

    # -- the single accepted-view pass ------------------------------------ #
    def _view_pass(self, release, obs, sealed_population, test_rows, conditions):
        source = self.source
        reader = release.source(source)
        count = len(obs)
        view_rows = np.asarray(sealed_population["view_row"], dtype=np.int64)
        profile_ids = obs["profile_id"].astype(str).to_numpy()
        roles = obs["role"].astype(str).to_numpy()
        if view_rows.min() < 0 or view_rows.max() >= count:
            raise HarnessError("the sealed view_row mapping leaves the accepted view")
        if not np.array_equal(profile_ids[view_rows],
                              sealed_population["profile_id"].astype(str).to_numpy()):
            raise HarnessError(
                "the sealed population does not name the accepted full-view rows")
        if not np.all(roles[view_rows] == "treated"):
            raise HarnessError("the sealed population names a non-treated accepted view row")
        treated = np.flatnonzero(roles == "treated")
        if not np.array_equal(np.sort(view_rows), treated):
            raise HarnessError(
                "the sealed population is not the accepted view's treated row set")

        # primary: destination slot per fold, in sealed split-row order
        destination = np.full(count, -1, dtype=np.int64)
        fold_of = np.full(count, -1, dtype=np.int64)
        for fold in self.folds:
            positions = view_rows[test_rows[fold]]
            destination[positions] = np.arange(len(positions), dtype=np.int64)
            fold_of[positions] = int(fold)

        # diagnostics: every member of every replicated condition of the source
        replicated_positions = np.flatnonzero(
            conditions["member_count"].to_numpy(dtype=np.int64) > 1)
        slots = np.full(count, -1, dtype=np.int64)
        slot = 0
        for position in replicated_positions:
            for member in conditions.iloc[int(position)]["members"]:
                slots[int(member)] = slot
                slot += 1
        members = slot

        tables = {}
        for fold in self.folds:
            directory = self.output / "primary" / f"fold-{fold}"
            directory.mkdir(parents=True, exist_ok=True)
            rows = len(test_rows[fold])
            tables[fold] = {
                "delta": np.lib.format.open_memmap(
                    directory / "delta.npy", mode="w+", dtype=np.float32,
                    shape=(rows, self.size)),
                "control": np.lib.format.open_memmap(
                    directory / "control.npy", mode="w+", dtype=np.float32,
                    shape=(rows, self.size)),
                "rows": rows,
            }
        staging = self.output / ".staging"
        staging.mkdir(parents=True, exist_ok=True)
        diagnostic_path = staging / "condition-delta.npy"
        diagnostic_table = np.lib.format.open_memmap(
            diagnostic_path, mode="w+", dtype=np.float32, shape=(members, self.size))

        written = 0
        delta_identity = None
        for start in range(0, count, BATCH):
            stop = min(start + BATCH, count)
            x = np.asarray(reader.X(start, stop), dtype=np.float32)
            control = np.asarray(reader.control(start, stop), dtype=np.float32)
            if x.shape != (stop - start, self.size) or control.shape != (
                    stop - start, self.size):
                raise HarnessError("the accepted view returned the wrong window shape")
            delta = x - control
            if delta_identity is None:
                stored = np.asarray(reader.delta(start, stop), dtype=np.float32)
                delta_identity = bool(np.array_equal(delta, stored))
                if not delta_identity:
                    raise HarnessError(
                        "the accepted view's delta layer differs from X - control")
            local_fold = fold_of[start:stop]
            local_destination = destination[start:stop]
            keep = local_destination >= 0
            for fold in self.folds:
                selected = keep & (local_fold == fold)
                if not selected.any():
                    continue
                rows = np.flatnonzero(selected)
                target = local_destination[rows]
                tables[fold]["delta"][target] = delta[rows]
                tables[fold]["control"][target] = control[rows]
                written += int(len(rows))
            local_slots = slots[start:stop]
            diagnostic = np.flatnonzero(local_slots >= 0)
            if len(diagnostic):
                diagnostic_table[local_slots[diagnostic]] = delta[diagnostic]
        expected = sum(table["rows"] for table in tables.values())
        if written != expected:
            raise HarnessError(
                f"the primary pass wrote {written} of {expected} cold-drug-dose test rows")
        for table in tables.values():
            table["delta"].flush()
            table["control"].flush()
            del table["delta"]
            del table["control"]
        diagnostic_table.flush()
        return diagnostic_path, diagnostic_table, replicated_positions, delta_identity

    # -- diagnostics ------------------------------------------------------ #
    def _diagnostics(self, conditions, diagnostic_path, diagnostic_table,
                     replicated_positions, partitions, condition_mean):
        member_counts = conditions["member_count"].to_numpy(dtype=np.int64)[
            replicated_positions]
        noise = np.empty(len(replicated_positions), dtype=np.float64)
        noise_population = np.empty(len(replicated_positions), dtype=np.float64)
        rhos = np.empty(len(replicated_positions), dtype=np.float64)
        pair_rhos = []
        rows = []
        offset = 0
        for position, condition_index in enumerate(replicated_positions):
            condition = conditions.iloc[int(condition_index)]
            count = int(member_counts[position])
            block = np.asarray(diagnostic_table[offset:offset + count], dtype=np.float64)
            offset += count
            metrics = condition_diagnostics(block)
            pairs = pair_correlations(block)
            noise[position] = metrics["noise_floor"]
            noise_population[position] = metrics["noise_floor_population"]
            rhos[position] = metrics["rho_mean"]
            pair_rhos.extend(pairs.tolist())
            partitions_by_fold = [str(partitions[fold][condition_index])
                                  for fold in self.folds]
            rows.append({
                "source": self.source,
                "condition_key": partition_key(condition),
                "condition_index": int(condition["condition_index"]),
                "cell_line": str(condition["cell_line"]),
                "drug_parent": str(condition["drug_parent_inchikey"]),
                "dose_um": float(condition["dose_um"]),
                "member_count": count,
                "n_pairs": int(count * (count - 1) // 2),
                "noise_floor": metrics["noise_floor"],
                "noise_floor_population": metrics["noise_floor_population"],
                "rho_mean": metrics["rho_mean"],
                "sqrt_rho": (float(np.sqrt(np.clip(metrics["rho_mean"], 0.0, 1.0)))
                             if np.isfinite(metrics["rho_mean"]) else float("nan")),
                "weight_profile": count,
                "fold0_partition": partitions_by_fold[0] if len(
                    partitions_by_fold) > 0 else "",
                "fold1_partition": partitions_by_fold[1] if len(
                    partitions_by_fold) > 1 else "",
                "fold2_partition": partitions_by_fold[2] if len(
                    partitions_by_fold) > 2 else "",
            })
        if offset != len(diagnostic_table):
            raise HarnessError("the diagnostic table does not cover its member rows")
        frame = pd.DataFrame(rows, columns=list(DIAGNOSTIC_COLUMNS))
        path = self.output / "noise_floor.csv"
        frame.to_csv(path, index=False, lineterminator="\n")
        del diagnostic_table
        diagnostic_path.unlink(missing_ok=True)

        pair = np.asarray(pair_rhos, dtype=np.float64)
        finite = pair[np.isfinite(pair)]
        ceiling = pcc_ceiling(rhos, member_counts)
        reference = ADR_REFERENCE[self.source]
        fold_counts = {}
        for fold in self.folds:
            selected = np.flatnonzero(np.isin(
                replicated_positions,
                np.asarray(condition_mean[fold], dtype=np.int64)))
            profiles = int(member_counts[selected].sum()) if len(selected) else 0
            fold_counts[f"fold-{fold}"] = {"conditions": int(len(selected)),
                                           "profiles": profiles}
        histogram = {str(int(count)): int((member_counts == count).sum())
                     for count in sorted({int(value) for value in member_counts})}
        diagnostics = {
            "schema": DIAGNOSTICS_SCHEMA,
            "source": self.source,
            "scenario": SCENARIO,
            "label": HEURISTIC,
            "condition_population": {
                "all_conditions": int(len(conditions)),
                "replicated_conditions": int(len(replicated_positions)),
                "replicated_profiles": int(member_counts.sum()),
                "member_count_histogram": histogram,
                "note": ("the replicated-condition population is not a random subset of "
                         "the source: repeats concentrate on multi-plate, multi-article "
                         "conditions"),
            },
            "noise_floor": {
                "definition": ("sqrt(mean over genes of the within-condition per-gene "
                               "variance of Delta, ddof=1)"),
                "unit": "same as Delta RMSE",
                "median": float(np.median(noise)),
                "mean": float(np.mean(noise)),
                "conditions": int(len(noise)),
                "adr_reference_median": reference["noise_floor_median"],
                "deviation": float(np.median(noise) - reference["noise_floor_median"]),
            },
            "noise_floor_population": {
                "definition": ("sqrt(mean over genes of the within-condition per-gene "
                               "variance of Delta, ddof=0) = RMS deviation from the "
                               "condition mean; the variant that reproduces ADR 0082"),
                "unit": "same as Delta RMSE",
                "median": float(np.median(noise_population)),
                "mean": float(np.mean(noise_population)),
                "conditions": int(len(noise_population)),
                "adr_reference_median": reference["noise_floor_median"],
                "deviation": float(np.median(noise_population)
                                   - reference["noise_floor_median"]),
            },
            "repeat_correlation": {
                "definition": ("per condition: mean pairwise Pearson correlation of the "
                               "member Delta vectors over genes; per source: the median"),
                "rho_median_condition": float(np.nanmedian(rhos)),
                "rho_median_pair": float(np.median(finite)) if len(finite) else float("nan"),
                "pairs": int(len(pair)),
                "conditions_with_defined_rho": int(np.isfinite(rhos).sum()),
                "adr_reference_median": reference["repeat_correlation_median"],
                "deviation": (float(np.median(finite))
                              - reference["repeat_correlation_median"]
                              if len(finite) else float("nan")),
            },
            "pcc_ceiling": {
                "label": HEURISTIC,
                "definition": ("profile-weighted mean of sqrt(rho_c) over replicated "
                               "conditions"),
                "value": ceiling,
                "weighting": "member-profile count per replicated condition",
                "input_condition_population": "every replicated condition of the source",
                "input_conditions": int(len(replicated_positions)),
                "input_profiles": int(member_counts.sum()),
                "unweighted_mean_sqrt_rho": float(np.nanmean(
                    np.sqrt(np.clip(rhos, 0.0, 1.0)))),
                "adr_reference": reference["pcc_ceiling"],
                "note": ("heuristic diagnostic, not a measurement-error model; the input "
                         "conditions are not a random subset of the source"),
            },
            "replicated_test_conditions": fold_counts,
            "delta_source": {
                "view": "accepted full view-5k",
                "definition": "delta = X - matched control (mean of normalized vehicles)",
                "accepted_view_sha256": self.context.source_manifests()[
                    "accepted_view_sha256"],
                "sealed": self.sealed_binding,
            },
            "file": {"name": "noise_floor.csv", "sha256": sha256(path),
                     "rows": int(len(frame)), "bytes": path.stat().st_size},
        }
        diagnostics_path = self.output / "diagnostics.json"
        write_json(diagnostics_path, diagnostics)
        return diagnostics, diagnostics_path

    # -- condition-mean truth --------------------------------------------- #
    def _condition_mean(self, conditions, condition_mean, sealed_splits):
        selected = []
        fold_slots = {}
        for fold in self.folds:
            indices = [int(index) for index in condition_mean[fold]]
            fold_slots[fold] = (len(selected), len(selected) + len(indices))
            selected.extend(indices)
        subset = conditions.iloc[np.asarray(selected, dtype=np.int64)].reset_index(
            drop=True)
        if not len(subset):
            raise HarnessError(f"{self.source}: no replicated test condition to pool")

        row_to_slot = np.full(self._view_row_count, -1, dtype=np.int64)
        for slot, condition in enumerate(subset.itertuples()):
            for member in condition.members:
                row_to_slot[int(member)] = slot

        genes = self.genes
        with self.context.corpus(self.source) as reader:
            if not np.array_equal(reader.obs.index.astype(str).to_numpy(),
                                  self._profile_ids):
                raise HarnessError(
                    "the integer corpus row order differs from the accepted view row order")
            pooled = pooled_condition_counts(reader, genes, row_to_slot, subset)
            group_sums = vehicle_group_sums(reader, genes)
        control, _retained = pooled_control(subset, group_sums)

        expected = np.asarray(
            [int(condition.library_totals.sum(dtype=np.int64))
             for condition in subset.itertuples()], dtype=np.int64)
        treated = np.asarray(pooled.sum(axis=1), dtype=np.int64)
        if not np.array_equal(treated, expected):
            raise HarnessError(
                f"{self.source}: a pooled retained total differs from its member sum")

        tables = {}
        for fold in self.folds:
            directory = self.output / "condition_mean" / f"fold-{fold}"
            directory.mkdir(parents=True, exist_ok=True)
            left, right = fold_slots[fold]
            tables[fold] = np.lib.format.open_memmap(
                directory / "delta.npy", mode="w+", dtype=np.float32,
                shape=(right - left, self.size))
        written = 0
        for start, block in pooled_condition_blocks(pooled, control):
            stop = start + len(block["delta"])
            for fold in self.folds:
                left, right = fold_slots[fold]
                begin = max(start, left)
                overlap = max(0, min(stop, right) - begin)
                if overlap <= 0:
                    continue
                tables[fold][begin - left:begin - left + overlap] = block["delta"][
                    begin - start:begin - start + overlap]
                written += overlap
        if written != len(subset):
            raise HarnessError(
                f"the condition-mean pass wrote {written} of {len(subset)} rows")
        for table in tables.values():
            table.flush()
            del table

        merge_check = self._merge_cross_check(subset)
        records = {}
        for fold in self.folds:
            directory = self.output / "condition_mean" / f"fold-{fold}"
            left, right = fold_slots[fold]
            block = subset.iloc[left:right]
            population = pd.DataFrame({
                "profile_id": [
                    merge_profile_id(self.source, condition.cell_line,
                                     condition.drug_parent_inchikey, condition.dose_um)
                    for condition in block.itertuples()],
                "cell_line": block["cell_line"].astype(str).to_numpy(),
                "drug_parent": block["drug_parent_inchikey"].astype(str).to_numpy(),
                "dose_um": block["dose_um"].astype(float).to_numpy(),
                "members": [";".join(str(member) for member in condition.profile_ids)
                            for condition in block.itertuples()],
                "count": block["member_count"].to_numpy(dtype=np.int64),
                "partition": "test",
            }, columns=list(CONDITION_COLUMNS))
            population.to_csv(directory / "population.csv", index=False,
                              lineterminator="\n")
            outputs = {"population.csv": _record(directory / "population.csv"),
                       "delta.npy": _record(directory / "delta.npy")}
            split_block = sealed_splits[(sealed_splits["scenario"] == SCENARIO)
                                        & (sealed_splits["fold"] == fold)]
            manifest = {
                "schema": FOLD_SCHEMA,
                "kind": "condition_mean",
                "source": self.source,
                "scenario": SCENARIO,
                "fold": int(fold),
                "partition": "test",
                "conditions": int(len(population)),
                "profiles": int(population["count"].sum()),
                "columns": list(CONDITION_COLUMNS),
                "pooling_rule": RECIPE["merge"],
                "control_rule": CONTROL_RULE["merge"],
                "delta_definition": (
                    "pooled X - pooled control, float32 elementwise; pooled X = "
                    "log1p_cp10k of the summed member retained counts on their own "
                    "selected-axis total"),
                "coverage": {
                    "fold_test_profiles": int((split_block["partition"] == "test").sum()),
                    "replicated_test_conditions": int(len(population)),
                    "replicated_test_profiles": int(population["count"].sum()),
                },
                "partition_purity": ("every member of every listed condition is in the "
                                     "test partition of this fold"),
                "split_provenance": {**self.sealed_binding, "scenario": SCENARIO,
                                     "fold": int(fold)},
                "merge_cross_check": {
                    "available": bool(merge_check["available"]),
                    "equal": merge_check["equal"],
                    "merge_product_manifest_sha256":
                        merge_check["merge_product_manifest_sha256"],
                    "checked_conditions": (int(len(population))
                                           if merge_check["available"] else 0),
                },
                "outputs": outputs,
                "generated_by": "stage2/tsg_view/execution/truth.py",
            }
            write_json(directory / "manifest.json", manifest)
            records[f"fold-{fold}"] = {
                "conditions": int(len(population)),
                "profiles": int(population["count"].sum()),
                "outputs": outputs,
                "manifest": _record(directory / "manifest.json"),
            }
        return records, merge_check

    def _merge_cross_check(self, subset):
        manifest_path = self.views_root / "merge" / f"{self.source}.json"
        record = {"available": False, "equal": None,
                  "merge_product_manifest_sha256": None, "checked_conditions": 0}
        if not manifest_path.is_file():
            return record
        manifest = json.loads(manifest_path.read_text())
        payload = manifest.get("payload") or {}
        payload_path = self.views_root / "merge" / str(payload.get("file", ""))
        profile_path = self.views_root / "merge" / str(
            (manifest.get("profile_ids") or {}).get("file", ""))
        if not payload_path.is_file() or not profile_path.is_file():
            return record
        table = pd.read_csv(profile_path, dtype={"profile_id": str})
        lookup = {str(profile): int(row) for profile, row in
                  zip(table["profile_id"], table["row"])}
        identities = [merge_profile_id(self.source, condition.cell_line,
                                       condition.drug_parent_inchikey, condition.dose_um)
                      for condition in subset.itertuples()]
        missing = [identity for identity in identities if identity not in lookup]
        if missing:
            raise HarnessError(
                f"the merge product does not cover a condition-mean row: {missing[:3]}")
        positions = np.asarray([lookup[identity] for identity in identities],
                               dtype=np.int64)
        slot = {int(position): index for index, position in enumerate(positions)}
        merged = np.empty((len(subset), self.size), dtype=np.float32)
        import h5py

        with h5py.File(payload_path, "r") as handle:
            dataset = handle["layers/delta"]
            for start in range(0, len(dataset), BATCH):
                stop = min(start + BATCH, len(dataset))
                rows = np.arange(start, stop, dtype=np.int64)
                selected = np.asarray([int(row) for row in rows if int(row) in slot],
                                      dtype=np.int64)
                if len(selected):
                    block = np.asarray(dataset[start:stop], dtype=np.float32)
                    merged[[slot[int(row)] for row in selected]] = block[
                        selected - start]
        ours = np.concatenate([
            np.asarray(np.load(self.output / "condition_mean" / f"fold-{fold}"
                               / "delta.npy", mmap_mode="r"))
            for fold in self.folds], axis=0)
        if not np.array_equal(ours, merged):
            raise HarnessError(
                "the condition-mean delta does not match the accepted merge product")
        return {"available": True, "equal": True,
                "merge_product_manifest_sha256": sha256(manifest_path),
                "checked_conditions": int(len(subset)),
                "merge_product": relative(payload_path, self.root)}

    # -- assemble --------------------------------------------------------- #
    def build(self):
        source = self.source
        output = self.output
        shutil.rmtree(output, ignore_errors=True)
        output.mkdir(parents=True, exist_ok=True)
        sealed_population = _read_population(self.sealed / "population.csv")
        sealed_splits = read_sealed_splits(self.sealed)
        self.folds, mapping = _fold_mapping(sealed_splits)
        if not self.folds:
            raise HarnessError(f"{source}: the sealed {SCENARIO} splits carry no fold")
        fingerprint = self._fingerprint_binding()
        parents_file, morgan_file = self._fingerprint_files()

        with self.context.accepted_view() as release:
            reader = release.source(source)
            obs = pd.DataFrame(reader.obs).copy()
            obs["profile_id"] = obs.index.astype(str)
            conditions, _ = condition_table(obs, source)
            conditions = conditions.reset_index(drop=True)
            self._profile_ids = obs["profile_id"].astype(str).to_numpy()
            self._view_row_count = len(obs)
            if not self.fixture:
                if len(conditions) != EXPECTED_CONDITIONS[source]:
                    raise HarnessError(
                        f"{source}: the accepted view carries {len(conditions)} "
                        f"conditions, expected {EXPECTED_CONDITIONS[source]}")
            test_rows, condition_mean, partitions = self._populations(
                sealed_population, sealed_splits, conditions, mapping)
            observed = {fold: len(test_rows[fold]) for fold in self.folds}
            if not self.fixture:
                got = tuple(observed[fold] for fold in self.folds)
                expected = EXPECTED_TEST_COUNTS[source]
                if got != expected:
                    raise HarnessError(
                        f"{source}: the sealed {SCENARIO} test population is {got}, "
                        f"expected {expected}")
            replicated = int((conditions["member_count"].to_numpy(dtype=np.int64)
                              > 1).sum())
            if not self.fixture and replicated != EXPECTED_REPLICATED[source]:
                raise HarnessError(
                    f"{source}: {replicated} replicated conditions, expected "
                    f"{EXPECTED_REPLICATED[source]}")
            self._partition_purity(conditions, condition_mean, sealed_splits,
                                   sealed_population, self.folds)
            diagnostic_path, diagnostic_table, replicated_positions, delta_identity = (
                self._view_pass(release, obs, sealed_population, test_rows, conditions))
            diagnostics, diagnostics_path = self._diagnostics(
                conditions, diagnostic_path, diagnostic_table, replicated_positions,
                partitions, condition_mean)

        primary_records = {}
        for fold in self.folds:
            directory = output / "primary" / f"fold-{fold}"
            test = sealed_splits[(sealed_splits["scenario"] == SCENARIO)
                                 & (sealed_splits["fold"] == fold)
                                 & (sealed_splits["partition"] == "test")]
            rows = sealed_population.iloc[np.asarray(test["row"], dtype=np.int64)]
            population = pd.DataFrame({
                "profile_id": rows["profile_id"].astype(str).to_numpy(),
                "cell_line": rows["cell_line"].astype(str).to_numpy(),
                "drug_parent": rows["drug_parent_inchikey"].astype(str).to_numpy(),
                "dose_um": rows["dose_um"].to_numpy(dtype=float),
                "partition": "test",
                "fold": int(fold),
                "scenario": SCENARIO,
            }, columns=list(PRIMARY_COLUMNS))
            population.to_csv(directory / "population.csv", index=False,
                              lineterminator="\n")
            shutil.copyfile(parents_file, directory / "parents.csv")
            shutil.copyfile(morgan_file, directory / "morgan.npy")
            outputs = {name: _record(directory / name) for name in
                       ("population.csv", "delta.npy", "control.npy", "parents.csv",
                        "morgan.npy")}
            manifest = {
                "schema": FOLD_SCHEMA,
                "kind": "primary",
                "source": source,
                "scenario": SCENARIO,
                "fold": int(fold),
                "partition": "test",
                "rows": int(len(population)),
                "conditions": int(population[["cell_line", "drug_parent", "dose_um"]]
                                  .drop_duplicates().shape[0]),
                "columns": list(PRIMARY_COLUMNS),
                "delta_definition": "delta = full-view X - matched control, float32",
                "control_definition": CONTROL_RULE["full"],
                "delta_equals_view_delta_layer": bool(delta_identity),
                "split_provenance": {**self.sealed_binding, "scenario": SCENARIO,
                                     "fold": int(fold), "partition": "test"},
                "fingerprints": fingerprint,
                "accepted_view": self.context.source_manifests(),
                "axis_size": self.size,
                "outputs": outputs,
                "generated_by": "stage2/tsg_view/execution/truth.py",
            }
            write_json(directory / "manifest.json", manifest)
            primary_records[f"fold-{fold}"] = {
                "rows": int(len(population)),
                "conditions": int(manifest["conditions"]),
                "outputs": outputs,
                "manifest": _record(directory / "manifest.json"),
            }

        condition_records, merge_check = self._condition_mean(
            conditions, condition_mean, sealed_splits)

        counts = {
            "primary": {key: value["rows"] for key, value in primary_records.items()},
            "condition_mean": {
                key: {"conditions": value["conditions"], "profiles": value["profiles"]}
                for key, value in condition_records.items()},
            "diagnostics": diagnostics["condition_population"],
        }
        manifest = {
            "schema": SCHEMA,
            "source": source,
            "study": "stage2/tsg_view",
            "scenario": SCENARIO,
            "view_independent": True,
            "folds": [int(fold) for fold in self.folds],
            "primary": primary_records,
            "condition_mean": condition_records,
            "counts": counts,
            "diagnostics": {
                "file": "diagnostics.json",
                "sha256": sha256(diagnostics_path),
                "noise_floor_median": diagnostics["noise_floor"]["median"],
                "noise_floor_population_median":
                    diagnostics["noise_floor_population"]["median"],
                "repeat_correlation_median":
                    diagnostics["repeat_correlation"]["rho_median_pair"],
                "pcc_ceiling": diagnostics["pcc_ceiling"]["value"],
                "label": HEURISTIC,
            },
            "noise_floor": _record(output / "noise_floor.csv"),
            "provenance": {
                "accepted_view": self.context.source_manifests(),
                "sealed": self.sealed_binding,
                "fingerprints": fingerprint,
            },
            "adr_reference": ADR_REFERENCE[source],
            "merge_cross_check": merge_check,
            "delta_equals_view_delta_layer": bool(delta_identity),
            "environment": environment(),
        }
        write_json(output / "manifest.json", manifest)
        self.receipt = {
            "schema": SCHEMA,
            "source": source,
            "fixture": self.fixture,
            "output": str(output),
            "manifest_sha256": sha256(output / "manifest.json"),
            "counts": counts,
            "diagnostics": {
                "noise_floor_median": diagnostics["noise_floor"]["median"],
                "noise_floor_population_median":
                    diagnostics["noise_floor_population"]["median"],
                "repeat_correlation_median_pair":
                    diagnostics["repeat_correlation"]["rho_median_pair"],
                "repeat_correlation_median_condition":
                    diagnostics["repeat_correlation"]["rho_median_condition"],
                "pcc_ceiling": diagnostics["pcc_ceiling"]["value"],
                "adr_reference": ADR_REFERENCE[source],
            },
            "merge_cross_check": merge_check,
            "delta_equals_view_delta_layer": bool(delta_identity),
            "files": {
                "primary": {key: value["outputs"] for key, value in
                            primary_records.items()},
                "condition_mean": {key: value["outputs"] for key, value in
                                   condition_records.items()},
            },
        }
        return self.receipt


# --------------------------------------------------------------------------- #
# public entry points
# --------------------------------------------------------------------------- #
def build_truth(*, source, output_dir, repo_root_path=None, size=SIZE, fixture=False,
                corpus_dir=None, axis_dir=None, accepted_dir=None, sealed_dir=None,
                views_root=None, full_freeze_dir=None, record=None):
    """Build one source's frozen truth product; return its receipt."""
    truth = _Truth(source=source, output_dir=output_dir, repo_root_path=repo_root_path,
                   size=size, fixture=fixture, corpus_dir=corpus_dir, axis_dir=axis_dir,
                   accepted_dir=accepted_dir, sealed_dir=sealed_dir,
                   views_root=views_root, full_freeze_dir=full_freeze_dir)
    receipt = truth.build()
    if record is not None:
        write_json(record, receipt)
    return receipt


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", choices=(*SOURCES, "all"), default="all")
    parser.add_argument("--output", type=Path, default=HERE.parent / "output/truth",
                        help="truth root; each source lands under <output>/<source>/")
    parser.add_argument("--repo-root", type=Path, default=None)
    parser.add_argument("--size", type=int, default=SIZE)
    parser.add_argument("--fixture", action="store_true")
    parser.add_argument("--corpus", type=Path, default=None)
    parser.add_argument("--axis", type=Path, default=None)
    parser.add_argument("--accepted", type=Path, default=None)
    parser.add_argument("--sealed", type=Path, default=None,
                        help="sealed Step 04 freeze directory (one source at a time)")
    parser.add_argument("--views", type=Path, default=None,
                        help="views root for the merge-product cross-check")
    parser.add_argument("--full-freeze", type=Path, default=None,
                        help="canonical full freeze providing parents.csv/morgan.npy")
    parser.add_argument("--record", type=Path)
    args = parser.parse_args()
    sources = tuple(SOURCES) if args.source == "all" else (args.source,)
    if args.sealed is not None and len(sources) > 1:
        raise HarnessError("--sealed names one source; run the sources separately")
    receipts = []
    for source in sources:
        started = time.monotonic()
        receipt = build_truth(
            source=source, output_dir=args.output, repo_root_path=args.repo_root,
            size=args.size, fixture=args.fixture, corpus_dir=args.corpus,
            axis_dir=args.axis, accepted_dir=args.accepted, sealed_dir=args.sealed,
            views_root=args.views, full_freeze_dir=args.full_freeze)
        receipt["runtime_seconds"] = round(time.monotonic() - started, 3)
        receipts.append(receipt)
        sys.stdout.write(f"{source}: primary={receipt['counts']['primary']} "
                         f"condition_mean={receipt['counts']['condition_mean']}\n")
    result = receipts[0] if len(receipts) == 1 else {"schema": SCHEMA,
                                                     "sources": receipts}
    if args.record is not None:
        write_json(args.record, result)
    sys.stdout.write(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
