"""Verified view-aware freeze for one `tsg_view` product.

Schema ``tsg-view-freeze-v1``. The file set is the Step 04 axis-local ``REQUIRED`` set
plus ``control_sources.csv``; the semantics differ deliberately (plan C7.6):

* ``expression.npy``/``control.npy`` are indexed by **population row** (0..n_pop-1) and
  hold each row's own expression and its paired control. There is no duplicated vehicle
  tail, so a derived view does not pay for 65k vehicle rows.
* ``all_obs.csv`` lists the population rows plus the distinct control-source vehicle rows
  (metadata only), and ``matched_controls.json`` maps a population row to those rows'
  offsets in ``all_obs`` for audit.
* ``control_sources.csv`` records each row's control recipe and the pooled retained total.
* ``splits.csv`` partitions come only from the sealed Step 04 ``(parent, dose)`` mapping;
  for ``full`` the file is byte-equal to the sealed ``splits.csv``.

Verification is fail-closed and re-reads the sealed freeze (read-only) to prove the fold
mapping, the parent population and the Morgan recipe did not change.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from views import CONDITION_KEYS, HarnessError, SOURCES, VIEWS, bind, modules, relative, sha256

FREEZE_SCHEMA = "tsg-view-freeze-v1"
REQUIRED = ("config.json", "population.csv", "all_obs.csv", "parents.csv", "genes.csv",
            "splits.csv", "matched_controls.json", "expression.npy", "control.npy",
            "morgan.npy")
COMPANION = "control_sources.csv"
OUTPUTS = REQUIRED + (COMPANION,)
POPULATION_COLUMNS = (
    "profile_id", "source_item", "source_item_manifest_sha256", "source_row", "role",
    "plate", "cell_line", "cell_line_id", "compound_name", "vehicle", "dose_um",
    "dose_original_value", "dose_original_unit", "drug_parent_inchikey",
    "drug_parent_smiles", "source_inchikey", "n_cells", "matched_control_group",
    "library_total", "view_row")
SPLIT_COLUMNS = ("row", "profile_id", "drug_parent", "dose_um", "scenario", "fold",
                 "partition")
DELTA_KINDS = ("expression_minus_control", "pooled")


def sealed_reference(root, source):
    """Resolve the sealed Step 04 axis-local freeze for one source."""
    if source not in SOURCES:
        raise HarnessError(f"unknown source: {source}")
    path = Path(root) / "stage2/04_tsg_dose/execution/checkpoints/round1-axis-local-inputs" / source
    if not (path / "splits.csv").is_file():
        raise HarnessError(
            f"the sealed Step 04 freeze is unavailable: {path}; run where the payload "
            "exists or pass --sealed")
    return path


def sealed_binding(path, root):
    path = Path(path)
    return {"file": relative(path, root),
            "splits_sha256": sha256(path / "splits.csv"),
            "parents_sha256": sha256(path / "parents.csv"),
            "morgan_sha256": sha256(path / "morgan.npy"),
            "population_sha256": sha256(path / "population.csv")}


def read_sealed_splits(path):
    splits = pd.read_csv(Path(path) / "splits.csv")
    return splits


def partition_map(splits):
    """The sealed ``(scenario, fold, parent, dose) -> partition`` table."""
    key = ["scenario", "fold", "drug_parent", "dose_um"]
    mapping = splits[key + ["partition"]].drop_duplicates()
    if mapping.duplicated(key).any():
        raise HarnessError("the sealed splits map one condition key to two partitions")
    return mapping


def check_splits(splits, population, sealed, source):
    """Cross-check a view's splits against the sealed fold definition and its population.

    Returns the sealed partition map. Raises ``HarnessError`` naming the failing rule.
    """
    if list(splits.columns) != list(SPLIT_COLUMNS):
        raise HarnessError(f"{source}: splits.csv columns differ from the Step 04 contract")
    if len(splits) != 6 * len(population):
        raise HarnessError(f"{source}: splits.csv does not cover six scenario x fold blocks")
    counts = np.bincount(splits["row"].to_numpy(dtype=int), minlength=len(population))
    if len(counts) != len(population) or not np.all(counts == 6):
        raise HarnessError(f"{source}: splits.csv row space is not the population row space")
    sealed_map = partition_map(sealed)
    key = ["scenario", "fold", "drug_parent", "dose_um"]
    view = splits[key + ["partition"]].drop_duplicates()
    if view.duplicated(key).any():
        raise HarnessError(
            f"{source}: a condition key maps to more than one partition (members split "
            "across folds)")
    joined = sealed_map.merge(view, on=key, how="outer", suffixes=("_sealed", "_view"),
                              indicator=True)
    if not (joined["_merge"] == "both").all():
        missing = joined[joined["_merge"] != "both"].head(5)
        raise HarnessError(
            f"{source}: the sealed (parent, dose) population is not covered: "
            f"{missing[key + ['_merge']].to_dict('records')}")
    if not (joined["partition_sealed"] == joined["partition_view"]).all():
        differ = joined[joined["partition_sealed"] != joined["partition_view"]].head(5)
        raise HarnessError(
            f"{source}: a (parent, dose) key changed partition: "
            f"{differ[key + ['partition_sealed', 'partition_view']].to_dict('records')}")
    sealed_keys = set(zip(sealed["drug_parent"], sealed["dose_um"]))
    if set(zip(joined["drug_parent"], joined["dose_um"])) != sealed_keys:
        raise HarnessError(f"{source}: the sealed parent/dose population differs")
    # Exactly one partition per condition per scenario x fold.
    conditions = population[list(CONDITION_KEYS)].astype(str)
    conditions["dose_um"] = population["dose_um"].to_numpy(dtype=float)
    condition_id = pd.factorize(pd.MultiIndex.from_frame(conditions), sort=True)[0]
    population = population.assign(_condition=condition_id)
    split_rows = splits.assign(_condition=condition_id[splits["row"].to_numpy(dtype=int)])
    grouped = split_rows.groupby(["scenario", "fold", "_condition"])["partition"].nunique()
    if (grouped > 1).any():
        raise HarnessError(
            f"{source}: a condition's members sit in different partitions: "
            f"{grouped[grouped > 1].head(3).to_dict()}")
    # Every partition is disjoint and the five partitions exactly cover each block.
    known = {"train", "validation", "test", "excluded-held-parent", "excluded-held-dose"}
    for (scenario, fold), block in split_rows.groupby(["scenario", "fold"]):
        members = block.groupby("partition")["row"].apply(set).to_dict()
        unknown = set(members) - known
        if unknown:
            raise HarnessError(
                f"{source}: unknown partition in {scenario} fold {fold}: {sorted(unknown)}")
        names = sorted(members)
        for left in range(len(names)):
            for right in range(left + 1, len(names)):
                if members[names[left]] & members[names[right]]:
                    raise HarnessError(
                        f"{source}: partitions overlap in {scenario} fold {fold}")
        if set().union(*members.values()) != set(range(len(population))):
            raise HarnessError(
                f"{source}: partitions do not cover {scenario} fold {fold}")
    return sealed_map


class Inputs:
    """Verified view freeze: identity, digests, fold mapping and paired expression."""

    def __init__(self, directory, sealed=None, repo_root=None):
        self.directory = Path(directory)
        self.root = bind(repo_root or _root(self.directory))["root"]
        self.manifest = self.read_json("manifest.json")
        if self.manifest.get("schema") != FREEZE_SCHEMA:
            raise HarnessError("view freeze manifest schema mismatch")
        if set(self.manifest.get("outputs", {})) != set(OUTPUTS):
            raise HarnessError("view freeze does not declare exactly its required files")
        for name, digest in self.manifest["outputs"].items():
            if sha256(self.directory / name) != digest:
                raise HarnessError(f"view freeze digest mismatch: {name}")
        self.config = self.read_json("config.json")
        self.source = self.config["source"]
        self.view = self.config["view"]
        if self.source not in SOURCES or self.view not in VIEWS:
            raise HarnessError("unknown view freeze identity")
        axis = self.config.get("view_generation", {}).get("axis", {})
        if int(axis.get("size", -1)) != 5000:
            raise HarnessError("a view freeze requires the G_5000 axis")
        if self.config["view_generation"]["normalization"].get("source") != (
                "ADR 0077 five-k-only axis-local normalization"):
            raise HarnessError("view freeze normalization differs from ADR 0077")
        if self.config.get("delta_kind") not in DELTA_KINDS:
            raise HarnessError("view freeze delta kind is unknown")
        if self.config.get("split_seed") != 226:
            raise HarnessError("view freeze does not use split seed 226")

        self.population = pd.read_csv(self.directory / "population.csv")
        if list(self.population.columns) != list(POPULATION_COLUMNS):
            raise HarnessError("view freeze population columns differ from the contract")
        self.population = self.population.set_index("profile_id", drop=False)
        if self.population.index.duplicated().any():
            raise HarnessError("view freeze population profile_id values are not unique")
        if not self.population["role"].eq("treated").all():
            raise HarnessError("view freeze population must be treated rows only")
        if len(self.population) != int(self.config["n_treated"]):
            raise HarnessError("view freeze population count differs from its config")
        self.all_obs = pd.read_csv(self.directory / "all_obs.csv", low_memory=False)
        if list(self.all_obs.columns) != list(POPULATION_COLUMNS):
            raise HarnessError("view freeze all_obs columns differ from the contract")
        self.all_obs = self.all_obs.set_index("profile_id", drop=False)
        if self.all_obs.index.duplicated().any():
            raise HarnessError("view freeze all_obs profile_id values are not unique")
        keys = ["profile_id", "view_row", "role", "drug_parent_inchikey", "dose_um",
                "cell_line", "matched_control_group"]
        if not self.all_obs.iloc[:len(self.population)][keys].equals(self.population[keys]):
            raise HarnessError("view freeze profile identity mismatch")
        self.controls_sources = pd.read_csv(self.directory / COMPANION)
        if list(self.controls_sources.columns) != [
                "profile_id", "control_rule", "control_groups", "pooled_retained_total"]:
            raise HarnessError("view freeze control_sources columns differ from the contract")
        if list(self.controls_sources["profile_id"].astype(str)) != list(
                self.population["profile_id"].astype(str)):
            raise HarnessError("view freeze control_sources does not cover the population")

        self.genes = pd.read_csv(self.directory / "genes.csv")
        if list(self.genes.columns) != ["gene_index", "symbol"] or len(self.genes) != 5000:
            raise HarnessError("view freeze genes.csv is not the G_5000 axis")
        if sha256(self.directory / "genes.csv") != self.config["view_generation"]["axis"][
                "columns_sha256"]:
            raise HarnessError("view freeze genes.csv differs from the released axis")
        self.parents = pd.read_csv(self.directory / "parents.csv")
        self.splits = pd.read_csv(self.directory / "splits.csv")
        self.controls = {int(key): value for key, value in
                         self.read_json("matched_controls.json").items()}
        for name in ("expression", "control"):
            setattr(self, name, np.load(self.directory / f"{name}.npy", mmap_mode="r"))
            matrix = getattr(self, name)
            if matrix.shape != (len(self.population), 5000) or matrix.dtype != np.float32:
                raise HarnessError(f"view freeze {name}.npy shape or dtype mismatch")
        self.morgan = np.load(self.directory / "morgan.npy", mmap_mode="r")
        if (self.morgan.shape != (len(self.parents), 2048)
                or self.morgan.dtype != np.uint8
                or self.parents["drug_parent_inchikey"].duplicated().any()
                or set(self.parents["drug_parent_inchikey"])
                != set(self.population["drug_parent_inchikey"])):
            raise HarnessError("view freeze Morgan identity mismatch")

        # -- sealed fold and recipe cross-checks ------------------------------
        self.sealed = Path(sealed) if sealed else sealed_reference(self.root, self.source)
        binding = sealed_binding(self.sealed, self.root)
        for name in ("splits_sha256", "parents_sha256", "morgan_sha256"):
            if binding[name] != self.config["sealed"][name]:
                raise HarnessError(f"the sealed Step 04 freeze changed: {name}")
        if sha256(self.directory / "morgan.npy") != binding["morgan_sha256"]:
            raise HarnessError("view freeze morgan.npy is not the sealed Morgan recipe")
        sealed_parents = pd.read_csv(self.sealed / "parents.csv")
        if not sealed_parents.equals(self.parents):
            raise HarnessError("view freeze parents.csv differs from the sealed parents")
        sealed = read_sealed_splits(self.sealed)
        check_splits(self.splits, self.population, sealed, self.source)
        if self.view == "full":
            if sha256(self.directory / "splits.csv") != binding["splits_sha256"]:
                raise HarnessError("full view splits.csv is not byte-equal to the sealed freeze")

        # -- paired control contract -----------------------------------------
        vehicles = self.all_obs.iloc[len(self.population):]
        if not vehicles["role"].eq("vehicle").all():
            raise HarnessError("view freeze all_obs tail must be vehicle metadata")
        if set(self.controls) != set(range(len(self.population))):
            raise HarnessError("view freeze matched_controls does not cover the population")
        body = {str(group): (positions + len(self.population)).tolist()
                for group, positions in vehicles.groupby("matched_control_group").indices.items()}
        for position, row in enumerate(self.controls_sources.itertuples(index=False)):
            groups = [group for group in str(row.control_groups).split(";") if group]
            expected = sorted({offset for group in groups
                               for offset in body.get(group, [])})
            if sorted(self.controls[position]) != expected:
                raise HarnessError(
                    "view freeze matched_controls does not resolve the recorded groups")
        self._check_control_composition()

    # -- helpers ---------------------------------------------------------- #
    def read_json(self, name):
        return json.loads((self.directory / name).read_text())

    def _check_control_composition(self):
        """Prove the control/expression recipe from the stored arrays alone.

        A normalized expression row is a CP10K composition (``sum(expm1(X)) == 10000``).
        A pooled-sum control is the same composition evaluated on the summed retained
        integer counts, so it recovers near-integer counts on the recorded total; a
        mean-of-normalized-vehicles control is a mean of log compositions and therefore
        only has to stay below the CP10K bound.
        """
        control = np.asarray(self.control, dtype=np.float64)
        if not np.isfinite(control).all() or (control < 0).any():
            raise HarnessError("a view freeze control value is negative or non-finite")
        composition = np.expm1(control).sum(axis=1)
        if (composition > 10000.0 * (1 + 1e-4)).any():
            raise HarnessError("a view freeze control exceeds a CP10K composition")
        pooled = self.controls_sources["control_rule"].eq("pooled-sum").to_numpy()
        if pooled.any():
            if not np.allclose(composition[pooled], 10000.0, rtol=1e-5, atol=1.0):
                raise HarnessError("a pooled-sum control is not a CP10K composition")
            totals = self.controls_sources["pooled_retained_total"].to_numpy(dtype=float)[
                pooled]
            if not np.isfinite(totals).all() or (totals <= 0).any():
                raise HarnessError("a pooled-sum control has no retained total")
            counts = np.expm1(control[pooled]) * totals[:, None] / 10000.0
            if (counts < 0).any() or not np.allclose(counts, np.round(counts), rtol=1e-5,
                                                      atol=1e-3):
                raise HarnessError(
                    "a pooled-sum control does not recover integer retained counts")
        if self.config["delta_kind"] == "expression_minus_control":
            expression = np.expm1(np.asarray(self.expression, dtype=np.float64)).sum(axis=1)
            if not np.allclose(expression, 10000.0, rtol=1e-5, atol=1.0):
                raise HarnessError(
                    "a view freeze expression row is not a CP10K composition")

    def rows(self, scenario, fold, partition):
        splits = self.splits
        return splits.loc[(splits.scenario == scenario) & (splits.fold == fold)
                          & (splits.partition == partition), "row"].to_numpy(dtype=int)

    def fingerprints(self, rows):
        positions = self.population.iloc[rows]["drug_parent_inchikey"].map(
            dict(zip(self.parents["drug_parent_inchikey"], range(len(self.parents)))))
        return np.asarray(self.morgan[positions.to_numpy()], dtype=np.float32)

    def doses(self, rows):
        return self.population.iloc[rows]["dose_um"].to_numpy(dtype=np.float32)

    def delta(self, rows):
        rows = np.asarray(rows, dtype=int)
        return (np.asarray(self.expression[rows], dtype=np.float64)
                - np.asarray(self.control[rows], dtype=np.float64))


def _root(start):
    for path in [Path(start).resolve(), *Path(start).resolve().parents]:
        if (path / "stage2/01_bases/output/api/bases/view.py").is_file():
            return path
    return modules()["root"] if modules() else None
