"""Readers for sealed raw-first input generations; retains study split and row APIs."""
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from expression_view import CONTROL, NORMALIZATION, SCHEMA as VIEW_SCHEMA, sha256
from prediction_axis import open_axes
from stage2_util.representations import WIDTHS
from stage2_util.seeds import read_seeds
from tsg_round1.inputs import REQUIRED, SOURCE_POPULATIONS
from tsg_round1.protocol import SPLIT_SEED, TRAINING_SEED, make_splits

PAIR_SCHEMA = "tsg-dose-raw-first-pair-v1"
INPUT_SCHEMA = "tsg-dose-raw-first-inputs-v1"
KPGT_SCHEMA = "tsg-dose-raw-first-kpgt-v1"
SOURCES = tuple(SOURCE_POPULATIONS)
GENE_COLUMNS = ["gene_index", "symbol"]
IDENTITY_COLUMNS = ["profile_id", "view_row", "role", "drug_parent_inchikey", "drug_parent_smiles",
                    "dose_um", "cell_line", "matched_control_group"]
PARENT_COLUMNS = ["drug_parent_inchikey", "drug_parent_smiles"]


def read_json(path):
    return json.loads(path.read_text())


def read_csv(path):
    return pd.read_csv(path, keep_default_na=False, float_precision="round_trip", low_memory=False)


def checked_files(directory, files, expected):
    if set(files) != set(expected):
        raise ValueError("raw-first manifest output coverage differs")
    for name, expected_digest in files.items():
        if sha256(directory / name) != expected_digest:
            raise ValueError(f"raw-first input digest differs: {name}")


def valid_digest(value):
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def check_generation(generation, axis):
    if (generation.get("schema") != VIEW_SCHEMA or generation.get("family") != "full"
            or generation.get("axis_size") != len(axis.gene_indices)
            or generation.get("axis_sha256") != axis.sha256
            or not all(valid_digest(generation.get(name)) for name in ("view_sha256", "corpus_sha256"))
            or generation.get("normalization") != NORMALIZATION or generation.get("control") != CONTROL):
        raise ValueError("raw-first view or prediction axis generation differs")


def match_controls(obs, vehicles):
    groups = {name: (np.asarray(positions) + len(obs)).tolist()
              for name, positions in vehicles.groupby("matched_control_group", observed=True, sort=False).indices.items()}
    controls = {index: groups.get(group) for index, group in enumerate(obs.matched_control_group)}
    if "" in groups or any(not members for members in controls.values()):
        raise ValueError("raw-first treated observation has no individual matched vehicles")
    return controls


class Inputs:
    """Open a raw-first source whose digest is pinned by a paired manifest."""

    def __init__(self, directory, *, expected_pair_sha256=None):
        self.directory = Path(directory)
        pair_path = self.directory.parent / "manifest.json"
        self.pair_sha256 = sha256(pair_path)
        if expected_pair_sha256 is not None and self.pair_sha256 != expected_pair_sha256:
            raise ValueError("raw-first pair digest differs")
        self.pair = read_json(pair_path)
        self.source = self.directory.name
        sources = self.pair.get("sources", {})
        if (self.pair.get("schema") != PAIR_SCHEMA or not set(sources) <= set(SOURCES)
                or self.source not in sources
                or sha256(self.directory / "manifest.json") != sources[self.source]):
            raise ValueError("raw-first input source is outside the pinned pair")
        self.manifest = self.read_json("manifest.json")
        generation = self.pair["view_generation"]
        if (self.manifest.get("schema") != INPUT_SCHEMA or self.manifest.get("source") != self.source
                or self.manifest.get("view_generation") != generation):
            raise ValueError("raw-first source input generation differs from the pair")
        axis = open_axes().axis(5000)
        check_generation(generation, axis)
        checked_files(self.directory, self.manifest["outputs"], REQUIRED)
        self.config = self.read_json("config.json")
        if (self.config.get("source") != self.source or self.config.get("view_generation") != generation
                or self.config.get("split_seed") != SPLIT_SEED
                or self.config.get("training_seed") != TRAINING_SEED
                or TRAINING_SEED not in read_seeds()
                or not all(valid_digest(self.config.get(name)) for name in
                           ("representation_sha256", "view_payload_sha256"))):
            raise ValueError("raw-first configuration provenance differs")
        self.obs = read_csv(self.directory / "population.csv").set_index("profile_id", drop=False)
        self.all_obs = read_csv(self.directory / "all_obs.csv").set_index("profile_id", drop=False)
        self.genes = read_csv(self.directory / "genes.csv")
        self.parents = read_csv(self.directory / "parents.csv")
        self.splits = read_csv(self.directory / "splits.csv")
        self.controls = {int(row): value for row, value in self.read_json("matched_controls.json").items()}
        self.expression = np.load(self.directory / "expression.npy", mmap_mode="r")
        self.control = np.load(self.directory / "control.npy", mmap_mode="r")
        self.morgan = np.load(self.directory / "morgan.npy", mmap_mode="r")
        self._validate(axis)
        self.parent_lookup = dict(zip(self.parents.drug_parent_inchikey, range(len(self.parents))))

    def _validate(self, axis):
        if (list(self.genes) != GENE_COLUMNS or tuple(self.genes.gene_index) != axis.gene_indices
                or tuple(self.genes.symbol) != axis.symbols):
            raise ValueError("raw-first ordered genes differ from the published axis")
        if (self.all_obs.index.has_duplicates or not self.obs.role.eq("treated").all()
                or not self.all_obs.iloc[:len(self.obs)][IDENTITY_COLUMNS].equals(self.obs[IDENTITY_COLUMNS])):
            raise ValueError("raw-first profile identities differ")
        vehicles = self.all_obs.iloc[len(self.obs):]
        if not vehicles.role.eq("vehicle").all() or self.controls != match_controls(self.obs, vehicles):
            raise ValueError("raw-first matched vehicle identities differ")
        shape = (len(self.all_obs), len(self.genes))
        if any(matrix.shape != shape or matrix.dtype != np.float32 for matrix in
               (self.expression, self.control)):
            raise ValueError("raw-first expression dimensions differ")
        if (not self.parents.equals(ordered_parents(self.obs))
                or self.morgan.shape != (len(self.parents), WIDTHS["ECFP4C"])
                or self.morgan.dtype != np.uint8 or (self.morgan > 1).any()):
            raise ValueError("raw-first parent identities or ECFP4 vectors differ")
        if not self.splits.equals(make_splits(self.obs)):
            raise ValueError("raw-first fold membership differs from the study protocol")

    def read_json(self, name):
        return read_json(self.directory / name)

    def rows(self, scenario, fold, partition):
        s = self.splits
        return s.loc[(s.scenario == scenario) & (s.fold == fold)
                     & (s.partition == partition), "row"].to_numpy(dtype=int)

    def fingerprints(self, rows):
        positions = self.obs.iloc[rows].drug_parent_inchikey.map(self.parent_lookup)
        return np.asarray(self.morgan[positions.to_numpy()], dtype=np.float32)

    def doses(self, rows):
        return self.obs.iloc[rows].dose_um.to_numpy(dtype=np.float32)


def ordered_parents(obs):
    parents = obs[PARENT_COLUMNS].drop_duplicates().sort_values(PARENT_COLUMNS[0], ignore_index=True)
    if parents.eq("").any().any() or parents.isna().any().any() or parents[PARENT_COLUMNS[0]].duplicated().any():
        raise ValueError("raw-first parent has missing or conflicting structures")
    return parents


class KPGTInputs(Inputs):
    """A tenv-only transport of shared KPGT vectors bound to the raw-first pair."""

    def __init__(self, directory, *, expected_pair_sha256=None):
        directory = Path(directory)
        receipt = read_json(directory / "manifest.json")
        if receipt.get("schema") != KPGT_SCHEMA:
            raise ValueError("not a raw-first KPGT input release")
        checked_files(directory, receipt["files"], (
            "pair/manifest.json", "pair/tenv/manifest.json", "KPGT.npy", "representation_parents.csv"
        ))
        super().__init__(directory / "pair/tenv", expected_pair_sha256=expected_pair_sha256)
        if receipt.get("representation_sha256") != self.config["representation_sha256"]:
            raise ValueError("KPGT representation provenance differs from the tenv input")
        if not read_csv(directory / "representation_parents.csv").equals(self.parents):
            raise ValueError("KPGT parent identities are not aligned to tenv")
        self.kpgt = np.load(directory / "KPGT.npy", mmap_mode="r")
        if (self.kpgt.shape != (len(self.parents), WIDTHS["KPGT"])
                or self.kpgt.dtype != np.float32 or not np.isfinite(self.kpgt).all()):
            raise ValueError("invalid KPGT vectors for raw-first tenv")
        self.directory = directory
        self.manifest = receipt

    def fingerprints(self, rows):
        positions = self.obs.iloc[rows].drug_parent_inchikey.map(self.parent_lookup)
        return self.kpgt[positions.to_numpy()]
