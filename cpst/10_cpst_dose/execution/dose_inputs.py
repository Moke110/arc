"""Bind Step 10 to the released cold-dose-drug partitions and same-plate controls."""

from contextlib import contextmanager
from dataclasses import dataclass
from hashlib import sha256

import numpy as np

from dose_model import model_features
from dose_protocol import AXIS, DEFAULT_STEP03, HELD_DOSES, STAGE, fold_numbers, protocol, representation
from training_data import open_tenv_axis_training
from tsg_axis.inputs import AxisInputs
from cpst_util.representations import RepresentationTable
from util.training.chemical import aligned_features

PARTITIONS = ("train", "valid", "test")
SMOKE_LIMITS = dict(zip(PARTITIONS, (64, 24, 24), strict=True))


def cold_dose_partitions(treatments, fold):
    """Project the frozen crossed fold into Step 03 treatment row order."""
    if type(fold) is not int or fold not in fold_numbers():
        raise ValueError("unknown cold-dose-drug fold")
    positions = protocol().folds[fold].positions_for_treatments(treatments)
    return {phase: np.asarray(rows, dtype=np.intp) for phase, rows in positions.items()}


class ColdDoseInputs(AxisInputs):
    def __init__(self, release):
        super().__init__(release)
        self.partitions = tuple(cold_dose_partitions(self.treatments, fold) for fold in fold_numbers())

    def rows(self, fold, partition):
        return self.partitions[fold][partition]

    def positions(self, fold, bounded):
        if bounded:
            return {phase: self.rows(fold, phase)[:limit] for phase, limit in SMOKE_LIMITS.items()}
        return {phase: self.rows(fold, phase) for phase in PARTITIONS}

    def fold_population(self, fold, bounded):
        selected = self.positions(fold, bounded)
        targets = self.model_input("test")
        split = protocol().folds[fold]
        return {
            "parent_group": split.parent_group,
            "held_dose_um": split.held_dose_um,
            "held_parents": sorted(split.held_parents),
            "full_counts": {phase: len(self.rows(fold, phase)) for phase in PARTITIONS},
            "excluded_strips": {reason: int(len(self.rows(fold, f"excluded-{reason}")))
                                for reason in ("held-parent", "held-dose")},
            "excluded_conditions": {reason: list(self.condition_ids(self.rows(fold, f"excluded-{reason}")))
                                    for reason in ("held-parent", "held-dose")},
            "control_sha256": {phase: sha256(targets.control[rows].tobytes()).hexdigest()
                               for phase, rows in selected.items()},
        }


@dataclass(frozen=True)
class ModelInputs:
    data: ColdDoseInputs
    features: np.ndarray

    def fingerprints(self, rows):
        return self.features[rows]


def features_for(data, vectors, arm):
    chemistry = aligned_features(vectors, representation(arm), data.parents)
    return ModelInputs(data, model_features(arm, chemistry, data.doses))


@contextmanager
def open_inputs(step03=DEFAULT_STEP03, representations=None):
    """Verify the frozen release, ordered axis, fold masks and chemical coverage."""
    with open_tenv_axis_training(AXIS, step03) as release:
        pinned = STAGE / "03_data_train/03.7_tenv_denoise-5k/output/manifest.json"
        axis = STAGE / "02_gene_axis/02.2_denoise-deg/output/release/gene_axis_5000.csv"
        if (release.manifest_sha256 != sha256(pinned.read_bytes()).hexdigest()
                or release.manifest["outputs"]["genes.csv"] != sha256(axis.read_bytes()).hexdigest()
                or len(release.genes) != 5000
                or tuple(release.manifest["inputs"]["folds"]["held_dose_um"]) != HELD_DOSES):
            raise ValueError("Step 10 requires the frozen Denoise-DEG 5k and cold-dose-drug releases")
        data = ColdDoseInputs(release)
        vectors = None if representations is None else RepresentationTable(representations, source="tenv")
        if vectors is not None and set(data.parents) != set(vectors.parents.drug_parent_inchikey):
            raise ValueError("frozen chemical parent population differs from the Step 10 training release")
        yield data, vectors
