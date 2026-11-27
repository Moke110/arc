"""Bind Step 10 to the released cold-dose-drug partitions and same-plate controls."""

from contextlib import contextmanager
from dataclasses import dataclass
from hashlib import sha256

import numpy as np

from dose_model import model_features
from dose_protocol import AXIS, DEFAULT_STEP03, HELD_DOSES, SCENARIO, STAGE, representation
from training_data import open_tenv_axis_training
from tsg_axis.inputs import AxisInputs
from stage2_util.representations import RepresentationTable
from util.training.chemical import aligned_features

PARTITIONS = ("train", "valid", "test")
SMOKE_LIMITS = dict(zip(PARTITIONS, (64, 24, 24), strict=True))


def cold_dose_partitions(treatments, fold):
    """Check the full intersection and both excluded strips against parent and dose."""
    if type(fold) is not int or fold not in range(len(HELD_DOSES)):
        raise ValueError("unknown cold-dose-drug fold")
    label = f"{SCENARIO}_{fold}"
    labels = treatments[label]
    parents = treatments.drug_parent
    held_parents = parents.loc[treatments[f"cold_drug_{fold}"].eq("test")].unique()
    held_parent = parents.isin(held_parents)
    held_dose = treatments.dose_um.eq(HELD_DOSES[fold])
    test = held_parent & held_dose
    development = ~held_parent & ~held_dose
    if not (labels.eq("test").equals(test) and labels.isin(("train", "valid")).equals(development)):
        raise ValueError(f"fold {fold} does not split the held-parent and held-dose intersection")
    reasons = treatments[f"{label}_reason"]
    if not (reasons.eq("held_parent").equals(held_parent & ~held_dose)
            and reasons.eq("held_dose").equals(~held_parent & held_dose)
            and labels.eq("").equals(held_parent ^ held_dose)):
        raise ValueError(f"fold {fold} has an incorrect single-holdout exclusion")
    result = {phase: np.flatnonzero(labels.eq(phase).to_numpy()) for phase in PARTITIONS}
    if any(not len(rows) for rows in result.values()):
        raise ValueError(f"fold {fold} has an empty partition")
    return result


class ColdDoseInputs(AxisInputs):
    def __init__(self, release):
        super().__init__(release)
        self.partitions = tuple(cold_dose_partitions(self.treatments, fold) for fold in range(3))

    def rows(self, fold, partition):
        return self.partitions[fold][partition]

    def positions(self, fold, bounded):
        return {phase: rows[:SMOKE_LIMITS[phase]] if bounded else rows
                for phase, rows in self.partitions[fold].items()}

    def fold_population(self, fold, bounded):
        selected = self.positions(fold, bounded)
        targets = self.model_input("test")
        return {
            "held_dose_um": HELD_DOSES[fold],
            "held_parents": sorted(set(self.parents[self.treatments[f"cold_drug_{fold}"].eq("test")])),
            "full_counts": {phase: len(rows) for phase, rows in self.partitions[fold].items()},
            "excluded_strips": {reason: int(self.treatments[f"{SCENARIO}_{fold}_reason"].eq(reason).sum())
                                for reason in ("held_parent", "held_dose")},
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
