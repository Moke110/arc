"""Prepare and evaluate the single-dose comparison on raw-first HVG-DEG views."""
import json
from pathlib import Path
import re
import sys
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

STAGE2 = Path(__file__).resolve().parents[4]
for relative in ("util/output/api", "03_data_train/03.1_hvg-deg-5k/output/api",
                 "02_gene_axis/02.1_hvg-deg/output/api"):
    sys.path.insert(0, str(STAGE2 / relative))

from expression_view import DEFAULT, PAYLOAD_SCHEMA, SCHEMA as VIEW_SCHEMA, open_view, sha256
from prediction_axis import open_axes
from stage2_util import metrics, palette, plotting
from stage2_util.representations import RepresentationTable, ordered_parents
from stage2_util.seeds import read_seeds
from tenv5k.protocol import FOLDS, make_splits

SCHEMA = "tenv5k-raw-first-inputs-v1"
SPLIT_SEED = 226
AXIS_SIZE = 5000
DOSE = 1
MODELS = ("PRnet", "TranSiGen", "ECFP4_MLP", "Mean", "Mean_cell", "Mean_drug")
TABLES = ("population", "all_obs", "genes", "parents", "splits")
ARRAYS = ("expression", "control", "delta", "morgan")
FILES = {*(f"{name}.csv" for name in TABLES), *(f"{name}.npy" for name in ARRAYS),
         "matched_controls.json", "config.json"}
LABELS = {"pcc": "Pearson", "spearman": "Spearman", "rmse": "RMSE", "r2": "R²"}


def read_json(path):
    return json.loads(path.read_text())


def write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")


def read_csv(path):
    return pd.read_csv(path, keep_default_na=False, float_precision="round_trip", low_memory=False)


def control_mapping(treated, vehicles):
    groups = vehicles.groupby("matched_control_group", sort=False, observed=True).indices
    if "" in groups or not set(treated.matched_control_group) <= set(groups):
        raise ValueError("single-dose observation lacks matching vehicle rows")
    return {i: (groups[group] + len(treated)).tolist()
            for i, group in enumerate(treated.matched_control_group)}


class Inputs:
    def __init__(self, directory):
        self.directory = Path(directory)
        self.manifest = read_json(self.directory / "manifest.json")
        if self.manifest.get("schema") != SCHEMA or set(self.manifest["outputs"]) != FILES:
            raise ValueError("raw-first single-dose input schema or outputs differ")
        for name, digest in self.manifest["outputs"].items():
            if sha256(self.directory / name) != digest:
                raise ValueError(f"raw-first single-dose payload differs: {name}")
        self.config = read_json(self.directory / "config.json")
        generation = self.config["view_generation"]
        axis = open_axes().axis(AXIS_SIZE)
        if (generation != self.manifest["view_generation"]
                or generation["schema"] not in (VIEW_SCHEMA, PAYLOAD_SCHEMA)
                or generation["family"] != "full" or generation["source"] != "tenv"
                or generation["axis_sha256"] != axis.sha256 or self.config["dose_um"] != DOSE
                or self.config["axis_size"] != AXIS_SIZE
                or any(re.fullmatch(r"[0-9a-f]{64}", generation.get(field, "")) is None
                       for field in ("view_sha256", "corpus_sha256", "payload_sha256"))
                or self.config["training_seeds"] != list(read_seeds())
                or self.config["split_seed"] != SPLIT_SEED or self.config["folds"] != FOLDS):
            raise ValueError("invalid raw-first single-dose view or study generation")
        self.obs = read_csv(self.directory / "population.csv").set_index("profile_id", drop=False)
        self.all_obs = read_csv(self.directory / "all_obs.csv").set_index("profile_id", drop=False)
        self.parents = read_csv(self.directory / "parents.csv")
        self.genes = read_csv(self.directory / "genes.csv")
        self.splits = read_csv(self.directory / "splits.csv")
        self.controls = {int(i): rows for i, rows in read_json(self.directory / "matched_controls.json").items()}
        for name in ARRAYS:
            setattr(self, name, np.load(self.directory / f"{name}.npy", mmap_mode="r"))
        self._validate(axis)
        self.parent_lookup = dict(zip(self.parents.drug_parent_inchikey, range(len(self.parents))))

    def _validate(self, axis):
        tested = ["profile_id", "view_row", "role", "matched_control_group",
                  "drug_parent_inchikey", "drug_parent_smiles", "dose_um", "cell_line"]
        if (self.all_obs.index.has_duplicates or not self.obs.role.eq("treated").all()
                or not self.obs.dose_um.eq(DOSE).all() or set(self.obs.cell_line) != {"hiPSC-CM"}
                or not self.obs[tested].equals(self.all_obs.iloc[:len(self.obs)][tested])):
            raise ValueError("raw-first single-dose population differs")
        vehicles = self.all_obs.iloc[len(self.obs):]
        if not vehicles.role.eq("vehicle").all() or self.controls != control_mapping(self.obs, vehicles):
            raise ValueError("raw-first single-dose vehicle rows differ")
        if (tuple(self.genes.gene_index) != axis.gene_indices or tuple(self.genes.symbol) != axis.symbols
                or any(array.shape != (len(self.all_obs), AXIS_SIZE) or array.dtype != np.float32
                       for array in (self.expression, self.control, self.delta))):
            raise ValueError("raw-first single-dose gene axis or array shape differs")
        if (not self.parents.equals(ordered_parents(self.obs)) or self.morgan.dtype != np.uint8
                or self.morgan.shape != (len(self.parents), 2048) or (self.morgan > 1).any()):
            raise ValueError("raw-first single-dose ECFP4 identities differ")
        if not self.splits.equals(make_splits(self.obs, seed=SPLIT_SEED)):
            raise ValueError("raw-first single-dose fold membership differs")

    def rows(self, scenario, fold, partition):
        splits = self.splits
        return splits.loc[(splits.scenario == scenario) & (splits.fold == fold)
                          & (splits.partition == partition), "row"].to_numpy(dtype=int)

    def fingerprints(self, rows):
        keys = self.obs.iloc[rows].drug_parent_inchikey.map(self.parent_lookup)
        return np.asarray(self.morgan[keys.to_numpy()], dtype=np.float32)


def _freeze(source, destination, axis, representations):
    if source.manifest.get("default_family") != "full":
        raise ValueError("single-dose preparation requires the default full view")
    obs = source.obs.copy()
    if not obs.index.astype(str).equals(pd.Index(obs.profile_id.astype(str))):
        raise ValueError("view profile identities differ from source observations")
    obs["view_row"] = np.arange(len(obs))
    selected = np.flatnonzero(obs.role.eq("treated") & obs.dose_um.eq(DOSE))
    population = obs.iloc[selected]
    vehicles = np.flatnonzero(obs.role.eq("vehicle") &
                              obs.matched_control_group.isin(population.matched_control_group))
    rows = np.r_[selected, vehicles]
    for name, layer in (("expression", "X"), ("control", "control"), ("delta", "delta")):
        matrix = np.lib.format.open_memmap(destination / f"{name}.npy", mode="w+",
                                           shape=(len(rows), AXIS_SIZE), dtype=np.float32)
        for start, part in source.batches(layer):
            mask = (rows >= start) & (rows < start + len(part))
            if not np.isfinite(part).all():
                raise ValueError(f"raw-first {layer} contains nonfinite values")
            matrix[mask] = part[rows[mask] - start]
        matrix.flush()
        del matrix
    parents = ordered_parents(population)
    keys = representations.verified_parent_keys(parents)
    np.save(destination / "morgan.npy", representations.vectors("ECFP4C", keys).astype(np.uint8))
    for name, frame in (("population", population), ("all_obs", obs.iloc[rows]), ("parents", parents),
                        ("splits", make_splits(population, seed=SPLIT_SEED)),
                        ("genes", source.var[["gene_index", "symbol"]])):
        frame.to_csv(destination / f"{name}.csv", index=False)
    write_json(destination / "matched_controls.json", control_mapping(population, obs.iloc[vehicles]))
    generation = {
        "schema": source.manifest["schema"], "source": source.source, "family": source.family,
        "view_sha256": source.manifest_sha256, "axis_sha256": axis.sha256,
        "corpus_sha256": source.manifest["corpus_sha256"],
        "payload_sha256": source.payload_sha256,
    }
    write_json(destination / "config.json", {
        "view_generation": generation, "representation_sha256": representations.manifest_sha256,
        "dose_um": DOSE, "split_seed": SPLIT_SEED, "training_seeds": list(read_seeds()),
        "folds": FOLDS, "axis_size": AXIS_SIZE,
    })
    write_json(destination / "manifest.json", {
        "schema": SCHEMA, "view_generation": generation,
        "outputs": {name: sha256(destination / name) for name in sorted(FILES)},
    })


def prepare(destination, *, view_path=None, representations=None, representation_sha256=None):
    destination = Path(destination)
    root = DEFAULT if view_path is None else Path(view_path)
    if destination.exists() or root.resolve() in destination.resolve().parents:
        raise ValueError("single-dose preparation needs a fresh destination outside the view")
    axis = open_axes().axis(AXIS_SIZE)
    table = RepresentationTable(representations, source="tenv", manifest_sha256=representation_sha256)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".single-dose-", dir=destination.parent) as temp:
        staging = Path(temp) / "inputs"
        staging.mkdir()
        with open_view(view_path, source="tenv", axis_sha256=axis.sha256) as source:
            _freeze(source, staging, axis, table)
        Inputs(staging)
        staging.rename(destination)
    return destination


def score_prediction(obs, truth, prediction, control):
    """Preserve well and parent aggregation while using shared vector formulas."""
    truth_delta = truth - control
    prediction_delta = prediction - control
    per_well, per_parent = [], []
    for index, parent in enumerate(obs.drug_parent_inchikey):
        expression = metrics.vector_metrics(truth[index], prediction[index])
        per_well.append({"profile_id": str(obs.index[index]), "drug_parent": parent,
                         **metrics.vector_metrics(truth_delta[index], prediction_delta[index]),
                         **{f"expression_{key}": value for key, value in expression.items()}})
    for parent, indices in obs.groupby("drug_parent_inchikey", sort=True).indices.items():
        per_parent.append({"drug_parent": parent, "n_wells": len(indices),
                           **metrics.vector_metrics(truth_delta[indices].mean(axis=0),
                                                    prediction_delta[indices].mean(axis=0))})
    return pd.DataFrame(per_well), pd.DataFrame(per_parent)


def evaluate(inputs, prediction_path, identity_path, output, *, scenario, fold, model, seed, run_id):
    data = Inputs(inputs)
    if (model not in MODELS or (model == "Mean_drug" and scenario == "cold-drug")
            or seed not in data.config["training_seeds"] or not run_id.strip()):
        raise ValueError("unregistered single-dose model, seed or run identity")
    identity = {"model": model, "seed": seed, "run_id": run_id}
    output = Path(output)
    rows = data.rows(scenario, fold, "test")
    ids = read_csv(identity_path)
    prediction = np.load(prediction_path, mmap_mode="r")
    if (not len(rows) or list(ids.columns) != ["profile_id"]
            or ids.profile_id.tolist() != data.obs.index[rows].tolist()):
        raise ValueError("prediction identity order mismatch")
    if prediction.shape != (len(rows), AXIS_SIZE):
        raise ValueError("prediction shape differs from the frozen test rows and genes")
    per_well, per_parent = score_prediction(data.obs.iloc[rows], data.expression[rows],
                                            prediction, data.control[rows])
    if output.exists():
        raise FileExistsError(f"evaluation output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".single-dose-score-", dir=output.parent) as temp:
        staging = Path(temp) / "report"
        staging.mkdir()
        per_well.assign(**identity).to_csv(staging / "per_well.csv", index=False)
        per_parent.assign(**identity).to_csv(staging / "per_parent.csv", index=False)
        provenance = data.config["view_generation"]
        pd.DataFrame([{**identity, "scenario": scenario, "fold": fold, "n_parents": len(per_parent),
                       "view_sha256": provenance["view_sha256"], "axis_sha256": provenance["axis_sha256"],
                       **{name: per_parent[name].mean() for name in LABELS}}]).to_csv(
                           staging / "summary.csv", index=False)
        staging.rename(output)


def render(summary, output):
    """Read recorded parent-level metrics from CSV and render with the shared palette."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    rows = plotting.read_rows(summary)
    if not rows:
        raise ValueError("single-dose summary CSV is empty")
    figure = Figure(figsize=(8, 5), layout="constrained")
    FigureCanvasAgg(figure)
    labels = [f"{row['model']} · {row['seed']} · {row['run_id']}\n{row['scenario']} · {row['fold']}"
              for row in rows]
    colors = [palette.color("blue" if row["scenario"] == "cold-drug" else "orange", palette.FILL_SHADE)
              for row in rows]
    for axis, (metric, label) in zip(figure.subplots(2, 2).flat, LABELS.items(), strict=True):
        axis.bar(range(len(rows)), [float(row[metric]) if row[metric] else np.nan for row in rows],
                 color=colors)
        axis.set_xticks(range(len(rows)), labels, rotation=25, ha="right")
        axis.set_title(f"Δ {label}")
        axis.grid(axis="y", color=palette.color("gray", 8))
        axis.set_axisbelow(True)
    plotting.save_figure(figure, output, dpi=180)
