"""Import native cardiac count tables in bounded gene blocks."""
from __future__ import annotations

import csv
import gzip
import re

import h5py
import numpy as np
import pandas as pd
from anndata.io import write_elem

from native_corpus import CorpusError

GENE_BATCH = 256


def _unique_axis(values):
    counts = pd.Series(values).value_counts()
    return set(counts[counts == 1].index)


def _coverage(gene_ids, namespace, axis):
    ids = pd.Index(gene_ids)
    if namespace == "ensembl_id":
        ids = ids.str.replace(r"\.\d+$", "", regex=True)
    unique = _unique_axis(ids)
    lookup = {gene: col for col, gene in enumerate(ids) if gene in unique}
    return {int(index): lookup[gene] for index, gene in axis[namespace].items()
            if gene in lookup}


def build_count_table(candidate, name, path, metadata, axis, namespace, sep="\t"):
    """Stage one immutable count table and return its measured common-axis coverage."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as header_file:
        header = next(csv.reader(header_file, delimiter=sep), [])
    observations = header[1:]
    if not observations or any(not name for name in observations):
        raise CorpusError(f"{name}: missing source observation identifiers")
    if len(set(observations)) != len(observations):
        raise CorpusError(f"{name}: duplicate source observation identifiers")
    gene_ids = []
    destination = candidate / f"{name}.h5ad"
    with (h5py.File(destination, "w") as handle,
          pd.read_csv(path, sep=sep, index_col=0, chunksize=GENE_BATCH) as frames):
        handle.attrs.update({"encoding-type": "anndata", "encoding-version": "0.1.0"})
        for frame in frames:
            if len(frame.columns) == 0:
                raise CorpusError(f"{name}: no source observations")
            if not gene_ids:
                obs = metadata(frame.columns)
                write_elem(handle, "obs", obs)
                matrix = handle.create_dataset(
                    "X", (len(obs), 0), maxshape=(len(obs), None), dtype="i4",
                    chunks=(min(64, len(obs)), GENE_BATCH),
                    compression="gzip", compression_opts=1, shuffle=True,
                )
                matrix.attrs.update({"encoding-type": "array", "encoding-version": "0.2.0"})
            values = frame.to_numpy()
            if values.dtype.kind not in "iu" or np.any(values < 0) or np.any(
                values > np.iinfo(np.int32).max
            ):
                raise CorpusError(f"{name}: source count is not a nonnegative int32")
            end = len(gene_ids) + len(frame)
            matrix.resize(end, axis=1)
            matrix[:, len(gene_ids):end] = values.T
            if frame.index.hasnans:
                raise CorpusError(f"{name}: empty native gene identifier")
            gene_ids.extend(frame.index.astype(str))
        if not gene_ids:
            raise CorpusError(f"{name}: no source genes")
        # Repeated publisher identifiers remain separate native columns.
        var = pd.DataFrame({"source_gene_id": gene_ids},
                           index=pd.Index([f"column-{i}" for i in range(len(gene_ids))],
                                          name="native_id"))
        write_elem(handle, "var", var)
    return _coverage(gene_ids, namespace, axis)


def _dtoxs_metadata(design_path, plates_are_separate):
    records = pd.read_csv(design_path, sep="\t", dtype=str, keep_default_na=False)
    if plates_are_separate:
        keys = records.Well
    else:
        fields = ["State", "Cell", "Experiment", "Dish", "Plate", "Well"]
        keys = records[fields].agg(".".join, axis=1)
    if keys.duplicated().any():
        raise CorpusError(f"{design_path.name}: repeated source observation key")
    records = records.set_index(keys)

    def align(columns):
        obs = records.reindex(columns).fillna("")
        obs.index = pd.Index(columns, name="profile_id")
        obs["role"] = np.where(obs.State.eq(""), "unknown", "treated")
        obs.loc[obs.State.isin(["CTRL", "CTRL-DMSO", "CTRL-EtOH"]), "role"] = "vehicle"
        obs["compound_code"] = obs.State
        # Distinct plates, cell lines and exposure times never share vehicles.
        obs["vehicle"] = obs.State.map({"CTRL": "DMSO" if not plates_are_separate else "",
                                       "CTRL-DMSO": "DMSO", "CTRL-EtOH": "EtOH"}).fillna("DMSO")
        obs["matched_control_group"] = ""
        resolved = obs.role.ne("unknown") & obs.vehicle.ne("")
        obs.loc[resolved, "matched_control_group"] = obs.loc[
            resolved, ["Plate", "Cell", "Time", "vehicle"]
        ].agg("::".join, axis=1)
        return obs

    return align


def _series_titles(path):
    with gzip.open(path, "rt", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if row and row[0] == "!Sample_title":
                return row[1:]
    raise CorpusError(f"{path.name}: series deposit has no sample titles")


def _screen_metadata(series):
    pattern = re.compile(
        r"(.+?), (Con\d+|NoTreatment), (Optimization|Screening) "
        r"plate (\d+), Well (\w+) \[(S[^\]]+)\]"
    )
    rows = []
    for title in _series_titles(series):
        match = pattern.fullmatch(title)
        if match is None:
            continue
        drug, concentration, assay, plate, well, identity = match.groups()
        vehicle = "NoTreatment" if concentration == "NoTreatment" else "DMSO"
        control = drug == "DMSO" or concentration in {"NoTreatment", "Con0"}
        rows.append({
            "profile_id": identity, "compound_name": drug, "dose_code": concentration,
            "plate": f"{assay} {plate}", "well": well, "vehicle": vehicle,
            "role": "vehicle" if control else "treated",
            "matched_control_group": f"{assay}::{plate}::{vehicle}",
        })
    records = pd.DataFrame(rows).set_index("profile_id")

    def align(columns):
        obs = records.reindex(columns).fillna("")
        obs.index = pd.Index(columns, name="profile_id")
        obs.loc[obs.role.eq(""), "role"] = "unknown"
        return obs

    return align


def _antiviral_metadata(series):
    pattern = re.compile(r"hCM, (.*), replicate (\d+), RNA-seq")
    rows = []
    for title in _series_titles(series):
        match = pattern.fullmatch(title)
        if match is None:
            raise CorpusError(f"unexpected GSE247541 sample title {title!r}")
        drug, replicate = match.groups()
        control = drug == "DMSO"
        identity = f"{'DMSO_control' if control else drug}_rep{replicate}"
        rows.append({
            "profile_id": identity, "compound_name": drug, "replicate": int(replicate),
            "role": "vehicle" if control else "treated",
            "matched_control_group": "GSE247541::hCM::DMSO",
        })
    records = pd.DataFrame(rows).set_index("profile_id")

    def align(columns):
        obs = records.loc[columns].copy()
        obs.index.name = "profile_id"
        return obs

    return align


def build_cardiac_sources(candidate, data, axis, pin_source, inputs):
    mappings = {}
    consumed = {}

    def register(name, item, table, design, metadata, namespace, sep="\t"):
        for path in (table, design):
            consumed[path.relative_to(data).as_posix()] = pin_source(data / item, path)
        mappings[name] = build_count_table(candidate, name, table, metadata(design),
                                axis, namespace, sep)

    for name, series, technology in (
        ("dtoxs_ipsc_conv", "GSE217421", "Conv"),
        ("dtoxs_ipsc_dge", "GSE217423", "DGE"),
    ):
        item = "DToxS-iPSC"
        source = data / item / "source"
        prefix = f"{series}_{technology}-RNAseq"
        register(name, item, source / f"{prefix}-Read-Counts.tsv.gz",
                 source / f"{prefix}-Experiment-Design.tsv.gz",
                 lambda path: _dtoxs_metadata(path, plates_are_separate=False), "symbol")

    item = "LINCS-LDS-1593"
    source = data / item / "source"
    tables = sorted((source / "Data/Level3").glob("**/*Read-Counts*.tsv"))
    if not tables:
        raise CorpusError("LINCS-LDS-1593 has no deposited read-count tables")
    for index, table in enumerate(tables, 1):
        design = table.with_name(table.name.replace("Read-Counts", "Experiment-Design"))
        register(f"dtoxs_primary_{index:02}", item, table, design,
                 lambda path: _dtoxs_metadata(path, plates_are_separate=True), "symbol")

    item = "GSE244740"
    source = data / item / "source"
    register("screen_ctox", item, source / "GSE244740_processed_data_counts.txt.gz",
             source / "GSE244740_series_matrix.txt.gz", _screen_metadata, "ensembl_id")

    item = "GSE247541"
    source = data / item / "source"
    register("gse247541", item, source / "GSE247541_raw_counts.csv.gz",
             source / "GSE247541_series_matrix.txt.gz", _antiviral_metadata, "ensembl_id", ",")

    inputs["additional_count_files"] = consumed
    return mappings
