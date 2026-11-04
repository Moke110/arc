"""Source loaders for the Stage 2 CM corpus.

Provenance
----------
The four raw-file readers are a **copy** of
``stage1/xpert-fda-miri-v3/03_cm_finetune_data/_datasets.py`` at commit
``c921f7d`` (the Stage 1 all-CM finetune loaders).  Stage 1 is read-only
reference and is never imported.  The copy carries the Stage 2 fixes:

* underscore / trailing-space compound names are normalized;
* DToxS code-only states resolve through the code table (``DIC`` -> diclofenac,
  ``SUN+`` -> sunitinib, ``AMI`` -> its own ZINC structure);
* antibodies / proteins / peptides are labelled instead of silently dropped;
* every row carries ``dose_um`` + ``dose_source``, ``vehicle`` and ``time_h``;
* ``temp_env`` (GSE262419) ``Blank`` wells are technical blanks, neither control nor
  perturbation;
* ``cell_type`` and ``donor`` are explicit covariates.

Ticket #15 added two row-level provenance columns without touching any resolved
value (ADR 0033; the corpus ``obs`` carries both):

* ``cell_line`` — the raw cell-line label.  For ``dtoxs_primary`` it is the
  ``Cell`` field of the composite well key (``A``/``B``/``D``/``E`` on plates
  1/2/5, ``PC-A``…``PC-E`` on plate 4), which is otherwise only recoverable by
  parsing ``well``; for ``dtoxs_ipsc`` it is the raw cell label
  (``MSN02-04R-CM`` / ``MSN09-04R-CM``).  ``screen_ctox`` / ``temp_env`` name no
  per-well line, so it is ``NA`` (same convention as ``donor``).
* ``n_raw_rows`` — how many raw design rows collapsed into one loaded well.
  ``dtoxs_primary`` has 956 design rows for 842 wells (36 keys reused across
  the experiment-5 sub-batches); every other source has one row per well.

Counts are raw (untransformed) and are returned as ``int32``.
"""

from __future__ import annotations

import glob
import gzip
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import paths
from .annotation import HgncIndex
from .genes import (SourceGeneMap, map_ensembl_space, map_probe_space,
                    map_symbol_space)
from .primary_doses import primary_dose
from .text import normalize_name

CONTROL_STATES = {"CTRL", "CTRL-DMSO", "CTRL-EtOH"}


def _dtoxs_key(design: pd.DataFrame) -> pd.Series:
    """Stable compound/well key shared by the DToxS design and counts files."""
    return (design["State"].astype(str) + "." + design["Cell"].astype(str)
            + "." + design["Experiment"].astype(str) + "."
            + design["Dish"].astype(str) + "." + design["Plate"].astype(str)
            + "." + design["Well"].astype(str))

# Non-small-molecule perturbagens in the DToxS panels.
DTOXS_BIOLOGICS = {
    "TRS": ("Trastuzumab", "antibody"),
    "END": ("Endothelin-1", "peptide"),
    "IGF": ("IGF-1", "protein"),
    "TNF": ("TNF-alpha", "protein"),
    "BEV": ("Bevacizumab", "antibody"),
    "CTX": ("Cetuximab", "antibody"),
    "RTX": ("Rituximab", "antibody"),
}

# The single diclofenac concentration recorded in the DToxS DCIC workbook
# (SR1 G5/G6 and SR2 C11: Drug.Conc.1 = 5.0 uM).  `sources.check_sources`
# re-reads the workbook and fails if this changes.
DIC_DOSE_UM = 5.0
DIC_CITATION = ("DToxS PromoCell-DGE-RNAseq-Experiment-Design.DCIC.xlsx, "
                "sheets SR1/SR2, Drug.Conc.1 = 5.0 uM (diclofenac)")

# Curated top concentration per compound from the GSE244740 paper (Table 1);
# `screen_ctox` derives its 8-point half-log series from this value.
SCREEN_CTOX_TOP_DOSE = {
    normalize_name(name): dose for name, dose in {
        "Acetaminophen": 4000.0, "Acetylsalicylic_acid": 1000.0,
        "Acyclovir": 1000.0, "Amoxicillin": 300.0, "Buspirone": 10.0,
        "Enalapril": 100.0, "Furosemide": 500.0, "Sildenafil": 40.0,
        "Tolbutamide": 6000.0, "Amitriptyline": 60.0, "Atenolol": 250.0,
        "Bepridil": 30.0, "Cisapride": 10.0, "Digoxin": 0.1, "Diltiazem": 50.0,
        "Dobutamine": 100.0, "Dopamine": 20.0, "Epinephrine": 1.0,
        "Levosimendan": 50.0, "Lidocaine": 1250.0, "Nifedipine": 50.0,
        "Propranolol": 200.0, "Sotalol": 1250.0, "Amiodarone": 40.0,
        "Amphotericin_B": 165.0, "Bortezomib": 50.0, "5_Fluorouracil": 500.0,
        "Rofecoxib": 25.0, "Rosiglitazone": 200.0, "Cyclophosphamide": 4000.0,
        "Dasatinib": 40.0, "Doxorubicin": 50.0, "Idarubicin": 20.0,
        "Imatinib": 75.0, "Isoproterenol": 10.0, "Lapatinib": 20.0,
        "Mitomycin_C": 1000.0, "Mitoxantrone": 10.0, "Verapamil": 50.0,
        "Sunitinib": 75.0, "Clozapine": 50.0,
        "3_Azido_3_Deoxythymidine": 2000.0,
    }.items()
}
SCREEN_CTOX_SUNITINIB_TOP = 75.0

# Stage 1 looked up these names *before* underscore normalization, so they
# silently fell back to the sunitinib top dose.  Stage 2 corrects them; the
# Stage 1 dose regression reports the correction instead of failing on it.
SCREEN_CTOX_STAGE1_UNDERSCORE_FIXES = {
    "Acetylsalicylic acid", "5 Fluorouracil", "Mitomycin C",
    "Amphotericin B", "3 Azido 3 Deoxythymidine",
}


class SourceError(RuntimeError):
    """A source file does not match the frozen build contract."""


@dataclass
class SourceData:
    dataset_id: str
    meta: pd.DataFrame
    counts: np.ndarray          # n_samples x n_source_genes, int32
    gene_ids: list[str]
    gene_map: SourceGeneMap
    mapped_counts: np.ndarray | None = None

    @property
    def n(self) -> int:
        return len(self.meta)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------
def _int_counts(counts: np.ndarray, dataset_id: str) -> np.ndarray:
    counts = np.asarray(counts, dtype=np.float64)
    if not np.isfinite(counts).all():
        raise SourceError(f"{dataset_id}: non-finite count values")
    if (counts < 0).any():
        raise SourceError(f"{dataset_id}: negative count values")
    rounded = np.rint(counts)
    if not np.allclose(counts, rounded):
        raise SourceError(f"{dataset_id}: non-integral raw counts")
    return rounded.astype(np.int32)


def _empty_frame(records: list[dict]) -> pd.DataFrame:
    meta = pd.DataFrame(records)
    ordered = [
        "well", "dataset_id", "compound_name", "compound_code", "dose_um",
        "dose_source", "time_h", "is_control", "is_blank", "vehicle", "plate",
        "cell_type", "donor", "cell_line", "n_raw_rows", "entity_type",
        "dataset_smiles", "dataset_inchikey", "annotation_note",
    ]
    for column in ordered:
        if column not in meta:
            meta[column] = 1 if column == "n_raw_rows" else ""
    meta["n_raw_rows"] = meta["n_raw_rows"].astype(np.int32)
    return meta[ordered]


# ---------------------------------------------------------------------------
# DToxS primary — LINCS LDS-1593 well-level package
# ---------------------------------------------------------------------------
def _load_dtoxs_lincs_raw() -> pd.DataFrame:
    """Union the per-batch read-count tables (symbol x well)."""
    parts = []
    for counts_path in sorted(glob.glob(
            str(paths.DTOXS_LEVEL3 / "**" / "PromoCell-DGE-Read-Counts-*.tsv"),
            recursive=True)):
        counts = pd.read_csv(counts_path, sep="\t", index_col=0)
        design_path = counts_path.replace("Read-Counts", "Experiment-Design")
        design = pd.read_csv(design_path, sep="\t")
        key = _dtoxs_key(design)
        well_to_key = dict(zip(design["Well"].astype(str), key))
        keep = [w for w in counts.columns if str(w) in well_to_key]
        counts = counts.loc[:, keep]
        counts.columns = [well_to_key[str(w)] for w in keep]
        parts.append(counts)
    if not parts:
        raise SourceError(f"no LINCS LDS-1593 read-count tables under "
                          f"{paths.DTOXS_LEVEL3}")
    raw = pd.concat(parts, axis=1)
    raw = raw.loc[:, ~raw.columns.duplicated()]
    raw = raw.groupby(raw.index).sum()
    raw.columns = raw.columns.astype(str)
    return raw


def _read_dtoxs_design() -> dict[str, tuple[str, str, str, str]]:
    out = {}
    for design_path in sorted(glob.glob(
            str(paths.DTOXS_LEVEL3 / "**"
                / "PromoCell-DGE-Experiment-Design-*.tsv"), recursive=True)):
        design = pd.read_csv(design_path, sep="\t")
        key = _dtoxs_key(design)
        for k, state, plate, well, cell in zip(
                key, design["State"], design["Plate"], design["Well"],
                design["Cell"]):
            out[k] = (str(state), str(plate), str(well), str(cell))
    return out


def _dtoxs_design_row_counts() -> dict[str, int]:
    """How many raw design rows collapsed into each composite well key.

    The eight experiment-design tables hold 956 rows for 842 wells: 36 keys are
    reused inside the experiment-5 sub-batches.  ``_read_dtoxs_design`` keeps
    one row per key (the loader's behaviour since ticket 03); this counter
    exposes the collapse so the corpus can record ``n_raw_rows``.
    """
    counts: dict[str, int] = {}
    for design_path in sorted(glob.glob(
            str(paths.DTOXS_LEVEL3 / "**"
                / "PromoCell-DGE-Experiment-Design-*.tsv"), recursive=True)):
        design = pd.read_csv(design_path, sep="\t")
        for key in _dtoxs_key(design):
            counts[key] = counts.get(key, 0) + 1
    return counts


def _dtoxs_small_molecules() -> dict[str, dict]:
    table = pd.read_csv(paths.DTOXS_SMALL_MOL, sep="\t", dtype=str)
    out = {}
    for row in table.to_dict("records"):
        code = str(row.get("sm_center_batch_id", "")).strip()
        if not code:
            continue
        out.setdefault(code, {
            "name": normalize_name(row.get("sm_name", "")),
            "smiles": row.get("sm_smiles_parent") or "",
            "inchikey": row.get("sm_inchi_key_parent") or "",
            "pubchem_cid": row.get("sm_pubchem_cid") or "",
        })
    return out


def _dtoxs_identity(code: str, small: dict[str, dict]) -> dict:
    if code in DTOXS_BIOLOGICS:
        name, entity = DTOXS_BIOLOGICS[code]
        return {"compound_name": name, "entity_type": entity, "smiles": "",
                "inchikey": "", "note": f"code {code}: {entity}"}
    if code == "SUN+":
        row = small.get("SUN", {})
        return {"compound_name": "Sunitinib", "entity_type": "small_molecule",
                "smiles": row.get("smiles", ""),
                "inchikey": row.get("inchikey", ""),
                "note": "code SUN+; same parent structure as SUN"}
    if code == "DIC":
        # No DToxS small-molecule metadata row; identity comes from the name
        # and the frozen external annotation.
        return {"compound_name": "Diclofenac", "entity_type": "small_molecule",
                "smiles": "", "inchikey": "",
                "note": f"code DIC; no DToxS metadata row; {DIC_CITATION}"}
    if code in small:
        row = small[code]
        # Structure wins over a provisional ZINC identifier: the row keeps the
        # source-visible name, never a drug-class relabel.
        return {"compound_name": row["name"], "entity_type": "small_molecule",
                "smiles": row["smiles"], "inchikey": row["inchikey"],
                "note": ""}
    raise SourceError(f"dtoxs: unresolved perturbation code {code!r}")


def _dtoxs_vehicle(state: str) -> tuple[str, bool]:
    if state == "CTRL-DMSO":
        return "DMSO", True
    if state == "CTRL-EtOH":
        return "EtOH", True
    if state == "CTRL":
        return "VEH", True
    return "DMSO", False


def _dcic_doses() -> pd.DataFrame:
    frames = []
    for sheet in ("SR1", "SR2", "SR4", "SR5"):
        table = pd.read_excel(paths.DTOXS_DCIC_XLSX, sheet_name=sheet)
        table["sheet"] = sheet
        frames.append(table)
    return pd.concat(frames, ignore_index=True)


def load_dtoxs_primary() -> SourceData:
    raw = _load_dtoxs_lincs_raw()
    design = _read_dtoxs_design()
    design_rows = _dtoxs_design_row_counts()
    small = _dtoxs_small_molecules()
    rows = []
    for column in raw.columns:
        if column not in design:
            raise SourceError(f"dtoxs_primary: well {column!r} has no design row")
        state, plate, _well, cell = design[column]
        vehicle, is_control = _dtoxs_vehicle(state)
        if is_control:
            identity = {"compound_name": vehicle, "entity_type": "control",
                        "smiles": "", "inchikey": "", "note": ""}
            dose, dose_source = 0.0, ""
        else:
            identity = _dtoxs_identity(state, small)
            if identity["entity_type"] != "small_molecule":
                # Antibodies / proteins / peptides are not dosed in the
                # transcriptomic design; the row says so rather than imputing.
                dose, dose_source = 0.0, ""
            else:
                dose = primary_dose(state)
                if dose is not None:
                    dose_source = "literature"
                elif state == "DIC":
                    dose, dose_source = DIC_DOSE_UM, "measured"
                else:
                    raise SourceError(
                        f"dtoxs_primary: no recorded dose for code {state!r}")
        entity = identity["entity_type"]
        name = identity["compound_name"]
        rows.append({
            "well": f"dtoxs_primary::{column}",
            "dataset_id": "dtoxs_primary",
            "compound_name": name,
            "compound_code": state,
            "dose_um": float(dose),
            "dose_source": dose_source,
            "time_h": 48.0,
            "is_control": is_control,
            "is_blank": False,
            "vehicle": vehicle,
            "plate": plate,
            "cell_type": "primary_cm",
            "donor": "NA",
            "cell_line": cell,
            "n_raw_rows": design_rows.get(column, 1),
            "entity_type": entity,
            "dataset_smiles": identity["smiles"] if not is_control else "",
            "dataset_inchikey": identity["inchikey"] if not is_control else "",
            "annotation_note": identity["note"],
        })
    meta = _empty_frame(rows)
    columns = [w.split("::", 1)[1] for w in meta["well"]]
    counts = raw.loc[:, columns].values.T
    return SourceData("dtoxs_primary", meta,
                      _int_counts(counts, "dtoxs_primary"),
                      list(raw.index), None)  # gene map attached later


# ---------------------------------------------------------------------------
# DToxS iPSC — GSE217423 DGE arm
# ---------------------------------------------------------------------------
def _ipsc_series_meta() -> dict[str, tuple[str, float]]:
    with gzip.open(paths.DTOXS_IPSC / "GSE217423_series_matrix.txt.gz", "rt") as fh:
        text = fh.read()
    codes = re.findall(r'"([A-Z0-9-]+) MSN[^"]*"', text)
    names = re.findall(r'"drug name: ([^"]*)"', text)
    concentrations = re.findall(r'"concentration: ([^"]*)"', text)
    out = {}
    for code, name, value in zip(codes, names, concentrations):
        out.setdefault(code, (normalize_name(name), _parse_conc_um(value)))
    return out


def _parse_conc_um(value: str) -> float:
    text = str(value).strip().lower().replace(" ", "")
    if text in ("0", "0um", "0nm"):
        return 0.0
    match = re.match(r"^([0-9.]+)(um|nm|mm)?$", text)
    if not match:
        match = re.match(r"^([0-9.]+)(um|nm|mm)?/", text)
    if not match:
        return 0.0
    amount, unit = float(match.group(1)), match.group(2) or "um"
    if unit == "nm":
        return amount / 1000.0
    if unit == "mm":
        return amount * 1000.0
    return amount


def load_dtoxs_ipsc() -> SourceData:
    design_path = sorted(glob.glob(
        str(paths.DTOXS_IPSC / "GSE217423*Design*.tsv.gz")))[0]
    counts_path = sorted(glob.glob(
        str(paths.DTOXS_IPSC / "GSE217423*Read-Counts*.tsv.gz")))[0]
    design = pd.read_csv(design_path, sep="\t")
    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    key = _dtoxs_key(design)
    valid = key.isin(counts.columns)
    counts = counts.loc[:, key[valid].values]
    counts.columns = [str(w) for w in key[valid].values]
    series = _ipsc_series_meta()
    small = _dtoxs_small_molecules()

    rows = []
    for column in counts.columns:
        code = column.split(".")[0]
        cell_line = column.split(".")[1]
        name, dose = series.get(code, (code, 0.0))
        is_control = code == "CTRL"
        if is_control:
            name, vehicle, entity = "DMSO", "DMSO", "control"
            dose = 0.0
        elif code in DTOXS_BIOLOGICS:
            name, entity = DTOXS_BIOLOGICS[code]
            vehicle = "DMSO"
        else:
            row = small.get(code, {})
            name = normalize_name(row.get("name") or name)
            vehicle = "DMSO"
            entity = "small_molecule"
        identity = small.get(code, {})
        donor_match = re.match(r"^(MSN\d+)", cell_line)
        rows.append({
            "well": f"dtoxs_ipsc::{column}",
            "dataset_id": "dtoxs_ipsc",
            "compound_name": name,
            "compound_code": code,
            "dose_um": float(dose),
            "dose_source": "series_matrix" if not is_control else "",
            "time_h": 48.0,
            "is_control": is_control,
            "is_blank": False,
            "vehicle": vehicle,
            "plate": column.split(".")[4],
            "cell_type": "ipsc_cm",
            "donor": donor_match.group(1) if donor_match else "NA",
            "cell_line": cell_line,
            "n_raw_rows": 1,
            "entity_type": entity,
            "dataset_smiles": identity.get("smiles", "") if entity == "small_molecule" else "",
            "dataset_inchikey": identity.get("inchikey", "") if entity == "small_molecule" else "",
            "annotation_note": "",
        })
    meta = _empty_frame(rows)
    ordered = [w.split("::", 1)[1] for w in meta["well"]]
    matrix = counts.loc[:, ordered].values.T
    return SourceData("dtoxs_ipsc", meta, _int_counts(matrix, "dtoxs_ipsc"),
                      list(counts.index), None)


# ---------------------------------------------------------------------------
# screen_ctox (GSE244740) — ScreenSeq hiPSC-CM
# ---------------------------------------------------------------------------
def load_screen_ctox() -> SourceData:
    with gzip.open(paths.SCREEN_CTOX_DIR / "GSE244740_series_matrix.txt.gz", "rt") as fh:
        text = fh.read()
    title_line = [line for line in text.splitlines()
                  if line.startswith("!Sample_title")][0]
    titles = re.findall(r'"([^"]*)"', title_line)
    treatment_line = [line for line in text.splitlines()
                      if line.startswith("!Sample_characteristics_ch1")
                      and "treatment" in line][0]
    treatments = re.findall(r'"treatment: ([^"]*)"', treatment_line)

    counts = pd.read_csv(
        paths.SCREEN_CTOX_DIR / "GSE244740_processed_data_counts.txt.gz", sep="\t")
    well_to_index = {column: i for i, column in enumerate(counts.columns)}
    rows = []
    for title, treatment in zip(titles, treatments):
        match = re.match(
            r"(.+?), (Con\d+|NoTreatment), (Optimization|Screening) "
            r"plate (\d+), Well (\w+) \[(S[^\]]+)\]", title)
        if not match:
            continue
        drug, concentration, assay, plate, _well, well_id = match.groups()
        if well_id not in well_to_index:
            continue
        drug = normalize_name(drug)
        if concentration in ("NoTreatment", "Con0") or drug == "DMSO":
            is_control = True
            dose, dose_source = 0.0, ""
            name = "NoTreatment" if concentration == "NoTreatment" else "DMSO"
            vehicle = name
        else:
            is_control = False
            index = int(concentration[3:])
            top = SCREEN_CTOX_TOP_DOSE.get(drug)
            if top is None:
                if drug == "Sunitinib":
                    top = SCREEN_CTOX_SUNITINIB_TOP
                else:
                    raise SourceError(
                        f"screen_ctox: no curated top dose for {drug!r}")
            dose = top * 10 ** (-0.5 * (8 - index))
            dose_source = "derived_top_dose"
            name = drug
            vehicle = "DMSO"
        rows.append({
            "well": f"screen_ctox::{well_id}",
            "dataset_id": "screen_ctox",
            "compound_name": name,
            "compound_code": "",
            "dose_um": float(dose),
            "dose_source": dose_source,
            "time_h": 24.0,
            "is_control": is_control,
            "is_blank": False,
            "vehicle": vehicle,
            # Optimization and Screening plates are numbered independently in
            # the GEO series, so the assay is part of the physical plate id.
            "plate": f"{assay.lower()}_plate{plate}",
            "cell_type": "ipsc_cm",
            "donor": "NA",
            "cell_line": "NA",
            "n_raw_rows": 1,
            "entity_type": "control" if is_control else "small_molecule",
            "dataset_smiles": "",
            "dataset_inchikey": "",
            "annotation_note": "",
        })
    meta = _empty_frame(rows)
    ordered = [w.split("::", 1)[1] for w in meta["well"]]
    matrix = counts[ordered].values.T
    return SourceData("screen_ctox", meta, _int_counts(matrix, "screen_ctox"),
                      list(counts["Ensembl"]), None)


# ---------------------------------------------------------------------------
# temp_env (GSE262419) — TempO-Seq hiPSC-CM
# ---------------------------------------------------------------------------
def _norm_well(value: str) -> str:
    return re.sub(r"^([A-Z]+)0*(\d+)$", r"\1\2", str(value))


def load_temp_env() -> SourceData:
    hash_path = paths.TEMP_ENV_DIR / "GSE262419_hash.csv.gz"
    hash_table = pd.read_csv(hash_path)
    hash_table["Well_ID"] = hash_table["Well_ID"].astype(str)
    plate_files = sorted(glob.glob(str(paths.TEMP_ENV_DIR / "GSE262419_Plate*.csv.gz")))
    gene_rows: dict[str, int] = {}
    for path in plate_files:
        frame = pd.read_csv(path)
        for gene in frame["Genes"].astype(str):
            gene_rows.setdefault(gene, len(gene_rows))
    gene_ids = list(gene_rows)
    matrix = np.zeros((len(gene_ids), 0), dtype=np.float64)
    sample_columns, sample_rows, plate_of_well = [], [], {}
    for path in plate_files:
        frame = pd.read_csv(path)
        frame.columns = [str(c) for c in frame.columns]
        plate_id = frame.columns[1].split("-")[0]
        block = np.zeros((len(gene_ids), frame.shape[1] - 1), dtype=np.float64)
        for i, gene in enumerate(frame["Genes"].astype(str)):
            block[gene_rows[gene], :] = frame.iloc[i, 1:].astype(float).values
        sample_columns.extend(frame.columns[1:])
        sample_rows.append(block)
        for column in frame.columns[1:]:
            plate_of_well[column] = plate_id
    matrix = np.concatenate(sample_rows, axis=1)
    well_info = {(int(row["Plate_ID"]), _norm_well(row["Well_ID"])): row
                 for _, row in hash_table.iterrows()}

    rows = []
    for well_id in sample_columns:
        plate_name, _, well_part = well_id.partition("-")
        short = _norm_well(well_part)
        info = well_info.get((int(plate_name.replace("Plate", "")), short))
        if info is None:
            continue
        source_name = normalize_name(info["Chemical_name"])
        concentration = float(info["Chemical_Concentration_uM"])
        if source_name in ("VEH", "Media"):
            is_control, is_blank = True, False
            name = "DMSO" if source_name == "VEH" else "Media"
            vehicle = "DMSO" if source_name == "VEH" else "Media"
            entity = "control"
            dose = 0.0
        elif source_name == "Blank":
            is_control, is_blank = False, True
            name, vehicle, entity, dose = "Blank", "DMSO", "blank", 0.0
        else:
            is_control, is_blank = False, False
            name = source_name
            vehicle = "DMSO"
            entity = "peptide" if source_name.upper() == "ET-1" else "small_molecule"
            dose = concentration
        rows.append({
            "well": f"temp_env::{well_id}",
            "dataset_id": "temp_env",
            "compound_name": name,
            "compound_code": str(info["Chemical_Index"]),
            "dose_um": dose,
            "dose_source": "measured" if not (is_control or is_blank) else "",
            "time_h": 24.0,
            "is_control": is_control,
            "is_blank": is_blank,
            "vehicle": vehicle,
            "plate": plate_of_well.get(well_id, ""),
            "cell_type": "ipsc_cm",
            "donor": "NA",
            "cell_line": "NA",
            "n_raw_rows": 1,
            "entity_type": entity,
            "dataset_smiles": "",
            "dataset_inchikey": "",
            "annotation_note": "",
        })
    meta = _empty_frame(rows)
    ordered = [w.split("::", 1)[1] for w in meta["well"]]
    matrix = matrix[:, [sample_columns.index(w) for w in ordered]].T
    return SourceData("temp_env", meta, _int_counts(matrix, "temp_env"),
                      gene_ids, None)


# ---------------------------------------------------------------------------
# Dispatch + gene mapping
# ---------------------------------------------------------------------------
LOADERS = {
    "dtoxs_primary": load_dtoxs_primary,
    "dtoxs_ipsc": load_dtoxs_ipsc,
    "screen_ctox": load_screen_ctox,
    "temp_env": load_temp_env,
}


def attach_gene_map(source: SourceData, index: HgncIndex,
                    entrez_map: dict[str, str] | None = None) -> SourceData:
    if source.dataset_id in ("dtoxs_primary", "dtoxs_ipsc"):
        source.gene_map = map_symbol_space(
            source.dataset_id, source.gene_ids, index, entrez_map)
    elif source.dataset_id == "screen_ctox":
        source.gene_map = map_ensembl_space(
            source.dataset_id, source.gene_ids, index)
    elif source.dataset_id == "temp_env":
        source.gene_map = map_probe_space(
            source.dataset_id, source.gene_ids, index)
    else:
        raise SourceError(f"unknown dataset {source.dataset_id!r}")
    return source


def load_source(dataset_id: str, index: HgncIndex,
                entrez_map: dict[str, str] | None = None) -> SourceData:
    source = LOADERS[dataset_id]()
    return attach_gene_map(source, index, entrez_map)


def load_dtoxs_entrez_map() -> dict[str, str]:
    table = pd.read_csv(paths.DTOXS_ENTREZ_MAP, sep="\t")
    out = {}
    for symbol, entrez in zip(table["HGNC.Symbol"].astype(str),
                              table["Entrez.ID"]):
        if pd.notna(entrez):
            out[symbol] = str(int(entrez))
    return out


# ---------------------------------------------------------------------------
# --check
# ---------------------------------------------------------------------------
def summarize(source: SourceData) -> dict[str, int]:
    meta = source.meta
    perturbed = ~meta["is_control"] & ~meta["is_blank"]
    return {
        "wells": len(meta),
        "controls": int(meta["is_control"].sum()),
        "blanks": int(meta["is_blank"].sum()),
        "perturbed": int(perturbed.sum()),
        "small_molecule_perturbed": int(
            (perturbed & (meta["entity_type"] == "small_molecule")).sum()),
        "distinct_compound_names": int(meta["compound_name"].nunique()),
    }


def check_sources() -> int:
    failures = 0
    for dataset_id in paths.SOURCE_DATASETS:
        source = LOADERS[dataset_id]()
        summary = summarize(source)
        expected = paths.EXPECTED_SOURCE_COUNTS[dataset_id]
        print(f"{dataset_id}: {summary}")
        for key, value in expected.items():
            actual = summary.get(key)
            if actual != value:
                print(f"  MISMATCH {key}: expected {value}, got {actual}")
                failures += 1
    # The DCIC workbook is the recorded source of the DIC dose.
    dcic = _dcic_doses()
    dic = dcic[dcic["State"].astype(str) == "DIC"]
    values = sorted(set(dic["Drug.Conc.1"].dropna().astype(float)))
    if values != [DIC_DOSE_UM]:
        print(f"  MISMATCH DIC dose: workbook {values}, frozen {DIC_DOSE_UM}")
        failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(check_sources())
