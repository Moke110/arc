"""
Raw whole-transcriptome CM dataset loaders for the Stage 2 v1 benchmark.

Each loader reads one dataset under `pert_data/` (read-only) and returns a
sample table plus a raw count matrix. No 978-landmark projection is applied
anywhere (ADR 0020): genes are kept in their native space and mapped to the
common universe later by `genes.py`.

Only compounds that carry **both** a resolvable SMILES and a real dose are
kept as treated samples; this yields the spec's ~79 drug-like compounds
(46 DToxS + 42 gse244740 − 9 shared).

Returns per dataset:
    meta:   DataFrame, columns
            well, dataset, drug, compound_code, smiles, dose_uM, time,
            is_control, plate
    counts: ndarray (n_samples, n_genes) raw counts
    gene_ids: list[str], length n_genes, in the dataset's native id space
"""

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import (DTOXS_CONTROL_STATES, INPUT_DIR, dataset_path)  # noqa: E402
from primary_concentrations import primary_dose  # noqa: E402

# gse244740 paper Table 1 "tested dose range" tops (uM); Con_i = top * 10^(-0.5*(8-i)).
GSE244740_TOP_DOSE = {
    "Acetaminophen": 4000.0, "Acetylsalicylic acid": 1000.0,
    "Acyclovir": 1000.0, "Amoxicillin": 300.0, "Buspirone": 10.0,
    "Enalapril": 100.0, "Furosemide": 500.0, "Sildenafil": 40.0,
    "Tolbutamide": 6000.0, "Amitriptyline": 60.0, "Atenolol": 250.0,
    "Bepridil": 30.0, "Cisapride": 10.0, "Digoxin": 0.1, "Diltiazem": 50.0,
    "Dobutamine": 100.0, "Dopamine": 20.0, "Epinephrine": 1.0,
    "Levosimendan": 50.0, "Lidocaine": 1250.0, "Nifedipine": 50.0,
    "Propranolol": 200.0, "Sotalol": 1250.0, "Amiodarone": 40.0,
    "Amphotericin B": 165.0, "Bortezomib": 50.0, "5-Fluorouracil": 500.0,
    "Rofecoxib": 25.0, "Rosiglitazone": 200.0, "Cyclophosphamide": 4000.0,
    "Dasatinib": 40.0, "Doxorubicin": 50.0, "Idarubicin": 20.0,
    "Imatinib": 75.0, "Isoproterenol": 10.0, "Lapatinib": 20.0,
    "Mitomycin C": 1000.0, "Mitoxantrone": 10.0, "Verapamil": 50.0,
    "Sunitinib": 75.0, "Clozapine": 50.0, "3_Azido_3_Deoxythymidine": 2000.0,
}
# Normalize the series-matrix token to the Table-1 spaced name.
_GSE_TO_TABLE = {
    "5_Fluorouracil": "5-Fluorouracil", "Mitomycin_C": "Mitomycin C",
    "Acetylsalicylic_acid": "Acetylsalicylic acid",
    "Amphotericin_B": "Amphotericin B",
    "3_Azido_3_Deoxythymidine": "3_Azido_3_Deoxythymidine",
}


def _canon(smiles):
    if not smiles:
        return None
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    return None if mol is None else Chem.MolToSmiles(mol, isomericSmiles=False)


def _dtox_code_meta():
    """DToxS small-molecule metadata: code -> (name, canonical SMILES)."""
    meta = pd.read_csv(dataset_path("dtoxs_primary") / "Metadata"
                       / "Small_Molecule_Metadata.txt", sep="\t")
    out = {}
    for _, r in meta.iterrows():
        code = r.get("sm_center_batch_id")
        if pd.isna(code) or pd.isna(r["sm_smiles_parent"]):
            continue
        out.setdefault(str(code), (str(r["sm_name"]), _canon(r["sm_smiles_parent"])))
    return out


# ---------------------------------------------------------------------------
# dtoxs_primary — LINCS LDS-1593 well-level Level-3 package
# ---------------------------------------------------------------------------
def load_dtoxs_primary():
    import glob

    level3 = dataset_path("dtoxs_primary") / "Data" / "Level3"
    parts, designs = [], []
    for counts_rel in sorted(glob.glob(
            str(level3 / "**" / "PromoCell-DGE-Read-Counts-*.tsv"), recursive=True)):
        c = pd.read_csv(counts_rel, sep="\t", index_col=0)
        d = pd.read_csv(counts_rel.replace("Read-Counts", "Experiment-Design"),
                        sep="\t")
        key = (d["State"].astype(str) + "." + d["Cell"].astype(str) + "."
               + d["Experiment"].astype(str) + "." + d["Dish"].astype(str) + "."
               + d["Plate"].astype(str) + "." + d["Well"].astype(str))
        well_to_key = dict(zip(d["Well"].astype(str), key))
        keep = [w for w in c.columns if str(w) in well_to_key]
        c = c.loc[:, keep]
        c.columns = [well_to_key[str(w)] for w in keep]
        parts.append(c)
        designs.append(d.assign(_key=key))
    raw = pd.concat(parts, axis=1)
    raw = raw.loc[:, ~raw.columns.duplicated()]
    raw = raw.groupby(raw.index).sum()  # duplicate gene rows -> sum
    design = pd.concat(designs, ignore_index=True)
    key_info = {k: (str(st), str(pl), str(wl)) for k, st, pl, wl in
                zip(design["_key"], design["State"], design["Plate"], design["Well"])}
    code_meta = _dtox_code_meta()
    rows = []
    for col in raw.columns:
        st, plate, well = key_info[col]
        name, smiles = code_meta.get(st, (st, None))
        is_ctl = st in DTOXS_CONTROL_STATES
        dose = 0.0 if is_ctl else primary_dose(st)
        if dose is None and not is_ctl:
            # compound outside the cited concentration map -> excluded below
            smiles = None
        rows.append({"well": f"dtoxs_primary::{col}", "dataset": "dtoxs_primary",
                     "drug": name, "compound_code": st, "smiles": smiles,
                     "dose_uM": dose, "time": 48.0, "is_control": is_ctl,
                     "plate": plate, "dose_level": 0})
    meta = pd.DataFrame(rows)
    counts = raw.loc[:, meta["well"].map(lambda w: w.split("::")[1])].values.T
    gene_ids = [str(g) for g in raw.index]
    return _finish(meta, counts.astype(float), gene_ids)


# ---------------------------------------------------------------------------
# dtoxs_ipsc — GSE217423 DGE arm
# ---------------------------------------------------------------------------
def load_dtoxs_ipsc():
    import glob

    root = dataset_path("dtoxs_ipsc")
    design = pd.read_csv(sorted(glob.glob(str(root / "GSE217423*Design*.tsv.gz")))[0],
                         sep="\t")
    counts = pd.read_csv(sorted(glob.glob(str(root / "GSE217423*Read-Counts*.tsv.gz")))[0],
                         sep="\t", index_col=0)
    key = (design["State"].astype(str) + "." + design["Cell"].astype(str) + "."
           + design["Experiment"].astype(str) + "." + design["Dish"].astype(str)
           + "." + design["Plate"].astype(str) + "." + design["Well"].astype(str))
    valid = key.isin(counts.columns)
    counts = counts.loc[:, key[valid].values]
    counts.columns = [str(w) for w in key[valid].values]
    code_meta = _dtox_code_meta()
    ipsc = _ipsc_code_meta()
    rows = []
    for col in counts.columns:
        st = col.split(".")[0]
        name, dose = ipsc.get(st, (st, None))
        _, smiles = code_meta.get(st, (name, None))
        is_ctl = st == "CTRL"
        rows.append({"well": f"dtoxs_ipsc::{col}", "dataset": "dtoxs_ipsc",
                     "drug": name, "compound_code": st, "smiles": smiles,
                     "dose_uM": 0.0 if is_ctl else dose, "time": 48.0,
                     "is_control": is_ctl, "plate": col.split(".")[4],
                     "dose_level": 0})
    meta = pd.DataFrame(rows)
    counts = counts.loc[:, meta["well"].map(lambda w: w.split("::")[1])]
    gene_ids = [str(g) for g in counts.index]
    return _finish(meta, counts.values.T.astype(float), gene_ids)


def _ipsc_code_meta():
    path = dataset_path("dtoxs_ipsc") / "GSE217423_series_matrix.txt.gz"
    text = gzip.open(path, "rt").read()
    titles = re.findall(r'"([A-Z0-9-]+) MSN[^"]*"', text)
    dn = re.findall(r'"drug name: ([^"]*)"', text)
    conc = re.findall(r'"concentration: ([^"]*)"', text)
    out = {}
    for code, name, c in zip(titles, dn, conc):
        out.setdefault(code, (name, _parse_conc_uM(c)))
    return out


def _parse_conc_uM(s):
    s = s.strip().lower().replace(" ", "")
    if s in ("0", "0um", "0nm"):
        return 0.0
    m = re.match(r"^([0-9.]+)(um|nm|mm)?", s)
    if not m:
        return None
    val, unit = float(m.group(1)), m.group(2) or "um"
    return val / 1000.0 if unit == "nm" else (val * 1000.0 if unit == "mm" else val)


# ---------------------------------------------------------------------------
# gse244740 — ScreenSeq
# ---------------------------------------------------------------------------
def load_gse244740_smiles():
    path = INPUT_DIR / "gse244740_smiles.tsv"
    df = pd.read_csv(path, sep="\t")
    return {r["drug_name"]: (r["smiles"] or None) for _, r in df.iterrows()}


def load_gse244740():
    path = dataset_path("gse244740") / "GSE244740_series_matrix.txt.gz"
    text = gzip.open(path, "rt").read()
    title_line = [l for l in text.splitlines() if l.startswith("!Sample_title")][0]
    titles = re.findall(r'"([^"]*)"', title_line)
    tline = [l for l in text.splitlines() if l.startswith("!Sample_characteristics_ch1")
             and "treatment" in l][0]
    treats = re.findall(r'"treatment: ([^"]*)"', tline)
    counts = pd.read_csv(dataset_path("gse244740")
                         / "GSE244740_processed_data_counts.txt.gz", sep="\t")
    idx = {c: i for i, c in enumerate(counts.columns)}
    smiles_map = load_gse244740_smiles()
    rows = []
    for title, tr in zip(titles, treats):
        m = re.match(r"(.+?), (Con\d+|NoTreatment), (?:Optimization|Screening) "
                     r"plate (\d+), Well (\w+) \[(S[^\]]+)\]", title)
        if not m:
            continue
        drug, conc, plate, well, wellid = m.groups()
        if wellid not in idx:
            continue
        is_ctl = conc in ("NoTreatment", "Con0") or drug == "DMSO"
        if is_ctl:
            dose = 0.0
        else:
            table = _GSE_TO_TABLE.get(drug, drug)
            top = GSE244740_TOP_DOSE.get(table)
            ci = int(conc[3:])
            dose = None if top is None else top * 10 ** (-0.5 * (8 - ci))
        rows.append({"well": f"gse244740::{wellid}", "dataset": "gse244740",
                     "drug": drug, "compound_code": "", "smiles": smiles_map.get(drug),
                     "dose_uM": dose, "time": 24.0, "is_control": is_ctl,
                     "plate": f"plate{plate}", "dose_level": 0 if is_ctl else ci})
    meta = pd.DataFrame(rows)
    cols = [w.split("::")[1] for w in meta["well"]]
    counts_m = counts[cols].values.T
    return _finish(meta, counts_m.astype(float), [str(e) for e in counts["Ensembl"]])


# ---------------------------------------------------------------------------
def _finish(meta, counts, gene_ids):
    """Drop treated wells lacking a SMILES or a real dose (controls kept)."""
    treated = ~meta["is_control"]
    bad = treated & (meta["smiles"].isna() | meta["dose_uM"].isna())
    if bad.any():
        meta = meta.loc[~bad].reset_index(drop=True)
        counts = counts[~bad.values]
    return meta, counts, gene_ids


LOADERS = {
    "dtoxs_primary": load_dtoxs_primary,
    "dtoxs_ipsc": load_dtoxs_ipsc,
    "gse244740": load_gse244740,
}


def load_dataset(name):
    return LOADERS[name]()
