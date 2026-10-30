"""
Per-dataset sample extraction for the CM finetune-data union (step 03).

Each loader reads the raw files of one union dataset under pert_data
(read-only) and normalizes it into (sample metadata, raw counts matrix, gene
index). Gene identity is resolved separately by _lib/gene_mapping.py; drug
identity by _lib/drug_resolution.py. No pert_data file is ever written.

Identifier spaces per dataset (matching _lib/gene_mapping.ID_SPACE):
  dtoxs_primary / dtoxs_ipsc   gene symbols
  gse244740 / gse247541        Ensembl IDs (version suffix stripped)
  gse262419                    SYMBOL_entrezid rows

Output contract of `extract_samples(dataset)`:
  (meta_df, counts, gene_ids)
    meta_df.columns: well, dataset_id, compound_name, compound_code, smiles,
                     dose, time, is_control, plate, cell_iname
    counts: ndarray (n_genes, n_samples) raw counts
    gene_ids: list of str (length n_genes) in the dataset's identifier space

Dose conventions (all in uM, log2(CPM+1) space is applied in build_union):
  dtoxs_primary                 real per-drug concentrations from the
                                primary_concentrations map (paper Table 1 KIs,
                                iPSC series matrix shared non-KIs, literature
                                Cmax primary-only drugs; ADR 0015) — the LINCS
                                LDS-1593 package records no dose; transcriptomes
                                from pert_data/LINCS-LDS-1593/source/Data/Level3
  dtoxs_ipsc                    per-drug concentration from the GSE217423
                                series matrix (uM), time 48 h
  gse244740                    8-point half-log series from the paper Table-1
                                top dose (Con_i = top * 10^(-0.5*(8-i))), 24 h
                                (finetune dataset in v3, ADR 0017)
  gse262419                    hash.csv Chemical_Concentration_uM, 24 h
                               (finetune dataset in v3, ADR 0017)
  gse247541                    single dose 1 uM (paper, clinically relevant),
                               144 h (6 d; excluded from the v3 finetune —
                                outside the author time map, 2 in-HG drugs)
Controls:
  dtoxs arms: CTRL / CTRL-DMSO / CTRL-EtOH (per-plate: plates 1/2/4/5 primary,
             plate 6 iPSC)
  gse244740:  DMSO (Con0) + NoTreatment (per Optimization-plate)
  gse262419:  VEH + Media (per-plate)
  gse247541:  DMSO_control (dataset-wide)
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_lib"))
from gene_mapping import load_dtoxs_lincs_raw  # noqa: E402
from drug_resolution import dtoxs_code_meta  # noqa: E402
from utils import CM_DATA_ROOT, DTOXS_LINCS_LEVEL3  # noqa: E402

# ---------------------------------------------------------------------------
# gse244740 paper Table 1 "tested dose range" tops (uM) — curated from the
# published Table 1 (Rosell-Hidalgo, Expert Opin Drug Metab Toxicol 2024).
# Only the top concentration is recorded; Con_i = top * 10^(-0.5*(8 - i))
# follows the paper's 8-point half-log dilution series.
# ---------------------------------------------------------------------------
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
    "Sunitinib": 75.0, "Clozapine": 50.0,
    "3_Azido_3_Deoxythymidine": 2000.0,
}
# Sunitinib is the positive control and is *not* in the 42-compound table; its
# tested range 0.03-75 uM is given in the same Table 1.
GSE244740_SUNITINIB_TOP = 75.0


def _read_csv(path, **kw):
    """pd.read_csv that transparently handles .gz files."""
    return pd.read_csv(path, **kw)


# ---------------------------------------------------------------------------
# dtoxs_primary — LINCS LDS-1593 well-level package (Level 3)
# ---------------------------------------------------------------------------
def extract_dtoxs_primary():
    raw, _ = load_dtoxs_lincs_raw()
    # The Experiment-5 Cor / Spec-Fixed / Spec-Var designs share control wells;
    # dedupe the raw columns so each physical well appears once.
    raw = raw.loc[:, ~raw.columns.duplicated()]
    code_meta = dtoxs_code_meta()
    # design map: compound column key -> (State, Plate, Well)
    design = _read_dtoxs_primary_design()
    rows = []
    for col in raw.columns:
        st, plate, well = design[col]
        name, smiles = code_meta.get(st, (st, None))
        rows.append({
            "well": f"dtoxs_primary::{col}",
            "dataset_id": "dtoxs_primary",
            "compound_name": name, "compound_code": st, "smiles": smiles,
            "dose": 0.0, "time": 48.0,
            "is_control": st in ("CTRL", "CTRL-DMSO", "CTRL-EtOH"),
            "plate": str(plate), "cell_iname": "cardio",
        })
    meta = pd.DataFrame(rows)
    counts = raw.loc[:, meta["well"].map(lambda w: w.split("::")[1])].values.T
    return meta, counts.astype(float), list(raw.index)


def _read_dtoxs_primary_design():
    import glob

    out = {}
    for design_rel in sorted(glob.glob(
            str(DTOXS_LINCS_LEVEL3 / "**" / "PromoCell-DGE-Experiment-Design-*.tsv"),
            recursive=True)):
        d = _read_csv(design_rel, sep="\t")
        key = (d["State"].astype(str) + "." + d["Cell"].astype(str) + "."
               + d["Experiment"].astype(str) + "." + d["Dish"].astype(str)
               + "." + d["Plate"].astype(str) + "." + d["Well"].astype(str))
        for k, st, pl, wl in zip(key, d["State"], d["Plate"], d["Well"]):
            out[k] = (st, pl, wl)
    return out


# ---------------------------------------------------------------------------
# dtoxs_ipsc — GSE217423 DGE arm
# ---------------------------------------------------------------------------
def extract_dtoxs_ipsc():
    import glob

    root = CM_DATA_ROOT / "DToxS-iPSC" / "source"
    design_path = sorted(glob.glob(str(root / "GSE217423*Design*.tsv.gz")))[0]
    counts_path = sorted(glob.glob(str(root / "GSE217423*Read-Counts*.tsv.gz")))[0]
    design = _read_csv(design_path, sep="\t")
    counts = _read_csv(counts_path, sep="\t", index_col=0)
    key = (design["State"].astype(str) + "." + design["Cell"].astype(str)
           + "." + design["Experiment"].astype(str) + "."
           + design["Dish"].astype(str) + "." + design["Plate"].astype(str)
           + "." + design["Well"].astype(str))
    valid = key.isin(counts.columns)
    counts = counts.loc[:, key[valid].values]
    counts.columns = [str(w) for w in key[valid].values]
    gene_ids = list(counts.index)
    code_meta = dtoxs_code_meta()
    ipsc_meta = _ipsc_code_meta()
    rows = []
    for col in counts.columns:
        st = col.split(".")[0]
        name, dose = ipsc_meta.get(st, (st, 0.0))
        _, smiles = code_meta.get(st, (name, None))
        rows.append({
            "well": f"dtoxs_ipsc::{col}",
            "dataset_id": "dtoxs_ipsc",
            "compound_name": name, "compound_code": st, "smiles": smiles,
            "dose": dose, "time": 48.0,
            "is_control": st == "CTRL",
            "plate": col.split(".")[4], "cell_iname": "cardio",
        })
    meta = pd.DataFrame(rows)
    counts = counts.loc[:, meta["well"].map(lambda w: w.split("::")[1])].values.T
    return meta, counts.astype(float), gene_ids


def _ipsc_code_meta():
    """code -> (drug name, dose uM) from the GSE217423 series matrix."""
    import gzip

    path = CM_DATA_ROOT / "DToxS-iPSC" / "source" / "GSE217423_series_matrix.txt.gz"
    with gzip.open(path, "rt") as f:
        text = f.read()
    titles = re.findall(r'\"([A-Z0-9-]+) MSN[^\"]*\"', text)
    dn = re.findall(r'\"drug name: ([^\"]*)\"', text)
    conc = re.findall(r'\"concentration: ([^\"]*)\"', text)
    out = {}
    for code, name, c in zip(titles, dn, conc):
        out.setdefault(code, (name, _parse_conc_uM(c)))
    return out


def _parse_conc_uM(s):
    s = s.strip().lower().replace(" ", "")
    if s in ("0", "0um", "0nm"):
        return 0.0
    m = re.match(r"^([0-9.]+)(um|nm|mm)?$", s)
    if not m:
        m = re.match(r"^([0-9.]+)(um|nm|mm)?/", s)  # '3uM /10uM'
    if not m:
        return 0.0
    val, unit = float(m.group(1)), m.group(2) or "um"
    return val / 1000.0 if unit == "nm" else (val * 1000.0 if unit == "mm" else val)


# ---------------------------------------------------------------------------
# gse244740 — ScreenSeq
# ---------------------------------------------------------------------------
def extract_gse244740():
    import gzip

    import scanpy as sc

    path = CM_DATA_ROOT / "GSE244740" / "source" / "GSE244740_series_matrix.txt.gz"
    with gzip.open(path, "rt") as f:
        text = f.read()
    title_line = [l for l in text.splitlines() if l.startswith("!Sample_title")][0]
    titles = re.findall(r'\"([^\"]*)\"', title_line)
    tline = [l for l in text.splitlines() if l.startswith("!Sample_characteristics_ch1")
             and "treatment" in l][0]
    treats = re.findall(r'\"treatment: ([^\"]*)\"', tline)

    counts = _read_csv(CM_DATA_ROOT / "GSE244740" / "source"
                       / "GSE244740_processed_data_counts.txt.gz",
                       sep="\t")
    well_to_idx = {c: i for i, c in enumerate(counts.columns)}

    rows = []
    for title, tr in zip(titles, treats):
        m = re.match(r"(.+?), (Con\d+|NoTreatment), (?:Optimization|Screening) "
                     r"plate (\d+), Well (\w+) \[(S[^\]]+)\]", title)
        if not m:
            continue
        drug, conc, plate, well, wellid = m.groups()
        if wellid not in well_to_idx:
            continue
        if conc == "NoTreatment" or conc == "Con0" or drug == "DMSO":
            is_ctl, dose = True, 0.0
        else:
            is_ctl = False
            ci = int(conc[3:])
            top = (GSE244740_TOP_DOSE.get(drug)
                   if drug in GSE244740_TOP_DOSE
                   else GSE244740_SUNITINIB_TOP)
            dose = top * 10 ** (-0.5 * (8 - ci))
        rows.append({
            "well": f"gse244740::{wellid}",
            "dataset_id": "gse244740",
            "compound_name": drug, "compound_code": "", "smiles": None,
            "dose": dose, "time": 24.0, "is_control": is_ctl,
            "plate": f"plate{plate}", "cell_iname": "cardio",
        })
    meta = pd.DataFrame(rows)
    cols = [w.split("::")[1] for w in meta["well"]]
    counts_m = counts[cols].values.T
    return meta, counts_m.astype(float), list(counts["Ensembl"])


# ---------------------------------------------------------------------------
# gse247541 — bulk RNA-seq (4 antivirals)
# ---------------------------------------------------------------------------
def extract_gse247541():
    counts = _read_csv(CM_DATA_ROOT / "GSE247541" / "source" / "GSE247541_raw_counts.csv.gz",
                       index_col=0)
    rows = []
    for col in counts.columns:
        drug = col.rsplit("_", 1)[0]
        is_ctl = drug == "DMSO_control" or col.startswith("DMSO")
        rows.append({
            "well": f"gse247541::{col}",
            "dataset_id": "gse247541",
            "compound_name": "DMSO" if is_ctl else drug, "compound_code": "",
            "smiles": None, "dose": 0.0 if is_ctl else 1.0, "time": 144.0,
            "is_control": is_ctl, "plate": "", "cell_iname": "cardio",
        })
    meta = pd.DataFrame(rows)
    counts_m = counts.loc[:, meta["well"].map(lambda w: w.split("::")[1])].values.T
    gene_ids = [str(e).split(".")[0] for e in counts.index]
    return meta, counts_m.astype(float), gene_ids


# ---------------------------------------------------------------------------
# gse262419 — TempO-Seq (16 plates + hash.csv)
# ---------------------------------------------------------------------------
def _norm_well(s):
    """Normalize a plate well id: 'B02' -> 'B2' (strip leading zeros)."""
    return re.sub(r"^([A-Z]+)0*(\d+)$", r"\1\2", str(s))


def extract_gse262419():
    import glob

    root = CM_DATA_ROOT / "GSE262419" / "source"
    hashdf = _read_csv(root / "GSE262419_hash.csv.gz")
    hashdf["Well_ID"] = hashdf["Well_ID"].astype(str)
    plate_files = sorted(glob.glob(str(root / "GSE262419_Plate*.csv.gz")))
    gene_rows = {}
    for f in plate_files:
        df = _read_csv(f)
        for g in df["Genes"].astype(str):
            gene_rows.setdefault(g, len(gene_rows))
    gene_ids = list(gene_rows)
    n_genes = len(gene_ids)
    sample_cols, sample_rows = [], []
    plate_of_well = {}
    for f in plate_files:
        df = _read_csv(f)
        df.columns = [str(c) for c in df.columns]
        plate_id = df.columns[1].split("-")[0]
        mat = np.zeros((n_genes, df.shape[1] - 1), dtype=float)
        for i, g in enumerate(df["Genes"].astype(str)):
            mat[gene_rows[g], :] = df.iloc[i, 1:].astype(float).values
        sample_cols.extend(df.columns[1:])
        sample_rows.append(mat)
        for c in df.columns[1:]:
            plate_of_well[c] = plate_id
    counts = np.concatenate(sample_rows, axis=1)
    well2info = {(int(r["Plate_ID"]), _norm_well(r["Well_ID"])): r
                 for _, r in hashdf.iterrows()}
    rows = []
    for wellid in sample_cols:
        plate_name, _, well_part = wellid.partition("-")
        short = _norm_well(well_part)
        info = well2info.get((int(plate_name.replace("Plate", "")), short))
        if info is None:
            continue
        name = str(info["Chemical_name"])
        is_ctl = name in ("VEH", "Media")
        rows.append({
            "well": f"gse262419::{wellid}",
            "dataset_id": "gse262419",
            "compound_name": "DMSO" if name == "VEH" else name,
            "compound_code": str(info["Chemical_Index"]), "smiles": None,
            "dose": float(info["Chemical_Concentration_uM"]),
            "time": 24.0, "is_control": is_ctl,
            "plate": plate_of_well.get(wellid, ""), "cell_iname": "cardio",
        })
    meta = pd.DataFrame(rows)
    # keep only wells that survived the meta join; counts is (n_genes, n_samples)
    well_order = [w.split("::")[1] for w in meta["well"]]
    counts = counts[:, [sample_cols.index(w) for w in well_order]].T
    return meta, counts, gene_ids


def extract_samples(dataset):
    """Dispatch to the per-dataset loader."""
    return {
        "dtoxs_primary": extract_dtoxs_primary,
        "dtoxs_ipsc": extract_dtoxs_ipsc,
        "gse244740": extract_gse244740,
        "gse247541": extract_gse247541,
        "gse262419": extract_gse262419,
    }[dataset]()
