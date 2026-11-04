"""Real uM concentrations for the DToxS-primary arm (copied from Stage 1).

Provenance: ``stage1/xpert-fda-miri-v3/03_cm_finetune_data/
primary_concentrations.py`` at commit ``c921f7d`` (ADR 0015).  Stage 1 is a
read-only reference; this module is a copy, not an import.

The map covers the 46 primary drugs that Stage 1 resolved to an L1000 drug
index.  Stage 2 additionally needs ``DIC`` (diclofenac, 3 wells), which has no
DToxS small-molecule metadata row; its 5 uM concentration is recorded in the
DToxS DCIC experiment-design workbook and is handled in ``sources.py`` as a
``measured`` dose rather than a literature value.

Sources
-------
1. DToxS-primary paper Table 1 concentrations (23 kinase inhibitors).
2. GSE217423 iPSC-CM series-matrix concentrations (15 shared non-KIs).
3. Literature-curated clinical Cmax for the 8 primary-only drugs.
"""

PAPER_TABLE1_CONC = {
    "AFA": 0.05, "AXI": 0.2, "BOS": 0.1, "CAB": 2.0, "CER": 1.0,
    "CRI": 0.25, "DAB": 2.5, "DAS": 0.1, "ERL": 3.0, "GEF": 1.0,
    "IMA": 5.0, "LAP": 2.0, "NIL": 3.0, "PAZ": 10.0, "PON": 0.1,
    "REG": 1.0, "RUX": 1.0, "SOR": 0.5, "SUN": 1.0, "TOF": 1.0,
    "TRA": 0.1, "VAN": 0.33, "VEM": 2.0,
}
IPSC_SERIES_CONC = {
    "AMI": 0.5, "CYC": 0.125, "DAU": 1.0, "DOB": 1.0, "DOX": 1.0,
    "EPI": 2.0, "EST": 0.005, "FLE": 0.333, "IDA": 0.2, "ISO": 2.0,
    "MIL": 2.0, "OLM": 1.0, "PHP": 0.1, "PRE": 0.5, "VRP": 0.2,
}
PRIMARY_ONLY_CONC = {
    "ALE": 0.24, "DTP": 0.05, "DOM": 0.022, "ENT": 0.030,
    "LOP": 0.0008, "PAR": 0.15, "SB-431542": 10.0, "URS": 8.0,
}
PAPER_TABLE1_CITATION = (
    "DToxS-primary paper Table 1 'Concentration (uM)' (clinical max plasma "
    "concentrations): Sharma A, et al. Nat Commun 11:4809 (2020), "
    "DOI 10.1038/s41467-020-18396-7")
IPSC_SERIES_CITATION = (
    "GSE217423 iPSC-CM series matrix treatment concentrations "
    "(Sharma A, et al. Sci Transl Med 2022, GSE217423)")
PRIMARY_ONLY_CITATION = (
    "literature-curated clinical Cmax; see Stage 1 primary_concentrations.py "
    "for the per-drug citation (ADR 0015)")

_PRIMARY_NAMES = {
    "AFA": "Afatinib", "AXI": "Axitinib", "BOS": "Bosutinib",
    "CAB": "Cabozantinib", "CER": "Ceritinib", "CRI": "Crizotinib",
    "DAB": "Dabrafenib", "DAS": "Dasatinib", "ERL": "Erlotinib",
    "GEF": "Gefitinib", "IMA": "Imatinib", "LAP": "Lapatinib",
    "NIL": "Nilotinib", "PAZ": "Pazopanib", "PON": "Ponatinib",
    "REG": "Regorafenib", "RUX": "Ruxolitinib", "SOR": "Sorafenib",
    "SUN": "Sunitinib", "TOF": "Tofacitinib", "TRA": "Trametinib",
    "VAN": "Vandetanib", "VEM": "Vemurafenib",
    "AMI": "Amiodarone", "CYC": "Cyclosporine", "DAU": "Daunorubicin",
    "DOB": "Dobutamine", "DOX": "Doxorubicin", "EPI": "Epirubicin",
    "EST": "Estradiol", "FLE": "Flecainide", "IDA": "Idarubicin",
    "ISO": "Isoproterenol", "MIL": "Milrinone", "OLM": "Olmesartan",
    "PHP": "Phenylephrine", "PRE": "Prednisolone", "VRP": "Verapamil",
    "ALE": "Alendronic acid", "DTP": "Diethylpropion",
    "DOM": "Domperidone", "ENT": "Entecavir", "LOP": "Loperamide",
    "PAR": "Paroxetine", "SB-431542": "SB-431542", "URS": "Ursodiol",
}
CODE_ALIASES = {"TGF": "SB-431542", "SUN+": "SUN"}

PRIMARY_CONC: dict[str, tuple[str, float, str]] = {}
for _code, _conc in PAPER_TABLE1_CONC.items():
    PRIMARY_CONC[_code] = (_PRIMARY_NAMES[_code], _conc, PAPER_TABLE1_CITATION)
for _code, _conc in IPSC_SERIES_CONC.items():
    PRIMARY_CONC[_code] = (_PRIMARY_NAMES[_code], _conc, IPSC_SERIES_CITATION)
for _code, _conc in PRIMARY_ONLY_CONC.items():
    PRIMARY_CONC[_code] = (_PRIMARY_NAMES[_code], _conc, PRIMARY_ONLY_CITATION)


def primary_dose(compound_code: str) -> float | None:
    code = CODE_ALIASES.get(compound_code, compound_code)
    row = PRIMARY_CONC.get(code)
    return None if row is None else row[1]
