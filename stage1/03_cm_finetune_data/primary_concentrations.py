#!/usr/bin/env python3
"""
Primary finetune concentration map (ADR 0015, step 03).

Every `dtoxs_primary` finetune drug is mapped to the real uM concentration its
transcriptomes were treated at, so the union finetune's dose bins reflect
actual clinical exposure instead of v1's single-dose (dose 0 / bin 0)
convention. Sources:

  1. the DToxS-primary paper Table 1 "Concentration (uM)" column — the 23
     kinase inhibitors (derived from clinically reported maximum plasma
     concentrations; the paper's methods: "drug exposures similar to
     clinically reported maximum plasma concentrations", not corrected for
     protein binding);
  2. the iPSC-CM series matrix (GSE217423) — the 15 non-KI drugs shared with
     the primary arm carry the iPSC study's treatment concentration;
  3. literature-curated clinical Cmax with citations — the 8 primary-only
     drugs. SB-431542 has no clinical Cmax and uses its standard in-vitro
     assay concentration (10 uM).

Every primary in-HG drug (46 compounds) has a concentration and a cited
source; no treated sample falls back to a placeholder dose. The intended SOR
divergence is documented: primary SOR 0.5 uM (paper Table 1) vs iPSC SOR
1.0 uM (series matrix). The bin assignment uses the author
`assign_pert_dose_corrected` bins (0-9) mirrored in `build_union.dose_idx_of`.

Usage (imported by build_union.py; no CLI):
    from primary_concentrations import PRIMARY_CONC, primary_dose
"""

# ---------------------------------------------------------------------------
# 1) The 23 KIs — DToxS-primary paper Table 1, concentration (uM) column.
#    Sharma et al., Nat Commun 11:4809 (2020), DOI 10.1038/s41467-020-18396-7.
#    Values are the clinically reported max total plasma concentrations the
#    study treated the primary CM cells at (SOR 0.5 uM, RUX 1.0 uM confirmed).
# ---------------------------------------------------------------------------
PAPER_TABLE1_CONC = {
    "AFA": 0.05, "AXI": 0.2, "BOS": 0.1, "CAB": 2.0, "CER": 1.0,
    "CRI": 0.25, "DAB": 2.5, "DAS": 0.1, "ERL": 3.0, "GEF": 1.0,
    "IMA": 5.0, "LAP": 2.0, "NIL": 3.0, "PAZ": 10.0, "PON": 0.1,
    "REG": 1.0, "RUX": 1.0, "SOR": 0.5, "SUN": 1.0, "TOF": 1.0,
    "TRA": 0.1, "VAN": 0.33, "VEM": 2.0,
}

# ---------------------------------------------------------------------------
# 2) The 15 shared non-KIs — GSE217423 iPSC-CM series matrix treatment
#    concentrations (the primary arm carries the same drugs; the primary
#    platform recorded no concentration, so the iPSC study's value is used).
# ---------------------------------------------------------------------------
IPSC_SERIES_CONC = {
    "AMI": 0.5,       # Amiodarone      0.5 uM
    "CYC": 0.125,     # Cyclosporine    125 nM
    "DAU": 1.0,       # Daunorubicin    1 uM
    "DOB": 1.0,       # Dobutamine      1 uM
    "DOX": 1.0,       # Doxorubicin     1 uM
    "EPI": 2.0,       # Epirubicin      2 uM
    "EST": 0.005,     # Estradiol       5 nM
    "FLE": 0.333,     # Flecainide      0.333 uM
    "IDA": 0.2,       # Idarubicin      0.2 uM
    "ISO": 2.0,       # Isoproterenol   2 uM
    "MIL": 2.0,       # Milrinone       2 uM
    "OLM": 1.0,       # Olmesartan      1 uM
    "PHP": 0.1,       # Phenylephrine   100 nM
    "PRE": 0.5,       # Prednisolone    0.5 uM
    "VRP": 0.2,       # Verapamil       0.2 uM
}

# ---------------------------------------------------------------------------
# 3) The 8 primary-only drugs — literature-curated clinical Cmax (uM), each
#    with a verifiable citation. SB-431542 has no clinical Cmax and uses its
#    standard in-vitro assay concentration (10 uM).
# ---------------------------------------------------------------------------
PRIMARY_ONLY_CONC = {
    "ALE": (0.24,
            "clinical Cmax (alendronic acid 70 mg oral, mean Cmax ~60 ng/mL)",
            "Zhang Y, et al. Drug Des Devel Ther 2017;11:2109-2119 "
            "(PMID 28744102)"),
    "DTP": (0.05,
            "approximate steady-state plateau estimate (75 mg SR; no discrete "
            "numeric Cmax published; plasma exposure is metabolite-dominated)",
            "Dangor CM, et al. Arzneimittelforschung 1987;37(6):736-739 "
            "(PMID 3663274)"),
    "DOM": (0.022,
            "clinical Cmax (domperidone 10 mg oral, ~9.4 ng/mL, dose-proportional "
            "from the 20 mg value)",
            "Huang YC, et al. J Clin Pharmacol 1986;26(8):628-632 "
            "(PMID 3793955)"),
    "ENT": (0.030,
            "clinical steady-state Cmax (entecavir 1 mg oral, 8.2 ng/mL)",
            "BARACLUDE (entecavir) US Prescribing Information, DailyMed, "
            "section 12.3 (SS Cmax 4.2 ng/mL at 0.5 mg; 8.2 ng/mL at 1 mg)"),
    "LOP": (0.0008,
            "clinical Cmax (loperamide 2 mg oral, ~0.4 ng/mL by modern LC-MS; "
            "historical RIA upper bound ~0.003 uM)",
            "Kamali F & Huang ML. Br J Clin Pharmacol 1996;41(2):125-128 "
            "(PMID 8838438)"),
    "PAR": (0.15,
            "clinical steady-state Cmax (paroxetine 20 mg/day, ~50 ng/mL, "
            "interpolated from the label's 30 mg/day value 61.7 ng/mL)",
            "PAXIL (paroxetine HCl) US Prescribing Information, DailyMed, "
            "section 12.3; Sindrup SH, et al. Clin Pharmacol Ther "
            "1992;51(3):278-287 (PMID 1531950)"),
    "SB-431542": (10.0,
                  "standard in-vitro assay concentration (TGF-beta/ALK5 "
                  "inhibitor; no clinical Cmax exists)",
                  "standard SB-431542 in-vitro assay concentration (10 uM)"),
    "URS": (8.0,
            "clinical Cmax (ursodiol 300 mg oral, ~3.1 ug/mL, dose-scaled "
            "from the 500 mg measured value 5.2 ug/mL)",
            "Huang L, et al. Front Pharmacol 2023;14:1168144 (PMID 37138846)"),
}

PAPER_TABLE1_CITATION = (
    "DToxS-primary paper Table 1 'Concentration (uM)' (clinical max plasma "
    "concentrations): Sharma A, et al. Nat Commun 11:4809 (2020), "
    "DOI 10.1038/s41467-020-18396-7")
IPSC_SERIES_CITATION = (
    "GSE217423 iPSC-CM series matrix treatment concentrations "
    "(Sharma A, et al. Sci Transl Med 2022, GSE217423)")
TGF_CITATION = PRIMARY_ONLY_CONC["SB-431542"][2]

# ---------------------------------------------------------------------------
# Full map: code -> (drug name, concentration uM, source label)
# ---------------------------------------------------------------------------
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

# The LINCS well-level package codes SB-431542 as TGF; the map keys on the
# primary platform's compound code.
CODE_ALIASES = {"TGF": "SB-431542"}

PRIMARY_CONC = {}
for _code, _conc in PAPER_TABLE1_CONC.items():
    PRIMARY_CONC[_code] = (_PRIMARY_NAMES[_code], _conc, PAPER_TABLE1_CITATION)
for _code, _conc in IPSC_SERIES_CONC.items():
    PRIMARY_CONC[_code] = (_PRIMARY_NAMES[_code], _conc, IPSC_SERIES_CITATION)
for _code, (_conc, _src, _cit) in PRIMARY_ONLY_CONC.items():
    PRIMARY_CONC[_code] = (_PRIMARY_NAMES[_code], _conc, f"{_src}; {_cit}")

# Dataset-level summary helpers used by the step-03 reports.
PRIMARY_KI_CODES = set(PAPER_TABLE1_CONC)
PRIMARY_SHARED_NONKI_CODES = set(IPSC_SERIES_CONC)
PRIMARY_ONLY_CODES = set(PRIMARY_ONLY_CONC)
PRIMARY_CODES = set(PRIMARY_CONC)


def primary_dose(compound_code):
    """uM concentration for a primary compound code (or None if unmapped)."""
    code = CODE_ALIASES.get(compound_code, compound_code)
    row = PRIMARY_CONC.get(code)
    return None if row is None else row[1]


if __name__ == "__main__":
    print(f"primary concentration map: {len(PRIMARY_CONC)} compounds "
          f"({len(PRIMARY_KI_CODES)} KIs, {len(PRIMARY_SHARED_NONKI_CODES)} "
          f"shared non-KIs, {len(PRIMARY_ONLY_CODES)} primary-only)")
    print("  dose-bin distribution (author assign_pert_dose_corrected):")
    from collections import Counter

    def dose_idx_of(dose):
        if 0 <= dose < 0.21:
            return 0
        if dose < 0.41:
            return 1
        if dose < 0.71:
            return 2
        if dose < 1.01:
            return 3
        if dose < 1.51:
            return 4
        if dose < 3.1:
            return 5
        if dose < 4.1:
            return 6
        if dose < 7.1:
            return 7
        if dose < 12.1:
            return 8
        return 9

    bins = Counter(dose_idx_of(v[1]) for v in PRIMARY_CONC.values())
    for b in sorted(bins):
        print(f"    bin {b}: {bins[b]} drugs")
