"""Compound identity, pharmacology tiers and exclusions.

Identity rules (ADR 0026 decisions 4-5)
---------------------------------------
* ``drug_key`` is the desalted parent InChIKey **with stereochemistry retained**
  (full 27-character key).  Connectivity-only 14-character keys are rejected:
  doxorubicin and epirubicin are 4'-epimers and share the connectivity block.
* ``inchikey`` is the exact registered structure key (salt forms keep their own
  key); ``drug_key`` is the parent they collapse to.
* Structure wins over name.  The DToxS ``AMI`` row carries a
  pentyl-benzofuran homolog structure, so it is a different compound from
  amiodarone even though the audit's first name hit was amiodarone.
* Tiers are boolean columns, never pre-applied filters:
    ``atc-moa``       ATC code or ChEMBL mechanism
    ``atc-moa-mesh``  that, or a PubChem MeSH class
    ``all-smiles``    any resolvable small-molecule structure
* The only manual supplement is SB-431542 on the curated DToxS / gse244740
  panels (Inman et al. 2002, Mol Pharmacol 62:65-74).  gse262419 stays
  automatic by construction.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from .annotation import MANUAL_SOURCE_PREFIX, MANUAL_SUPPLEMENT, CompoundAnnotation
from .text import normalize_name

SMALL_MOLECULE = "small_molecule"
NON_SMALL_MOLECULES = {"antibody", "protein", "peptide"}

TIER_LABELS = {"atc-moa": "pharm_atc_moa",
               "atc-moa-mesh": "pharm_mesh",
               "all-smiles": "pharm_all_smiles"}


def dataset_members(value: object) -> set[str]:
    """Parse the ``datasets`` column (``;``-joined source ids)."""
    return {part for part in str(value).split(";") if part}

CURATED_PANELS = {"dtoxs_primary", "dtoxs_ipsc", "gse244740"}


class CompoundError(RuntimeError):
    """A compound violated the frozen identity contract."""


@dataclass
class Compound:
    drug_key: str
    inchikey: str
    compound_name: str
    smiles: str
    pubchem_cid: str = ""
    chembl_id: str = ""
    atc: str = ""
    moa: str = ""
    mesh_classes: str = ""
    annotation_source: str = ""
    structure_source: str = ""
    datasets: set[str] = field(default_factory=set)
    n_wells: int = 0
    notes: str = ""

    @property
    def pharm_atc_moa(self) -> bool:
        return bool(self.atc or self.moa)

    @property
    def pharm_mesh(self) -> bool:
        return bool(self.pharm_atc_moa or self.mesh_classes)

    @property
    def pharm_all_smiles(self) -> bool:
        return bool(self.drug_key)


@dataclass
class CompoundTable:
    compounds: pd.DataFrame          # one row per (drug_key, compound_name)
    exclusions: pd.DataFrame

    def members(self, tier: str) -> set[str]:
        column = TIER_LABELS[tier]
        return set(self.compounds.loc[self.compounds[column], "drug_key"])


def canonical_smiles(smiles: str) -> str:
    if not smiles:
        return ""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def _validate_key(key: str, name: str) -> str:
    key = str(key or "").strip()
    if not key:
        raise CompoundError(f"{name}: no resolvable structure")
    if len(key) != 27 or key.count("-") != 2:
        raise CompoundError(
            f"{name}: {key!r} is not a full InChIKey (connectivity-only keys "
            "are forbidden because they merge 4'-epimers)")
    return key


def _small_molecule_records(source_meta: pd.DataFrame) -> pd.DataFrame:
    mask = ((~source_meta["is_control"]) & (~source_meta["is_blank"])
            & (source_meta["entity_type"] == SMALL_MOLECULE))
    return source_meta.loc[mask].copy()


def _group_records(frame: pd.DataFrame) -> pd.DataFrame:
    return (frame.groupby(
        ["dataset_id", "compound_name", "compound_code", "dataset_inchikey",
         "dataset_smiles", "entity_type"], dropna=False, sort=True)
        .size().rename("n_wells").reset_index())


def build_compound_table(source_metas: dict[str, pd.DataFrame],
                         annotation: CompoundAnnotation) -> CompoundTable:
    frames = []
    exclusions = []
    for dataset_id, meta in source_metas.items():
        small = _small_molecule_records(meta)
        frames.append(small)
        excluded = meta.loc[
            (~meta["is_control"]) & (~meta["is_blank"])
            & (meta["entity_type"].isin(NON_SMALL_MOLECULES))]
        for name, group in excluded.groupby("compound_name", sort=True):
            exclusions.append({
                "dataset_id": dataset_id, "compound_name": name,
                "compound_code": ",".join(sorted(set(group["compound_code"]))),
                "entity_type": group["entity_type"].iloc[0],
                "n_wells": len(group),
                "reason": f"non-small-molecule perturbagen ({group['entity_type'].iloc[0]})",
            })
        blanks = meta.loc[meta["is_blank"]]
        if len(blanks):
            exclusions.append({
                "dataset_id": dataset_id, "compound_name": "Blank",
                "compound_code": "",
                "entity_type": "technical_blank", "n_wells": len(blanks),
                "reason": "technical blank (concentration 0; neither control nor perturbation)",
            })
    small = pd.concat(frames, ignore_index=True)
    grouped = _group_records(small)

    resolved: dict[tuple[str, str], Compound] = {}
    unresolved = []
    for record in grouped.to_dict("records"):
        name = normalize_name(record["compound_name"])
        dataset_key = str(record["dataset_inchikey"] or "")
        dataset_smiles = str(record["dataset_smiles"] or "")
        ann = annotation.get(name, dataset_key)
        if dataset_key:
            # Structure-first: the source's own structure is the identity.
            drug_key = _validate_key(dataset_key, name)
            exact = (ann or {}).get("inchikey", "") or drug_key
            smiles = canonical_smiles(dataset_smiles) or (ann or {}).get("smiles", "")
            structure_source = "DToxS_metadata"
        else:
            if not ann or not ann.get("parent_inchikey"):
                unresolved.append({
                    "dataset_id": record["dataset_id"],
                    "compound_name": name,
                    "compound_code": record["compound_code"],
                    "entity_type": record["entity_type"],
                    "n_wells": record["n_wells"],
                    "reason": "no resolvable structure in the frozen external annotation",
                })
                continue
            drug_key = _validate_key(ann["parent_inchikey"], name)
            exact = ann.get("inchikey", "") or drug_key
            smiles = canonical_smiles(ann.get("smiles", "")) or ann.get("smiles", "")
            structure_source = "pubchem"

        key = (drug_key, name)
        if key not in resolved:
            atc = (ann or {}).get("atc", "") or ""
            moa = (ann or {}).get("moa", "") or ""
            mesh = (ann or {}).get("mesh_classes", "") or ""
            annotation_source = ("chembl" if (atc or moa)
                                 else ("pubchem" if mesh else ""))
            resolved[key] = Compound(
                drug_key=drug_key, inchikey=exact, compound_name=name,
                smiles=smiles,
                pubchem_cid=str((ann or {}).get("pubchem_cid", "") or ""),
                chembl_id=(ann or {}).get("chembl_id", "") or "",
                atc=atc, moa=moa, mesh_classes=mesh,
                annotation_source=annotation_source,
                structure_source=structure_source,
            )
        compound = resolved[key]
        compound.datasets.add(record["dataset_id"])
        compound.n_wells += int(record["n_wells"])

    # Manual supplement: curated panels only.
    for compound in resolved.values():
        manual = MANUAL_SUPPLEMENT.get(compound.drug_key)
        if manual is None:
            continue
        if not compound.datasets.issubset(CURATED_PANELS):
            raise CompoundError(
                f"{compound.compound_name}: manual supplement is not allowed "
                f"on datasets {sorted(compound.datasets - CURATED_PANELS)}")
        if not (compound.atc or compound.moa):
            compound.moa = manual["moa"]
            compound.annotation_source = MANUAL_SOURCE_PREFIX + manual["citation"]
            compound.notes = "manual pharmacology supplement"

    if not resolved:
        raise CompoundError("no compound resolved")

    rows = []
    for compound in sorted(resolved.values(),
                           key=lambda c: (c.drug_key, c.compound_name)):
        rows.append({
            "drug_key": compound.drug_key,
            "inchikey": compound.inchikey,
            "compound_name": compound.compound_name,
            "smiles": compound.smiles,
            "pubchem_cid": compound.pubchem_cid,
            "chembl_id": compound.chembl_id,
            "atc": compound.atc,
            "moa": compound.moa,
            "mesh_classes": compound.mesh_classes,
            "pharm_atc_moa": compound.pharm_atc_moa,
            "pharm_mesh": compound.pharm_mesh,
            "pharm_all_smiles": compound.pharm_all_smiles,
            "annotation_source": compound.annotation_source,
            "structure_source": compound.structure_source,
            "datasets": ";".join(sorted(compound.datasets)),
            "n_wells": compound.n_wells,
            "notes": compound.notes,
        })
    compounds = pd.DataFrame(rows)

    exclusion_rows = exclusions + unresolved
    exclusions_frame = pd.DataFrame(exclusion_rows)
    if exclusions_frame.empty:
        exclusions_frame = pd.DataFrame(
            columns=["dataset_id", "compound_name", "compound_code",
                     "entity_type", "n_wells", "reason"])
    else:
        exclusions_frame = exclusions_frame.sort_values(
            ["dataset_id", "compound_name"], kind="stable").reset_index(drop=True)

    _identity_regressions(compounds)
    return CompoundTable(compounds, exclusions_frame)


def _identity_regressions(compounds: pd.DataFrame) -> None:
    """Hard identity checks that must never regress."""
    by_name = {}
    for row in compounds.to_dict("records"):
        by_name.setdefault(row["compound_name"], set()).add(row["drug_key"])
    doxo = by_name.get("Doxorubicin", set())
    epi = by_name.get("Epirubicin", set())
    if doxo and epi and doxo == epi:
        raise CompoundError("doxorubicin and epirubicin collapsed to one drug_key "
                            "(connectivity-only key regression)")
    # Salt/free-base pairs must share a key while remaining separate names.
    pairs = [("Tamoxifen", "Tamoxifen citrate"),
             ("Diltiazem", "Diltiazem hydrochloride"),
             ("Estradiol", "17beta-Estradiol"),
             ("Perfluorooctanoic acid", "Ammonium perfluorooctanoate")]
    for free, salt in pairs:
        if free in by_name and salt in by_name:
            if by_name[free] != by_name[salt]:
                raise CompoundError(
                    f"{free} / {salt} do not share a drug_key: "
                    f"{by_name[free]} vs {by_name[salt]}")


def annotate_source_meta(meta: pd.DataFrame,
                         table: CompoundTable) -> pd.DataFrame:
    """Attach drug identity + tier membership to every source row."""
    out = meta.copy()
    by_name: dict[str, list[dict]] = {}
    for row in table.compounds.to_dict("records"):
        by_name.setdefault(normalize_name(row["compound_name"]), []).append(row)
    for column in ("drug_key", "inchikey", "smiles", "atc", "moa",
                   "mesh_classes", "pharm_atc_moa", "pharm_mesh",
                   "pharm_all_smiles", "annotation_source"):
        out[column] = ""
    for index, row in out.iterrows():
        candidates = by_name.get(normalize_name(row["compound_name"]), [])
        compound = None
        if len(candidates) == 1:
            compound = candidates[0]
        elif candidates:
            # Disambiguate a shared display name by the source structure.  The
            # source key is either the parent (`drug_key`) or the exact salt
            # (`inchikey`).
            source_key = str(row.get("dataset_inchikey") or "")
            for candidate in candidates:
                if source_key and source_key in (candidate["drug_key"],
                                                 candidate["inchikey"]):
                    compound = candidate
                    break
            if compound is None and len({c["drug_key"] for c in candidates}) == 1:
                compound = candidates[0]
        if compound is None:
            out.at[index, "drug_key"] = ""
            out.at[index, "pharm_atc_moa"] = False
            out.at[index, "pharm_mesh"] = False
            out.at[index, "pharm_all_smiles"] = False
            continue
        out.at[index, "drug_key"] = compound["drug_key"]
        out.at[index, "inchikey"] = compound["inchikey"]
        out.at[index, "smiles"] = compound["smiles"]
        out.at[index, "atc"] = compound["atc"]
        out.at[index, "moa"] = compound["moa"]
        out.at[index, "mesh_classes"] = compound["mesh_classes"]
        out.at[index, "pharm_atc_moa"] = compound["pharm_atc_moa"]
        out.at[index, "pharm_mesh"] = compound["pharm_mesh"]
        out.at[index, "pharm_all_smiles"] = compound["pharm_all_smiles"]
        out.at[index, "annotation_source"] = compound["annotation_source"]
    return out


def tier_counts(table: CompoundTable,
                datasets: list[str] | None = None) -> dict:
    compounds = table.compounds
    if datasets:
        wanted = set(datasets)
        compounds = compounds[compounds["datasets"].apply(
            lambda value: bool(dataset_members(value) & wanted))]
    counts = {}
    for tier, column in TIER_LABELS.items():
        counts[tier] = compounds.loc[compounds[column], "drug_key"].nunique()
    return counts
