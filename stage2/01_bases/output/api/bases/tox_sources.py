"""Measurement adapters for the independently admitted TOX sources (ADR 0072).

``prepare_source`` exposes source membership, scale and matched controls before
any fitting. ``source_contrasts`` applies the shared source-appropriate DEG
worker; callers collapse its strata once, giving one vote per dataset.
"""

from pathlib import Path
import re
import gzip
import tarfile
import io
import numpy as np
import pandas as pd
from . import bulk_deg, disease_deg
from .release import ReleaseError

HGNC = "pert_data/Annotation/source/hgnc_complete_set.txt"


def _matrix(frame, symbols=None):
    result = frame.copy()
    if symbols is None:
        symbols = frame.attrs.get("gene_symbol_rows")
    result.insert(
        0,
        "gene_symbol",
        list(
            [
                frame.attrs.get("gene_symbols", {}).get(str(g), str(g))
                for g in frame.index
            ]
            if symbols is None
            else symbols
        ),
    )
    result.insert(0, "source_gene_id", frame.index.astype(str))
    return result.reset_index(drop=True)


def _part(frame, samples, kind="counts", formula="~0 + group"):
    return {
        "matrix": _matrix(frame),
        "samples": samples,
        "design": {
            "kind": kind,
            "formula": formula,
            "factor": "group",
            "contrasts": [
                {"label": "TOX-vehicle", "numerator": "TOX", "denominator": "vehicle"}
            ],
        },
    }


def prepare_source(repo, source, admitted_drugs):
    """Prepare deposited measurements and exact admitted exposure/control audit."""
    repo = Path(repo)
    if (
        "disposition" in admitted_drugs
        and not admitted_drugs.disposition.eq("admitted").all()
    ):
        raise ReleaseError("source adapter received non-admitted drug records")
    if admitted_drugs.empty:
        raise ReleaseError(f"{source}: no independently admitted drugs")
    if source in {
        "gse106297",
        "gse157282",
        "gse242692",
        "gse255505",
        "gse313746",
    } and "DOX" not in set(admitted_drugs.source_id.astype(str)):
        raise ReleaseError(f"{source}: no admitted DOX identity")
    if source != "gse313746":
        return _attach_metadata(
            repo, source, _prepare_other(repo, source, admitted_drugs)
        )
    path = "disease_data/GSE313746/source/GSE313746_Counts_RNA_EMPfortmiller.txt.gz"
    counts = pd.read_csv(repo / path, sep=r"\s+", index_col=0)
    counts.index = counts.index.astype(str)
    hgnc = pd.read_csv(repo / HGNC, sep="\t", dtype=str, keep_default_na=False)
    approved = hgnc[hgnc.status.eq("Approved") & hgnc.entrez_id.ne("")].copy()
    approved["locus_rank"] = approved.locus_group.map(
        {"protein-coding gene": 0, "non-coding RNA": 1, "pseudogene": 2}
    ).fillna(99)
    approved["hgnc_num"] = pd.to_numeric(
        approved.hgnc_id.str.replace("HGNC:", "", regex=False)
    )
    approved = approved.sort_values(
        ["entrez_id", "locus_rank", "hgnc_num"], kind="stable"
    )
    mapping = (
        approved.drop_duplicates("entrez_id").set_index("entrez_id").symbol.to_dict()
    )
    unmapped = int(sum(g not in mapping for g in counts.index))
    counts = counts.loc[[g for g in counts.index if g in mapping]].copy()
    counts.index = [mapping[g] for g in counts.index]
    counts = counts.groupby(level=0, sort=True).sum()
    metadata = []
    for sample in counts:
        match = re.fullmatch(r"(DOX|VEH)_tx(\d+)_(\d+)(R?)", sample)
        if not match:
            raise ReleaseError(f"cannot parse GSE313746 sample {sample}")
        drug, time, donor, repeat = match.groups()
        metadata.append(
            {
                "sample": sample,
                "source_id": drug,
                "donor": donor,
                "stratum": f"tx{time}",
                "time": time,
                "role": "treated" if drug == "DOX" else "vehicle",
                "analysis_sample": f"{drug}_tx{time}_{donor}",
                "technical_repeat": bool(repeat),
            }
        )
    audit = []
    strata = []
    for stratum in sorted({m["stratum"] for m in metadata}):
        rows = [m for m in metadata if m["stratum"] == stratum]
        grouped = {}
        for row in rows:
            grouped.setdefault(row["analysis_sample"], []).append(row["sample"])
        matrix = pd.DataFrame(
            {name: counts[wells].sum(axis=1) for name, wells in grouped.items()}
        )
        matrix.attrs = counts.attrs.copy()
        samples = pd.DataFrame(
            [
                {
                    "sample": name,
                    "group": "TOX" if name.startswith("DOX") else "vehicle",
                    "donor": "donor_"
                    + next(m["donor"] for m in rows if m["analysis_sample"] == name),
                }
                for name in matrix
            ]
        )
        for row in rows:
            controls = [
                m["sample"]
                for m in rows
                if m["role"] == "vehicle" and m["donor"] == row["donor"]
            ]
            if not controls:
                raise ReleaseError("GSE313746 donor lacks matched vehicle")
            audit.append(dict(row, matched_controls=controls))
        part = _part(matrix, samples, formula="~0 + group + donor")
        part["stratum"] = stratum
        strata.append(part)
    return _attach_metadata(
        repo,
        source,
        {
            "strata": strata,
            "audit": audit,
            "inputs": [path, HGNC],
            "method": {
                "name": "limma-voom",
                "scale": "raw integer counts; technical repeat counts summed per donor/arm/recovery",
                "design": "within recovery, donor fixed effect; admitted DOX pooled vs matched vehicle",
                "entrez_mapping": "approved HGNC; protein-coding < non-coding RNA < pseudogene < other, then lowest HGNC integer; counts summed per symbol",
                "unmapped_entrez_rows": unmapped,
            },
        },
    )


def _assemble(
    frame,
    rows,
    inputs,
    *,
    kind="counts",
    formula="~0 + group",
    scale="raw integer counts",
    conversion=None,
):
    """Bind a source measurement column to its declared stratum and control."""
    audit = []
    strata = []
    for row in rows:
        if row["sample"] not in frame:
            raise ReleaseError(f"sample absent from deposited matrix: {row['sample']}")
    for stratum in sorted({r["stratum"] for r in rows if r["role"] == "treated"}):
        selected = [
            dict(r)
            for r in rows
            if r["stratum"] == stratum and r["role"] in ("treated", "vehicle")
        ]
        controls = [r for r in selected if r["role"] == "vehicle"]
        if not controls:
            raise ReleaseError(f"{stratum}: no matched vehicle")
        for row in selected:
            matched = [
                r["sample"]
                for r in controls
                if not row.get("donor") or r.get("donor") == row["donor"]
            ]
            if not matched:
                raise ReleaseError(f"{row['sample']}: no donor-matched vehicle")
            row["matched_controls"] = matched
            row["analysis_sample"] = row["sample"]
            audit.append(row)
        samples = pd.DataFrame(
            [
                dict(
                    sample=r["sample"],
                    group="TOX" if r["role"] == "treated" else "vehicle",
                    **{k: r[k] for k in ("donor", "batch", "family") if k in r},
                )
                for r in selected
            ]
        )
        part = _part(
            frame[samples["sample"].tolist()], samples, kind=kind, formula=formula
        )
        part["stratum"] = stratum
        strata.append(part)
    if not strata:
        raise ReleaseError("source has no admitted exposure with matched vehicle")
    audit.extend(
        dict(
            r, matched_controls=[], reason=r.get("reason", "drug identity not admitted")
        )
        for r in rows
        if r["role"] == "excluded"
    )
    result = {
        "strata": strata,
        "audit": audit,
        "inputs": inputs,
        "method": {
            "name": "limma-voom" if kind == "counts" else "limma lmFit/eBayes",
            "scale": scale,
            "design": formula,
            "collapse": "maximum-effect across declared strata; one dataset vote",
        },
    }
    if conversion is not None:
        result["method"]["integerisation"] = conversion
    return result


def _archive(repo, path, source):
    """Read only the compact gene-by-library source tables, never model payloads."""
    vectors = {}
    symbols = {}
    with tarfile.open(repo / path) as archive:
        for member in archive:
            if not member.isfile():
                continue
            name = (
                Path(member.name)
                .name.removesuffix(".gz")
                .removesuffix(".txt")
                .removesuffix(".fpkm_tracking")
            )
            raw = archive.extractfile(member).read()
            if member.name.endswith(".gz"):
                raw = gzip.decompress(raw)
            if source == "gse157282":
                table = pd.read_csv(io.BytesIO(raw), sep="\t", header=None)
                numeric = pd.to_numeric(table.iloc[:, 3], errors="coerce")
                table = table[numeric.notna()]
                vector = pd.Series(
                    pd.to_numeric(table.iloc[:, 3]).to_numpy(float),
                    index=table.iloc[:, 0].astype(str),
                )
                symbols.update({g: g for g in table.iloc[:, 0].astype(str)})
            elif source == "gse106297":
                table = pd.read_csv(io.BytesIO(raw), sep="\t")
                vector = pd.Series(
                    table.FPKM.to_numpy(float), index=table.tracking_id.astype(str)
                )
                symbols.update(
                    {
                        str(g): str(v) if str(v) not in ("-", "", "nan") else str(g)
                        for g, v in zip(table.tracking_id, table.gene_short_name)
                    }
                )
            else:
                table = pd.read_csv(io.BytesIO(raw), sep="\t", index_col=0)
                vector = table.iloc[:, 0]
            if vector.index.duplicated().any():
                # Deposits can repeat the same symbol over gene IDs. Keep IDs for
                # alignment, rather than inventing independent observations.
                vector = vector.groupby(level=0, sort=False).sum()
            vectors[name] = vector
    result = pd.DataFrame(vectors).fillna(0)
    result.attrs["gene_symbols"] = symbols
    return result


def _prepare_other(repo, source, drugs):
    accession = source.upper()
    category = "pert_data" if source in ("gse247541", "gse262419") else "disease_data"
    base = f"{category}/{accession}/source"
    root = repo / base
    admitted = set(drugs.source_id.astype(str))
    rows = []
    if source in ("gse157282", "gse106297", "gse255505"):
        path = f"{base}/{accession}_RAW.tar"
        frame = _archive(repo, path, source)
        inputs = [path]
        series_path = f"{base}/{accession}_series_matrix.txt.gz"
        series = (
            bulk_deg.read_series_matrix(repo / series_path)
            if (repo / series_path).exists()
            else None
        )
        for sample in frame:
            role = "excluded"
            drug = ""
            family = ""
            if source == "gse157282":
                role = (
                    "vehicle"
                    if re.search(r"_Control[0-9]+_", sample)
                    else "treated"
                    if re.search(r"_Dox[0-9]+_", sample)
                    else "excluded"
                )
            elif source == "gse255505":
                role = (
                    "vehicle"
                    if "_wo_DOX_" in sample
                    else "treated"
                    if "_wi_DOX_" in sample
                    else "excluded"
                )
            else:
                gsm = sample.split("_")[0]
                title = str(series.loc[gsm, "title"])
                dose = float(re.search(r"DOXO ([0-9.]+)", title).group(1))
                role = "vehicle" if dose == 0 else "treated"
                family = re.search(r"sample ([AB])", title).group(1)
            row = dict(
                sample=sample,
                role=role,
                source_id="DOX" if role == "treated" else "vehicle",
                stratum="all",
                reason=(
                    "co-exposure or non-admitted intervention"
                    if role == "excluded"
                    else ""
                ),
            )
            if family:
                row["family"] = family
            rows.append(row)
        if series is not None:
            inputs.append(series_path)
        if source == "gse255505":
            frame, stats = bulk_deg.integerise_counts(frame)
            return _assemble(
                frame,
                rows,
                inputs,
                conversion=stats,
                scale="fractional estimated counts; shared nearest-even integerisation",
            )
        expression = np.log2(frame + 1)
        expression.attrs = frame.attrs.copy()
        return _assemble(
            expression,
            rows,
            inputs,
            kind="expression",
            formula="~0 + group + family" if source == "gse106297" else "~0 + group",
            scale="log2(FPKM + 1)",
        )
    if source == "gse232331":
        path = f"{base}/{accession}_hg38-counts.tsv.gz"
        frame = pd.read_csv(repo / path, sep="\t", index_col=0)
        for sample in frame:
            treatment = (
                "DMSO"
                if "_DMSO_" in sample
                else (
                    "DOX"
                    if "DOXO_" in sample
                    else (
                        "ACLA"
                        if "ACLA_" in sample
                        else "AMR"
                        if "AMR_" in sample
                        else "other"
                    )
                )
            )
            rows.append(
                dict(
                    sample=sample,
                    source_id=treatment,
                    role=(
                        "vehicle"
                        if treatment == "DMSO"
                        else "treated"
                        if treatment in admitted
                        else "excluded"
                    ),
                    stratum="all",
                    batch="date_" + re.search(r"_([0-9]{6})_", sample).group(1),
                    dose_um=0 if treatment == "DMSO" else 5,
                    time_h=24,
                )
            )
        return _assemble(frame, rows, [path], formula="~0 + group + batch")
    if source == "gse235470":
        path = f"{base}/{accession}_rawcounts.counts.txt.gz"
        table = pd.read_csv(repo / path, sep="\t", index_col=0)
        symbols = table.pop("symbol")
        frame = table
        # Preserve row-wise symbols even when the deposit repeats a feature ID.
        frame.attrs["gene_symbol_rows"] = symbols.astype(str).tolist()
        for sample in frame:
            token = re.sub(r"[0-9]+$", "", sample)
            role = (
                "vehicle"
                if token == "DMSO"
                else "treated"
                if token in admitted
                else "excluded"
            )
            rows.append(
                dict(
                    sample=sample,
                    source_id=token,
                    stratum="all",
                    role=role,
                    dose_um=5 if role == "treated" else 0,
                    time_h=48,
                )
            )
        return _assemble(frame, rows, [path])
    if source == "gse242692":
        path = f"{base}/{accession}_countsMatrix_raw.txt.gz"
        frame = pd.read_csv(repo / path, sep="\t", index_col=0)
        series_path = f"{base}/{accession}_series_matrix.txt.gz"
        series = bulk_deg.read_series_matrix(repo / series_path)
        for sample in frame:
            donor, dose = re.split(r"_DOX_", sample, flags=re.IGNORECASE)
            matched = series[series.description.eq(donor)]
            concentration = (
                0
                if dose == "DMSO"
                else 0.3
                if dose == "03uM"
                else float(dose.replace("uM", ""))
            )
            title_doses = matched.title.str.extract(r" ([0-9.]+) ")[0].astype(float)
            matches = matched.loc[title_doses.eq(concentration)]
            flows = set(matches.flowcell)
            if len(flows) != 1:
                raise ReleaseError(f"{sample}: ambiguous flowcell metadata")
            role = "vehicle" if dose == "DMSO" else "treated"
            rows.append(
                dict(
                    sample=sample,
                    source_id="DOX" if role == "treated" else "DMSO",
                    role=role,
                    donor=donor,
                    stratum="all",
                    flowcell=next(iter(flows)),
                    source_title=str(matches.title.iloc[0]),
                    dose_um=concentration,
                )
            )
        result = _assemble(
            frame, rows, [path, series_path], formula="~0 + group + donor"
        )
        part = result["strata"][0]
        samples = part["samples"]
        flow_by_sample = {r["sample"]: r["flowcell"] for r in rows}
        base_design = pd.concat(
            [
                pd.get_dummies(samples.group, dtype=float),
                pd.get_dummies(samples.donor, dtype=float).iloc[:, 1:],
            ],
            axis=1,
        )
        design = base_design.to_numpy()
        rank = int(np.linalg.matrix_rank(design))
        retained = []
        omitted = []
        for flow in sorted(set(flow_by_sample.values())):
            column = np.array(
                [int(flow_by_sample[name] == flow) for name in samples["sample"]]
            )
            candidate = np.column_stack([design, column])
            new_rank = int(np.linalg.matrix_rank(candidate))
            if new_rank > rank:
                field = "flow_" + re.sub("[^a-zA-Z0-9]", "_", flow)
                samples[field] = column
                retained.append(field)
                design = candidate
                rank = new_rank
            else:
                omitted.append(flow)
        formula = "~0 + group + donor" + "".join(" + " + field for field in retained)
        part["design"]["formula"] = formula
        result["method"].update(
            design=formula,
            design_rank=rank,
            design_columns=int(design.shape[1]),
            retained_flow_covariates=retained,
            omitted_collinear_flow_effects=omitted,
            flow_policy="sequencing batch adjustment only; matched vehicle requires donor identity, not same sequencing run",
        )
        return result
    if source == "gse247541":
        path = f"{base}/{accession}_raw_counts.csv.gz"
        frame = pd.read_csv(repo / path, index_col=0)
        for sample in frame:
            drug = sample.split("_rep")[0]
            role = (
                "vehicle"
                if drug == "DMSO_control"
                else "treated"
                if drug in admitted
                else "excluded"
            )
            rows.append(dict(sample=sample, source_id=drug, stratum="all", role=role))
        return _assemble(frame, rows, [path])
    if source == "gse262419":
        return _prepare_262419(repo, base, drugs)
    raise ReleaseError(f"unsupported TOX source: {source}")


def _prepare_262419(repo, base, drugs):
    path = f"{base}/GSE262419_hash.csv.gz"
    panel = pd.read_csv(repo / path)
    admitted = set(drugs.source_id.astype(str))
    parts = []
    audit = []
    inputs = [path]
    for plate in sorted(set(panel.loc[panel.Chemical_name.isin(admitted), "Plate_ID"])):
        plate = int(plate)
        metadata = panel[panel.Plate_ID.eq(plate)].copy()
        metadata["sample"] = [
            f"Plate{plate}-{re.match(r'[A-Za-z]+', w).group()}{int(re.search(r'[0-9]+', w).group()):02d}"
            for w in metadata.Well_ID.astype(str)
        ]
        file = f"{base}/GSE262419_Plate{plate:02d}.csv.gz"
        inputs.append(file)
        available = set(pd.read_csv(repo / file, nrows=0).columns) - {"Genes"}
        rows = []
        for row in metadata.itertuples(index=False):
            role = (
                "vehicle"
                if row.Chemical_name == "VEH"
                else "treated"
                if row.Chemical_name in admitted
                else "excluded"
            )
            entry = dict(
                sample=row.sample,
                role=role,
                source_id=row.Chemical_name,
                stratum=f"plate{plate}",
                plate=plate,
                dose_um=float(row.Chemical_Concentration_uM),
            )
            if row.sample not in available:
                audit.append(
                    dict(
                        entry,
                        role="excluded",
                        reason="source well has no deposited count column",
                        matched_controls=[],
                    )
                )
            else:
                rows.append(entry)
        selected = [r["sample"] for r in rows if r["role"] != "excluded"]
        if not any(r["role"] == "treated" for r in rows):
            continue
        table = pd.read_csv(repo / file, usecols=["Genes"] + selected, index_col=0)
        table.attrs["gene_symbols"] = {
            str(g): re.sub(r"_[0-9]+$", "", str(g)) for g in table.index
        }
        # Blank/media and unadmitted drugs are audited but never read into the fit.
        included = [r for r in rows if r["role"] != "excluded"]
        value = _assemble(table, included, [file])
        parts.extend(value["strata"])
        audit.extend(value["audit"])
        audit.extend(
            dict(
                r,
                reason="technical blank/media or drug identity not admitted",
                matched_controls=[],
            )
            for r in rows
            if r["role"] == "excluded"
        )
    if not parts:
        raise ReleaseError("GSE262419 has no measured independently admitted exposure")
    return {
        "strata": parts,
        "audit": audit,
        "inputs": inputs,
        "method": {
            "name": "limma-voom",
            "scale": "raw integer counts",
            "design": "same-plate admitted exposures pooled vs VEH; plate-specific contrasts",
            "collapse": "maximum effect across plate contrasts; one dataset vote",
        },
    }


def source_contrasts(repo, source_id, admitted_drugs, workdir):
    """Fit prepared strata, retaining the complete sample/control audit."""
    prepared = prepare_source(repo, source_id, admitted_drugs)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    frames = []
    contrasts = []
    summaries = []
    for index, part in enumerate(prepared["strata"]):
        matrix = workdir / f"{source_id}_{index}.matrix.tsv"
        samples = workdir / f"{source_id}_{index}.samples.csv"
        out = workdir / f"{source_id}_{index}.contrasts.csv"
        part["matrix"].to_csv(matrix, sep="\t", index=False)
        part["samples"].to_csv(samples, index=False)
        summary = disease_deg.run_deg_worker(matrix, samples, part["design"], out)
        frame = pd.read_csv(out)
        label = f"TOX|{part['stratum']}"
        frame["contrast"] = label
        n_case = int(part["samples"].group.eq("TOX").sum())
        n_control = int(part["samples"].group.eq("vehicle").sum())
        frame["n_case"] = n_case
        frame["n_control"] = n_control
        frames.append(frame)
        contrasts.append(
            {
                "contrast": label,
                "stratum": part["stratum"],
                "n_case": n_case,
                "n_control": n_control,
            }
        )
        summaries.append(summary)
    prepared.pop("strata")
    prepared["method"]["fits"] = summaries
    return pd.concat(frames, ignore_index=True), contrasts, prepared


def _attach_metadata(repo, source, result):
    category = "pert_data" if source in ("gse247541", "gse262419") else "disease_data"
    path = f"{category}/{source.upper()}/source/{source.upper()}_series_matrix.txt.gz"
    if not (repo / path).is_file():
        return result
    series = bulk_deg.read_series_matrix(repo / path)
    if path not in result["inputs"]:
        result["inputs"].append(path)
    descriptions = {
        re.sub(r"^Library name: ", "", str(row.get("description", ""))): row
        for _, row in series.iterrows()
    }
    titles = {str(row.title): row for _, row in series.iterrows()}
    for row in result["audit"]:
        sample = row["sample"]
        gsm = sample.split("_")[0]
        record = None
        if gsm in series.index:
            record = series.loc[gsm]
        elif sample in descriptions:
            record = descriptions[sample]
        elif source == "gse242692":
            record = titles.get(row.get("source_title"))
        elif source == "gse262419":
            record = titles.get(sample.replace("-", "_"))
        elif source == "gse235470":
            match = re.fullmatch(r"(DMSO|DOX|SPEDOX)([0-9]+)", sample)
            if match:
                token, replicate = match.groups()
                record = next(
                    (
                        v
                        for k, v in titles.items()
                        if f" {token} " in k and k.endswith("replicate " + replicate)
                    ),
                    None,
                )
        if record is not None:
            row["source_metadata"] = {
                str(k): str(v) for k, v in record.to_dict().items()
            }
    return result
