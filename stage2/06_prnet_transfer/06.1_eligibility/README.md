# eligibility — 460-structure pharmacological eligibility audit

**Status: v1.3 built and verified (2026-11-02).** Ticket #1 of the PRnet
formal-transfer effort, revised after the ticket #1 review. The audit covers every
structure identity of the frozen temp_env panel and replaces the old automatic
`atc-moa` (113 keys) / `atc-moa-mesh` (252 keys) tiers as the basis for the main
training set. Nothing here trains a model; every number is an evidence count, not
a performance result.

Rule v1.3 fixes two rule-consistency defects the review found in the locked v1.2
rule: the mechanism-tool clause admitted agrochemical/industrial compounds on a
mechanism text alone (Ethoxyquin, Dibutyl phthalate, 2-Phenylphenol, Daminozide,
Trifloxystrobin, Cyfluthrin, Paraoxon), and the approved/clinical-versus-other-use
precedence was not stated and was applied inconsistently to natural products with
registered human evidence (Quercetin was `other_use` while Genistein was
`clinical`). v1.3 states the precedence explicitly, encodes it and re-decides
every one of the 460 keys with per-row evidence.

**The v1.2 record is preserved byte-identical** in
[`output/v1.2/`](output/v1.2) (ADR 0030: no silent overwrite). v1.3 is a separate
versioned artifact set in [`output/v1.3/`](output/v1.3), and
`output/current.json` is the machine-readable pointer to the current record.

| artifact set | rule | include | exclude | uncertain |
|---|---|---:|---:|---:|
| `output/v1.2/` (frozen) | `eligibility-rule-v1.2` sha256 `e1af8cdb…` | 126 | 193 | 141 |
| `output/v1.3/` (current) | `eligibility-rule-v1.3` sha256 `504b8e42…` | **118** | **314** | **28** |

## What this delivers

- **460/460 rows.** Every `drug_key` with `temp_env_all-smiles=True` in the frozen
  [`compounds.csv`](../../02_cm_data/output/datasets-1/compounds.csv) (467 name
  records collapse to 460 keys; no structure is merged or dropped). Each row
  carries all names, InChIKey, SMILES, decision, decision class, **decision
  precedence label**, the per-row human-use or other-use evidence string,
  documented phase, MoA, targets, ATC codes, MeSH, evidence URLs with retrieval
  timestamp, rationale, the three old-tier flags and the **v1.2 decision** it
  replaces.
- **A rule locked before the frozen input is built** (`execution/rules.py`;
  version and sha256 in every artifact). v1.0–v1.2 are earlier revisions of the
  same rule; no model result existed at any revision.
- **A change report** [`output/v1.3/change_report.csv`](output/v1.3/change_report.csv):
  one row for each of the 127 decisions that moved from v1.2, with the evidence
  that decided it.
- **Coverage counts** recomputed from the data: decision classes, precedence
  labels, human-use-marker coverage, ATC normalisation in both readings, MoA
  annotation shape, identity resolution, old-tier false negatives/positives.
- **An include list, an unresolved list, and handoff v2** — the ticket #3
  interface.

## The locked decision rule (v1.3)

`RULE_TEXT` in [`execution/rules.py`](execution/rules.py) is the canonical rule;
`rule_sha256()` hashes it and the hash is recorded in the table, summary, handoff,
change report, manifest and `current.json`. **`eligibility-rule-v1.3`**, sha256
`504b8e4215b6d3c4403124146f2ef63b16fe45dd08f9e829d75215f8ec459ce6`.

**Precedence, first match wins:**

> documented human drug use > documented environmental/industrial use > mechanism text

| step | clause | decision |
|---|---|---|
| **D0** | identity unresolved (frozen InChIKey maps to no PubChem/Broad/BindingDB record) | uncertain |
| **D1** | documented **human drug use** (approval or clinical phase, corroborated when an environmental/industrial use is also documented) | include — `approved` / `clinical` |
| **D2** | documented **environmental/industrial use** and no D1 match | exclude — `other_use` |
| **D3** | mechanism/target evidence **and** research-tool or human-use evidence | include — `mechanism_tool` |
| **D4** | nothing above documented | uncertain |

**Documented human drug use** is `H1` approval (PubChem Max Phase `Approved` or
Broad `clinical_phase` `Launched`), `H2` a documented clinical phase (Early
Phase 1, Phase 1/2/3/4, Phase 1/Phase 2, Phase 2/Phase 3, Withdrawn), or `H3` a
human therapeutic indication (Broad `indication`). It is **corroborated by a
human-use marker**: a 7-character human ATC code (not the Q veterinary group), a
human therapeutic indication, a PubChem `FDA Approved Drugs` listing, or a
PubChem `WHO Essential Medicines` listing.

**Why corroboration matters.** A bare status string is *not* human drug use when
the same record documents an environmental/industrial use: PubChem Max Phase is
ChEMBL-derived and marks agrochemicals `Approved`/`Phase 2` (Carbaryl, Dieldrin),
and it carries spurious phases for industrial chemicals (triglycidyl
isocyanurate, `Phase 2`) and food additives (sucrose, `Phase 3`). In v1.2 the
other-use test ran first and excluded those rows; v1.3 keeps them excluded for
the stated reason instead of by test order. When **no** environmental/industrial
use is documented, H1/H2 are sufficient on their own — that keeps approved salt
forms with sparse records (moxifloxacin hydrochloride, raloxifene hydrochloride,
citalopram hydrobromide, …) included by status, as in v1.2.

**Documented environmental/industrial use** is a non-pharmacological use in
PubChem Use Classification, PubChem Industry Uses, or the bare PubChem
**Use and Manufacturing → Uses** statements (HSDB/EPA "Uses" text), naming a
pesticide/agrochemical, an industrial chemical, a cosmetic/fragrance use, or a
food additive. The bare Uses statements are new in v1.3: v1.2 never read them, so
a pesticide with a mechanism text and no Use Classification string looked
undocumented. When several uses are documented the reported category follows the
declaration order pesticide/agrochemical, industrial, cosmetic/fragrance, food
additive, and the matched string is kept in the row.

**Mechanism/target evidence** is Broad moa or target, PubChem Mechanism of Action
text, PubChem Drug Targets, or a BindingDB ligand-target pair. **Research-tool use
evidence** is Broad Repurposing Hub membership. **Mechanism text alone never
satisfies eligibility:** D2 runs before D3, so a documented environmental or
industrial use excludes a compound even when a mechanism text and Broad Hub
membership are present, and an ATC code alone does not escape D2.

The other clauses that fix the old tier's failure modes are unchanged from v1.2:
status is never inferred from a non-empty ATC/MoA string; "unknown" is never an
exclude reason; an unresolved identity is never included or excluded.

### What is *not* evidence

- **An ATC code alone.** It is annotation, not a use: it does not escape D2 and is
  not a substitute for H1/H2/H3 (2-Phenylphenol, Cyfluthrin, Paraoxon and Dibutyl
  phthalate all carry human ATC codes and are now excluded by their documented
  pesticide/industrial use).
- **The generic PubChem "Pharmaceuticals" use-classification string.** It is not
  a human-use marker (Daminozide carries `Pharmaceuticals` beside a documented
  pesticide use).
- **PubChem Drug Indication / Therapeutic Uses text.** It is recorded in the table
  (columns `drug_indications`, `therapeutic_uses`) but not used as human-use
  evidence: PubChem carries diagnostic patch-test reagents ("approved for use
  within allergenic epicutaneous patch tests"), mis-mapped entries (the herbicide
  fenoxaprop-ethyl carries a human angina indication) and historical notes
  (paraoxon "formerly used for glaucoma"). Such a string does not corroborate a
  status claim.
- **Historical or experimental treatment notes.** 2,4-Dinitrophenol's 1930s
  anti-obesity use and paraoxon's former glaucoma use are recorded but are not a
  current, registered human indication.
- **The frozen `compounds.csv` ATC/MoA/MeSH strings** — annotation under audit,
  never eligibility evidence.
- **ChEMBL** — attempted, unreachable (HTTP 000/500); the attempt is recorded in
  `output/cache/chembl_probe.jsonl`.

### The review's named compounds, re-decided

| compound | v1.2 | v1.3 | deciding evidence (in the change report) |
|---|---|---|---|
| Ethoxyquin | include, mechanism_tool | **exclude**, other_use | documented `Farming (Pesticides)` use + food-feed antioxidant use; a mechanism text is not eligibility |
| Dibutyl phthalate | include, mechanism_tool | **exclude**, other_use | documented plasticiser / USEPA/OPP pesticide-code use |
| 2-Phenylphenol | include, mechanism_tool | **exclude**, other_use | documented `Fungicides` use |
| Daminozide | include, mechanism_tool | **exclude**, other_use | documented plant-growth-regulator/pesticide use |
| Trifloxystrobin | include, mechanism_tool | **exclude**, other_use | bare Uses: "Used as broad-spectrum fungicide …" |
| Cyfluthrin | include, mechanism_tool | **exclude**, other_use | documented agricultural insecticide use |
| Paraoxon | include, mechanism_tool | **exclude**, other_use | documented insecticide use; former glaucoma note is not a current indication |
| Quercetin | exclude, other_use | **include**, approved | PubChem Phase 3 / Broad `Launched` + human indications (fatigue, allergic rhinitis, drowsiness); the food-additive use is outranked |
| Genistein | include, clinical | include, clinical | unchanged: documented clinical phase (Phase 2 / Phase 2-Phase 3) |
| Carbaryl | exclude, other_use | exclude, other_use | `Approved` status is not corroborated as human drug use; documented insecticide use decides |
| Lindane | exclude, other_use | **include**, approved | PubChem `Approved` + Broad `Launched` + human indication `scabies` |

## Evidence sources (network reality on 2026-11-02)

| source | endpoint / file | status at build | used as |
|---|---|---|---|
| PubChem PUG-REST | `rest/pug/compound/inchikey/<key>/cids/JSON` | working (HTTP 200) | structure identity (CID) |
| PubChem PUG-View | `rest/pug_view/data/compound/<cid>/JSON` | working (HTTP 200) | Max Phase, ATC, Drug Classes, indications, therapeutic uses, MoA, targets, use classification, industry uses, bare Uses, synonyms |
| Broad Repurposing Hub | local Biomni parquet (`clinical_phase`, `moa`, `target`, `indication`, `disease_area`) | local file | documented phase + human indication + mechanism + tool-library membership |
| BindingDB 2024-09 | local 6.2 GB TSV, streamed by `awk` on `Ligand InChI Key` | local file | **recorded, not used** — the InChIKey join produced no usable per-key target pairs (deviation D5) |
| ChEMBL | `ebi.ac.uk/chembl/api/data/...` | **unreachable** (HTTP 000 / 500) | **not used** — the attempt is recorded in `output/cache/chembl_probe.jsonl` |
| frozen annotation | `compounds.csv` `atc`/`moa`/`mesh_classes` | frozen | annotation **under audit**, never eligibility evidence |

Every PubChem response is cached under `output/cache/` (a JSONL index with URL,
HTTP status and retrieval timestamp, plus one gzip body per URL), so `--offline`
rebuilds and re-verifies the record without the network. Local-source retrieval
timestamps come from the file mtime and `generated_at` is the latest retrieval in
the record, so an offline rebuild is byte-identical (verified: two consecutive
`--offline` builds produce the same sha256 for every artifact).

## Results (v1.3)

### Decisions

| | keys |
|---|---:|
| total structure identities | 460 |
| total name records | 467 |
| include | **118** |
| exclude | **314** |
| uncertain | **28** |

| decision class | keys |
|---|---:|
| approved | 107 |
| clinical | 9 |
| mechanism_tool | 2 |
| other_use | 314 |
| unknown | 28 |

| precedence label | keys |
|---|---:|
| human_drug_use | 116 |
| mechanism_tool | 2 |
| environmental_or_industrial_use | 314 |
| insufficient_evidence | 28 |

The include set is the candidate main training set; QC (ticket #2) still decides
which observations survive, so this is an eligibility count, not a well count.

The only mechanism tools admitted by D3 are **Cyclopamine** and **PHA-00568487**
(documented mechanism + Broad Hub research-library membership, no documented
environmental/industrial use). The seven v1.2 mechanism tools that were
agrochemical/industrial are excluded; see the change report.

Human-use coverage: 107 keys carry a human-use marker (97 human ATC, 82 human
indication; some carry both). 15 includes are approved/clinical **without** an
independent human ATC/indication marker — mostly salt forms such as moxifloxacin
hydrochloride, raloxifene hydrochloride, citalopram hydrobromide, plus
sematilide, azimilide, PK 11195, pirinixic acid, candoxatril, genistein. They have
no documented environmental/industrial use, so their status is not contradicted;
they are flagged `status_without_human_use_marker` rather than silently
reclassified. 24 excluded keys carry a documented phase (v1.2: 23; the revision
adds four status-bearing exclusions and re-includes three status-bearing natural
products).

### Change report v1.2 → v1.3 (127 rows)

| movement | keys | what happened |
|---|---:|---|
| include → exclude | **11** | 7 review-named agrochemicals/industrial compounds + 2,4-Dinitrophenol, Icaridin, Phenothiazine, Triglycidyl isocyanurate — documented other use, no corroborated human drug use |
| exclude → include | **3** | Quercetin, Lindane, Eugenol — documented human drug use outranks the documented food/fragrance/pesticide use |
| uncertain → exclude | **113** | the bare PubChem Uses statements now document a pesticide use (85) or an industrial use (28); v1.2 lacked that evidence |
| other cells | 0 | — |

Every moved row is in [`output/v1.3/change_report.csv`](output/v1.3/change_report.csv)
with its v1.2 decision, its v1.3 decision, the precedence label, the per-row
deciding evidence (human-use basis or other-use basis), the v1.3 rationale and the
evidence URLs. The 11 include→exclude rows:

| compound | v1.2 class | v1.3 | matched use evidence |
|---|---|---|---|
| Ethoxyquin | mechanism_tool | exclude | `Farming (Pesticides)` |
| Dibutyl phthalate | mechanism_tool | exclude | plasticiser / USEPA/OPP pesticide code |
| 2-Phenylphenol | mechanism_tool | exclude | `Fungicides` |
| Daminozide | mechanism_tool | exclude | plant growth regulator / pesticides |
| Trifloxystrobin | mechanism_tool | exclude | broad-spectrum fungicide |
| Cyfluthrin | mechanism_tool | exclude | agricultural insecticide |
| Paraoxon | mechanism_tool | exclude | insecticide (Merck Index) |
| 2,4-Dinitrophenol | clinical | exclude | pesticide / wood preservative / dye intermediate |
| Icaridin | clinical | exclude | USEPA/OPP pesticide code (insect repellent) |
| Phenothiazine | clinical | exclude | dye intermediate / insect control |
| Triglycidyl isocyanurate | clinical | exclude | powder-coating crosslinker |

The 3 exclude→include rows: **Quercetin** (approved; PubChem Phase 3 vs Broad
`Launched` phase conflict, human indications), **Lindane** (approved; human
indication `scabies`), **Eugenol** (approved; Broad `Launched`, human indication
`toothache`). Each row keeps both source phases and the `phase_conflict` flag.

### Uncertain / unresolved (28)

None is silently included or excluded; all are listed in `handoff.json` under
`uncertain` and `summary.json` `uncertain_keys`. 2 carry mechanism evidence but
no research-tool or human-use marker (4-Androstene-3,17-dione, CP-122721); 3
carry a phase string that is not an approval/clinical phase (two `Unknown`, one
`Preclinical`); the rest have no documented status, mechanism or use in the
sources reachable today. Several are pesticides whose PubChem record contains no
Use section (Thiophanate-methyl, Propamocarb hydrochloride,
Trifloxysulfuron-sodium): the rule does not fabricate an exclusion reason, so
they stay uncertain and out of the training set. The ticket-#3 rule is that an
uncertain key stays out of the main training set unless a new, recorded evidence
decision admits it.

### ATC normalisation — both readings, recomputed (unchanged from v1.2)

The frozen `atc` field is annotation under audit; these numbers are recomputed from
it (complete = `\b[A-Z][0-9]{2}[A-Z]{2}[0-9]{2}\b`; partial = 5-character codes that
are not the prefix of a complete code in the same string). Levels are the audit's
L1–L5 = ATC prefix lengths 1/3/4/5/7.

| reading | keys with ATC | complete codes | partial codes | L1 | L2 | L3 | L4 | L5 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| complete 7-character only | 106 | 172 | — | 14 | 58 | 88 | 124 | 172 |
| complete + partial | 106 | 172 | 19 | 14 | 60 | 92 | 134 | 172 |

The audit's prefix claims (`14/58/88/124/172` and `14/60/92/134/172`) **reproduce
exactly**. The audit's "22 partial codes" does not: the frozen field holds 19
distinct 5-character codes (9 keys carry one). ATC is a multi-label annotation —
Aspirin alone carries five complete codes — so a code count is not a compound count
and L3 is not a mutually exclusive class set. Separately, `atc_pubchem` holds the
substance-level codes documented in PubChem (97 keys, 152 codes); PubChem and the
frozen annotation share a complete code for 93 keys.

### MoA is a multi-label annotation string, not "canonical mechanisms" (unchanged)

At structure level the frozen field has **82 keys with a MoA**, 59 distinct complete
strings and 71 semicolon-separated fragments; at the original name-record level the
same field gives 87 rows, **60 strings and 72 fragments**. The table's `moa` column
is *documented* mechanism (Broad `moa` first, else a PubChem Mechanism-of-Action
excerpt) with `moa_sources` naming which; 259 keys have documented mechanism
evidence (117 Broad, 140 PubChem). The frozen string is kept separately as
`moa_frozen`.

### Structure identity (unchanged)

- PubChem CID resolved for **460/460** keys; identity resolved (PubChem, Broad or
  BindingDB) for all 460; **0 unresolved**.
- 451 frozen names appear among PubChem synonyms, 2 do not, 7 have no synonym list to
  check against.
- 27 keys have a registered structure that is a salt/parent of `drug_key`; all names
  of the 7 multi-name structures are kept in `all_names`, and no two structures are
  merged.
- 3 rows have a frozen SMILES whose recomputed InChIKey differs from the frozen
  `inchikey`; recorded in `identity_notes`, not "fixed".
- MeSH: 240 keys carry a frozen MeSH class; 52 of them are corroborated by a PubChem
  class/mechanism string.

### Old tier comparison (v1.3 decisions)

Old tiers are kept as flags; nothing was rewritten in the frozen corpus.

| old tier | old keys | still include | now exclude | now uncertain | false negatives | false positives |
|---|---:|---:|---:|---:|---:|---:|
| `atc-moa` | 113 | 105 | 7 | 1 | 13 | 7 |
| `atc-moa-mesh` | 252 | 114 | 134 | 4 | 4 | 134 |
| `all-smiles` | 460 | 118 | 314 | 28 | 0 | 314 |

`atc-moa` false negatives (13): Lindane, N-acetyl-procainamide, Sematilide,
PHA-00568487, Azimilide, Oxytetracycline hydrochloride, Cyclopamine, PK 11195,
Quercetin, Eugenol, Pirinixic acid, Genistein, Candoxatril — real omissions from
the old `bool(atc or moa)` tier. Eight v1.2 false negatives (2,4-Dinitrophenol,
2-Phenylphenol, Daminozide, Ethoxyquin, Icaridin, Phenothiazine, Trifloxystrobin,
Triglycidyl isocyanurate) are now correctly excluded, while Quercetin, Lindane and
Eugenol are new false negatives (v1.2 had 18, v1.3 has 13). `atc-moa` false
positives are 7 (v1.2: 2): Dibutyl phthalate, L-Tartaric acid, Dimethyl phthalate,
Ethylparaben, Cyfluthrin, Indoxacarb, Paraoxon — the revision adds five compounds
the old tier admitted and the evidence now excludes. The one "now uncertain" key
is Chlorhexidine diacetate. The `atc-moa-mesh` false positives (134) are the
environmental/other-use compounds the MeSH tier let in — that tier is not "more
drugs", it is more annotations.

## Handoff for ticket #3 (v2)

`output/v1.3/handoff.json`, schema **`prnet-eligibility-handoff-v2`**, sha256
`e52cac242978ca4f936024efaf24339d85a2cec687effe117ed212abd2eb2087`. It is the
interface ticket #3 consumes; `output/current.json` points at it (path, schema,
sha256) so ticket #3 need not hard-code the version directory:

```json
{"schema": "prnet-eligibility-handoff-v2",
 "artifact_version": "v1.3",
 "rule": {"version": "eligibility-rule-v1.3", "sha256": "504b8e42…"},
 "predecessor": {"artifact_version": "v1.2", "rule": {...}, "files": [...]},
 "counts": {"total_keys": 460, "include": 118, "exclude": 314, "uncertain": 28},
 "precedence_counts": {...},
 "include": [{"drug_key", "decision", "decision_class", "decision_precedence",
              "compound_name", "all_names", "inchikey", "pubchem_cid",
              "clinical_phase", "moa", "targets", "atc_codes",
              "human_use_basis", "other_use_basis", "evidence_urls",
              "retrieved_at", "rationale"}],
 "uncertain": [...], "exclude": [...]}
```

The include, uncertain and exclude lists partition all 460 keys; their sha256 and
the table's are in `output/v1.3/manifest.json`. Ticket #3 consumes the `include`
list plus the per-key metadata and keeps `uncertain` keys out of the main training
set.

### Artifact digests (v1.3)

| file | sha256 |
|---|---|
| `eligibility_table.csv` | `37485af6587d577ab24ff22e8a6222f3f9d4695d60fd1007d0c24d7f5fc4ba85` |
| `handoff.json` | `e52cac242978ca4f936024efaf24339d85a2cec687effe117ed212abd2eb2087` |
| `change_report.csv` | `46488e485cfe93d8aaef274d50a18262f7a72742b671623fbcb3473d6f1dfcc4` |
| `old_tier_comparison.csv` | `8135df438629329560627ac98f6fd096eadbbb24be6eeb777cbcfa2c167ee606` |
| `evidence_manifest.csv` | `fde46afdd4b6e4829cd9b9ee26d0b0518c0782cf852effc8742e0098bba77ecf` |
| `summary.json` | `19e1e93b2b7ce4f9100dc3f9aa00c52ae58d55cf8f0fd7ec980eb19e2bf6ffd9` |
| `manifest.json` | `f3243576d6f745588e143a6c75dbc7d5a74f8108bec4995d08c0f0808685c7f0` |
| `current.json` | `592ba6aeefd39f1214524094f4393766f556d1717dce857ef47c49b002e1b0fc` |

The frozen v1.2 record keeps its original digests (table `0663d725…`, handoff
`3aa2282b…`, manifest `9dd864d0…`); the v1.3 manifest's `predecessor` block
carries the v1.2 rule and file hashes, and the `v12_preserved` gate recomputes
them.

## Gates

`python execution/build_eligibility.py --check` re-runs on the committed record:

| gate | expectation |
|---|---|
| row_count / key_uniqueness / key_set_matches_frozen | 460 rows, keys equal the frozen key set |
| name_records | 467 frozen name records accounted for |
| decision_domain | decision ∈ {include, exclude, uncertain}; class ∈ the five |
| rule_identity_v13 | summary, handoff, manifest and current.json carry `eligibility-rule-v1.3` and its sha256 |
| precedence_recorded | every row's precedence label is consistent with its decision |
| include_records_human_use | every include has a recorded human-use basis |
| exclude_records_other_use | every exclude has a recorded other-use basis |
| mechanism_alone_never_includes | no mechanism tool without a use/tool marker |
| traceable_decisions | every include/exclude has a rationale and ≥1 HTTP-200/local evidence URL |
| unknown_never_excluded | no exclude with class `unknown` |
| exclude_requires_other_use | every exclude has a documented other-use flag |
| unresolved_identity_uncertain | unresolved identity never included/excluded |
| old_tier_flags | boolean old-tier flags on every row |
| retrieval_timestamps | every row has a retrieval timestamp |
| v12_preserved | the frozen v1.2 files still match their manifest hashes |
| change_report_complete | every moved row is reported once, with deciding evidence and its v1.2/v1.3 decisions |
| manifest_hashes | every tracked artifact sha256 matches the v1.3 manifest; the predecessor block matches v1.2 |
| handoff_schema_v2 | handoff schema is `prnet-eligibility-handoff-v2`; counts equal the summary |
| handoff_partition | include + uncertain + exclude = 460 |
| current_pointer | `current.json` points at the v1.3 record and matches its digests |

All **22 PASS** on the committed record (2026-11-02). The unittest suite
(`execution/tests/`, **66 tests**: 32 rule, 16 evidence, 18 record/gates) covers the
rule engine, ATC normalisation, evidence extraction, identity mapping and the same
completeness/precedence gates, all PASS.

## Reproduce

From the repository root:

```bash
bash stage2/06_prnet_transfer/06.1_eligibility/execution/run.sh            # fetch (cached), build, gates
bash stage2/06_prnet_transfer/06.1_eligibility/execution/run.sh --offline  # rebuild from cache (byte-identical)
bash stage2/06_prnet_transfer/06.1_eligibility/execution/run.sh --check    # gates only
make -C stage2/06_prnet_transfer/06.1_eligibility/execution test           # unittest
```

`--out-root` overrides the step directory that owns `output/`; `--bindingdb` points
at the BindingDB TSV (default: the local Biomni data lake). `--jobs` sets the
parallel PubChem fetch count (default 3).

## Environment and resources

`/home/chang/miniconda3/envs/arc/bin/python` (Python 3.13, pandas 2.3.3, rdkit
2026.03.2). No GPU, no training. `curl` is used for HTTP because this environment
cannot verify the site certificate chain through Python's `urllib`.

Measured for the offline build (2026-11-02, WSL2 envelope in
[`docs/agents/execution-environment.md`](../../../docs/agents/execution-environment.md)):
run under `systemd-run --user --scope -p MemoryMax=8G -p MemorySwapMax=8G` with
`/usr/bin/time -v`, **peak RSS 313,140 KiB (≈306 MiB)** on the final timed
rebuild (298,832 KiB on an earlier run), wall time ≈4 s, exit 0.
The build is single-process, reads BindingDB streamed (never whole) and materialises
no observation × gene matrix.

## Tracked vs local

Tracked (force-added per ADR 0030): everything under `output/v1.2/` and
`output/v1.3/` and `output/current.json`. Local only: `output/cache/` (PubChem
bodies + JSONL index + BindingDB match file), logs.

## Version history and deviations (D0–D6)

| version | date | change | status |
|---|---|---|---|
| v1.0 | 2026-11-02 | took a documented clinical phase as sufficient | rejected pre-final |
| v1.1 | 2026-11-02 | other-use test first + human-use escape | rejected pre-final |
| v1.2 | 2026-11-02 | dropped PubChem `Drug Classes` as tool evidence; locked; first built record (126/193/141) | **preserved in `output/v1.2/`** |
| v1.3 | 2026-11-02 | explicit precedence (human drug use > environmental/industrial use > mechanism text), bare Uses statements as use evidence, human-indication corroboration, re-decided all 460 keys | **current (118/314/28)** |

No model result existed at any revision: the table is built only from
documentation and tickets #4–#11 have not run, so the revision is still
pre-registration development, not result-driven tuning.

| id | deviation | why |
|---|---|---|
| D0 | ChEMBL is attempted and recorded, never used as evidence | the EBI API answered HTTP 000/500 on 2026-11-02; no substitute was invented |
| D1 | frozen ATC/MoA/MeSH are annotation under audit, not evidence; frozen and PubChem ATC are reported separately | the old `bool(atc or moa)` admitted `Unknown` and agrochemical mechanisms; v1.3 also refuses a bare ATC code as human-use evidence |
| D2 | three rule revisions before the locked v1.2 table (v1.0/v1.1/v1.2), then the v1.3 revision after review | no model result existed at any revision; the hash of each rule is in its artifacts |
| D3 | the audit's "22 partial ATC codes" is not reproduced (19 found) | recomputation from the frozen field; the discrepancy is reported rather than copied |
| D4 | PubChem Max Phase is the primary status source even where it conflicts with Broad (34 `phase_conflict` rows) | both phases and the conflict flag are kept; D1 requires a human-use marker when an environmental/industrial use is documented, which neutralises the ChEMBL-derived approvals of agrochemicals |
| D5 | BindingDB is attempted and recorded but supplied **no** target evidence in this build; PubChem Drug Indication / Therapeutic Uses are recorded but not used as human-use evidence | the BindingDB InChIKey join yielded no usable per-key pairs (`summary.json` sources), and PubChem indications include patch-test reagents and mis-mapped entries (fenoxaprop-ethyl) |
| D6 | v1.3 splits the record into `output/v1.2/` (frozen) and `output/v1.3/` (current) with `current.json` as the pointer, and the table schema is `…-v2` | ADR 0030: a new rule version gets a new artifact identity; the v1.2 record is never overwritten, and ticket #3 reads the pointer instead of a moving path |

## Out of scope

- No per-observation QC, no dose/vehicle/control decision, no well counts (ticket #2);
  eligibility ≠ measured response quality.
- No fold, split, input build or model claim (tickets #3–#11).
- No rewrite of the frozen corpus, `datasets-1` views or the old tiers; the new
  eligibility layer is a derived training input only.
- No claim that an included compound has a measurable CM response; a low-response but
  eligible compound is still eligible.
