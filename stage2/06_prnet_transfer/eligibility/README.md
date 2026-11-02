# eligibility — 460-structure pharmacological eligibility audit

**Status: built and verified (2026-11-02).** Ticket #1 of the PRnet formal-transfer
effort. The audit covers every structure identity of the frozen temp_env panel and
replaces the old automatic `atc-moa` (113 keys) / `atc-moa-mesh` (252 keys) tiers as
the basis for the main training set. Nothing here trains a model; every number is
an evidence count, not a performance result.

The durable artifacts are the tracked record in [`output/`](output):
`eligibility_table.csv` (one row per structure), `old_tier_comparison.csv`,
`evidence_manifest.csv`, `summary.json`, `manifest.json`, and `handoff.json` — the
interface ticket #3 consumes.

## What this delivers

- **460/460 rows.** Every `drug_key` with `temp_env_all-smiles=True` in the frozen
  [`compounds.csv`](../../02_cm_data/output/datasets-1/compounds.csv) (467 name
  records collapse to 460 keys; no structure is merged or dropped). Each row carries
  all names, InChIKey, SMILES, decision, decision class, documented clinical phase,
  MoA, targets, ATC codes, MeSH, evidence URLs with retrieval timestamp, rationale,
  and the three old-tier flags.
- **A rule locked before the final table** (`execution/rules.py`; version and sha256
  in every artifact).
- **Coverage counts** recomputed from the data: class split, ATC normalisation in
  both readings, MoA annotation shape, identity resolution, old-tier false
  negatives/positives.
- **An include list, an unresolved list, and a machine-readable handoff** for
  ticket #3.

## The locked decision rule

`RULE_TEXT` in [`execution/rules.py`](execution/rules.py) is the canonical rule;
`rule_sha256()` hashes it and the hash is recorded in the table, summary, handoff
and manifest. **`eligibility-rule-v1.2`**, sha256
`e1af8cdb370bb1a8adc59ab49a5d658029f18e372e67248034fcccc2fd545bfa`.

The rule cannot be influenced by model results: no model is trained in this ticket,
the table is built only from documentation, and tickets #4–#11 have not run. The
rule text was finalized before the final table build; earlier drafts (v1.0, v1.1)
were rejected during evidence review, before any model result existed (D2).

| decision class | when it applies | decision |
|---|---|---|
| `approved` | PubChem Max Phase `Approved` or Broad `clinical_phase` `Launched` | include |
| `clinical` | documented phase (`Early Phase 1`, `Phase 1/2/3/4`, `Phase 1/Phase 2`, `Phase 2/Phase 3`, `Withdrawn`) in either source | include |
| `mechanism_tool` | mechanism/target evidence **and** research-tool or human-use evidence | include |
| `other_use` | documented non-pharmacological use (pesticide/agrochemical, industrial, cosmetic/fragrance/flavour, food additive) **and** no documented pharmaceutical use | exclude |
| `unknown` | nothing documented | uncertain |

Decision procedure, first match wins: **D0** identity gate (unresolved ⇒ uncertain),
**D1** other_use ⇒ exclude, **D2** approved ⇒ include, **D3** clinical ⇒ include,
**D4** mechanism_tool ⇒ include, **D5** unknown ⇒ uncertain.

The clauses that fix the old tier's failure mode:

- **Status is never inferred from a non-empty ATC/MoA string.** The old tier was
  `bool(atc or moa)`, which admitted `moa = Unknown`; here approval/clinical status
  comes only from PubChem Max Phase or Broad `clinical_phase`.
- **A pesticide/herbicide mechanism does not qualify a compound.** D1 runs before
  D4: a documented agrochemical/industrial use with no documented pharmaceutical use
  is an exclude even when an aggregator records an approval for it (Carbaryl,
  Lindane) or when PubChem has a mechanism text (Indoxacarb, Fipronil).
- **Two markers, deliberately different.** The D1 escape uses a 7-character human
  ATC code or a PubChem use string naming pharmaceutical use (this keeps a real drug
  that also carries a fragrance/food use, e.g. methyl salicylate, out of the
  exclude). The D4 qualifier uses an ATC code only, because the generic
  "Pharmaceuticals" use string would promote agrochemicals with a mechanism text.
- **"unknown" is never an exclude reason.** Keys with no documented status,
  mechanism or use stay `uncertain` and are listed separately.

Identity gate: a key is resolved when the frozen InChIKey maps to a PubChem CID or
to a Broad/BindingDB InChIKey record. Unresolved identity ⇒ `uncertain`, whatever
else is documented.

## Evidence sources (network reality on 2026-11-02)

| source | endpoint / file | status at build | used as |
|---|---|---|---|
| PubChem PUG-REST | `rest/pug/compound/inchikey/<key>/cids/JSON` | working (HTTP 200) | structure identity (CID) |
| PubChem PUG-View | `rest/pug_view/data/compound/<cid>/JSON` | working (HTTP 200) | Max Phase, ATC, Drug Classes, indications, MoA, targets, use classification, synonyms |
| Broad Repurposing Hub | local Biomni parquet (`clinical_phase`, `moa`, `target`) | local file | documented phase + mechanism + tool-library membership |
| BindingDB 2024-09 | local 6.2 GB TSV, streamed by `awk` on `Ligand InChI Key` | local file | supplementary target evidence |
| ChEMBL | `ebi.ac.uk/chembl/api/data/...` | **unreachable** (HTTP 000 / 500) | **not used** — the attempt is recorded in `output/cache/chembl_probe.jsonl` |
| frozen annotation | `compounds.csv` `atc`/`moa`/`mesh_classes` | frozen | annotation **under audit**, never eligibility evidence |

Every PubChem response is cached under `output/cache/` (a JSONL index with URL, HTTP
status and retrieval timestamp, plus one gzip body per URL), so `--offline` rebuilds
and re-verifies the record without the network.

## Results

### Decisions

| | keys |
|---|---:|
| total structure identities | 460 |
| total name records | 467 |
| include | **126** |
| exclude | **193** |
| uncertain | **141** |

| decision class | keys |
|---|---:|
| approved | 104 |
| clinical | 13 |
| mechanism_tool | 9 |
| other_use | 193 |
| unknown | 141 |

The include set is the candidate main training set; QC (ticket #2) still decides
which observations survive, so this is an eligibility count, not a well count.

Status detail: 97 keys are `Approved` by PubChem Max Phase and 97 by Broad
`Launched` (union 104 approved keys); 34 rows have a PubChem–Broad phase conflict
(kept as `phase_conflict`, plus both phase columns). 23 keys carry a documented
phase but are excluded by D1 (Carbaryl, Lindane, Dieldrin, Bisphenol A, Sucrose,
Quercetin, …) — the documented agrochemical/industrial use is the reason and is in
the row. 26 keys are approved/clinical without an independent ATC marker (mostly
salt forms such as moxifloxacin hydrochloride, plus PK 11195, sematilide,
candoxatril); they are flagged `status_without_human_use_marker` rather than
silently reclassified.

### ATC normalisation — both readings, recomputed

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
frozen annotation share a complete code for 93 keys. PubChem also lists ATC
hierarchy levels as separate entries; those are not counted, so the PubChem "partial"
reading is reported as 0 rather than inflated.

### MoA is a multi-label annotation string, not "canonical mechanisms"

At structure level the frozen field has **82 keys with a MoA**, 59 distinct complete
strings and 71 semicolon-separated fragments; at the original name-record level the
same field gives 87 rows, **60 strings and 72 fragments** — the audit's 60/72 are
record-level counts, and neither is a count of canonical mechanisms (semicolons also
join target/action phrases). The table's `moa` column is *documented* mechanism
(Broad `moa` first, else a PubChem Mechanism-of-Action excerpt) with `moa_sources`
naming which; 259 keys have documented mechanism evidence (117 Broad, 140 PubChem).
The frozen string is kept separately as `moa_frozen` (e.g. Coumarin's frozen `Unknown`
is replaced by the documented "vitamin K antagonist").

### Structure identity

- PubChem CID resolved for **460/460** keys; identity resolved (PubChem, Broad or
  BindingDB) for all 460; **0 unresolved**.
- 451 frozen names appear among PubChem synonyms, 2 do not, 7 have no synonym list to
  check against.
- 27 keys have a registered structure that is a salt/parent of `drug_key` (the
  `drug_key` is the desalted parent by frozen contract); all names of the 7
  multi-name structures are kept in `all_names`, and no two structures are merged.
- 3 rows have a frozen SMILES whose recomputed InChIKey differs from the frozen
  `inchikey`; recorded in `identity_notes`, not "fixed".
- MeSH: 240 keys carry a frozen MeSH class; 52 of them are corroborated by a PubChem
  class/mechanism string. The rest stay annotation under audit.

### Old tier comparison

Old tiers are kept as flags; nothing was rewritten in the frozen corpus.

| old tier | old keys | still include | now exclude | now uncertain | false negatives | false positives |
|---|---:|---:|---:|---:|---:|---:|
| `atc-moa` | 113 | 108 | 2 | 3 | 18 | 2 |
| `atc-moa-mesh` | 252 | 121 | 94 | 37 | 5 | 94 |
| `all-smiles` | 460 | 126 | 193 | 141 | 0 | 193 |

- `atc-moa` false negatives (18) are the real omissions: Cyclopamine, Genistein,
  Candoxatril, Azimilide, Sematilide, Oxytetracycline hydrochloride, PK 11195,
  Pirinixic acid, and others with no ATC/MoA string in the frozen annotation.
- `atc-moa` false positives (2) and the 3 "now uncertain" keys are Indoxacarb
  (moa `External insecticide`), Methyl salicylate and Coumarin added by the old
  `bool(atc or moa)`; Chlorhexidine diacetate, Dimethyl phthalate and Ethylparaben
  were admitted by ATC alone but have no documented mechanism or use in the sources
  reachable today.
- `atc-moa-mesh` false positives (94) are the environmental/other-use compounds the
  MeSH tier let in; that tier is not "more drugs", it is more annotations.

Named audit candidates and their outcomes (evidence in the row):

| candidate | outcome | basis |
|---|---|---|
| Cyclopamine | include, mechanism_tool | Broad `Preclinical`, hedgehog pathway inhibitor |
| Genistein | include, clinical | Broad `Phase 2/Phase 3` |
| Candoxatril | include, clinical | PubChem Max Phase 2 |
| Azimilide | include, clinical | PubChem/Broad Phase 3 |
| Sematilide | include, clinical | PubChem Max Phase 2 |
| Oxytetracycline hydrochloride | include, approved | PubChem Max Phase `Approved` |
| PK 11195 | include, clinical | Broad `Phase 1` |
| Pirinixic acid | include, clinical | PubChem Max Phase 2 |
| Methyl salicylate | include, approved | PubChem Max Phase `Approved` |
| Coumarin | include, approved | PubChem Max Phase `Approved` |
| Indoxacarb | exclude, other_use | documented `Insecticides` use, no ATC/indication |

### Uncertain / unresolved (141)

- 46 have mechanism evidence but no research-tool or human-use marker (4-nitrophenol,
  imazapyr, lactofen, tebuthiuron, nitrobenzene, …).
- 52 carry only use strings (often both "Pharmaceuticals" and a pesticide entry), and
  3 carry an ATC code but no mechanism or status (Chlorhexidine diacetate, Dimethyl
  phthalate, Ethylparaben).
- 64 have no documented status, mechanism or use in the reachable sources.

None is silently included or excluded; all are listed in `handoff.json` under
`uncertain`. The ticket-#3 rule is that an uncertain key stays out of the main
training set unless a new, recorded evidence decision admits it.

### Residual includes a reviewer may want to revisit

9 mechanism tools are admitted by D4: Cyclopamine, PHA-00568487 (research tools),
2-Phenylphenol (ATC D08AE06), Cyfluthrin (P03BA01), Paraoxon (S01EB10) — human ATC
markers — and Ethoxyquin, Dibutyl phthalate, Daminozide, Trifloxystrobin, which are
Broad Repurposing Hub `Preclinical` entries with a documented MoA. The Broad entries
for the last four look like library/tool entries for agrochemical or industrial
compounds; they are included because the locked rule names Broad membership as
research-tool evidence, and they are named here so a reclassification (if wanted) is
a recorded change, not a silent one.

## Handoff for ticket #3

`output/handoff.json` (`schema: prnet-eligibility-handoff-v1`) is the interface:

```json
{"rule": {"version": "...", "sha256": "..."},
 "counts": {"total_keys": 460, "include": 126, "exclude": 193, "uncertain": 141},
 "include": [{"drug_key", "decision", "decision_class", "clinical_phase", "moa",
              "targets", "atc_codes", "evidence_urls", "retrieved_at", "rationale"}],
 "uncertain": [...], "exclude": [...]}
```

The include, uncertain and exclude lists partition all 460 keys; their sha256 and
the table's are in `output/manifest.json`. Ticket #3 consumes the `include` list plus
the per-key metadata and keeps `uncertain` keys out of the main training set.

## Gates

`python execution/build_eligibility.py --check` re-runs on the committed record:

| gate | expectation |
|---|---|
| row_count / key_uniqueness / key_set_matches_frozen | 460 rows, keys equal the frozen key set |
| name_records | 467 frozen name records accounted for |
| decision_domain | decision ∈ {include, exclude, uncertain}; class ∈ the five |
| traceable_decisions | every include/exclude has a rationale and ≥1 HTTP-200/local evidence URL |
| unknown_never_excluded | no exclude with class `unknown` |
| exclude_requires_other_use | every exclude has a documented other-use flag |
| unresolved_identity_uncertain | unresolved identity never included/excluded |
| old_tier_flags | boolean old-tier flags on every row |
| retrieval_timestamps | every row has a retrieval timestamp |
| handoff_hash / table_hash | sha256 matches `manifest.json` |
| handoff_partition | include + uncertain + exclude = 460 |

All PASS on the committed record (2026-11-02). The unittest suite
(`execution/tests/`) covers the rule engine, ATC normalisation, identity mapping and
the same completeness gates.

## Reproduce

From the repository root:

```bash
bash stage2/06_prnet_transfer/eligibility/execution/run.sh            # fetch (cached), build, gates
bash stage2/06_prnet_transfer/eligibility/execution/run.sh --offline  # rebuild from cache
bash stage2/06_prnet_transfer/eligibility/execution/run.sh --check    # gates only
make -C stage2/06_prnet_transfer/eligibility/execution test           # unittest
```

`--out-root` overrides the step directory that owns `output/`; `--bindingdb` points
at the BindingDB TSV (default: the local Biomni data lake). `--jobs` sets the
parallel PubChem fetch count (default 3).

## Environment

`/home/chang/miniconda3/envs/arc/bin/python` (Python 3.13, pandas 2.3.3, rdkit
2026.03.2). No GPU, no training. `curl` is used for HTTP because this environment
cannot verify the site certificate chain through Python's `urllib`.

## Tracked vs local

Tracked (force-added per ADR 0030 precedent): `output/eligibility_table.csv`,
`output/old_tier_comparison.csv`, `output/evidence_manifest.csv`,
`output/summary.json`, `output/handoff.json`, `output/manifest.json`.
Local only: `output/cache/` (PubChem bodies + JSONL index + BindingDB match file),
logs.

## Deviations (D0–D4)

| id | deviation | why |
|---|---|---|
| D0 | ChEMBL is attempted and recorded, never used as evidence | the EBI API answered HTTP 000/500 on 2026-11-02; no substitute was invented, so keys whose only phase/mechanism source would have been ChEMBL stay `uncertain` |
| D1 | frozen ATC/MoA/MeSH are annotation under audit, not evidence; frozen and PubChem ATC are reported separately | the old `bool(atc or moa)` admitted `Unknown` and agrochemical mechanisms |
| D2 | three rule revisions before the final table: v1.0 trusted Max Phase (Carbaryl/Dieldrin admitted), v1.1 used PubChem Drug Classes as tool evidence (agrochemical class strings admitted), v1.2 is the locked rule | no model result existed at any revision, so this is pre-registration development, not result-driven tuning; the hash of the final rule is in every artifact |
| D3 | the audit's "22 partial ATC codes" is not reproduced (19 found) | recomputation from the frozen field; the discrepancy is reported rather than copied |
| D4 | PubChem Max Phase is the primary status source even where it conflicts with Broad (34 `phase_conflict` rows) and even where an aggregator marks an agrochemical approved | it is the only documented machine-readable approval field reachable without ChEMBL; both phases and the conflict flag are kept, and D1 overrides it for documented agrochemical/industrial use |

## Out of scope

- No per-observation QC, no dose/vehicle/control decision, no well counts (ticket #2);
  eligibility ≠ measured response quality.
- No fold, split, input build or model claim (tickets #3–#11).
- No rewrite of the frozen corpus, `datasets-1` views or the old tiers; the new
  eligibility layer is a derived training input only.
- No claim that an included compound has a measurable CM response; a low-response but
  eligible compound is still eligible.
