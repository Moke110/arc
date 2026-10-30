# AGENTS.md — drug_emb_model

Drug biological-embedding model libraries as two external 项 (ADR 0030). Unlike the other
external categories these are **code libraries**: their `source/` trees stay tracked; only
the downloaded feature tensors and pretrained weights are gitignored
(`source/feature_extract/`, `source/models/`).

- `HG/` — HeteroGNN drug biological-embedding library, author-provided (ADR 0006), version `unknown`; `script/download.md`, `script/verify.sh`, `author_log.json`
- `KPGT/` — vendored KPGT knowledge-guided drug-embedding model, Apache-2.0, upstream URL/commit `unknown`; downloaded `source/models/` (894 MB weights) + `source/feature_extract/` (34 MB tensors); `script/download.md`, `script/verify.sh`, `author_log.json`

`source/` holds the author-provided/vendored library as given; `local/` holds declared local
additions (HG's `AGENTS.md` routing note; KPGT has none). Every item carries `manifest.sha256`
(442 GB hashed, 27/27 verified); `script/verify.sh` checks it offline. Category overview and
the Stage 1 HG-embedding role: `README.md`. Layout contract: ADR 0030.
