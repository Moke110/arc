# Heliste source-first corpus implementation — #52 / #53

The current step replaces the endpoint-gated workflow under ADR 0044. Baseline and
H-R each contain 2,732 admitted source values and 24 pending values: all 5,512
numerical Dataset S1 cells are accounted for. There are 689 source-label rows,
683 resolved labels and 682 unique reference canonical structures. Both Amiloride
hydrochloride labels remain separate despite resolving to one reference structure.

Pending chemical identities are Verteporfin (Visudyne), Sevelamer HCl, Quinine
hydrochloride dihydrate, Methscopolamine (Pamine), Formoterol hemifumarate and
Vitamin D2. The reports retain statuses, reasons, original cells/values and raw
lookup evidence. Sevelamer has no confirmed unique structure; the other five have
ambiguous reference results. No structure was fabricated to admit these values.

Chemical acquisition is separate from the offline corpus command. The declared
`pheno_data/Heliste-2021/local/chemical_resolution_v1.json` freezes raw PubChem
property/synonym responses, URLs, CIDs, retrieval times and deterministic RDKit
2026.03.2 canonicalization. Validation requires one CID, exact supplied primary
name/synonym confirmation and preservation of explicitly supplied salt and
stereochemistry. A full-label lookup for Sulindac initially returned an impurity;
the primary-name guard rejected it and the final lookup uses the confirmed
Sulindac reference. Author-supplied labels remain unchanged in the corpus.

Source limitations persist: author-normalized CTG is not raw RLU; the exact
normalization denominator, drug persistence and physical replicate identity are
unknown. The endpoint follows approximately 96 h elapsed, including 24 h
pretreatment, 48 h hypoxia and 24 h reoxygenation; normoxia is duration matched.
No vehicle values, replicate indices, SD, n, p-values or endpoint labels were
invented. These are source-measurement corpora, with no training readiness claim.

Validation evidence is in `validation/`: public build and offline reconstructive
verify pass; two public acceptance tests pass, including literal source-cell
values, CSV/Parquet equivalence, complete accounting, idempotence, rehashed output
tampering and rehashed raw identity-provenance tampering. Tests were first red
before the builder existed. The final measured peak RSS is 208,340 KiB for build,
256,088 KiB for verify and 259,816 KiB for the test command. Heliste source
verification and repository path/tracking checks pass. The full repository
integrity sweep is an integration-run responsibility because unrelated training
may still be active; no whole-repository integrity success is implied here.

Run commands from the repository root with the `arc` Python environment:

```bash
python stage3/01_corpus/execution/corpus.py build
python stage3/01_corpus/execution/corpus.py verify
python -m unittest discover -s stage3/01_corpus/execution/tests -v
bash pheno_data/Heliste-2021/script/verify.sh
bash script/check_paths.sh
bash script/check_tracked.sh
```

The unchanged prior release evidence lives under `execution/historical/`;
retired code is recoverable at commit `74f7e7a`. Standards/Spec review and tracker
resolution are recorded by the integration run after this implementation commit.
