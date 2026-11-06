# 02_base_morgan — structure and dose baseline

This step consumes the frozen [01_corpus](../01_corpus/README.md) release under
[ADR 0045](../../docs/adr/0045-stage3-morgan-direction-baseline.md). Implementation is in
progress under issue #59; the current public tracer bullet prepares the three
cold-drug folds. It does not constitute a completed benchmark.

```bash
/home/chang/miniconda3/envs/arc/bin/python stage3/02_base_morgan/execution/workflow.py prepare --input stage3/01_corpus/output --output /tmp/morgan-prepared
/home/chang/miniconda3/envs/arc/bin/python -m unittest discover -s stage3/02_base_morgan/execution/tests -p test_workflow.py -v
```

`--input` accepts the corpus release directory (condition subdirectories containing
`corpus.csv`) or one CSV with the same fields for synthetic acceptance. Outputs
are deterministic under input row reordering. `prepared.json` retains individual
source measurements, canonical parent structures, InChIKey parent identities,
source identities, directions and doses. Its roles align with the saved row order;
all aliases and doses of a parent receive the same role across conditions.

The split seed is 20261105. Three outer folds test each parent once; validation
holds out the ceiling of 20% of each remaining parent pool (at least one parent).
Preparation fails when any condition–phenotype target lacks a role. Input and
split fingerprints use SHA256 of canonical sorted-key JSON. Reusing an output
requires byte-identical preparation; choose a new directory for changed inputs.
Morgan metadata freezes radius 2, 2048 binary bits, chirality and RDKit version.
No fingerprint cache or training is produced by prepare yet.
