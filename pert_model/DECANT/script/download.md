# Download — DECANT

**method: manual** — deployed model, but there is **no public data file**.

## Repository (automatic)

```bash
git clone https://github.com/bliulab/DECANT pert_model/DECANT/source
git -C pert_model/DECANT/source checkout e50642a980cbdf9fa73ec152dbe81a6047ab396b
# vendored read-only: the nested .git directory was removed on 2026-10-19
```

## 4. DECANT — no public data file

DECANT (deployed) explicitly excludes processed cell bags, teacher prototypes,
and checkpoints from its repository. There is no public download. The setup smoke
uses a synthetic 2,000-gene forward; reconstructing matched treated/control cell
bags is a separate data-engineering task.

Paper: <https://doi.org/10.1093/bioinformatics/btag662>. Deployed as a
representation baseline (dose is context to suppress, not a generative response
axis). Run card `stage2/01_setup/execution/decant/`.
