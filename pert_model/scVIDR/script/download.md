# Download — scVIDR

**method: manual** — deployed model. The repository has code only; the
figure/notebook data and pretrained single-/multi-dose models are in the author's
Google Drive archive.

## Repository (automatic)

```bash
git clone https://github.com/BhattacharyaLab/scVIDR pert_model/scVIDR/source
git -C pert_model/scVIDR/source checkout bcb6a96c85cddc4eef052be739c122efb737fe0e
# vendored read-only: the nested .git directory was removed on 2026-10-19
```

## 3. scVIDR — `data.zip` (deployed; includes pretrained models)

The repository has code only; the figure/notebook data and pretrained
single-/multi-dose models are in the author's Google Drive archive.

**Source:** <https://drive.google.com/file/d/11fzDbp0B19Dy47MtD742Jl4Hz2bdSiiq/view?usp=sharing>

1. Download `data.zip`.
2. Put it in `pert_model/scVIDR/source/` and extract so the notebooks see
   `source/data/`.
3. The archive contains `nault2021_singleDose.h5ad`,
   `nault2021_multiDose.h5ad`, and the reported `.pt` model directories.
4. Record the archive SHA256 in `author_log.json` sources[].

The 2,000-gene synthetic smoke in `stage2/01_setup/scvidr/` does **not** need
this archive; it is only required for the published liver dose-response results.
The archive was not downloaded on this host, so `source/scvidr_data/` currently
holds only the empty upstream placeholder notebook.
