# Download — ExPO

**method: repo only** — upstream clone (MIT). No pretrained weights are released.

## Repository

```bash
git clone https://github.com/MLBC-lab/ExPO pert_model/ExPO/source
git -C pert_model/ExPO/source checkout 7e88cbae17b44397ff7d75983475334b40488d65
# vendored read-only: the nested .git directory was removed on 2026-10-19
```

## Notes

Reference-only: 978 landmark genes is the published scope (fails revised G1), so
it is not a Stage 2 entrant. Upstream packaging bug: `expo/models/expo_model.py`
has a cp1252 byte (0x96) at line 19 that breaks `import expo` under Python 3's
UTF-8 default; re-encode that file in a scratch copy to run the tests.
