# Pretraining review corrections (2026-11-16)

Addresses all five findings from the standards/spec review of 9e4de08..5b8b477.
No formal training, publication or old-step retirement was performed for this fix.

Behavioral checks (arc Python 3.13, two threads):
- Four existing comparison tests pass, preserving OOF/seed/pairing behavior.
- Four matrix tests pass, including a new red→green interruption/replay test:
  interrupted bytes survive in attempts; awaiting replay invokes replay only;
  complete corrupt/smoke candidates still fail closed.
- One CSV-only figure test passes, including all-method N/A and constant/undefined
  annotations from CSV and byte-identical fresh-subprocess PNG redraw.
- Four new boundary tests pass: pristine copied native trees (including TranSiGen
  utils) retain every file byte/name after two fresh-process native loads;
  publication rejects all destinations before changes and survives staging-copy,
  second-rename and later PNG-rename failure with old bytes intact and retry
  success; a simulated journal interruption restores only journaled destinations;
  fake free-g available=11 prevents the heavy command from running.

The publication and memory tests first demonstrated failures in the previous
implementation; missing recovery/annotation seams first failed before their
implementation. Native-import test setup explicitly redirects third-party Numba
and Matplotlib caches to /tmp, then checks actual native source snapshots with
PYTHONDONTWRITEBYTECODE unset in the child environment. No scientific/native
configuration, seed, inference definition or stopping rule changed.

Source cache inventory (read-only): main native trees already contain PRnet
models/__pycache__/{PRnet.cpython-39,__init__.cpython-39,PRnet.cpython-313}.pyc and
TranSiGen src/__pycache__/{model.cpython-313,utils.cpython-313}.pyc. No files were
removed. Python 3.9 cache files may be upstream payload bytes; any later cleanup
must first distinguish downloaded provenance from generated 3.13 cache bytes.
