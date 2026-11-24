# Initial matrix implementation checks (2026-11-16)

Base dev b8616338f354277c6d4a95f2a575e1cef0cd7e23; #230 integrated.
No formal matrix execution or replacement acceptance is claimed here.

Approved #226 seam TDD: comparison tests first failed because comparison API did
not exist; matrix tests first failed because run_matrix.py did not exist; figure
test first failed because CSV renderer did not exist. After each minimal slice,
8 new behavioral tests pass with arc Python 3.13: four comparison tests (OOF
identity/order, missing seed rejection, paired constant-difference CI and no
forced winner, nonfinite refusal), three matrix planning/resume tests (90 unique
identities and refusal of incomplete/duplicate/smoke/failed candidates), and one
subprocess CSV-only exact PNG redraw test without any model/input payloads.

`execution/run_matrix.py --plan-only` emitted 90 identities without model imports.
GPU visibility and systemd MemoryMax=10G/MemorySwapMax=8G were checked under
approved escalation: RTX4070TiSUPER, 16376MiB total, 1731MiB occupied at check.
Source trees, Stage3, old02 and ignored inputs were not modified. Full scientific
matrix acceptance, complete suite, resource/profile evidence and review remain
required before publication/retirement.
