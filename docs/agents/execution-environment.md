# Execution environment (WSL2)

Measured inside the development WSL2 distro on 2026-11-26 (Asia/Shanghai). Every agent that
loads data, trains, or runs the test suite must respect this envelope: one VM with one shared
memory pool backs every worktree and subagent, so one unconstrained process can take the whole
session — and every other agent's work — down with it. This file records the measured limits
and the rules that go with them; re-measure after any WSL restart and update the numbers.

## Measured resources

| Resource | Measured value | Notes |
|---|---|---|
| CPU in WSL | 8 vCPU | `nproc`; `.wslconfig` sets `processors=8` |
| RAM, 2026-11-22 01:39 boot | **23.47 GiB** | `MemTotal: 24,612,504 KiB`. Use live `MemAvailable` for scheduling. |
| Swap | 64 GiB | `SwapTotal: 67,108,864 KiB`. Heavy swap activity can slow the VM. |
| GPU | NVIDIA GeForce RTX 4070 Ti SUPER, 16,376 MiB | 2,446 MiB used at measurement time. Query GPU usage before training. |
| Disk | 1.2 TiB, **66 GiB free** | ext4 on `/dev/sdd`, 95% used at measurement time. |
| Windows host | 63.77 GiB RAM total | 25.55 GiB available at measurement time. |

Use live WSL and Windows host memory readings when budgeting jobs. Changes to
`.wslconfig` take effect after WSL restarts.

## Incident: 2026-11-02 WSL crash during `/implement-spec`

The run had six implementer/merge subagents working concurrently (ticket #1 rule revision,
ticket #3 input build, ticket #4 pipeline, ticket #8 baselines). The journal for boot `-2`
ends abruptly at 12:33:50 with no shutdown sequence; that boot started 1,086
`dsh-subprocess` scopes and a single scope peaked at **9.7 GiB** (`4min 33s` CPU). That
scope ran `script/verify_all.sh` over the 52 GB of frozen external data, so most of the
peak is page cache charged to the scope rather than anonymous memory. No `oom-kill` line
was captured: the VM died rather than killing a userspace process. The host survived.

`.wslconfig` was then edited to give WSL 48 GB and a 64 GB swap file; the distro was restarted
cleanly at 14:31 and again at 14:35 (systemd `exit.target` in the journal, no crash). The
12:33 crash occurred with 31.3 GiB of WSL RAM. The 14:35 boot reported 47.0 GiB of WSL RAM
and `hv_balloon: Max. dynamic memory size: 49152 MB`. The 12:33:50 boot end is the only
one in this incident without a shutdown sequence.

## Rules

1. **At most two heavy jobs repo-wide.** All worktrees and subagents share this VM.
   Compare the jobs' measured memory peaks with available WSL and host RAM, leaving
   headroom for page cache and the agent runtime. Use one heavy job if the peaks
   are unknown or there is insufficient headroom for both. Run long training jobs alone.
2. **Cap every heavy process.** Preferred: `systemd-run --user --scope -p MemoryMax=10G
   -p MemorySwapMax=8G -- <cmd>`; use `MemoryMax=18G` for a single formal pretrain (measured
   15.1 GiB) and `MemoryMax=12G` for the 06.3/06.4 data jobs. Portable fallback:
   `ulimit -v` in the launching shell. A cap that trips is a signal to restructure the job,
   not to raise the cap.
3. **Pre-flight check.** `free -g` before launching; abort if `available` is below 12 GiB.
   After any suspicious exit, check `free -h` and `dmesg | tail` before relaunching.
4. **Bounded data access.** Never materialise an observation × 18,164 dense matrix. Never
   open a large `.h5ad` with `backed="r"` (see the trap below): read headers with `h5py` and
   slice through `read_elem`/CSR row runs, use `float32`, iterate per fold / dose / plate, and
   memmap large caches to `.npy`. Never `n_jobs=-1`; use `n_jobs=2`.
5. **Thread caps when parallel.** With two jobs running:
   `OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2`, and
   `torch.set_num_threads(2)`.
6. **Measure and record.** Run anything above ~1 GiB under `/usr/bin/time -v` (or rely on
   the systemd-run scope accounting) and put the peak RSS in the ticket record. A rising
   peak in a test suite or build gate is a defect.
7. **Tests are not free.** Run test files individually (`-p 'test_x.py'`) instead of the
   full suite alongside another heavy job; the 12:33 crash happened with six implementers
   each running suites, not with two.

## After a WSL restart

```bash
wsl.exe --shutdown          # from Windows, not from inside the distro
free -h
rg '^btime' /proc/stat
nvidia-smi --query-gpu=memory.total,memory.used --format=csv
df -h /home/chang/arc
```

Re-measure WSL RAM, swap, GPU, disk, and Windows host RAM. Update the table
for the new boot. If the RAM differs from the intended `.wslconfig` limit,
check that WSL restarted and that the file has the expected section and encoding.

## Known heavy jobs (measured)

| job | measured peak | consequence |
|---|---|---|
| `06.1_eligibility` offline build | 306 MiB | cheap |
| `05_prnet_pretrain` 68-epoch GPU pretrain on the 67,018 × 18,164 view | **15,150 MiB host RSS** | a full pretrain must run alone; leave at least 16 GiB free before starting |
| `script/verify_all.sh` (27 offline verifiers, 52 GB hashed) | 9.7 GiB charged, `4min 33s` CPU | mostly page cache; run it when no heavy job is running, not concurrently |
| `06.4_pipeline` real-data smoke (ticket #4 Phase B) | measure per run and record here | the Tahoe view is 16.2 GB; the loader must be header/memmap-based |
| `06.8_baselines --all` on the synthetic fixture (post-fix) | 1.2-1.7 GiB | dominated by imports and the shared evaluator, not the data |
| `06.7_cm_tahoe` bounded wrapper profile (1 CM unit, 2 epochs, 2 workers) | **2,293 MiB host process tree**, 1,583 MiB scope cgroup, 1,325 MiB run manifest | historical measurement beside the #5 pretrain; a 45-unit matrix is expected to match #6's ~3.0 GiB |
| full test suites of the 06 substeps | < 1 GiB each | two suites may run together, a third heavy job may not |

The pretrain number is the reason the training tickets (#5–#7) are serialised: two
concurrent 15 GiB trainers would consume the whole VM.

## Measured trap: `anndata.read_h5ad(backed="r")` is not bounded

`backed="r"` keeps `X` lazy but **still materialises every `layers` entry**. Measured on the
06.3 inputs: 8.1 GB for the 16.2 GB Tahoe view and 12.9 GB for the 25.8 GB CM corpus — the
symptom is a SIGKILL (exit 137) that looks like an OOM. Use `h5py` header reads plus
`read_elem`/CSR row slices for metadata and per-row access; ticket #3's
`contract.read_h5ad_header()` is the pattern, with regression tests that poison
`anndata.read_h5ad`. The frozen 06.3 output is 16.4 GB on disk; do not load it whole.
