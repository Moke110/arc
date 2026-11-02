# Execution environment (WSL2)

Measured inside the development WSL2 distro on 2026-11-02 (Asia/Shanghai). Every agent that
loads data, trains, or runs the test suite must respect this envelope: the VM has no cgroup
memory cap, so one unconstrained process can take the whole session — and every other
agent's work — down with it. This file records the measured limits and the rules that go
with them; re-measure after any WSL restart and update the numbers.

## Measured resources

| Resource | Measured value | Notes |
|---|---|---|
| CPU in WSL | 8 vCPU | `nproc`; `.wslconfig` sets `processors=8` |
| RAM, this boot | **31.3 GiB** | `MemTotal: 32,870,068 kB`; `memory.max = max` (no per-cgroup limit) |
| RAM, requested for next boot | 48 GiB | `.wslconfig` `memory=48GB`, edited 2026-11-02 12:40; **only applies after `wsl --shutdown`** |
| Swap | 64 GiB | `swapFile=D:\temp\wsl-swap.vhdx` — on a Windows drive, so heavy swap thrash is slow |
| GPU | NVIDIA RTX 4070 Ti SUPER, 16,376 MiB | ~1.3 GiB already in use at idle |
| Disk | 1.2 TB, **222 GB free** | ext4 on `/dev/sdd`, 81% used |
| Windows host | 63.8 GiB RAM total | WSL can never exceed the host, whatever `.wslconfig` says |

The current boot still runs the old default envelope (50% of host RAM = 32 GiB). Do not
budget against the requested 48 GiB until `free -h` actually reports it.

## Incident: 2026-11-02 WSL crash during `/implement-spec`

The run had six implementer/merge subagents working concurrently (ticket #1 rule revision,
ticket #3 input build, ticket #4 pipeline, ticket #8 baselines). The journal for boot `-2`
ends abruptly at 12:33:50 with no shutdown sequence; that boot started 1,086
`dsh-subprocess` scopes and a single scope peaked at **9.7 GiB** (`4min 33s` CPU). That
scope ran `script/verify_all.sh` over the 52 GB of frozen external data, so most of the
peak is page cache charged to the scope rather than anonymous memory. No `oom-kill` line
was captured: the VM died rather than killing a userspace process. The host survived.

`.wslconfig` was then edited to give WSL 48 GB and a 64 GB swap file. The distro was not
restarted cleanly, so this session still runs on 31.3 GiB. Treat the numbers in the table as
the budget until `free -h` says otherwise.

## Rules

1. **Two heavy jobs repo-wide, not per agent.** All worktrees and subagents share this VM.
   At most two jobs that each load datasets or train may run at once; keep their combined
   peak RSS under ~20 GiB. When a job is heavy and long (a training ticket), run it alone.
2. **Cap every heavy process.** Preferred: `systemd-run --user --scope -p MemoryMax=10G
   -p MemorySwapMax=8G -- <cmd>`; portable fallback: `ulimit -v 10485760` in the shell that
   launches it. A cap that trips is a signal to restructure the job, not to raise the cap.
3. **Pre-flight check.** `free -g` before launching; abort if `available` is below 10 GiB.
   After any suspicious exit, check `free -h` and `dmesg | tail` before relaunching.
4. **Bounded data access.** Never materialise an observation × 18,164 dense matrix. Open
   `.h5ad` with `backed="r"` or slice with `read_elem`, use `float32`, iterate per fold /
   dose / plate, and memmap large caches to `.npy`. Never `n_jobs=-1`; use `n_jobs=2`.
5. **Thread caps when parallel.** With two jobs running:
   `OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2`, and
   `torch.set_num_threads(2)`.
6. **Measure and record.** Run anything above ~1 GiB under `/usr/bin/time -v` (or rely on
   the systemd-run scope accounting) and put the peak RSS in the ticket record. A rising
   peak in a test suite or build gate is a defect.
7. **Tests are not free.** Run test files individually (`-p 'test_x.py'`) instead of the
   full suite alongside another heavy job; a second concurrent full suite is what the crash
   record shows is unsafe.

## After a WSL restart

```bash
wsl.exe --shutdown          # from Windows, not from inside the distro
free -h                     # expect ~48 GiB; if not, the .wslconfig did not apply
nvidia-smi --query-gpu=memory.total,memory.used --format=csv
```

If `free -h` still reports ~32 GiB after a clean restart, check `.wslconfig` for a UTF-8 BOM
or a misspelled section before trusting the 48 GB setting.

## Known heavy jobs (measured)

| job | measured peak | consequence |
|---|---|---|
| `06.1_eligibility` offline build | 306 MiB | cheap |
| `05_prnet_pretrain` 68-epoch GPU pretrain on the 67,018 × 18,164 view | **15,150 MiB host RSS** | a full pretrain must run alone; leave at least 16 GiB free before starting |
| `script/verify_all.sh` (27 offline verifiers, 52 GB hashed) | 9.7 GiB charged, `4min 33s` CPU | mostly page cache; run it when no heavy job is running, not concurrently |
| `06.8_baselines --all` on the synthetic fixture (post-fix) | 1.2-1.7 GiB | dominated by imports and the shared evaluator, not the data |
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
