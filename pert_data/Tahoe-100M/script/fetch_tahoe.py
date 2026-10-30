#!/usr/bin/env python3
"""Fetch ``tahoebio/Tahoe-100M`` into ``pert_data/Tahoe-100M/source/`` behind one scope switch.

Why this script exists
----------------------
Tahoe-100M is a 100M-cell drug-perturbation atlas served from HuggingFace as
3,388 expression parquet shards (337.6 GB) plus 1,026 pseudobulk-DE shards
(88.9 GB) plus the metadata tables.  The full snapshot is 428.8 GB; this host
has 287.2 GB free.  A *full* fetch therefore cannot be the default, but the
data home has to be the same one a full fetch would use, so that later work
never has to re-download or re-point anything.

The switch
----------
``--scope subset`` (default)
    8 expression shards, picked for plate diversity, plus the **complete**
    metadata tables.  ~3.1 GB.  Enough to exercise every data-shaped code path
    (all 62,710 genes, all 1,138 drug-dose labels, all 14 plates, all 50 cell
    lines, plate-matched DMSO controls) without moving the atlas.

``--scope full``
    all 3,388 expression shards + the complete metadata tables.  ~339.9 GB.
    (``--include-pseudobulk`` adds the 1,026 DE shards: 428.8 GB total.)

Both scopes write into the *same* directory and the *same* ledger, so switching
scope is an incremental top-up — files already present are never fetched twice.
The disk preflight refuses to start rather than dying at shard 2,900, because a
scope switch changes the file set, not the size of the disk.

Layout produced
---------------
::

    pert_data/Tahoe-100M/
      source/                             the download (upstream layout, bytes unchanged)
        data/train-NNNNN-of-03388.parquet expression shards (subset or full)
        metadata/*.parquet                the five metadata tables + vocabularies
        README.md  LICENSE.md  tutorials/ upstream docs (provenance)
        MANIFEST.jsonl                    path, bytes, sha256, source per file
        AXES.json                         the scope-independent category axes
        DATASET.json                      scope + counts + axis digests (the seam)
      local/_state/                       ledger, status, plate probes, HF cache
      local/Tahoe_summary.md              the local provenance note

``AXES.json``/``DATASET.json`` are what a training script should read: the axes
are derived from the metadata tables, so they describe the *whole* dataset even
when only a subset of shards is present.  A consumer can therefore tell whether
it is looking at ``scope: subset`` or ``scope: full`` without hard-coding
shard counts.

Interruption / resume
---------------------
* ``_state/ledger.jsonl`` — append-only, one line per completed file
  (``path``, ``bytes``, ``sha256``, ``revision``, ``completed_at``).
* ``_state/status.json`` — atomically rewritten progress snapshot (done/total,
  bytes, rate, ETA, failures) so an external watcher can read it safely.
* ``_state/lock`` — ``flock``; a second concurrent run exits with status 3.
* in-flight bytes live in ``hf_hub_download``'s ``.incomplete`` files under
  ``_state/hf`` and resume across runs; SIGINT stops new work, drains in-flight
  files, then exits cleanly leaving a resumable state.
* per-file retries: 3 attempts with 4/12/36 s backoff; a permanently failed
  file is recorded and the run continues (exit status 2).

Why the env vars are rewritten
------------------------------
Only the session workspace is writable on this host: ``/home/chang/.cache`` is
mounted read-only and ``hf_xet`` writes its chunk cache there by default, dying
with ``Read-only file system (os error 30)``.  ``HF_HOME``/``HF_XET_CACHE`` are
therefore redirected under ``_state/hf`` before ``huggingface_hub`` is imported.

Usage
-----
::

    # default: subset + full metadata, verify, write MANIFEST/AXES/DATASET
    python pert_data/Tahoe-100M/script/fetch_tahoe.py

    # show the plan and the disk arithmetic, touch nothing
    python pert_data/Tahoe-100M/script/fetch_tahoe.py --dry-run

    # top up to the full expression atlas (needs 287+ GB free)
    python pert_data/Tahoe-100M/script/fetch_tahoe.py --scope full

    # re-check an existing tree without touching the network for payload
    python pert_data/Tahoe-100M/script/fetch_tahoe.py --verify-only deep
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import io
import json
import os
import shutil
import signal
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path

def find_repo_root(start: Path) -> Path:
    """Walk up from ``start`` to the directory containing .git (ADR 0030)."""
    for candidate in [start, *start.parents]:
        if (candidate / ".git").exists():
            return candidate
    raise SystemExit(f"cannot find the repository root (.git) above {start}")


HERE = Path(__file__).resolve().parent
REPO_ROOT = find_repo_root(HERE)
ITEM_DIR = REPO_ROOT / "pert_data" / "Tahoe-100M"
SOURCE_DIR = ITEM_DIR / "source"      # the download: unmodified, upstream layout
LOCAL_DIR = ITEM_DIR / "local"        # local additions, declared in author_log.json


def state_dir(root: Path) -> Path:
    """Downloader state (ledger/status/caches) lives under the item's local/."""
    return LOCAL_DIR if root.resolve() == SOURCE_DIR.resolve() else root.parent / "local"

REPO_ID = "tahoebio/Tahoe-100M"
# Pinned so that a later `--scope full` fetches the exact revision the subset
# was validated against.  Update deliberately, never implicitly.
PINNED_REVISION = "349e2a55f5538a3e73b32a06c6bdd8e9e87011e8"

N_SHARDS = 3388
DEFAULT_SUBSET_SHARDS = 8
DEFAULT_WORKERS = 4          # measured: HF's CDN path saturates at ~4.5 MB/s
HEADROOM_BYTES = 2 * 1024**3  # refuse to start with less slack than this
# Metadata tables are loaded into RAM only if they are small.  obs_metadata is
# 2.29 GB compressed / 100.6M rows: loading it costs more RAM than this host has,
# and its row count is already checked from the parquet footer.
MAX_IN_MEMORY_META_BYTES = 512 * 1024**2
ASSUMED_BYTES_PER_SEC = 4.5e6

META_FILES = [
    "metadata/sample_metadata.parquet",
    "metadata/gene_metadata.parquet",
    "metadata/drug_metadata.parquet",
    "metadata/cell_line_metadata.parquet",
    "metadata/obs_metadata.parquet",
    "metadata/gene_vocabulary.json",
    "metadata/gene_vocabulary.jsonl",
]
DOC_FILES = [
    "README.md",
    "LICENSE.md",
    "tutorials/loading_data.ipynb",
    "tutorials/hackathon_setup.md",
]
PSEUDOBULK_PREFIX = "metadata/pseudobulk_differential_expression/"
# The upstream repo ships its own top-level README.md (the dataset card); our
# provenance note therefore follows the repo-wide dataset-summary convention
# (`DToxS_summary.md`, `L1000_summary.md`) instead of shadowing it.
SUMMARY_NAME = "Tahoe_summary.md"

# Frozen expectations, used as hard assertions by `--verify deep`.  These come
# from the dataset card and from the released `ST-HVG-Tahoe/var_dims.pkl`, which
# was cross-checked against these exact tables.
EXPECTED = {
    "metadata/gene_metadata.parquet": {"rows": 62710,
                                       "columns": ["gene_symbol", "ensembl_id", "token_id"]},
    "metadata/sample_metadata.parquet": {"rows": 1344,
                                         "columns": ["sample", "plate", "mean_gene_count",
                                                     "mean_tscp_count", "mean_mread_count",
                                                     "mean_pcnt_mito", "drug",
                                                     "drugname_drugconc"]},
    "metadata/drug_metadata.parquet": {"rows": 379},
    "metadata/cell_line_metadata.parquet": {"rows": 1000},
    "metadata/obs_metadata.parquet": {"rows": 100648790},
}
CONTROL_PERT = "[('DMSO_TF', 0.0, 'uM')]"
N_GENES, N_PERT, N_PLATES = 62710, 1138, 14


def root_dir(override: str | Path | None) -> Path:
    return Path(override).resolve() if override else SOURCE_DIR


# ---------------------------------------------------------------------------
# HF env must be redirected before huggingface_hub is imported.
# ---------------------------------------------------------------------------
def configure_hf_env(root: Path) -> None:
    hf_home = state_dir(root) / "hf"
    hf_home.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("HF_HOME", str(hf_home))
    os.environ.setdefault("HF_HUB_CACHE", str(hf_home / "hub"))
    os.environ.setdefault("HF_XET_CACHE", str(hf_home / "xet"))
    os.environ.setdefault("HF_HUB_DISABLE_TELEMETRY", "1")
    os.environ.setdefault("HF_HUB_DISABLE_PROGRESS_BARS", "1")
    os.environ.setdefault("HF_XET_HIGH_PERFORMANCE", "1")


# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------
def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def human(n: float) -> str:
    """1024-based sizes, labelled honestly (KiB/MiB/GiB) with the decimal value."""
    scaled = float(n)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if abs(scaled) < 1024 or unit == "TiB":
            if unit == "B":
                return f"{int(scaled)} B"
            return f"{scaled:.2f} {unit} ({n/1e9:.2f} GB)"
        scaled /= 1024
    return f"{scaled:.2f} TiB"


def sha256_of(path: Path, chunk: int = 8 * 1024**2) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def write_json_atomic(path: Path, payload) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    os.replace(tmp, path)


def append_jsonl(path: Path, record: dict) -> None:
    with path.open("a") as fh:
        fh.write(json.dumps(record, ensure_ascii=False) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def read_ledger(path: Path) -> dict[str, dict]:
    done: dict[str, dict] = {}
    if not path.exists():
        return done
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            rec = json.loads(line)
        except json.JSONDecodeError:
            continue  # tolerate a torn tail from an interrupted append
        if rec.get("status") == "ok":
            done[rec["path"]] = rec
    return done


# ---------------------------------------------------------------------------
# Shard selection: one plate per shard, so diversity == distinct plates.
# Probing a shard's plate costs two HTTP Range reads (~65 KB, ~2-3 s).
# ---------------------------------------------------------------------------
class _RangeFile(io.RawIOBase):
    """Minimal seekable HTTP-range reader, enough for a parquet footer/row group."""

    def __init__(self, url: str, size: int, session):
        super().__init__()
        self.url, self.size, self._pos, self._session = url, size, 0, session

    def readable(self):
        return True

    def seekable(self):
        return True

    def seek(self, off, whence=0):
        self._pos = off if whence == 0 else self._pos + off if whence == 1 else self.size + off
        return self._pos

    def tell(self):
        return self._pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self._pos
        if n <= 0 or self._pos >= self.size:
            return b""
        end = min(self._pos + n, self.size) - 1
        r = self._session.get(self.url, headers={"Range": f"bytes={self._pos}-{end}"},
                              timeout=120)
        r.raise_for_status()
        self._pos += len(r.content)
        return r.content

    def readinto(self, b):
        data = self.read(len(b))
        b[:len(data)] = data
        return len(data)


def shard_name(idx: int) -> str:
    return f"data/train-{idx:05d}-of-{N_SHARDS:05d}.parquet"


def probe_shard_plate(idx: int, revision: str, session, url_base: str) -> str | None:
    """Return the (single) plate that expression shard ``idx`` is built from."""
    import pyarrow.parquet as pq

    url = f"{url_base}/{shard_name(idx)}"
    head = session.head(url, allow_redirects=True, timeout=120)
    head.raise_for_status()
    size = int(head.headers["Content-Length"])
    rf = _RangeFile(head.url, size, session)
    pf = pq.ParquetFile(rf)
    table = pf.read_row_group(0, columns=["plate"])
    plates = sorted(set(table.column("plate").to_pylist()))
    return ",".join(plates) if plates else None


def probe_shard_dmso(idx: int, revision: str, session, url_base: str) -> dict:
    """Return ``{plate, dmso_cells, n_drugs}`` for expression shard ``idx``.

    Costs ~1 column read (~16-21 s, 29 row groups) — far more than the plate-only
    finger, so it is only run on shards the selection actually wants.  The
    control arm is mandatory: a shard of a plate that carries no ``DMSO_TF``
    cannot supply plate-matched vehicle controls for that plate.
    """
    import pyarrow.parquet as pq

    url = f"{url_base}/{shard_name(idx)}"
    head = session.head(url, allow_redirects=True, timeout=120)
    head.raise_for_status()
    size = int(head.headers["Content-Length"])
    pf = pq.ParquetFile(_RangeFile(head.url, size, session))
    plate = sorted(set(pf.read_row_group(0, columns=["plate"]).column("plate").to_pylist()))
    drugs = pf.read(columns=["drug"]).column("drug").to_pylist()
    counts: dict[str, int] = {}
    for name in drugs:
        if name is not None:
            counts[name] = counts.get(name, 0) + 1
    return {"plate": ",".join(plate) if plate else None,
            "dmso_cells": sum(v for k, v in counts.items() if "DMSO" in k),
            "n_drugs": len(counts)}


def _dmso_of(idx: int, revision: str, session, url_base: str,
             cache: dict) -> dict | None:
    key = str(idx)
    if key in cache and isinstance(cache[key], dict):
        return cache[key]
    try:
        record = probe_shard_dmso(idx, revision, session, url_base)
    except Exception as exc:
        print(f"  ! control probe failed for shard {idx}: {exc}", file=sys.stderr)
        return None
    cache[key] = record
    return record


def pick_shards(n: int, revision: str, cache_path: Path, dmso_cache_path: Path,
                use_probe: bool = True):
    """Choose ``n`` plate-diverse shards that also carry a vehicle-control arm.

    Step 1 (cheap, ~2 requests / 65 KB per probe) finds shards from distinct
    plates.  Step 2 (one column read each) checks that a picked shard actually
    contains ``DMSO_TF`` cells; shards without controls — about one in three —
    are topped up with a neighbouring shard from the same plate.
    """
    import requests

    grid = sorted({round(i * (N_SHARDS - 1) / max(1, 4 * n - 1)) for i in range(4 * n)})
    plate_cache: dict[str, str] = {}
    if cache_path.exists():
        try:
            plate_cache = json.loads(cache_path.read_text())
        except json.JSONDecodeError:
            plate_cache = {}
    dmso_cache: dict[str, dict] = {}
    if dmso_cache_path.exists():
        try:
            dmso_cache = json.loads(dmso_cache_path.read_text())
        except json.JSONDecodeError:
            dmso_cache = {}

    plate_of: dict[int, str] = {}
    session = requests.Session()
    session.headers.update({"User-Agent": "arc-tahoe-fetch/1.0"})
    url_base = f"https://huggingface.co/datasets/{REPO_ID}/resolve/{revision}"

    if use_probe:
        for idx in grid:
            key = str(idx)
            if key in plate_cache:
                plate_of[idx] = plate_cache[key]
                continue
            try:
                plate = probe_shard_plate(idx, revision, session, url_base)
            except Exception as exc:
                print(f"  ! plate probe failed for shard {idx}: {exc}", file=sys.stderr)
                plate = None
            if plate:
                plate_of[idx] = plate
                plate_cache[key] = plate
        write_json_atomic(cache_path, plate_cache)

    picked: list[int] = []
    seen_plates: set[str] = set()
    for idx in grid:
        plate = plate_of.get(idx)
        if plate and plate in seen_plates:
            continue
        picked.append(idx)
        if plate:
            seen_plates.add(plate)
        if len(picked) == n:
            break
    for idx in grid:
        if len(picked) == n:
            break
        if idx not in picked:
            picked.append(idx)

    rows: list[dict] = []
    strategy_parts = [f"probed {len(plate_of)} grid shards, {len(seen_plates)} distinct plates"]
    if not use_probe:
        rows = [{"shard": i, "plate": None, "dmso_cells": None, "role": "primary"}
                for i in picked[:n]]
        strategy_parts.append("no control probe (--no-probe): plate-matched DMSO NOT guaranteed")
        return sorted(r["shard"] for r in rows), rows, "; ".join(strategy_parts)

    topups = 0
    for idx in picked[:n]:
        record = _dmso_of(idx, revision, session, url_base, dmso_cache)
        if record is None:  # probe unusable: keep the shard, flag it
            rows.append({"shard": idx, "plate": None, "dmso_cells": None, "role": "primary"})
            continue
        if record["dmso_cells"] > 0:
            rows.append({"shard": idx, "plate": record["plate"],
                         "dmso_cells": record["dmso_cells"], "role": "primary"})
            continue
        # No vehicle control in this shard: the (plate, sublib) block spans ~1.5
        # shards, so a neighbour of the same plate almost always carries it.
        print(f"  · shard {idx} ({record['plate']}) has no DMSO_TF; looking for a "
              f"same-plate neighbour with controls", file=sys.stderr)
        rows.append({"shard": idx, "plate": record["plate"], "dmso_cells": 0,
                     "role": "primary-no-control"})
        for delta in (1, -1, 2, -2):
            cand = idx + delta
            if not 0 <= cand < N_SHARDS or cand in picked:
                continue
            nrec = _dmso_of(cand, revision, session, url_base, dmso_cache)
            if nrec and nrec["dmso_cells"] > 0 and nrec["plate"] == record["plate"]:
                rows.append({"shard": cand, "plate": nrec["plate"],
                             "dmso_cells": nrec["dmso_cells"], "role": "control-topup"})
                topups += 1
                print(f"  · + shard {cand} supplies {nrec['dmso_cells']} DMSO_TF "
                      f"cells for {nrec['plate']}", file=sys.stderr)
                break
    write_json_atomic(dmso_cache_path, dmso_cache)
    if topups:
        strategy_parts.append(f"+{topups} control-topup shard(s)")
    return (sorted(r["shard"] for r in rows), rows, "; ".join(strategy_parts))


# ---------------------------------------------------------------------------
# Planning
# ---------------------------------------------------------------------------
@dataclass
class FilePlan:
    path: str
    size: int
    sha256: str


def build_plan(api, revision: str, scope: str, n_shards: int,
               include_pseudobulk: bool, probe: bool, root: Path,
               reselect: bool = False):
    info = api.repo_info(REPO_ID, repo_type="dataset", revision=revision,
                         files_metadata=True)
    sizes = {s.rfilename: s for s in info.siblings}

    selection_path = state_dir(root) / "selection.json"
    if scope == "full":
        shards = list(range(N_SHARDS))
        strategy = "all expression shards"
    else:
        # Reuse a previously recorded selection: the fetched shard set must not
        # drift between runs (probing is best-effort and could otherwise change
        # which shards "subset" means).
        reused = None
        if selection_path.exists() and not reselect:
            try:
                record = json.loads(selection_path.read_text())
            except json.JSONDecodeError:
                record = None
            # schema 2 records the control arm; schema 1 selections had shards
            # without DMSO_TF and must not be reused.
            if (record and record.get("schema") == 2
                    and record.get("revision") == revision
                    and len(record.get("shards", [])) >= n_shards):
                reused = record
        if reused:
            shards = reused["shards"]
            strategy = (f"reused _state/selection.json (recorded "
                        f"{reused.get('recorded_at')}); {reused.get('strategy')}")
        else:
            shards, rows, how = pick_shards(
                n_shards, revision, state_dir(root) / "plate_probe.json",
                state_dir(root) / "shard_probe.json", use_probe=probe)
            strategy = f"{how}; shards={shards}"
            write_json_atomic(selection_path, {
                "schema": 2, "revision": revision, "shards": shards,
                "shards_detail": rows, "strategy": how, "recorded_at": now_iso()})

    wanted = [shard_name(i) for i in shards] + META_FILES + DOC_FILES
    if include_pseudobulk:
        wanted += sorted(p for p in sizes if p.startswith(PSEUDOBULK_PREFIX))
    # A shard list may name a file the revision does not contain (guard against
    # a wrong N_SHARDS for an older revision).
    missing = [p for p in wanted if p not in sizes]
    if missing:
        raise SystemExit(f"revision {revision} is missing {len(missing)} planned "
                         f"files, e.g. {missing[:3]}")

    plan = []
    for p in wanted:
        meta = sizes[p]
        sha = (meta.lfs.sha256 if getattr(meta, "lfs", None) else None)
        plan.append(FilePlan(path=p, size=int(meta.size or 0), sha256=sha or ""))
    return plan, strategy


def preflight(plan: list[FilePlan], root: Path) -> tuple[int, int, list[FilePlan]]:
    """Return (bytes_to_fetch, free_bytes, pending_files); refuse if it cannot fit."""
    root.mkdir(parents=True, exist_ok=True)
    pending, to_fetch = [], 0
    for item in plan:
        dest = root / item.path
        if dest.exists() and dest.stat().st_size == item.size:
            continue
        pending.append(item)
        # A stale partial file will be rewritten by hf_hub_download; count its
        # full size because the completed file replaces it.
        to_fetch += item.size
    free = shutil.disk_usage(root).free
    return to_fetch, free, pending


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------
class Progress:
    def __init__(self, total_files: int, total_bytes: int, status_path: Path):
        self.total_files, self.total_bytes = total_files, total_bytes
        self.done_files = self.failed = 0
        self.done_bytes = 0
        self.started = time.time()
        self.status_path = status_path
        self._lock = threading.Lock()

    def record(self, size: int, path: str, ok: bool = True) -> None:
        with self._lock:
            self.done_files += 1
            if ok:
                self.done_bytes += size
            else:
                self.failed += 1
            self._write()

    def _write(self) -> None:
        elapsed = max(time.time() - self.started, 1e-6)
        rate = self.done_bytes / elapsed
        remaining = max(self.total_bytes - self.done_bytes, 0)
        write_json_atomic(self.status_path, {
            "updated_at": now_iso(),
            "done_files": self.done_files,
            "total_files": self.total_files,
            "failed_files": self.failed,
            "done_bytes": self.done_bytes,
            "total_bytes": self.total_bytes,
            "rate_bytes_per_sec": round(rate, 1),
            "eta_seconds": round(remaining / rate, 1) if rate > 0 else None,
            "elapsed_seconds": round(elapsed, 1),
        })


def download_all(plan: list[FilePlan], pending: list[FilePlan], root: Path,
                 revision: str, workers: int, ledger_path: Path,
                 status_path: Path, stop: threading.Event) -> list[dict]:
    from huggingface_hub import hf_hub_download

    done_now = read_ledger(ledger_path)
    progress = Progress(len(pending), sum(p.size for p in pending), status_path)
    progress._write()
    failures: list[dict] = []

    def fetch(item: FilePlan) -> tuple[FilePlan, str | None]:
        last_err = None
        for attempt, delay in enumerate((0, 4, 12, 36)):
            if delay:
                time.sleep(delay)
            if stop.is_set():
                return item, "interrupted"
            try:
                path = hf_hub_download(repo_id=REPO_ID, filename=item.path,
                                       repo_type="dataset", revision=revision,
                                       local_dir=str(root))
                return item, None if Path(path).exists() else "missing after download"
            except Exception as exc:  # network, TLS, xet reconstruction…
                last_err = f"{type(exc).__name__}: {exc}"
        return item, last_err

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {}
        for item in pending:
            if stop.is_set():
                break
            futures[pool.submit(fetch, item)] = item
        for fut in as_completed(futures):
            item, err = fut.result()
            if err == "interrupted":
                # Paused via SIGTERM: queued files were never attempted.
                # Do not pollute the ledger/status with fake failures — the
                # resume run recomputes pending from disk.
                continue
            if err:
                rec = {"path": item.path, "bytes": item.size, "status": "failed",
                       "error": err, "completed_at": now_iso()}
                failures.append(rec)
                append_jsonl(ledger_path, rec)
                progress.record(item.size, item.path, ok=False)
                print(f"  ✗ {item.path}: {err}", file=sys.stderr)
            else:
                local = root / item.path
                rec = {"path": item.path, "bytes": item.size,
                       "sha256": item.sha256 or None, "revision": revision,
                       "status": "ok", "completed_at": now_iso()}
                done_now[item.path] = rec
                append_jsonl(ledger_path, rec)
                progress.record(item.size, item.path, ok=True)
                print(f"  ✓ {item.path} ({human(item.size)})", flush=True)
    return failures


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------
def verify(root: Path, plan: list[FilePlan], mode: str) -> dict:
    """Check the tree against ``plan``.  ``fast`` = presence+size, ``deep`` also
    recomputes sha256, reads every parquet footer, and asserts the frozen
    metadata expectations."""
    report = {"mode": mode, "checked_at": now_iso(), "files": 0, "problems": []}
    problems: list[str] = report["problems"]

    for item in plan:
        dest = root / item.path
        report["files"] += 1
        if not dest.exists():
            problems.append(f"missing: {item.path}")
            continue
        size = dest.stat().st_size
        if size != item.size:
            problems.append(f"size mismatch: {item.path} {size} != {item.size}")
            continue
        if mode == "deep" and item.sha256:
            got = sha256_of(dest)
            if got != item.sha256:
                problems.append(f"sha256 mismatch: {item.path} {got[:12]}… != "
                                f"{item.sha256[:12]}…")

    if mode != "deep":
        return report

    import pyarrow.parquet as pq
    import pandas as pd

    parquet_rows: dict[str, int] = {}
    for item in plan:
        if not item.path.endswith(".parquet"):
            continue
        dest = root / item.path
        if not dest.exists():
            continue
        try:
            meta = pq.ParquetFile(dest).metadata
            rg_rows = sum(meta.row_group(i).num_rows for i in range(meta.num_row_groups))
            parquet_rows[item.path] = meta.num_rows
            if rg_rows != meta.num_rows:
                problems.append(f"row-group sum mismatch: {item.path} "
                                f"{rg_rows} != {meta.num_rows}")
        except Exception as exc:
            problems.append(f"unreadable parquet: {item.path} ({exc})")

    frames = {}
    for rel in META_FILES:
        dest = root / rel
        if not rel.endswith(".parquet") or not dest.exists():
            continue
        if dest.stat().st_size > MAX_IN_MEMORY_META_BYTES:
            continue  # row count already asserted from the footer
        frames[rel] = pd.read_parquet(dest)
    for rel, expect in EXPECTED.items():
        if "rows" in expect:
            got = parquet_rows.get(rel, len(frames[rel]) if rel in frames else None)
            if got is not None and got != expect["rows"]:
                problems.append(f"{rel}: rows {got:,} != {expect['rows']:,}")
        if rel in frames and "columns" in expect:
            if list(frames[rel].columns) != expect["columns"]:
                problems.append(f"{rel}: columns {list(frames[rel].columns)} "
                                f"!= {expect['columns']}")

    genes = cells = samples = None
    if "metadata/gene_metadata.parquet" in frames:
        genes = frames["metadata/gene_metadata.parquet"]
        if len(genes) != N_GENES:
            problems.append(f"gene axis {len(genes)} != {N_GENES}")
        if genes["gene_symbol"].iloc[0] != "TSPAN6":
            problems.append(f"gene order: first symbol {genes['gene_symbol'].iloc[0]!r} != 'TSPAN6'")
        if genes["token_id"].iloc[0] != 3:
            problems.append(f"token_id offset: {genes['token_id'].iloc[0]} != 3")
    if "metadata/sample_metadata.parquet" in frames:
        samples = frames["metadata/sample_metadata.parquet"]
        if samples["drugname_drugconc"].nunique() != N_PERT:
            problems.append(f"pert axis {samples['drugname_drugconc'].nunique()} != {N_PERT}")
        if samples["plate"].nunique() != N_PLATES:
            problems.append(f"plate axis {samples['plate'].nunique()} != {N_PLATES}")
        if CONTROL_PERT not in set(samples["drugname_drugconc"]):
            problems.append(f"control label {CONTROL_PERT!r} absent from the pert axis")

    report["plates_observed"] = sorted(samples["plate"].unique().tolist()) if samples is not None else None
    report["problems"] = problems
    report["ok"] = not problems
    return report


def subset_coverage(root: Path) -> dict:
    """What the *downloaded* shards actually cover (plate/drug/cell-line counts)."""
    import pyarrow.parquet as pq

    shards = sorted((root / "data").glob("train-*.parquet"))
    if not shards:
        return {}
    plates, drugs, cell_lines, rows = set(), set(), set(), 0
    dmso_per_plate: dict[str, int] = {}
    for path in shards:
        table = pq.read_table(path, columns=["plate", "drug", "cell_line_id"])
        rows += table.num_rows
        pl = table.column("plate").to_pylist()
        dr = table.column("drug").to_pylist()
        plates |= set(pl)
        drugs |= set(dr)
        cell_lines |= set(table.column("cell_line_id").to_pylist())
        for plate, drug in zip(pl, dr):
            if drug and "DMSO" in drug:
                dmso_per_plate[plate] = dmso_per_plate.get(plate, 0) + 1
    return {"shards": len(shards), "rows": rows,
            "plates": sorted(plates), "n_drugs": len(drugs),
            "n_cell_lines": len(cell_lines), "cell_lines": sorted(cell_lines),
            "dmso_cells_per_plate": dict(sorted(dmso_per_plate.items())),
            "plates_without_controls": sorted(p for p in plates if not dmso_per_plate.get(p))}


# ---------------------------------------------------------------------------
# Artifacts: MANIFEST.jsonl, AXES.json, DATASET.json, README.md
# ---------------------------------------------------------------------------
def write_artifacts(root: Path, plan: list[FilePlan], revision: str, scope: str,
                    strategy: str, coverage: dict, verify_report: dict) -> None:
    import pandas as pd

    present = [p for p in plan if (root / p.path).exists()]
    manifest_path = root / "MANIFEST.jsonl"
    lines = []
    for item in sorted(present, key=lambda p: p.path):
        lines.append(json.dumps({
            "path": item.path,
            "bytes": item.size,
            "sha256": item.sha256 or sha256_of(root / item.path),
            "source": f"https://huggingface.co/datasets/{REPO_ID}/resolve/{revision}/{item.path}",
        }, ensure_ascii=False))
    manifest_path.write_text("\n".join(lines) + "\n")

    axes = {}
    if (root / "metadata/gene_metadata.parquet").exists():
        genes = pd.read_parquet(root / "metadata/gene_metadata.parquet")
        axes["genes"] = {
            "count": len(genes),
            "symbols": genes["gene_symbol"].tolist(),
            "ensembl_ids": genes["ensembl_id"].tolist(),
            "note": "Tahoe gene order; a row's genes/expressions lists index into this",
        }
    if (root / "metadata/sample_metadata.parquet").exists():
        samples = pd.read_parquet(root / "metadata/sample_metadata.parquet")
        axes["perturbations"] = {
            "count": int(samples["drugname_drugconc"].nunique()),
            # sorted, because cell_load builds its one-hot map as sorted(set(...))
            "labels": sorted(samples["drugname_drugconc"].unique().tolist()),
            "note": "drugname_drugconc strings, kept byte-exact (some carry a "
                    "trailing space); the control label is CONTROL_PERT.",
        }
        axes["plates"] = {"count": int(samples["plate"].nunique()),
                          "labels": sorted(samples["plate"].unique().tolist())}
        axes["drugs"] = {"count": int(samples["drug"].nunique()),
                         "labels": sorted(samples["drug"].unique().tolist())}
    if (root / "metadata/cell_line_metadata.parquet").exists():
        lines_df = pd.read_parquet(root / "metadata/cell_line_metadata.parquet")
        axes["cell_lines_catalog"] = {
            "count": int(lines_df["cell_name"].nunique()),
            "cellosaurus_ids": sorted(lines_df["Cell_ID_Cellosaur"].unique().tolist()),
            "note": "the catalog covers every line; the 50 realised in the atlas "
                    "are the subset with expression data",
        }
    axes["control"] = CONTROL_PERT
    axes["revision"] = revision
    write_json_atomic(root / "AXES.json", axes)

    shards_present = sorted(p.path for p in present if p.path.startswith("data/"))
    dataset = {
        "name": "Tahoe-100M",
        "repo_id": REPO_ID,
        "revision": revision,
        "scope": scope,
        "generated_at": now_iso(),
        "subset_strategy": strategy,
        "shards_present": len(shards_present),
        "shards_in_full_dataset": N_SHARDS,
        "files_present": len(present),
        "bytes_present": sum(p.size for p in present),
        "complete_metadata": True,
        "pseudobulk_included": any(p.path.startswith(PSEUDOBULK_PREFIX) for p in present),
        "axes_file": "AXES.json",
        "axes_digest": {k: (v.get("count") if isinstance(v, dict) else v)
                        for k, v in axes.items() if isinstance(v, dict)},
        "subset_coverage": coverage,
        "verification": {"mode": verify_report.get("mode"),
                         "ok": verify_report.get("ok"),
                         "problems": verify_report.get("problems", [])[:20]},
        "switch_to_full": ("re-run with --scope full; the same directory and "
                           "ledger are reused, so only missing shards are fetched"),
        "license": "CC0-1.0",
    }
    write_json_atomic(root / "DATASET.json", dataset)

    readme = [
        "# Tahoe-100M (vendored subset)",
        "",
        f"Fetched by `pert_data/Tahoe-100M/script/fetch_tahoe.py`.",
        "",
        f"- source: <https://huggingface.co/datasets/{REPO_ID}>",
        f"- pinned revision: `{revision}`",
        f"- license: CC0-1.0",
        f"- scope: **{scope}** — {dataset['shards_present']} of {N_SHARDS} expression shards",
        f"- complete metadata tables: yes (5 tables + 2 gene vocabularies)",
        f"- pseudobulk DE included: {dataset['pseudobulk_included']}",
        f"- bytes on disk: {human(dataset['bytes_present'])}",
        f"- verified: {verify_report.get('mode')} — ok={verify_report.get('ok')}",
        (f"- subset coverage: {coverage['rows']:,} cells, "
         f"{len(coverage['plates'])} plates, {coverage['n_cell_lines']} cell lines, "
         f"{coverage['n_drugs']} drugs" if coverage else "- subset coverage: n/a"),
        (f"- plate-matched vehicle control: {coverage['dmso_cells_per_plate']}"
         if coverage else "- plate-matched vehicle control: n/a"),
        "",
        "## What is here",
        "",
        "```",
        "source/data/train-NNNNN-of-03388.parquet  expression shards",
        "source/metadata/                          sample/gene/drug/cell_line/obs tables",
        "source/README.md LICENSE.md tutorials/    upstream documentation",
        "source/MANIFEST.jsonl                     path, bytes, sha256 per file",
        "local/Tahoe_summary.md                    this note (the upstream README.md is the dataset card)",
        "source/AXES.json                          category axes (scope-independent)",
        "source/DATASET.json                       scope + counts + axis digests",
        "local/_state/                             ledger, status, fetch caches",
        "```",
        "",
        "## Reading it",
        "",
        "`genes`/`expressions` in each row are parallel lists; the first entry is a",
        "marker/CLS token to drop.  Index into `AXES.json:genes.symbols` for gene",
        "identity.  `DMSO_TF` is the vehicle control and must be paired with `plate`",
        "for plate-matched controls.",
        "",
        "## Switching to the full atlas",
        "",
        "```",
        "python pert_data/Tahoe-100M/script/fetch_tahoe.py --scope full",
        "```",
        "",
        "The ledger and directory are shared, so this fetches only the shards that",
        "are not present yet.  The full expression atlas is ~339.9 GB; the disk",
        "preflight refuses to start if it cannot fit.",
        "",
    ]
    LOCAL_DIR.mkdir(parents=True, exist_ok=True)
    (LOCAL_DIR / SUMMARY_NAME).write_text("\n".join(readme))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Fetch Tahoe-100M into pert_data/Tahoe-100M/source/ (subset|full, resumable).")
    parser.add_argument("--scope", choices=["subset", "full"], default="subset",
                        help="file set to fetch (default: subset)")
    parser.add_argument("--shards", type=int, default=DEFAULT_SUBSET_SHARDS,
                        help=f"expression shards in subset scope (default {DEFAULT_SUBSET_SHARDS})")
    parser.add_argument("--revision", default=PINNED_REVISION,
                        help="dataset revision; defaults to the pinned commit")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                        help=f"parallel downloads (default {DEFAULT_WORKERS})")
    parser.add_argument("--include-pseudobulk", action="store_true",
                        help="also fetch metadata/pseudobulk_differential_expression/ (+88.9 GB)")
    parser.add_argument("--verify", dest="verify_after", choices=["none", "fast", "deep"],
                        default="deep", help="verification after the fetch (default deep)")
    parser.add_argument("--verify-only", choices=["fast", "deep"], default=None,
                        help="skip downloading, only verify what is on disk")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the plan and the disk arithmetic, fetch nothing")
    parser.add_argument("--no-probe", action="store_true",
                        help="skip plate probing; pick subset shards from a plain grid")
    parser.add_argument("--reselect", action="store_true",
                        help="ignore the recorded subset selection and pick fresh shards")
    parser.add_argument("--root", default=None, help="override the data home")
    args = parser.parse_args(argv)

    root = root_dir(args.root)
    configure_hf_env(root)
    from huggingface_hub import HfApi

    print(f"data home   : {root}")
    print(f"revision    : {args.revision}")
    print(f"scope       : {args.scope}"
          + (f" ({args.shards} shards)" if args.scope == "subset" else "")
          + (" + pseudobulk DE" if args.include_pseudobulk else ""))

    api = HfApi()
    plan, strategy = build_plan(api, args.revision, args.scope, args.shards,
                                args.include_pseudobulk, not args.no_probe, root,
                                reselect=args.reselect)
    print(f"selection   : {strategy}")

    total_bytes = sum(p.size for p in plan)
    to_fetch, free, pending = preflight(plan, root)
    print(f"planned     : {len(plan)} files, {human(total_bytes)}")
    print(f"pending     : {len(pending)} files, {human(to_fetch)}")
    print(f"free disk   : {human(free)}")
    if to_fetch:
        eta_h = to_fetch / ASSUMED_BYTES_PER_SEC / 3600
        print(f"estimated   : {eta_h:.1f} h at {ASSUMED_BYTES_PER_SEC/1e6:.1f} MB/s")

    if args.verify_only:
        report = verify(root, plan, args.verify_only)
        print(json.dumps({k: v for k, v in report.items() if k != "problems"}, indent=2))
        for problem in report["problems"][:40]:
            print(f"  ! {problem}")
        coverage = subset_coverage(root) if args.scope == "subset" else {}
        if coverage:
            print(f"subset covers: {len(coverage['plates'])} plates, "
                  f"{coverage['n_cell_lines']} cell lines, {coverage['n_drugs']} drugs, "
                  f"{coverage['rows']:,} cells")
        write_artifacts(root, plan, args.revision, args.scope, strategy, coverage, report)
        return 0 if report["ok"] else 1

    fits = to_fetch + HEADROOM_BYTES <= free
    if not fits:
        shortfall = to_fetch + HEADROOM_BYTES - free
        print(f"\nWOULD REFUSE: needs {human(to_fetch)} + {human(HEADROOM_BYTES)} "
              f"headroom, has {human(free)} (short by {human(shortfall)}).",
              file=sys.stderr)
        print("Scopes change the file set, not the disk. Free space, attach storage, "
              "or use a smaller scope.", file=sys.stderr)

    if args.dry_run:
        print(f"\ndry run: nothing fetched ({'fits' if fits else 'WOULD NOT FIT'})")
        return 0 if fits else 4

    if not fits:
        return 4

    state = state_dir(root)
    state.mkdir(parents=True, exist_ok=True)
    lock_path = state / "lock"
    lock_fh = lock_path.open("w")
    try:
        fcntl.flock(lock_fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        print(f"another fetch is running ({lock_path})", file=sys.stderr)
        return 3

    if pending:
        stop = threading.Event()

        def _handler(signum, _frame):
            print(f"\nsignal {signum}: draining in-flight files, then stopping",
                  file=sys.stderr)
            stop.set()

        signal.signal(signal.SIGINT, _handler)
        signal.signal(signal.SIGTERM, _handler)
        started = time.time()
        failures = download_all(plan, pending, root, args.revision, args.workers,
                                state / "ledger.jsonl", state / "status.json", stop)
        elapsed = time.time() - started
        fetched = sum(p.size for p in pending) - sum(f["bytes"] for f in failures)
        print(f"\nfetched {len(pending) - len(failures)}/{len(pending)} files, "
              f"{human(fetched)} in {elapsed/60:.1f} min "
              f"({fetched/max(elapsed,1e-6)/1e6:.2f} MB/s)")
        if failures:
            print(f"{len(failures)} file(s) failed; re-run to retry them", file=sys.stderr)

    coverage = subset_coverage(root) if args.scope == "subset" else {}
    if coverage:
        print(f"subset covers: {len(coverage['plates'])} plates, "
              f"{coverage['n_cell_lines']} cell lines, {coverage['n_drugs']} drugs, "
              f"{coverage['rows']:,} cells")
        print(f"control arm  : DMSO_TF cells per plate "
              f"{coverage['dmso_cells_per_plate']}")
        if coverage["plates_without_controls"]:
            print(f"  ! plates WITHOUT plate-matched vehicle control: "
                  f"{coverage['plates_without_controls']}", file=sys.stderr)

    report = {"mode": "skipped", "ok": None, "problems": []}
    if coverage.get("plates_without_controls"):
        # The subset must be loadable as a control-bearing dataset: every plate
        # needs its own vehicle wells (the repo's per-plate normalization
        # contract, ADR 0016).
        report["problems"] = [f"plate without vehicle control: {p}"
                              for p in coverage["plates_without_controls"]]
        report["ok"] = False
    if args.verify_after != "none":
        report = verify(root, plan, args.verify_after)
        if coverage.get("plates_without_controls"):
            report["problems"] = list(report.get("problems", [])) + [
                f"plate without vehicle control: {p}"
                for p in coverage["plates_without_controls"]]
            report["ok"] = False
        print(f"\nverify({args.verify_after}): {report['files']} files, "
              f"{len(report['problems'])} problem(s), ok={report['ok']}")
        for problem in report["problems"][:40]:
            print(f"  ! {problem}")

    write_artifacts(root, plan, args.revision, args.scope, strategy, coverage, report)
    print(f"\nwrote {root/'MANIFEST.jsonl'}\nwrote {root/'AXES.json'}\n"
          f"wrote {root/'DATASET.json'}\nwrote {LOCAL_DIR/SUMMARY_NAME}")
    return 0 if report.get("ok", True) else 1


if __name__ == "__main__":
    raise SystemExit(main())
