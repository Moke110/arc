"""Verify formal CUDA phase measurements and ten-second coordinator coverage."""

from chem_runs import PHASES, axis_run
from matrix_measured import GIB, MEMORY_CAP, SWAP_CAP
from report_shape import gaps_summary, read_samples


def audit_resources(activity, coordinator, fresh, *, bracket=False):
    sample_path = activity / "resource-samples.jsonl"
    matrix = read_samples(sample_path)
    scope_files = sorted(activity.glob("scope-*.json"))
    if not scope_files or not any(row["phase"] == "train" for row in matrix):
        raise ValueError("formal acceptance requires bounded CUDA training resource evidence")
    scopes = [axis_run.read_json(path) for path in scope_files]
    if any((scope["memory_max_bytes"], scope["memory_swap_max_bytes"]) != (MEMORY_CAP, SWAP_CAP)
           for scope in scopes):
        raise ValueError("formal memory cgroup differs from the accepted gate")
    start = matrix[0]["monotonic_seconds"]
    latest = matrix[-1]
    completed = int(scope_files[-1].stem.removeprefix("scope-")) / 1e9
    end = latest["monotonic_seconds"] + max(0, completed - latest["unix_seconds"])
    coverage = read_samples(coordinator)
    before = [row for row in coverage if row["monotonic_seconds"] <= start]
    running = [row for row in coverage if start < row["monotonic_seconds"] < end]
    after = [row for row in coverage if row["monotonic_seconds"] >= end]
    if not before or not after:
        raise ValueError("coordinator readings must bracket the entire formal matrix")
    selected = [before[-1], *running, after[0]]
    if not bracket and coverage != selected:
        raise ValueError("sealed coordinator record must contain only the formal window")
    coverage_gaps = gaps_summary(selected)
    windows = [row["windows_host_total_and_available_kib"] for row in selected]
    if coverage_gaps["gaps_over_12_seconds"] or any(not isinstance(value, str) for value in windows):
        raise ValueError("ten-second GPU, Windows and WSL sampling is required for the formal matrix")
    available = min(row["wsl_available_kib"] for row in matrix)
    coordinator_available = min(row["wsl"]["MemAvailable_kib"] for row in selected)
    host_available = min(int(value.split(",")[1]) for value in windows)
    gpu_available = min(row["gpu_total_mib"] - row["gpu_used_mib"] for row in matrix)
    disk_available = min(row["disk_free_bytes"] for row in matrix)
    scope_peak = max(scope["memory_peak_bytes"] for scope in scopes)
    if (min(available, coordinator_available) < 12 * GIB // 1024 or host_available <= 0
            or gpu_available < 512 or disk_available < GIB or scope_peak > MEMORY_CAP):
        raise ValueError("formal matrix violated its GPU, RAM or disk limits")

    peaks, files = {}, list(scope_files)
    for arm, fold, seed in fresh:
        directory = activity / "runs" / arm / f"fold-{fold}" / f"seed-{seed}"
        for phase in ("register", *PHASES):
            records = sorted(directory.glob(f"{phase}-*.resources"))
            if not records:
                raise ValueError(f"missing bounded {phase} measurement for {(arm, fold, seed)}")
            for record in records:
                output = record.with_suffix(".log")
                measure = next((line.rsplit(": ", 1)[1] for line in record.read_text().splitlines()
                                if "Maximum resident set size (kbytes):" in line), None)
                if measure is None or not output.is_file():
                    raise ValueError(f"incomplete phase resource measurement at {record}")
                files.extend((record, output))
                peaks[phase] = max(peaks.get(phase, 0), int(measure))
    return selected, {
        "status": "PASS", "matrix": gaps_summary(matrix), "coordinator": coverage_gaps,
        "scope_memory_max_bytes": MEMORY_CAP, "scope_swap_max_bytes": SWAP_CAP,
        "peak_scope_memory_bytes": scope_peak, "phase_peak_rss_kib": peaks,
        "minimum_wsl_available_kib": min(available, coordinator_available),
        "minimum_windows_available_kib": host_available, "minimum_gpu_available_mib": gpu_available,
        "maximum_gpu_used_mib": max(int(row["gpu_total_used_mib_and_utilization_percent"].split(",")[1])
                                    for row in selected),
        "minimum_disk_free_bytes": disk_available, "matrix_samples_sha256": axis_run.digest(sample_path),
        "phase_audit_sha256": axis_run.file_hashes(activity, [path.relative_to(activity).as_posix() for path in files]),
    }
