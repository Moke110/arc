"""Record the coordinator's ten-second GPU and RAM samples until signaled."""

import argparse
from datetime import datetime
import json
from pathlib import Path
import shutil
import signal
import subprocess
import time

INTERVAL_SECONDS = 10
MEMORY_FIELDS = {"MemTotal", "MemAvailable", "SwapTotal", "SwapFree"}
GPU_COMMAND = (
    "nvidia-smi", "--id=0", "--query-gpu=memory.total,memory.used,utilization.gpu",
    "--format=csv,noheader,nounits",
)
WINDOWS_MEMORY_COMMAND = (
    "powershell.exe", "-NoProfile", "-NonInteractive", "-Command",
    '$m=Get-CimInstance Win32_OperatingSystem; '
    'Write-Output ($m.TotalVisibleMemorySize.ToString()+","+$m.FreePhysicalMemory.ToString())',
)


def read_sensor(command):
    return subprocess.run(command, capture_output=True, text=True, timeout=4, check=True).stdout.strip()


def windows_memory():
    if shutil.which(WINDOWS_MEMORY_COMMAND[0]) is None:
        return None
    try:
        return read_sensor(WINDOWS_MEMORY_COMMAND)
    except (OSError, subprocess.SubprocessError) as error:
        return {"error": str(error)}


def sample(deadline):
    started = time.monotonic()
    memory = {}
    for line in Path("/proc/meminfo").read_text().splitlines():
        key, value = line.split(":", maxsplit=1)
        if key in MEMORY_FIELDS:
            memory[f"{key}_kib"] = int(value.split()[0])
    gpu = read_sensor(GPU_COMMAND)
    total, used, utilization = map(int, gpu.split(","))
    if not (0 <= used <= total and total > 0 and 0 <= utilization <= 100):
        raise ValueError(f"invalid GPU sample: {gpu}")
    return {
        "wall_time": datetime.now().astimezone().isoformat(),
        "monotonic_seconds": started,
        "schedule_lag_seconds": started - deadline,
        "wsl": memory,
        "windows_host_total_and_available_kib": windows_memory(),
        "gpu_total_used_mib_and_utilization_percent": gpu,
        "workspace_free_bytes": shutil.disk_usage(Path.cwd()).free,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="path to a new JSONL file")
    args = parser.parse_args()
    stopped = False

    def request_stop(_signum, _frame):
        nonlocal stopped
        stopped = True

    signal.signal(signal.SIGINT, request_stop)
    signal.signal(signal.SIGTERM, request_stop)
    try:
        with args.output.open("x", encoding="utf-8") as stream:
            deadline = time.monotonic()
            while not stopped:
                print(json.dumps(sample(deadline)), file=stream, flush=True)
                deadline += INTERVAL_SECONDS
                while deadline <= time.monotonic():
                    deadline += INTERVAL_SECONDS
                delay = deadline - time.monotonic()
                while not stopped and delay > 0:
                    time.sleep(min(1, delay))
                    delay = deadline - time.monotonic()
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        parser.exit(1, f"Resource sampling failed: {error}\n")


if __name__ == "__main__":
    main()
