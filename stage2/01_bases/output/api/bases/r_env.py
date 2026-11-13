"""Resolve the R interpreter and ``R_HOME`` from configuration, not host paths.

Two launchers consume this: :mod:`bases.cardiotox` runs ``Rscript`` (the
documented ``conda:biomni_e1`` build with edgeR/DESeq2) and
:mod:`bases.disease_deg` runs ``<R_HOME>/bin/exec/R`` (the ``arc`` env's R with
limma). Which R that is, is configuration rather than science, so the order is
declared once here: an explicit argument, then ``ARC_RSCRIPT`` / ``ARC_R_HOME``,
then the R that belongs to the conda env running the build, then the ``Rscript``
on ``PATH`` with ``R_HOME`` derived from it. Pinned manifests record only the
environment name and version, never a host path.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

RSCRIPT_VARIABLE = "ARC_RSCRIPT"
R_HOME_VARIABLE = "ARC_R_HOME"


class REnvironmentError(RuntimeError):
    """The R interpreter cannot be resolved from the current configuration."""


def rscript(explicit=None):
    """Return the ``Rscript`` command: explicit argument, ``ARC_RSCRIPT``, then ``PATH``."""
    override = explicit or os.environ.get(RSCRIPT_VARIABLE)
    if override:
        path = Path(override)
        if not path.is_file():
            raise REnvironmentError(
                f"Rscript not found: {path} (from {RSCRIPT_VARIABLE})")
        return str(path)
    found = shutil.which("Rscript")
    if not found:
        raise REnvironmentError(
            "no Rscript found on PATH; set ARC_RSCRIPT to the Rscript of an R "
            "environment that provides the required packages (the documented "
            "cardiotox build is conda:biomni_e1)")
    return found


def r_home(explicit=None, rscript_command=None):
    """Return ``R_HOME``: explicit, ``ARC_R_HOME``, the running env's R, then ``Rscript``."""
    if explicit:
        return Path(explicit)
    override = os.environ.get(R_HOME_VARIABLE)
    if override:
        return Path(override)
    bundled = _environment_r_home()
    if bundled is not None:
        return bundled
    command = rscript(rscript_command)
    derived = _derive_r_home(command)
    if derived is None:
        raise REnvironmentError(
            f"cannot derive R_HOME from {command}; set {R_HOME_VARIABLE}")
    return derived


def _environment_r_home():
    """The R of the conda env running this build (``<prefix>/lib/R``), if it ships one."""
    candidate = Path(sys.prefix) / "lib" / "R"
    if (candidate / "bin" / "exec" / "R").is_file():
        return candidate
    return None


def _derive_r_home(command):
    """Ask ``Rscript`` (then ``R``) for the home of the R it runs; ``None`` if neither answers."""
    answers = [[command, "-e", "cat(R.home())"]]
    r = shutil.which("R")
    if r:
        answers.append([r, "RHOME"])
    for answer in answers:
        try:
            completed = subprocess.run(answer, capture_output=True, text=True,
                                       check=False)
        except OSError:
            continue
        lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
        if completed.returncode == 0 and lines:
            return Path(lines[-1])
    return None


def entry_point(r_home_path=None, rscript_command=None):
    """Return the ``<R_HOME>/bin/exec/R`` entry point and its home."""
    home = r_home(r_home_path, rscript_command)
    executable = home / "bin" / "exec" / "R"
    if not executable.is_file():
        raise FileNotFoundError(f"pinned R entry point not found: {executable}")
    return executable, home
