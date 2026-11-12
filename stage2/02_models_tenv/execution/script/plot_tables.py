"""ADR 0064 plot-value tables for the tenv comparison figures (#169).

The drawing API in ``stage2/01_bases/output/api/bases/plots.py`` is the single
schema authority: this module binds the shared panel constants, computes the
direct and ``@d`` reliability plotted values, and serialises them in the shared
long-form :data:`PLOT_VALUE_COLUMNS` schema.  It never duplicates a panel
constant or recomputes a scientific metric.

The two figures are ``metrics.png`` (one row by four direct-metric panels) and
``metrics_indirect.png`` (two rows by four: the four ``@d`` curves over ``d'``
with each run's ``d*`` marker and the ``F1@d`` reference on the Reliability
panel, then ``Reliability@d*`` and ``d*`` bars, the macro three-class confusion
heat map at ``d*`` and the null/changed Reliability split).  No retired ``@k``
metric name appears here.
"""
import sys
from pathlib import Path


def _api_root():
    """Locate the public implementation without assuming a fixed worktree."""
    for candidate in Path(__file__).resolve().parents:
        root = candidate / "stage2/01_bases/output/api"
        if (root / "bases/plots.py").is_file():
            return root
    raise ImportError("shared plotting implementation not found")


API = _api_root()
if str(API) not in sys.path:
    sys.path.insert(0, str(API))

from bases import plots  # noqa: E402  (public implementation bound above)

METRICS_SCHEMA = plots.METRICS_SCHEMA
METRICS_INDIRECT_SCHEMA = plots.INDIRECT_SCHEMA
LEARNING_CURVE_SCHEMA = plots.LEARNING_CURVE_SCHEMA
METRICS_TITLE = plots.METRICS_TITLE
METRICS_INDIRECT_TITLE = plots.INDIRECT_TITLE
LEARNING_CURVE_TITLE = plots.LEARNING_CURVE_TITLE
PLOT_VALUE_COLUMNS = plots.PLOT_VALUE_COLUMNS
DIRECT_LAYOUT = plots.METRICS_LAYOUT
DIRECT_PANEL_ORDER = plots.PANEL_ORDER
INDIRECT_LAYOUT = plots.INDIRECT_LAYOUT
INDIRECT_PANEL_ORDER = plots.INDIRECT_PANEL_ORDER
CURVE_METRICS = plots.CURVE_METRICS
CONFUSION_CLASSES = plots.CONFUSION_CLASSES
LINE_SHADE = plots.palette.LINE_SHADE


def palette_entry_hue(index):
    """The entry hue assigned by the shared stage palette."""
    return plots.palette.entry_hue(index)


def metrics_values(entries, *, title=METRICS_TITLE):
    """Complete direct-metric plotted values (one row by four columns)."""
    return plots.metrics_values(entries, title=title)


def indirect_metrics_values(entries, *, title=METRICS_INDIRECT_TITLE):
    """Complete ``@d`` reliability plotted values (two rows by four columns)."""
    return plots.indirect_metrics_values(entries, title=title)


def learning_curve_values(models, *, title=LEARNING_CURVE_TITLE):
    """Complete learning-curve plotted values through the shared API."""
    return plots.learning_curve_values(models, title=title)


def plot_rows(values):
    """Flatten the plotted values into the recorded long-form CSV rows."""
    return plots.plot_rows(values)


def write_plot_values(values, path):
    """Write the complete plotted values as the recorded long-form CSV."""
    return plots.write_plot_values(values, path)
