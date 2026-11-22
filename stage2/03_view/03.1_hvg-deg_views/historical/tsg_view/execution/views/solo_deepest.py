"""The `solo-deepest` sensitivity view: one profile per condition, deepest retained total.

This is the same member-choice recipe as :mod:`views.solo` with the draw replaced by
``argmax(library_total)`` and a deterministic tie-break to the smallest ``profile_id``
(``views.solo_choice`` owns both rules). A single-member condition keeps its only member, so
the view equals ``solo`` exactly there; every row still reuses the chosen member's full-view
``X``/control, so the freeze sees the same splits and the same vehicle tail.
"""
from . import CONTROL_RULE, RECIPE
from .solo import SELECTION, materialize


def build(context):
    return materialize(context, rule="deepest", selection=SELECTION["deepest"],
                       seed=None, recipe=RECIPE["solo-deepest"],
                       control_rule=CONTROL_RULE["solo-deepest"])
