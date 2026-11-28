"""Compare the model's normal loss call against a worked example."""
import unittest

import torch

from module_objective import configure_module_objective
from util.training.measured_roster import CONTROL, CONTROL_ONLY, DIRECT_DELTA, MEASURED_SEPARATE, MEASURED_SHARED


class ModuleObjectiveTest(unittest.TestCase):
    def test_every_readout_has_the_registered_full_loss(self):
        tensor = lambda value: torch.tensor([[value]], dtype=torch.float64)
        for arm in (CONTROL, DIRECT_DELTA, CONTROL_ONLY, MEASURED_SHARED, MEASURED_SEPARATE):
            with self.subTest(arm=arm):
                class Model:
                    def reconstruct_measured_delta(self, delta):
                        torch.testing.assert_close(delta, tensor(3))
                        return tensor(1), tensor(2), tensor(0)

                model = Model()
                configure_module_objective(model, arm)
                treated_decoder = arm in (CONTROL, DIRECT_DELTA)
                control_only = arm == CONTROL_ONLY
                terms = model.loss(
                    tensor(4), tensor(3), tensor(1), tensor(0), tensor(7),
                    tensor(5) if treated_decoder else None,
                    tensor(2) if treated_decoder else None,
                    tensor(0) if treated_decoder else None,
                    tensor(6), tensor(3), tensor(0),
                )
                expected = (2.05, 1, 0, 1, .5, 0, 0, 0) if control_only else (
                    6.3, 1, 4, 1, .5, 2, .5, 0)
                torch.testing.assert_close(torch.stack(terms), torch.tensor(expected, dtype=torch.float64))
