"""Documentation invariants for the sealed `tsg_view` report.

These tests read the real study documents (no synthetic fixture): the pre-registration must
be pinned verbatim and appear before any result, and the report's product links must resolve.
"""
import re
import unittest
from pathlib import Path

STEP = Path(__file__).resolve().parents[2]
README = STEP / "README.md"
PREREGISTRATION = STEP / "log" / "acceptance" / "pre_registration.md"

RESULTS_HEADING = "## Results"
PREREGISTRATION_HEADING = "## Pre-registration"
VERDICT_HEADING = "## Verdict against the pre-registration"


class PreRegistrationTest(unittest.TestCase):
    def setUp(self):
        self.readme = README.read_text()
        self.preregistration = PREREGISTRATION.read_text()

    def test_the_block_is_pinned_as_its_own_file(self):
        self.assertTrue(PREREGISTRATION.is_file())
        self.assertIn("tenv Δ PCC gain", self.preregistration)
        self.assertIn("81.5%", self.preregistration)
        # the file holds the block only: no heading, no result, no verdict
        self.assertNotIn("##", self.preregistration)
        self.assertNotIn("Not supported", self.preregistration)

    def test_the_file_is_contained_verbatim_in_the_readme(self):
        self.assertIn(self.preregistration.strip(), self.readme)

    def test_the_pre_registration_precedes_the_first_result(self):
        headings = [line for line in self.readme.splitlines() if line.startswith("## ")]
        self.assertIn(PREREGISTRATION_HEADING, headings)
        self.assertIn(RESULTS_HEADING, headings)
        self.assertLess(headings.index(PREREGISTRATION_HEADING),
                        headings.index(RESULTS_HEADING))
        block_at = self.readme.index(self.preregistration.strip())
        results_at = self.readme.index(RESULTS_HEADING)
        self.assertLess(block_at, results_at)

    def test_the_verdict_follows_the_results(self):
        self.assertGreater(self.readme.index(VERDICT_HEADING),
                           self.readme.index(RESULTS_HEADING))
        verdict = self.readme[self.readme.index(VERDICT_HEADING):]
        self.assertIn("**Not supported.**", verdict)


class ReportLinksTest(unittest.TestCase):
    def setUp(self):
        self.readme = README.read_text()

    def test_every_figure_and_report_link_resolves(self):
        targets = re.findall(r"\]\((figure/[^)]+|output/report/[^)]+)\)", self.readme)
        self.assertGreaterEqual(len(targets), 15)
        for target in targets:
            with self.subTest(target=target):
                self.assertTrue((STEP / target).is_file(), target)

    def test_the_report_publishes_no_weights(self):
        self.assertIn("no model weights", self.readme.lower().replace("**", ""))
        self.assertIn("No weights are published", self.readme)
        self.assertFalse((STEP / "output" / "model").exists())


class FigureCoverageTest(unittest.TestCase):
    def setUp(self):
        self.readme = README.read_text()

    def test_every_published_png_is_linked(self):
        pngs = sorted(path.name for path in (STEP / "figure").glob("*.png"))
        self.assertTrue(pngs, "the figure suite is empty")
        for name in pngs:
            with self.subTest(name=name):
                self.assertIn(f"figure/{name}", self.readme)


class AgentsHermeticityTest(unittest.TestCase):
    def setUp(self):
        self.agents = (STEP / "AGENTS.md").read_text()

    def test_no_unconditional_hermeticity_claim(self):
        self.assertNotIn("The tests are hermetic (synthetic fixtures, no payloads", self.agents)
        self.assertIn("two-tier", self.agents.lower())

    def test_the_two_tier_contract_names_the_override_variables(self):
        for name in ("TSG_NATIVE_SOURCE", "TSG_PAYLOAD_ROOT"):
            with self.subTest(name=name):
                self.assertIn(name, self.agents)
        self.assertIn("skipUnless", self.agents)
        self.assertIn("28 tests skip", self.agents)


if __name__ == "__main__":
    unittest.main()
