"""The #165/#170 matrix runner: declared plan, resource parsing, fail-closed
publish and the no-training acceptance/replay path."""
import contextlib
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

STEP_ROOT = Path(__file__).resolve().parents[1].parent
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_matrix  # noqa: E402

BASES = Path("/home/chang/arc/stage2/01_bases/output")


class PlanTests(unittest.TestCase):
    def test_plan_declares_the_two_by_two_matrix_in_order(self):
        cells = run_matrix.declared_identities(BASES / "fold_n_seed")
        self.assertEqual([(cell["fold"], cell["seed"]) for cell in cells],
                         [(0, 17), (0, 42), (1, 17), (1, 42)])
        self.assertEqual([cell["fold_row"] for cell in cells], [0, 0, 1, 1])
        self.assertEqual([cell["seed_row"] for cell in cells], [0, 1, 0, 1])
        self.assertEqual(run_matrix.ENTRIES, ("prnet", "mean"))

    def test_role_limits_are_the_full_fold_cohorts(self):
        import pandas as pd
        cohorts = pd.read_csv(BASES / "fold_n_seed/cohorts.csv", keep_default_na=False)
        limits = run_matrix.role_limits(cohorts, 0)
        self.assertEqual(limits, {"train": 2508, "validation": 297, "test": 1401})
        self.assertEqual(sum(limits.values()), len(cohorts[cohorts.fold.eq(0)]))


class TimeParsingTests(unittest.TestCase):
    SAMPLE = """
\tCommand being timed: "python"
\tUser time (seconds): 12.34
\tSystem time (seconds): 1.23
\tElapsed (wall clock) time (m:ss): 1:02.50
\tMaximum resident set size (kbytes): 1048576
"""

    def test_parses_wall_time_and_peak_rss(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "run.time"
            path.write_text(self.SAMPLE)
            value = run_matrix.parse_time(path)
        self.assertAlmostEqual(value["elapsed_seconds"], 62.5, places=3)
        self.assertEqual(value["peak_rss_kib"], 1048576)
        self.assertAlmostEqual(value["user_seconds"], 12.34)

    def test_missing_resource_fields_fail_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "run.time"
            path.write_text("no measurement here")
            with self.assertRaisesRegex(run_matrix.MatrixError, "resident"):
                run_matrix.parse_time(path)


class FailClosedTests(unittest.TestCase):
    def test_only_a_complete_terminal_status_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp) / "run"
            run.mkdir()
            for status in ("failed", "timeout", "over_budget", "interrupted"):
                (run / "terminal.json").write_text(json.dumps({"status": status}))
                with self.assertRaisesRegex(run_matrix.MatrixError, "not complete"):
                    run_matrix.check_terminal(run)
            (run / "terminal.json").write_text(json.dumps({"status": "complete"}))
            self.assertEqual(run_matrix.check_terminal(run)["status"], "complete")

    def test_publish_tree_replaces_a_stable_path_only_from_a_candidate(self):
        with tempfile.TemporaryDirectory() as tmp:
            candidate = Path(tmp) / "candidate/inputs"
            candidate.mkdir(parents=True)
            (candidate / "manifest.json").write_text('{"ok": true}')
            stable = Path(tmp) / "output/inputs"
            stable.mkdir(parents=True)
            (stable / "old.txt").write_text("superseded")
            run_matrix.publish_tree(candidate, stable)
            self.assertEqual((stable / "manifest.json").read_text(), '{"ok": true}')
            self.assertFalse((stable / "old.txt").exists())
            self.assertFalse((Path(tmp) / "output/inputs.publishing").exists())
            self.assertFalse((Path(tmp) / "output/inputs.superseded").exists())

    def test_missing_candidate_product_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(run_matrix.MatrixError, "missing"):
                run_matrix.publish_tree(Path(tmp) / "absent", Path(tmp) / "output/inputs")


class FigureValueContractTests(unittest.TestCase):
    def values(self):
        sys.path.insert(0, str(run_matrix.repo_root() / "stage2/01_bases/output/api"))
        from bases import plots
        data = STEP_ROOT / "figure/data"
        return plots, plots.read_plot_values(data / "metrics_values.csv"), \
            plots.read_plot_values(data / "metrics_indirect_values.csv"), \
            plots.read_plot_values(data / "learning_curve_values.csv")

    def test_recorded_tables_follow_the_adr_layout(self):
        plots, metrics, indirect, curve = self.values()
        self.assertEqual(metrics["schema"], plots.METRICS_SCHEMA)
        self.assertEqual((metrics["layout"]["rows"], metrics["layout"]["columns"]), (1, 4))
        self.assertEqual([panel["panel"] for panel in metrics["panels"]],
                         list(plots.PANEL_ORDER))
        self.assertEqual([panel["kind"] for panel in metrics["panels"]], ["bar"] * 4)
        self.assertEqual(len(metrics["panels"][0]["entries"]), 2)
        self.assertEqual({entry["series"] for entry in metrics["panels"][0]["entries"]},
                         {"PRnet", "Mean"})
        self.assertEqual(indirect["schema"], plots.INDIRECT_SCHEMA)
        self.assertEqual((indirect["layout"]["rows"], indirect["layout"]["columns"]), (2, 4))
        self.assertEqual([panel["panel"] for panel in indirect["panels"]],
                         list(plots.INDIRECT_PANEL_ORDER))
        self.assertEqual([panel["kind"] for panel in indirect["panels"]],
                         ["curve"] * 4 + ["bar", "bar", "heatmap", "bar"])
        reliability = indirect["panels"][3]
        self.assertEqual(reliability["panel"], "Reliability@d")
        self.assertEqual(len(reliability["entries"]), 2)
        for entry in reliability["entries"]:
            self.assertIn("reference", entry)
            self.assertEqual(len(entry["d_star"]), 4)
        confusion = indirect["panels"][6]
        self.assertEqual(confusion["classes"], ("up", "not_significant", "down"))
        self.assertEqual([entry["series"] for entry in confusion["entries"]],
                         ["PRnet", "Mean"])
        self.assertTrue(all(len(entry["matrix"]) == 3 for entry in confusion["entries"]))
        self.assertEqual(curve["schema"], plots.LEARNING_CURVE_SCHEMA)
        self.assertEqual(len(curve["panels"]), 1)
        self.assertEqual(curve["panels"][0]["panel"], "PRnet")
        self.assertEqual(len(curve["panels"][0]["runs"]), 4)

    def test_three_figures_are_published_and_verified(self):
        self.assertEqual(list(run_matrix.FIGURE_FILES),
                         ["metrics", "metrics_indirect", "learning_curve"])
        metrics, indirect, curve = run_matrix.verify_figure_values(STEP_ROOT / "figure/data")
        self.assertEqual(len(metrics["panels"]), 4)
        self.assertEqual(len(indirect["panels"]), 8)
        self.assertEqual(indirect["layout"]["rows"], 2)
        self.assertEqual(indirect["layout"]["columns"], 4)
        self.assertEqual(curve["schema"], "stage2-tenv-learning-curve-values-v1")


def _fake_run(root, entry, fold, seed, *, status="complete", digest_ok=True):
    run = Path(root) / "runs" / entry / f"fold_{fold}" / f"seed_{seed}"
    run.mkdir(parents=True)
    (run / "prediction.npy").write_bytes(b"prediction-bytes")
    (run / "case_ids.json").write_text('["case-a", "case-b"]\n')
    prediction = run_matrix.sha256(run / "prediction.npy") if digest_ok else "0" * 64
    (run / "identity.json").write_text(json.dumps({
        "identity": {"model": entry},
        "prediction_sha256": prediction,
        "case_ids_sha256": run_matrix.sha256(run / "case_ids.json"),
    }))
    (run / "terminal.json").write_text(json.dumps({"status": status}))
    return run


class AcceptedProductTests(unittest.TestCase):
    def test_recorded_accepted_hashes_match_the_published_figures_and_values(self):
        self.assertEqual(sorted(run_matrix.ACCEPTED_DATA), sorted(run_matrix.FIGURE_FILES))
        self.assertEqual(sorted(run_matrix.ACCEPTED_FIGURES), sorted(run_matrix.FIGURE_FILES))
        for name, digest in run_matrix.ACCEPTED_DATA.items():
            path = STEP_ROOT / "figure/data" / run_matrix.FIGURE_FILES[name][0]
            self.assertEqual(run_matrix.sha256(path), digest)
        for name, digest in run_matrix.ACCEPTED_FIGURES.items():
            path = STEP_ROOT / "figure" / run_matrix.FIGURE_FILES[name][1]
            self.assertEqual(run_matrix.sha256(path), digest)

    def test_digest_guard_refuses_changed_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "metrics_values.csv"
            path.write_text("changed")
            with self.assertRaisesRegex(run_matrix.MatrixError, "accepted"):
                run_matrix.verify_accepted_digests(
                    {"metrics": path}, {"metrics": "0" * 64})

    def test_replay_comparison_refuses_a_changed_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            first, second = Path(tmp) / "first", Path(tmp) / "second"
            for root in (first, second):
                root.mkdir()
                (root / "per_run.csv").write_text("a,b\n1,2\n")
            run_matrix.compare_product_trees(first, second, ("per_run",))
            (second / "per_run.csv").write_text("a,b\n1,3\n")
            with self.assertRaisesRegex(run_matrix.MatrixError, "replay"):
                run_matrix.compare_product_trees(first, second, ("per_run",))


class AcceptCommandTests(unittest.TestCase):
    def test_accept_mode_runs_only_the_rescore_and_renderer(self):
        commands = run_matrix.accept_commands(
            bases=BASES, output_root=Path("/tmp/tenv-output"), python="python")
        scripts = [Path(command[1]).name for command in commands]
        self.assertEqual(scripts, ["evaluation_driver.py", "evaluation_driver.py",
                                   "render_figures.py"])
        joined = " ".join(" ".join(command) for command in commands)
        for forbidden in ("prnet_fit_predict.py", "prepare_inputs.py", "mean_reference.py"):
            self.assertNotIn(forbidden, joined)


class PublishedRunTests(unittest.TestCase):
    def test_absent_run_set_fails_closed_before_any_command(self):
        with tempfile.TemporaryDirectory() as tmp:
            original = run_matrix.heavy
            run_matrix.heavy = lambda *args, **kwargs: self.fail(
                "accept mode invoked a command before verifying the run set")
            try:
                with self.assertRaisesRegex(run_matrix.MatrixError, "missing"):
                    run_matrix.run_accept(
                        bases=BASES, output_root=Path(tmp), python="python",
                        figure_dir=Path(tmp) / "figure",
                        checkpoint_root=Path(tmp) / "checkpoints",
                        record_path=Path(tmp) / "acceptance.json")
            finally:
                run_matrix.heavy = original

    def test_verify_published_runs_requires_every_complete_run(self):
        cells = [{"fold": 0, "seed": 17}]
        with tempfile.TemporaryDirectory() as tmp:
            _fake_run(tmp, "prnet", 0, 17)
            with self.assertRaisesRegex(run_matrix.MatrixError, "missing"):
                run_matrix.verify_published_runs(Path(tmp), cells)
            _fake_run(tmp, "mean", 0, 17, status="failed")
            with self.assertRaisesRegex(run_matrix.MatrixError, "not complete"):
                run_matrix.verify_published_runs(Path(tmp), cells)

    def test_verify_published_runs_binds_prediction_digest(self):
        cells = [{"fold": 0, "seed": 17}]
        with tempfile.TemporaryDirectory() as tmp:
            _fake_run(tmp, "prnet", 0, 17)
            _fake_run(tmp, "mean", 0, 17, digest_ok=False)
            with self.assertRaisesRegex(run_matrix.MatrixError, "digest"):
                run_matrix.verify_published_runs(Path(tmp), cells)


class SupersededPointerTests(unittest.TestCase):
    def test_pointer_is_added_without_reordering_the_historical_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "165-acceptance.json"
            path.write_text('{\n  "z_first": 1,\n  "a_second": 2\n}\n')
            pointer = {"record": "log/166-acceptance.json",
                       "schema": "stage2-tenv-acceptance-v2", "ticket": 170}
            run_matrix.mark_superseded(path, pointer)
            text = path.read_text()
            self.assertLess(text.index('"z_first"'), text.index('"a_second"'))
            self.assertEqual(json.loads(text)["superseded_by"], pointer)
            run_matrix.mark_superseded(path, pointer)
            self.assertEqual(path.read_text(), text)


class AcceptModeGuardTests(unittest.TestCase):
    def test_accept_only_refuses_training_plan_and_resume(self):
        for argv in (["--accept-only", "--all"],
                     ["--accept-only", "--plan-only"],
                     ["--accept-only", "--resume"]):
            with contextlib.redirect_stderr(io.StringIO()) as error:
                self.assertEqual(run_matrix.main(argv), 2)
            self.assertIn("accept-only", error.getvalue())

    def test_accept_only_refuses_an_unpublished_run_set(self):
        with tempfile.TemporaryDirectory() as tmp:
            with contextlib.redirect_stderr(io.StringIO()) as error:
                code = run_matrix.main(["--accept-only", "--output-root", tmp,
                                        "--record", str(Path(tmp) / "record.json")])
            self.assertEqual(code, 2)
            self.assertIn("missing", error.getvalue())
            self.assertFalse((Path(tmp) / "record.json").exists())


if __name__ == "__main__":
    unittest.main()
