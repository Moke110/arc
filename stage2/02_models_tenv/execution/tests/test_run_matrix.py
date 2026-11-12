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
            self.assertTrue(run_matrix.compare_product_trees(first, second, ("per_run",)))
            (second / "per_run.csv").write_text("a,b\n1,3\n")
            with self.assertRaisesRegex(run_matrix.MatrixError, "replay"):
                run_matrix.compare_product_trees(first, second, ("per_run",))


class AcceptCommandTests(unittest.TestCase):
    def test_accept_mode_runs_only_the_two_rescore_commands_it_records(self):
        commands = run_matrix.accept_commands(
            bases=BASES, output_root=Path("/tmp/tenv-output"), python="python")
        scripts = [Path(command[1]).name for command in commands]
        self.assertEqual(scripts, ["evaluation_driver.py", "evaluation_driver.py"])
        joined = " ".join(" ".join(command) for command in commands)
        self.assertNotIn("render_figures.py", joined)
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


class RestoreTests(unittest.TestCase):
    SAMPLE = ("\tUser time (seconds): 1.00\n"
              "\tSystem time (seconds): 0.50\n"
              "\tElapsed (wall clock) time (m:ss): 0:02.00\n"
              "\tMaximum resident set size (kbytes): 2048\n")

    def test_restore_removes_products_that_did_not_preexist(self):
        with tempfile.TemporaryDirectory() as tmp:
            output, backup = Path(tmp) / "output", Path(tmp) / "backup"
            backup.mkdir()
            for name in ("evaluation", "figure_data"):
                (output / name).mkdir(parents=True)
                (output / name / "new.txt").write_text("first publish")
            run_matrix._restore_product(output, backup, ("evaluation", "figure_data"))
            self.assertFalse((output / "evaluation").exists())
            self.assertFalse((output / "figure_data").exists())

    def test_restore_puts_a_previous_publish_back(self):
        with tempfile.TemporaryDirectory() as tmp:
            output, backup = Path(tmp) / "output", Path(tmp) / "backup"
            (backup / "evaluation").mkdir(parents=True)
            (backup / "evaluation" / "old.txt").write_text("accepted")
            (output / "evaluation").mkdir(parents=True)
            (output / "evaluation" / "new.txt").write_text("failed publish")
            run_matrix._restore_product(output, backup, ("evaluation",))
            self.assertEqual((output / "evaluation/old.txt").read_text(), "accepted")
            self.assertFalse((output / "evaluation/new.txt").exists())

    def test_failed_second_command_removes_a_first_publish(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            output, figure = tmp / "output", tmp / "figure"
            calls = []

            def fake_heavy(command, *, stdout_path, time_path, **kwargs):
                calls.append(command)
                path = Path(time_path)
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(self.SAMPLE)
                if len(calls) == 1:
                    return run_matrix.parse_time(path)
                for name in ("evaluation", "figure_data"):
                    (output / name).mkdir(parents=True, exist_ok=True)
                    (output / name / "partial.txt").write_text("wrong bytes")
                raise run_matrix.MatrixError("driver failed after writing products")

            originals = {name: getattr(run_matrix, name) for name in (
                "heavy", "verify_published_runs", "verify_figure_values",
                "verify_accepted_digests")}
            run_matrix.heavy = fake_heavy
            run_matrix.verify_published_runs = lambda output_root, cells: []
            run_matrix.verify_figure_values = lambda path: None
            run_matrix.verify_accepted_digests = lambda paths, expected: True
            original_resources = run_matrix.RESOURCE_DIR
            run_matrix.RESOURCE_DIR = tmp / "resources"
            try:
                with self.assertRaisesRegex(run_matrix.MatrixError, "driver failed"):
                    run_matrix.run_accept(
                        bases=BASES, output_root=output, python="python",
                        figure_dir=figure, checkpoint_root=tmp / "checkpoints",
                        record_path=tmp / "record.json")
            finally:
                for name, value in originals.items():
                    setattr(run_matrix, name, value)
                run_matrix.RESOURCE_DIR = original_resources
            self.assertEqual(len(calls), 2)
            self.assertFalse((output / "evaluation").exists())
            self.assertFalse((output / "figure_data").exists())
            self.assertFalse((tmp / "record.json").exists())


class AllRecordTargetTests(unittest.TestCase):
    def test_all_writes_the_v2_record_to_166_and_only_marks_165(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            output, figure, candidate = tmp / "output", tmp / "figure", tmp / "candidate"
            for fold, seed in ((0, 17), (0, 42), (1, 17), (1, 42)):
                cell = {"fold": fold, "seed": seed}
                for entry in ("prnet", "mean"):
                    run = (output / "runs" / entry
                           / f"fold_{cell['fold']}" / f"seed_{cell['seed']}")
                    run.mkdir(parents=True)
                    (run / "resources.json").write_text(json.dumps(
                        {"resources": {"record": str(run / "r.time"),
                                       "peak_rss_kib": 1, "elapsed_seconds": 1.0}}))
            for root in (candidate / "evaluation", output / "evaluation"):
                root.mkdir(parents=True)
                (root / "evaluation.json").write_text(json.dumps({"schema": "x"}))
                for table in run_matrix.REPLAY_TABLES:
                    (root / f"{table}.csv").write_text(f"{table}\n")
            accepted_data = {}
            for name, (csv_name, _) in run_matrix.FIGURE_FILES.items():
                path = output / "figure_data" / csv_name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(name)
                accepted_data[name] = run_matrix.sha256(path)
            accepted_figures = {}
            for name, (_, png_name) in run_matrix.FIGURE_FILES.items():
                path = figure / png_name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(name)
                accepted_figures[name] = run_matrix.sha256(path)
            captured = {}

            def fake_builder(**kwargs):
                captured.update(kwargs)
                return {"schema": run_matrix.ACCEPT_SCHEMA, "ticket": 170}

            historical = tmp / "log/165-acceptance.json"
            historical.parent.mkdir(parents=True)
            historical.write_text('{\n  "z_first": 1,\n  "a_second": 2\n}\n')
            target = tmp / "log/166-acceptance.json"
            originals = {name: getattr(run_matrix, name) for name in (
                "STEP_ROOT", "ACCEPT_RECORD", "HISTORICAL_RECORD", "ACCEPTED_DATA",
                "ACCEPTED_FIGURES", "build_release", "train_prnet", "build_mean",
                "score", "verify_figure_values", "publish", "acceptance_record_v2",
                "environment_record", "verify_published_runs")}
            run_matrix.STEP_ROOT = tmp
            run_matrix.ACCEPT_RECORD = target
            run_matrix.HISTORICAL_RECORD = historical
            run_matrix.ACCEPTED_DATA = accepted_data
            run_matrix.ACCEPTED_FIGURES = accepted_figures
            run_matrix.build_release = lambda *a, **k: (None, {})
            run_matrix.train_prnet = lambda *a, **k: None
            run_matrix.build_mean = lambda *a, **k: None
            run_matrix.score = lambda *a, **k: None
            run_matrix.verify_figure_values = lambda *a, **k: None
            run_matrix.publish = lambda *a, **k: ({}, {}, {})
            run_matrix.environment_record = lambda python, commands: {
                "git": {"branch": "test", "head_sha": "0" * 40},
                "commands": []}
            run_matrix.acceptance_record_v2 = fake_builder
            run_matrix.verify_published_runs = lambda output_root, cells: []
            try:
                run_matrix.run_all(
                    bases=BASES, candidate=candidate, python="python",
                    prepare_python="python", train_python="python",
                    config=Path(tmp) / "config.json", figure_dir=figure,
                    resume=True)
            finally:
                for name, value in originals.items():
                    setattr(run_matrix, name, value)
            self.assertEqual(json.loads(target.read_text())["schema"],
                             run_matrix.ACCEPT_SCHEMA)
            self.assertTrue(captured["replay"]["tables_byte_identical"])
            self.assertTrue(captured["replay"]["evaluation_record_equal_except_paths"])
            self.assertTrue(captured["replay"]["figure_data_byte_identical_to_accepted"])
            self.assertTrue(captured["replay"]["figures_byte_identical_to_accepted"])
            historical_record = json.loads(historical.read_text())
            self.assertEqual(historical_record["z_first"], 1)
            self.assertEqual(historical_record["a_second"], 2)
            self.assertEqual(historical_record["superseded_by"]["record"],
                             "log/166-acceptance.json")
            text = historical.read_text()
            self.assertLess(text.index('"z_first"'), text.index('"a_second"'))


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
