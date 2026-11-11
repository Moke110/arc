"""External process and CLI encoding contracts without launching an R process."""
from contextlib import redirect_stdout
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bases.truth_successor import runtime_seal


class RuntimeProbeContract(unittest.TestCase):
    def test_r_argument_uses_one_escaped_newline_and_reads_two_versions(self):
        with tempfile.TemporaryDirectory() as temporary:
            executable = Path(temporary) / 'Rscript'
            executable.write_bytes(b'fixture executable, never launched')

            def process(arguments, **options):
                self.assertEqual(arguments[:3], [str(executable), '--vanilla', '-e'])
                self.assertIn(r'sep="\n"', arguments[3])
                self.assertNotIn(r'sep="\\n"', arguments[3])
                self.assertEqual(arguments[3].count('\\'), 1)
                self.assertEqual(options, {'capture_output': True, 'text': True, 'check': True})
                return SimpleNamespace(stdout='4.5.3\n1.50.2\n', stderr='', returncode=0)

            with patch('bases.truth_successor.subprocess.run', side_effect=process):
                report = runtime_seal(executable)
            self.assertEqual(report['r_version'], '4.5.3')
            self.assertEqual(report['deseq2_version'], '1.50.2')
            self.assertEqual(report['rscript_sha256'], hashlib.sha256(executable.read_bytes()).hexdigest())

    def test_cli_record_is_valid_json_with_a_real_trailing_newline(self):
        script = Path(__file__).resolve().parents[1] / 'bind_truth_successor.py'
        spec = importlib.util.spec_from_file_location('successor_cli_contract', script)
        cli = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cli)
        with tempfile.TemporaryDirectory() as temporary:
            record = Path(temporary) / 'report.json'
            arguments = ['bind_truth_successor.py', '--accepted-verification',
                         str(Path(temporary) / 'accepted.json'), '--record', str(record)]
            stdout = io.StringIO()
            with patch.object(sys, 'argv', arguments), \
                    patch.object(cli, 'create_successor', return_value={'result':'PASS'}), \
                    redirect_stdout(stdout):
                cli.main()
            self.assertEqual(json.loads(record.read_text()), {'result':'PASS'})
            self.assertEqual(record.read_bytes()[-1:], b'\n')
            self.assertEqual(json.loads(stdout.getvalue()), {'result':'PASS'})


if __name__ == '__main__':
    unittest.main()
