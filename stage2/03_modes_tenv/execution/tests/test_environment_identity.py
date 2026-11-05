"""Installed external dependencies must invalidate stale artifact CLI resumes."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from test_pipeline import fixture, SCRIPT


class EnvironmentIdentity(unittest.TestCase):
    def test_installed_package_change_rejects_resume_without_lock_changes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture(root/'input')
            package = root/'metadata/tenv_identity_fixture-1.0.dist-info'
            package.mkdir(parents=True)
            metadata = package/'METADATA'
            metadata.write_text('Metadata-Version: 2.1\nName: tenv-identity-fixture\nVersion: 1.0\n')
            adapter = root/'adapter.py'
            adapter.write_text("""import argparse,hashlib,json,numpy as np
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--request');r=json.loads(Path(p.parse_args().request).read_text());o=Path(r['output_directory'])
np.save(o/'prediction.npy',np.zeros((2,4),dtype='float32'))
(o/'case_ids.json').write_text(json.dumps(['test1','test2']))
(o/'model.bin').write_bytes(b'fixture')
(o/'checkpoint.json').write_text(json.dumps({'identity':r['identity'],'blobs':{'model.bin':hashlib.sha256(b'fixture').hexdigest()}}))
""")
            (root/'env.lock').write_text('intentionally unchanged recorded environment\n')
            (root/'adapter.identity.json').write_text(json.dumps({'source_files':['adapter.py'],
                'environment_files':['env.lock'],'config_files':['adapter.py'],
                'environment_probe_paths':['metadata']}))
            (root/'commands.json').write_text(json.dumps({'fixture':[sys.executable,str(adapter)]}))
            env = dict(os.environ)  # Only the declared interpreter overlay exposes fixture metadata.
            command = [sys.executable,str(SCRIPT),'predict','--input',str(root/'input'),
                       '--output',str(root/'output'),'--model','fixture','--commands',str(root/'commands.json')]
            first = subprocess.run(command,env=env,capture_output=True,text=True)
            self.assertEqual(first.returncode,0,first.stderr)
            unchanged = subprocess.run(command,env=env,capture_output=True,text=True)
            self.assertEqual(unchanged.returncode,0,unchanged.stderr)
            metadata.write_text(metadata.read_text().replace('Version: 1.0','Version: 2.0'))
            package.rename(package.with_name('tenv_identity_fixture-2.0.dist-info'))
            changed = subprocess.run(command,env=env,capture_output=True,text=True)
            self.assertNotEqual(changed.returncode,0,'Changed installed package was accepted with an unchanged lock file')
            self.assertIn('stale run identity',changed.stderr)

    def test_overlay_file_change_rejects_resume_with_unchanged_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture(root/'input')
            overlay = root/'overlay'; overlay.mkdir()
            dependency = root/'dependency'; dependency.mkdir()
            target = dependency/'native.so'; target.write_bytes(b'original')
            library = overlay/'native.so'; library.symlink_to(target)
            import hashlib
            (root/'overlay.sha256').write_text(hashlib.sha256(b'original').hexdigest()+'  native.so\n')
            adapter = root/'adapter.py'
            adapter.write_text("""import argparse,hashlib,json,numpy as np
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--request');r=json.loads(Path(p.parse_args().request).read_text());o=Path(r['output_directory'])
np.save(o/'prediction.npy',np.zeros((2,4),dtype='float32'))
(o/'case_ids.json').write_text(json.dumps(['test1','test2']))
(o/'model.bin').write_bytes(b'fixture')
(o/'checkpoint.json').write_text(json.dumps({'identity':r['identity'],'blobs':{'model.bin':hashlib.sha256(b'fixture').hexdigest()}}))
""")
            (root/'adapter.identity.json').write_text(json.dumps({'source_files':['adapter.py'],
                'environment_files':['overlay.sha256'],'config_files':['adapter.py'],
                'environment_file_manifests':[{'manifest':'overlay.sha256','root':'overlay',
                    'symlink_roots':['dependency']}]}))
            (root/'commands.json').write_text(json.dumps({'fixture':[sys.executable,str(adapter)]}))
            command = [sys.executable,str(SCRIPT),'predict','--input',str(root/'input'),
                       '--output',str(root/'output'),'--model','fixture','--commands',str(root/'commands.json')]
            first = subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(first.returncode,0,first.stderr)
            target.write_bytes(b'changed')
            changed = subprocess.run(command,capture_output=True,text=True)
            self.assertNotEqual(changed.returncode,0,'Changed overlay library accepted with unchanged manifest')
            self.assertIn('environment manifest file mismatch',changed.stderr)
            target.write_bytes(b'original')
            declaration = json.loads((root/'adapter.identity.json').read_text())
            declaration['environment_file_manifests'][0]['symlink_roots'] = []
            (root/'adapter.identity.json').write_text(json.dumps(declaration))
            denied = subprocess.run(command,capture_output=True,text=True)
            self.assertNotEqual(denied.returncode,0)
            self.assertIn('symlink escapes allowed roots',denied.stderr)
            declaration['environment_file_manifests'][0]['symlink_roots'] = ['dependency']
            (root/'adapter.identity.json').write_text(json.dumps(declaration))
            (root/'overlay.sha256').write_text(hashlib.sha256(b'original').hexdigest()+'  ../dependency/native.so\n')
            traversal = subprocess.run(command,capture_output=True,text=True)
            self.assertNotEqual(traversal.returncode,0)
            self.assertIn('manifest path escapes root',traversal.stderr)


    def test_failed_external_environment_probe_is_explicit(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture(root/'input')
            interpreter = root/'unavailable-python'
            interpreter.write_text('#!/bin/sh\nexit 1\n')
            interpreter.chmod(0o755)
            adapter = root/'adapter.py'; adapter.write_text('pass\n')
            (root/'adapter.identity.json').write_text(json.dumps({'source_files':['adapter.py'],
                'environment_files':['adapter.py'],'config_files':['adapter.py']}))
            (root/'commands.json').write_text(json.dumps({'fixture':[str(interpreter),str(adapter)]}))
            result = subprocess.run([sys.executable,str(SCRIPT),'predict','--input',str(root/'input'),
                '--output',str(root/'output'),'--model','fixture','--commands',str(root/'commands.json')],
                capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('external environment probe failed',result.stderr)


if __name__ == '__main__':
    unittest.main()
