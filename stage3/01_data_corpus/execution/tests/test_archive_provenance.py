"""Source archive verifier rejects aggregate provenance drift offline."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'.git').exists())

class ArchiveProvenance(unittest.TestCase):
    def test_source_tree_summary_must_match_originals_and_manifest_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            (root/'.git').mkdir()
            (root/'script').mkdir()
            shutil.copy(ROOT/'script/phenotype_archive.py',root/'script/phenotype_archive.py')
            original=ROOT/'pheno_data/Liu-2021'
            item=root/'pheno_data/Liu-2021'
            shutil.copytree(original,item)
            log_path=item/'author_log.json'
            log=json.loads(log_path.read_text())
            expected={'files':2,'bytes':263362,'tree_sha256':hashlib.sha256((item/'manifest.sha256').read_bytes()).hexdigest()}
            log['source_tree']=expected.copy()
            log_path.write_text(json.dumps(log))
            def verify(success):
                result=subprocess.run([sys.executable,str(item/'script/download.py'),'--verify'],capture_output=True,text=True)
                self.assertEqual(result.returncode==0,success,result.stdout+result.stderr)
            verify(True)
            for field,value in [('files',3),('bytes',263363),('tree_sha256','0'*64)]:
                log['source_tree']=dict(expected,**{field:value})
                log_path.write_text(json.dumps(log))
                verify(False)
            del log['source_tree']
            log_path.write_text(json.dumps(log))
            verify(False)

if __name__=='__main__': unittest.main()
