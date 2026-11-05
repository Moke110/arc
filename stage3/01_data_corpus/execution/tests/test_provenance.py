"""Offline production CLI regression using a temporary, isolated source inventory."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import test_workflow as workflow_helpers

class FrozenLocalProvenance(unittest.TestCase):
    def test_consumed_local_evidence_must_match_author_log_hash_and_bytes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            script=root/'corpus.py'
            shutil.copy(workflow_helpers.SCRIPT,script)
            subprocess.run(['git','init','-q',tmp],check=True)
            subprocess.run(['git','-C',tmp,'add','corpus.py'],check=True)
            subprocess.run(['git','-C',tmp,'-c','user.name=Test','-c','user.email=test@example.invalid',
                            'commit','-qm','isolated public CLI fixture'],check=True)
            item=root/'pheno_data/TemporarySource'
            (item/'source').mkdir(parents=True)
            (item/'local').mkdir()
            original=b'original evidence'
            curation=b'curation version 1'
            sha=lambda content:hashlib.sha256(content).hexdigest()
            (item/'source/source.txt').write_bytes(original)
            local=item/'local/curation.json'
            local.write_bytes(curation)
            (item/'manifest.sha256').write_text(sha(original)+'  source.txt\n')
            log={'dataset_id':'temporary_source','local_additions':[{'path':'local/curation.json',
                 'sha256':sha(curation),'bytes':len(curation),'origin':'frozen v1'}]}
            log_path=item/'author_log.json'
            log_path.write_text(json.dumps(log))
            (root/'pheno_data/registry.json').write_text(json.dumps({'schema_version':1,'items':[
                {'item':'TemporarySource','dataset_id':'temporary_source','status':'active'}]}))
            helper=workflow_helpers.Workflow()
            helper.asset=local
            record=helper.record('a',1.1,None)
            record['source']='temporary_source'
            record['measurements'][0]['asset']='pheno_data/TemporarySource/local/curation.json'
            bundle={'schema_version':1,'sources':[{'source':'temporary_source','item':'pheno_data/TemporarySource',
                'status':'active','test_only':False,'assets':[
                    {'path':'pheno_data/TemporarySource/source/source.txt','sha256':sha(original)},
                    {'path':'pheno_data/TemporarySource/local/curation.json','sha256':sha(curation)}]}],
                'records':[record]}
            def build(name,success):
                path=root/'bundle.json'
                path.write_text(json.dumps(bundle))
                result=subprocess.run([sys.executable,str(script),'build','--input',str(path),
                    '--output',str(root/name)],capture_output=True,text=True)
                self.assertEqual(result.returncode==0,success,result.stdout+result.stderr)
            build('valid',True)
            changed=b'curation changed without author-log revision'
            local.write_bytes(changed)
            bundle['sources'][0]['assets'][1]['sha256']=sha(changed)
            build('changed',False)
            local.write_bytes(curation)
            bundle['sources'][0]['assets'][1]['sha256']=sha(curation)
            log['local_additions'][0]['bytes']=len(curation)+1
            log_path.write_text(json.dumps(log))
            build('wrong-bytes',False)
            log['local_additions'][0]['bytes']=len(curation)
            del log['local_additions'][0]['sha256']
            log_path.write_text(json.dumps(log))
            build('missing-declared-hash',False)

if __name__=='__main__':
    unittest.main()
