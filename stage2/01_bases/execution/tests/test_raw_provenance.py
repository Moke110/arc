"""Raw release provenance resolves the producer's real manifest filenames."""
import hashlib
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import public_api  # noqa: F401  (binds the single public implementation in output/api)
from bases.raw import default_bindings, input_files, source_provenance


class RawProvenanceTest(unittest.TestCase):
    def test_production_defaults_bind_each_producer_manifest(self):
        with tempfile.TemporaryDirectory() as directory:
            bindings = default_bindings(Path(directory))
            for source, name in (("tenv", "tenv_manifest.json"), ("tahoe", "manifest.json")):
                raw = bindings[f"{source}.h5ad"]
                raw.parent.mkdir(parents=True, exist_ok=True)
                (raw.parent / name).write_text('{"source":"' + source + '"}\n')
            evidence = source_provenance(input_files(bindings), ("tenv", "tahoe"), fixture=False)
            tenv_manifest = bindings["tenv.h5ad"].parent / "tenv_manifest.json"
            self.assertEqual(evidence["tenv"]["release_manifest"], str(tenv_manifest.resolve()))
            self.assertEqual(evidence["tenv"]["release_manifest_sha256"],
                             hashlib.sha256(tenv_manifest.read_bytes()).hexdigest())
            tenv_manifest.unlink()
            with self.assertRaisesRegex(ValueError, "tenv_manifest.json"):
                source_provenance(input_files(bindings), ("tenv",), fixture=False)


if __name__ == "__main__":
    unittest.main()
