"""Input-byte validation, including deliberate damage to disposable copies."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class InputIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="figure1-input-test-")
        self.root = Path(self.temp.name)
        shutil.copytree(ROOT / "data", self.root / "data")
        (self.root / "scripts").mkdir()
        shutil.copy2(ROOT / "scripts/verify_inputs.py", self.root / "scripts/verify_inputs.py")

    def tearDown(self):
        self.temp.cleanup()

    def run_check(self):
        return subprocess.run(
            [sys.executable, str(self.root / "scripts/verify_inputs.py")],
            cwd=self.root.parent, capture_output=True, text=True)

    def test_unchanged_inputs_pass(self):
        result = self.run_check()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("7 input files", result.stdout)

    def test_modified_matrix_fails(self):
        matrix = self.root / "data/ani40_matrix.tsv"
        matrix.write_bytes(matrix.read_bytes() + b"\n")
        result = self.run_check()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("checksum mismatch", result.stderr)

    def test_missing_input_fails(self):
        (self.root / "data/tree_41tips_source.newick").unlink()
        self.assertNotEqual(self.run_check().returncode, 0)

    def test_manifest_cannot_escape_input_directory(self):
        path = self.root / "data/input_checksums.json"
        path.write_text(json.dumps({"../../outside.txt": "0" * 64}))
        result = self.run_check()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Invalid input path", result.stderr)


if __name__ == "__main__":
    unittest.main()
