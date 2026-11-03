"""Regression and rejection checks for the full-panel repeat comparator."""
import importlib.util
import json
import shutil
import tempfile
import unittest
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('repeat_comparator', BASE/'bin/compare_repeat.py')
comparator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparator)


class FullRepeat(unittest.TestCase):
    def test_all_calls_and_saved_comparison(self):
        actual = comparator.compare(BASE/'evidence', BASE/'full_panel_repeat/results')
        self.assertEqual(actual, json.loads((BASE/'full_panel_repeat/comparison.json').read_text()))

    def reject(self, mutate):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            for name in ('rbh_calls.json','presence_matrix.tsv','query_id_check.json','proteome_manifest.json'):
                shutil.copy2(BASE/'full_panel_repeat/results'/name, out/name)
            mutate(out)
            with self.assertRaises(ValueError):
                comparator.compare(BASE/'evidence', out)

    def test_changed_target_rejected(self):
        def mutate(p):
            f=p/'rbh_calls.json'; data=json.loads(f.read_text())
            data[0]['forward_hit']='invalid_target'; f.write_text(json.dumps(data))
        self.reject(mutate)

    def test_duplicate_call_rejected(self):
        def mutate(p):
            f=p/'rbh_calls.json'; data=json.loads(f.read_text())
            data[1]=data[0]; f.write_text(json.dumps(data))
        self.reject(mutate)

    def test_matrix_change_rejected(self):
        def mutate(p):
            f=p/'presence_matrix.tsv'; rows=f.read_text().splitlines()
            cells=rows[1].split('\t'); cells[1]=str(1-int(cells[1]))
            rows[1]='\t'.join(cells); f.write_text('\n'.join(rows)+'\n')
        self.reject(mutate)

    def test_protein_hash_change_rejected(self):
        def mutate(p):
            f=p/'proteome_manifest.json'; data=json.loads(f.read_text())
            data[next(iter(data))]['proteome_sha256']='0'*64
            f.write_text(json.dumps(data))
        self.reject(mutate)


if __name__=='__main__':
    unittest.main()
