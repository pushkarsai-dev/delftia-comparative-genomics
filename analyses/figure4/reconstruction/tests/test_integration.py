"""Cross-check the packaged reconstruction against its figure inputs."""
import csv
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
MODULE = BASE.parent


class IntegrationTests(unittest.TestCase):
    def test_reconstruction_matches_all_archived_gene_counts(self):
        with (BASE / "evidence/presence_matrix.tsv").open() as stream:
            rows = list(csv.reader(stream, delimiter="\t"))
        actual = {r[0]: sum(map(int, r[1:])) for r in rows[1:]}
        with (MODULE / "data/region002_gene_conservation.csv").open() as stream:
            archived = list(csv.DictReader(stream))
        expected = {r["gene"]: int(r["genomes_with_ortholog_of_40"]) for r in archived}
        self.assertEqual(actual, expected)

    def test_reconstruction_matches_all_archived_genome_counts(self):
        calls = json.loads((BASE / "evidence/rbh_calls.json").read_text())
        with (MODULE / "data/region002_rbh_matrix.csv").open() as stream:
            archived = list(csv.DictReader(stream))
        self.assertEqual(len(archived), 40)
        for row in archived:
            selected = [c for c in calls if c["genome"] == row["genome"] and c["reciprocal"]]
            self.assertEqual(len(selected), int(row["rbh_of_45"]))
            self.assertEqual(len({c["contig"] for c in selected}), int(row["n_contigs"]))

    def test_six_repeat_records_match_primary_raw_evidence(self):
        evidence = BASE / "evidence"
        calls = json.loads((evidence / "rbh_calls.json").read_text())
        forward = json.loads((evidence / "forward_hits.json").read_text())
        reverse = json.loads((evidence / "reverse_hits.json").read_text())
        manifest = json.loads((BASE / "repeat_evidence/bounded_repeat.json").read_text())
        self.assertEqual(len(manifest["genomes"]), 6)
        for entry in manifest["per_genome"]:
            genome = entry["genome"]
            data = (BASE / "repeat_evidence" / (genome + ".repeat.json")).read_bytes()
            self.assertEqual(hashlib.sha256(data).hexdigest(), entry["repeat_sha256"])
            obj = json.loads(data)
            self.assertEqual(obj["calls"], [c for c in calls if c["genome"] == genome])
            self.assertEqual(obj["forward"], forward[genome])
            self.assertEqual(obj["reverse"], reverse[genome])

    def test_prepare_inputs_checks_bytes_and_refuses_reuse(self):
        spec = importlib.util.spec_from_file_location("prepare", BASE / "bin/prepare_inputs.py")
        prep = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(prep)
        with tempfile.TemporaryDirectory(prefix="figure4-input-test-") as tmp:
            dest = Path(tmp) / "inputs"
            self.assertEqual(prep.prepare(dest, BASE), 41)
            with self.assertRaisesRegex(ValueError, "already exists"):
                prep.prepare(dest, BASE)


if __name__ == "__main__":
    unittest.main()
