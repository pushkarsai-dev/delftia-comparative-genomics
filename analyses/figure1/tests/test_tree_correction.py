#!/usr/bin/env python3
"""Regression tests for the distance-preserving Figure 1 tree correction.

Run with ``python -m unittest discover -s tests -p 'test_tree_correction.py'``.
All comparison inputs are bundled and read-only; no external paths or optional
environment variables are needed. The first 17 tests cover the original source,
correction, and historical-output checks. Additional tests exercise portability,
CLI reproducibility, and preservation of every bundled Newick input.
"""

import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
SOURCE = DATA / "tree_41tips_source.newick"
HISTORICAL = DATA / "tree_40tips_historical.newick"
CORRECTED = DATA / "tree_40tips_corrected.newick"
SCRIPT = ROOT / "scripts" / "rebuild_tree.py"
spec = importlib.util.spec_from_file_location("tree_rebuild", SCRIPT)
rb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rb)
RAW = SOURCE.read_text(encoding="utf-8")
SOURCE_TREE = rb.parse_newick(RAW)
FIXED = rb.prune_preserving(rb.parse_newick(RAW), rb.DROP_DEFAULT)
ACID3 = {"user14", "user39", "user5"}


class TestSource(unittest.TestCase):
    def test_source_is_41_tips(self):
        self.assertEqual(len(rb.leaves(SOURCE_TREE)), 41)

    def test_extra_tip_present(self):
        self.assertIn(rb.DROP_DEFAULT, [leaf.name for leaf in rb.leaves(SOURCE_TREE)])

    def test_shared_stem_is_not_sph1_only(self):
        nodes = [n for n in rb.walk(SOURCE_TREE)
                 if n.children and abs(n.length - 0.00643) < 1e-9]
        self.assertEqual(len(nodes), 1)
        self.assertEqual({leaf.name for leaf in rb.leaves(nodes[0])},
                         ACID3 | {rb.DROP_DEFAULT})
        self.assertEqual(nodes[0].name, "1.000")

    def test_sph1_root_distance_is_not_the_maximum(self):
        def rootdist(n, distance=0.0, out=None):
            out = {} if out is None else out
            distance += n.length
            if not n.children:
                out[n.name] = round(distance, 10)
            for child in n.children:
                rootdist(child, distance, out)
            return out

        distances = rootdist(SOURCE_TREE)
        self.assertAlmostEqual(distances[rb.DROP_DEFAULT], 0.01354, places=8)
        self.assertLess(distances[rb.DROP_DEFAULT], max(distances.values()))
        farthest = sorted(k for k, v in distances.items()
                          if abs(v - max(distances.values())) < 1e-12)
        self.assertEqual(farthest, ["user30", "user5"])


class TestCorrection(unittest.TestCase):
    def test_40_tips_out(self):
        self.assertEqual(len(rb.leaves(FIXED)), 40)
        self.assertNotIn(rb.DROP_DEFAULT, [leaf.name for leaf in rb.leaves(FIXED)])

    def test_all_780_distances_preserved(self):
        source, fixed = rb.pairwise(SOURCE_TREE), rb.pairwise(FIXED)
        self.assertEqual(len(fixed), 780)
        self.assertEqual([k for k in fixed if abs(fixed[k] - source[k]) > 1e-12], [])

    def test_unrooted_topology_unchanged(self):
        keep = {leaf.name for leaf in rb.leaves(FIXED)}

        def restrict(splits):
            out = set()
            for split in splits:
                first, second = tuple(split)
                a, b = first & keep, second & keep
                if a and b:
                    out.add(frozenset([a, b]))
            return out

        self.assertEqual(restrict(rb.splits(SOURCE_TREE)), rb.splits(FIXED))

    def test_root_not_moved(self):
        self.assertEqual([round(c.length, 8) for c in FIXED.children],
                         [round(c.length, 8) for c in SOURCE_TREE.children])
        biggest = max(SOURCE_TREE.children, key=lambda c: len(rb.leaves(c)))
        self.assertEqual(biggest.name, "0.769")

    def test_root_sides_rebalanced_without_moving_the_root(self):
        def maxd(node):
            def walk_distances(x, distance=0.0):
                distance += x.length
                if not x.children:
                    return [distance]
                return [v for child in x.children for v in walk_distances(child, distance)]
            return round(max(walk_distances(node)), 8)

        self.assertEqual(sorted(maxd(c) for c in FIXED.children), [0.01449, 0.01449])

    def test_missing_support_stays_missing(self):
        count = sum(1 for x in rb.walk(FIXED)
                    if x.children and x is not FIXED and not x.labelled)
        self.assertEqual(count, 10)

    def test_552_retained_and_1000_not_transferred(self):
        nodes = [x for x in rb.walk(FIXED)
                 if x.children and {leaf.name for leaf in rb.leaves(x)} == ACID3]
        self.assertEqual(len(nodes), 1)
        self.assertEqual(nodes[0].name, "0.552")
        self.assertAlmostEqual(nodes[0].length, 0.00698, places=8)

    def test_serializer_round_trips(self):
        again = rb.parse_newick(rb.write_newick(FIXED))
        self.assertEqual(rb.pairwise(again), rb.pairwise(FIXED))
        count = sum(1 for x in rb.walk(again)
                    if x.children and x is not again and not x.labelled)
        self.assertEqual(count, 10)

    def test_cli_refuses_a_nonempty_outdir(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory)
            sentinel = destination / "existing_result.txt"
            sentinel.write_bytes(b"existing result\n")
            result = subprocess.run(
                [sys.executable, str(SCRIPT), "--out-dir", str(destination)],
                capture_output=True, text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("not empty", result.stderr)
            self.assertEqual(sentinel.read_bytes(), b"existing result\n")
            self.assertEqual(list(destination.iterdir()), [sentinel])


class TestAgainstArchive(unittest.TestCase):
    def setUp(self):
        self.historical = rb.parse_newick(HISTORICAL.read_text(encoding="utf-8"))

    def test_archive_has_the_111_shortened_distances(self):
        source, archived = rb.pairwise(SOURCE_TREE), rb.pairwise(self.historical)
        differences = [round(source[k] - archived[k], 8) for k in archived]
        self.assertEqual(sorted(set(differences)), [0.0, 0.00643])
        self.assertEqual(differences.count(0.00643), 111)

    def test_every_shortened_pair_crosses_the_species_boundary(self):
        source, archived = rb.pairwise(SOURCE_TREE), rb.pairwise(self.historical)
        changed = [k for k in archived if abs(source[k] - archived[k]) > 1e-12]
        self.assertEqual(len(changed), 111)
        self.assertTrue(all((k[0] in ACID3) != (k[1] in ACID3) for k in changed))

    def test_archive_labels_ten_unsupported_clades_as_1(self):
        keep = {leaf.name for leaf in rb.leaves(self.historical)}
        source = {frozenset(leaf.name for leaf in rb.leaves(x)) & keep: x.labelled
                  for x in rb.walk(SOURCE_TREE) if x.children and x is not SOURCE_TREE}
        incorrect = [frozenset(leaf.name for leaf in rb.leaves(x))
                     for x in rb.walk(self.historical)
                     if x.children and x is not self.historical and x.name == "1"
                     and source.get(frozenset(leaf.name for leaf in rb.leaves(x))) is False]
        self.assertEqual(len(incorrect), 10)

    def test_archive_is_never_written(self):
        before = HISTORICAL.read_bytes()
        rb.prune_preserving(rb.parse_newick(RAW), rb.DROP_DEFAULT)
        self.assertEqual(HISTORICAL.read_bytes(), before)


class TestPortableCLI(unittest.TestCase):
    def test_rebuild_matches_accepted_corrected_tree_exactly(self):
        sources = [SOURCE, HISTORICAL, CORRECTED]
        before = {path: path.read_bytes() for path in sources}
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "rebuild"
            result = subprocess.run(
                [sys.executable, str(SCRIPT), "--out-dir", str(destination)],
                cwd=directory, capture_output=True, text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            output = (destination / "tree_40tips_corrected.newick").read_bytes()
            self.assertEqual(output, before[CORRECTED])
            report = json.loads((destination / "rebuild_report.json").read_text())
            self.assertEqual(report["tips_in"], 41)
            self.assertEqual(report["tips_out"], 40)
            self.assertEqual(report["retained_pairs"], 780)
            self.assertEqual(report["distances_equal_to_source"], 780)
            self.assertEqual(report["unlabelled_internal_nodes"], 10)
            self.assertEqual(report["source_file_sha256"],
                             hashlib.sha256(before[SOURCE]).hexdigest())
            self.assertEqual(report["output_file_sha256"],
                             hashlib.sha256(output).hexdigest())
        self.assertEqual({path: path.read_bytes() for path in sources}, before)

    def test_repository_input_directory_cannot_be_used_as_output(self):
        before = {path: path.read_bytes() for path in DATA.glob("*.newick")}
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--out-dir", str(DATA)],
            capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("not empty", result.stderr)
        self.assertEqual({path: path.read_bytes() for path in DATA.glob("*.newick")}, before)

    def test_missing_source_does_not_create_output(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "rebuild"
            missing = Path(directory) / "missing.newick"
            result = subprocess.run(
                [sys.executable, str(SCRIPT), "--source", str(missing),
                 "--out-dir", str(destination)],
                capture_output=True, text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(destination.exists())

    def test_missing_tip_does_not_create_output(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "rebuild"
            result = subprocess.run(
                [sys.executable, str(SCRIPT), "--drop", "absent-tip",
                 "--out-dir", str(destination)],
                capture_output=True, text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("not found", result.stderr)
            self.assertFalse(destination.exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
