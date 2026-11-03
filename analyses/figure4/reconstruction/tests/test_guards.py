#!/usr/bin/env python3
"""Guard tests for the hardened runner. Small fixtures only; no full searches.

    python3 tests/test_guards.py
"""
import copy, importlib.util, json, os, shutil, tempfile, unittest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location("rr", os.path.join(PKG, "bin", "run_rbh.py"))
rr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rr)

MAN = {"runner_revision": "2.0.0", "versions": {"pyhmmer": "0.12.2"},
       "search": {"reverse_db": "proteome", "e_threshold": 1e-5},
       "inputs": {"g1": "a" * 64, "g2": "b" * 64}}
FP = rr.manifest_fingerprint(MAN)
QSET = {"q1", "q2"}


def _hit(target, e=1e-20, score=90.0):
    return dict(target=target, evalue=e, score=score, bias=0.0,
                dom_evalue=e, dom_score=score, env_from=1, env_to=10,
                ali_from=1, ali_to=10, hmm_from=1, hmm_to=10)


def _call(genome, q, hit):
    """A minimal but fully consistent call record."""
    return dict(genome=genome, query=q, forward_hit=hit, forward_evalue=1e-20,
                forward_score=90.0, forward_ambiguous=False, forward_tied_with=[],
                reverse_hit=q, reverse_evalue=1e-20, reverse_score=90.0,
                reverse_ambiguous=False, reverse_tied_with=[], reciprocal=True,
                contig="c1", target_begin=1, target_end=30, target_strand=1,
                query_len=10, target_len=10, query_coverage=1.0,
                target_coverage=1.0, n_forward_hits=1)


# A small checkpoint that is genuinely self-consistent, so the deep validation
# in revision 3 accepts it. Revision 2's fixture had empty raw evidence, which is
# exactly the hole the review found.
GOOD = {"calls": [_call("g1", "q1", "t1"), _call("g1", "q2", "t2")],
        "forward": {"q1": [_hit("t1")], "q2": [_hit("t2")]},
        "reverse": {"t1": [_hit("q1")], "t2": [_hit("q2")]},
        "fingerprint": FP, "n_forward_db": 10, "n_reverse_db": 20}


class TestFingerprint(unittest.TestCase):
    def test_fingerprint_changes_with_the_reverse_database(self):
        other = copy.deepcopy(MAN); other["search"]["reverse_db"] = "region45"
        self.assertNotEqual(rr.manifest_fingerprint(other), FP)

    def test_fingerprint_changes_with_an_input_hash(self):
        other = copy.deepcopy(MAN); other["inputs"]["g1"] = "c" * 64
        self.assertNotEqual(rr.manifest_fingerprint(other), FP)

    def test_fingerprint_changes_with_the_software_version(self):
        other = copy.deepcopy(MAN); other["versions"]["pyhmmer"] = "0.11.0"
        self.assertNotEqual(rr.manifest_fingerprint(other), FP)

    def test_fingerprint_is_stable_under_key_reordering(self):
        other = json.loads(json.dumps(MAN))
        self.assertEqual(rr.manifest_fingerprint(other), FP)


class TestCheckpointGuards(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="guard_")
        self.p = os.path.join(self.tmp, "g1.json")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def write(self, obj, raw=None):
        with open(self.p, "w") as fh:
            fh.write(raw if raw is not None else json.dumps(obj))

    def test_valid_checkpoint_is_accepted(self):
        self.write(GOOD)
        self.assertEqual(rr.load_checkpoint(self.p, FP, "g1", QSET)["n_forward_db"], 10)

    def test_foreign_fingerprint_is_refused(self):
        bad = copy.deepcopy(GOOD); bad["fingerprint"] = "z" * 64
        self.write(bad)
        with self.assertRaises(rr.Refused):
            rr.load_checkpoint(self.p, FP, "g1", QSET)

    def test_truncated_checkpoint_is_refused(self):
        self.write(None, raw=json.dumps(GOOD)[: len(json.dumps(GOOD)) // 2])
        with self.assertRaises(rr.Refused):
            rr.load_checkpoint(self.p, FP, "g1", QSET)

    def test_wrong_schema_is_refused(self):
        self.write({"calls": []})
        with self.assertRaises(rr.Refused):
            rr.load_checkpoint(self.p, FP, "g1", QSET)

    def test_wrong_query_set_is_refused(self):
        bad = copy.deepcopy(GOOD); bad["calls"] = [bad["calls"][0]]
        self.write(bad)
        with self.assertRaises(rr.Refused):
            rr.load_checkpoint(self.p, FP, "g1", QSET)

    def test_refusal_leaves_the_checkpoint_unchanged(self):
        bad = copy.deepcopy(GOOD); bad["fingerprint"] = "z" * 64
        self.write(bad)
        before = open(self.p, "rb").read()
        with self.assertRaises(rr.Refused):
            rr.load_checkpoint(self.p, FP, "g1", QSET)
        self.assertEqual(open(self.p, "rb").read(), before)


class TestDeepCheckpoint(unittest.TestCase):
    """Review probe 3: query NAMES match but the evidence is hollow."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="deep_")
        src = os.path.join(PKG, "repeat_evidence", "ALG19_illumina.hardened_runner.json")
        if not os.path.exists(src):
            self.fail("required packaged checkpoint evidence is missing")
        with open(src) as fh:
            self.good = json.load(fh)
        self.qs = {c["query"] for c in self.good["calls"]}
        self.p = os.path.join(self.tmp, "ALG19_illumina.json")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def write(self, obj):
        with open(self.p, "w") as fh:
            json.dump(obj, fh)

    def test_genuine_checkpoint_is_accepted(self):
        self.write(self.good)
        rr.load_checkpoint(self.p, self.good["fingerprint"], "ALG19_illumina", self.qs)

    def test_reviewers_probe_duplicate_call_and_empty_evidence(self):
        bad = json.loads(json.dumps(self.good))
        bad["calls"].append(dict(bad["calls"][0]))
        bad["forward"] = {}
        bad["reverse"] = {}
        self.write(bad)
        with self.assertRaises(rr.Refused):
            rr.load_checkpoint(self.p, bad["fingerprint"], "ALG19_illumina", self.qs)

    def test_duplicate_call_alone_is_refused(self):
        bad = json.loads(json.dumps(self.good))
        bad["calls"].append(dict(bad["calls"][0]))
        self.write(bad)
        with self.assertRaises(rr.Refused):
            rr.load_checkpoint(self.p, bad["fingerprint"], "ALG19_illumina", self.qs)

    def test_call_belonging_to_another_genome_is_refused(self):
        bad = json.loads(json.dumps(self.good))
        bad["calls"][0]["genome"] = "SOME_OTHER_GENOME"
        self.write(bad)
        with self.assertRaises(rr.Refused):
            rr.load_checkpoint(self.p, bad["fingerprint"], "ALG19_illumina", self.qs)

    def test_incomplete_forward_coverage_is_refused(self):
        bad = json.loads(json.dumps(self.good))
        bad["forward"].pop(sorted(bad["forward"])[0])
        self.write(bad)
        with self.assertRaises(rr.Refused):
            rr.load_checkpoint(self.p, bad["fingerprint"], "ALG19_illumina", self.qs)

    def test_invalid_database_size_is_refused(self):
        bad = json.loads(json.dumps(self.good))
        bad["n_forward_db"] = 0
        self.write(bad)
        with self.assertRaises(rr.Refused):
            rr.load_checkpoint(self.p, bad["fingerprint"], "ALG19_illumina", self.qs)

    def test_tampered_call_inconsistent_with_raw_hits_is_refused(self):
        bad = json.loads(json.dumps(self.good))
        for c in bad["calls"]:
            if c["forward_hit"]:
                c["forward_hit"] = "GHOST_PROTEIN"
                break
        self.write(bad)
        with self.assertRaises(rr.Refused):
            rr.load_checkpoint(self.p, bad["fingerprint"], "ALG19_illumina", self.qs)


class TestSettingsInFingerprint(unittest.TestCase):
    def test_effective_settings_are_part_of_the_manifest(self):
        obs, dec = rr.effective_settings()
        self.assertTrue(obs)
        man = dict(MAN)
        man["effective_settings_observed"] = obs
        man["effective_settings_declared"] = dec
        base = rr.manifest_fingerprint(man)
        for key, newval in (("builder.popen", 0.05), ("pipeline.F2", 0.01),
                            ("builder.score_matrix", "PAM30")):
            alt = json.loads(json.dumps(man))
            alt["effective_settings_observed"][key] = newval
            self.assertNotEqual(rr.manifest_fingerprint(alt), base,
                                f"changing {key} did not change the fingerprint")


class TestValidateBeforeWrite(unittest.TestCase):
    def test_populated_directory_without_a_manifest_is_refused_unchanged(self):
        tmp = tempfile.mkdtemp(prefix="vbw_")
        try:
            stray = os.path.join(tmp, "stray.txt")
            with open(stray, "w") as fh:
                fh.write("preexisting")
            before = sorted(os.listdir(tmp))
            with self.assertRaises(rr.Refused):
                rr.audit_output_dir(tmp, "f" * 64, {"g1"}, {"q1"})
            self.assertEqual(sorted(os.listdir(tmp)), before,
                             "a refused directory was written to")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_empty_directory_is_allowed(self):
        tmp = tempfile.mkdtemp(prefix="vbw2_")
        try:
            rr.audit_output_dir(tmp, "f" * 64, {"g1"}, {"q1"})
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_unknown_checkpoint_file_is_refused(self):
        tmp = tempfile.mkdtemp(prefix="vbw3_")
        try:
            man = {"fingerprint": "f" * 64, "manifest": MAN}
            man["fingerprint"] = rr.manifest_fingerprint(MAN)
            with open(os.path.join(tmp, "compatibility_manifest.json"), "w") as fh:
                json.dump(man, fh)
            ck = os.path.join(tmp, "per_genome")
            os.makedirs(ck)
            with open(os.path.join(ck, "NOT_A_PANEL_GENOME.json"), "w") as fh:
                json.dump({}, fh)
            with self.assertRaises(rr.Refused):
                rr.audit_output_dir(tmp, rr.manifest_fingerprint(MAN), {"g1", "g2"}, {"q1"})
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


class TestFastaGuards(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="fasta_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def write(self, text):
        p = os.path.join(self.tmp, "x.faa")
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def test_duplicate_fasta_identifier_is_refused(self):
        p = self.write(">a\nMKV\n>a\nMKW\n")
        with self.assertRaises(rr.Refused):
            list(rr.read_fasta(p))

    def test_unique_identifiers_are_accepted(self):
        p = self.write(">a\nMKV\n>b\nMKW\n")
        self.assertEqual(dict(rr.read_fasta(p)), {"a": "MKV", "b": "MKW"})


class TestAtomicWrite(unittest.TestCase):
    def test_no_temporary_file_survives(self):
        tmp = tempfile.mkdtemp(prefix="atomic_")
        try:
            p = os.path.join(tmp, "out.json")
            rr.atomic_write_json(p, {"a": 1})
            self.assertEqual(json.load(open(p)), {"a": 1})
            self.assertEqual([f for f in os.listdir(tmp) if f.endswith(".tmp")], [])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    def test_failed_write_leaves_the_original_intact(self):
        tmp = tempfile.mkdtemp(prefix="atomic2_")
        try:
            p = os.path.join(tmp, "out.json")
            rr.atomic_write_json(p, {"good": True})
            before = open(p, "rb").read()
            class Unserialisable:
                pass
            with self.assertRaises(TypeError):
                rr.atomic_write_json(p, {"bad": Unserialisable()})
            self.assertEqual(open(p, "rb").read(), before)
            self.assertEqual([f for f in os.listdir(tmp) if f.endswith(".tmp")], [])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


class TestSettingsReadback(unittest.TestCase):
    def test_observed_settings_come_from_the_installed_api(self):
        obs, dec = rr.effective_settings()
        self.assertEqual(obs["builder.seed"], 42)
        self.assertEqual(obs["builder.score_matrix"], "BLOSUM62")
        self.assertEqual(obs["background.L"], 350)
        self.assertIn("builder.symfrac", dec, "unreadable settings must be declared separately")
        self.assertIn("pipeline.Z", dec)


if __name__ == "__main__":
    unittest.main(verbosity=2)
