#!/usr/bin/env python3
"""Mutation tests: every corruption must be REJECTED by the shared validator.

Each test copies the evidence to a temporary directory, mutates one thing, and
requires validate_evidence.validate() to report at least one problem. The
original package is never written to; one test asserts that explicitly.

    RBH_RESULTS=<path to results dir> python3 tests/test_mutations.py
"""
import copy, csv, hashlib, importlib.util, json, os, shutil, tempfile, unittest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
SRC = os.environ.get("RBH_RESULTS", os.path.join(PKG, "evidence"))

spec = importlib.util.spec_from_file_location(
    "ve", os.path.join(PKG, "bin", "validate_evidence.py"))
ve = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ve)


def dir_hash(path):
    h = hashlib.sha256()
    for root, _, files in os.walk(path):
        for f in sorted(files):
            p = os.path.join(root, f)
            h.update(os.path.relpath(p, path).encode())
            with open(p, "rb") as fh:
                h.update(fh.read())
    return h.hexdigest()


class MutationCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="rbh_mut_")
        self.res = os.path.join(self.tmp, "results")
        shutil.copytree(SRC, self.res)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def load(self, n):
        with open(os.path.join(self.res, n)) as fh:
            return json.load(fh)

    def save(self, n, obj):
        with open(os.path.join(self.res, n), "w") as fh:
            json.dump(obj, fh)

    def assertRejected(self, why):
        try:
            problems, _ = ve.validate(self.res, verbose=False)
        except ve.Failure:
            return
        self.assertTrue(problems, f"validator ACCEPTED corrupted evidence: {why}")

    def assertAccepted(self):
        problems, n = ve.validate(self.res, verbose=False)
        self.assertEqual(problems, [], "validator rejected the unmodified evidence")
        return n


class TestBaseline(MutationCase):
    def test_unmodified_evidence_is_accepted(self):
        self.assertEqual(self.assertAccepted(), 1774)


class TestWrongTargetIdentifier(MutationCase):
    def test_reviewers_exact_probe_is_rejected(self):
        """ALG19_illumina_contig_2_7 -> ..._999999999, raw hits left intact.

        This passed all 21 tests in revision 1. It must fail now.
        """
        calls = self.load("rbh_calls.json")
        hit = [c for c in calls
               if c["forward_hit"] == "ALG19_illumina_contig_2_7"]
        self.assertTrue(hit, "the probe target is not in this evidence")
        for c in calls:
            if c["forward_hit"] == "ALG19_illumina_contig_2_7":
                c["forward_hit"] = "ALG19_illumina_contig_2_999999999"
        self.save("rbh_calls.json", calls)
        self.assertRejected("called target renamed to a nonexistent protein")

    def test_arbitrary_wrong_target_is_rejected(self):
        calls = self.load("rbh_calls.json")
        for c in calls:
            if c["forward_hit"]:
                c["forward_hit"] = "NOT_A_REAL_PROTEIN"
                break
        self.save("rbh_calls.json", calls)
        self.assertRejected("a called target that is absent from the raw hits")


class TestRecordIntegrity(MutationCase):
    def test_missing_call_record_is_rejected(self):
        calls = self.load("rbh_calls.json")
        self.save("rbh_calls.json", calls[:-1])
        self.assertRejected("a call record was deleted")

    def test_duplicate_call_record_is_rejected(self):
        calls = self.load("rbh_calls.json")
        self.save("rbh_calls.json", calls + [copy.deepcopy(calls[0])])
        self.assertRejected("a call record was duplicated")

    def test_missing_genome_in_raw_evidence_is_rejected(self):
        fwd = self.load("forward_hits.json")
        del fwd[sorted(fwd)[0]]
        self.save("forward_hits.json", fwd)
        self.assertRejected("a genome was removed from the forward evidence")

    def test_truncated_json_is_rejected(self):
        p = os.path.join(self.res, "rbh_calls.json")
        with open(p) as fh:
            txt = fh.read()
        with open(p, "w") as fh:
            fh.write(txt[:len(txt)//2])
        self.assertRejected("rbh_calls.json was truncated")


class TestNumericIntegrity(MutationCase):
    def test_nonfinite_evalue_is_rejected(self):
        fwd = self.load("forward_hits.json")
        g = sorted(fwd)[0]; q = sorted(fwd[g])[0]
        if not fwd[g][q]:
            self.skipTest("no hits for this pair")
        fwd[g][q][0]["evalue"] = float("nan")
        with open(os.path.join(self.res, "forward_hits.json"), "w") as fh:
            fh.write(json.dumps(fwd))
        self.assertRejected("a raw E-value was set to NaN")

    def test_negative_evalue_is_rejected(self):
        fwd = self.load("forward_hits.json")
        g = sorted(fwd)[0]; q = sorted(fwd[g])[0]
        if not fwd[g][q]:
            self.skipTest("no hits for this pair")
        fwd[g][q][0]["evalue"] = -1.0
        self.save("forward_hits.json", fwd)
        self.assertRejected("a raw E-value was set negative")

    def test_threshold_violation_is_rejected(self):
        calls = self.load("rbh_calls.json")
        for c in calls:
            if c["reciprocal"]:
                c["forward_evalue"] = 1.0
                break
        self.save("rbh_calls.json", calls)
        self.assertRejected("a reciprocal call above the E-value threshold")


class TestDerivedOutputs(MutationCase):
    def test_flipped_matrix_cell_is_rejected(self):
        p = os.path.join(self.res, "presence_matrix.tsv")
        with open(p) as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        for r in rows[1:]:
            for i in range(1, len(r)):
                if r[i] == "1":
                    r[i] = "0"
                    break
            else:
                continue
            break
        with open(p, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerows(rows)
        self.assertRejected("one matrix cell was flipped")

    def test_nonreciprocal_call_marked_reciprocal_is_rejected(self):
        calls = self.load("rbh_calls.json")
        for c in calls:
            if not c["reciprocal"] and c["forward_hit"]:
                c["reciprocal"] = True
                break
        self.save("rbh_calls.json", calls)
        self.assertRejected("a non-reciprocal pair was marked reciprocal")

    def test_altered_rbh_total_is_rejected(self):
        calls = self.load("rbh_calls.json")
        n = 0
        for c in calls:
            if c["reciprocal"]:
                c["reciprocal"] = False; n += 1
                if n == 3: break
        self.save("rbh_calls.json", calls)
        self.assertRejected("three reciprocal calls were switched off")


class TestSyntheticTies(MutationCase):
    def test_synthetic_exact_tie_must_not_count_as_unique(self):
        """Duplicate the best forward hit under another target id, same E and
        score. The pair becomes ambiguous and must not remain a unique RBH."""
        fwd = self.load("forward_hits.json")
        calls = self.load("rbh_calls.json")
        target = next(c for c in calls if c["reciprocal"])
        g, q = target["genome"], target["query"]
        rows = fwd[g][q]
        clone = copy.deepcopy(rows[0]); clone["target"] = rows[0]["target"] + "_TIE"
        fwd[g][q] = [rows[0], clone] + rows[1:]
        self.save("forward_hits.json", fwd)
        self.assertRejected("an exact tie was introduced but the call still claims a unique best hit")


class TestReverseDirection(MutationCase):
    """Revision 2 compared the forward direction only. These are the review's
    probe 1 and its neighbours."""

    def test_reviewers_probe_negative_reverse_evalue(self):
        calls = self.load("rbh_calls.json")
        next(c for c in calls if c["reciprocal"])["reverse_evalue"] = -1.0
        self.save("rbh_calls.json", calls)
        self.assertRejected("reverse_evalue set to -1.0 with raw hits untouched")

    def test_altered_reverse_evalue_is_rejected(self):
        calls = self.load("rbh_calls.json")
        next(c for c in calls if c["reciprocal"])["reverse_evalue"] = 1e-9
        self.save("rbh_calls.json", calls)
        self.assertRejected("reverse_evalue silently changed")

    def test_altered_reverse_score_is_rejected(self):
        calls = self.load("rbh_calls.json")
        next(c for c in calls if c["reciprocal"])["reverse_score"] = 1.0
        self.save("rbh_calls.json", calls)
        self.assertRejected("reverse_score silently changed")

    def test_nonfinite_reverse_score_is_rejected(self):
        calls = self.load("rbh_calls.json")
        next(c for c in calls if c["reciprocal"])["reverse_score"] = float("inf")
        with open(os.path.join(self.res, "rbh_calls.json"), "w") as fh:
            fh.write(json.dumps(calls))
        self.assertRejected("reverse_score set to infinity")

    def test_flipped_reverse_ambiguity_is_rejected(self):
        calls = self.load("rbh_calls.json")
        next(c for c in calls if c["reciprocal"])["reverse_ambiguous"] = True
        self.save("rbh_calls.json", calls)
        self.assertRejected("reverse_ambiguous flipped without evidence")

    def test_invented_reverse_tie_alternatives_are_rejected(self):
        calls = self.load("rbh_calls.json")
        next(c for c in calls if c["reciprocal"])["reverse_tied_with"] = ["SOMETHING"]
        self.save("rbh_calls.json", calls)
        self.assertRejected("reverse_tied_with invented")

    def test_invented_forward_tie_alternatives_are_rejected(self):
        calls = self.load("rbh_calls.json")
        next(c for c in calls if c["reciprocal"])["forward_tied_with"] = ["SOMETHING"]
        self.save("rbh_calls.json", calls)
        self.assertRejected("forward_tied_with invented")

    def test_negative_bit_scores_in_raw_hits_are_allowed(self):
        """A negative bit score is legitimate and must not be rejected."""
        fwd = self.load("forward_hits.json")
        g = sorted(fwd)[0]
        q = next(q for q in sorted(fwd[g]) if len(fwd[g][q]) > 1)
        fwd[g][q][-1]["score"] = -3.5
        self.save("forward_hits.json", fwd)
        problems, _ = ve.validate(self.res, verbose=False)
        self.assertEqual([p for p in problems if "bit score" in p], [],
                         "a legitimate negative bit score was rejected")


class TestMatrixIdentifiers(MutationCase):
    """Review probe 2: shape preserved, identifiers corrupted."""

    def _universal_rows(self, rows):
        return [i for i, r in enumerate(rows[1:], 1) if all(v == "1" for v in r[1:])]

    def test_reviewers_probe_duplicated_universal_row(self):
        p = os.path.join(self.res, "presence_matrix.tsv")
        with open(p) as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        uni = self._universal_rows(rows)
        rows[uni[1]] = list(rows[uni[0]])
        with open(p, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerows(rows)
        self.assertRejected("one query row duplicated, another lost, shape unchanged")

    def test_duplicate_genome_column_is_rejected(self):
        p = os.path.join(self.res, "presence_matrix.tsv")
        with open(p) as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        rows[0][2] = rows[0][1]
        with open(p, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerows(rows)
        self.assertRejected("a genome column header duplicated")

    def test_renamed_query_row_is_rejected(self):
        p = os.path.join(self.res, "presence_matrix.tsv")
        with open(p) as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        rows[1][0] = "NOT_A_QUERY"
        with open(p, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerows(rows)
        self.assertRejected("a query row renamed to an unknown identifier")

    def test_short_row_is_rejected(self):
        p = os.path.join(self.res, "presence_matrix.tsv")
        with open(p) as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        rows[1] = rows[1][:-1]
        with open(p, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerows(rows)
        self.assertRejected("a matrix row is one field short")

    def test_non_binary_cell_is_rejected(self):
        p = os.path.join(self.res, "presence_matrix.tsv")
        with open(p) as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        rows[1][1] = "2"
        with open(p, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerows(rows)
        self.assertRejected("a matrix cell holds a non-binary value")


class TestNonDestructive(MutationCase):
    def test_source_evidence_is_never_written(self):
        before = dir_hash(SRC)
        calls = self.load("rbh_calls.json")
        calls[0]["forward_hit"] = "MUTATED"
        self.save("rbh_calls.json", calls)
        try:
            ve.validate(self.res, verbose=False)
        except ve.Failure:
            pass
        self.assertEqual(dir_hash(SRC), before,
                         "the original evidence directory changed during mutation testing")


if __name__ == "__main__":
    unittest.main(verbosity=2)
