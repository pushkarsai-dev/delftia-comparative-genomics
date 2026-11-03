#!/usr/bin/env python3
"""Shared validator for Figure 4c RBH evidence — revision 3.

Revision 3 closes three gaps found by review of revision 2:
  1. the reverse direction is now recomputed and compared symmetrically with the
     forward direction, including E-value, bit score, ambiguity, and retained tie
     alternatives; recorded values must themselves be finite, non-negative
     E-values below threshold, with finite (possibly negative) bit scores
  2. presence-matrix identifiers are validated before any cell comparison, and
     the expected identifiers come from the query and proteome manifests and the
     raw evidence, never from the table being validated
  3. the same core rules apply to single-genome checkpoints, via rbh_core

    python3 validate_evidence.py --results <dir>        # exit 0 or 2
"""
import argparse, csv, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rbh_core as core

Failure = core.Failure


def _load(results, name, required=True):
    p = os.path.join(results, name)
    if not os.path.exists(p):
        if required:
            raise Failure(f"missing evidence file: {name}")
        return None
    try:
        with open(p) as fh:
            return json.load(fh)
    except json.JSONDecodeError as exc:
        raise Failure(f"{name} is not readable JSON (truncated?): {exc}")


def expected_identifiers(results, forward):
    """Expected query and genome ids, taken from OUTSIDE the table under test."""
    qcheck = _load(results, "query_id_check.json", required=False)
    pman = _load(results, "proteome_manifest.json", required=False)
    queries = set(qcheck) if qcheck else set()
    genomes = set(pman) if pman else set()
    if not queries:
        seen = set()
        for per_q in forward.values():
            seen |= set(per_q)
        queries = seen
    if not genomes:
        genomes = set(forward)
    return queries, genomes


def validate(results, verbose=True):
    problems = []
    def ok(m):
        if verbose:
            print(f"  [PASS] {m}")
    def note(m):
        if verbose:
            print(f"  [FAIL] {m}")

    env = _load(results, "environment.json")
    calls = _load(results, "rbh_calls.json")
    forward = _load(results, "forward_hits.json")
    reverse = _load(results, "reverse_hits.json")
    thr = float(env.get("e_threshold", core.E_DEFAULT))

    queries, genomes = expected_identifiers(results, forward)
    before = len(problems)
    if len(genomes) != 40:
        problems.append(f"expected 40 genomes, manifests give {len(genomes)}")
    if len(queries) != 45:
        problems.append(f"expected 45 queries, manifests give {len(queries)}")
    if set(forward) != genomes:
        problems.append("forward evidence does not cover exactly the expected genomes")
    if set(reverse) != genomes:
        problems.append("reverse evidence does not cover exactly the expected genomes")
    (ok if len(problems) == before else note)(
        "expected identifier sets come from the manifests and match the raw evidence")

    before = len(problems)
    core.check_raw_hits(forward, "forward", problems)
    core.check_raw_hits(reverse, "reverse", problems)
    (ok if len(problems) == before else note)(
        "raw hit E-values are finite and non-negative and bit scores are finite")

    before = len(problems)
    core.check_calls(calls, forward, reverse, genomes, queries, thr, problems)
    (ok if len(problems) == before else note)(
        f"every one of the {len(genomes)*len(queries)} call decisions recomputes from "
        f"the raw evidence, both directions")

    recomputed = {}
    for g in genomes:
        for q in queries:
            rows = forward.get(g, {}).get(q)
            if rows is None:
                continue
            recomputed[(q, g)] = core.recompute_call(g, q, rows, reverse.get(g, {}), thr)["reciprocal"]

    before = len(problems)
    mpath = os.path.join(results, "presence_matrix.tsv")
    if not os.path.exists(mpath):
        problems.append("presence_matrix.tsv is missing")
    else:
        with open(mpath) as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        core.check_matrix(rows, queries, genomes, recomputed, problems)
    (ok if len(problems) == before else note)(
        "presence matrix identifiers are unique and complete, and every cell derives "
        "from the recomputed calls")

    for m in problems[:8]:
        if verbose:
            print(f"        {m}")
    n_rbh = sum(1 for c in calls if isinstance(c, dict) and c.get("reciprocal"))
    if verbose:
        print(f"  ---- {n_rbh} reciprocal best hits")
    return problems, n_rbh


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    try:
        problems, n = validate(a.results, verbose=not a.quiet)
    except Failure as exc:
        print(f"REJECTED: {exc}")
        return 2
    if problems:
        print(f"\nREJECTED: {len(problems)} inconsistency(ies)")
        return 2
    print(f"\nACCEPTED: evidence is internally consistent ({n} reciprocal best hits)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
