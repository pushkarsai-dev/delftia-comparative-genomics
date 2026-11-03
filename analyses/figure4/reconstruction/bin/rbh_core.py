#!/usr/bin/env python3
"""Shared rules for Figure 4c RBH evidence.

One module so that a single-genome checkpoint and the final aggregated evidence
are judged by exactly the same rules. Revision 3.
"""
import math

E_DEFAULT = 1e-5
CALL_KEYS = {"genome", "query", "forward_hit", "forward_evalue", "forward_score",
             "forward_ambiguous", "forward_tied_with", "reverse_hit", "reverse_evalue",
             "reverse_score", "reverse_ambiguous", "reverse_tied_with", "reciprocal"}


class Failure(Exception):
    pass


def finite_number(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(float(v))


def valid_evalue(v):
    """Finite, a real number, and non-negative. Negative E-values are impossible."""
    return finite_number(v) and float(v) >= 0.0


def valid_score(v):
    """Bit scores must be finite. They may legitimately be negative."""
    return finite_number(v)


def pick_best(rows):
    """Lowest E, then highest bit score. Exact remaining ties are ambiguous."""
    if not rows:
        return None, False, []
    ordered = sorted(rows, key=lambda r: (r["evalue"], -r["score"]))
    best = ordered[0]
    tied = [r for r in ordered
            if r["evalue"] == best["evalue"] and r["score"] == best["score"]]
    return best, (len(tied) > 1), tied


def check_raw_hits(store, side, problems):
    """Type and range checks on a {genome: {query: [hit, ...]}} store."""
    for g, per_q in store.items():
        if not isinstance(per_q, dict):
            problems.append(f"{side} evidence for {g} is not an object")
            continue
        for qid, rows in per_q.items():
            if not isinstance(rows, list):
                problems.append(f"{side} evidence {g}/{qid} is not a list")
                continue
            for r in rows:
                if not isinstance(r, dict) or "target" not in r:
                    problems.append(f"{side} hit in {g}/{qid} is malformed")
                    continue
                if not valid_evalue(r.get("evalue")):
                    problems.append(f"{side} hit {g}/{qid}->{r['target']} has an "
                                    f"invalid E-value {r.get('evalue')!r}")
                if not valid_score(r.get("score")):
                    problems.append(f"{side} hit {g}/{qid}->{r['target']} has a "
                                    f"non-finite bit score {r.get('score')!r}")


def recompute_call(genome, qid, fwd_rows, rev_store, threshold):
    """The single source of truth for what a call record must contain."""
    fbest, fambig, ftied = pick_best(fwd_rows)
    exp = dict(genome=genome, query=qid,
               forward_hit=None, forward_evalue=None, forward_score=None,
               forward_ambiguous=fambig,
               forward_tied_with=[r["target"] for r in ftied[1:]] if fambig else [],
               reverse_hit=None, reverse_evalue=None, reverse_score=None,
               reverse_ambiguous=False, reverse_tied_with=[], reciprocal=False)
    if fbest is None:
        return exp
    exp.update(forward_hit=fbest["target"], forward_evalue=fbest["evalue"],
               forward_score=fbest["score"])
    rbest, rambig, rtied = pick_best(rev_store.get(fbest["target"], []))
    if rbest is not None:
        exp.update(reverse_hit=rbest["target"], reverse_evalue=rbest["evalue"],
                   reverse_score=rbest["score"], reverse_ambiguous=rambig,
                   reverse_tied_with=[r["target"] for r in rtied[1:]] if rambig else [])
        exp["reciprocal"] = (rbest["target"] == qid) and not fambig and not rambig
    return exp


def check_calls(calls, forward, reverse, expected_genomes, expected_queries,
                threshold, problems):
    """Both directions, symmetrically, against the raw best hits.

    `expected_queries` and `expected_genomes` must come from OUTSIDE the table
    being validated.
    """
    seen = {}
    for c in calls:
        if not isinstance(c, dict) or not CALL_KEYS <= set(c):
            problems.append("a call record is missing required fields")
            continue
        key = (c["genome"], c["query"])
        seen[key] = seen.get(key, 0) + 1
        if c["genome"] not in expected_genomes:
            problems.append(f"call names an unexpected genome {c['genome']!r}")
            continue
        if c["query"] not in expected_queries:
            problems.append(f"call names an unexpected query {c['query']!r}")
            continue
        rows = forward.get(c["genome"], {}).get(c["query"])
        if rows is None:
            problems.append(f"no forward evidence for {c['genome']}/{c['query']}")
            continue
        exp = recompute_call(c["genome"], c["query"],
                             rows, reverse.get(c["genome"], {}), threshold)
        for f in sorted(CALL_KEYS):
            if c.get(f) != exp[f]:
                problems.append(f"{c['genome']}/{c['query']} {f}: recorded "
                                f"{c.get(f)!r}, recomputed {exp[f]!r}")
        # recorded numbers must themselves be well formed and below threshold
        for side in ("forward", "reverse"):
            ev, sc = c.get(f"{side}_evalue"), c.get(f"{side}_score")
            if ev is not None:
                if not valid_evalue(ev):
                    problems.append(f"{c['genome']}/{c['query']} {side}_evalue "
                                    f"{ev!r} is not a valid E-value")
                elif ev >= threshold:
                    problems.append(f"{c['genome']}/{c['query']} {side}_evalue "
                                    f"{ev} is not below {threshold}")
            if sc is not None and not valid_score(sc):
                problems.append(f"{c['genome']}/{c['query']} {side}_score "
                                f"{sc!r} is not finite")
        if c["reciprocal"]:
            if c["reverse_hit"] != c["query"]:
                problems.append(f"{c['genome']}/{c['query']} marked reciprocal but "
                                f"the reverse hit is not the query")
            if c["reverse_evalue"] is None:
                problems.append(f"{c['genome']}/{c['query']} reciprocal without a "
                                f"reverse hit below threshold")
            if c["forward_ambiguous"] or c["reverse_ambiguous"]:
                problems.append(f"{c['genome']}/{c['query']} counts an ambiguous tie "
                                f"as a unique reciprocal best hit")
        if c["forward_hit"] is not None:
            targets = {r["target"] for r in rows}
            if c["forward_hit"] not in targets:
                problems.append(f"{c['genome']}/{c['query']} calls target "
                                f"{c['forward_hit']!r}, absent from the raw forward hits")
    dups = [k for k, v in seen.items() if v > 1]
    if dups:
        problems.append(f"{len(dups)} duplicate (genome, query) call records, "
                        f"e.g. {dups[0]}")
    missing = {(g, q) for g in expected_genomes for q in expected_queries} - set(seen)
    if missing:
        problems.append(f"{len(missing)} (genome, query) pairs have no call record, "
                        f"e.g. {sorted(missing)[0]}")
    return {k: v for k, v in seen.items()}


def check_matrix(rows, expected_queries, expected_genomes, recomputed, problems):
    """Identifiers first, cells second. Expected ids come from outside the table."""
    if not rows:
        problems.append("presence matrix is empty")
        return
    header = rows[0]
    if header[0] != "query":
        problems.append("presence matrix header does not start with 'query'")
    cols = header[1:]
    if len(cols) != len(set(cols)):
        dup = sorted({c for c in cols if cols.count(c) > 1})
        problems.append(f"presence matrix has duplicate genome columns: {dup[:3]}")
    if set(cols) != set(expected_genomes):
        problems.append("presence matrix columns are not exactly the expected genomes")
    qids = [r[0] for r in rows[1:]]
    if len(qids) != len(set(qids)):
        dup = sorted({q for q in qids if qids.count(q) > 1})
        problems.append(f"presence matrix has duplicate query rows: {dup[:3]}")
    if set(qids) != set(expected_queries):
        missing = sorted(set(expected_queries) - set(qids))
        extra = sorted(set(qids) - set(expected_queries))
        problems.append(f"presence matrix rows are not exactly the expected queries; "
                        f"missing {missing[:3]}, unexpected {extra[:3]}")
    for r in rows[1:]:
        if len(r) != len(header):
            problems.append(f"presence matrix row {r[0]!r} has {len(r)} fields, "
                            f"expected {len(header)}")
            continue
        for v in r[1:]:
            if v not in ("0", "1"):
                problems.append(f"presence matrix row {r[0]!r} holds a non-binary "
                                f"value {v!r}")
                break
    bad = 0
    for r in rows[1:]:
        if len(r) != len(header):
            continue
        for g, v in zip(cols, r[1:]):
            if v not in ("0", "1"):
                continue
            key = (r[0], g)
            if key not in recomputed or int(v) != int(bool(recomputed[key])):
                bad += 1
    if bad:
        problems.append(f"{bad} presence matrix cells disagree with the recomputed calls")


def check_checkpoint(obj, genome, expected_queries, problems):
    """Deep validation of a single-genome checkpoint, same rules as the whole."""
    required = {"calls", "forward", "reverse", "fingerprint", "n_forward_db", "n_reverse_db"}
    if not isinstance(obj, dict) or not required <= set(obj):
        problems.append(f"checkpoint for {genome} does not match the expected schema")
        return
    for k in ("n_forward_db", "n_reverse_db"):
        v = obj[k]
        if not isinstance(v, int) or isinstance(v, bool) or v <= 0:
            problems.append(f"checkpoint for {genome} has an invalid {k}: {v!r}")
    if not isinstance(obj["forward"], dict) or not obj["forward"]:
        problems.append(f"checkpoint for {genome} has empty or malformed forward evidence")
    if not isinstance(obj["reverse"], dict):
        problems.append(f"checkpoint for {genome} has malformed reverse evidence")
    fwd_q = set(obj["forward"]) if isinstance(obj["forward"], dict) else set()
    if fwd_q != set(expected_queries):
        problems.append(f"checkpoint for {genome} forward evidence covers {len(fwd_q)} "
                        f"queries, expected {len(expected_queries)}")
    wrong = [c for c in obj["calls"] if isinstance(c, dict) and c.get("genome") != genome]
    if wrong:
        problems.append(f"checkpoint for {genome} holds {len(wrong)} calls belonging "
                        f"to another genome")
    check_raw_hits({genome: obj["forward"]}, "forward", problems)
    check_raw_hits({genome: obj["reverse"]}, "reverse", problems)
    check_calls(obj["calls"], {genome: obj["forward"]}, {genome: obj["reverse"]},
                {genome}, set(expected_queries), E_DEFAULT, problems)
