#!/usr/bin/env python3
"""Figure 4c whole-proteome RBH reconstruction - hardened runner (revision 3).

An explicitly documented RECONSTRUCTION, not a recovered execution. See README.

Revision 2 hardening:
  * a compatibility manifest binds every cache and checkpoint to the exact
    inputs, versions, and settings that produced it; reuse is refused when any
    of them changes
  * the proteome cache is JSON, never a pickle
  * per-genome checkpoints are written atomically and validated on load
  * the run ABORTS if any regenerated query differs from its archived sequence
  * duplicate FASTA or protein identifiers are rejected, never silently merged
  * effective builder/pipeline settings are read back from the installed API and
    separated from declarations that the API cannot confirm
  * actual per-genome forward and reverse database sizes are recorded

Revision 3:
  * the effective builder/pipeline configuration is part of the fingerprint, so a
    changed setting invalidates caches and checkpoints. Revision 2's README
    claimed this; the code did not do it.
  * checkpoints are validated deeply, by the same rbh_core rules as the final
    evidence: unique call pairs, genome identity, complete forward-query
    coverage, reverse evidence present, valid database sizes, and every call
    consistent with the raw hits
  * the checkpoint directory must hold only expected genomes; unknown JSON is
    refused rather than aggregated
  * existing-output, checkpoint, and cache compatibility checks precede output
    initialization. A populated directory without a valid manifest is refused
    unchanged. Query identity is checked later, before searches.
"""
import argparse, hashlib, json, os, platform, sys, tempfile, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rbh_core as core

RUNNER_REVISION = "3.0.0"
E_THRESHOLD = 1e-5
SEED = 42
QUERY_GENOME = "BAA554_hybrid"
SUPPORTED = {"pyhmmer": "0.12.2", "pyrodigal": "3.7.1"}


class Refused(Exception):
    pass


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_fasta(path):
    """Yield (id, seq). Raises on duplicate identifiers."""
    seen, name, buf = set(), None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name = line[1:].split()[0]
                if name in seen:
                    raise Refused(f"duplicate FASTA identifier {name!r} in {path}")
                seen.add(name)
                buf = []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def predict_proteome(fna):
    import pyrodigal
    recs = list(read_fasta(fna))
    finder = pyrodigal.GeneFinder(meta=False)
    finder.train(*[s.encode() for _, s in recs])
    out = {}
    for contig, seq in recs:
        for j, gene in enumerate(finder.find_genes(seq.encode()), start=1):
            pid = f"{contig}_{j}"
            if pid in out:
                raise Refused(f"duplicate protein identifier {pid!r} in {fna}")
            out[pid] = dict(seq=gene.translate().rstrip("*"), contig=contig,
                            begin=gene.begin, end=gene.end, strand=gene.strand)
    return out


def effective_settings():
    """Read back what the installed API actually reports; separate the rest."""
    import pyhmmer
    from pyhmmer.easel import Alphabet
    from pyhmmer.plan7 import Builder, Pipeline, Background
    al = Alphabet.amino()
    b = Builder(al, seed=SEED)
    bg = Background(al)
    pl = Pipeline(al, background=bg, seed=SEED)
    observed, declared = {}, {}
    for k in ("seed", "architecture", "weighting", "effective_number", "prior_scheme",
              "popen", "pextend", "score_matrix", "window_length", "window_beta"):
        v = getattr(b, k, None)
        (observed if v is not None else declared)[f"builder.{k}"] = v
    for k in ("bias_filter", "null2", "seed", "F1", "F2", "F3"):
        v = getattr(pl, k, None)
        (observed if v is not None else declared)[f"pipeline.{k}"] = v
    observed["background.uniform"] = bg.uniform
    observed["background.L"] = bg.L
    declared["builder.symfrac"] = "0.5 (constructor default; not readable back)"
    declared["builder.fragthres"] = "0.5 (constructor default; not readable back)"
    declared["pipeline.Z"] = "set per search to the actual database size"
    return observed, declared


def compatibility_manifest(panel_dir, fnas, qfaa, reverse_db):
    import pyhmmer, pyrodigal
    observed, declared = effective_settings()
    return {
        "effective_settings_observed": observed,
        "effective_settings_declared": declared,
        "runner_revision": RUNNER_REVISION,
        "versions": {"pyhmmer": pyhmmer.__version__, "pyrodigal": pyrodigal.__version__,
                     "python": platform.python_version()},
        "prediction": {"tool": "pyrodigal", "meta": False, "training": "per-genome",
                       "id_scheme": "{contig}_{1-based ordinal}"},
        "search": {"profiles": "single-sequence, pyhmmer Builder",
                   "directions": "hmmsearch both ways", "seed": SEED,
                   "e_threshold": E_THRESHOLD, "reverse_db": reverse_db,
                   "ranking": "lowest sequence E, then highest bit score; "
                              "exact ties ambiguous and excluded"},
        "query_fasta": {"path": os.path.relpath(qfaa, os.path.dirname(panel_dir)),
                        "sha256": sha256_file(qfaa)},
        "inputs": {f[:-4]: sha256_file(os.path.join(panel_dir, f)) for f in sorted(fnas)},
    }


def manifest_fingerprint(m):
    return hashlib.sha256(json.dumps(m, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def atomic_write_json(path, obj):
    d = os.path.dirname(path)
    fd, tmp = tempfile.mkstemp(dir=d, suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(obj, fh)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise


CHECKPOINT_KEYS = {"calls", "forward", "reverse", "fingerprint", "n_forward_db", "n_reverse_db"}


def load_checkpoint(path, fingerprint, genome, expect_queries):
    """Return the checkpoint or raise Refused. Never repairs, never overwrites.

    Revision 3 applies the same rbh_core rules used for the final evidence, so a
    checkpoint whose query NAMES happen to match but whose raw evidence is empty
    or whose calls are duplicated is refused.
    """
    try:
        with open(path) as fh:
            obj = json.load(fh)
    except json.JSONDecodeError as exc:
        raise Refused(f"checkpoint for {genome} is truncated or malformed: {exc}")
    if not isinstance(obj, dict) or not CHECKPOINT_KEYS <= set(obj):
        raise Refused(f"checkpoint for {genome} does not match the expected schema")
    if obj["fingerprint"] != fingerprint:
        raise Refused(
            f"checkpoint for {genome} was produced under different inputs or settings.\n"
            f"    Refusing to mix runs. Use a fresh --outdir.")
    problems = []
    core.check_checkpoint(obj, genome, expect_queries, problems)
    if problems:
        raise Refused(f"checkpoint for {genome} failed validation: {problems[0]}"
                      + (f" (and {len(problems)-1} more)" if len(problems) > 1 else ""))
    return obj


def audit_output_dir(outdir, fingerprint, expected_genomes, expect_queries):
    """Validate an existing output directory WITHOUT writing anything."""
    if not os.path.isdir(outdir) or not os.listdir(outdir):
        return
    mpath = os.path.join(outdir, "compatibility_manifest.json")
    if not os.path.exists(mpath):
        raise Refused(f"--outdir is populated but holds no compatibility_manifest.json: "
                      f"{outdir}\n    Refusing to write into it. Use a fresh --outdir.")
    try:
        with open(mpath) as fh:
            prev = json.load(fh)
    except json.JSONDecodeError as exc:
        raise Refused(f"compatibility_manifest.json is truncated or malformed: {exc}")
    if not isinstance(prev, dict) or "manifest" not in prev:
        raise Refused("compatibility_manifest.json does not match the expected schema")
    if manifest_fingerprint(prev["manifest"]) != fingerprint:
        raise Refused(
            "this output directory was produced under different inputs, versions, or "
            "settings.\n    Refusing to mix runs. Use a fresh --outdir.")
    ckpt = os.path.join(outdir, "per_genome")
    if os.path.isdir(ckpt):
        names = sorted(os.listdir(ckpt))
        unknown = [n for n in names if not n.endswith(".json") or n[:-5] not in expected_genomes]
        if unknown:
            raise Refused(f"checkpoint directory holds unexpected entries: {unknown[:3]}")
        for n in names:
            load_checkpoint(os.path.join(ckpt, n), fingerprint, n[:-5], expect_queries)


def main():
    from pyhmmer.easel import Alphabet, TextSequence, DigitalSequenceBlock
    from pyhmmer.plan7 import Builder, Background, Pipeline
    import pyhmmer, pyrodigal

    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--reverse-db", choices=["proteome", "region45"], default="proteome")
    ap.add_argument("--proteome-cache", default=None, help="JSON cache; never a pickle")
    ap.add_argument("--time-budget", type=float, default=0.0)
    a = ap.parse_args()

    for pkg, want in SUPPORTED.items():
        got = {"pyhmmer": pyhmmer.__version__, "pyrodigal": pyrodigal.__version__}[pkg]
        if got != want:
            raise Refused(f"{pkg} {got} is not the supported version {want}")

    R = a.project_root
    panel = os.path.join(R, "A3 antiSMASH Local (2 Aug 2026)", "panel40_fasta")
    qfaa = os.path.join(R, "A11 Region Synteny (19 Aug 2026)", "region002.faa")
    fnas = sorted(f for f in os.listdir(panel) if f.endswith(".fna"))
    if len(fnas) != 40:
        raise Refused(f"expected 40 panel FASTAs, found {len(fnas)}")

    man = compatibility_manifest(panel, fnas, qfaa, a.reverse_db)
    fp = manifest_fingerprint(man)
    expected_genomes = {f[:-4] for f in fnas}

    # The archived queries are read here so the checkpoint audit can use them.
    archived_probe = dict(read_fasta(qfaa))

    # ---- audit existing outputs and cache compatibility before initialization ----
    audit_output_dir(a.outdir, fp, expected_genomes, set(archived_probe))
    if a.proteome_cache and os.path.exists(a.proteome_cache):
        try:
            with open(a.proteome_cache) as fh:
                _probe = json.load(fh)
        except json.JSONDecodeError as exc:
            raise Refused(f"proteome cache is truncated or malformed: {exc}")
        if not isinstance(_probe, dict) or _probe.get("fingerprint") != fp:
            raise Refused("proteome cache was built from different inputs or settings")
        if sorted(_probe.get("proteomes", {})) != sorted(expected_genomes):
            raise Refused("proteome cache covers a different genome set")

    os.makedirs(a.outdir, exist_ok=True)
    mpath = os.path.join(a.outdir, "compatibility_manifest.json")
    if not os.path.exists(mpath):
        atomic_write_json(mpath, {"fingerprint": fp, "manifest": man})
    ckpt = os.path.join(a.outdir, "per_genome")
    os.makedirs(ckpt, exist_ok=True)

    # ---- proteomes, JSON cache bound to the fingerprint ---------------------
    proteomes = None
    if a.proteome_cache and os.path.exists(a.proteome_cache):
        try:
            with open(a.proteome_cache) as fh:
                blob = json.load(fh)
        except json.JSONDecodeError as exc:
            raise Refused(f"proteome cache is truncated or malformed: {exc}")
        if blob.get("fingerprint") != fp:
            raise Refused("proteome cache was built from different inputs or settings")
        if sorted(blob["proteomes"]) != [f[:-4] for f in fnas]:
            raise Refused("proteome cache covers a different genome set")
        proteomes = blob["proteomes"]
        print(f"  loaded {len(proteomes)} cached proteomes (fingerprint verified)", flush=True)
    if proteomes is None:
        proteomes = {}
        for f in fnas:
            g = f[:-4]
            proteomes[g] = predict_proteome(os.path.join(panel, f))
            print(f"  predicted {g}: {len(proteomes[g])}", flush=True)
        if a.proteome_cache:
            atomic_write_json(a.proteome_cache, {"fingerprint": fp, "proteomes": proteomes})

    # ---- queries: verify by identifier, ABORT on any difference -------------
    archived = dict(read_fasta(qfaa))
    qprot = proteomes[QUERY_GENOME]
    qcheck, bad = {}, []
    for qid, arch in archived.items():
        regen = qprot.get(qid)
        same = regen is not None and regen["seq"] == arch.rstrip("*")
        qcheck[qid] = dict(present_in_regeneration=regen is not None, sequence_identical=same)
        if not same:
            bad.append(qid)
    if bad:
        raise Refused(f"{len(bad)} regenerated queries differ from the archived sequences: "
                      f"{bad[:5]}\n    Refusing to search on mismatched queries.")
    queries = {q: qprot[q]["seq"] for q in archived}
    if len(queries) != 45:
        raise Refused(f"recovered {len(queries)} of 45 query proteins")

    alphabet = Alphabet.amino()
    background = Background(alphabet)

    def block(prots):
        return DigitalSequenceBlock(alphabet, [
            TextSequence(name=k.encode(), sequence=v["seq"]).digitize(alphabet)
            for k, v in prots.items()])

    def profiles(items):
        b = Builder(alphabet, seed=SEED)
        out = {}
        for pid, seq in items:
            hmm, _, _ = b.build(
                TextSequence(name=pid.encode(), sequence=seq).digitize(alphabet), background)
            out[pid] = hmm
        return out

    def search(profs, blk, n):
        pl = Pipeline(alphabet, background=background, seed=SEED, Z=n)
        res = {}
        for pid, hmm in profs.items():
            pl.clear(); pl.Z = n
            rows = []
            for h in pl.search_hmm(hmm, blk):
                if h.evalue >= E_THRESHOLD:
                    continue
                d = h.best_domain
                nm = h.name.decode() if isinstance(h.name, bytes) else str(h.name)
                rows.append(dict(target=nm, evalue=float(h.evalue), score=float(h.score),
                                 bias=float(h.bias), dom_evalue=float(d.i_evalue),
                                 dom_score=float(d.score),
                                 env_from=int(d.env_from), env_to=int(d.env_to),
                                 ali_from=int(d.alignment.target_from),
                                 ali_to=int(d.alignment.target_to),
                                 hmm_from=int(d.alignment.hmm_from),
                                 hmm_to=int(d.alignment.hmm_to)))
            rows.sort(key=lambda r: (r["evalue"], -r["score"]))
            res[pid] = rows
        return res

    def pick(rows):
        if not rows:
            return None, False, []
        best = rows[0]
        tied = [r for r in rows if r["evalue"] == best["evalue"] and r["score"] == best["score"]]
        return best, len(tied) > 1, tied

    qprofiles = profiles(list(queries.items()))
    hdir = os.path.join(a.outdir, "profiles"); os.makedirs(hdir, exist_ok=True)
    with open(os.path.join(hdir, "queries.hmm"), "wb") as fh:
        for pid in sorted(qprofiles):
            qprofiles[pid].write(fh)

    rev_pool = proteomes[QUERY_GENOME] if a.reverse_db == "proteome" else {k: qprot[k] for k in queries}
    rev_block = block(rev_pool); n_rev = len(rev_pool)
    expect_q = set(queries)
    t0 = time.time()
    todo = [g for g in sorted(proteomes) if not os.path.exists(os.path.join(ckpt, g + ".json"))]
    for g in sorted(proteomes):
        p = os.path.join(ckpt, g + ".json")
        if os.path.exists(p):
            load_checkpoint(p, fp, g, expect_q)          # validates or refuses
            continue
        if a.time_budget and (time.time() - t0) > a.time_budget:
            break
        blk = block(proteomes[g]); n_fwd = len(proteomes[g])
        fwd = search(qprofiles, blk, n_fwd)
        bt = {}
        for qid, rows in fwd.items():
            b, _, _ = pick(rows)
            if b:
                bt[b["target"]] = proteomes[g][b["target"]]["seq"]
        rev = search(profiles(list(bt.items())), rev_block, n_rev)
        calls = []
        for qid, rows in fwd.items():
            fb, famb, ft = pick(rows)
            rec = dict(genome=g, query=qid, forward_hit=None, forward_evalue=None,
                       forward_score=None, forward_ambiguous=famb,
                       forward_tied_with=[r["target"] for r in ft[1:]] if famb else [],
                       reverse_hit=None, reverse_evalue=None, reverse_score=None,
                       reverse_ambiguous=False, reverse_tied_with=[], reciprocal=False,
                       contig=None, target_begin=None, target_end=None, target_strand=None,
                       query_len=len(queries[qid]), target_len=None,
                       query_coverage=None, target_coverage=None, n_forward_hits=len(rows))
            if fb:
                t = proteomes[g][fb["target"]]
                rec.update(forward_hit=fb["target"], forward_evalue=fb["evalue"],
                           forward_score=fb["score"], contig=t["contig"],
                           target_begin=t["begin"], target_end=t["end"], target_strand=t["strand"],
                           target_len=len(t["seq"]),
                           query_coverage=round((fb["hmm_to"]-fb["hmm_from"]+1)/len(queries[qid]), 4),
                           target_coverage=round((fb["ali_to"]-fb["ali_from"]+1)/len(t["seq"]), 4))
                rb, ramb, rt = pick(rev.get(fb["target"], []))
                if rb:
                    rec.update(reverse_hit=rb["target"], reverse_evalue=rb["evalue"],
                               reverse_score=rb["score"], reverse_ambiguous=ramb,
                               reverse_tied_with=[r["target"] for r in rt[1:]] if ramb else [])
                    rec["reciprocal"] = (rb["target"] == qid) and not famb and not ramb
            calls.append(rec)
        atomic_write_json(p, dict(calls=calls, forward=fwd, reverse=rev, fingerprint=fp,
                                  n_forward_db=n_fwd, n_reverse_db=n_rev))
        print(f"  searched {g}: {sum(1 for c in calls if c['reciprocal'])} RBH", flush=True)

    names = sorted(os.listdir(ckpt))
    unknown = [n for n in names if not n.endswith(".json") or n[:-5] not in expected_genomes]
    if unknown:
        raise Refused(f"checkpoint directory holds unexpected entries: {unknown[:3]}")
    done = sorted(n[:-5] for n in names)
    if len(done) < 40:
        print(f"\nCHECKPOINT: {len(done)}/40 complete. Re-run the same command to resume.")
        return 0

    forward, reverse, calls, dbsizes = {}, {}, [], {}
    for g in done:
        obj = load_checkpoint(os.path.join(ckpt, g + ".json"), fp, g, expect_q)
        forward[g] = obj["forward"]; reverse[g] = obj["reverse"]; calls.extend(obj["calls"])
        dbsizes[g] = dict(forward_db=obj["n_forward_db"], reverse_db=obj["n_reverse_db"])

    observed, declared = effective_settings()
    env = dict(
        runner_revision=RUNNER_REVISION,
        versions={"pyhmmer": pyhmmer.__version__, "pyrodigal": pyrodigal.__version__,
                  "python": platform.python_version(), "platform": platform.platform()},
        seed=SEED, e_threshold=E_THRESHOLD, reverse_db=a.reverse_db,
        n_reverse_db_sequences=n_rev,
        settings_observed_from_api=observed,
        settings_declared_not_readable_back=declared,
        per_genome_database_sizes=dbsizes,
        note=("search_wall_seconds_this_invocation covers only the genomes searched in "
              "THIS invocation; the run is checkpointed, so it is not total compute time"),
        search_wall_seconds_this_invocation=round(time.time()-t0, 1),
        profiles_deposited="the 45 forward query profiles only; reverse profiles are "
                           "rebuilt deterministically from the forward best hits")
    def dump(n, o):
        atomic_write_json(os.path.join(a.outdir, n), o)
    dump("environment.json", env)
    dump("query_id_check.json", qcheck)
    dump("proteome_manifest.json",
         {g: dict(n_proteins=len(p),
                  proteome_sha256=hashlib.sha256(
                      "".join(f">{k}\n{v['seq']}\n" for k, v in p.items()).encode()).hexdigest())
          for g, p in proteomes.items()})
    dump("forward_hits.json", forward)
    dump("reverse_hits.json", reverse)
    dump("rbh_calls.json", calls)
    import csv as _csv
    with open(os.path.join(a.outdir, "rbh_calls.tsv"), "w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=list(calls[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(calls)
    genomes = sorted(proteomes)
    with open(os.path.join(a.outdir, "presence_matrix.tsv"), "w", newline="") as fh:
        w = _csv.writer(fh, delimiter="\t"); w.writerow(["query"] + genomes)
        for qid in sorted(queries):
            row = {c["genome"]: c["reciprocal"] for c in calls if c["query"] == qid}
            w.writerow([qid] + [int(bool(row.get(g))) for g in genomes])
    print(f"\ntotal reciprocal best hits: {sum(1 for c in calls if c['reciprocal'])}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Refused as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        sys.exit(2)
