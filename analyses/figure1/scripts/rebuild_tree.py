#!/usr/bin/env python3
"""Produce the 40-tip Figure 1 tree by distance-preserving pruning.

The input is the recovered 41-tip KBase Newick. This is a correction of the
export/pruning step, not a new phylogenetic inference. It removes SPH-1
(GCF_000018665.1), retains the original root, and preserves both branch distances
and the distinction between a reported support value and a missing label.

The historical KBase pruning discarded a 0.00643 stem shared by SPH-1 and the
three retained D. acidovorans genomes. Consequently, 111 of the 780 retained-tip
distances were shortened by 0.00643. Its ETE serialization also replaced ten
missing support labels with 1. This script preserves distances when suppressing
unary nodes and never transfers support from a removed clade to its survivor.

Only a fresh or empty output directory is accepted. Repository inputs are read
only. The accompanying tests retain the historical output for comparison.
"""

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_DEFAULT = ROOT / "data" / "tree_41tips_source.newick"
DROP_DEFAULT = "GCF_000018665.1"


class Node:
    __slots__ = ("name", "length", "children", "parent", "labelled")

    def __init__(self):
        self.name = None
        self.length = 0.0
        self.children = []
        self.parent = None
        self.labelled = False


def parse_newick(text):
    """Parse this numeric-support Newick format, preserving label presence.

    The archived input uses unquoted tip names, numeric internal labels, and
    branch lengths. This is not intended as a general-purpose Newick parser.
    """
    s = text.strip().rstrip(";")
    i = 0

    def node():
        nonlocal i
        n = Node()
        if s[i] == "(":
            i += 1
            while True:
                c = node()
                c.parent = n
                n.children.append(c)
                if s[i] == ",":
                    i += 1
                else:
                    break
            assert s[i] == ")", s[i]
            i += 1
        j = i
        while i < len(s) and s[i] not in ",():":
            i += 1
        lab = s[j:i]
        n.labelled = lab != ""
        n.name = lab or None
        if i < len(s) and s[i] == ":":
            i += 1
            j = i
            while i < len(s) and (s[i].isdigit() or s[i] in ".-+eE"):
                i += 1
            n.length = float(s[j:i])
        return n

    root = node()
    root.length = 0.0
    return root


def walk(n, acc=None):
    acc = [] if acc is None else acc
    acc.append(n)
    for c in n.children:
        walk(c, acc)
    return acc


def leaves(n):
    return [x for x in walk(n) if not x.children]


def fmt_len(v):
    """Fixed decimal with trailing zeros trimmed; no scientific notation."""
    s = f"{v:.10f}".rstrip("0")
    return s + "0" if s.endswith(".") else s


def write_newick(n):
    """Serialize without inventing labels for unlabelled internal nodes."""
    if not n.children:
        return f"{n.name}:{fmt_len(n.length)}"
    inner = ",".join(write_newick(c) for c in n.children)
    lab = n.name if n.labelled and n.name else ""
    if n.parent is None:
        return f"({inner}){lab};"
    return f"({inner}){lab}:{fmt_len(n.length)}"


def prune_preserving(root, drop):
    """Remove one tip, retaining distances, the root, and surviving labels.

    Suppressed unary nodes contribute their branch lengths to the surviving
    child, but never their support labels: those described a different clade.
    """
    targets = [leaf for leaf in leaves(root) if leaf.name == drop]
    if not targets:
        raise SystemExit(f"leaf {drop!r} not found")
    n = targets[0]
    p = n.parent
    p.children.remove(n)
    while p is not None and len(p.children) == 1 and p.parent is not None:
        only = p.children[0]
        only.length += p.length
        only.parent = p.parent
        idx = p.parent.children.index(p)
        p.parent.children[idx] = only
        p = only.parent
    if root.children and len(root.children) == 1:
        raise SystemExit("root became unary; refusing to guess a new root")
    return root


def pairwise(root):
    """Return unordered tip-pair distances, rounded to ten decimal places."""
    adj = {}
    for x in walk(root):
        adj.setdefault(id(x), [])
        for c in x.children:
            adj[id(x)].append((id(c), c.length))
            adj.setdefault(id(c), []).append((id(x), c.length))
    tips = leaves(root)
    out = {}
    for a in tips:
        dist = {id(a): 0.0}
        stack = [id(a)]
        while stack:
            u = stack.pop()
            for v, weight in adj[u]:
                if v not in dist:
                    dist[v] = dist[u] + weight
                    stack.append(v)
        for b in tips:
            if a.name < b.name:
                out[(a.name, b.name)] = round(dist[id(b)], 10)
    return out


def splits(root):
    all_tips = frozenset(leaf.name for leaf in leaves(root))
    out = set()
    for x in walk(root):
        if x is root:
            continue
        clade = frozenset(leaf.name for leaf in leaves(x))
        if 0 < len(clade) < len(all_tips):
            out.add(frozenset([clade, all_tips - clade]))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", type=Path, default=SOURCE_DEFAULT,
        help="41-tip source Newick (default: module data input)",
    )
    parser.add_argument(
        "--out-dir", type=Path, required=True,
        help="fresh or empty output directory; existing results are never overwritten",
    )
    parser.add_argument("--drop", default=DROP_DEFAULT)
    args = parser.parse_args()
    if args.out_dir.exists() and (
        not args.out_dir.is_dir() or any(args.out_dir.iterdir())
    ):
        parser.error(f"--out-dir is not empty or is not a directory: {args.out_dir}")

    # Read and validate before creating the destination directory.
    raw_bytes = args.source.read_bytes()
    raw = raw_bytes.decode("utf-8")
    source = parse_newick(raw)
    source_distances = pairwise(source)
    fixed = prune_preserving(parse_newick(raw), args.drop)
    out_newick = write_newick(fixed)
    output_bytes = (out_newick + "\n").encode("utf-8")
    retained_distances = pairwise(fixed)
    equal_pairs = sum(
        abs(distance - source_distances[pair]) < 1e-12
        for pair, distance in retained_distances.items()
    )
    if equal_pairs != len(retained_distances):
        raise SystemExit("distance preservation failed; no output written")
    report = {
        "source_file_sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "source_string_sha256": hashlib.sha256(raw.strip().encode()).hexdigest(),
        "dropped": args.drop,
        "tips_in": len(leaves(source)),
        "tips_out": len(leaves(fixed)),
        "retained_pairs": len(retained_distances),
        "distances_equal_to_source": equal_pairs,
        "unlabelled_internal_nodes": sum(
            1
            for x in walk(fixed)
            if x.children and x is not fixed and not x.labelled
        ),
        "root_children_lengths": [c.length for c in fixed.children],
        "output_string_sha256": hashlib.sha256(out_newick.encode()).hexdigest(),
        "output_file_sha256": hashlib.sha256(output_bytes).hexdigest(),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "tree_40tips_corrected.newick").write_bytes(output_bytes)
    (args.out_dir / "rebuild_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
