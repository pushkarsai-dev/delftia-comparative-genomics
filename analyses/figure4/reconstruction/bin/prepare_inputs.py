#!/usr/bin/env python3
"""Restore the hash-verified FASTA inputs into a fresh staging directory."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def prepare(out, base=None):
    base = base or Path(__file__).resolve().parents[1]
    out = Path(out)
    if out.exists():
        raise ValueError("Output already exists; choose a new directory.")
    inventory = json.loads((base / "evidence/input_inventory.json").read_text())
    payloads = []
    for genome, info in sorted(inventory["assemblies"].items()):
        if Path(genome).name != genome:
            raise ValueError("Invalid genome identifier")
        src = base / "inputs/assemblies" / (genome + ".fna.gz")
        data = gzip.decompress(src.read_bytes())
        if len(data) != info["bytes"] or hashlib.sha256(data).hexdigest() != info["sha256"]:
            raise ValueError(f"Assembly checksum mismatch: {genome}")
        payloads.append((Path(inventory["panel_dir"]) / (genome + ".fna"), data))
    query = (base.parent / "data/region002.faa").read_bytes()
    if hashlib.sha256(query).hexdigest() != inventory["query_fasta"]["sha256"]:
        raise ValueError("Query FASTA checksum mismatch")
    payloads.append((Path(inventory["query_fasta"]["path"]), query))
    for rel, _ in payloads:
        if rel.is_absolute() or ".." in rel.parts:
            raise ValueError("Unsafe inventory path")
    out.mkdir(parents=True, exist_ok=False)
    for rel, data in payloads:
        path = out / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(data)
    return len(payloads)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        count = prepare(args.output_dir)
    except (ValueError, OSError) as error:
        parser.exit(2, f"REFUSED: {error}\n")
    print(f"Restored {count} verified FASTA files in {args.output_dir}")
