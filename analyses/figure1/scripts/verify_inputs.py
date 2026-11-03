#!/usr/bin/env python3
"""Verify the bundled Figure 1 input bytes without modifying any files."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/input_checksums.json"


def verify():
    expected = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if not isinstance(expected, dict) or not expected:
        raise ValueError("Input manifest must be a nonempty object")
    for name, digest in expected.items():
        path = (ROOT / name).resolve()
        if not path.is_relative_to(ROOT / "data"):
            raise ValueError(f"Invalid input path: {name}")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != digest:
            raise ValueError(f"Input checksum mismatch: {name}")
    return len(expected)


if __name__ == "__main__":
    try:
        print(f"Verified {verify()} input files; no files changed.")
    except (OSError, ValueError) as error:
        raise SystemExit(f"Input verification failed: {error}")
