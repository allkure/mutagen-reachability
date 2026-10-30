#!/usr/bin/env python3
"""Generate SHA-256 checksums for every file in data/derived/.

Writes data/derived/CHECKSUMS.sha256 with one line per file in the
`sha256sum` output format: `<hex>  <relpath>`.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DERIVED = ROOT / "data" / "derived"
OUT = DERIVED / "CHECKSUMS.sha256"


def sha256_of(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    files = sorted(
        p for p in DERIVED.rglob("*")
        if p.is_file() and p != OUT
        and ".git" not in p.parts
        and ".pre_" not in p.name     # exclude per-revision backups
        and ".bak" not in p.name
    )
    lines = []
    for p in files:
        rel = p.relative_to(DERIVED)
        lines.append(f"{sha256_of(p)}  {rel}")
    OUT.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT}: {len(lines)} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
