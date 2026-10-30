#!/usr/bin/env python3
"""Check that every registered number actually appears in the manuscript.

The golden rule says each quantitative claim in paper.tex comes from the
registry. compute_paper_numbers.py enforces one direction: every registered
number is computed from source. This script enforces the other: a number that
is registered as appearing in the paper must actually be found there, in the
form the registry says.

That catches the failure mode this project hit repeatedly during revision,
where a recomputation moved a value and the prose kept the old one.

Entries whose `where` names no paper location are skipped, as are those whose
formatted value is too short to search for without matching unrelated digits.

Exit status is non-zero if anything is missing, so it can gate a build.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "data" / "derived" / "paper_numbers.json"
PAPER = ROOT / "paper" / "paper.tex"

# Locations in `where` that mean "this appears in the manuscript body".
PAPER_HINTS = ("abstract", "§", "table", "methods", "discussion",
               "background", "results", "conclusion")
# "Fig 2A" alone means the value is drawn inside a rendered figure, which is
# produced by its own script and is not text in paper.tex.

# Values too generic to search for literally.
MIN_LEN = 3
SKIP_KEYS = {"software_versions"}


def normalise(s: str) -> str:
    """Strip LaTeX noise so prose and registry strings compare fairly."""
    s = s.replace("\\%", "%").replace("\\,", "").replace("~", " ")
    s = s.replace("{,}", ",")
    # scientific notation: $4.7 \times 10^{-5}$ -> 4.7e-05
    s = re.sub(r"\$?([\d.]+)\s*\\times\s*10\^\{(-?\d+)\}\$?",
               lambda m: f"{float(m.group(1)):g}e{int(m.group(2)):+03d}", s)
    # fold ratios: $30.6\times$ -> 30.6x
    s = re.sub(r"([\d.]+)\s*\\times", r"\1x", s)
    s = s.replace("\\", "").replace("$", "")
    # thousands separators inside numbers: 1,069 -> 1069
    s = re.sub(r"(?<=\d),(?=\d\d\d\b)", "", s)
    s = re.sub(r"\s+", " ", s)
    return s.lower()


def variants(val: str) -> list[str]:
    """Forms a registry value may legitimately take in the prose."""
    out = {val}
    # en-dash range vs "to"
    out.add(val.replace("--", " to "))
    out.add(re.sub(r"\s+to\s+", "--", val))
    # trailing zero in a percentage: 0.20% == 0.2%
    m = re.fullmatch(r"(\d+)\.(\d*?)0+%", val)
    if m:
        out.add(f"{m.group(1)}.{m.group(2)}%".replace(".%", "%"))
    # scientific notation written either way
    m = re.fullmatch(r"([\d.]+)e([+-]\d+)", val)
    if m:
        out.add(f"{float(m.group(1)):g}e{int(m.group(2)):+03d}")
        out.add(f"{float(m.group(1)):g}e{int(m.group(2)):d}")
    # a recovery value reported with its interval, e.g. "63% [53%--72%]",
    # is written in the prose as the point value or as a range
    m = re.fullmatch(r"([\d.]+%) \[([\d.]+%)--([\d.]+%)\]", val)
    if m:
        out.update({m.group(1), f"{m.group(2)}--{m.group(3)}"})
    return sorted(out)


def split_paper(raw: str) -> tuple[str, str]:
    """Return (abstract, body) as normalised text.

    A value registered for the Abstract has to be in the abstract, not merely
    somewhere in the file. Checking the whole file let two wrong values stand
    in the Results for several revisions, because each happened to coincide
    with an unrelated figure elsewhere in the manuscript.
    """
    a0 = raw.index(r"\section*{Abstract}")
    a1 = raw.index(r"\section{Background}")
    return normalise(raw[a0:a1]), normalise(raw[a1:])


def main() -> int:
    if not REGISTRY.exists() or not PAPER.exists():
        sys.exit("run compute_paper_numbers.py first")
    reg = json.loads(REGISTRY.read_text())
    raw = PAPER.read_text()
    paper = normalise(raw)
    abstract, body = split_paper(raw)

    misplaced: list[tuple[str, str, list, str]] = []
    missing, checked, skipped = [], 0, 0
    for key, entry in sorted(reg.items()):
        if key in SKIP_KEYS:
            skipped += 1
            continue
        locs = [str(w) for w in entry.get("where", [])]
        # entries explicitly marked as not quoted, or quoted only in a
        # figure or table built outside paper.tex, are not checked here
        if any(l.lower().startswith(("not reported", "supporting")) for l in locs):
            skipped += 1
            continue
        where = " ".join(locs).lower()
        if not any(h in where for h in PAPER_HINTS):
            skipped += 1
            continue
        val = normalise(str(entry.get("formatted", "")))
        if len(val) < MIN_LEN:
            skipped += 1
            continue
        checked += 1
        forms = [normalise(v) for v in variants(val)]
        if not any(f in paper for f in forms):
            missing.append((key, entry["formatted"], entry.get("where", [])))
            continue

        # Location check. "Abstract" is enforced exactly; a section or table
        # reference is enforced only as "somewhere in the body", because the
        # §x.y numbers in `where` predate the current section numbering and
        # renumbering them is a separate job.
        wants_abstract = any("abstract" in l.lower() for l in locs)
        wants_body = any(h in where for h in PAPER_HINTS if h != "abstract")
        if wants_abstract and not any(f in abstract for f in forms):
            misplaced.append((key, entry["formatted"], entry.get("where", []),
                              "registered for the Abstract, not found there"))
        elif wants_body and not any(f in body for f in forms):
            misplaced.append((key, entry["formatted"], entry.get("where", []),
                              "registered for the body, found only in the Abstract"))

    print(f"registered: {len(reg)}   checked against paper.tex: {checked}   "
          f"skipped: {skipped}")
    rc = 0
    if missing:
        print(f"\n{len(missing)} registered value(s) not found in the manuscript:")
        for key, val, where in missing:
            print(f"  {key:40s} {val!r:30s} registered for {where}")
        print("\nEither the prose is stale or the `where` field is wrong.")
        rc = 1
    if misplaced:
        print(f"\n{len(misplaced)} registered value(s) in the wrong place:")
        for key, val, where, why in misplaced:
            print(f"  {key:40s} {val!r:24s} {why}")
            print(f"  {'':40s} registered for {where}")
        print("\nFix the `where` field, or put the value where it claims to be.")
        rc = 1
    if not rc:
        print("\nAll checked values appear in the manuscript, in the places "
              "the registry claims.")
    return rc


if __name__ == "__main__":
    sys.exit(main())
