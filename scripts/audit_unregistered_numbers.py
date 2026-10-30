#!/usr/bin/env python3
"""Find quantitative claims in the prose that no registry entry accounts for.

audit_paper_numbers.py runs registry to paper: every registered value must
appear where the registry says. This script runs the other way: every number
written in the prose must be traceable to a registered value.

That is the direction nothing checked, and it is how a hand-typed "5.2 clones"
survived next to a registry that said 5.1.

Numbers are matched loosely, because the prose legitimately writes a value in
several forms: 1,427 or 1427, 0.20% or 0.2%, 4.7e-05 or $4.7 \\times 10^{-5}$.
A prose number counts as accounted for when it appears in some registered
value's formatted string, in its raw computed value, or in the small set of
literals below that are structural rather than empirical.

Exit status is non-zero if anything is unaccounted for, so it can gate a build.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "data" / "derived" / "paper_numbers.json"
PAPER = ROOT / "paper" / "paper.tex"

# Literals that are structural, definitional, or bibliographic rather than
# results: counts of things the paper names, model constants, years, and the
# six substitution classes' own arithmetic.
ALLOWED = {
    # Small integers used for counting sections, panels, classes and steps.
    "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
    "1.0", "2.0", "0.5", "100", "1000",
    # Years.
    *(str(y) for y in range(1990, 2031)),
    # Method constants fixed in Methods, not measured here.
    "0.25", "4", "95", "97.5", "5000", "20261013", "1e-12",
    "200",            # the third top-N cut of the N2 null, N = 50, 100, 200
    "16.7",           # 1/6, the uniform chemistry-free baseline
    "39",             # 189 - 151 + 1, the KRAS-4A hypervariable window
    "0.91",           # sqrt(0.82), derived in the sentence that states it
    # Figures quoted from the cited literature rather than computed here.
    "22500", "1210",  # Oh et al. atlas size, from that paper
    "587",            # size of the Arnold et al. germline ENU compilation
    # Values that appear in a table whose cells tests/test_manuscript_tables.py
    # checks against the derived tables directly, and are repeated in prose.
    "0.78", "0.44", "0.57", "493", "45",
    # Derived in the sentence that reports them, from registered inputs.
    "9.7", "2.6", "26.2", "7.0", "49.7",   # ENU C>G inversion, see §3.5
    "33", "1.1", "1.5", "16",               # sub-library enrichment, §3.6
    "2.1",                                  # 10^0.32, stated inline
}


# Lines that are not prose claims: comments, labels, citations, file paths,
# the bibliography command, and the figure include statements.
SKIP_LINE = re.compile(
    r"^\s*%|\\label\{|\\includegraphics|\\bibliography|\\cite|"
    r"\\usepackage|\\documentclass|\\newcommand|\\texttt\{[^}]*\d")

NUMBER = re.compile(r"(?<![\w.-])(\d[\d,]*\.?\d*)(?![\w-])")

# Compound identifiers (NVL-655, BLU-945, TPX-0131, SHA-256), software
# versions (0.11.10) and transcript accessions carry digits that are names,
# not measurements.
IDENTIFIER = re.compile(
    r"[A-Za-z]{2,}[-_]?\d[\w.-]*|\d+\.\d+\.\d+|NM\\?_\d+|\d+e[-+]?\d+"
    # ordinals such as "the 97.5th percentile": the number is a method
    # constant, and the suffix would otherwise truncate it to "97"
    r"|\d+(?:\.\d+)?(?:st|nd|rd|th)\b")


def norm(s: str) -> str:
    """Collapse a number to a comparable form."""
    s = s.replace(",", "").rstrip(".")
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s or "0"


def registry_forms(reg: dict) -> set[str]:
    """Every number any registered entry could legitimately license."""
    out: set[str] = set()
    for entry in reg.values():
        for field in ("formatted", "value"):
            v = entry.get(field)
            for piece in re.findall(r"\d[\d,]*\.?\d*", str(v)):
                out.add(norm(piece))
                # a value quoted to fewer decimals than it is stored with
                n = norm(piece)
                if "." in n:
                    whole, frac = n.split(".")
                    for k in range(1, len(frac)):
                        out.add(f"{whole}.{frac[:k]}")
                    out.add(whole)
        # inputs often carry the raw operands a sentence may quote
        for v in (entry.get("inputs") or {}).values():
            for piece in re.findall(r"\d[\d,]*\.?\d*", str(v)):
                out.add(norm(piece))
    return out


def body_sections(raw: str) -> list[tuple[str, int, str]]:
    """(section name, line number, line) for the prose we audit."""
    start = raw.index(r"\section{Background}")
    end = raw.index(r"\section*{Availability of data and materials}")
    offset = raw[:start].count("\n") + 1
    section = "Background"
    out = []
    in_table = False
    for i, line in enumerate(raw[start:end].split("\n")):
        m = re.match(r"\\(sub)?section\*?\{(.+?)\}", line)
        if m:
            section = m.group(2)[:46]
        # Table bodies are hand-typed LaTeX, but tests/test_manuscript_tables.py
        # checks their cells against the derived tables directly, which is a
        # stronger check than matching digits against the register.
        if re.search(r"\\begin\{tabular", line):
            in_table = True
        if re.search(r"\\end\{tabular", line):
            in_table = False
            continue
        if in_table or SKIP_LINE.search(line):
            continue
        out.append((section, offset + i, line))
    return out


def main() -> int:
    reg = json.loads(REGISTRY.read_text())
    raw = PAPER.read_text()
    known = registry_forms(reg) | {norm(a) for a in ALLOWED}

    unaccounted: list[tuple[str, int, str, str]] = []
    checked = 0
    for section, lineno, line in body_sections(raw):
        # strip math-mode scaffolding that carries no claim
        probe = line.replace("{,}", "")          # 1{,}069 is one number
        probe = re.sub(r"\\(ref|label|cite\w*)\{[^}]*\}", " ", probe)
        probe = re.sub(r"\\[a-zA-Z]+", " ", probe)
        probe = IDENTIFIER.sub(" ", probe)
        for m in NUMBER.finditer(probe):
            tok = norm(m.group(1))
            checked += 1
            if tok not in known:
                unaccounted.append((section, lineno, m.group(1), line.strip()))

    print(f"numbers found in prose: {checked}   "
          f"unaccounted for: {len(unaccounted)}")
    if unaccounted:
        print()
        for section, lineno, tok, line in unaccounted:
            print(f"  line {lineno:>5}  {tok:>12}   [{section}]")
            print(f"              {line[:96]}")
        print("\nEach of these is either an unregistered claim or a literal "
              "that belongs in ALLOWED.")
        return 1
    print("\nEvery number in the prose is accounted for by the registry.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
