#!/usr/bin/env python3
"""Rewrite the abstract block of paper/biorxiv_submission.md from paper.tex.

The submission metadata drifted from the manuscript twice during revision.
Extracting it mechanically removes that failure mode: paper.tex is the
single source, this script is the only writer of that block.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEX = ROOT / "paper" / "paper.tex"
MD = ROOT / "paper" / "biorxiv_submission.md"

SECTIONS = ["Background", "Methods", "Results", "Conclusions"]

REPL = [
    (r"\\textit\{([^}]*)\}", r"*\1*"), (r"\\emph\{([^}]*)\}", r"*\1*"),
    (r"\\citep\{[^}]*\}", ""), (r"\\citeyearpar\{[^}]*\}", ""),
    (r"\\S\\ref\{[^}]*\}", ""), (r"\\ref\{[^}]*\}", ""),
    (r"\$P\(v\) \\propto S_\{\\mathrm\{mutagen\}\}\(v\) \\cdot \\exp\(\\beta\\sigma\(v\)\)\$",
     "P(v) ∝ S_mutagen(v)·exp(βσ(v))"),
    (r"\$\\log S \+ \\sigma\$", "log S + σ"), (r"\$\\log S\$", "log S"),
    (r"\$\\sigma\(v\)\$", "σ(v)"), (r"\$\\sigma\$", "σ"), (r"\$\\beta\$", "β"),
    (r"\$S_\{\\mathrm\{mutagen\}\}\$", "S_mutagen"),
    (r"\$P = 3\.2 \\times 10\^\{-5\}\$", "P = 3.2 × 10⁻⁵"),
    (r"\$n = 4\$--11", "n = 4–11"), (r"\$K\$", "K"), (r"top-\$K\$", "top-K"),
    (r"C\$>\$G", "C>G"), (r"\\%", "%"), (r"--", "–"),
    (r"\\ ", " "), (r"\\\\", " "), (r"~", " "),
    (r"\{|\}", ""), (r"\$", ""),
]


def tex_to_md(block: str) -> str:
    for pat, rep in REPL:
        block = re.sub(pat, rep, block)
    return re.sub(r"\s+", " ", block).strip()


def main() -> int:
    tex = TEX.read_text()
    abs_start = tex.index(r"\section*{Abstract}")
    abs_end = tex.index(r"\paragraph{Keywords.}")
    abstract = tex[abs_start:abs_end]

    parts = []
    for i, name in enumerate(SECTIONS):
        start = abstract.index("\\paragraph{%s.}" % name) + len("\\paragraph{%s.}" % name)
        end = (abstract.index("\\paragraph{%s.}" % SECTIONS[i + 1])
               if i + 1 < len(SECTIONS) else len(abstract))
        parts.append(f"> **{name}.** {tex_to_md(abstract[start:end])}")

    new_block = "\n>\n".join(parts) + "\n"

    md = MD.read_text()
    head = md.index("> **Background.**")
    tail = md.index("## 5. Keywords")
    MD.write_text(md[:head] + new_block + "\n" + md[tail:])
    words = sum(len(p.split()) for p in parts)
    print(f"Synced abstract into {MD.name} ({words} words).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
