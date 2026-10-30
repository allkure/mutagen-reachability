#!/usr/bin/env python3
"""Generate notebooks/reproduce_figures.ipynb.

The notebook is the reader's entry point: it redraws the five main figures
and prints the manuscript's headline numbers beside the register that
produced them, so a reader can check the paper against the code without
reading the pipeline.

It is generated rather than hand-edited so that it cannot drift from the
scripts it calls. Edit this file, not the .ipynb.
"""
from __future__ import annotations

import json
from pathlib import Path

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "notebooks" / "reproduce_figures.ipynb"

MD = [
    # 0
    """# Reproducing the figures

This notebook redraws every main figure in *Chemical reachability limits
mutagen-based TKI-resistance screens* and prints the headline numbers next
to the register entry each one comes from.

It calls the same scripts the manuscript was built with. Nothing here
recomputes anything a different way, because a second implementation that
agreed would prove nothing and one that disagreed would not say which was
right. The checks that do that work live in `tests/` and `scripts/audit_*`.

**Runtime** is about a minute. **Prerequisite** is `pip install -r
requirements.txt` from the repository root.""",
    # 1
    """## 1. The register

Every quantitative claim in the manuscript is computed by a named script
and recorded in `data/derived/paper_numbers.json`, with its source file,
the function that produced it, and where it appears in the paper. The
printed table below is that register, filtered to the values the abstract
quotes.""",
    # 2
    """## 2. The decomposition

A screen's expected output is modelled as a supply term times a selection
term,

$$P(v) \\;\\propto\\; S_{\\mathrm{mutagen}}(v)\\;\\cdot\\;\\exp(\\beta\\,\\sigma(v))$$

with $S$ from codon tables and the mutagen's published six-class spectrum,
and $\\sigma$ from a chemistry-free atlas. The six substitution classes and
the four chemistries are below; note that a base-editor pair gives four of
the six classes a rate of exactly zero, which is the whole of the hard
reachability claim.""",
    # 3
    """## 3. Figures 1 to 5

Each cell runs the figure's own script, which reads only from
`data/derived/` and writes a PDF and a PNG into `paper/`. The scripts share
`paper/_style.py`, so the five figures use one palette and one type scale.""",
    # 4
    """## 4. The two audits

The register is checked against the manuscript in both directions.

`audit_paper_numbers.py` asks whether every registered value appears in
`paper.tex` where the register says it does. `audit_unregistered_numbers.py`
asks the reverse: whether every number written in the prose is accounted
for by a registered value.

The second one is the one that earns its keep. It is how six stale recovery
figures and four stale variant ranks were found before submission, all of
them left behind by a codon-table rebuild that moved the ALK universe from
1,259 to 1,427 variants.""",
    # 5
    """## 5. What the figures show

**Figure 1** states the model and shows that the variants a chemistry
cannot reach sit at the top of the resistance ranking, not the bottom.

**Figure 2** gives the four chemistries' six-class spectra and the coverage
tiers they fall into.

**Figure 3** is the central one. Panel A is the 24-variant resistance panel
against each chemistry; panel B is top-$K$ recovery; panel C is the
per-residue reachability inventory across all three kinases, where ENU and
Mutazyme II block no residue in any gene and a base-editor pair blocks
63%–66% of them.

**Figure 4** is the retrospective validation across three ALK screens.

**Figure 5** is the MET audit, the one pipeline whose screen reports its
hits at nucleotide resolution, so the mutagen's output is observed rather
than reconstructed.""",
]

CODE = [
    # after MD[1]
    """import json, subprocess, sys
from pathlib import Path

ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
sys.path.insert(0, str(ROOT))

register = json.loads((ROOT / "data/derived/paper_numbers.json").read_text())
print(f"{len(register)} registered values\\n")

for key, entry in sorted(register.items()):
    if any("abstract" in str(w).lower() for w in entry.get("where", [])):
        print(f"  {key:34s} {entry['formatted']:>22s}   {entry['code_fn']}")""",
    # after MD[2]
    """import yaml
import pandas as pd

spectra = yaml.safe_load((ROOT / "data/reference/spectra.yaml").read_text())["spectra"]
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

rows = {
    "ENU":         dict(spectra["ENU_Arnold2012_nominal"]["values"]),
    "Mutazyme II": dict(spectra["MutazymeII_DoiSupp2022"]["values"]),
    "CBE+ABE":     {"C>A": 0.0, "C>G": 0.0, "C>T": 0.5,
                    "T>A": 0.0, "T>C": 0.5, "T>G": 0.0},
    "Twist / PE":  {c: 1 / 6 for c in CLASSES},
}
table = pd.DataFrame(rows).T[CLASSES]
display(table.style.format("{:.3f}").background_gradient(cmap="Blues", axis=None))

blocked = {m: [c for c in CLASSES if v[c] == 0] for m, v in rows.items()}
for m, cs in blocked.items():
    print(f"{m:12s} supplies {6 - len(cs)} of 6 classes"
          + (f"; zero rate for {', '.join(cs)}" if cs else ""))""",
    # after MD[3]
    """for n in range(1, 6):
    script = ROOT / "paper" / f"figure_{n}.py"
    r = subprocess.run([sys.executable, str(script)],
                       capture_output=True, text=True, cwd=ROOT)
    status = "ok" if r.returncode == 0 else "FAILED"
    print(f"figure_{n}.py  {status}")
    if r.returncode:
        print(r.stderr[-800:])""",
    # display cell
    """from IPython.display import Image, display as show

for n in range(1, 6):
    png = ROOT / "paper" / f"figure_{n}.png"
    print(f"\\nFigure {n}")
    show(Image(filename=str(png), width=900))""",
    # after MD[4]
    """for audit in ["audit_paper_numbers.py", "audit_unregistered_numbers.py"]:
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / audit)],
                       capture_output=True, text=True, cwd=ROOT)
    print(f"--- {audit} (exit {r.returncode}) ---")
    print(r.stdout.strip() or r.stderr.strip())
    print()""",
]


def main() -> None:
    nb = nbf.v4.new_notebook()
    cells = [
        nbf.v4.new_markdown_cell(MD[0]),
        nbf.v4.new_markdown_cell(MD[1]),
        nbf.v4.new_code_cell(CODE[0]),
        nbf.v4.new_markdown_cell(MD[2]),
        nbf.v4.new_code_cell(CODE[1]),
        nbf.v4.new_markdown_cell(MD[3]),
        nbf.v4.new_code_cell(CODE[2]),
        nbf.v4.new_code_cell(CODE[3]),
        nbf.v4.new_markdown_cell(MD[4]),
        nbf.v4.new_code_cell(CODE[4]),
        nbf.v4.new_markdown_cell(MD[5]),
    ]
    nb["cells"] = cells
    nb["metadata"] = {
        "kernelspec": {"display_name": "Python 3", "language": "python",
                       "name": "python3"},
        "language_info": {"name": "python", "version": "3.13.16"},
    }
    OUT.parent.mkdir(exist_ok=True)
    nbf.write(nb, OUT)
    print(f"wrote {OUT.relative_to(ROOT)}  ({len(cells)} cells)")


if __name__ == "__main__":
    main()
