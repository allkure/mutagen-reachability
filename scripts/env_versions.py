#!/usr/bin/env python3
"""Emit the actual interpreter and package versions used for the analysis.

The manuscript previously listed hand-typed versions, one of which (pyyaml 7.0)
does not exist.  In a paper whose claim is auditability, the software stanza
must be generated, not remembered.  Writes a LaTeX-ready sentence plus a JSON
record for NUMBERS.md.
"""
from __future__ import annotations

import json
import platform
import sys
from pathlib import Path

import matplotlib
import numpy
import pandas
import scipy
import yaml

ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = ROOT / "data" / "derived" / "env_versions.json"
OUT_TEX = ROOT / "data" / "derived" / "env_versions.tex"

PKGS = {
    "pandas": pandas.__version__,
    "numpy": numpy.__version__,
    "scipy": scipy.__version__,
    "matplotlib": matplotlib.__version__,
    "pyyaml": yaml.__version__,
}


def main() -> None:
    rec = {
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "packages": PKGS,
    }
    OUT_JSON.write_text(json.dumps(rec, indent=2))

    pkg_list = ", ".join(f"{k} {v}" for k, v in list(PKGS.items())[:-1])
    sentence = (f"All analysis is in Python {rec['python']} with {pkg_list}, "
                f"and pyyaml {PKGS['pyyaml']}.")
    OUT_TEX.write_text(sentence + "\n")
    print(sentence)


if __name__ == "__main__":
    main()
