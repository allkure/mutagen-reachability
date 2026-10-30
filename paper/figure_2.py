"""Figure 2 v2 — mutagen substitution spectra used in preclinical practice.

Two panels only (the previous Panel B "current practice" table moves to
a LaTeX table in paper.tex per MR 2026-10-14 feedback).

Panel A: six-class spectra for ENU, Mutazyme II, CBE+ABE, Twist/PE.
Panel B: class coverage score (number of classes supplied at ≥ 1 %).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

NOTEBOOK_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = NOTEBOOK_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from paper._style import PALETTE, apply_rc
from src.spectrum import load_spectra

apply_rc()

SPECTRA = load_spectra(PROJECT_ROOT / "data" / "reference" / "spectra.yaml")
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

MUTAGEN_PANEL = [
    ("ENU",         {k: SPECTRA["ENU_Arnold2012_nominal"][k] for k in CLASSES},
     PALETTE["ENU"]),
    ("Mutazyme II", {k: SPECTRA["MutazymeII_DoiSupp2022"][k] for k in CLASSES},
     PALETTE["MutII"]),
    ("CBE + ABE",   {"C>A": 0.0, "C>G": 0.0, "C>T": 0.5, "T>A": 0.0,
                     "T>C": 0.5, "T>G": 0.0}, PALETTE["BE"]),
    ("Twist / PE",  {k: 1 / 6 for k in CLASSES}, PALETTE["chemfree"]),
]

fig = plt.figure(figsize=(19.0, 10.5))
gs  = fig.add_gridspec(1, 2, width_ratios=[1.9, 1.0], wspace=0.22,
                       left=0.07, right=0.98, top=0.94, bottom=0.20)
axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])

# ===========================================================================
# A — spectrum bar chart
# ===========================================================================
axA.set_title("A   Six-class substitution spectra", loc="left",
              fontweight="bold", pad=8)

n_mut, n_cls = len(MUTAGEN_PANEL), len(CLASSES)
bar_w = 0.19
x_pos = np.arange(n_cls)
for i, (name, spec, colour) in enumerate(MUTAGEN_PANEL):
    offs = (i - (n_mut - 1) / 2) * bar_w
    vals = [spec[c] for c in CLASSES]
    axA.bar(x_pos + offs, vals, width=bar_w, label=name,
            color=colour, edgecolor="#1a1a1a", linewidth=0.5)

axA.set_xticks(x_pos)
axA.set_xticklabels(CLASSES, fontsize=24, fontweight="bold")
axA.set_ylabel("Probability mass per class")
axA.set_ylim(0, 0.65)
# "1/6 uniform" baseline — thicker dotted line + inline label placed in the
# gap between C>G and C>T, well clear of y-axis tick labels and of bars.
axA.axhline(1 / 6, color=PALETTE["neutral"], lw=1.4, ls=":")
axA.text(1.5, 1 / 6 + 0.014, "1/6  (uniform baseline)", fontsize=18,
         color=PALETTE["neutral"], ha="center", va="bottom",
         fontstyle="italic")

# Legend below the panel
axA.legend(loc="upper center", frameon=False, fontsize=24,
           ncol=4, columnspacing=2.0, bbox_to_anchor=(0.5, -0.22),
           handletextpad=0.8, handlelength=1.8)

# Light background band for the two base-editing classes
for cls_idx, cls in enumerate(CLASSES):
    if cls in ("C>T", "T>C"):
        axA.axvspan(cls_idx - 0.45, cls_idx + 0.45,
                    color=PALETTE["BE"], alpha=0.12, zorder=0)

# ===========================================================================
# B — class coverage
# ===========================================================================
axB.set_title("B   Class coverage (≥ 1 % supply)", loc="left",
              fontweight="bold", pad=8)

names    = [m[0] for m in MUTAGEN_PANEL]
colours  = [m[2] for m in MUTAGEN_PANEL]
coverage = [sum(1 for c in CLASSES if m[1][c] >= 0.01) for m in MUTAGEN_PANEL]
y_pos    = np.arange(len(names))

axB.barh(y_pos, coverage, color=colours, edgecolor="#1a1a1a", linewidth=1.2)
axB.set_yticks(y_pos)
axB.set_yticklabels(names, fontsize=24)
axB.invert_yaxis()
axB.set_xlabel("Classes supplied (of 6)")
axB.set_xlim(0, 7.5)
axB.set_xticks(range(7))
# "all 6" baseline — rotated label hanging off the dashed line to the right
axB.axvline(6, color=PALETTE["neutral"], lw=1.4, ls=":")
axB.text(6.25, 1.5, "all 6", fontsize=22,
         color=PALETTE["neutral"], ha="left", va="center", rotation=0)

# Numbers at bar end
for yi, c in enumerate(coverage):
    axB.text(c + 0.28, yi, str(c), fontsize=26,
             va="center", fontweight="bold")

out_png = PROJECT_ROOT / "paper" / "figure_2.png"
out_pdf = PROJECT_ROOT / "paper" / "figure_2.pdf"
fig.savefig(out_png, dpi=300, bbox_inches="tight")
fig.savefig(out_pdf,              bbox_inches="tight")
print(f"Wrote {out_png}")
