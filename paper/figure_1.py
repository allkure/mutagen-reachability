"""Figure 1 v3 — per MR 2026-10-14 pm feedback:

  - Panel A redesigned (no grey, no redundant "drug selection" box, cleaner flow)
  - Panel B has shaded background bands + dashed separators for the 3 groups
  - Panels C/D use log-log axes, generic axis labels (no author names),
    per-panel legends below each panel
  - Overall font bump (base 15)
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

NOTEBOOK_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = NOTEBOOK_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from paper._style import PALETTE, apply_rc

apply_rc()

# ---------------------------------------------------------------------------
kras_results = pd.read_csv(PROJECT_ROOT / "data" / "derived" / "feng_decomposition.csv")
alk_results  = pd.read_csv(PROJECT_ROOT / "data" / "derived" / "oh_decomposition_all_drugs.csv")

AWAD_PANEL = [
    ("R68L", "unreachable", 3.47, 4.34),
    ("V8M",  "unreachable", 1.02, 1.10),
    ("V8L",  "multi",       3.08, 1.38),
    ("R68S", "multi",       1.26, 2.67),
    ("H95Q", "multi",       0.33, 2.36),
    ("C12F", "forced",      2.76, 2.53),
    ("Y96C", "forced",      3.32, 3.30),
    ("Y96D", "forced",      2.32, 1.82),
    ("Y96H", "forced",      2.71, 2.10),
    ("H95D", "forced",     -0.05, 2.50),
    ("H95R", "forced",      0.42, 3.77),
    ("R68M", "forced",      0.92, 0.90),
]

CAT_COLORS = {"unreachable": PALETTE["unreachable"],
              "multi":       PALETTE["undercounted"],
              "forced":      PALETTE["accessible"]}
CAT_LABEL  = {"unreachable": "SNV-unreachable (2+ nt from WT codon)",
              "multi":       "Multi-class SNV (ambiguous chemistry)",
              "forced":      "Class-forced SNV (chemistry predicts rank)"}

# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(18.0, 16.5))
gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.26,
                      left=0.07, right=0.98, top=0.96, bottom=0.08)

axA = fig.add_subplot(gs[0, 0])
axB = fig.add_subplot(gs[0, 1])
axC = fig.add_subplot(gs[1, 0])
axD = fig.add_subplot(gs[1, 1])

# ===========================================================================
# A — schematic (redesigned: cream library at top, two routes, equation)
# ===========================================================================
axA.set_xlim(0, 10)
axA.set_ylim(-0.3, 10.2)
axA.axis("off")
axA.set_title("A   Supply × selection decomposition", loc="left",
              fontweight="bold", pad=14)

# ---- modern geometry ----
# Common box style: softly rounded, subtle shadow via a lighter rectangle
# drawn below each visible rectangle
def _shadow(x, y, w, h, offset=0.08):
    axA.add_patch(mpatches.FancyBboxPatch(
        (x + offset, y - offset), w, h,
        boxstyle="round,pad=0.12,rounding_size=0.25",
        facecolor="#E6E4DE", edgecolor="none", zorder=0, alpha=0.65))

def _box(x, y, w, h, fc, ec=None):
    _shadow(x, y, w, h)
    axA.add_patch(mpatches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.12,rounding_size=0.25",
        facecolor=fc, edgecolor=ec or "#1a1a1a",
        linewidth=1.4, alpha=0.85, zorder=1))

# ---- Top: library of variants ----
_box(2.0, 8.1, 6.0, 1.3, PALETTE["cream"], ec=PALETTE["neutral"])
axA.text(5.0, 8.75, "Library of variants",
         ha="center", va="center", fontsize=26, fontweight="bold")

# ---- Two routes ----
_box(0.2, 4.9, 4.3, 2.0, PALETTE["ENU"])
axA.text(2.35, 5.90, "Chemistry-biased\nmutagen",
         ha="center", va="center", fontsize=26, fontweight="bold")

_box(5.5, 4.9, 4.3, 2.0, PALETTE["chemfree"])
axA.text(7.65, 5.90, "Chemistry-free\nlibrary",
         ha="center", va="center", fontsize=26, fontweight="bold")

# ---- Modern arrows: thick lines with filled triangle heads ----
arrow_kw = dict(arrowstyle="-|>,head_width=0.5,head_length=0.9",
                color=PALETTE["neutral"], linewidth=3.0,
                zorder=4, mutation_scale=18,
                connectionstyle="arc3,rad=0")
axA.annotate("", xy=(2.9, 6.95), xytext=(4.2, 8.1), arrowprops=arrow_kw)
axA.annotate("", xy=(7.1, 6.95), xytext=(5.8, 8.1), arrowprops=arrow_kw)
axA.annotate("", xy=(4.1, 3.35), xytext=(2.9, 4.85), arrowprops=arrow_kw)
axA.annotate("", xy=(5.9, 3.35), xytext=(7.1, 4.85), arrowprops=arrow_kw)

# ---- Equation band ----
_box(1.4, 1.9, 7.2, 1.4, "#FFFFFF", ec=PALETTE["neutral"])
# Must match Eq. 1 in the Methods exactly: proportionality, not equality,
# and the exponential of beta*sigma. sigma is a log-scale score that goes
# negative, so a bare product would give negative probabilities.
axA.text(5.0, 2.6,
         r"$P(v) \;\propto\; S_{mutagen}(v) \;\cdot\; \exp(\beta\,\sigma(v))$",
         ha="center", va="center", fontsize=27, fontweight="bold")

# (Panel B category legend is placed directly under Panel B — see Panel B
# code — so Panel A stays free of overlapping legends.)

# ===========================================================================
# B — KRAS G12C resistance panel categorised (composite: Awad 2021 clinical
# cohort + Tanaka 2021 + Awad 2021 in vitro deep mutational scan — see
# paper Supp Table S1 for per-variant provenance)
# ===========================================================================
axB.set_title("B   KRAS G12C resistance panel", loc="left",
              fontweight="bold", pad=10)

# Group backgrounds (shaded) + dashed separators
BAND_ALPHA = 0.14
axB.axvspan(-0.3, 1.25, facecolor=PALETTE["unreachable"], alpha=BAND_ALPHA, zorder=0)
axB.axvspan(1.25, 2.75, facecolor=PALETTE["undercounted"], alpha=BAND_ALPHA, zorder=0)
axB.axvspan(2.75, 4.6,  facecolor=PALETTE["accessible"],   alpha=BAND_ALPHA, zorder=0)
axB.axvline(1.25, color=PALETTE["neutral"], lw=1.4, ls=(0, (6, 4)), zorder=1)
axB.axvline(2.75, color=PALETTE["neutral"], lw=1.4, ls=(0, (6, 4)), zorder=1)

MANUAL_XY = {
    "R68L":   (0.40, 4.34),
    "V8M":    (0.60, 1.10),
    "V8L":    (1.75, 3.08),
    "R68S":   (1.85, 2.67),
    "H95Q":   (2.40, 2.00),
    "H95R":   (3.05, 3.77),
    "Y96C":   (3.55, 3.30),
    "C12F":   (4.00, 2.90),
    "Y96H":   (3.05, 2.71),
    "H95D":   (3.75, 2.45),
    "Y96D":   (3.30, 1.95),
    "R68M":   (3.60, 0.95),
}
LABEL_OFFSET = {
    "R68L":  (0.12, 0.00),  "V8M":   (0.12, 0.10),
    "V8L":   (0.12, 0.12),  "R68S":  (0.12, 0.12),
    "H95Q":  (0.12, 0.12),  "H95R":  (-0.45, 0.18),
    "Y96C":  (0.12, 0.10),  "C12F":  (0.12, 0.10),
    "Y96H":  (-0.45, 0.10), "H95D":  (0.12, 0.10),
    "Y96D":  (0.12, 0.10),  "R68M":  (0.12, 0.10),
}
for v, cat, sz, az in AWAD_PANEL:
    xj, y = MANUAL_XY[v]
    axB.scatter(xj, y, s=360, c=CAT_COLORS[cat],
                edgecolors="#1a1a1a", linewidths=1.4, zorder=3)
    ox, oy = LABEL_OFFSET[v]
    axB.annotate(v, (xj, y), xytext=(xj + ox, y + oy),
                 fontsize=24)

axB.set_xticks([0.5, 2.0, 3.5])
axB.set_xticklabels(["Unreachable\n(n = 2)",
                     "Multi-class\n(n = 3)",
                     "Class-forced\n(n = 7)"], fontsize=24)
axB.set_xlim(-0.3, 4.6)
axB.set_ylim(-0.6, 5.6)
axB.set_ylabel("Max chemistry-free Z-score across drugs")
axB.axhline(0, color=PALETTE["neutral"], lw=0.5, zorder=0)

# R68L callout — compact
axB.annotate("rank-1 variant\nunder both drugs",
             xy=(0.42, 4.34), xytext=(0.95, 4.95),
             fontsize=20, ha="left",
             arrowprops=dict(arrowstyle="->", color=PALETTE["neutral"], lw=1.4))

# Panel B category legend — horizontal row directly under Panel B
legend_handles_B = [
    mpatches.Patch(facecolor=PALETTE["unreachable"], edgecolor="#1a1a1a",
                   linewidth=0.9, label="SNV-unreachable"),
    mpatches.Patch(facecolor=PALETTE["undercounted"], edgecolor="#1a1a1a",
                   linewidth=0.9, label="Multi-class SNV"),
    mpatches.Patch(facecolor=PALETTE["accessible"], edgecolor="#1a1a1a",
                   linewidth=0.9, label="Class-forced SNV"),
]
axB.legend(handles=legend_handles_B, loc="upper center",
           bbox_to_anchor=(0.5, -0.19), ncol=3, frameon=False,
           fontsize=26, columnspacing=2.0, handletextpad=0.8,
           handlelength=1.8, handleheight=1.6)

# ===========================================================================
# C — KRAS scatter (ENU × adagrasib), log-log
# ===========================================================================
sub = kras_results[(kras_results["mutagen"] == "ENU") &
                   (kras_results["drug"] == "adagrasib")].copy()

axC.set_title("C   KRAS G12C · adagrasib · ENU",
              loc="left", fontweight="bold", pad=10)

# Diagonal (log-log)
axC.plot([1, 2000], [1, 2000], ls="--", color=PALETTE["neutral"], lw=1.0)
axC.scatter(sub["rank_chemfree"], sub["rank_simulated"],
            s=22, c=PALETTE["grid"], alpha=0.65, edgecolors="none")

AWAD_FORCED = ["C12F", "Y96C", "Y96D", "Y96H", "H95D", "H95R", "R68M"]
for v in AWAD_FORCED:
    r = sub[sub["aa_change"] == v]
    if len(r) == 0:
        continue
    rf = int(r["rank_chemfree"].iloc[0])
    rs = int(r["rank_simulated"].iloc[0])
    axC.scatter(rf, rs, s=200, c=PALETTE["clinical"],
                edgecolors="#1a1a1a", linewidths=1.2, zorder=4)

# Per-variant labels — explicit (x, y) in log coords, separated so the
# bottom-left cluster and the two top outliers each get breathing room.
LABEL_C = {
    # (xytext_x, xytext_y, ha)
    "H95R": (7,    1.1,  "left"),
    "Y96C": (7,    2.2,  "left"),
    "Y96H": (33,   6.2,  "left"),
    "C12F": (20,   30,   "left"),
    "Y96D": (55,   310,  "left"),
    "R68M": (170,  420,  "left"),
    "H95D": (28,   1050, "left"),
}
for v, (tx, ty, ha) in LABEL_C.items():
    r = sub[sub["aa_change"] == v]
    if len(r) == 0:
        continue
    rf = int(r["rank_chemfree"].iloc[0])
    rs = int(r["rank_simulated"].iloc[0])
    axC.annotate(v, xy=(rf, rs), xytext=(tx, ty),
                 fontsize=22, ha=ha,
                 arrowprops=dict(arrowstyle="-",
                                 color=PALETTE["neutral_soft"], lw=1.0))

axC.set_xscale("log")
axC.set_yscale("log")
axC.set_xlim(0.8, 1500)
axC.set_ylim(0.8, 1500)
axC.set_xlabel("Chemistry-free rank")
axC.set_ylabel("ENU-simulated preclinical rank")

# Legend below panel
handle_C = plt.Line2D([0], [0], marker="o", linestyle="",
                       markerfacecolor=PALETTE["clinical"],
                       markeredgecolor="#1a1a1a", markersize=22,
                       label="KRAS G12C resistance-panel variant (class-forced)")
handle_C_diag = plt.Line2D([0, 1], [0, 1], ls="--",
                            color=PALETTE["neutral"], lw=1.6,
                            label="no shift (diagonal)")
axC.legend(handles=[handle_C, handle_C_diag], loc="upper center",
           bbox_to_anchor=(0.5, -0.20), ncol=1, frameon=False, fontsize=26)

# ===========================================================================
# D — ALK scatter (MutII × alectinib), log-log
# ===========================================================================
sub_alk = alk_results[(alk_results["mutagen"] == "MutazymeII") &
                      (alk_results["drug"] == "alectinib")].copy()

axD.set_title("D   ALK · alectinib · Mutazyme II",
              loc="left", fontweight="bold", pad=10)

axD.plot([1, 2000], [1, 2000], ls="--", color=PALETTE["neutral"], lw=1.0)
axD.scatter(sub_alk["rank_chemfree"], sub_alk["rank_simulated"],
            s=22, c=PALETTE["grid"], alpha=0.65, edgecolors="none")

ALK_HILITE = ["V1180L", "L1196M", "G1269A", "F1174C",
              "I1171T", "I1171N", "C1156Y", "L1152R"]
for v in ALK_HILITE:
    r = sub_alk[sub_alk["aa_change"] == v]
    if len(r) == 0:
        continue
    rf = int(r["rank_chemfree"].iloc[0])
    rs = int(r["rank_simulated"].iloc[0])
    axD.scatter(rf, rs, s=200, c=PALETTE["clinical"],
                edgecolors="#1a1a1a", linewidths=1.2, zorder=4)

# Explicit label positions in log coords; all labels INSIDE the plot area.
LABEL_D = {
    "V1180L": (1.5,  2.3,  "left"),
    "L1196M": (6.5,  2.0,  "left"),
    "G1269A": (15,   13,   "left"),
    "I1171T": (90,   35,   "left"),
    "F1174C": (50,   90,   "left"),
    "C1156Y": (35,   60,   "left"),
    "I1171N": (300,  160,  "left"),
    "L1152R": (450,  900,  "left"),
}
for v, (tx, ty, ha) in LABEL_D.items():
    r = sub_alk[sub_alk["aa_change"] == v]
    if len(r) == 0:
        continue
    rf = int(r["rank_chemfree"].iloc[0])
    rs = int(r["rank_simulated"].iloc[0])
    axD.annotate(v, xy=(rf, rs), xytext=(tx, ty),
                 fontsize=22, ha=ha,
                 arrowprops=dict(arrowstyle="-",
                                 color=PALETTE["neutral_soft"], lw=1.0))

axD.set_xscale("log")
axD.set_yscale("log")
axD.set_xlim(0.8, 1500)
axD.set_ylim(0.8, 1500)
axD.set_xlabel("Chemistry-free rank")
axD.set_ylabel("Mutazyme-simulated preclinical rank")

handle_D = plt.Line2D([0], [0], marker="o", linestyle="",
                       markerfacecolor=PALETTE["clinical"],
                       markeredgecolor="#1a1a1a", markersize=22,
                       label="ALK resistance-panel variant (class-forced)")
axD.legend(handles=[handle_D, handle_C_diag], loc="upper center",
           bbox_to_anchor=(0.5, -0.20), ncol=1, frameon=False, fontsize=26)

out_png = PROJECT_ROOT / "paper" / "figure_1.png"
out_pdf = PROJECT_ROOT / "paper" / "figure_1.pdf"
fig.savefig(out_png, dpi=300, bbox_inches="tight")
fig.savefig(out_pdf,              bbox_inches="tight")
print(f"Wrote {out_png}")
