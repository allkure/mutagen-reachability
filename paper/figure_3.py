"""Figure 3 v4 — Two panels.

Panel A (left):  Clinical variant × mutagen chemistry reachability matrix.
                 Rows: composite KRAS G12C resistance panel (top 12 — Awad
                       2021 clinical + in vitro DMS, Tanaka 2021 clinical)
                       and ALK resistance panel (bottom 12 — Doi 2022 +
                       Takei 2025). Per-variant provenance in Table S1.
                 Cols: four mutagen chemistries.

Panel B (right): Top-K recovery small-multiples. One subpanel per
                 (gene × drug) combination; one curve per mutagen chemistry.
                 Y-axis: fraction of the chemistry-free top-K recovered by
                 the chemistry-biased screen's top-K. X-axis: K on log scale,
                 1 to full saturation-atlas size (1,205 KRAS or 1,259 ALK).
                 Twist/PE (chemistry-free) sits flat at y=1 as the ideal
                 reference. CBE+ABE plateaus at the fraction of the atlas
                 it can supply, directly showing the coverage ceiling.

Style notes (per MR 2026-10-14 feedback):
  - Pastel colour-blind-friendly palette (locked across all 3 figures).
  - No figure suptitle; no coloured text.
  - Larger fonts throughout.
  - 5th distinct graphical primitive across the paper (schematic, scatter,
    bar, matrix, curve) — no visual element is reused.
"""

from __future__ import annotations

import sys
from pathlib import Path

import json

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import yaml

NOTEBOOK_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = NOTEBOOK_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from paper._style import PALETTE, apply_rc
from src.harmonized_cohort import classify_amino_acid_change
from src.spectrum import load_spectra

apply_rc()
SPECTRA = load_spectra(PROJECT_ROOT / "data" / "reference" / "spectra.yaml")

MUTAGENS = {
    "ENU":           SPECTRA["ENU_Arnold2012_nominal"].values,
    "Mutazyme\nII":  SPECTRA["MutazymeII_DoiSupp2022"].values,
    "CBE\n+ ABE":    {"C>A": 0.0, "C>G": 0.0, "C>T": 0.5, "T>A": 0.0,
                      "T>C": 0.5, "T>G": 0.0},
    "Twist\n/ PE":   {k: 1 / 6 for k in ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]},
}

with open(PROJECT_ROOT / "data" / "reference" / "codons_kras.yaml") as f:
    KRAS_CODONS = {int(r): v["codon"]
                   for r, v in yaml.safe_load(f)["kras_4a_nm_033360_4"]["codons"].items()}
KRAS_CODONS[12] = "TGT"

with open(PROJECT_ROOT / "data" / "reference" / "codons.yaml") as f:
    ALK_CODONS = {int(r): v["codon"]
                  for r, v in yaml.safe_load(f)["alk_nm_004304_5"]["codons"].items()
                  if "codon" in v}

AWAD_PANEL = [
    ("C12F", "KRAS G12C", KRAS_CODONS[12]),
    ("R68L", "KRAS G12C", KRAS_CODONS[68]),
    ("Y96C", "KRAS G12C", KRAS_CODONS[96]),
    ("V8L",  "KRAS G12C", KRAS_CODONS[8]),
    ("V8M",  "KRAS G12C", KRAS_CODONS[8]),
    ("Y96D", "KRAS G12C", KRAS_CODONS[96]),
    ("Y96H", "KRAS G12C", KRAS_CODONS[96]),
    ("H95D", "KRAS G12C", KRAS_CODONS[95]),
    ("H95Q", "KRAS G12C", KRAS_CODONS[95]),
    ("H95R", "KRAS G12C", KRAS_CODONS[95]),
    ("R68S", "KRAS G12C", KRAS_CODONS[68]),
    ("R68M", "KRAS G12C", KRAS_CODONS[68]),
]
ALK_PANEL = [
    ("V1180L", "ALK", ALK_CODONS.get(1180)),
    ("L1196M", "ALK", ALK_CODONS.get(1196)),
    ("I1171N", "ALK", ALK_CODONS.get(1171)),
    ("I1171T", "ALK", ALK_CODONS.get(1171)),
    ("G1202R", "ALK", ALK_CODONS.get(1202)),
    ("G1269A", "ALK", ALK_CODONS.get(1269)),
    ("F1174C", "ALK", ALK_CODONS.get(1174)),
    ("F1174L", "ALK", ALK_CODONS.get(1174)),
    ("L1196Q", "ALK", ALK_CODONS.get(1196)),
    ("L1198F", "ALK", ALK_CODONS.get(1198)),
    ("C1156Y", "ALK", ALK_CODONS.get(1156)),
    ("L1152R", "ALK", ALK_CODONS.get(1152)),
]

UNDERCOUNT = 0.10

def cell_status(wt_aa, mut_aa, wt_codon, spec):
    if wt_codon is None:
        return "no_codon", None, 0.0
    info = classify_amino_acid_change(wt_aa, mut_aa, wt_codon=wt_codon)
    if not info["reachable_1nt"]:
        return "unreachable", None, 0.0
    classes = info["classes"]
    max_sup = max(spec[c] for c in classes)
    return ("accessible" if max_sup >= UNDERCOUNT else "undercounted",
            classes, max_sup)


rows = []
for v, gene, wt_codon in AWAD_PANEL + ALK_PANEL:
    wt_aa = v[0]; mut_aa = v[-1]
    r = {"variant": v, "gene": gene}
    for mn, spec in MUTAGENS.items():
        s, _, sup = cell_status(wt_aa, mut_aa, wt_codon, spec)
        r[mn] = {"status": s, "sup": sup}
    rows.append(r)

# ===========================================================================
# Top-K recovery curves for Panel B
# ===========================================================================
recov_df = pd.read_csv(PROJECT_ROOT / "data" / "derived" / "topk_recovery.csv")

# ---------------------------------------------------------------------------
# Figure layout
# ---------------------------------------------------------------------------
STATUS_COLOURS = {
    "accessible":    PALETTE["accessible"],
    "undercounted":  PALETTE["undercounted"],
    "unreachable":   PALETTE["unreachable"],
    "no_codon":      PALETTE["grid"],
}
# Palette for the four mutagens in Panel B (matches Fig 2 bar colours)
MUT_COLOUR = {
    "ENU":   PALETTE["ENU"],
    "MutII": PALETTE["MutII"],
    "BE":    PALETTE["BE"],
    "Twist": PALETTE["chemfree"],
}
MUT_LABEL = {
    "ENU":   "ENU",
    "MutII": "Mutazyme II",
    "BE":    "CBE + ABE",
    "Twist": "Twist / PE",
}

fig = plt.figure(figsize=(30.0, 25.0))
outer = fig.add_gridspec(2, 2, width_ratios=[0.88, 1.35],
                         height_ratios=[1.0, 0.40],
                         wspace=0.26, hspace=0.34,
                         left=0.055, right=0.985, top=0.94, bottom=0.055)
axA = fig.add_subplot(outer[0, 0])
inner = outer[0, 1].subgridspec(2, 3, wspace=0.48, hspace=0.48)
# Panel C spans the full width under A and B. Column widths are proportional
# to each gene's residue count so one residue is the same width throughout.
_trk = json.loads((PROJECT_ROOT / "data" / "derived"
                   / "reachability_track.json").read_text())
C_GENES = [("ALK", "alectinib"), ("KRAS G12C", "adagrasib"),
           ("MET", "NVP-BVU972")]
_cwid = [_trk["by_gene"][g]["n_residues"] for g, _ in C_GENES]
cgrid = outer[1, :].subgridspec(1, len(C_GENES) + 1,
                               width_ratios=_cwid + [max(_cwid) * 0.30],
                               wspace=0.045)
axC = [fig.add_subplot(cgrid[0, i]) for i in range(len(C_GENES))]
axC_cbar = fig.add_subplot(cgrid[0, len(C_GENES)])
# Row 1: KRAS G12C (sotorasib, adagrasib) + legend cell
# Row 2: ALK       (alectinib, lorlatinib, zotizalkib)
axB_sot  = fig.add_subplot(inner[0, 0])
axB_ada  = fig.add_subplot(inner[0, 1])
axB_leg  = fig.add_subplot(inner[0, 2])
axB_alec = fig.add_subplot(inner[1, 0])
axB_lor  = fig.add_subplot(inner[1, 1])
axB_zot  = fig.add_subplot(inner[1, 2])
B_AXES = {
    ("KRAS G12C","sotorasib"):  axB_sot,
    ("KRAS G12C","adagrasib"):  axB_ada,
    ("ALK","alectinib"):        axB_alec,
    ("ALK","lorlatinib"):       axB_lor,
    ("ALK","zotizalkib"):       axB_zot,
}

# ===========================================================================
# A — reachability matrix
# ===========================================================================
# Title drawn as fig.text later so it aligns exactly with Panel B's title
# (set_title + custom pad is fragile across figure sizes).

mutagen_names = list(MUTAGENS.keys())
n_cols = len(mutagen_names)
n_rows = len(rows)

CELL_W = 3.0
CELL_H = 1.05
GAP    = 0.07
matrix_right = n_cols * CELL_W
total_h      = n_rows * (CELL_H + GAP)

for ri, row in enumerate(reversed(rows)):
    y = ri * (CELL_H + GAP)
    for ci, mut in enumerate(mutagen_names):
        x = ci * CELL_W
        cell = row[mut]
        colour = STATUS_COLOURS[cell["status"]]
        axA.add_patch(mpatches.Rectangle((x, y), CELL_W, CELL_H,
                                         facecolor=colour,
                                         edgecolor="white", linewidth=1.0))
        if cell["status"] in ("accessible", "undercounted"):
            pct = cell["sup"] * 100
            axA.text(x + CELL_W / 2, y + CELL_H / 2, f"{pct:.1f}%",
                     ha="center", va="center",
                     fontsize=24, color="#1a1a1a")
        else:
            axA.text(x + CELL_W / 2, y + CELL_H / 2, "✗",
                     ha="center", va="center",
                     fontsize=34, color="#1a1a1a")

# Variant labels (left side)
for ri, row in enumerate(reversed(rows)):
    y = ri * (CELL_H + GAP) + CELL_H / 2
    axA.text(-0.30, y, row["variant"],
             ha="right", va="center", fontsize=24)

# Column labels (top)
y_top = total_h + 0.5
for ci, mut in enumerate(mutagen_names):
    axA.text(ci * CELL_W + CELL_W / 2, y_top, mut,
             ha="center", va="bottom", fontsize=26, fontweight="bold")

# Divider between KRAS G12C (12 top) and ALK (12 bottom)
y_divider = (n_rows - 12) * (CELL_H + GAP) - GAP / 2
axA.plot([-0.4, matrix_right + 0.15], [y_divider, y_divider],
         color=PALETTE["neutral"], lw=1.4)

# Right-side brackets
BRACKET_X  = matrix_right + 0.30
CAP_WIDTH  = 0.35
BRACKET_LW = 2.5
_br_kw = dict(color=PALETTE["neutral"], lw=BRACKET_LW,
              solid_capstyle="round", clip_on=False)

y_alk_bot = 0.08
y_alk_top = y_divider - GAP / 2 - 0.05
axA.plot([BRACKET_X, BRACKET_X], [y_alk_bot, y_alk_top], **_br_kw)
axA.plot([BRACKET_X - CAP_WIDTH, BRACKET_X], [y_alk_bot, y_alk_bot], **_br_kw)
axA.plot([BRACKET_X - CAP_WIDTH, BRACKET_X], [y_alk_top, y_alk_top], **_br_kw)

y_kras_bot = y_divider + GAP / 2 + 0.05
y_kras_top = total_h - 0.08
axA.plot([BRACKET_X, BRACKET_X], [y_kras_bot, y_kras_top], **_br_kw)
axA.plot([BRACKET_X - CAP_WIDTH, BRACKET_X], [y_kras_bot, y_kras_bot], **_br_kw)
axA.plot([BRACKET_X - CAP_WIDTH, BRACKET_X], [y_kras_top, y_kras_top], **_br_kw)

# Gene group text labels: moved OUT of Panel A — they now live in the
# margin strip between A and B (serving both panels), to avoid duplication.
LABEL_X = BRACKET_X + 0.35

axA.set_xlim(-3.9, LABEL_X + 0.3)
axA.set_ylim(-0.6, total_h + 1.6)
axA.axis("off")

# Legend for Panel A — only the 3 cell-category colours
legend_items = [
    mpatches.Patch(facecolor=PALETTE["accessible"], edgecolor="#1a1a1a",
                   linewidth=0.8,
                   label="Chemistry-accessible (class supply ≥ 10 %)"),
    mpatches.Patch(facecolor=PALETTE["undercounted"], edgecolor="#1a1a1a",
                   linewidth=0.8,
                   label="Chemistry-under-counted (class supply < 10 %)"),
    mpatches.Patch(facecolor=PALETTE["unreachable"], edgecolor="#1a1a1a",
                   linewidth=0.8,
                   label="SNV-unreachable (2+ nt from WT codon)"),
]
axA.legend(handles=legend_items, loc="upper left",
           bbox_to_anchor=(0.0, -0.012), ncol=1, frameon=False, fontsize=23,
           handlelength=1.7, handleheight=1.5, labelspacing=0.45)

# ===========================================================================
# B — Top-K recovery, small multiples (one subpanel per gene × drug)
# ===========================================================================
B_MUTAGENS = ["Twist", "ENU", "MutII", "BE"]  # Twist first (reference)

def _smooth(y, w=5):
    """Rolling-median smoothing to tame small-K noise without erasing it."""
    if w <= 1: return y
    import numpy as _np
    pad = w // 2
    yp = _np.pad(y, (pad, pad), mode="edge")
    out = _np.empty_like(y, dtype=float)
    for i in range(len(y)):
        out[i] = _np.median(yp[i:i + w])
    return out

for (gene, drug), ax in B_AXES.items():
    sub = recov_df[(recov_df["gene"] == gene) & (recov_df["drug"] == drug)]
    N = int(sub["N"].iloc[0])
    for mg in B_MUTAGENS:
        s = sub[sub["mutagen"] == mg].sort_values("K")
        colour = MUT_COLOUR[mg]
        K   = s["K"].values
        rec = _smooth(s["recovery"].values, w=5) if mg != "Twist" else \
              s["recovery"].values
        # Wilson 95% CI band (R9) — smooth the same way, then shade
        if mg != "Twist" and "ci_lo" in s.columns:
            lo = _smooth(s["ci_lo"].values, w=5)
            hi = _smooth(s["ci_hi"].values, w=5)
            ax.fill_between(K, lo, hi, color=colour, alpha=0.14,
                            linewidth=0, zorder=1)
        lw     = 2.8 if mg != "Twist" else 2.0
        ls     = ":" if mg == "Twist" else "-"
        alpha  = 0.75 if mg == "Twist" else 0.95
        ax.plot(K, rec, color=colour, lw=lw, ls=ls, alpha=alpha,
                zorder=3 if mg != "Twist" else 2)

    # Dashed reference at y=1 (ideal)
    ax.axhline(1.0, color=PALETTE["neutral_soft"], lw=0.9, ls="--", zorder=1)

    # Axes cosmetics
    ax.set_xscale("log")
    ax.set_xlim(1, N)
    ax.set_ylim(-0.03, 1.08)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", ".25", ".5", ".75", "1"], fontsize=22)
    ax.set_xticks([1, 10, 100, 1000])
    ax.set_xticklabels(["1", "10", "100", "1,000"], fontsize=22)
    ax.tick_params(axis="both", which="major", length=6, width=1.4)
    ax.tick_params(axis="x", which="minor", length=4, width=0.8,
                   color=PALETTE["grid"])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_linewidth(1.4)
    ax.spines["bottom"].set_linewidth(1.4)
    ax.grid(which="major", axis="y", color=PALETTE["grid"], lw=0.8,
            ls=":", alpha=0.6, zorder=0)

    # Subpanel title (drug name)
    ax.set_title(drug, fontsize=26, fontweight="bold", pad=10, loc="center")

# ---- Shared panel labels
# Both panel titles aligned to the same y so A and B feel like siblings.
# Shortened to avoid A overflowing into B at this font size.
TITLE_Y = 0.955
_apos   = axA.get_position()
_bpos   = axB_sot.get_position()
fig.text(_apos.x0, TITLE_Y,
         "A   Resistance panel × chemistry reachability",
         ha="left", va="bottom", fontsize=32, fontweight="bold",
         color="#1a1a1a")
fig.text(_bpos.x0 - 0.025, TITLE_Y,
         "B   Top-K recovery  (biased vs chemistry-free)",
         ha="left", va="bottom", fontsize=32, fontweight="bold",
         color="#1a1a1a")

# X-axis label — one, centred under the ALK (bottom) row
axB_lor.set_xlabel("Top-K rank threshold  (chemistry-biased screen)",
                   fontsize=26, labelpad=18)

# Gene group labels — one big label per group in the clean margin strip
# between Panel A and Panel B, spanning BOTH rows of each gene group so a
# single visual anchor serves both panels. No y-axis labels compete for
# this space now, so the labels can breathe.
_sot_pos  = axB_sot.get_position()
_alec_pos = axB_alec.get_position()
# Place labels close to the matrix's right edge (clear of Panel B's y-ticks)
_margin_x = axA.get_position().x1 + 0.025
_kras_y   = (_sot_pos.y1 + _sot_pos.y0) / 2
_alk_y    = (_alec_pos.y1 + _alec_pos.y0) / 2
fig.text(_margin_x, _kras_y, "KRAS\nG12C",
         ha="left", va="center", fontsize=36, fontweight="bold",
         color="#1a1a1a", linespacing=1.05)
fig.text(_margin_x, _alk_y, "ALK",
         ha="left", va="center", fontsize=36, fontweight="bold",
         color="#1a1a1a")

# ---- Legend cell (top-right of Panel B grid)
axB_leg.axis("off")
legend_lines = [
    plt.Line2D([0],[0], color=MUT_COLOUR["ENU"],   lw=3.4,
               label="ENU"),
    plt.Line2D([0],[0], color=MUT_COLOUR["MutII"], lw=3.4,
               label="Mutazyme II"),
    plt.Line2D([0],[0], color=MUT_COLOUR["BE"],    lw=3.4,
               label="CBE + ABE"),
    plt.Line2D([0],[0], color=MUT_COLOUR["Twist"], lw=2.4, ls=":",
               label="Twist / PE  (reference, y = 1)"),
]
leg = axB_leg.legend(handles=legend_lines, loc="center",
                     frameon=False, fontsize=24,
                     handlelength=2.6, handletextpad=0.9, labelspacing=1.2,
                     borderpad=0.5, title="Mutagen chemistry",
                     title_fontsize=26, alignment="left")
leg.get_title().set_fontweight("bold")

# ===========================================================================
# C — per-residue reachability inventory, all three kinases
# ===========================================================================
# For each residue, the substitution the chemistry-free atlas ranks highest
# there, shaded by the one-step supply each chemistry gives that substitution.
# Scope: class-forced substitutions only, so a chemistry is blocked here only
# where its spectrum gives that class zero. Variants needing two or more
# nucleotide changes lie outside this view by construction.
C_LANES = ["Twist/PE", "Mutazyme II", "ENU", "CBE+ABE"]
C_CMAP = plt.get_cmap("Blues")
C_SMAX = 0.5                      # supply at which the ramp saturates
C_BLOCKED = PALETTE["unreachable"]

_blocked_share = {}

for ax, (gene, drug) in zip(axC, C_GENES):
    G = _trk["by_gene"][gene]
    res = sorted({r["residue"] for r in G["track"]})
    index = {r: i for i, r in enumerate(res)}
    rgba = np.zeros((len(C_LANES), len(res), 4), dtype=float)

    for li, lane in enumerate(C_LANES):
        rows_l = [r for r in G["track"]
                  if r["mutagen"] == lane and r["drug"] == drug]
        assert len(rows_l) == len(res), (gene, lane, len(rows_l), len(res))
        for r in rows_l:
            ci = index[r["residue"]]
            if r["S"] > 0:
                rgba[li, ci] = C_CMAP(min(r["S"] / C_SMAX, 1.0))
            else:
                rgba[li, ci] = mcolors.to_rgba(C_BLOCKED)
        if lane == "CBE+ABE":
            n_out = sum(1 for r in rows_l if r["S"] <= 0)
            _blocked_share[gene] = (n_out, len(res))

    ax.imshow(rgba, aspect="auto", interpolation="nearest",
              extent=(res[0] - 0.5, res[-1] + 0.5, len(C_LANES) - 0.5, -0.5))
    ax.set_yticks(range(len(C_LANES)))
    ax.set_yticklabels(C_LANES if ax is axC[0] else [""] * len(C_LANES),
                       fontsize=24)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", labelsize=22)
    for spine in ax.spines.values():
        spine.set_visible(False)

    n_out, n_res = _blocked_share[gene]
    ax.set_title(f"{gene}  ·  {drug}",
                 fontsize=26, fontweight="bold", pad=12, loc="left",
                 color="#1a1a1a")
    ax.set_xlabel(f"{n_out} of {n_res} residues blocked under CBE + ABE",
                  fontsize=21, labelpad=12, color=PALETTE["neutral"])

# Shared supply ramp, plus the one colour that is not on the ramp.
axC_cbar.axis("off")
_cb = fig.colorbar(
    plt.cm.ScalarMappable(norm=mcolors.Normalize(0, C_SMAX), cmap=C_CMAP),
    ax=axC_cbar, fraction=0.30, pad=0.02, aspect=11)
_cb.set_label("S, one-step supply", fontsize=22, labelpad=14)
_cb.ax.tick_params(labelsize=20)
_cbpos = _cb.ax.get_position()
fig.patches.append(mpatches.Rectangle(
    (_cbpos.x0, _cbpos.y0 - 0.032), _cbpos.width, 0.013,
    facecolor=C_BLOCKED, edgecolor="white", linewidth=1.0,
    transform=fig.transFigure, figure=fig))
fig.text(_cbpos.x1 + 0.004, _cbpos.y0 - 0.0255, "class not\nsupplied",
         ha="left", va="center", fontsize=20, linespacing=1.15,
         color="#1a1a1a")

fig.text(axA.get_position().x0, axC[0].get_position().y1 + 0.028,
         "C   Per-residue reachability of the top-ranked substitution",
         ha="left", va="bottom", fontsize=32, fontweight="bold",
         color="#1a1a1a")

out_png = PROJECT_ROOT / "paper" / "figure_3.png"
out_pdf = PROJECT_ROOT / "paper" / "figure_3.pdf"
fig.savefig(out_png, dpi=300, bbox_inches="tight")
fig.savefig(out_pdf,              bbox_inches="tight")
print(f"Wrote {out_png}")
