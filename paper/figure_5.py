#!/usr/bin/env python3
"""Figure 5 — the MET audit.

A  Same mutagen, two drugs. Observed six-class composition of the clones each
   screen recovered, against the germline ENU reference. If chemistry set the
   composition the two screens would agree with each other and with the
   reference; they do neither.
B  The Y1230 single-codon experiment. Every single-nucleotide route out of the
   wild-type codon, its substitution class, and the clones recovered against
   what the ENU class supply alone predicts.
C  Reachability of the chemistry-free top-100 MET priorities.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent
sys.path.insert(0, str(PROJECT_ROOT))

from paper._style import CLASS_COLORS, PALETTE, apply_rc  # noqa: E402

AUDIT = PROJECT_ROOT / "data" / "derived" / "met_audit.json"
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]


def main() -> None:
    apply_rc()
    a = json.loads(AUDIT.read_text())
    enu = a["germline_enu"]
    screens = list(a["screens"])

    fig = plt.figure(figsize=(19.0, 15.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.08],
                          hspace=0.46, wspace=0.34,
                          left=0.085, right=0.975, top=0.875, bottom=0.095)

    # ---------------------------------------------------------------- A ----
    axA = fig.add_subplot(gs[0, :])
    x = np.arange(len(CLASSES))
    w = 0.33
    for i, s in enumerate(screens):
        frac = [a["screens"][s]["observed_class_fraction"][c] for c in CLASSES]
        axA.bar(x + (i - 0.5) * w, frac, w,
                color=[CLASS_COLORS[j] for j in range(6)],
                edgecolor="black", linewidth=1.5,
                hatch="" if i == 0 else "///")
    axA.plot(x, [enu[c] for c in CLASSES], "o--", color=PALETTE["neutral"],
             markersize=15, linewidth=2.6, markerfacecolor="white",
             markeredgewidth=2.4, zorder=5)

    axA.set_xticks(x)
    axA.set_xticklabels(CLASSES, fontsize=26, fontweight="bold")
    axA.set_ylabel("fraction of clones", fontsize=25)
    axA.set_ylim(0, 1.02)
    axA.grid(axis="y", color=PALETTE["grid"], linewidth=0.9, alpha=0.7)
    axA.set_axisbelow(True)
    axA.tick_params(labelsize=23)

    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor="#CCCCCC", edgecolor="black",
                      linewidth=1.5,
                      label=f"{screens[0]}  ({a['screens'][screens[0]]['n_clones']} clones)"),
        plt.Rectangle((0, 0), 1, 1, facecolor="#CCCCCC", edgecolor="black",
                      linewidth=1.5, hatch="///",
                      label=f"{screens[1]}  ({a['screens'][screens[1]]['n_clones']} clones)"),
        plt.Line2D([0], [0], color=PALETTE["neutral"], marker="o", linestyle="--",
                   markersize=13, markerfacecolor="white", markeredgewidth=2.2,
                   linewidth=2.4, label="germline ENU reference"),
    ]
    axA.legend(handles=handles, loc="upper left", frameon=False, fontsize=21,
               handlelength=2.0, borderaxespad=0.3, ncol=1)

    h = a["two_screen_homogeneity"]
    axA.text(0.995, 0.97,
             "same mutagen, same cells, same lab\n"
             "compositions differ between the two screens\n"
             rf"$\chi^2$ = {h['chi2']:.1f},  df = {h['dof']},  $P$ = {h['p']:.3f}",
             transform=axA.transAxes, ha="right", va="top", fontsize=21,
             linespacing=1.45,
             bbox=dict(boxstyle="round,pad=0.5", facecolor=PALETTE["cream"],
                       edgecolor=PALETTE["neutral_soft"], linewidth=1.2))
    axA.set_title("A    One mutagen, two drugs: the class composition follows the drug",
                  loc="left", fontsize=27, fontweight="bold", pad=16)

    # ---------------------------------------------------------------- B ----
    axB = fig.add_subplot(gs[1, 0])
    y = a["y1230"]
    items = sorted(y["per_variant"].items(), key=lambda kv: -kv[1]["clones_observed"])
    names = [k for k, _ in items]
    obs = [v["clones_observed"] for _, v in items]
    exp = [y["expected_clones_under_supply"][k] for k in names]
    cls = ["/".join(r["class"] for r in v["routes"]) for _, v in items]
    yy = np.arange(len(names))[::-1].astype(float)

    axB.barh(yy + 0.20, obs, 0.38, color=PALETTE["ENU"], edgecolor="black",
             linewidth=1.4, label="clones recovered", zorder=3)
    axB.barh(yy - 0.20, exp, 0.38, color=PALETTE["cream"], edgecolor=PALETTE["neutral"],
             linewidth=1.4, hatch="///", label="expected from ENU supply alone", zorder=3)
    for i, (o, c) in enumerate(zip(obs, cls)):
        axB.text(-0.8, yy[i], c, ha="right", va="center", fontsize=20,
                 color=PALETTE["neutral"])
        if o == 0:
            axB.text(0.8, yy[i] + 0.20, "not recovered", ha="left", va="center",
                     fontsize=19, style="italic", color=PALETTE["clinical"], zorder=4)

    axB.set_yticks(yy)
    axB.set_yticklabels(names, fontsize=23)
    axB.set_xlabel("clones", fontsize=25)
    axB.set_xlim(-9.5, 36)
    axB.set_ylim(-0.75, len(names) - 0.15)
    axB.tick_params(labelsize=23)
    axB.grid(axis="x", color=PALETTE["grid"], linewidth=0.9, alpha=0.7)
    axB.set_axisbelow(True)
    axB.legend(loc="lower center", ncol=2, frameon=False, fontsize=19,
               bbox_to_anchor=(0.5, -0.30), columnspacing=1.6, handlelength=1.4)
    axB.set_title(f"B    One-step routes out of Y1230 ({y['wt_codon']})",
                  loc="left", fontsize=27, fontweight="bold", pad=16)
    span = (y["supply_range_of_those_classes"][1]
            / y["supply_range_of_those_classes"][0])
    axB.text(0.0, -0.42,
             f"{y['n_recovered']} of {y['n_single_nt_missense_routes']} routes recovered, "
             f"spanning classes ENU supplies over a {span:.1f}-fold range",
             transform=axB.transAxes, ha="left", va="top", fontsize=20,
             color=PALETTE["neutral"])

    # ---------------------------------------------------------------- C ----
    axC = fig.add_subplot(gs[1, 1])
    order = ["class-forced", "multi-class SNV", "two-step"]
    colors = [PALETTE["accessible"], PALETTE["undercounted"], PALETTE["unreachable"]]
    bottom = np.zeros(len(screens))
    for cat, col in zip(order, colors):
        vals = np.array([a["screens"][s]["top100_reachability"].get(cat, 0)
                         for s in screens], dtype=float)
        axC.bar(np.arange(len(screens)), vals, 0.5, bottom=bottom, color=col,
                edgecolor="black", linewidth=1.4, label=cat, zorder=3)
        for i, (v, b) in enumerate(zip(vals, bottom)):
            if v >= 8:
                axC.text(i, b + v / 2, f"{int(v)}", ha="center", va="center",
                         fontsize=25, fontweight="bold", zorder=4)
        bottom += vals

    axC.set_xticks(np.arange(len(screens)))
    axC.set_xticklabels(["NVP-BVU972", "AMG 458"], fontsize=23)
    axC.set_ylabel("variants in the top 100", fontsize=24)
    axC.set_ylim(0, 108)
    axC.tick_params(labelsize=23)
    axC.grid(axis="y", color=PALETTE["grid"], linewidth=0.9, alpha=0.7)
    axC.set_axisbelow(True)
    axC.legend(loc="lower center", ncol=3, frameon=False, fontsize=19,
               bbox_to_anchor=(0.5, -0.20), columnspacing=1.4, handlelength=1.4)
    axC.set_title("C    What a one-step mutagen could reach", loc="left",
                  fontsize=27, fontweight="bold", pad=16)

    fig.suptitle("MET: auditing an ENU screen against a chemistry-free atlas",
                 fontsize=31, fontweight="bold", y=0.962)

    for ext in ("pdf", "png"):
        fig.savefig(HERE / f"figure_5.{ext}", dpi=300, bbox_inches="tight")
    print(f"wrote {HERE / 'figure_5.pdf'} and .png")


if __name__ == "__main__":
    main()
