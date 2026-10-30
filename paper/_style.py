"""Shared matplotlib style for all paper figures.

Pastel, colour-blind-friendly palette derived from Okabe & Ito (2008)
with ~25% desaturation and no coloured text anywhere — all labels in
black or neutral grey; colour lives in bar/point fills and legend swatches.
"""

import matplotlib.pyplot as plt


# Okabe-Ito-derived pastels (safe for deuteranopia + protanopia + tritanopia)
PALETTE = {
    # Four mutagen chemistries — hues at ~90° so distinguishable in
    # deuteranopia/protanopia/tritanopia AND grayscale. Light enough to
    # read as a pastel palette but with enough chroma that they do not
    # blur into neighbours.
    "chemfree":      "#7FA8D1",  # soft cyan-blue
    "ENU":           "#E68B88",  # salmon-red
    "MutII":         "#E8C45E",  # warm gold (not amber — further from ENU)
    "BE":            "#B07CC6",  # purple (not lavender — further from blue)

    # Three Panel-B categories (lighter so markers' black edges stand out)
    "accessible":    "#B8D4EC",  # light cyan-blue
    "undercounted":  "#F3D894",  # light gold
    "unreachable":   "#EDA19F",  # light salmon

    "clinical":      "#365F8F",  # darker blue — highlighted clinical point
    "neutral":       "#5C6066",  # neutral grey — dividers, swatches
    "neutral_soft":  "#9EA3AA",  # softer grey — dashed separators
    "cream":         "#F5EFE0",  # warm cream — neutral box fill, not grey
    "grid":          "#C8CCD1",  # pale grey — scatter background points
    "gene_kras":     "#C78F15",  # muted gold — KRAS gene swatch
    "gene_alk":      "#3B6FA8",  # muted blue — ALK gene swatch
}

# Six-class COSMIC-adjacent pastels (desaturated for the spectrum bar chart)
CLASS_COLORS = ["#9ECCE4", "#9FA3A8", "#E89C9B", "#DDDAD6", "#AED188", "#D7B7BD"]
# order: C>A  C>G  C>T  T>A  T>C  T>G


def apply_rc(scale=1.0):
    """Apply matplotlib rcParams for all paper figures."""
    s = scale
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size":          int(26 * s),
        "axes.titlesize":     int(28 * s),
        "axes.labelsize":     int(26 * s),
        "xtick.labelsize":    int(24 * s),
        "ytick.labelsize":    int(24 * s),
        "legend.fontsize":    int(24 * s),
        "axes.spines.top":    False,
        "axes.spines.right":  False,
        "axes.linewidth":     1.4,
        "xtick.major.width":  1.4,
        "ytick.major.width":  1.4,
        "axes.labelcolor":    "#1a1a1a",
        "text.color":         "#1a1a1a",
    })
