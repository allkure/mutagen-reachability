"""Figure 4 — Retrospective validation of P(v) = S·σ across three preclinical ALK screens.

Panel A: Mutazyme II × alectinib Ba/F3 screen (single-mutation, primary anchor).
Panel B: ePCR × neladalkib on an I1171N fixed-V1 background (compound-mutation).
Panel C: ENU × multi-drug industry screen (aggregate hits across arms).

Design notes (v2, 2026-10-15):
  - All text sizes increased for printed-column legibility.
  - Grey subtitles (that overlaid other elements) removed; the context is
    now in the figure caption, not in the panel.
  - Study names removed from panel titles; descriptive labels only.
"""
from pathlib import Path
import sys
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parent
sys.path.insert(0, str(PROJECT_ROOT))
from paper._style import *  # noqa

DATA = PROJECT_ROOT / "data"

df = pd.read_csv(DATA / "derived" / "retrospective_validation.csv")

# Pyrimidine class → color (consistent with other figures)
CLASS_COLORS = {
    "C>A": "#4e79a7", "C>G": "#e15759", "C>T": "#59a14f",
    "T>A": "#f28e2b", "T>C": "#76b7b2", "T>G": "#b07aa1",
}

# ---- Figure layout (large; prints at ~A3 landscape to keep text legible).
# Supertitle gets its own ~7% of the vertical space at the top; gridspec
# occupies the rest so panel A's title does not collide with the supertitle.
fig = plt.figure(figsize=(24, 22))
gs = fig.add_gridspec(nrows=3, ncols=1, left=0.055, right=0.985, top=0.905, bottom=0.025,
                      hspace=0.30, height_ratios=[6.5, 5.5, 13.5])
ax_doi = fig.add_subplot(gs[0, 0])
ax_tak = fig.add_subplot(gs[1, 0])
ax_nuv = fig.add_subplot(gs[2, 0])

# Font sizes — ONE place to tune. All bumped ~1.4x from v2.
TITLE_FS      = 32
HEADER_FS     = 24
VARIANT_FS    = 24
CODON_FS      = 22
REACH_FS      = 21
CLASS_FS      = 22
PCT_FS        = 23
RANK_FS       = 23
SUPTITLE_FS   = 30

# ---- Helper: draw one screen as a horizontal table
def draw_panel(ax, screen_df, title_text):
    screen_df = screen_df.reset_index(drop=True).copy()
    n = len(screen_df)
    ax.set_xlim(-0.3, 6.9)
    ax.set_ylim(-1.0, n + 1.1)
    ax.invert_yaxis()
    ax.axis('off')

    # Panel title (bold, large, dark). Left-anchored, well above the header row.
    ax.text(-0.25, -0.65, title_text, fontsize=TITLE_FS, fontweight='bold',
            color='#1a1a1a', ha='left', va='bottom')

    # Column headers
    cols_x = [0.0, 1.3, 2.5, 3.55, 4.55, 5.70, 6.85]
    hdrs   = ["Variant", "WT codon", "Reachability", "Class",
              "Rate (%)", "Chem-free rank", "Simulated rank"]
    for i, (x, h) in enumerate(zip(cols_x, hdrs)):
        align = 'left' if i < 2 else 'center'
        ax.text(x, 0.3, h, fontsize=HEADER_FS, fontweight='bold', color='#1a1a1a',
                ha=align, va='bottom')
    ax.hlines(y=0.55, xmin=-0.25, xmax=6.95, colors='#555', linewidth=1.6)

    # Rows
    for i, row in screen_df.iterrows():
        y = i + 1.15
        # Variant name
        ax.text(cols_x[0], y, row['variant'], fontsize=VARIANT_FS, fontweight='bold',
                family='monospace', color='#1a1a1a', ha='left', va='center')
        # WT codon
        ax.text(cols_x[1], y, str(row.get('wt_codon', '?')), fontsize=CODON_FS,
                family='monospace', color='#333', ha='left', va='center')
        # Reachability
        reach_text = str(row.get('reach', '?'))
        if reach_text in ('nan', '?', ''):
            reach_text, reach_color = 'unanchored', '#888'
        else:
            reach_color = '#1a1a1a' if reach_text == 'class-forced' else '#a05200'
        ax.text(cols_x[2], y, reach_text, fontsize=REACH_FS, color=reach_color,
                ha='center', va='center')
        # Forced class — colored chip (sized for bigger fonts)
        cls = row.get('class', None)
        if cls and str(cls) != 'nan' and cls in CLASS_COLORS:
            chip = Rectangle((cols_x[3] - 0.52, y - 0.40), 1.04, 0.80,
                             facecolor=CLASS_COLORS[cls], edgecolor='none')
            ax.add_patch(chip)
            ax.text(cols_x[3], y, cls, fontsize=CLASS_FS, color='white', fontweight='bold',
                    ha='center', va='center', family='monospace')
        else:
            ax.text(cols_x[3], y, '—', fontsize=CLASS_FS, color='#aaa',
                    ha='center', va='center')
        # Mutagen rate
        rate = row.get('mutagen_rate', None)
        if pd.notna(rate):
            pct = rate * 100
            txt = f"{pct:.1f}" if pct >= 1.0 else f"{pct:.2f}"
            color = '#c0392b' if pct < 1.0 else ('#e67e22' if pct < 10 else '#1a1a1a')
            weight = 'bold' if pct < 1.0 else 'normal'
            ax.text(cols_x[4], y, txt, fontsize=PCT_FS, color=color, fontweight=weight,
                    ha='center', va='center', family='monospace')
        else:
            ax.text(cols_x[4], y, '—', fontsize=PCT_FS, color='#aaa',
                    ha='center', va='center')
        # Chem-free rank
        rc = row.get('rank_chemfree', None)
        if pd.notna(rc):
            ax.text(cols_x[5], y, f"{int(rc)}", fontsize=RANK_FS,
                    family='monospace', color='#1a1a1a', ha='center', va='center')
        else:
            ax.text(cols_x[5], y, '—', fontsize=RANK_FS, color='#aaa',
                    ha='center', va='center')
        # Simulated rank
        rs = row.get('rank_simulated', None)
        if pd.notna(rs):
            ax.text(cols_x[6], y, f"{int(rs)}", fontsize=RANK_FS,
                    family='monospace', color='#1a1a1a', ha='center', va='center')
        else:
            ax.text(cols_x[6], y, '—', fontsize=RANK_FS, color='#aaa',
                    ha='center', va='center')

# Panel A — Mutazyme II × alectinib (single-mutation)
doi_df = df[df['screen'] == 'Doi 2022'].copy()
draw_panel(ax_doi, doi_df,
           "A   Mutazyme II  ×  alectinib   (single-mutation anchor)")

# Panel B — ePCR × neladalkib on I1171N fixed-V1
tak_df = df[df['screen'] == 'Takei 2025'].copy()
draw_panel(ax_tak, tak_df,
           "B   ePCR  ×  neladalkib   (compound-mutation, I1171N fixed-V1)")

# Panel C — ENU × multi-drug industry screen
nuv_df = df[df['screen'] == 'Nuvalent 2025'].copy()
draw_panel(ax_nuv, nuv_df,
           "C   ENU  ×  neladalkib + alectinib + lorlatinib   (industry screen, aggregate)")

# Supertitle — placed above the gridspec, with its own clear band at the top.
fig.suptitle(r"Retrospective validation of  $P(v) = S_{\mathrm{mutagen}}(v) \cdot \sigma(v)$  across three preclinical ALK TKI resistance screens",
             fontsize=SUPTITLE_FS, fontweight='bold', y=0.965)

fig.savefig(HERE / "figure_4.png", dpi=300, bbox_inches='tight')
fig.savefig(HERE / "figure_4.pdf", bbox_inches='tight')
print(f"Wrote {HERE / 'figure_4.png'}")
print("Wrote figure_4.pdf and figure_4.png")
