# Title & abstract — draft for review (2026-10-14)

**Target journal**: Genome Biology (Research article).
**Rationale**: methods-driven, data-reuse, actionable for screen designers.
**Framing per MR 2026-10-14**: keep both case studies to show the framework
generalises; lead with KRAS because the inhibitor pipeline is still being
developed and the predictions are forward-looking.

---

## Title (locked, 2026-10-14 late pm)

> **Quantifying the preclinical blind spots of mutagen-based TKI-resistance screens from chemistry-free saturation mutagenesis atlases**

Rationale: "quantifying" names the primary scientific act of the paper
(magnitudes of rank shifts, recovery curves, class-supply percentages,
Jeffreys intervals on spectrum estimates); predictions are derivative and
appear in the Discussion's Prospective paragraph. Drops specific kinase
names in favour of the framework-level scope ("mutagen-based TKI-resistance
screens") since the decomposition $P(v) = S_{\text{mutagen}}(v)\cdot\sigma(v)$
is drug-agnostic and the ALK + KRAS G12C case studies demonstrate it; names
the data (chemistry-free atlases) as the asset the paper actually brings;
avoids "miss" (the well-known half of the story). No colon. Keywords and
abstract still carry "ALK" and "KRAS G12C" so the paper remains searchable
by anyone looking for those specific kinases. Earlier draft used
"Predicting"; switched to "Quantifying" on 2026-10-14 pm after the paper's
caveat structure (Jeffreys intervals, "sign and order of magnitude are
robust but point estimates are not", neladalkib-atlas gap) made the
quantitative framing more defensible than the predictive one.

---

## Abstract (280 words)

**Background.** Ba/F3-style preclinical mutagenesis screens remain the first
line of evidence for predicting kinase-inhibitor resistance variants, with
chemistry-biased mutagens — ENU and Mutazyme II — supplying variants whose
pyrimidine substitution classes vary by two orders of magnitude in
probability. Whether and by how much this chemistry bias distorts
preclinical rankings relative to clinically observed resistance has not been
quantified against a chemistry-free ground truth.

**Methods.** We decompose preclinical screen output as a product of
mutagen supply and drug selection, P(v) = S_mutagen(v) · σ(v), with σ(v)
read from two recently published chemistry-free saturation mutagenesis
atlases: the Oh et al. (2026) prime-editing atlas for ALK (3,208 SNVs under
alectinib, lorlatinib, zotizalkib) and the Feng et al. (2022) Twist
Biosciences synthetic library for KRAS G12C (3,567 variants under sotorasib
and adagrasib). Each variant is classified by the pyrimidine substitution
class forced by its wild-type codon.

**Results.** For KRAS G12C, two of twelve clinically observed Awad 2021
resistance variants — including R68L, the top-ranked variant in Feng's
atlas under both drugs — require two or more nucleotide changes from the
WT codon and are structurally unreachable by any single-nucleotide-variant
mutagen. Among the seven class-forced clinical variants, chemistry
restructures rankings by up to 1,048 positions (H95D under adagrasib + ENU,
where C>G accounts for 0.2 % of ENU's output). The clinically "novel"
variants R68L and Y96D — described as previously unidentified by Tanaka
et al. (2021) — are precisely those our framework predicts should be
absent or deeply buried in standard ENU preclinical screens. The ALK case
reproduces the Doi et al. (2022) Mutazyme × alectinib headline observations
(V1180L rank 1, L1196M rank 2) with Spearman ρ = 0.88 between simulated
and chemistry-free rankings.

**Conclusions.** Chemistry bias in preclinical TKI-resistance screens is
both quantitative (ranking shift) and structural (SNV-unreachability).
Screen designs for the next generation of KRAS G12C inhibitors should
either use chemistry-free libraries (Twist, prime-editing) or acknowledge
the specific clinical variants their mutagen cannot produce.

---

## Figure 1 — concept + money shot

Four-panel figure, A4 landscape at 2 × 2 grid.

- **Panel A (concept)**: Schematic. Horizontal arrow from "KRAS G12C library
  (uniform)" through two routes — "Mutazyme II / ENU (chemistry-biased)"
  vs "Twist synthetic (chemistry-free)" — into "Drug-treated pool". A box
  below each route lists class-forced variants that would be observed. Key
  text: P_screen(v) = S_mutagen(v) · σ(v).

- **Panel B (KRAS Awad 2021 panel categorised)**: One row per clinical
  variant; three x-axis bins — "SNV-unreachable", "multi-class SNV",
  "class-forced SNV". R68L anchored at the top of the SNV-unreachable bin
  with Feng Z = +4.34 (adagrasib). Point colour = Feng Z-score magnitude.
  Labels on each point.

- **Panel C (KRAS rank-shift scatter)**: y = rank under ENU simulation,
  x = rank under chemistry-free Feng. One point per class-forced variant
  (n = 981). Diagonal = no-shift line. Highlight H95D and Y96D as off-
  diagonal outliers. Zoom inset at top-50 fit shows how clinical variants
  cluster vs non-clinical hits.

- **Panel D (ALK parallel)**: Same scatter structure for Oh × MutazymeII
  × alectinib, demonstrating the framework generalises. Highlight V1180L
  and L1196M (Doi 2022 top-2) and I1171N as off-diagonal clinical outlier.

Rendered as `paper/figure_1.png` with the plotting script
`paper/figure_1.py`.
