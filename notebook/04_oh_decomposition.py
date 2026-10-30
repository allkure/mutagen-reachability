# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
# ---

# %% [markdown]
# # Two-term decomposition using Oh 2026 as σ reference
#
# Model:
#     log P_screen(v) = α · log S_mutagen(v) + β · OhResistanceScore(v, drug) + const
#
# where
#     - P_screen(v) is observed variant frequency in a mutagenesis screen
#     - S_mutagen(v) is class-forced supply weight for mutagen's spectrum
#     - OhResistanceScore(v, drug) is Oh 2026 Table S2's AA-level resistance score,
#       which is a z-score-standardized log2-fold-change of (drug-treated / untreated)
#       variant frequency, measured by prime editing in H3122 cells.
#       Monotonic surrogate for log σ_drug(v).
#
# Oh methodology (from Oh 2026 Methods, p.5): LFCs of SNV frequencies between
# drug-treated and untreated conditions, standardized to the distribution of
# synonymous-SNV LFCs, then averaged across replicates. For protein-level scores
# (SAAVs), the mean of the SNV-level scores is used.
#
# Caveat: Oh used H3122 (ALK+ NSCLC) cells, not Ba/F3. The σ we're borrowing is
# "H3122 selection coefficient under drug," not Ba/F3. For Doi 2022 (Ba/F3 Mutazyme
# screen of alectinib), we are implicitly assuming σ_H3122 ≈ σ_BaF3 up to a scale.
# This is an approximation flagged in the preregistration §9.

# %%
from __future__ import annotations

import sys
import re
from pathlib import Path

import numpy as np
import pandas as pd

NOTEBOOK_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = NOTEBOOK_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.harmonized_cohort import classify_amino_acid_change, load_alk_codons
from src.spectrum import load_spectra
from src.codon_tools import SUB_CLASSES

ALK_CODONS = load_alk_codons()
SPECTRA = load_spectra(PROJECT_ROOT / "data" / "reference" / "spectra.yaml")

# %%
# 1. Load Oh 2026 Table S2.
oh_raw = pd.read_csv(PROJECT_ROOT / "data" / "raw" / "oh_2026" / "table_s2_resistance_scores.csv")
oh_m = oh_raw[oh_raw["type"] == "Missense"].copy()
oh_m["residue"] = pd.to_numeric(oh_m["residue"], errors="coerce")
oh_m = oh_m[oh_m["residue"].between(1116, 1392)].copy()
print(f"Oh atlas missense in ALK KD: {len(oh_m)} SNV-level rows")

# 2. Collapse to AA-level using MEAN across SNVs encoding the same AA change
# (Oh 2026 Methods, p.5: "we used the mean of SNV resistance scores to calculate
# the resistance score of each protein variant").
oh_aa = oh_m.groupby("aa_change").agg(
    alec_score=("alec_score", "mean"),
    lor_score=("lor_score", "mean"),
    zot_score=("zot_score", "mean"),
    fitness_score=("fitness_score", "mean"),
    wt_aa=("wt_aa", "first"),
    mut_aa=("mut_aa", "first"),
    residue=("residue", "first"),
    n_snvs=("aa_change", "size"),
).reset_index()
print(f"After AA-level mean aggregation: {len(oh_aa)} unique AA changes")
print(f"  Distribution of n_snvs per AA: "
      f"min={oh_aa['n_snvs'].min()}, median={oh_aa['n_snvs'].median()}, "
      f"max={oh_aa['n_snvs'].max()}")

# %%
# 3. Compute class-forced S_Mutazyme(v) for every AA change using the ALK WT codon
# lookup (data/reference/codons.yaml) where available, else enumerating over
# synonymous codons. This matches the methodology in src.harmonized_cohort.
mut_spectrum = SPECTRA["MutazymeII_DoiSupp2022"]
enu_spectrum = SPECTRA["ENU_Arnold2012_nominal"]


def annotate(row, spectrum, label):
    """Return dict with class-forced status + supply weight for this AA change."""
    wt = row["wt_aa"]
    res = int(row["residue"])
    mut = row["mut_aa"]
    info = classify_amino_acid_change(wt, mut, wt_codon=ALK_CODONS.get(res))
    return pd.Series({
        f"class_forced_{label}": info["class_forced"],
        f"forced_class_{label}": info["forced_class"],
        f"reachable_1nt_{label}": info["reachable_1nt"],
        f"S_{label}": (spectrum[info["forced_class"]] if info["class_forced"] else None),
    })


oh_aa = pd.concat([oh_aa, oh_aa.apply(lambda r: annotate(r, mut_spectrum, "mut"), axis=1)], axis=1)
oh_aa = pd.concat([oh_aa, oh_aa.apply(lambda r: annotate(r, enu_spectrum, "enu"), axis=1)], axis=1)
oh_cf_mut = oh_aa[oh_aa["class_forced_mut"] == True].copy()
print(f"Class-forced AA changes (Mutazyme): {len(oh_cf_mut)}")

# %%
# 4. Load Doi 2022 Fig 3a observed alectinib-resistance variants.
doi = pd.read_csv(PROJECT_ROOT / "data" / "raw" / "doi_2022" / "fig3a_alectinib_counts.csv", comment="#")
doi_r = doi[doi["is_resistance_mutation"] == True].copy()
doi_r["p_doi"] = doi_r["peak_freq_approx"] / doi_r["peak_freq_approx"].sum()
print(f"\nDoi 2022 Fig 3a resistance variants: {len(doi_r)} (visible labels only)")

# %%
# 5. Merge Doi observations onto Oh atlas (class-forced, Mutazyme).
merged = oh_cf_mut.merge(
    doi_r[["variant", "p_doi"]].rename(columns={"variant": "aa_change"}),
    on="aa_change", how="left",
)
merged["observed_in_doi"] = merged["p_doi"].notna()

# %%
# 6. TEST 1 — Does Oh σ + chemistry predict Doi's observations?
#    Multiple regression: log(p_doi) ~ log(S_Mut) + alec_score.
#    Only on Doi-observed class-forced variants (k = small).
reg_data = merged[merged["observed_in_doi"]].dropna(subset=["S_mut", "alec_score"]).copy()
if len(reg_data) >= 3:
    y = np.log(reg_data["p_doi"].values + 1e-6)
    X = np.column_stack([
        np.ones(len(reg_data)),
        np.log(reg_data["S_mut"].values),
        reg_data["alec_score"].values,
    ])
    coef, residuals, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    y_pred = X @ coef
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    print(f"\n=== TEST 1 — Two-term fit on Doi 2022 class-forced observations ===")
    print(f"  k = {len(reg_data)} variants, 3 parameters, {len(reg_data)-3} residual df")
    print(f"  log(p_doi) = {coef[0]:+.3f} + {coef[1]:+.3f}·log(S_Mut) + {coef[2]:+.3f}·alec_score")
    print(f"  R² = {r2:+.3f}")
    # Univariate predictors
    for name, x in [("log(S_Mut)", np.log(reg_data["S_mut"].values)),
                     ("alec_score", reg_data["alec_score"].values)]:
        r = np.corrcoef(x, y)[0, 1]
        print(f"  Pearson(log(p_doi), {name}) = {r:+.3f}")

    print("\n  Per-variant residuals:")
    reg_data["log_p_obs"] = y
    reg_data["log_p_pred"] = y_pred
    reg_data["residual"] = y - y_pred
    print(reg_data[["aa_change", "forced_class_mut", "p_doi", "S_mut", "alec_score",
                     "log_p_obs", "log_p_pred", "residual"]].to_string(index=False))

# %%
# 7. TEST 2 — Simulated Mutazyme-alectinib ranking vs chemistry-free ranking.
# For every class-forced Mutazyme ALK KD missense variant, compute:
#   P_hat_mut(v)  ∝ S_Mut(v) · exp(alec_score(v))       # with chemistry
#   P_hat_free(v) ∝                exp(alec_score(v))   # chemistry-free (prime editing)
#
# Compare the two rankings. The rank-gap IS the chemistry contribution.
sim = oh_cf_mut.dropna(subset=["S_mut", "alec_score"]).copy()
sim["log_P_mut_alec"] = np.log(sim["S_mut"]) + sim["alec_score"]
sim["log_P_free_alec"] = sim["alec_score"]
sim["rank_mut"] = sim["log_P_mut_alec"].rank(ascending=False, method="min")
sim["rank_free"] = sim["log_P_free_alec"].rank(ascending=False, method="min")
spearman_mut_vs_free = sim[["log_P_mut_alec", "log_P_free_alec"]].corr(method="spearman").iloc[0, 1]
pearson_mut_vs_free = sim[["log_P_mut_alec", "log_P_free_alec"]].corr(method="pearson").iloc[0, 1]
print(f"\n=== TEST 2 — Simulated preclinical vs chemistry-free ranking (alectinib) ===")
print(f"  N variants: {len(sim)}  (all class-forced ALK KD missense in Oh × Mutazyme)")
print(f"  Spearman ρ (simulated-Mutazyme-rank vs chemistry-free-rank): {spearman_mut_vs_free:+.3f}")
print(f"  Pearson r on log scores:                                      {pearson_mut_vs_free:+.3f}")

# %%
# 8. TEST 3 — Doi's observed top-7 variants: where do they rank in P_hat_mut_alec?
# If the model is right, Doi's 7 observed resistance variants should rank in the
# top-N of P_hat_mut_alec, where N is small compared to total.
sim_top20 = sim.sort_values("rank_mut").head(20)
doi_in_sim = sim[sim["aa_change"].isin(doi_r["variant"])].sort_values("rank_mut")

print(f"\n=== TEST 3 — Simulated top-20 ranking (Mutazyme × alectinib) ===")
print(sim_top20[["aa_change", "forced_class_mut", "S_mut", "alec_score",
                 "log_P_mut_alec", "rank_mut", "rank_free"]].to_string(index=False))
print(f"\nDoi 2022 observed resistance variants and their simulated ranks:")
print(doi_in_sim[["aa_change", "forced_class_mut", "S_mut", "alec_score",
                   "log_P_mut_alec", "rank_mut", "rank_free"]].to_string(index=False))

# %%
# 9. Chemistry-under-counted and chemistry-over-counted variants.
# A variant is "chemistry-under-counted" if its chemistry-free rank is high
# (fit) but its simulated-Mutazyme rank is poor (low S_Mut). These are
# predictions: "a Mutazyme screen will miss them even though they are fit."
sim["under_count"] = sim["rank_mut"] - sim["rank_free"]
sim_top_free = sim.sort_values("rank_free").head(50).copy()
sim_top_free["chemistry_delta"] = sim_top_free["log_P_mut_alec"] - sim_top_free["log_P_free_alec"]
print(f"\n=== Chemistry-under-counted variants (fit per Oh, chemistry-rare in Mutazyme) ===")
print(f"  Looking at the top-50 chemistry-free (Oh σ_alec) variants:")
most_under = sim_top_free.sort_values("under_count", ascending=False).head(10)
print(most_under[["aa_change", "forced_class_mut", "S_mut", "alec_score",
                   "rank_free", "rank_mut", "under_count"]].to_string(index=False))

# %%
# 10. Save derived table for downstream analyses.
out = PROJECT_ROOT / "data" / "derived" / "oh_decomposition_alectinib.csv"
sim.to_csv(out, index=False)
print(f"\nWrote {out}")
