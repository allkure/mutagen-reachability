# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
# ---

# %% [markdown]
# # Prospective predictions across three ALK TKI generations
#
# Extends `notebook/04_oh_decomposition.py` from alectinib (2nd gen) to
# lorlatinib (3rd gen) and zotizalkib (4th gen), using Oh 2026 Table S2's
# per-variant drug-selection scores for all three drugs.
#
# For each (mutagen, drug) pair we compute:
#   - log P̂_simulated(v)   = log S_mutagen(v) + (ResistanceScore - FitnessScore)
#   - log P̂_chemfree(v)    =                     (ResistanceScore - FitnessScore)
#
# where ResistanceScore and FitnessScore come from Oh Table S2 at AA-level
# (mean across SNVs encoding the same AA change).
#
# Chemistry-under-counted variants: those with high P̂_chemfree rank (fit per
# Oh) but low P̂_simulated rank (chemistry-rare under the mutagen). These are
# falsifiable predictions: "preclinical screens using this mutagen will miss
# these variants even though they are fit under this drug."
#
# Caveat (same as notebook/04): Oh's σ is measured in H3122 (ALK+ NSCLC human
# cells), not Ba/F3. For Doi-style Ba/F3 Mutazyme screens we assume
# σ_H3122 ≈ σ_BaF3 up to a scale — this is validated at the top of the
# alectinib ranking (V1180L #1, L1196M #2 match Doi 2022 Fig 3a observations).

# %%
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

NOTEBOOK_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = NOTEBOOK_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.harmonized_cohort import classify_amino_acid_change, load_alk_codons
from src.spectrum import load_spectra

ALK_CODONS = load_alk_codons()
SPECTRA = load_spectra(PROJECT_ROOT / "data" / "reference" / "spectra.yaml")

DRUGS = {
    "alectinib":  ("alec_score",  "2nd gen"),
    "lorlatinib": ("lor_score",   "3rd gen"),
    "zotizalkib": ("zot_score",   "4th gen"),
}
MUTAGENS = {
    "MutazymeII": SPECTRA["MutazymeII_DoiSupp2022"],
    "ENU":        SPECTRA["ENU_Arnold2012_nominal"],
}

# %%
# 1. Oh atlas AA-level aggregation (as in notebook/04).
oh = pd.read_csv(PROJECT_ROOT / "data" / "raw" / "oh_2026" / "table_s2_resistance_scores.csv")
oh_m = oh[oh["type"] == "Missense"].copy()
oh_m["residue"] = pd.to_numeric(oh_m["residue"], errors="coerce")
oh_m = oh_m[oh_m["residue"].between(1116, 1392)].copy()
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
print(f"Oh AA-level in ALK KD (1116-1392): {len(oh_aa)} unique AA changes")

# %%
# 2. For each mutagen, annotate class-forced supply weight.
rows = []
for mutagen_name, spectrum in MUTAGENS.items():
    for _, r in oh_aa.iterrows():
        info = classify_amino_acid_change(
            r["wt_aa"], r["mut_aa"], wt_codon=ALK_CODONS.get(int(r["residue"]))
        )
        if not info["class_forced"]:
            continue
        rows.append({
            "mutagen": mutagen_name,
            "aa_change": r["aa_change"],
            "residue": int(r["residue"]),
            "forced_class": info["forced_class"],
            "S_mutagen": spectrum[info["forced_class"]],
            "alec_score": r["alec_score"],
            "lor_score": r["lor_score"],
            "zot_score": r["zot_score"],
            "fitness_score": r["fitness_score"],
        })
annotated = pd.DataFrame(rows)
print(f"Class-forced per mutagen: {annotated.groupby('mutagen').size().to_dict()}")

# %%
# 3. For each (mutagen × drug) compute drug-selected scores and rankings.
summary_rows = []
for drug, (score_col, gen) in DRUGS.items():
    for mutagen in MUTAGENS:
        sub = annotated[annotated["mutagen"] == mutagen].copy()
        # σ_drug(v) = Resistance(v). Per Oh et al. 2026 Methods, the atlas's
        # Resistance score is already a same-timepoint drug-vs-vehicle contrast
        # (log(drug_D10 / untreated_D10)), so intrinsic growth is already
        # divided out by the common untreated-D10 denominator. Fitness
        # (log(untreated_D10 / D0)) is an orthogonal time-course growth axis;
        # subtracting it from Resistance would double-correct and systematically
        # bias the ranking toward drug-addiction-like variants (R>>0, F<<0) and
        # against driver-like variants (R>>0, F>>0). Confirmed empirically:
        # r(resistance, fitness) ≈ 0 across all three drugs. The two scores are
        # explicitly used by Oh as the x- and y-axes of a 2D variant classifier.
        sub["drug_selec"] = sub[score_col]
        sub["log_P_simulated"] = np.log(sub["S_mutagen"]) + sub["drug_selec"]
        sub["log_P_chemfree"]  = sub["drug_selec"]
        sub["rank_simulated"] = sub["log_P_simulated"].rank(ascending=False, method="min")
        sub["rank_chemfree"]  = sub["log_P_chemfree"].rank(ascending=False, method="min")
        sub["rank_shift"] = sub["rank_simulated"] - sub["rank_chemfree"]
        sub["drug"] = drug
        sub["drug_gen"] = gen
        summary_rows.append(sub)
all_results = pd.concat(summary_rows, ignore_index=True)

# %%
# 4. Rank-correlation (chemistry impact) per (mutagen × drug).
print("\n=== Spearman ρ (simulated-preclinical rank vs chemistry-free rank) ===")
corr_table = []
for drug, (score_col, gen) in DRUGS.items():
    for mutagen in MUTAGENS:
        sub = all_results[(all_results["mutagen"] == mutagen) & (all_results["drug"] == drug)]
        rho = sub[["log_P_simulated", "log_P_chemfree"]].corr(method="spearman").iloc[0, 1]
        pear = sub[["log_P_simulated", "log_P_chemfree"]].corr(method="pearson").iloc[0, 1]
        corr_table.append({"mutagen": mutagen, "drug": drug, "gen": gen,
                           "n_variants": len(sub), "spearman_rho": rho, "pearson_r": pear})
        print(f"  {mutagen:10s} × {drug:11s} ({gen}): N={len(sub):4d}  ρ={rho:+.3f}  r={pear:+.3f}")
pd.DataFrame(corr_table).to_csv(
    PROJECT_ROOT / "data" / "derived" / "chemistry_impact_by_mutagen_drug.csv", index=False
)

# %%
# 5. Top-N chemistry-under-counted variants per (mutagen × drug).
# "Chemistry-under-counted" = high rank_chemfree (fit), low rank_simulated (chemistry-rare).
# Rank_shift > 0 means mutagen chemistry PUSHES the variant DOWN the ranking.
# Report top 15 under-counted variants per combination that are in top-50 chemistry-free.
print("\n\n=== CHEMISTRY-UNDER-COUNTED PREDICTIONS per (mutagen × drug) ===")
all_predictions = []
for drug, (score_col, gen) in DRUGS.items():
    for mutagen in MUTAGENS:
        sub = all_results[(all_results["mutagen"] == mutagen) & (all_results["drug"] == drug)].copy()
        top_fit = sub[sub["rank_chemfree"] <= 50].copy()
        top_fit = top_fit.sort_values("rank_shift", ascending=False).head(15)
        if len(top_fit) == 0:
            continue
        print(f"\n  {mutagen} × {drug} ({gen}) — top-50-fit variants, ranked by chemistry under-counting:")
        print(top_fit[["aa_change", "forced_class", "S_mutagen", "drug_selec",
                       "rank_chemfree", "rank_simulated", "rank_shift"]].to_string(index=False))
        top_fit["mutagen"] = mutagen
        top_fit["drug"] = drug
        top_fit["drug_gen"] = gen
        all_predictions.append(top_fit)
pd.concat(all_predictions, ignore_index=True).to_csv(
    PROJECT_ROOT / "data" / "derived" / "chemistry_undercount_predictions.csv", index=False
)

# %%
# 6. Save a per-variant table of simulated rankings for the paper's main figure.
out = PROJECT_ROOT / "data" / "derived" / "oh_decomposition_all_drugs.csv"
all_results.to_csv(out, index=False)
print(f"\nWrote {out}")
print(f"Wrote data/derived/chemistry_impact_by_mutagen_drug.csv")
print(f"Wrote data/derived/chemistry_undercount_predictions.csv")
