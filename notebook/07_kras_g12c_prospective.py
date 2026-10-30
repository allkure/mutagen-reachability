# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
# ---

# %% [markdown]
# # KRAS G12C chemistry decomposition — Feng 2022 atlas × preclinical mutagens
#
# Replicates the ALK framework (notebook/05) in a KRAS G12C context, using
# Feng et al. 2022 PNAS (DOI 10.1073/pnas.2120512119) as the chemistry-free
# σ atlas for sotorasib and adagrasib. Feng used TWIST Biosciences synthetic
# DNA — a saturation mutagenesis library whose variants are oligo-synthesised
# (not PCR- or ENU-generated), so every single-AA substitution of KRAS-G12C
# enters the pool at ~uniform baseline. This is the KRAS analogue of the
# prime-editing library in Oh 2026.
#
# The DECOMPOSITION is slightly different from the ALK case:
#   - Oh 2026 reports a Resistance score AND a Fitness score separately;
#     drug_selection = ResistanceScore − FitnessScore.
#   - Feng 2022 reports a single Z-score per (variant, drug), normalised to
#     the untreated control (DMSO Day 0). The Z-score IS the drug-selection
#     signal, so we use it directly.
#
# For each class-forced variant v, mutagen m, drug d:
#     log_P_simulated(v)  = log S_m(v) + Z_d(v)
#     log_P_chemfree(v)   =             Z_d(v)
#
# Chemistry-under-counted variants: variants fit per the Feng atlas but
# chemistry-rare per the mutagen. These are variants that a KRAS G12C
# preclinical screen using Mutazyme II or ENU would MISS even though they
# are drug-selected.
#
# Validation anchor: Awad 2021 (NEJM 384:2382) reported clinical sotorasib/
# adagrasib resistance mutations in patients: Y96D, R68S/M, H95D/Q/R,
# V8L/M, C12F. These are the real-world "escape routes" and provide a
# ground-truth list to cross-reference our predictions.

# %%
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

NOTEBOOK_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = NOTEBOOK_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.harmonized_cohort import classify_amino_acid_change
from src.spectrum import load_spectra

SPECTRA = load_spectra(PROJECT_ROOT / "data" / "reference" / "spectra.yaml")


def load_kras_codons() -> dict[int, str]:
    """Return {residue -> WT codon} for KRAS-4A NM_033360.4."""
    import yaml
    with open(PROJECT_ROOT / "data" / "reference" / "codons_kras.yaml") as f:
        raw = yaml.safe_load(f)
    return {int(r): v["codon"] for r, v in raw["kras_4a_nm_033360_4"]["codons"].items()}


KRAS_CODONS = load_kras_codons()
# The Feng library has G12C baked in: position 12 codon changes GGT → TGT.
KRAS_CODONS[12] = "TGT"

MUTAGENS = {
    "MutazymeII": SPECTRA["MutazymeII_DoiSupp2022"],
    "ENU":        SPECTRA["ENU_Arnold2012_nominal"],
}
DRUGS = ["sotorasib", "adagrasib"]

# %%
# 1. Load Feng atlas (3568 scored variants per drug).
feng = pd.read_csv(PROJECT_ROOT / "data" / "raw" / "feng_2022" / "table_s1_zscores.csv")
# aa_change column mirrors the Oh table format for compatibility.
feng["aa_change"] = feng["wt_aa"] + feng["position"].astype(str) + feng["mut_aa"]
print(f"Feng atlas: {len(feng)} rows ({len(feng.dropna(subset=['sotorasib_z']))} sotorasib, "
      f"{len(feng.dropna(subset=['adagrasib_z']))} adagrasib)")
print(f"  positions covered: {feng['position'].min()} - {feng['position'].max()}")
print(f"  unique AA changes: {feng['aa_change'].nunique()}")

# %%
# 2. Codon-locked class-forcing (replaces AA-grammar only).
# The Feng library has G12C baked in, so position 12 codon = TGT not GGT.
rows = []
for _, r in feng.iterrows():
    pos = int(r["position"])
    wt_codon = KRAS_CODONS.get(pos)
    info = classify_amino_acid_change(r["wt_aa"], r["mut_aa"], wt_codon=wt_codon)
    rows.append({
        **r.to_dict(),
        "wt_codon": wt_codon,
        "reachable_1nt": info["reachable_1nt"],
        "class_forced": info["class_forced"],
        "forced_class": info["forced_class"],
        "n_classes": len(info["classes"]),
    })
classified = pd.DataFrame(rows)
print(f"\nCodon-locked classification of Feng atlas (KRAS NM_033360.4):")
print(f"  Reachable in 1 nt (SNV):      {classified['reachable_1nt'].sum()} / {len(classified)} "
      f"({classified['reachable_1nt'].mean()*100:.1f}%)")
print(f"  Class-forced (1 class only):  {classified['class_forced'].sum()} "
      f"({classified['class_forced'].mean()*100:.1f}%)")
print(f"\n  Forced-class breakdown:")
print(classified[classified["class_forced"]]["forced_class"].value_counts().to_string())

# %%
# 3. For each (mutagen × drug) compute P_simulated and P_chemfree.
#    Keep ONLY class-forced variants with Z-scores (drug-selection known).
cf = classified[classified["class_forced"]].copy()

all_results = []
for drug in DRUGS:
    z_col = f"{drug}_z"
    for mutagen_name, spectrum in MUTAGENS.items():
        sub = cf[cf[z_col].notna()].copy()
        sub["mutagen"] = mutagen_name
        sub["drug"] = drug
        sub["S_mutagen"] = sub["forced_class"].map(spectrum.values)
        sub["drug_selec"] = sub[z_col]  # Feng Z IS the drug-selection signal
        sub["log_P_simulated"] = np.log(sub["S_mutagen"]) + sub["drug_selec"]
        sub["log_P_chemfree"]  = sub["drug_selec"]
        sub["rank_simulated"] = sub["log_P_simulated"].rank(ascending=False, method="min")
        sub["rank_chemfree"]  = sub["log_P_chemfree"].rank(ascending=False, method="min")
        sub["rank_shift"] = sub["rank_simulated"] - sub["rank_chemfree"]
        all_results.append(sub)
results = pd.concat(all_results, ignore_index=True)
print(f"\nTotal (mutagen × drug) rows: {len(results)}")

# %%
# 4. Rank correlation per combination.
print("\n=== Spearman ρ (simulated vs chemistry-free) by (mutagen × drug) ===")
corr_rows = []
for drug in DRUGS:
    for mutagen in MUTAGENS:
        sub = results[(results["mutagen"] == mutagen) & (results["drug"] == drug)]
        rho = sub[["log_P_simulated", "log_P_chemfree"]].corr(method="spearman").iloc[0, 1]
        pear = sub[["log_P_simulated", "log_P_chemfree"]].corr(method="pearson").iloc[0, 1]
        corr_rows.append({"mutagen": mutagen, "drug": drug, "n": len(sub),
                          "spearman_rho": rho, "pearson_r": pear})
        print(f"  {mutagen:10s} × {drug:10s}: N={len(sub):4d}  ρ={rho:+.3f}  r={pear:+.3f}")
pd.DataFrame(corr_rows).to_csv(
    PROJECT_ROOT / "data" / "derived" / "feng_chemistry_impact_by_mutagen_drug.csv", index=False
)

# %%
# 5. Spot check: Awad 2021 clinical variants, where do they rank?
print("\n\n=== CLINICAL VARIANTS (Awad 2021 NEJM): ranks per (mutagen × drug) ===")
awad_variants = ["R68L", "R68S", "R68M", "H95D", "H95Q", "H95R",
                 "Y96D", "Y96C", "Y96H", "V8L", "V8M", "C12F"]
clinical_rows = []
for v in awad_variants:
    for drug in DRUGS:
        for mutagen in MUTAGENS:
            sub = results[(results["mutagen"] == mutagen) & (results["drug"] == drug)]
            row = sub[sub["aa_change"] == v]
            if len(row) == 0:
                clinical_rows.append({"variant": v, "drug": drug, "mutagen": mutagen,
                                      "status": "NOT CLASS-FORCED"})
            else:
                r = row.iloc[0]
                clinical_rows.append({
                    "variant": v, "drug": drug, "mutagen": mutagen,
                    "forced_class": r["forced_class"],
                    "S_mutagen": r["S_mutagen"],
                    "drug_z": r["drug_selec"],
                    "rank_simulated": int(r["rank_simulated"]),
                    "rank_chemfree":  int(r["rank_chemfree"]),
                    "rank_shift": int(r["rank_shift"]),
                })
cl = pd.DataFrame(clinical_rows)
print(cl.to_string(index=False))
cl.to_csv(PROJECT_ROOT / "data" / "derived" / "feng_awad2021_variant_ranks.csv", index=False)

# %%
# 6. Chemistry-under-counted predictions (top-50 fit × highest rank_shift).
#    Interpretation: "a KRAS G12C preclinical screen using Mutagen M against
#    Drug D will miss these variants, even though Feng says they are fit."
print("\n\n=== CHEMISTRY-UNDER-COUNTED PREDICTIONS per (mutagen × drug) ===")
all_predictions = []
for drug in DRUGS:
    for mutagen in MUTAGENS:
        sub = results[(results["mutagen"] == mutagen) & (results["drug"] == drug)].copy()
        top_fit = sub[sub["rank_chemfree"] <= 50]
        top_fit = top_fit.sort_values("rank_shift", ascending=False).head(15)
        if len(top_fit) == 0:
            continue
        print(f"\n  {mutagen} × {drug} — top-50-fit variants, ranked by chemistry under-counting:")
        print(top_fit[["aa_change", "forced_class", "S_mutagen", "drug_selec",
                       "rank_chemfree", "rank_simulated", "rank_shift"]].to_string(index=False))
        all_predictions.append(top_fit)

pd.concat(all_predictions, ignore_index=True).to_csv(
    PROJECT_ROOT / "data" / "derived" / "feng_chemistry_undercount_predictions.csv", index=False
)

# %%
# 7. Save the full per-variant results.
results.to_csv(PROJECT_ROOT / "data" / "derived" / "feng_decomposition.csv", index=False)
print(f"\nWrote data/derived/feng_decomposition.csv")
print(f"Wrote data/derived/feng_chemistry_impact_by_mutagen_drug.csv")
print(f"Wrote data/derived/feng_chemistry_undercount_predictions.csv")
print(f"Wrote data/derived/feng_awad2021_variant_ranks.csv")
