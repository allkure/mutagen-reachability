# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
# ---

# %% [markdown]
# # Compound-mutation predictions for current 4th-gen ALK TKI preclinical screens
#
# For each published preclinical compound screen in the literature, enumerate
# the specific second-site (V2) mutations that will be chemistry-under-counted
# by the mutagen used, independent of their fitness. For drugs with Oh 2026 σ
# available (alectinib, lorlatinib, zotizalkib), the predictions are further
# ranked by Oh-measured fitness so the top candidates are the ones that
# matter: fit per Oh atlas, chemistry-rare per mutagen.
#
# **Setup for compound screens.** In a preclinical compound screen with fixed
# V1 backbone and chemistry-generated V2:
#     S_screen(V1 + V2) = S_mutagen(V2)      [V1 fixed by library, no chemistry cost]
#
# So the chemistry-undercount prediction for V2 is identical to the single-
# mutation prediction we computed in notebook/05. This script re-casts those
# predictions in compound-mutation form for each published screen.
#
# **Caveat on fitness.** For compound fitness σ(V1 + V2), we use Oh's single-
# mutation σ(V2) as a proxy. This is NOT equivalent to measured compound
# fitness and ignores epistasis — e.g. G1202R+L1196M has super-additive
# lorlatinib resistance that would not be captured. The chemistry prediction
# is EXACT; the fitness interpretation is approximate. For drugs not in Oh
# 2026 (repotrectinib, ensartinib, gilteritinib, neladalkib) we fall back to
# chemistry-only predictions, which are weaker but still useful.
#
# **Output:** `data/derived/compound_predictions_by_screen.csv` with one row
# per (screen × V1 backbone × V2 candidate) prediction, ranked by chemistry
# under-count magnitude.

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

# %%
# 1. Published preclinical compound screens of ALK TKIs.
# For each: (study_tag, mutagen_name, V1_backbone, drug, drug_generation,
#            oh_drug_column, proxy_note, data_source_citation)
#
# σ proxy choices (used when the screened drug is not in Oh 2026's set
# {alectinib, lorlatinib, zotizalkib}):
#   - ensartinib → alectinib proxy: both 2nd-gen ATP-pocket binders, overlapping
#     resistance variants (L1196M, G1202R, I1171N, F1174)
#   - repotrectinib → zotizalkib proxy: both 4th-gen designed to overcome
#     compound mutations, similar binding conformation (shared target design logic)
#   - neladalkib (NVL-655) → zotizalkib proxy: both 4th-gen ALK-selective
#     designed for post-lorlatinib resistance landscape
#   - gilteritinib → NO proxy: FLT3 inhibitor repurposed off-label for I1171N+
#     compound screens; binding mode is NOT ALK-optimized, so Oh ALK σ is
#     scientifically inappropriate as a proxy. Chemistry-only prediction with
#     an explicit flag that fitness is unknown.
#
# These proxy choices are approximations; the predictions carry a
# `proxy_note` column so downstream users can filter or re-weight.
SCREENS = [
    ("Doi 2022 Fig 3b",  "MutazymeII", "G1202R", "repotrectinib", "4th gen", "zot_score",
     "proxy: zotizalkib (both 4th-gen compound-targeting ALK-TKIs)",
     "Doi Y et al. Target Oncol 2022;17:695-707 Fig 3b"),
    ("Doi 2022 Fig 3c",  "MutazymeII", "I1171N", "ensartinib",    "2nd gen", "alec_score",
     "proxy: alectinib (both 2nd-gen ATP-pocket ALK-TKIs)",
     "Doi Y et al. Target Oncol 2022;17:695-707 Fig 3c"),
    ("Takei 2025 zot G1202R+", "MutazymeII", "G1202R", "zotizalkib",  "4th gen", "zot_score",
     "direct σ from Oh 2026",
     "Takei Y et al. bioRxiv 2025.05.26.656085 (Semba lab)"),
    ("Takei 2025 gil I1171N+", "MutazymeII", "I1171N", "gilteritinib","4th gen (off-target)", None,
     "no proxy: FLT3-first drug, off-label on ALK I1171N+ screen",
     "Takei Y et al. bioRxiv 2025.05.26.656085 (Semba lab)"),
    ("Takei 2025 nel G1202R+", "MutazymeII", "G1202R", "neladalkib",  "4th gen", "zot_score",
     "proxy: zotizalkib (both 4th-gen ALK-selective)",
     "Takei Y et al. bioRxiv 2025.05.26.656085 (Semba lab)"),
    ("Takei 2025 nel I1171N+", "MutazymeII", "I1171N", "neladalkib",  "4th gen", "zot_score",
     "proxy: zotizalkib (both 4th-gen ALK-selective)",
     "Takei Y et al. bioRxiv 2025.05.26.656085 (Semba lab)"),
]

# %%
# 2. Load Oh 2026 AA-level resistance scores (as in notebook/05).
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
).reset_index()
print(f"Oh atlas AA-level ALK KD: {len(oh_aa)} variants")

# %%
# 3. For each ALK KD missense variant, compute class-forced S_mutagen for
# each mutagen we care about.
MUTAGENS_USED = {"MutazymeII": SPECTRA["MutazymeII_DoiSupp2022"]}
# Also include ENU as a counter-factual: "what if someone had used ENU
# instead of Mutazyme II on these backbones?" so the paper can also
# recommend complementary-chemistry screens.
MUTAGENS_USED["ENU"] = SPECTRA["ENU_Arnold2012_nominal"]

def annotate_chem(row, spectrum):
    info = classify_amino_acid_change(
        row["wt_aa"], row["mut_aa"], wt_codon=ALK_CODONS.get(int(row["residue"]))
    )
    return pd.Series({
        "class_forced": info["class_forced"],
        "forced_class": info["forced_class"],
        "S_mut": spectrum[info["forced_class"]] if info["class_forced"] else None,
    })

chem_by_mutagen = {}
for name, sp in MUTAGENS_USED.items():
    tmp = oh_aa.copy()
    tmp = pd.concat([tmp, tmp.apply(lambda r: annotate_chem(r, sp), axis=1)], axis=1)
    tmp = tmp[tmp["class_forced"] == True].copy()
    chem_by_mutagen[name] = tmp

# %%
# 4. Build the compound-prediction table.
all_rows = []
for study_tag, mutagen, v1, drug, gen, oh_col, proxy_note, citation in SCREENS:
    chem = chem_by_mutagen[mutagen].copy()
    chem["study"] = study_tag
    chem["mutagen"] = mutagen
    chem["v1_backbone"] = v1
    chem["drug"] = drug
    chem["drug_gen"] = gen
    chem["citation"] = citation
    chem["sigma_source"] = proxy_note

    # Drop V2 variants that are at the same residue as V1 (can't compound
    # with themselves and typically the preclinical library excludes them).
    try:
        v1_res = int(''.join(c for c in v1 if c.isdigit()))
        chem = chem[chem["residue"] != v1_res].copy()
    except Exception:
        pass

    # Chemistry score: higher = more chemistry supply for V2 → more likely to arise
    # under this mutagen. Pure chemistry rank (independent of fitness).
    chem["chem_rank"] = chem["S_mut"].rank(ascending=False, method="min")

    if oh_col is not None:
        # Drug σ available: compute Oh-anchored predictions.
        chem["drug_selec_Oh"] = chem[oh_col] - chem["fitness_score"]
        chem["log_P_sim"] = np.log(chem["S_mut"]) + chem["drug_selec_Oh"]
        chem["log_P_chemfree"] = chem["drug_selec_Oh"]
        chem["rank_simulated"] = chem["log_P_sim"].rank(ascending=False, method="min")
        chem["rank_chemfree"]  = chem["log_P_chemfree"].rank(ascending=False, method="min")
        chem["under_count_shift"] = chem["rank_simulated"] - chem["rank_chemfree"]
        chem["oh_sigma_available"] = True
        # Compound-mutation prediction: top-N most under-counted among top-50 Oh-fit
        top_fit = chem[chem["rank_chemfree"] <= 50].copy()
        preds = top_fit.sort_values("under_count_shift", ascending=False).head(15)
        preds["compound_variant"] = v1 + " + " + preds["aa_change"]
    else:
        # No Oh σ for this drug — pure chemistry prediction only.
        chem["drug_selec_Oh"] = None
        chem["log_P_sim"] = None
        chem["log_P_chemfree"] = None
        chem["rank_simulated"] = None
        chem["rank_chemfree"] = None
        chem["under_count_shift"] = None
        chem["oh_sigma_available"] = False
        # Chemistry-only prediction: V2 variants with lowest S_mut per class.
        # These will be systematically under-represented regardless of fitness.
        # To make the list specific and short, take the top-15 chemistry-poorest
        # variants in classes that fall below the mutagen's median class rate.
        median_class_rate = np.median(list(MUTAGENS_USED[mutagen].values.values()))
        chem_poor = chem[chem["S_mut"] < median_class_rate].copy()
        preds = chem_poor.sort_values("S_mut").head(15)
        preds["compound_variant"] = v1 + " + " + preds["aa_change"]

    preds["note"] = (
        "Chemistry under-count: fit per Oh σ (top-50 chemistry-free), S_mut low under this mutagen."
        if oh_col is not None else
        "Chemistry-only prediction (no Oh σ for this drug): V2 in a mutagen-poor substitution class."
    )
    all_rows.append(preds)

out = pd.concat(all_rows, ignore_index=True)
# Keep only the columns that make sense for the deliverable table.
cols = [
    "study", "mutagen", "v1_backbone", "drug", "drug_gen", "compound_variant",
    "forced_class", "S_mut", "oh_sigma_available", "sigma_source",
    "drug_selec_Oh", "rank_chemfree", "rank_simulated", "under_count_shift",
    "note", "citation",
]
out = out[cols]
out_path = PROJECT_ROOT / "data" / "derived" / "compound_predictions_by_screen.csv"
out.to_csv(out_path, index=False)
print(f"\nWrote {out_path}")
print(f"Total compound predictions: {len(out)} rows across {out['study'].nunique()} screens")

# %%
# 5. Compact view per screen.
print("\n" + "="*90)
print("COMPOUND-MUTATION CHEMISTRY-UNDERCOUNT PREDICTIONS")
print("="*90)
for (study, v1, drug, gen), g in out.groupby(["study", "v1_backbone", "drug", "drug_gen"], sort=False):
    print(f"\n{study} ({drug}, {gen}, mutagen={g['mutagen'].iloc[0]}, V1={v1}):")
    if g["oh_sigma_available"].iloc[0]:
        print(f"  Oh σ available → chemistry-undercount predictions (fit per Oh, chemistry-rare):")
        display_cols = ["compound_variant", "forced_class", "S_mut",
                        "drug_selec_Oh", "rank_chemfree", "rank_simulated", "under_count_shift"]
    else:
        print(f"  No Oh σ for this drug → chemistry-only predictions (mutagen-poor V2):")
        display_cols = ["compound_variant", "forced_class", "S_mut"]
    print(g[display_cols].to_string(index=False))
