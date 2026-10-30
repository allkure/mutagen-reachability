# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
# ---

# %% [markdown]
# # KRAS G12C summary stats for the paper (codon-locked)
#
# With the KRAS-4A NM_033360.4 codon table (data/reference/codons_kras.yaml,
# positions 1-150 ClinVar-verified at anchor residues; HVR 151-189 consensus),
# we can split the Awad 2021 clinical panel into three categories:
#
#   A. NOT SNV-REACHABLE from the KRAS G12C WT codon: variants requiring 2+ nt
#      changes. The Feng-Twist (synthetic DNA) library captures them, but
#      NO SNV mutagen (ENU, Mutazyme, EMS, ...) can produce them from the
#      G12C-mutant parent plasmid.
#   B. MULTI-CLASS: reachable in 1 nt but via multiple pyrimidine classes.
#      Chemistry acts, but the forced-class attribution is ambiguous.
#   C. CLASS-FORCED: reachable in 1 nt via a single pyrimidine class. These
#      are the clean cases where we can compute P̂_simulated(v) and
#      compare it to Feng's chemistry-free P̂_chemfree(v).

# %%
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

NOTEBOOK_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = NOTEBOOK_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.harmonized_cohort import classify_amino_acid_change

# %%
# Load the full decomposition (contains only class-forced variants).
results = pd.read_csv(PROJECT_ROOT / "data" / "derived" / "feng_decomposition.csv")
feng = pd.read_csv(PROJECT_ROOT / "data" / "raw" / "feng_2022" / "table_s1_zscores.csv")
feng["aa_change"] = feng["wt_aa"] + feng["position"].astype(str) + feng["mut_aa"]

# Load KRAS codons
import yaml
with open(PROJECT_ROOT / "data" / "reference" / "codons_kras.yaml") as f:
    raw = yaml.safe_load(f)
KRAS_CODONS = {int(r): v["codon"] for r, v in raw["kras_4a_nm_033360_4"]["codons"].items()}
KRAS_CODONS[12] = "TGT"  # G12C library

# %%
# Awad 2021 clinical panel (NEJM 384:2382, Table 1 + Fig 2)
# Sotorasib-emergent: V8L, V8M, Y96D, Y96H (+ shared)
# Adagrasib-emergent: H95D, H95Q, H95R, R68S, R68M (+ shared)
# Shared/either: C12F, R68L, Y96C
AWAD_SHARED = ["C12F", "R68L", "Y96C"]
AWAD_SOT    = ["V8L", "V8M", "Y96D", "Y96H"]
AWAD_ADA    = ["H95D", "H95Q", "H95R", "R68S", "R68M"]
AWAD_ALL = AWAD_SHARED + AWAD_SOT + AWAD_ADA

# %%
# Categorise each Awad variant
print("="*80)
print("TABLE A — Awad 2021 clinical panel: SNV reachability from KRAS WT codons")
print("="*80)
print()
print(f"{'Variant':10s}  {'WT codon':9s}  {'Reachable_1nt':15s}  {'Class_forced':14s}  {'Classes':10s}  Comment")
reach_by_var = {}
for v in AWAD_ALL:
    wt_aa = v[0]
    pos = int(v[1:-1])
    mut_aa = v[-1]
    wt_codon = KRAS_CODONS.get(pos, "?")
    info = classify_amino_acid_change(wt_aa, mut_aa, wt_codon=wt_codon)
    reach_by_var[v] = info
    cls = ",".join(info["classes"]) if info["classes"] else "-"
    cmt = ""
    if not info["reachable_1nt"]:
        cmt = "requires 2+ nt — NO SNV mutagen reaches this from WT"
    elif not info["class_forced"]:
        cmt = f"multi-class ({len(info['classes'])} classes)"
    else:
        cmt = "clean case"
    print(f"  {v:8s}  {wt_codon:9s}  {str(info['reachable_1nt']):15s}  "
          f"{str(info['class_forced']):14s}  {cls:10s}  {cmt}")

n_unreach = sum(1 for v in AWAD_ALL if not reach_by_var[v]["reachable_1nt"])
n_multi   = sum(1 for v in AWAD_ALL if reach_by_var[v]["reachable_1nt"] and not reach_by_var[v]["class_forced"])
n_forced  = sum(1 for v in AWAD_ALL if reach_by_var[v]["class_forced"])
print(f"\nAwad panel breakdown: {n_unreach} not-SNV-reachable, {n_multi} multi-class, "
      f"{n_forced} class-forced  (total {len(AWAD_ALL)})")

# %%
print("\n" + "="*80)
print("TABLE B — Chemistry effect on class-forced clinical variants")
print("="*80)
print()
print(f"{'Variant':8s}  {'Drug':10s}  {'Class':6s}  {'Z':>6s}  "
      f"{'rk_free':>8s}  {'MutII_rk':>9s}  {'ENU_rk':>7s}  {'MutII_sh':>9s}  {'ENU_sh':>7s}")
for v in AWAD_ALL:
    if not reach_by_var[v]["class_forced"]:
        continue
    for drug in ("sotorasib", "adagrasib"):
        mr = results[(results["mutagen"] == "MutazymeII") & (results["drug"] == drug) &
                     (results["aa_change"] == v)]
        er = results[(results["mutagen"] == "ENU") & (results["drug"] == drug) &
                     (results["aa_change"] == v)]
        if len(mr) == 0 or len(er) == 0:
            continue
        m, e = mr.iloc[0], er.iloc[0]
        print(f"  {v:6s}  {drug:10s}  {m['forced_class']:6s}  "
              f"{m['drug_selec']:+6.2f}  "
              f"{int(m['rank_chemfree']):>8d}  "
              f"{int(m['rank_simulated']):>9d}  "
              f"{int(e['rank_simulated']):>7d}  "
              f"{int(m['rank_shift']):>+9d}  "
              f"{int(e['rank_shift']):>+7d}")

# %%
print("\n" + "="*80)
print("TABLE C — Headline predictions (clinical variants chemistry under-counts)")
print("="*80)
print()
print("Falsifiable claim per entry: a Mutagen-M-based KRAS G12C preclinical")
print("screen will MISS this clinical variant in a top-N ranking.")
print()
# For each clinical class-forced variant, show the MOST severe shift under
# a specific chemistry (|shift| > 100 AND rank_free in top-50).
for v in AWAD_ALL:
    if not reach_by_var[v]["class_forced"]:
        continue
    best_hits = []
    for drug in ("sotorasib", "adagrasib"):
        for mutagen in ("MutazymeII", "ENU"):
            sub = results[(results["drug"] == drug) & (results["mutagen"] == mutagen) &
                          (results["aa_change"] == v)]
            if len(sub) == 0:
                continue
            r = sub.iloc[0]
            if r["rank_chemfree"] <= 50 and r["rank_shift"] > 100:
                best_hits.append((drug, mutagen, int(r["rank_chemfree"]),
                                 int(r["rank_simulated"]), int(r["rank_shift"])))
    for h in best_hits:
        d, m, rf, rs, sh = h
        print(f"  {v:5s} × {d:10s} × {m:11s}  rank_free={rf:3d}  rank_sim={rs:4d}  shift=+{sh:4d}")

# %%
print("\n" + "="*80)
print("TABLE D — The 'SNV-unreachable' clinical variants (chemistry-INDEPENDENT undercount)")
print("="*80)
print()
print("Variants no SNV mutagen can produce from WT KRAS G12C (regardless of spectrum).")
print("These variants are visible in Feng-Twist because Twist is a synthetic library,")
print("but are STRUCTURALLY INVISIBLE to any ENU, Mutazyme, EMS, or similar screen.")
print()
for v in AWAD_ALL:
    if reach_by_var[v]["reachable_1nt"]:
        continue
    wt_aa = v[0]
    pos = int(v[1:-1])
    mut_aa = v[-1]
    wt_codon = KRAS_CODONS.get(pos, "?")
    sot_z = feng[feng["aa_change"] == v]["sotorasib_z"].values
    ada_z = feng[feng["aa_change"] == v]["adagrasib_z"].values
    sot_s = f"{sot_z[0]:+5.2f}" if len(sot_z) and pd.notna(sot_z[0]) else "   NA"
    ada_s = f"{ada_z[0]:+5.2f}" if len(ada_z) and pd.notna(ada_z[0]) else "   NA"
    print(f"  {v:6s}  WT codon {wt_codon} → {mut_aa} needs 2+ nt  "
          f"Feng-Z sot={sot_s}  ada={ada_s}")

print("\n(This category is a stronger statement than chemistry under-counting:")
print("it is a structural prediction that no mutagen chemistry can repair.)")
