#!/usr/bin/env python3
"""Recompute Top-K recovery curves with 95% Wilson CIs on recall@K.

Writes a replacement data/derived/topk_recovery.csv with extra columns
`ci_lo` and `ci_hi`. Also writes recall@K CIs for the three fixed K values
quoted in the abstract (K=10, 20, 50).

Covers all four mutagens: ENU, Mutazyme II, CBE+ABE, Twist/PE. The ENU and
Mutazyme II rankings are read from the per-mutagen decomposition CSVs;
CBE+ABE and Twist/PE are computed in-place from the raw drug-selection
scores (and the same class-forced universe), because the earlier
decomposition notebook only handles the two SNV mutagens.

Mathematical model.  recall@K = k/K where k is the number of chemistry-free
top-K variants that are also in the chemistry-biased top-K.  This is a
proportion over K trials, so the 95 % Wilson score interval is the right
finite-sample uncertainty bound (R9 of revision plan).

NB: a bootstrap over variants is not appropriate: the underlying rank
statistic itself (sigma, S*sigma) is deterministic given the atlas; the
only sampling uncertainty is in the finite K draws from the chemistry-free
priority set, which the Wilson interval captures exactly.
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scipy.stats import norm

from src.harmonized_cohort import classify_amino_acid_change, load_alk_codons
from src.spectrum import load_spectra

Z = norm.ppf(0.975)

CLASSES = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")


def wilson_ci(k: int, n: int, z: float = Z) -> tuple[float, float]:
    """Two-sided Wilson score interval for binomial proportion k/n."""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    denom = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / denom
    half = (z / denom) * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))
    return max(0.0, centre - half), min(1.0, centre + half)


def compute_curve(df: pd.DataFrame, K_values: list[int],
                  sim_col: str = "rank_simulated",
                  cf_col: str = "rank_chemfree") -> list[dict]:
    """For a single (gene, drug, mutagen) atlas, compute recovery and CIs at each K."""
    df = df.copy().reset_index(drop=True)
    rows = []
    N = len(df)
    for K in K_values:
        K = min(K, N)
        sim = set(df.nsmallest(K, sim_col).index)
        cf  = set(df.nsmallest(K, cf_col).index)
        k = len(sim & cf)
        recov = k / K
        lo, hi = wilson_ci(k, K)
        rows.append({"K": K, "recovery": recov, "k": k, "ci_lo": lo, "ci_hi": hi})
    return rows


MUT_RENAME = {"MutazymeII": "MutII", "CBE_ABE": "BE", "Twist_PE": "Twist", "ENU": "ENU"}


def load_oh_class_forced_scored() -> pd.DataFrame:
    """ALK atlas as a class-forced long table with all four mutagens.

    Returns a DataFrame with columns:
      drug, mutagen, aa_change, residue, forced_class, S_mutagen,
      drug_selec, log_P_simulated, log_P_chemfree,
      rank_simulated, rank_chemfree
    where mutagen in {ENU, MutII, BE, Twist} per class-forced variant.
    """
    SPECTRA = load_spectra(ROOT / "data" / "reference" / "spectra.yaml")
    # One source of truth for the ALK codon table (see recompute_benchmark_r3).
    ALK_CODONS = load_alk_codons(ROOT / "data" / "reference" / "codons_alk_full.yaml")
    MUTAGENS = {
        "ENU":   dict(SPECTRA["ENU_Arnold2012_nominal"].values),
        "MutII": dict(SPECTRA["MutazymeII_DoiSupp2022"].values),
        "BE":    {"C>A": 0.0, "C>G": 0.0, "C>T": 0.5, "T>A": 0.0,
                  "T>C": 0.5, "T>G": 0.0},
        "Twist": {c: 1.0 / 6.0 for c in CLASSES},
    }
    oh = pd.read_csv(ROOT / "data" / "raw" / "oh_2026" / "table_s2_resistance_scores.csv")
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
    DRUGS = {"alectinib": "alec_score", "lorlatinib": "lor_score", "zotizalkib": "zot_score"}

    rows = []
    for _, r in oh_aa.iterrows():
        info = classify_amino_acid_change(
            r["wt_aa"], r["mut_aa"], wt_codon=ALK_CODONS.get(int(r["residue"]))
        )
        if not info["class_forced"]:
            continue
        for mg, spec in MUTAGENS.items():
            rows.append({
                "mutagen": mg, "aa_change": r["aa_change"],
                "residue": int(r["residue"]),
                "forced_class": info["forced_class"],
                "S_mutagen": spec[info["forced_class"]],
                "alec_score": r["alec_score"],
                "lor_score": r["lor_score"],
                "zot_score": r["zot_score"],
                "fitness_score": r["fitness_score"],
            })
    annotated = pd.DataFrame(rows)
    out = []
    for drug, score_col in DRUGS.items():
        for mg in MUTAGENS:
            sub = annotated[annotated["mutagen"] == mg].copy()
            # σ_drug(v) = Resistance(v) per Oh 2026 Methods (see notebook/05 header
            # comment). The Resistance score is already a same-timepoint drug/vehicle
            # contrast; subtracting Fitness would double-correct.
            sub["drug_selec"] = sub[score_col]
            sub["log_P_simulated"] = np.where(
                sub["S_mutagen"] > 0,
                np.log(sub["S_mutagen"].clip(lower=1e-12)) + sub["drug_selec"],
                -1e18,  # sentinel: zero supply pushes out of top-K
            )
            sub["log_P_chemfree"] = sub["drug_selec"]
            sub["rank_simulated"] = sub["log_P_simulated"].rank(ascending=False, method="min")
            sub["rank_chemfree"]  = sub["log_P_chemfree"].rank(ascending=False, method="min")
            sub["drug"] = drug
            sub["gene"] = "ALK"
            out.append(sub)
    return pd.concat(out, ignore_index=True)


def load_feng_class_forced_scored() -> pd.DataFrame:
    """KRAS G12C atlas as a class-forced long table with all four mutagens."""
    SPECTRA = load_spectra(ROOT / "data" / "reference" / "spectra.yaml")
    MUTAGENS = {
        "ENU":   dict(SPECTRA["ENU_Arnold2012_nominal"].values),
        "MutII": dict(SPECTRA["MutazymeII_DoiSupp2022"].values),
        "BE":    {"C>A": 0.0, "C>G": 0.0, "C>T": 0.5, "T>A": 0.0,
                  "T>C": 0.5, "T>G": 0.0},
        "Twist": {c: 1.0 / 6.0 for c in CLASSES},
    }
    feng = pd.read_csv(ROOT / "data" / "derived" / "feng_decomposition.csv")
    cf = feng[feng["class_forced"] == True].drop_duplicates(["aa_change", "drug", "mutagen"], keep="first")
    # The ENU+MutII rows are already there. We add BE and Twist for every (drug, aa_change).
    base_cols = ["position", "wt_aa", "mut_aa", "sotorasib_z", "adagrasib_z",
                 "aa_change", "forced_class", "drug"]
    # unique variants per drug (already class-forced)
    uniq = cf[base_cols].drop_duplicates().reset_index(drop=True)
    rows = []
    for drug, sel_col in [("sotorasib", "sotorasib_z"), ("adagrasib", "adagrasib_z")]:
        drug_uniq = uniq[uniq["drug"] == drug].copy()
        for mg, spec in MUTAGENS.items():
            sub = drug_uniq.copy()
            sub["mutagen"] = mg
            sub["S_mutagen"] = sub["forced_class"].map(spec)
            sub["drug_selec"] = sub[sel_col]
            sub["log_P_simulated"] = np.where(
                sub["S_mutagen"] > 0,
                np.log(sub["S_mutagen"].clip(lower=1e-12)) + sub["drug_selec"],
                -1e18,
            )
            sub["log_P_chemfree"] = sub["drug_selec"]
            sub["rank_simulated"] = sub["log_P_simulated"].rank(ascending=False, method="min")
            sub["rank_chemfree"]  = sub["log_P_chemfree"].rank(ascending=False, method="min")
            sub["gene"] = "KRAS G12C"
            rows.append(sub)
    return pd.concat(rows, ignore_index=True)


def main():
    oh = load_oh_class_forced_scored()
    feng = load_feng_class_forced_scored()
    both = pd.concat([
        oh[["gene", "drug", "mutagen", "aa_change", "rank_simulated", "rank_chemfree"]],
        feng[["gene", "drug", "mutagen", "aa_change", "rank_simulated", "rank_chemfree"]],
    ], ignore_index=True)

    out_rows = []
    for (gene, drug), group in both.groupby(["gene", "drug"]):
        N = int(group.groupby("mutagen").size().max())
        # Dense K sampling for the figure band
        K_values = sorted(set(
            list(range(1, min(50, N) + 1))
            + list(range(50, min(200, N) + 1, 5))
            + list(range(200, N + 1, max(1, (N - 200) // 60)))
            + [N]
        ))
        for mg, sub in group.groupby("mutagen"):
            sub = sub.dropna(subset=["rank_simulated", "rank_chemfree"]).reset_index(drop=True)
            if sub.empty:
                continue
            curve = compute_curve(sub, K_values)
            for row in curve:
                out_rows.append({
                    "gene": gene, "drug": drug, "mutagen": mg,
                    "K": row["K"], "recovery": row["recovery"],
                    "ci_lo": row["ci_lo"], "ci_hi": row["ci_hi"],
                    "N": len(sub),
                })
    out = pd.DataFrame(out_rows)
    out_path = ROOT / "data" / "derived" / "topk_recovery.csv"
    bak = out_path.with_suffix(".csv.pre_r9")
    if out_path.exists() and not bak.exists():
        bak.write_bytes(out_path.read_bytes())
    out.to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({len(out)} rows, with ci_lo/ci_hi)")

    # Also emit abstract-quoted CIs at K=10, 20, 50, 100 to a summary CSV
    summary = out[out["K"].isin([10, 20, 50, 100])].copy()
    summary.to_csv(ROOT / "data" / "derived" / "topk_recovery_abstract_CIs.csv", index=False)
    print(f"Wrote abstract CI summary: {len(summary)} rows at K in {{10, 20, 50, 100}}")

    # Diagnostic: V1180L and L1196M ranks under Mutazyme / chemfree for alec
    alec_mut = oh[(oh.drug == "alectinib") & (oh.mutagen == "MutII")]
    for v in ("V1180L", "L1196M", "C1156Y", "L1152R", "I1171N", "G1202R"):
        rows = alec_mut[alec_mut.aa_change == v]
        if not rows.empty:
            r = rows.iloc[0]
            print(f"  alec Mut {v}: rank_cf={r.rank_chemfree:.0f} / {len(alec_mut)} "
                  f"rank_sim={r.rank_simulated:.0f}")


if __name__ == "__main__":
    main()
