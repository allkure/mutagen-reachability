#!/usr/bin/env python3
"""R3 benchmark + R10 proxy sensitivity.

Produces data/derived/benchmark_r3.csv — recall@K and AUROC for σ-alone,
log S + σ (framework at β=1), and log S alone, on four retrospective anchors.
Permutation p-values (5,000 class-label shuffles) for the recall@100 gain of
the framework over σ-alone.

R10 (proxy sensitivity): for the industry ENU anchor the Oh atlas does not
include neladalkib, so σ is proxied by alectinib's σ in the headline.  Here
we rerun the same anchor with lorlatinib σ as an alternative proxy and emit
both results so §2.5 can state whether conclusions are proxy-sensitive.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.harmonized_cohort import classify_amino_acid_change, load_alk_codons
from src.spectrum import load_spectra

RNG = np.random.default_rng(20261013)
CLASSES = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

SPECTRA = load_spectra(ROOT / "data" / "reference" / "spectra.yaml")
MUT_MUT = dict(SPECTRA["MutazymeII_DoiSupp2022"].values)
MUT_ENU = dict(SPECTRA["ENU_Arnold2012_nominal"].values)

# Empirical in vitro ENU C>G rate inverted from Nuvalent anchor (R4)
ENU_EMP = dict(MUT_ENU)
# Keep the five non-C>G classes proportional, inflate C>G to posterior upper (0.23 median)
_old_cg = MUT_ENU["C>G"]
_new_cg = 0.23  # median of 7-50% Jeffreys posterior
_rest = {k: v for k, v in MUT_ENU.items() if k != "C>G"}
_rest_sum = sum(_rest.values())
_scale = (1 - _new_cg) / _rest_sum
ENU_EMP = {k: (_new_cg if k == "C>G" else v * _scale) for k, v in MUT_ENU.items()}

# One source of truth for the ALK codon table. This previously inlined the
# 23-codon curated slice, which silently gave this script a smaller
# class-forced universe than the decomposition notebook was using.
ALK_CODONS = load_alk_codons(ROOT / "data" / "reference" / "codons_alk_full.yaml")


def load_oh_scored() -> pd.DataFrame:
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
    rows = []
    for _, r in oh_aa.iterrows():
        info = classify_amino_acid_change(
            r["wt_aa"], r["mut_aa"], wt_codon=ALK_CODONS.get(int(r["residue"]))
        )
        if not info["class_forced"]:
            continue
        row = {"aa_change": r["aa_change"],
               "residue": int(r["residue"]),
               "forced_class": info["forced_class"],
               "alec_score": r["alec_score"], "lor_score": r["lor_score"],
               "zot_score": r["zot_score"], "fitness_score": r["fitness_score"]}
        rows.append(row)
    return pd.DataFrame(rows)


def load_feng_scored() -> pd.DataFrame:
    feng = pd.read_csv(ROOT / "data" / "derived" / "feng_decomposition.csv")
    cf = feng[feng["class_forced"] == True].drop_duplicates(["aa_change", "drug"], keep="first")
    return cf[["aa_change", "position", "forced_class", "sotorasib_z",
               "adagrasib_z", "drug"]]


def recall_at_k(ranks: pd.Series, hits: set[str], aa_index: pd.Series, k: int) -> tuple[float, int]:
    sub = ranks.reset_index(drop=True)
    aa_sub = aa_index.reset_index(drop=True)
    topk_aas = set(aa_sub.iloc[sub.nsmallest(k).index])
    recall = len(topk_aas & hits) / len(hits) if hits else 0.0
    return recall, len(topk_aas & hits)


def recall_count_at_k(ranks: np.ndarray, hits_mask: np.ndarray, k: int) -> int:
    """Return # of hits in the top-k entries (ranks are 1..N)."""
    order = np.argsort(ranks, kind="stable")[:k]
    return int(hits_mask[order].sum())


def auroc(scores: np.ndarray, hits_mask: np.ndarray) -> float:
    """ROC-AUC with scores (higher = more likely to be a hit)."""
    pos_scores = scores[hits_mask == 1]
    neg_scores = scores[hits_mask == 0]
    n_pos, n_neg = len(pos_scores), len(neg_scores)
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    all_scores = np.concatenate([pos_scores, neg_scores])
    ranks = pd.Series(all_scores).rank(ascending=True, method="average").values
    pos_ranks = ranks[:n_pos]
    u = pos_ranks.sum() - n_pos * (n_pos + 1) / 2
    return float(u / (n_pos * n_neg))


def permutation_p(sigma: np.ndarray, log_S: np.ndarray, forced_class: np.ndarray,
                  hits_mask: np.ndarray, obs_gain: float, k: int = 100,
                  n_iter: int = 5000) -> dict[str, float]:
    """Class-label permutation: shuffle forced_class across variants,
    recompute log_S from a fixed class→S mapping, recompute framework recall@K
    gain over σ-alone, return empirical two-sided p-value for the observed gain.

    The gain is expressed as (hits_captured_framework - hits_captured_sigma_alone)
    / n_hits — the SAME normalization used in the output rows, so obs_gain and
    simulated gains live on the same scale.

    Also returns the test's resolution, because these anchors are small:

      mdg      the smallest |gain| that would reach two-sided p < 0.05, read off
               the permutation null as its 97.5th percentile of |gain|.  A
               non-significant p above this threshold is uninformative rather
               than evidence of no effect.
      quantum  1 / n_hits.  Recall at these n can only move in whole hits, so
               the gain is quantized and no effect smaller than one hit is
               observable at all.
    """
    rng = np.random.default_rng(20261013)
    uniq_cls = np.unique(forced_class)
    cls_to_logS = {}
    for c in uniq_cls:
        vals = log_S[forced_class == c]
        if len(vals) > 0:
            cls_to_logS[c] = float(vals[0])
    gains = []
    n_hits = max(1, int(hits_mask.sum()))
    base_rank_sigma = pd.Series(-sigma).rank(method="min").values
    base_hits_sigma = recall_count_at_k(base_rank_sigma, hits_mask, k)
    for _ in range(n_iter):
        shuffled_cls = rng.permutation(forced_class)
        new_logS = np.array([cls_to_logS[c] for c in shuffled_cls])
        new_score = sigma + new_logS
        new_rank = pd.Series(-new_score).rank(method="min").values
        new_hits = recall_count_at_k(new_rank, hits_mask, k)
        gains.append((new_hits - base_hits_sigma) / n_hits)
    gains = np.array(gains)
    return {
        "p": float((np.sum(np.abs(gains) >= np.abs(obs_gain)) + 1) / (n_iter + 1)),
        "mdg": float(np.quantile(np.abs(gains), 0.975)),
        "quantum": 1.0 / n_hits,
    }


def run_anchor(name: str, N_df: pd.DataFrame, sigma: np.ndarray, S_mut: np.ndarray,
               hits: set[str], aa_index: np.ndarray,
               ranking_label: str, extra: dict | None = None) -> list[dict]:
    """Run σ-alone, log S + σ, log S alone for one anchor and emit rows.

    sigma, S_mut, aa_index: length-N arrays aligned to the atlas (class-forced).
    hits: set of AA-change strings.
    """
    hits_mask = np.array([1 if a in hits else 0 for a in aa_index], dtype=int)

    # Every named hit must be present in the ranked universe. The recall rows
    # below divide by len(hits) while the permutation null divides by
    # hits_mask.sum(); the two agree only when no named hit is missing. A hit
    # absent from the universe is an unreachable variant, which belongs in the
    # reachability inventory, not in a ranking denominator.
    missing = sorted(hits - set(aa_index))
    assert not missing, (
        f"{name}: named hits absent from the ranked universe: {missing}. "
        "Either they are not class-forced and belong in the two-step section, "
        "or the universe is wrong.")

    log_S = np.log(np.clip(S_mut, 1e-12, None))
    log_S = np.where(S_mut > 0, log_S, -1e18)

    rows = []
    for label, score in [
        ("σ alone",                 sigma),
        (f"Framework (log S + σ{', '+extra['note'] if extra and extra.get('note') else ''})",  sigma + log_S),
        (f"Supply alone (log S{', '+extra['note'] if extra and extra.get('note') else ''})",   log_S),
    ]:
        rank = pd.Series(-score).rank(method="min").values
        r10  = recall_count_at_k(rank, hits_mask, 10)  / max(1, len(hits))
        r50  = recall_count_at_k(rank, hits_mask, 50)  / max(1, len(hits))
        r100 = recall_count_at_k(rank, hits_mask, 100) / max(1, len(hits))
        au   = auroc(score, hits_mask)
        rows.append({
            "anchor": name, "ranking": label, "N": len(sigma),
            "n_hits": int(hits_mask.sum()),
            "recall@10": r10, "recall@50": r50, "recall@100": r100, "AUROC": au,
        })
    # permutation p-value for framework over σ-alone at K=100
    obs_gain = rows[1]["recall@100"] - rows[0]["recall@100"]
    forced = np.array([f.replace("$", "") for f in N_df["forced_class"].values])
    res = permutation_p(sigma, log_S, forced, hits_mask, obs_gain, k=100, n_iter=5000)
    rows[1]["perm_p_at_100"] = res["p"]
    rows[1]["min_detectable_gain"] = res["mdg"]
    rows[1]["gain_quantum"] = res["quantum"]
    rows[1]["observed_gain"] = obs_gain
    for i in (0, 2):
        rows[i]["perm_p_at_100"] = np.nan
        rows[i]["min_detectable_gain"] = np.nan
        rows[i]["gain_quantum"] = np.nan
        rows[i]["observed_gain"] = np.nan
    return rows


def main():
    oh = load_oh_scored()
    feng = load_feng_scored()

    # Hit sets
    DOI_HITS = {"V1180L", "L1196M", "C1156Y", "L1152R", "I1171N"}  # class-forced subset of Doi 2022
    # Nuvalent industry ENU: 11 class-forced (excluding R1192P outside atlas, G1202R multi-class)
    NUVALENT_HITS = {"I1171N", "I1171S", "I1171T", "F1174C", "F1174V", "V1180L",
                     "T1151K", "T1151R", "E1129K", "L1196M", "G1269A"}
    # Awad KRAS G12C, class-forced subset of the patient-observed hits.
    # R68L is a sotorasib hit but is NOT class-forced: AGG (Arg) is two
    # nucleotide changes from every Leu codon, so it has no one-step route and
    # cannot appear in a class-forced universe. It is reported in the two-step
    # section instead, with V8M. Including it here would put a hit in the
    # recall denominator that no ranking over this universe can ever retrieve,
    # which is a reachability statement and not a ranking one.
    AWAD_SOT_HITS = {"Y96D", "Y96H", "H95R"}
    AWAD_ADA_HITS = {"Y96D", "Y96H", "H95D", "H95R"}

    all_rows = []

    # σ_drug(v) = Resistance(v) per Oh 2026 Methods — see notebook/05 for rationale
    # ---- Anchor 1: Doi × alec × MutII
    oh_alec = oh.copy()
    oh_alec["sigma"] = oh_alec["alec_score"]
    oh_alec["S_mut"] = oh_alec["forced_class"].map(MUT_MUT)
    all_rows += run_anchor("Doi 2022 (ALK Mutazyme × alec)",
                           oh_alec, oh_alec["sigma"].values, oh_alec["S_mut"].values,
                           DOI_HITS, oh_alec["aa_change"].values, "mut")

    # ---- Anchor 2a: Nuvalent industry × alec proxy × ENU empirical
    oh_prox_alec = oh.copy()
    oh_prox_alec["sigma"] = oh_prox_alec["alec_score"]
    oh_prox_alec["S_mut"] = oh_prox_alec["forced_class"].map(ENU_EMP)
    all_rows += run_anchor("Nuvalent (ALK ENU × alec proxy)",
                           oh_prox_alec, oh_prox_alec["sigma"].values,
                           oh_prox_alec["S_mut"].values, NUVALENT_HITS,
                           oh_prox_alec["aa_change"].values, "ENU",
                           extra={"note": "ENU emp."})

    # ---- Anchor 2b (R10): Nuvalent industry × lor proxy × ENU empirical
    oh_prox_lor = oh.copy()
    oh_prox_lor["sigma"] = oh_prox_lor["lor_score"]
    oh_prox_lor["S_mut"] = oh_prox_lor["forced_class"].map(ENU_EMP)
    all_rows += run_anchor("Nuvalent (ALK ENU × lor proxy)",
                           oh_prox_lor, oh_prox_lor["sigma"].values,
                           oh_prox_lor["S_mut"].values, NUVALENT_HITS,
                           oh_prox_lor["aa_change"].values, "ENU",
                           extra={"note": "ENU emp."})

    # ---- Anchor 2c (R10): Nuvalent × zot proxy × ENU empirical
    oh_prox_zot = oh.copy()
    oh_prox_zot["sigma"] = oh_prox_zot["zot_score"]
    oh_prox_zot["S_mut"] = oh_prox_zot["forced_class"].map(ENU_EMP)
    all_rows += run_anchor("Nuvalent (ALK ENU × zot proxy)",
                           oh_prox_zot, oh_prox_zot["sigma"].values,
                           oh_prox_zot["S_mut"].values, NUVALENT_HITS,
                           oh_prox_zot["aa_change"].values, "ENU",
                           extra={"note": "ENU emp."})

    # ---- Anchor 3: Awad KRAS × sotorasib × MutII
    feng_sot = feng[feng["drug"] == "sotorasib"].copy()
    feng_sot["sigma"] = feng_sot["sotorasib_z"]
    feng_sot["S_mut"] = feng_sot["forced_class"].map(MUT_MUT)
    all_rows += run_anchor("Awad clinical (KRAS Mut × sot)",
                           feng_sot, feng_sot["sigma"].values,
                           feng_sot["S_mut"].values, AWAD_SOT_HITS,
                           feng_sot["aa_change"].values, "mut")

    # ---- Anchor 4: Awad KRAS × adagrasib × MutII
    feng_ada = feng[feng["drug"] == "adagrasib"].copy()
    feng_ada["sigma"] = feng_ada["adagrasib_z"]
    feng_ada["S_mut"] = feng_ada["forced_class"].map(MUT_MUT)
    all_rows += run_anchor("Awad clinical (KRAS Mut × ada)",
                           feng_ada, feng_ada["sigma"].values,
                           feng_ada["S_mut"].values, AWAD_ADA_HITS,
                           feng_ada["aa_change"].values, "mut")

    out = pd.DataFrame(all_rows)
    # Back up original
    out_path = ROOT / "data" / "derived" / "benchmark_r3.csv"
    bak = out_path.with_suffix(".csv.pre_r10")
    if out_path.exists() and not bak.exists():
        bak.write_bytes(out_path.read_bytes())
    out.to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({len(out)} rows)")

    # Also summarize the R10 proxy sensitivity
    sens = out[out["anchor"].str.contains("proxy")]
    print("\n=== R10 proxy sensitivity (Industry ENU, three proxies) ===")
    for _, r in sens.iterrows():
        if "σ alone" in r["ranking"]:
            print(f"  {r['anchor']:50s} | σ-alone r@100 = {r['recall@100']:.2f}  AUROC = {r['AUROC']:.3f}")
        elif "Framework" in r["ranking"]:
            print(f"  {r['anchor']:50s} | framework r@100 = {r['recall@100']:.2f}  AUROC = {r['AUROC']:.3f}  perm p = {r['perm_p_at_100']:.3f}")


if __name__ == "__main__":
    main()
