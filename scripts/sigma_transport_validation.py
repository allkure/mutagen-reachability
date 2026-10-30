#!/usr/bin/env python3
"""External validation of the sigma-transport assumption.

Referee point addressed
-----------------------
The benchmark's load-bearing assumption is that sigma, measured in the Oh
atlas at mild selection in H3122 (20 nM alectinib, 5 nM lorlatinib, 300 nM
zotizalkib), orders variants the same way a Ba/F3 screen at 5-15x IC50 does.
Until now this was asserted and flagged as a limitation.

This script tests it directly against measured Ba/F3 cellular IC50 values
compiled from the published ALK TKI literature (the author's harmonised
potency dataset).  That is a stronger test than the atlas's own internal
dose series, because it uses a different cell background, a different assay
and different laboratories.

Method
------
Fold resistance is computed WITHIN each (study, drug) stratum against that
study's own wild-type EML4-ALK measurement, so inter-laboratory differences
in absolute IC50 cancel.  Variants measured in more than one study are
summarised by the median of log10 fold change, with the inter-study spread
reported.  Agreement with the atlas is Spearman rank correlation, which is
the relevant statistic because every downstream use of sigma in this paper
is a ranking.

A second question the same data answers: the atlas has no neladalkib, so
the paper proxies it with alectinib, lorlatinib or zotizalkib.  Measured
neladalkib potency says which proxy is defensible.

Outputs
-------
data/derived/sigma_transport.csv    per-variant merged table
data/derived/sigma_transport.json   correlations and proxy ranking
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
POTENCY = ROOT / "data" / "raw" / "potency" / "alk_tki_potency_dataset_v0.csv"
DECOMP = ROOT / "data" / "derived" / "oh_decomposition_all_drugs.csv"
OUT_CSV = ROOT / "data" / "derived" / "sigma_transport.csv"
OUT_JSON = ROOT / "data" / "derived" / "sigma_transport.json"

# atlas sigma column for each drug the atlas actually measured
ATLAS_SIGMA = {"alectinib": "alec_score",
               "lorlatinib": "lor_score",
               "zotizalkib": "zot_score"}
PROXY_TARGET = "neladalkib"          # drug the atlas does not contain
N_BOOT = 10000
SEED = 20261016


def load_potency() -> pd.DataFrame:
    d = pd.read_csv(POTENCY)
    d = d[(d.gene == "ALK")
          & (d.value_type == "Cellular_IC50")
          & (d.cell_line == "Ba/F3")
          & (d.expression_construct == "EML4-ALK_v1")
          & d.value_nM.notna()].copy()
    # single amino-acid substitutions only: drop compounds and indels
    pv = d.protein_variant.fillna("WT")
    single = (~pv.str.contains(r"\+", na=False)) & (~pv.str.contains("del|ins|dup", case=False, na=False))
    d = d[single].copy()
    d["protein_variant"] = pv[single]
    # upper bounds ("> x nM") censor the value; keep but flag.
    # The column is the strings "yes"/"no", so astype(bool) would mark every
    # row True.  Map explicitly.
    d["censored"] = (d.value_is_upper_bound.astype(str).str.strip().str.lower()
                     .isin({"yes", "true", "1"}))
    return d


def fold_change_vs_wt(d: pd.DataFrame) -> pd.DataFrame:
    """log10 fold resistance within each (study, drug) stratum."""
    rows = []
    for (paper, drug), g in d.groupby(["paper_id", "drug_canonical"]):
        wt = g[g.protein_variant == "WT"]
        if wt.empty:
            continue
        wt_nM = float(wt.value_nM.median())
        if not np.isfinite(wt_nM) or wt_nM <= 0:
            continue
        for r in g[g.protein_variant != "WT"].itertuples():
            if r.value_nM <= 0:
                continue
            rows.append(dict(paper_id=paper, drug=drug, variant=r.protein_variant,
                             ic50_nM=float(r.value_nM), wt_nM=wt_nM,
                             log10_fc=float(np.log10(r.value_nM / wt_nM)),
                             censored=bool(r.censored)))
    return pd.DataFrame(rows)


def boot_spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    n = len(x)
    out = []
    for _ in range(N_BOOT):
        i = rng.integers(0, n, n)
        if len(np.unique(x[i])) < 3 or len(np.unique(y[i])) < 3:
            continue
        out.append(stats.spearmanr(x[i], y[i]).statistic)
    if not out:
        return (float("nan"), float("nan"))
    return (float(np.quantile(out, 0.025)), float(np.quantile(out, 0.975)))


def main() -> None:
    if not POTENCY.exists():
        sys.exit(f"missing input: {POTENCY}")

    pot = load_potency()
    fc = fold_change_vs_wt(pot)

    # median across studies, with spread
    agg = (fc.groupby(["drug", "variant"])
             .agg(log10_fc=("log10_fc", "median"),
                  n_studies=("paper_id", "nunique"),
                  n_meas=("log10_fc", "size"),
                  spread=("log10_fc", lambda s: float(s.max() - s.min())),
                  any_censored=("censored", "any"))
             .reset_index())

    oh = pd.read_csv(DECOMP)
    atlas = (oh[oh.drug == "alectinib"]
             .drop_duplicates(subset="aa_change")
             [["aa_change", "forced_class", "alec_score", "lor_score", "zot_score"]]
             .rename(columns={"aa_change": "variant"}))
    m = agg.merge(atlas, on="variant", how="inner")
    m.to_csv(OUT_CSV, index=False)

    result: dict = {"n_potency_rows": int(len(pot)),
                    "n_studies": int(pot.paper_id.nunique()),
                    "transport": {}, "proxy": {}}

    # ---- transport: atlas sigma vs measured Ba/F3 fold resistance ---------
    for drug, col in ATLAS_SIGMA.items():
        sub = m[(m.drug == drug) & m[col].notna()]
        if len(sub) < 5:
            result["transport"][drug] = {"n": int(len(sub)), "note": "too few overlapping variants"}
            continue
        r = stats.spearmanr(sub[col], sub.log10_fc)
        lo, hi = boot_spearman(sub[col].to_numpy(float), sub.log10_fc.to_numpy(float))
        result["transport"][drug] = {
            "n": int(len(sub)),
            "spearman_rho": float(r.statistic), "p": float(r.pvalue),
            "boot95": [lo, hi],
            "pearson_r": float(stats.pearsonr(sub[col], sub.log10_fc).statistic),
            "variants": sorted(sub.variant),
            "median_inter_study_spread_log10": float(sub.spread.median()),
        }

    # ---- proxy choice for the drug the atlas lacks ------------------------
    tgt = m[m.drug == PROXY_TARGET]
    for drug, col in ATLAS_SIGMA.items():
        sub = tgt[tgt[col].notna()]
        if len(sub) < 5:
            result["proxy"][drug] = {"n": int(len(sub)), "note": "too few overlapping variants"}
            continue
        r = stats.spearmanr(sub[col], sub.log10_fc)
        lo, hi = boot_spearman(sub[col].to_numpy(float), sub.log10_fc.to_numpy(float))
        result["proxy"][drug] = {"n": int(len(sub)), "spearman_rho": float(r.statistic),
                                 "p": float(r.pvalue), "boot95": [lo, hi]}
    ranked = [(k, v["spearman_rho"]) for k, v in result["proxy"].items() if "spearman_rho" in v]
    if ranked:
        result["best_proxy"] = max(ranked, key=lambda kv: kv[1])[0]

    # ---- selection-pressure context --------------------------------------
    wt_only = pot[pot.protein_variant == "WT"]
    result["wt_ic50_nM_by_drug"] = {
        d: float(g.value_nM.median()) for d, g in wt_only.groupby("drug_canonical")}
    result["atlas_dose_nM"] = {"alectinib": 20, "lorlatinib": 5, "zotizalkib": 300}
    result["atlas_dose_over_wt_ic50"] = {
        d: result["atlas_dose_nM"][d] / result["wt_ic50_nM_by_drug"][d]
        for d in ATLAS_SIGMA if d in result["wt_ic50_nM_by_drug"]}

    OUT_JSON.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
