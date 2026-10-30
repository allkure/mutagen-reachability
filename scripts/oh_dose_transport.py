#!/usr/bin/env python3
"""Internal dose-transport analysis of the Oh 2026 atlas.

Referee point addressed (second half)
-------------------------------------
Section 2.6 tests sigma against external Ba/F3 potency.  This script asks the
complementary internal question: within the atlas's own assay, how much does
the sigma ranking move when the drug concentration changes?

Inputs
------
Table S5  exon 23 sub-library, three drugs at three doses each  (primary)
Table S3  exon 23 sub-library, lorlatinib at three doses        (replication)
Table S9  second sub-library, three drugs at two doses each     (replication)
Table S2  the main single-dose atlas                            (comparator)

Controls, because a bare dose-to-dose correlation is not interpretable
---------------------------------------------------------------------
A  Which dose does the main atlas resemble?  Establishes which column
   corresponds to the condition the rest of the paper uses as sigma.
B  Cross-drug correlation at matched versus unmatched dose rank.  If matched
   pairs agree better, part of the dose-to-dose agreement is shared batch
   structure rather than drug response.
C  Does dispersion grow with dose?  A compressed low-dose column would
   attenuate its own correlations for reasons that are not about ranking.
D  Is the sub-library range-restricted relative to the full atlas?  Range
   restriction attenuates correlation and would be a confound.

Outputs
-------
data/derived/oh_dose_transport.csv
data/derived/oh_dose_transport.json
"""
from __future__ import annotations

import itertools
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "oh_2026"
OUT_CSV = ROOT / "data" / "derived" / "oh_dose_transport.csv"
OUT_JSON = ROOT / "data" / "derived" / "oh_dose_transport.json"

DOSE_RE = re.compile(r"_([\d.]+)uM_")
TOP_K = 25
WINDOW = (1172, 1215)          # exon 23 sub-library residue span

# The workbook spells lorlatinib "Loraltinib" in the dose tables.
SERIES = {
    "S5": ("oh_s5_three_drug_doses.csv",
           {"alectinib": "Alectinib", "lorlatinib": "Loraltinib", "zotizalkib": "Zotizalkib"}),
    "S3": ("oh_s3_lorlatinib_doses.csv", {"lorlatinib": "Lorlatinib"}),
    "S9": ("oh_s9_two_doses.csv",
           {"alectinib": "Alectinib", "lorlatinib": "Loraltinib", "zotizalkib": "Zotizalkib"}),
}
ATLAS_COL = {"alectinib": "Alectinib_Resistance_Score",
             "lorlatinib": "Lorlatinib_Resistance_Score",
             "zotizalkib": "Zotizalkib_Resistance_Score"}


def dose_of(col: str) -> float:
    m = DOSE_RE.search(col)
    if not m:
        raise ValueError(f"no dose in column name: {col}")
    return float(m.group(1))


def load_missense(fname: str) -> pd.DataFrame:
    """Missense rows, one per amino-acid change (duplicate pegRNAs averaged)."""
    d = pd.read_csv(RAW / fname)
    d = d[d.Type == "Missense"]
    return d.groupby("AA_Change").mean(numeric_only=True).reset_index()


def dose_columns(d: pd.DataFrame, prefix: str) -> list[str]:
    cols = [c for c in d.columns if c.startswith(prefix + "_") and c.endswith("Resistance_Score")]
    return sorted(cols, key=dose_of)


def rho(a: pd.Series, b: pd.Series) -> float:
    return float(stats.spearmanr(a, b).statistic)


def main() -> None:
    res: dict = {"top_k": TOP_K, "window": list(WINDOW), "series": {}, "controls": {}}
    rows = []

    # ------------------------------------------------------------------
    # dose-to-dose agreement within each (table, drug) series
    # ------------------------------------------------------------------
    for tag, (fname, drugs) in SERIES.items():
        d = load_missense(fname)
        for drug, prefix in drugs.items():
            cols = dose_columns(d, prefix)
            if len(cols) < 2:
                continue
            sub = d[["AA_Change"] + cols].dropna()
            lo, hi = cols[0], cols[-1]
            top_lo = set(sub.nlargest(TOP_K, lo).AA_Change)
            top_hi = set(sub.nlargest(TOP_K, hi).AA_Change)
            adj = [rho(sub[cols[i]], sub[cols[i + 1]]) for i in range(len(cols) - 1)]
            entry = {
                "n_variants": int(len(sub)),
                "doses": [dose_of(c) for c in cols],
                "fold_span": dose_of(hi) / dose_of(lo),
                "rho_low_vs_high": rho(sub[lo], sub[hi]),
                "rho_adjacent_steps": adj,
                "topk_overlap": int(len(top_lo & top_hi)),
                "sd_by_dose": [float(sub[c].std()) for c in cols],
            }
            res["series"][f"{tag} {drug}"] = entry
            rows.append(dict(table=tag, drug=drug, n=entry["n_variants"],
                             fold_span=entry["fold_span"],
                             rho_low_high=entry["rho_low_vs_high"],
                             rho_adjacent_median=float(np.median(adj)),
                             topk_overlap=entry["topk_overlap"]))

    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)

    # ------------------------------------------------------------------
    # control A: which S5 dose does the main atlas resemble?
    # ------------------------------------------------------------------
    s2 = pd.read_csv(RAW / "oh_s2_main.csv")
    s2 = s2[s2.Type == "Missense"]
    s5 = load_missense(SERIES["S5"][0])
    j = s2.merge(s5, on="AA_Change", how="inner")
    ctrl_a = {}
    for drug, prefix in SERIES["S5"][1].items():
        per = {}
        for c in dose_columns(s5, prefix):
            pair = j[[ATLAS_COL[drug], c]].dropna()
            per[str(dose_of(c))] = rho(pair[ATLAS_COL[drug]], pair[c])
        ctrl_a[drug] = {"rho_by_dose": per,
                        "best_match_dose": max(per, key=per.get),
                        "n": int(len(j))}
    res["controls"]["A_atlas_vs_dose"] = ctrl_a

    # ------------------------------------------------------------------
    # control B: cross-drug agreement, matched vs unmatched dose rank
    # ------------------------------------------------------------------
    all_cols = [c for c in s5.columns if c.endswith("Resistance_Score")]
    rank = {}
    for drug, prefix in SERIES["S5"][1].items():
        for i, c in enumerate(dose_columns(s5, prefix)):
            rank[c] = i
    dd = s5[all_cols].dropna()
    matched, unmatched = [], []
    for c1, c2 in itertools.combinations(all_cols, 2):
        if c1.split("_")[0] == c2.split("_")[0]:
            continue
        (matched if rank[c1] == rank[c2] else unmatched).append(rho(dd[c1], dd[c2]))
    res["controls"]["B_cross_drug"] = {
        "matched_dose_rank": {"n_pairs": len(matched), "median_rho": float(np.median(matched))},
        "unmatched_dose_rank": {"n_pairs": len(unmatched), "median_rho": float(np.median(unmatched))},
    }

    # ------------------------------------------------------------------
    # control C: dispersion by dose (already stored per series; summarise)
    # ------------------------------------------------------------------
    res["controls"]["C_sd_monotone_increasing"] = {
        k: bool(np.all(np.diff(v["sd_by_dose"]) > 0)) for k, v in res["series"].items()}

    # ------------------------------------------------------------------
    # control D: is the sub-library range-restricted vs the full atlas?
    # ------------------------------------------------------------------
    s2r = s2.dropna(subset=["residue"]).copy()
    s2r["residue"] = s2r.residue.astype(int)
    win = s2r[(s2r.residue >= WINDOW[0]) & (s2r.residue <= WINDOW[1])]
    ctrl_d = {}
    for drug, col in ATLAS_COL.items():
        f, w = s2r[col].dropna(), win[col].dropna()
        ctrl_d[drug] = {"full_sd": float(f.std()), "window_sd": float(w.std()),
                        "sd_ratio_window_over_full": float(w.std() / f.std())}
    res["controls"]["D_range_restriction"] = {
        "n_full_missense": int(len(s2r)), "n_window_missense": int(len(win)),
        "by_drug": ctrl_d,
        "range_restricted": bool(all(v["sd_ratio_window_over_full"] < 1 for v in ctrl_d.values())),
    }

    # headline summary
    lows = [v["rho_low_vs_high"] for v in res["series"].values()]
    adjs = [r for v in res["series"].values() for r in v["rho_adjacent_steps"]]
    ovl = [v["topk_overlap"] for v in res["series"].values()]
    res["headline"] = {
        "n_series": len(res["series"]),
        "rho_low_vs_high_range": [float(min(lows)), float(max(lows))],
        "rho_low_vs_high_median": float(np.median(lows)),
        "rho_adjacent_median": float(np.median(adjs)),
        "topk_overlap_range": [int(min(ovl)), int(max(ovl))],
        "cross_drug_matched_median_rho":
            res["controls"]["B_cross_drug"]["matched_dose_rank"]["median_rho"],
    }

    OUT_JSON.write_text(json.dumps(res, indent=2))
    print(json.dumps(res["headline"], indent=2))
    print("\ncontrol A (atlas resembles which dose):")
    for k, v in ctrl_a.items():
        print(f"  {k:11s} best match {v['best_match_dose']:>5s}  rho by dose "
              + ", ".join(f"{d}={r:+.3f}" for d, r in v["rho_by_dose"].items()))
    print("\ncontrol B:", json.dumps(res["controls"]["B_cross_drug"]))
    print("control C (sd rises with dose):", res["controls"]["C_sd_monotone_increasing"])
    print("control D range-restricted:", res["controls"]["D_range_restriction"]["range_restricted"],
          {k: round(v["sd_ratio_window_over_full"], 2) for k, v in ctrl_d.items()})


if __name__ == "__main__":
    main()
