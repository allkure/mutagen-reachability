#!/usr/bin/env python3
"""Post-selection null for the industry ENU C>G class-excess test.

Referee point addressed
-----------------------
The thirteen variants reported by the industry ENU screen
(Tangpeerachaikul 2025) are *post-selection*: they appear in the hit list
because they survived drug selection.  Testing their class distribution
against the raw pre-selection germline ENU spectrum attributes to the
supply term S what the paper's own model assigns to the selection term
sigma.  A C>G excess among hits is exactly what the model predicts if
C>G-reachable variants at those codons carry high sigma.

This script builds the correct nulls and re-runs the test.

Nulls
-----
N0  raw spectrum            pi(c) = S_ENU(c)                      [the old, confounded test]
N1  atlas-composition       pi(c) propto S_ENU(c) * n_c           [supply x how many variants exist in c]
N2  sigma-ranked top-N      pi(c) = fraction of the top-N sigma-ranked
                                    class-forced variants in class c
                                    (the referee's suggested null; supply-free)
N3  model-implied           pi(c) propto sum_{v in c} S_ENU(c) * exp(beta * sigma(v))
                                    (the paper's own model used as its own null)

N3 is the principled one: under H0 "the germline ENU spectrum describes the
in vitro library", N3 *is* the hit-class distribution the framework predicts.

Inversion
---------
The published Jeffreys interval on the observed hit proportion (3/13) is
mapped back through the N3 model to an implied in vitro C>G *supply* rate,
by solving pi_3(C>G; s) = p for the spectrum in which C>G carries mass s and
the remaining five classes keep their germline ratios scaled to 1 - s.
This replaces the previous inversion, which read the hit proportion directly
as a supply rate and therefore absorbed the whole sigma tilt into S.

Outputs
-------
data/derived/enu_class_test.csv      one row per (null, beta, drug-arm)
data/derived/enu_class_test.json     headline numbers for NUMBERS.md
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.spectrum import load_spectra  # noqa: E402

DECOMP = ROOT / "data" / "derived" / "oh_decomposition_all_drugs.csv"
OUT_CSV = ROOT / "data" / "derived" / "enu_class_test.csv"
OUT_JSON = ROOT / "data" / "derived" / "enu_class_test.json"

CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

# The industry ENU screen's aggregate single-amino-acid hit list
# (Tangpeerachaikul 2025).  Twelve of thirteen are class-forced; R1192P lies
# outside the atlas residue window and is excluded from the class test.
N_HITS = 13
N_CG_HITS = 3                      # V1180L, T1151R, G1269A
CG_HITS = ["V1180L", "T1151R", "G1269A"]

# sigma column per drug arm.  The screen runs alectinib, lorlatinib and
# neladalkib; the atlas has no neladalkib, so zotizalkib stands in as the
# fourth-generation proxy, matching the convention used elsewhere in the paper.
ARM_SIGMA = {"alectinib": "alec_score", "lorlatinib": "lor_score",
             "neladalkib (zot proxy)": "zot_score"}

BETAS = [0.25, 0.5, 1.0, 2.0, 4.0]
BETA_PRIMARY = 1.0
TOP_N = [50, 100, 200]
Z = 1.959963984540054


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def rescaled_spectrum(base: dict[str, float], s_cg: float) -> dict[str, float]:
    """Germline spectrum with C>G reset to s_cg, other five rescaled to 1-s_cg."""
    rest = {k: v for k, v in base.items() if k != "C>G"}
    scale = (1.0 - s_cg) / sum(rest.values())
    out = {k: v * scale for k, v in rest.items()}
    out["C>G"] = s_cg
    return out


def pi_model(var: pd.DataFrame, spec: dict[str, float], sigma_col: str,
             beta: float) -> dict[str, float]:
    """N3: class probabilities implied by P(v) propto S(c(v)) * exp(beta*sigma(v))."""
    w = var[sigma_col].to_numpy(dtype=float)
    w = np.exp(beta * (w - w.max()))                 # shift for numerical stability
    s = var["forced_class"].map(spec).to_numpy(dtype=float)
    p = s * w
    tot = p.sum()
    return {c: float(p[(var["forced_class"] == c).to_numpy()].sum() / tot) for c in CLASSES}


def pi_composition(var: pd.DataFrame, spec: dict[str, float]) -> dict[str, float]:
    """N1: supply weighted by how many class-forced variants exist in each class."""
    n = var["forced_class"].value_counts()
    p = {c: spec[c] * int(n.get(c, 0)) for c in CLASSES}
    tot = sum(p.values())
    return {c: p[c] / tot for c in CLASSES}


def pi_sigma_topn(var: pd.DataFrame, sigma_col: str, n_top: int) -> dict[str, float]:
    """N2: class composition of the top-n sigma-ranked class-forced variants."""
    top = var.nlargest(n_top, sigma_col)
    vc = top["forced_class"].value_counts()
    return {c: float(vc.get(c, 0)) / n_top for c in CLASSES}


def binom_tail(k: int, n: int, p: float) -> float:
    """P(X >= k) under Binomial(n, p)."""
    if p <= 0.0:
        return 0.0 if k > 0 else 1.0
    return float(stats.binom.sf(k - 1, n, p))


def jeffreys(k: int, n: int) -> tuple[float, float, float]:
    a, b = k + 0.5, n - k + 0.5
    return (float(stats.beta.ppf(0.025, a, b)),
            float(stats.beta.mean(a, b)),
            float(stats.beta.ppf(0.975, a, b)))


# --------------------------------------------------------------------------
def main() -> None:
    spectra = load_spectra(ROOT / "data" / "reference" / "spectra.yaml")
    enu = dict(spectra["ENU_Arnold2012_nominal"].values)

    df = pd.read_csv(DECOMP)
    # One row per class-forced atlas variant; sigma columns are mutagen-invariant.
    var = (df[df.mutagen == "MutazymeII"]
           .drop_duplicates(subset="aa_change")
           [["aa_change", "forced_class", "alec_score", "lor_score", "zot_score"]]
           .reset_index(drop=True))
    # Universe size is whatever the decomposition holds; it changed from 1,259 to
    # 1,427 when the ALK codon table was rebuilt from the Oh Table S11 construct.
    assert 1000 < len(var) < 2200, len(var)
    assert var.aa_change.is_unique and var.forced_class.notna().all()

    # Pooled sigma: a variant enters the aggregate hit list if it wins in any
    # arm, so the pooled selection weight is the best score across arms.
    var["pooled_score"] = var[["alec_score", "lor_score", "zot_score"]].max(axis=1)

    rows = []

    # --- N0: raw spectrum (the confounded test, retained for comparison) ----
    rows.append(dict(null="N0 raw germline spectrum", arm="--", beta=np.nan,
                     n_top=np.nan, pi_cg=enu["C>G"],
                     p_ge3=binom_tail(N_CG_HITS, N_HITS, enu["C>G"])))

    # --- N1: supply x atlas composition ------------------------------------
    pi1 = pi_composition(var, enu)
    rows.append(dict(null="N1 supply x atlas composition", arm="--", beta=np.nan,
                     n_top=np.nan, pi_cg=pi1["C>G"],
                     p_ge3=binom_tail(N_CG_HITS, N_HITS, pi1["C>G"])))

    # --- N2: sigma-ranked top-N (referee's null, supply-free) --------------
    for arm, col in list(ARM_SIGMA.items()) + [("pooled", "pooled_score")]:
        for n_top in TOP_N:
            pi2 = pi_sigma_topn(var, col, n_top)
            rows.append(dict(null="N2 sigma-ranked top-N", arm=arm, beta=np.nan,
                             n_top=n_top, pi_cg=pi2["C>G"],
                             p_ge3=binom_tail(N_CG_HITS, N_HITS, pi2["C>G"])))

    # --- N4: atlas composition only (no supply, no selection) --------------
    pi4_cg = float((var.forced_class == "C>G").sum()) / len(var)
    rows.append(dict(null="N4 atlas composition only", arm="--", beta=np.nan,
                     n_top=np.nan, pi_cg=pi4_cg,
                     p_ge3=binom_tail(N_CG_HITS, N_HITS, pi4_cg)))

    # --- N3: model-implied --------------------------------------------------
    for arm, col in list(ARM_SIGMA.items()) + [("pooled", "pooled_score")]:
        for beta in BETAS:
            pi3 = pi_model(var, enu, col, beta)
            rows.append(dict(null="N3 model-implied", arm=arm, beta=beta,
                             n_top=np.nan, pi_cg=pi3["C>G"],
                             p_ge3=binom_tail(N_CG_HITS, N_HITS, pi3["C>G"])))

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)

    # --- inversion: implied in vitro C>G SUPPLY rate ------------------------
    # Map the Jeffreys interval on the hit proportion back through N3.
    lo_p, mean_p, hi_p = jeffreys(N_CG_HITS, N_HITS)

    def solve_s(target_p: float, col: str, beta: float) -> float:
        f = lambda s: pi_model(var, rescaled_spectrum(enu, s), col, beta)["C>G"] - target_p
        if f(1e-9) > 0:
            return 0.0
        if f(1 - 1e-9) < 0:
            return 1.0
        return float(optimize.brentq(f, 1e-9, 1 - 1e-9, xtol=1e-10))

    inv = {}
    for label, col in (("pooled", "pooled_score"), ("alectinib", "alec_score")):
        inv[label] = {
            "s_lo": solve_s(lo_p, col, BETA_PRIMARY),
            "s_point": solve_s(N_CG_HITS / N_HITS, col, BETA_PRIMARY),
            "s_hi": solve_s(hi_p, col, BETA_PRIMARY),
        }

    pooled_pi3 = pi_model(var, enu, "pooled_score", BETA_PRIMARY)
    headline = {
        "n_hits": N_HITS,
        "n_cg_hits": N_CG_HITS,
        "cg_hits": CG_HITS,
        "n_classforced_atlas": int(len(var)),
        "n_cg_classforced_atlas": int((var.forced_class == "C>G").sum()),
        "germline_cg_rate": enu["C>G"],
        "jeffreys_hit_proportion": {"lo": lo_p, "mean": mean_p, "hi": hi_p},
        "pi_cg": {
            "N0_raw": enu["C>G"],
            "N1_composition": pi1["C>G"],
            "N2_top100_pooled": pi_sigma_topn(var, "pooled_score", 100)["C>G"],
            "N3_pooled_beta1": pooled_pi3["C>G"],
        },
        "p_ge3": {
            "N0_raw": binom_tail(N_CG_HITS, N_HITS, enu["C>G"]),
            "N1_composition": binom_tail(N_CG_HITS, N_HITS, pi1["C>G"]),
            "N2_top100_pooled": binom_tail(
                N_CG_HITS, N_HITS, pi_sigma_topn(var, "pooled_score", 100)["C>G"]),
            "N3_pooled_beta1": binom_tail(N_CG_HITS, N_HITS, pooled_pi3["C>G"]),
        },
        "p_ge3_N3_beta_range": {
            str(b): binom_tail(N_CG_HITS, N_HITS,
                               pi_model(var, enu, "pooled_score", b)["C>G"])
            for b in BETAS
        },
        "implied_cg_supply_rate": inv,
        "n_cg_in_top100_pooled": int(
            (var.nlargest(100, "pooled_score").forced_class == "C>G").sum()),
        "pi_cg_N4_composition": pi4_cg,
        "p_ge3_N4_composition": binom_tail(N_CG_HITS, N_HITS, pi4_cg),
        # Is C>G actually supplied at practical library depth under the
        # germline rate?  P_detect = 1 - exp(-N*mu*S*r) with r = 1 route.
        "p_detect_cg_germline": {
            "N": 1e7, "mu": 1e-4, "S": enu["C>G"], "routes": 1,
            "expected_carriers": 1e7 * 1e-4 * enu["C>G"],
            "p_detect": float(1.0 - np.exp(-1e7 * 1e-4 * enu["C>G"])),
        },
    }
    OUT_JSON.write_text(json.dumps(headline, indent=2))

    pd.set_option("display.width", 160)
    print(out.to_string(index=False, float_format=lambda v: f"{v:.6g}"))
    print("\nHEADLINE")
    print(json.dumps(headline, indent=2))


if __name__ == "__main__":
    main()
