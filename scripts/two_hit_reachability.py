#!/usr/bin/env python3
"""Quantify how strongly two-nucleotide variants are suppressed, rather than
asserting that they are impossible.

Referee point addressed
-----------------------
The draft claims R68L and V8M are "unreachable to any single-nucleotide
mutagen at any library depth".  That is overstated.  A mutagen applies a
per-base rate to every base, so two substitutions in one codon are rare but
not forbidden, and the paper's own compound-mutation anchor (Takei 2025) is
an error-prone PCR screen that exists *because* multi-nucleotide events
happen.  This script replaces the impossibility claim with the arithmetic.

Model
-----
A specific substitution at a specific base occurs at rate mu * S(class)
per base per round, where mu is the per-base rate and S(class) the mutagen's
six-class pyrimidine weight for that substitution.

One-nucleotide variant, route with class c:
    lambda_1 = N * mu * S(c)                                   (expected carriers)

Two-nucleotide variant, route requiring substitutions of classes c1 and c2:
    simultaneous, one round   lambda_2 = N * mu^2 * S(c1) * S(c2)
    sequential over R rounds  an intermediate single mutant arises and then
                              acquires the second hit, which over R rounds of
                              mutagenesis-plus-expansion contributes a further
                              factor of order R(R+1)/2.
Both orderings of the two substitutions are counted.

Route multiplicity is enumerated from the genetic code rather than assumed.

Outputs
-------
data/derived/two_hit_reachability.csv
data/derived/two_hit_reachability.json
"""
from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.codon_tools import sub_class, translate  # noqa: E402
from src.spectrum import load_spectra  # noqa: E402

OUT_CSV = ROOT / "data" / "derived" / "two_hit_reachability.csv"
OUT_JSON = ROOT / "data" / "derived" / "two_hit_reachability.json"

BASES = "ACGT"

# The two structurally blind-spot variants in the composite KRAS G12C panel.
# Both are reported in the accompanying in vitro deep mutational scan of
# Awad 2021, not in that paper's clinical cohort.
TARGETS = [
    ("R68L", "AGG", "L"),
    ("V8M", "GTA", "M"),
]

N_LIB = 1e7          # Ba/F3 library size used throughout the paper
MU = 1e-4            # per base per round (Methods, library-depth subsection)
ROUNDS = [1, 5, 20]  # mutagenesis rounds considered
MU_SENS = [1e-5, 1e-4, 1e-3]


def two_nt_routes(wt: str, target_aa: str) -> list[tuple[str, tuple[str, str]]]:
    """Every target codon exactly two substitutions away, with the two classes."""
    out = []
    for codon in ("".join(c) for c in itertools.product(BASES, repeat=3)):
        if translate(codon) != target_aa:
            continue
        diff = [i for i in range(3) if codon[i] != wt[i]]
        if len(diff) != 2:
            continue
        classes = tuple(sub_class(wt[i], codon[i]) for i in diff)
        out.append((codon, classes))
    return out


def one_nt_routes(wt: str, target_aa: str) -> list[tuple[str, str]]:
    out = []
    for codon in ("".join(c) for c in itertools.product(BASES, repeat=3)):
        if translate(codon) != target_aa:
            continue
        diff = [i for i in range(3) if codon[i] != wt[i]]
        if len(diff) == 1:
            out.append((codon, sub_class(wt[diff[0]], codon[diff[0]])))
    return out


def main() -> None:
    spectra = load_spectra(ROOT / "data" / "reference" / "spectra.yaml")
    mutagens = {
        "ENU (germline ref)": dict(spectra["ENU_Arnold2012_nominal"].values),
        "Mutazyme II": dict(spectra["MutazymeII_DoiSupp2022"].values),
    }

    rows = []
    detail = {}
    for variant, wt, aa in TARGETS:
        r2 = two_nt_routes(wt, aa)
        r1 = one_nt_routes(wt, aa)
        detail[variant] = {
            "wt_codon": wt, "target_aa": aa,
            "one_nt_routes": r1,
            "two_nt_routes": [{"codon": c, "classes": list(cl)} for c, cl in r2],
            "n_two_nt_routes": len(r2),
        }
        for mname, spec in mutagens.items():
            for mu in MU_SENS:
                # simultaneous, single round, summed over routes and both orderings
                lam_simul = sum(
                    N_LIB * mu * spec[c1] * mu * spec[c2] * 2.0 for _, (c1, c2) in r2)
                for R in ROUNDS:
                    lam_seq = lam_simul * (R * (R + 1) / 2.0)
                    rows.append(dict(
                        variant=variant, mutagen=mname, mu=mu, rounds=R,
                        n_routes=len(r2),
                        expected_carriers=lam_seq,
                        p_at_least_one=1.0 - np.exp(-lam_seq)))

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)

    # Headline comparison at the paper's own parameters: a typical one-nucleotide
    # variant against the two-nucleotide variants, same library, same mutagen.
    enu = mutagens["ENU (germline ref)"]
    mut = mutagens["Mutazyme II"]

    def lam1(spec: dict[str, float], cls: str, mu: float = MU) -> float:
        return N_LIB * mu * spec[cls]

    headline = {
        "N": N_LIB, "mu_primary": MU, "rounds_primary": 20,
        "routes": {v: detail[v]["n_two_nt_routes"] for v, _, _ in TARGETS},
        "two_nt_detail": detail,
    }
    for mname, spec in mutagens.items():
        key = mname.split(" (")[0].replace(" ", "")
        sel = out[(out.mutagen == mname) & (out.mu == MU) & (out.rounds == 20)]
        per_variant = {r.variant: {"expected_carriers": r.expected_carriers,
                                   "p_at_least_one": r.p_at_least_one}
                       for r in sel.itertuples()}
        # reference one-nucleotide variant: median-supply class for this mutagen
        med_cls = sorted(spec, key=lambda c: spec[c])[len(spec) // 2]
        ref = lam1(spec, med_cls)
        per_variant["_reference_one_nt"] = {
            "class": med_cls, "S": spec[med_cls], "expected_carriers": ref}
        for v, _, _ in TARGETS:
            per_variant[v]["fold_below_one_nt"] = ref / per_variant[v]["expected_carriers"]
        headline[key] = per_variant

    OUT_JSON.write_text(json.dumps(headline, indent=2, default=float))

    pd.set_option("display.width", 160)
    print(out.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
    print("\nHEADLINE")
    print(json.dumps(headline, indent=2, default=float))


if __name__ == "__main__":
    main()
