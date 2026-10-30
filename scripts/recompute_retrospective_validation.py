#!/usr/bin/env python3
"""Regenerate `data/derived/retrospective_validation.csv` under σ = R.

Rebuilds the per-variant chemistry-free and simulated ranks for the three
retrospective validation anchors (Doi 2022, Takei 2025, Nuvalent 2025)
using the up-to-date σ=Resistance(v) convention and the current state of
`data/derived/oh_decomposition_all_drugs.csv`. Writes a schema-compatible
replacement so Figure 4 renders without changes.

Each anchor pairs its reported variants with a (mutagen, drug) chemistry
pair:
  - Doi 2022     — Mutazyme II × alectinib
  - Takei 2025   — Mutazyme II × zotizalkib (as proxy for neladalkib)
  - Nuvalent 2025 — ENU (empirical) × alectinib (as proxy for neladalkib)

The Nuvalent rows use the ENU-empirical class rates (ENU with the
inverted in vitro C>G rate), matching the paper's primary industry-anchor
column; the mutagen_rate column still reports that empirical rate so the
figure can colour by class and surface the chemistry rate inline.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.spectrum import load_spectra


DOI_HITS = [
    ("V1180L", "GTT", "C>G"),
    ("L1196M", "CTG", "C>A"),
    ("I1171N", "ATC", "T>A"),
    ("I1171T", "ATC", "T>C"),
    ("L1198F", "CTC", "C>T"),
    ("G1202R", "GGA", None),  # multi-class
]

TAKEI_HITS = [
    ("L1256F", "CTC", "C>T"),
    ("E1129V", "GAG", "T>A"),
    ("E1129K", "GAG", "C>T"),
    ("C1156Y", "TGC", "C>T"),
    ("T1151R", "ACG", "C>G"),
]

# Nuvalent industry ENU hits (aggregate across arms).  R1192P is outside the
# atlas residue window and remains unanchored; G1202R is multi-class.
# ENU-empirical class rates for the C>G class use the median of the Jeffreys
# posterior from §2.3 (23 %); the other classes keep the Arnold 2012 nominal
# rates, renormalised.
NUVALENT_HITS = [
    ("I1171N", "ATC", "T>A"),
    ("I1171S", "ATC", "T>G"),
    ("I1171T", "ATC", "T>C"),
    ("F1174C", "TTC", "T>G"),  # Oh atlas aa-level; multi-codon handled below
    ("F1174V", "TTC", "T>G"),
    ("V1180L", "GTT", "C>G"),
    ("T1151K", "ACG", "C>A"),
    ("T1151R", "ACG", "C>G"),
    ("E1129K", "GAG", "C>T"),
    ("R1192P", "?",   None),  # outside atlas window
    ("L1196M", "CTG", "C>A"),
    ("G1269A", "GGA", "C>G"),
    ("G1202R", "GGA", None),  # multi-class
]

# Spectrum class rates for display
SPECTRA = load_spectra(ROOT / "data" / "reference" / "spectra.yaml")
MUT_MUT = dict(SPECTRA["MutazymeII_DoiSupp2022"].values)
MUT_ENU_NOMINAL = dict(SPECTRA["ENU_Arnold2012_nominal"].values)

# ENU empirical: C>G at 0.23 (median Jeffreys), other classes proportional.
_rest = {k: v for k, v in MUT_ENU_NOMINAL.items() if k != "C>G"}
_scale = (1 - 0.23) / sum(_rest.values())
MUT_ENU_EMP = {
    k: (0.23 if k == "C>G" else v * _scale)
    for k, v in MUT_ENU_NOMINAL.items()
}


def load_oh_ranks() -> pd.DataFrame:
    """Load the (drug, mutagen, aa_change) → (rank_chemfree, rank_simulated) table
    produced by notebook/05_prospective_predictions.py under the current σ=R
    convention."""
    df = pd.read_csv(ROOT / "data" / "derived" / "oh_decomposition_all_drugs.csv")
    return df[[
        "drug", "mutagen", "aa_change",
        "rank_chemfree", "rank_simulated",
    ]]


def rows_for_anchor(
    hits, screen_label: str, drug: str, mutagen_label: str,
    spectrum: dict, oh_ranks: pd.DataFrame,
) -> list[dict]:
    """One row per hit; look up ranks from Oh decomposition under (drug, mutagen)."""
    sub = oh_ranks[(oh_ranks["drug"] == drug) & (oh_ranks["mutagen"] == mutagen_label)]
    out = []
    for variant, wt_codon, cls in hits:
        row = {
            "variant": variant, "wt_codon": wt_codon,
            "screen": screen_label,
        }
        if cls is None:
            row["reach"] = "multi-class" if wt_codon != "?" else ""
            row["class"] = ""
            row["mutagen_rate"] = ""
            row["rank_chemfree"] = ""
            row["rank_simulated"] = ""
        else:
            row["reach"] = "class-forced"
            row["class"] = cls
            row["mutagen_rate"] = round(spectrum[cls], 3)
            m = sub[sub["aa_change"] == variant]
            if not m.empty:
                r = m.iloc[0]
                row["rank_chemfree"]  = int(r["rank_chemfree"])
                row["rank_simulated"] = int(r["rank_simulated"])
            else:
                row["rank_chemfree"] = ""
                row["rank_simulated"] = ""
        out.append(row)
    return out


def main() -> int:
    oh = load_oh_ranks()
    rows: list[dict] = []
    rows += rows_for_anchor(DOI_HITS,      "Doi 2022",      "alectinib",  "MutazymeII", MUT_MUT, oh)
    rows += rows_for_anchor(TAKEI_HITS,    "Takei 2025",    "zotizalkib", "MutazymeII", MUT_MUT, oh)
    # Nuvalent anchor: display the NOMINAL ENU rate so the figure's rate column
    # and the ranks are the same spectrum (rank_sim in oh_decomp uses nominal).
    rows += rows_for_anchor(NUVALENT_HITS, "Nuvalent 2025", "alectinib",  "ENU",        MUT_ENU_NOMINAL, oh)
    out = pd.DataFrame(rows, columns=[
        "variant", "wt_codon", "reach", "class", "mutagen_rate",
        "rank_chemfree", "rank_simulated", "screen",
    ])
    out_path = ROOT / "data" / "derived" / "retrospective_validation.csv"
    bak = out_path.with_suffix(".csv.pre_sigma_fix")
    if out_path.exists() and not bak.exists():
        bak.write_bytes(out_path.read_bytes())
    out.to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({len(out)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
