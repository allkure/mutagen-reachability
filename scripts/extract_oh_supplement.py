#!/usr/bin/env python3
"""Extract the Oh 2026 supplementary workbook tables to long-format CSV.

Source: Genome Biology 27:72 Additional file 2 (13059_2026_3977_MOESM2_ESM.xlsx).

Tables extracted
----------------
S2  main atlas, single dose per drug, 3,208 variants      -> oh_s2_main.csv
S3  exon 23 sub-library, lorlatinib at three doses        -> oh_s3_lorlatinib_doses.csv
S5  exon 23 sub-library, three drugs at three doses each  -> oh_s5_three_drug_doses.csv
S9  second sub-library, three drugs at two doses each     -> oh_s9_two_doses.csv

Every sheet has a title in row 1 and the header in row 2, so the body starts
at row 3.  Scores are written as text in some cells and as floats in others;
they are coerced and the coercion failures are counted rather than dropped
silently.

Integrity checks run at extraction time (any failure aborts):
  - header matches the expected column list for that table
  - the Type vocabulary is one of the two spellings the workbook uses
  - every Missense AA_Change parses as <WT aa><position><mutant aa>
  - no score column is entirely empty
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "raw" / "oh_2026" / "13059_2026_3977_MOESM2_ESM.xlsx"
OUT = ROOT / "data" / "raw" / "oh_2026"

AA_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY](\d+)[ACDEFGHIKLMNPQRSTVWY*]$")
TYPES_OK = {"Missense", "Synonymous", "Syn", "Intron", "Nonsense", "Splice Acceptor/Donor"}

SPECS = {
    "Table S2": ("oh_s2_main.csv", [
        "ID", "AA_Change", "Type",
        "Alectinib_Classification", "Alectinib_Resistance_Score",
        "Lorlatinib_Classification", "Lorlatinib_Resistance_Score",
        "Zotizalkib_Classification", "Zotizalkib_Resistance_Score",
        "Fitness_Score"]),
    "Table S3": ("oh_s3_lorlatinib_doses.csv", [
        "ID", "AA_Change", "Type",
        "Lorlatinib_2uM_Resistance_Score", "Lorlatinib_5uM_Resistance_Score",
        "Lorlatinib_10uM_Resistance_Score"]),
    "Table S5": ("oh_s5_three_drug_doses.csv", [
        "ID", "AA_Change", "Type",
        "Alectinib_20uM_Resistance_Score", "Alectinib_40uM_Resistance_Score",
        "Alectinib_80uM_Resistance_Score",
        "Loraltinib_5uM_Resistance_Score", "Loraltinib_10uM_Resistance_Score",
        "Loraltinib_20uM_Resistance_Score",
        "Zotizalkib_0.3uM_Resistance_Score", "Zotizalkib_0.9uM_Resistance_Score",
        "Zotizalkib_2.7uM_Resistance_Score"]),
    "Table S9": ("oh_s9_two_doses.csv", [
        "ID", "AA_Change", "Type",
        "Alectinib_3uM_Resistance_Score", "Alectinib_6uM_Resistance_Score",
        "Loraltinib_20uM_Resistance_Score", "Loraltinib_40uM_Resistance_Score",
        "Zotizalkib_0.3uM_Resistance_Score", "Zotizalkib_1uM_Resistance_Score",
        "Fitness_score"]),
}


def main() -> int:
    if not SRC.exists():
        sys.exit(f"missing workbook: {SRC}")

    for sheet, (fname, expected) in SPECS.items():
        raw = pd.read_excel(SRC, sheet_name=sheet, header=1, dtype=object)
        cols = [c for c in raw.columns if not str(c).startswith("Unnamed")]
        got = [str(c) for c in cols]
        if got != expected:
            sys.exit(f"{sheet}: header mismatch\n  expected {expected}\n  got      {got}")

        df = raw[cols].copy()
        df.columns = expected
        df = df[df["ID"].notna()].reset_index(drop=True)

        bad_types = set(df["Type"].dropna().unique()) - TYPES_OK
        if bad_types:
            sys.exit(f"{sheet}: unexpected Type values {bad_types}")

        mis = df[df["Type"] == "Missense"]
        unparsed = [a for a in mis["AA_Change"] if not AA_RE.match(str(a))]
        if unparsed:
            sys.exit(f"{sheet}: {len(unparsed)} missense AA_Change unparsed, e.g. {unparsed[:5]}")

        score_cols = [c for c in expected if c.endswith(("Score", "score"))]
        for c in score_cols:
            before = df[c].notna().sum()
            df[c] = pd.to_numeric(df[c], errors="coerce")
            after = df[c].notna().sum()
            if after == 0:
                sys.exit(f"{sheet}: score column {c} is entirely empty after coercion")
            if after < before:
                print(f"  note {sheet}.{c}: {before - after} values would not coerce to float")

        df["residue"] = [int(m.group(1)) if (m := AA_RE.match(str(a))) else pd.NA
                         for a in df["AA_Change"]]
        df.to_csv(OUT / fname, index=False)
        print(f"{sheet:9s} -> {fname:32s} {len(df):5d} rows, "
              f"{int((df.Type == 'Missense').sum()):4d} missense, "
              f"residues {df.residue.min()}-{df.residue.max()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
