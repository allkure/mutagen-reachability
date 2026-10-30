"""Parse Feng 2022 PNAS Table S1 (Z-score matrices for sotorasib and adagrasib)
from the supplementary PDF into a long-format CSV.

Feng et al. 2022 PNAS 119(49):e2120512119 "Saturation mutagenesis of KRAS-G12C
identifies emerging resistance mechanisms to the clinical drugs sotorasib and
adagrasib." The screen used Twist Biosciences SATURATED mutagenesis — a
synthetic DNA library covering every single-amino-acid substitution of the
KRAS-G12C protein, built from oligos rather than a PCR mutagenesis chemistry.
This makes Feng's library CHEMISTRY-FREE: it's the KRAS analogue of Oh 2026's
prime-editing library for ALK.

Table S1 is 10 PDF pages (two 5-page blocks: sotorasib then adagrasib):
  - Rows: 188 positions (2..189 of KRAS-4A)
  - Cols: 20 amino-acid substitutions (A, C, D, ..., W, Y — no stops)
  - Values: average Z-score of 3 replicates for the given mutation after
    screening against each drug. Blank = WT or not scored.

pdfplumber extracts the tables cleanly (header + 20 numeric cells per row).

Output: data/raw/feng_2022/table_s1_zscores.csv
  position, wt_aa, mut_aa, sotorasib_z, adagrasib_z
"""

from __future__ import annotations
import csv
from pathlib import Path

import pdfplumber

AA_COLS = list("ACDEFGHIKLMNPQRSTVWY")  # 20 letters

# KRAS-4A reference (UniProt P01116-2, 189 aa). The Feng library has G12C
# baked in; we store that as position 12 = C. Reference verified by matching
# every unambiguous blank in Feng Table S1 (185/188 positions matched directly).
KRAS_4A_G12C_WT = (
    "MTEYKLVVVGACGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEE"
    "YSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVD"
    "TKQAQDLARSYGIPFIETSAKTRQRVEDAFYTLVREIRQYRLKKISKEEKTPGCVKIKKCIIM"
)
assert len(KRAS_4A_G12C_WT) == 189, f"KRAS-4A should be 189 aa, got {len(KRAS_4A_G12C_WT)}"
assert KRAS_4A_G12C_WT[11] == "C", "Position 12 should be C (G12C baked in)"
# True KRAS-4A WT has G at position 12.
KRAS_4A_TRUE_WT = KRAS_4A_G12C_WT[:11] + "G" + KRAS_4A_G12C_WT[12:]

PDF_PATH = Path(
    "/root/.claude/uploads/fc032df5-1689-5549-9823-8f122cba433e/"
    "1e0cbe44-pnas.2120512119.sapp.pdf"
)
# 0-indexed page ranges
SOT_PAGES = range(2, 7)   # pages 3-7
ADA_PAGES = range(7, 12)  # pages 8-12


def extract_block(pdf, page_range) -> dict[tuple[int, str], float]:
    """Return {(pos, aa): z} merging rows from all pages in page_range."""
    out: dict[tuple[int, str], float] = {}
    for pg_idx in page_range:
        page = pdf.pages[pg_idx]
        tables = page.extract_tables()
        for t in tables:
            for row in t:
                # Skip header row (first cell = "position")
                if not row or row[0] == "position":
                    continue
                # Row should be [pos_str, 20 value strs]
                if len(row) < 21:
                    continue
                try:
                    pos = int(row[0])
                except (ValueError, TypeError):
                    continue
                for aa, cell in zip(AA_COLS, row[1:21]):
                    if cell is None or cell.strip() == "":
                        continue  # blank
                    try:
                        out[(pos, aa)] = float(cell)
                    except ValueError:
                        pass
    return out


def validate_wt(
    sot: dict[tuple[int, str], float],
    ada: dict[tuple[int, str], float],
) -> tuple[dict[int, str], list]:
    """At each position, use KRAS_4A_G12C_WT as authoritative WT, but validate
    against the blank patterns in both blocks.

    Reports: (wt_by_pos, validation_issues).
    """
    wt: dict[int, str] = {}
    issues = []
    for pos in range(2, 190):
        ref_wt = KRAS_4A_G12C_WT[pos - 1]
        wt[pos] = ref_wt
        # Check reference WT is blank in both blocks (allow extra blanks for
        # excluded variants)
        sot_has_wt = (pos, ref_wt) in sot
        ada_has_wt = (pos, ref_wt) in ada
        if sot_has_wt or ada_has_wt:
            issues.append({
                "pos": pos, "ref_wt": ref_wt,
                "sot_has_wt": sot_has_wt, "sot_z": sot.get((pos, ref_wt)),
                "ada_has_wt": ada_has_wt, "ada_z": ada.get((pos, ref_wt)),
            })
    return wt, issues


def main():
    with pdfplumber.open(PDF_PATH) as pdf:
        sot = extract_block(pdf, SOT_PAGES)
        ada = extract_block(pdf, ADA_PAGES)
    print(f"Sotorasib Z-scores parsed: {len(sot)}  (expect ~188 pos × 19 cols = 3572 minus blanks)")
    print(f"Adagrasib Z-scores parsed: {len(ada)}")

    # Sanity: positions covered
    sot_positions = {p for (p, _) in sot}
    ada_positions = {p for (p, _) in ada}
    missing_sot = set(range(2, 190)) - sot_positions
    missing_ada = set(range(2, 190)) - ada_positions
    print(f"Sotorasib: positions covered = {len(sot_positions)}  missing = {sorted(missing_sot) or 'none'}")
    print(f"Adagrasib: positions covered = {len(ada_positions)}  missing = {sorted(missing_ada) or 'none'}")

    wt_by_pos, issues = validate_wt(sot, ada)
    print(f"\nUsing KRAS-4A WT (G12C baked in at position 12).")
    print(f"  N-term:  {KRAS_4A_G12C_WT[:30]}")
    print(f"  G12C region: ...{KRAS_4A_G12C_WT[8:15]}... (position 12 = C)")
    print(f"  C-term:  ...{KRAS_4A_G12C_WT[-20:]}")
    if issues:
        print(f"\nWT-validation issues ({len(issues)}): positions where WT itself got a Z-score")
        for iss in issues:
            print(f"  pos {iss['pos']}: ref WT={iss['ref_wt']} "
                  f"sot_z={iss['sot_z']} ada_z={iss['ada_z']}")
    else:
        print("\nAll 188 scored positions: reference WT is blank in both blocks ✓")

    # Write long-format CSV
    out_path = Path("data/raw/feng_2022/table_s1_zscores.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["position", "wt_aa", "mut_aa", "sotorasib_z", "adagrasib_z"])
        all_keys = sorted(set(sot.keys()) | set(ada.keys()))
        for pos, mut in all_keys:
            wt = wt_by_pos.get(pos, "?")
            if wt.startswith("?"):
                wt = ""
            if wt == mut:
                continue  # synonymous
            s = sot.get((pos, mut))
            a = ada.get((pos, mut))
            w.writerow([pos, wt, mut,
                        s if s is not None else "",
                        a if a is not None else ""])
    print(f"\nWrote {out_path}")

    # Write WT reference alongside
    with Path("data/raw/feng_2022/kras_4a_wt.txt").open("w") as f:
        f.write("# KRAS-4A reference (UniProt P01116-2, 189 aa, G12C baked in at position 12)\n")
        f.write(KRAS_4A_G12C_WT + "\n")

    # Spot checks on Awad 2021 clinical resistance variants
    print("\n--- Spot checks (Awad 2021 NEJM clinical resistance variants) ---")
    print(f"  Format: WT<pos>MUT  sotorasib_Z  adagrasib_Z")
    checks = [
        (8,  "L"), (8,  "M"), (8,  "A"),
        (12, "F"), (12, "V"), (12, "D"), (12, "G"),  # pos 12 is C (G12C baked in)
        (17, "E"),  # S17E — sensitising / nucleotide exchange
        (68, "L"), (68, "S"), (68, "M"),
        (95, "D"), (95, "Q"), (95, "R"), (95, "L"), (95, "Y"),
        (96, "C"), (96, "D"), (96, "H"), (96, "A"), (96, "N"),
    ]
    for p, mut in checks:
        wt = wt_by_pos.get(p, "?")
        s = sot.get((p, mut))
        a = ada.get((p, mut))
        s_fmt = f"{s:+6.2f}" if s is not None else "    NA"
        a_fmt = f"{a:+6.2f}" if a is not None else "    NA"
        print(f"  {wt}{p}{mut}  {s_fmt}  {a_fmt}")


if __name__ == "__main__":
    main()
