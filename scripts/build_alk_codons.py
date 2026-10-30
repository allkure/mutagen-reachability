#!/usr/bin/env python3
"""Build the ALK kinase-domain codon table from the Oh 2026 atlas construct.

Source
------
Oh et al. 2026, Genome Biology 27:72, Additional file 2, Table S11, row
"EML4_ALK variant1_WT": the wild-type coding sequence of the construct the
prime-editing atlas was built in.  Using the construct's own sequence means
the codon table describes exactly the molecules the atlas assayed.

Why this replaces the previous reference
----------------------------------------
Earlier work used a 23-codon hand-curated slice plus a conservative fallback
that enumerated routes over every synonymous codon for the wild-type residue
and assigned a forced class only when all of them agreed.  That fallback is
sound but leaves many variants unassigned, and the curated slice contained
two transcription errors that this table exposes:

    residue 1181   recorded CGG, actually CGC   (same residue, wrong 3rd base)
    residue 1268   recorded GAC, actually ATT   (wrong residue entirely)

Verification
------------
The construct is searched in all three reading frames for the wild-type
peptide the atlas itself implies, and the resulting table must reproduce the
wild-type residue at every atlas position.  Nothing is written unless that
check is total.

Output
------
data/reference/codons_alk_full.yaml
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.codon_tools import translate  # noqa: E402

XLSX = ROOT / "data" / "raw" / "oh_2026" / "13059_2026_3977_MOESM2_ESM.xlsx"
ATLAS = ROOT / "data" / "raw" / "oh_2026" / "oh_s2_main.csv"
CURATED = ROOT / "data" / "reference" / "codons.yaml"
OUT = ROOT / "data" / "reference" / "codons_alk_full.yaml"

PROBE_LEN = 12


def main() -> None:
    s11 = pd.read_excel(XLSX, sheet_name="Table S11", header=1)
    hit = s11[s11.iloc[:, 1].astype(str).str.contains("variant1_WT", na=False)]
    if hit.empty:
        sys.exit("Table S11 has no 'variant1_WT' row")
    seq = re.sub(r"[^acgtACGT]", "", str(hit.iloc[0, 2])).upper()

    atlas = pd.read_csv(ATLAS)
    atlas = atlas[atlas.Type == "Missense"].dropna(subset=["residue"]).copy()
    atlas["residue"] = atlas.residue.astype(int)
    wt = (atlas.drop_duplicates("residue").set_index("residue")
                .AA_Change.str[0].sort_index())
    lo, hi = int(wt.index.min()), int(wt.index.max())
    if set(range(lo, hi + 1)) - set(wt.index):
        sys.exit("atlas wild-type residues are not contiguous; alignment unsafe")
    peptide = "".join(wt[r] for r in range(lo, hi + 1))

    placed = None
    for frame in (0, 1, 2):
        prot = "".join(translate(seq[i:i + 3]) for i in range(frame, len(seq) - 2, 3))
        idx = prot.find(peptide[:PROBE_LEN])
        if idx >= 0:
            if prot.count(peptide[:PROBE_LEN]) > 1:
                sys.exit("probe is not unique in the construct; alignment unsafe")
            placed = frame + 3 * idx
            break
    if placed is None:
        sys.exit("could not locate the atlas peptide in the construct")

    table, checked = {}, 0
    for k, r in enumerate(range(lo, hi + 1)):
        a, b = placed + 3 * k, placed + 3 * k + 3
        if b > len(seq):
            sys.exit(f"construct ends before residue {r}")
        codon = seq[a:b]
        if translate(codon) != wt[r]:
            sys.exit(f"residue {r}: codon {codon} encodes {translate(codon)}, "
                     f"atlas says {wt[r]}")
        table[r] = codon
        checked += 1
    print(f"aligned at construct nt {placed + 1}; "
          f"{checked}/{checked} atlas positions reproduce their wild-type residue")

    curated = yaml.safe_load(CURATED.read_text())["alk_nm_004304_5"]["codons"]
    agree, disagree = 0, []
    for r, v in curated.items():
        c = v["codon"] if isinstance(v, dict) else v
        r = int(r)
        if r in table:
            if table[r] == c:
                agree += 1
            else:
                disagree.append((r, c, table[r]))
    print(f"curated slice: {agree} agree, {len(disagree)} disagree")
    for r, old, new in disagree:
        print(f"  residue {r}: curated {old} ({translate(old)}) -> "
              f"construct {new} ({translate(new)})")

    OUT.write_text(yaml.safe_dump({
        "source": ("Oh et al. 2026 Genome Biology 27:72, Additional file 2, Table S11, "
                   "row 'EML4_ALK variant1_WT' (the construct the prime-editing atlas "
                   "was built in)."),
        "alignment": {"construct_nt_of_first_residue": placed + 1,
                      "residue_range": [lo, hi]},
        "verification": {"atlas_positions_checked": checked,
                         "all_reproduce_atlas_wildtype": True,
                         "curated_slice_agreements": agree,
                         "curated_slice_disagreements":
                             [{"residue": r, "curated": o, "construct": n}
                              for r, o, n in disagree]},
        "codons": {int(r): c for r, c in sorted(table.items())},
    }, sort_keys=False, width=100))
    print(f"\nwrote {OUT.relative_to(ROOT)} with {len(table)} codons ({lo}-{hi})")


if __name__ == "__main__":
    main()
