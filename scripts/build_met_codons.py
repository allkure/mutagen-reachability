#!/usr/bin/env python3
"""Build the MET kinase-domain codon table from RefSeq NM_000245.2.

Why this exact record
---------------------
Tiedt et al. 2011 Supplementary Tables S1 and S2 state that their amino-acid
and base-pair numbering refers to NM_000245.2, and they report every hit at
the nucleotide level.  Using the same record means the audit reads their
substitution classes directly off their own coordinates instead of inferring
them.

Three independent verifications, all of which must pass
-------------------------------------------------------
1. The CDS translated from the parsed sequence must equal the /translation
   field recorded in the same GenBank entry.
2. The codon at every one of the 287 kinase-domain positions must encode the
   wild-type residue that the Estevam DMS reports there.
3. Every nucleotide alteration Tiedt reports must fall inside the codon of
   the residue they assign it to, must start from the base this table holds,
   and must produce the amino acid they report.

Coordinates
-----------
The GenBank CDS is 188..4360 in mRNA coordinates, so CDS position p is mRNA
position 187 + p, and residue r occupies CDS positions 3r-2 .. 3r.  Tiedt's
"a3689g" style coordinates are CDS positions, which verification 3 confirms.

Outputs
-------
data/reference/codons_met.yaml
data/raw/met_refseq/met_cds.txt
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

GB = ROOT / "data" / "raw" / "met_refseq" / "NM_000245.2.gb"
DMS = ROOT / "data" / "raw" / "met_estevam2024" / "met_dms_long.csv"
OUT_YAML = ROOT / "data" / "reference" / "codons_met.yaml"
OUT_CDS = ROOT / "data" / "raw" / "met_refseq" / "met_cds.txt"

DOMAIN = (1059, 1345)          # kinase-domain span covered by the DMS

# Tiedt 2011 Supplementary Tables S1 and S2, transcribed verbatim.
# (residue, aa_change, cds_pos, ref, alt, clones)
TIEDT_S1 = [          # NVP-BVU972, 85 clones, 10 alterations in 4 residues
    (1230, "Y1230C", 3689, "a", "g", 31),
    (1230, "Y1230H", 3688, "t", "c", 20),
    (1230, "Y1230N", 3688, "t", "a", 19),
    (1230, "Y1230S", 3689, "a", "c", 4),
    (1230, "Y1230D", 3688, "t", "g", 1),
    (1228, "D1228N", 3682, "g", "a", 5),
    (1228, "D1228E", 3684, "c", "a", 2),
    (1228, "D1228A", 3683, "a", "c", 1),
    (1155, "V1155L", 3463, "g", "t", 1),
    (1200, "F1200I", 3598, "t", "a", 1),
]
TIEDT_S2 = [          # AMG 458, 60 clones, 15 alterations / 13 aa changes / 10 residues
    (1200, "F1200L", 3598, "t", "c", 17),
    (1200, "F1200L", 3600, "t", "a", 8),
    (1200, "F1200L", 3600, "t", "g", 3),
    (1200, "F1200I", 3598, "t", "a", 8),
    (1250, "M1250T", 3749, "t", "c", 6),
    (1211, "M1211L", 3631, "a", "c", 5),
    (1195, "L1195I", 3583, "c", "a", 3),
    (1195, "L1195V", 3583, "c", "g", 1),
    (1163, "G1163V", 3488, "g", "t", 2),
    (1163, "G1163E", 3488, "g", "a", 1),
    (1155, "V1155L", 3463, "g", "c", 1),
    (1159, "Y1159N", 3475, "t", "a", 2),
    (1193, "K1193I", 3578, "a", "t", 1),
    (1228, "D1228A", 3683, "a", "c", 1),
    (1258, "Q1258L", 3773, "a", "t", 1),
]


def parse_genbank(text: str) -> tuple[str, str, tuple[int, int]]:
    m = re.search(r"^\s+CDS\s+(\d+)\.\.(\d+)", text, re.M)
    if not m:
        sys.exit("no CDS feature found")
    cds_span = (int(m.group(1)), int(m.group(2)))

    tm = re.search(r'/translation="([^"]+)"', text, re.S)
    if not tm:
        sys.exit("no /translation field found")
    protein = re.sub(r"\s+", "", tm.group(1))

    om = re.search(r"^ORIGIN\s*$(.*?)^//", text, re.S | re.M)
    if not om:
        sys.exit("no ORIGIN block found")
    seq = re.sub(r"[^acgtn]", "", om.group(1).lower())
    return seq, protein, cds_span


def main() -> None:
    text = GB.read_text()
    mrna, protein_ref, (c0, c1) = parse_genbank(text)
    print(f"mRNA length {len(mrna):,}   CDS {c0}..{c1}   /translation {len(protein_ref)} aa")

    cds = mrna[c0 - 1:c1]
    if len(cds) % 3 != 0:
        sys.exit(f"CDS length {len(cds)} is not a multiple of 3")

    codons = [cds[i:i + 3].upper() for i in range(0, len(cds), 3)]
    protein = "".join(translate(c) for c in codons)
    if not protein.endswith("*"):
        sys.exit("CDS does not end in a stop codon")
    protein = protein[:-1]

    # --- verification 1: translation matches the record's own protein -------
    if protein != protein_ref:
        for i, (a, b) in enumerate(zip(protein, protein_ref), start=1):
            if a != b:
                sys.exit(f"translation mismatch at residue {i}: got {a}, record has {b}")
        sys.exit(f"translation length mismatch: {len(protein)} vs {len(protein_ref)}")
    print(f"  [1/3] translation matches the record's /translation exactly ({len(protein)} aa)")

    table = {r: codons[r - 1] for r in range(DOMAIN[0], DOMAIN[1] + 1)}

    # --- verification 2: codons encode the DMS wild-type residues -----------
    dms = pd.read_csv(DMS).drop_duplicates("residue")[["residue", "wildtype"]]
    dms = dms[dms.residue.between(*DOMAIN)]
    bad = [(int(r.residue), r.wildtype, translate(table[int(r.residue)]))
           for r in dms.itertuples() if translate(table[int(r.residue)]) != r.wildtype]
    if bad:
        sys.exit(f"{len(bad)} codons disagree with the DMS wild-type residue, e.g. {bad[:5]}")
    print(f"  [2/3] all {len(dms)} DMS kinase-domain positions agree with the codon table")

    # --- verification 3: every Tiedt alteration is consistent ---------------
    checked = 0
    for label, rows in (("S1", TIEDT_S1), ("S2", TIEDT_S2)):
        for res, aa_change, pos, ref, alt, _ in rows:
            lo = 3 * res - 2
            if not lo <= pos <= lo + 2:
                sys.exit(f"{label} {aa_change}: CDS {pos} is outside codon {res} ({lo}-{lo+2})")
            codon = table[res]
            idx = pos - lo
            if codon[idx].lower() != ref:
                sys.exit(f"{label} {aa_change}: table has {codon[idx]} at CDS {pos}, "
                         f"Tiedt reports a {ref}")
            mutated = codon[:idx] + alt.upper() + codon[idx + 1:]
            got = translate(mutated)
            want = aa_change[-1]
            if got != want:
                sys.exit(f"{label} {aa_change}: {codon}->{mutated} gives {got}, not {want}")
            checked += 1
    print(f"  [3/3] all {checked} Tiedt nucleotide alterations are consistent with the table")

    OUT_CDS.write_text(cds.upper() + "\n")
    payload = {
        "source": ("RefSeq NM_000245.2 (MET transcript variant 2), GenBank record as "
                   "retrieved by the author; CDS 188..4360. This is the record Tiedt "
                   "et al. 2011 Supplementary Tables S1 and S2 state their amino-acid "
                   "and base-pair numbering refers to."),
        "verification": {
            "translation_matches_record": True,
            "dms_wildtype_positions_checked": int(len(dms)),
            "tiedt_alterations_checked": checked,
        },
        "cds_span_mrna": [c0, c1],
        "domain": list(DOMAIN),
        "codons": {int(r): c for r, c in sorted(table.items())},
    }
    OUT_YAML.write_text(yaml.safe_dump(payload, sort_keys=False, width=100))
    print(f"\nwrote {OUT_YAML.relative_to(ROOT)} with {len(table)} codons "
          f"({DOMAIN[0]}-{DOMAIN[1]})")
    print("spot check:", {r: table[r] for r in (1155, 1163, 1195, 1200, 1228, 1230, 1250)})


if __name__ == "__main__":
    main()
