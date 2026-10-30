"""Derive substitution classes for Azam 2003 (BCR-ABL1 + XL1-Red) variants.

Same methodology as the Bradeen 2006 derivation:
 - For each distinct amino-acid change v, enumerate every WT codon encoding
   the WT residue, find all one-nt routes to the mutant residue, and report
   the set of pyrimidine-referenced substitution classes reachable.
 - A variant whose class set has size 1 is "class-forced by the amino-acid
   change" — no RefSeq CDS lookup needed to label it.
 - Multi-class variants (size > 1) need the WT codon identity; we flag them
   for RefSeq CDS verification (NM_005157.6).

Azam 2003 reports variants at 5 μM and 10 μM STI-571. For the T2 primary
pairing (Azam 5 μM vs Bradeen 2+4 μM pooled, both low-stringency broad-
permissive-set regimes), we focus on Table 1 (kinase domain) and expand
compound notation "A/B/C" into individual sub-variants.

Expansion rules:
 - "Y253H/C/F" at a given concentration contributes 3 sub-variants, each
   inheriting the row's total isolate count divided equally if the row
   doesn't say otherwise. If the paper notes a specific IC50 for one of
   them (e.g., "17.7 H"), the dominant sub-variant is marked but the raw
   per-sub-variant count is still ambiguous without the paper's own
   per-isolate breakdown.
 - For class-forced variants, this doesn't change the supply weight; for
   multi-class variants, we propagate the ambiguity.

For now, this script enumerates only the DISTINCT sub-variants and their
class forced-ness status. Count assignment is handled separately.

Shared with Bradeen (for ratio test):
    M244V, L248R, L248V, G250E, Q252H, Y253H, E255K, D276G, F311L, T315I,
    F317L, M351T, V379I, L384M, H396R.
We include all Azam Table 1 kinase-domain variants, not only shared ones,
because the full S_XL1R spectrum is computable from Azam's full output.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from codon_tools import GENETIC_CODE, enumerate_snv_routes  # noqa: E402

AA_TO_CODONS: dict[str, list[str]] = {}
for codon, aa in GENETIC_CODE.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)


# ----------------------------------------------------------------------------
# Azam 2003 variants (Table 1, kinase domain, converted to type Ia).
# Compound notations expanded to individual sub-variants.
# ClinVar/cascade-forced WT codons noted where available.
# ----------------------------------------------------------------------------
AZAM_VARIANTS = [
    # P-loop
    ("M244V", "M", "V", "ATG"),      # Met only ATG
    ("M244I", "M", "I", "ATG"),
    ("L248R", "L", "R", "CTN"),      # L248R forces WT = CTN (same as Bradeen)
    ("L248V", "L", "V", "CTN"),
    ("G250E", "G", "E", None),
    ("G250R", "G", "R", None),
    ("Q252H", "Q", "H", None),       # multi-class
    ("Y253H", "Y", "H", "TAC"),      # ClinVar anchor
    ("Y253C", "Y", "C", "TAC"),      # TAC->TGC pos2 A>C (= T>G pyr)
    ("Y253F", "Y", "F", "TAC"),      # TAC->TTC pos2 A>T (= T>A pyr)
    ("E255K", "E", "K", None),
    # C-helix
    ("Y257C", "Y", "C", None),
    ("E258D", "E", "D", None),
    ("S265T", "S", "T", None),       # Ser has 6 codons: TCN + AGC/AGT
    ("S265I", "S", "I", None),
    ("L266M", "L", "M", None),
    ("L266V", "L", "V", None),
    ("V268A", "V", "A", None),
    ("V270A", "V", "A", None),
    ("A269V", "A", "V", None),
    ("E275K", "E", "K", None),
    ("D276V", "D", "V", None),
    ("M278L", "M", "L", "ATG"),
    ("E279K", "E", "K", None),
    ("E281K", "E", "K", None),
    ("E282D", "E", "D", None),
    ("F283L", "F", "L", None),       # multi-class, like F317L
    ("L284F", "L", "F", None),
    ("V289S", "V", "S", None),
    ("M290L", "M", "L", "ATG"),
    ("M290T", "M", "T", "ATG"),
    ("K291E", "K", "E", None),
    ("K291R", "K", "R", None),
    ("E292Q", "E", "Q", None),
    # SH3-contact
    ("K294R", "K", "R", None),
    ("Q300H", "Q", "H", None),       # multi-class like Q252H
    ("L301F", "L", "F", None),
    # Drug-contact
    ("F311L", "F", "L", None),       # multi-class
    ("F311V", "F", "V", None),
    ("T315I", "T", "I", "ACT"),      # ClinVar anchor
    ("T315S", "T", "S", "ACT"),      # ACT->TCT pos1 A>T (= T>A pyr); or AGT via pos2
    ("T315G", "T", "G", "ACT"),      # ACT->GCT not G; ACT->GGT 2 changes; no 1-nt route to G
    ("E316D", "E", "D", None),
    ("F317L", "F", "L", None),       # multi-class
    # SH2-contact
    ("G321W", "G", "W", None),
    ("N331S", "N", "S", None),
    ("V338G", "V", "G", None),
    ("V339A", "V", "A", None),
    ("V339G", "V", "G", None),
    ("M343T", "M", "T", "ATG"),
    ("A344V", "A", "V", None),
    ("Q346H", "Q", "H", None),       # multi-class
    ("M351T", "M", "T", "ATG"),
    ("M351I", "M", "I", "ATG"),
    ("E352K", "E", "K", None),
    ("F359C", "F", "C", None),
    ("I360F", "I", "F", None),
    ("A366D", "A", "D", None),
    ("G372R", "G", "R", None),
    ("E373K", "E", "K", None),
    # Activation loop
    ("V379E", "V", "E", None),
    ("V379A", "V", "A", None),
    ("L384M", "L", "M", "CTG"),
    ("M388I", "M", "I", "ATG"),
    ("H396R", "H", "R", None),
    ("H396P", "H", "P", None),
    ("G398R", "G", "R", None),
    # SH2-contact (extended)
    ("Y440C", "Y", "C", None),
    ("E450K", "E", "K", None),
    ("L451M", "L", "M", "CTG"),
    ("G463D", "G", "D", None),
    ("M472I", "M", "I", "ATG"),
    ("R473L", "R", "L", None),
    ("F486S", "F", "S", None),
    ("E494A", "E", "A", None),
    ("E499K", "E", "K", None),
    ("E499I", "E", "I", None),
    ("I502M", "I", "M", None),
    ("E509D", "E", "D", None),
]


def expand_shorthand(codon_shorthand: str, wt_aa: str) -> list[str]:
    """Expand CTN/GTN-type shorthand to the specific codons encoding wt_aa."""
    if codon_shorthand is None:
        return []
    if "N" in codon_shorthand:
        pos = codon_shorthand.index("N")
        out = []
        for b in "ACGT":
            codon = codon_shorthand[:pos] + b + codon_shorthand[pos + 1 :]
            if GENETIC_CODE[codon] == wt_aa:
                out.append(codon)
        return out
    return [codon_shorthand]


def main() -> None:
    rows = []
    forced_count = 0
    unreachable: list[str] = []
    multiclass: list[str] = []

    for variant, wt_aa, mut_aa, forced_codon in AZAM_VARIANTS:
        if forced_codon is not None:
            candidate_codons = expand_shorthand(forced_codon, wt_aa)
        else:
            candidate_codons = []
            for codon in AA_TO_CODONS[wt_aa]:
                if enumerate_snv_routes(codon, mut_aa):
                    candidate_codons.append(codon)

        all_classes: set[str] = set()
        n_routes_per_codon: list[int] = []
        for codon in candidate_codons:
            codon_routes = enumerate_snv_routes(codon, mut_aa)
            n_routes_per_codon.append(len(codon_routes))
            for route in codon_routes:
                all_classes.add(route.sub_class)

        if not any(n_routes_per_codon):
            unreachable.append(variant)
            continue

        if n_routes_per_codon and (min(n_routes_per_codon) != max(n_routes_per_codon)):
            routes_per_codon_note = (
                f"VARIES {min(n_routes_per_codon)}..{max(n_routes_per_codon)}"
            )
        else:
            routes_per_codon_note = (
                str(n_routes_per_codon[0]) if n_routes_per_codon else "0"
            )

        if len(all_classes) == 1:
            forced_count += 1
            forced_class = next(iter(all_classes))
        else:
            multiclass.append(variant)
            forced_class = ""

        rows.append(
            {
                "variant_ia": variant,
                "wt_aa": wt_aa,
                "mut_aa": mut_aa,
                "candidate_wt_codons": "|".join(candidate_codons),
                "n_routes_total_across_candidates": sum(n_routes_per_codon),
                "n_routes_per_wt_codon": routes_per_codon_note,
                "substitution_classes_set": "|".join(sorted(all_classes)),
                "class_forced_by_aa_change": "yes" if len(all_classes) == 1 else "no",
                "forced_class": forced_class,
                "clinvar_or_cascade_anchor": "yes" if forced_codon else "no",
            }
        )

    out_path = SCRIPT_DIR / "azam_variant_routes.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Reachable variants: {len(rows)} / {len(AZAM_VARIANTS)} enumerated")
    print(f"Class-forced:       {forced_count}")
    print(f"Multi-class:        {len(multiclass)}: {', '.join(multiclass)}")
    if unreachable:
        print(f"Unreachable in 1 nt (filter these out of count data): {unreachable}")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
