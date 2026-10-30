"""Derive substitution classes for Bradeen 2006 imatinib variants.

For each amino-acid change v in the Bradeen imatinib arm (Table 1), this script
enumerates every WT codon that encodes the WT amino acid, finds all one-nt
routes to the mutant amino acid, and reports the SET OF SUBSTITUTION CLASSES
reachable.

For 18 of the 19 non-zero Bradeen imatinib variants, the substitution class is
FORCED by the amino-acid change itself — the only one-nt routes from any WT
codon to the mutant amino acid all share the same pyrimidine-referenced class.
These do not need RefSeq-CDS lookup to label correctly.

For Q252H and F317L, multiple substitution classes are reachable and the WT
codon identity matters. We enumerate both scenarios and leave disambiguation
to NM_005157.6 CDS verification (or Bradeen Table S3, which gives direct
nucleotide calls per isolate).

Four ClinVar anchors pin specific WT codons:
    Y253 WT = TAC    (ClinVar c.757T>C -> Y253H)
    E255 WT = GAG|GAA (pos2 A; both consistent with c.764A>T -> E255V)
    T315 WT = ACT    (ClinVar c.944C>T -> T315I)
    H396 WT = CAC|CAT (pos2 A; both consistent with c.1187A>G -> H396R)

Output:
    - stdout: human-readable table of (variant, forced_class?, routes...)
    - CSV at bradeen_imatinib_variant_routes.csv
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

# Add src/ to the import path.
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from codon_tools import GENETIC_CODE, enumerate_snv_routes  # noqa: E402

# ----------------------------------------------------------------------------
# All codons for each amino acid, from the standard genetic code
# ----------------------------------------------------------------------------
AA_TO_CODONS: dict[str, list[str]] = {}
for codon, aa in GENETIC_CODE.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)

# ----------------------------------------------------------------------------
# Bradeen 2006 Table 1 imatinib arm variants.
# Only entries with incidence_n > 0 are useful for the supply-weight ratio test.
# Codons constrained by ClinVar anchors are listed; otherwise all synonymous
# codons for the WT amino acid are enumerated below.
# ----------------------------------------------------------------------------
BRADEEN_VARIANTS = [
    ("M244V", "M", "V", None),       # Met only has ATG -> locked
    ("L248R", "L", "R", None),       # L248R only via CTN -> locks L248 WT = CTN
    # Cascade: L248R forces L248 WT codon to CTN (TTA/TTG cannot reach R in 1 nt).
    # So L248V WT codon = CTN -> class forced to C>G (pos1).
    ("L248V", "L", "V", "CTN"),
    ("G250E", "G", "E", None),
    ("Q252H", "Q", "H", None),
    ("Y253H", "Y", "H", "TAC"),      # ClinVar anchor
    ("E255K", "E", "K", None),       # WT codon GAG or GAA — pos3 ambiguous
    ("E255V", "E", "V", None),       # same ambiguity
    ("D276G", "D", "G", None),
    ("F311I", "F", "I", None),
    ("F311V", "F", "V", None),
    ("T315I", "T", "I", "ACT"),      # ClinVar anchor
    ("F317L", "F", "L", None),       # KNOWN MULTI-ROUTE
    ("M351T", "M", "T", "ATG"),      # Met only ATG
    ("E355G", "E", "G", None),
    ("F359C", "F", "C", None),
    ("F359V", "F", "V", None),
    ("V379I", "V", "I", None),
    ("L384M", "L", "M", "CTG"),      # Only CTG -> ATG is one-nt
    ("H396R", "H", "R", None),       # WT CAC or CAT — pos3 ambiguous, no effect on class
]


def main() -> None:
    rows = []
    print(
        f"{'variant':<8}  {'WT_codon(s)':<16}  "
        f"{'routes':<48}  {'classes_set':<16}  {'forced?':<8}"
    )
    print("-" * 110)

    for variant, wt_aa, mut_aa, forced_codon in BRADEEN_VARIANTS:
        if forced_codon is not None:
            # Expand CTN/GTN-style shorthand to the specific codons it covers.
            if "N" in forced_codon:
                prefix = forced_codon.replace("N", "")
                position_of_N = forced_codon.index("N")
                candidate_codons = []
                for b in "ACGT":
                    codon = (
                        forced_codon[:position_of_N] + b + forced_codon[position_of_N + 1 :]
                    )
                    if GENETIC_CODE[codon] == wt_aa:
                        candidate_codons.append(codon)
            else:
                candidate_codons = [forced_codon]
        else:
            # Only WT codons that CAN reach mut_aa via a single-nt change
            # are informative. A Val codon that cannot reach Ile in one step
            # tells us the WT codon identity, not the supply, so we drop it.
            candidate_codons = []
            for codon in AA_TO_CODONS[wt_aa]:
                if enumerate_snv_routes(codon, mut_aa):
                    candidate_codons.append(codon)

        all_routes: list[tuple[str, str, int, str]] = []
        all_classes: set[str] = set()
        n_routes_per_codon: list[int] = []
        for codon in candidate_codons:
            codon_routes = enumerate_snv_routes(codon, mut_aa)
            n_routes_per_codon.append(len(codon_routes))
            for route in codon_routes:
                all_routes.append(
                    (
                        route.wt_codon,
                        route.mut_codon,
                        route.position + 1,  # 1-indexed for human reading
                        route.sub_class,
                    )
                )
                all_classes.add(route.sub_class)

        # Per-WT-codon route count. For class-forced variants this is typically 1
        # across all synonymous candidates (and MUST be consistent for the supply
        # calculation to be well-defined without CDS lookup). We flag a mismatch
        # across candidates as a data integrity warning.
        if n_routes_per_codon:
            min_r = min(n_routes_per_codon)
            max_r = max(n_routes_per_codon)
            if min_r != max_r:
                routes_per_codon_note = f"VARIES {min_r}..{max_r}"
            else:
                routes_per_codon_note = f"{min_r}"
        else:
            routes_per_codon_note = "0"

        route_desc = "; ".join(
            f"{wc}->{mc}@p{p}:{cls}" for wc, mc, p, cls in all_routes
        )
        forced = "yes" if len(all_classes) == 1 else "no"

        print(
            f"{variant:<8}  {','.join(candidate_codons):<16}  "
            f"{route_desc[:48]:<48}  {','.join(sorted(all_classes)):<16}  {forced:<8}"
        )

        rows.append(
            {
                "variant_ia": variant,
                "wt_aa": wt_aa,
                "mut_aa": mut_aa,
                "candidate_wt_codons": "|".join(candidate_codons),
                "n_routes_total_across_candidates": len(all_routes),
                "n_routes_per_wt_codon": routes_per_codon_note,
                "substitution_classes_set": "|".join(sorted(all_classes)),
                "class_forced_by_aa_change": "yes" if len(all_classes) == 1 else "no",
                "forced_class": next(iter(all_classes)) if len(all_classes) == 1 else "",
                "routes_detail": route_desc,
                "clinvar_anchor": "yes" if forced_codon else "no",
            }
        )

    out_path = SCRIPT_DIR / "bradeen_imatinib_variant_routes.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {out_path}")

    forced_count = sum(1 for r in rows if r["class_forced_by_aa_change"] == "yes")
    multiroute = [r for r in rows if r["class_forced_by_aa_change"] == "no"]
    print(f"\n{forced_count}/{len(rows)} variants have class forced by amino-acid change")
    print(f"{len(multiroute)} multi-class variants: " + ", ".join(r["variant_ia"] for r in multiroute))


if __name__ == "__main__":
    main()
