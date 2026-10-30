"""Codon tools — genetic code, SNV route enumeration, pyrimidine-referenced classes.

All WT codons used in this project come from RefSeq CDS of the relevant gene,
NOT from inferred fingerprints of any downstream dataset. The RefSeq lookup is
performed once (data/reference/codons.yaml) and spot-checked against ClinVar
records that fix the third base of a codon via a known c.N>N variant.

Pyrimidine-referenced substitution class notation (COSMIC/standard):

- any mutation at a C or T is reported directly: e.g. "C>A", "T>G"
- any mutation at a G or A is reported after complementing to the pyrimidine
  reading-frame-agnostic equivalent: G>A on sense = C>T on antisense.

Note: this convention makes six substitution classes exhaustive for a
double-stranded context.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

__all__ = [
    "GENETIC_CODE",
    "translate",
    "sub_class",
    "enumerate_snv_routes",
    "VARIANT_ROUTES",
    "SnvRoute",
]

# ----------------------------------------------------------------------------
# Standard genetic code
# ----------------------------------------------------------------------------

BASES: tuple[str, ...] = ("A", "C", "G", "T")

_CODON_ORDER = (
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
)
# Index order TCAG × TCAG × TCAG, canonical BioPython layout
_BASE_INDEX = {b: i for i, b in enumerate("TCAG")}


def translate(codon: str) -> str:
    """Translate a three-base DNA codon to its amino-acid single-letter code.

    Parameters
    ----------
    codon
        Three-base string over {A, C, G, T}. Case-insensitive.

    Returns
    -------
    One-letter amino-acid code (or '*' for stop).
    """
    codon = codon.upper()
    if len(codon) != 3 or any(b not in BASES for b in codon):
        raise ValueError(f"invalid codon: {codon!r}")
    i = _BASE_INDEX[codon[0]]
    j = _BASE_INDEX[codon[1]]
    k = _BASE_INDEX[codon[2]]
    return _CODON_ORDER[i * 16 + j * 4 + k]


# Expose a dict form for callers that prefer it
GENETIC_CODE: dict[str, str] = {
    a + b + c: translate(a + b + c) for a in BASES for b in BASES for c in BASES
}


# ----------------------------------------------------------------------------
# Pyrimidine-referenced substitution classes
# ----------------------------------------------------------------------------

_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def sub_class(ref: str, alt: str) -> str:
    """Return the pyrimidine-referenced 6-class label for a ref -> alt substitution.

    Examples
    --------
    >>> sub_class("C", "T")
    'C>T'
    >>> sub_class("G", "A")  # G>A equals C>T on antisense
    'C>T'
    """
    ref = ref.upper()
    alt = alt.upper()
    if ref == alt or ref not in BASES or alt not in BASES:
        raise ValueError(f"invalid substitution: {ref}>{alt}")
    if ref in "CT":
        return f"{ref}>{alt}"
    return f"{_COMP[ref]}>{_COMP[alt]}"


SUB_CLASSES: tuple[str, ...] = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")


# ----------------------------------------------------------------------------
# SNV route enumeration
# ----------------------------------------------------------------------------


@dataclass(frozen=True)
class SnvRoute:
    """A single-nucleotide route from a wild-type codon to a mutant codon.

    Attributes
    ----------
    wt_codon : str
        The three-base wild-type codon.
    mut_codon : str
        The three-base mutant codon (differs from `wt_codon` at exactly one position).
    position : int
        0-indexed position within the codon (0, 1, or 2) that is changed.
    ref_base : str
        The wild-type base at `position`.
    alt_base : str
        The mutant base at `position`.
    sub_class : str
        Pyrimidine-referenced 6-class label for the substitution.
    origin : str
        "A:T" if the modified base is A or T; "C:G" if C or G.
    """

    wt_codon: str
    mut_codon: str
    position: int
    ref_base: str
    alt_base: str
    sub_class: str
    origin: str

    @property
    def wt_aa(self) -> str:
        return translate(self.wt_codon)

    @property
    def mut_aa(self) -> str:
        return translate(self.mut_codon)


def enumerate_snv_routes(wt_codon: str, target_aa: str) -> list[SnvRoute]:
    """Return all one-nucleotide routes from `wt_codon` to a codon encoding `target_aa`.

    Multi-route variants (list length > 1) are those for which more than one SNV
    reaches an alternative codon encoding the target amino acid.

    Parameters
    ----------
    wt_codon : str
        Wild-type codon.
    target_aa : str
        Single-letter target amino acid.

    Returns
    -------
    list[SnvRoute]
        One entry per distinct single-nucleotide route.
    """
    wt_codon = wt_codon.upper()
    target_aa = target_aa.upper()
    routes: list[SnvRoute] = []
    for pos in range(3):
        ref = wt_codon[pos]
        for alt in BASES:
            if alt == ref:
                continue
            mut_codon = wt_codon[:pos] + alt + wt_codon[pos + 1 :]
            if translate(mut_codon) != target_aa:
                continue
            routes.append(
                SnvRoute(
                    wt_codon=wt_codon,
                    mut_codon=mut_codon,
                    position=pos,
                    ref_base=ref,
                    alt_base=alt,
                    sub_class=sub_class(ref, alt),
                    origin="A:T" if ref in "AT" else "C:G",
                )
            )
    return routes


def VARIANT_ROUTES(wt_codon: str, target_aa: str) -> Iterator[SnvRoute]:
    """Alias of enumerate_snv_routes returning an iterator."""
    yield from enumerate_snv_routes(wt_codon, target_aa)


# ----------------------------------------------------------------------------
# Minimal self-check (runs on import for safety; cheap)
# ----------------------------------------------------------------------------

def _self_check() -> None:
    # Known codons we rely on throughout the paper
    assert translate("CTG") == "L", "CTG should encode Leu"
    assert translate("TTC") == "F", "TTC should encode Phe"
    assert translate("ATC") == "I", "ATC should encode Ile"
    assert translate("GAG") == "E", "GAG should encode Glu"
    assert translate("GGA") == "G", "GGA should encode Gly"
    assert translate("TGC") == "C", "TGC should encode Cys"

    # Pyrimidine-referenced class consistency
    assert sub_class("C", "T") == "C>T"
    assert sub_class("G", "A") == "C>T"  # complement
    assert sub_class("T", "C") == "T>C"
    assert sub_class("A", "G") == "T>C"  # complement

    # L1196M from CTG: single route (pos-1 C>A). L1196L's silent routes excluded.
    r_m = enumerate_snv_routes("CTG", "M")
    assert len(r_m) == 1, f"L1196M expected 1 route, got {len(r_m)}"
    assert r_m[0].position == 0 and r_m[0].sub_class == "C>A"

    # F1174L from TTC: three routes (pos-1 T>C, pos-3 C>A, pos-3 C>G)
    r_fl = enumerate_snv_routes("TTC", "L")
    assert len(r_fl) == 3, f"F1174L (TTC->L) expected 3 routes, got {len(r_fl)}"
    sub_classes_seen = {r.sub_class for r in r_fl}
    assert sub_classes_seen == {"T>C", "C>A", "C>G"}, f"got {sub_classes_seen}"

    # F1174L from TTT: three routes but different (pos-1 T>C, pos-3 T>A, pos-3 T>G)
    # This is what the v0.18 paper assumed. Keeping the test to document the
    # difference and prevent regression if someone passes TTT.
    r_fl_wrong = enumerate_snv_routes("TTT", "L")
    assert len(r_fl_wrong) == 3
    wrong_classes = {r.sub_class for r in r_fl_wrong}
    assert wrong_classes == {"T>C", "T>A", "T>G"}


_self_check()
