"""Per-variant supply weight computation.

For a resistance mutation v = (residue r, mutant amino acid a) in a gene whose
wild-type codon at r is c, and a mutational process d with 6-class spectrum S_d,
the per-variant supply weight is:

    S_d(v) = sum over one-nt routes from c to any codon encoding a : S_d(class of route)

If v is reachable by only one route, S_d(v) is a single class frequency. If v is
reachable by multiple routes (e.g., F1174L from TTC via T>C at pos 1, C>A at pos 3,
or C>G at pos 3), S_d(v) sums over all three.

Variants with multi-route reachability and mixed-origin routes (A:T + C:G) are
flagged "mixed" by convention and MUST be excluded from the strict A:T versus C:G
test by the pre-registered rules.
"""

from __future__ import annotations

from dataclasses import dataclass

from .codon_tools import enumerate_snv_routes, SnvRoute
from .spectrum import Spectrum

__all__ = ["SupplyRecord", "supply_weight", "variant_origin"]


@dataclass(frozen=True)
class SupplyRecord:
    """Per-variant supply result with full provenance of the routes used."""

    residue: int
    wt_aa: str
    mut_aa: str
    wt_codon: str
    routes: tuple[SnvRoute, ...]
    weight: float
    origin: str          # "A:T", "C:G", "mixed", or "none" if no route exists
    spectrum_name: str

    @property
    def variant(self) -> str:
        return f"{self.wt_aa}{self.residue}{self.mut_aa}"


def supply_weight(
    residue: int,
    wt_codon: str,
    mut_aa: str,
    spectrum: Spectrum,
) -> SupplyRecord:
    """Compute S_d(v) for one variant under one spectrum.

    Parameters
    ----------
    residue : int
        Residue number (1-indexed, in the convention of the primary literature).
    wt_codon : str
        Wild-type codon at that residue (from RefSeq).
    mut_aa : str
        Target amino acid in one-letter code.
    spectrum : Spectrum
        6-class pyrimidine-referenced spectrum.

    Returns
    -------
    SupplyRecord
    """
    routes = tuple(enumerate_snv_routes(wt_codon, mut_aa))
    weight = sum(spectrum[r.sub_class] for r in routes)
    origin = variant_origin(routes)
    from .codon_tools import translate
    return SupplyRecord(
        residue=residue,
        wt_aa=translate(wt_codon),
        mut_aa=mut_aa,
        wt_codon=wt_codon,
        routes=routes,
        weight=weight,
        origin=origin,
        spectrum_name=spectrum.name,
    )


def variant_origin(routes: tuple[SnvRoute, ...]) -> str:
    """Return 'A:T', 'C:G', 'mixed', or 'none' based on routes' base-pair origins."""
    if not routes:
        return "none"
    origins = {r.origin for r in routes}
    if origins == {"A:T"}:
        return "A:T"
    if origins == {"C:G"}:
        return "C:G"
    return "mixed"
