"""Tests for src.codon_tools.

These tests mirror the known codon/variant pairs from the pre-registration
and act as a regression suite against the kind of input errors (TTT vs TTC,
ATC vs ATT, GAA vs GAG) that broke our earlier submission.
"""

from __future__ import annotations

import pytest

from src.codon_tools import (
    enumerate_snv_routes,
    sub_class,
    translate,
)


class TestTranslate:
    """Standard genetic code at the codons we rely on."""

    @pytest.mark.parametrize(
        "codon,aa",
        [
            # ALK kinase-domain WT codons (RefSeq NM_004304.5)
            ("GAG", "E"),  # E1129, E1210
            ("CCT", "P"),  # P1153
            ("TGC", "C"),  # C1156
            ("ATC", "I"),  # I1171
            ("TTC", "F"),  # F1174, F1245
            ("CTG", "L"),  # L1196 (gatekeeper)
            ("CTC", "L"),  # L1198, L1256
            ("GGA", "G"),  # G1202, G1269
            ("GAC", "D"),  # D1203
            # Known-wrong assumptions from v0.18 (kept as regression):
            ("TTT", "F"),  # also Phe but wrong codon for F1174
            ("ATT", "I"),  # also Ile but wrong codon for I1171
            ("GAA", "E"),  # also Glu but wrong codon for E1129/E1210
        ],
    )
    def test_translate(self, codon: str, aa: str) -> None:
        assert translate(codon) == aa


class TestSubClass:
    """Pyrimidine-referenced class convention."""

    def test_direct_pyrimidine(self) -> None:
        assert sub_class("C", "T") == "C>T"
        assert sub_class("T", "C") == "T>C"

    def test_complemented_purine(self) -> None:
        assert sub_class("G", "A") == "C>T"  # == complement
        assert sub_class("A", "G") == "T>C"

    def test_invalid(self) -> None:
        with pytest.raises(ValueError):
            sub_class("C", "C")
        with pytest.raises(ValueError):
            sub_class("N", "A")


class TestL1196Routes:
    """L1196 at codon CTG — three primary variants in the paper."""

    def test_l1196m_single_route(self) -> None:
        routes = enumerate_snv_routes("CTG", "M")
        assert len(routes) == 1
        assert routes[0].position == 0
        assert routes[0].sub_class == "C>A"
        assert routes[0].origin == "C:G"

    def test_l1196q_single_route(self) -> None:
        routes = enumerate_snv_routes("CTG", "Q")
        assert len(routes) == 1
        assert routes[0].position == 1
        assert routes[0].sub_class == "T>A"
        assert routes[0].origin == "A:T"

    def test_l1196p_single_route(self) -> None:
        routes = enumerate_snv_routes("CTG", "P")
        assert len(routes) == 1
        assert routes[0].position == 1
        assert routes[0].sub_class == "T>C"
        assert routes[0].origin == "A:T"


class TestF1174Routes:
    """F1174 critical: wrong codon (TTT) in v0.18 gave F1174L a different route set."""

    def test_f1174l_mixed_origin_from_ttc(self) -> None:
        """From TTC, F1174L has three routes with MIXED origins."""
        routes = enumerate_snv_routes("TTC", "L")
        assert len(routes) == 3
        classes = {r.sub_class for r in routes}
        origins = {r.origin for r in routes}
        assert classes == {"T>C", "C>A", "C>G"}
        assert origins == {"A:T", "C:G"}

    def test_f1174l_wrong_codon_regression(self) -> None:
        """From TTT (v0.18 assumption), F1174L has three different routes.

        Keeping this test as a regression against anyone passing TTT by mistake.
        """
        routes = enumerate_snv_routes("TTT", "L")
        classes = {r.sub_class for r in routes}
        origins = {r.origin for r in routes}
        assert classes == {"T>C", "T>A", "T>G"}
        assert origins == {"A:T"}

    def test_f1174c_single_route_rare_class(self) -> None:
        """F1174C from TTC is a single pos-2 T>G — rarest Mutazyme II class."""
        routes = enumerate_snv_routes("TTC", "C")
        assert len(routes) == 1
        assert routes[0].position == 1
        assert routes[0].sub_class == "T>G"
        assert routes[0].origin == "A:T"


class TestG1202Routes:
    """G1202R — two routes from GGA (pos-1 G>A and pos-1 G>C) — C:G origin."""

    def test_g1202r_two_routes_cg_origin(self) -> None:
        routes = enumerate_snv_routes("GGA", "R")
        classes = {r.sub_class for r in routes}
        origins = {r.origin for r in routes}
        assert len(routes) == 2
        assert classes == {"C>T", "C>G"}
        assert origins == {"C:G"}
