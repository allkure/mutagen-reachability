"""Regression tests for the Estevam 2024 MET DMS extraction.

The MET audit is only as good as this extraction, and two things about the
source repository make it easy to get wrong: the kinase domain is numbered
locally 1-287, and the repository's own scripts disagree about which file
holds which background.  These tests pin down both, plus the sigma
definition.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
MET = ROOT / "data" / "raw" / "met_estevam2024"

# Residues whose wild-type identity is fixed by their own variant nomenclature
# in the MET resistance literature.
ANCHORS = {1092: "V", 1094: "H", 1163: "G", 1195: "L", 1200: "F",
           1211: "M", 1228: "D", 1230: "Y", 1234: "Y", 1235: "Y"}


def _need(p: Path):
    if not p.exists():
        pytest.skip("run scripts/extract_met_dms.py first")


@pytest.fixture(scope="module")
def long_df() -> pd.DataFrame:
    _need(MET / "met_dms_long.csv")
    return pd.read_csv(MET / "met_dms_long.csv")


@pytest.fixture(scope="module")
def wide() -> pd.DataFrame:
    _need(MET / "met_dms_wide.csv")
    return pd.read_csv(MET / "met_dms_wide.csv")


@pytest.fixture(scope="module")
def meta() -> dict:
    _need(MET / "met_dms_meta.json")
    return json.loads((MET / "met_dms_meta.json").read_text())


class TestNumbering:
    """The +1058 offset must be re-derivable, not assumed."""

    def test_offset_is_1058_for_both_backgrounds(self, meta):
        assert set(meta["offset_applied"].values()) == {1058}

    def test_canonical_range(self, meta):
        assert meta["canonical_residue_range"] == [1059, 1345]
        assert meta["n_positions"] == 287

    def test_every_anchor_residue_has_the_expected_wildtype(self, long_df):
        wt = (long_df.drop_duplicates("residue")
                     .set_index("residue").wildtype.to_dict())
        for res, aa in ANCHORS.items():
            assert wt.get(res) == aa, f"{aa}{res} got {wt.get(res)}"

    def test_a_shifted_offset_would_fail(self, long_df):
        """Guards against an off-by-one slipping through unnoticed."""
        wt = (long_df.drop_duplicates("residue")
                     .set_index("residue").wildtype.to_dict())
        for shift in (-1, 1):
            hits = sum(1 for r, aa in ANCHORS.items() if wt.get(r + shift) == aa)
            assert hits < len(ANCHORS) / 2


class TestStructure:
    def test_both_backgrounds_present(self, long_df):
        assert set(long_df.background.unique()) == {"WT", "exon14del"}

    def test_inhibitor_set(self, meta):
        assert set(meta["inhibitors"]) == {
            "A458", "Cabo", "Camp", "Crizo", "DMSO", "Gle",
            "Glu", "Mere", "NVP", "Savo", "Tepo", "Tiv"}

    def test_tiedt_inhibitors_are_in_the_atlas(self, meta):
        """The MET audit only works because Estevam assayed Tiedt's two drugs."""
        assert meta["tiedt_2011_inhibitors_present"] == {
            "NVP": "NVP-BVU972", "A458": "AMG 458"}

    def test_aa_change_is_consistent_with_its_parts(self, long_df):
        s = long_df.sample(2000, random_state=0)
        built = s.wildtype + s.residue.astype(str) + s.mutation
        assert (built == s.aa_change).all()

    def test_one_wildtype_per_residue(self, long_df):
        assert long_df.groupby("residue").wildtype.nunique().max() == 1


class TestSigmaDefinition:
    """sigma = drug - DMSO, the MET analogue of Oh's Resistance score."""

    def test_delta_columns_exist_for_every_drug(self, wide):
        drugs = ["A458", "Cabo", "Camp", "Crizo", "Gle", "Glu",
                 "Mere", "NVP", "Savo", "Tepo", "Tiv"]
        for d in drugs:
            assert f"delta_{d}" in wide.columns

    def test_delta_is_exactly_drug_minus_dmso(self, wide):
        s = wide.dropna(subset=["Crizo", "DMSO"]).head(500)
        assert ((s.delta_Crizo - (s.Crizo - s.DMSO)).abs() < 1e-9).all()

    def test_dmso_is_not_orthogonal_to_the_raw_drug_score(self, wide):
        """The reason the subtraction is needed at all. If this ever drops to
        ~0, the MET arms have changed and the sigma choice must be revisited."""
        s = wide.dropna(subset=["DMSO", "Crizo"])
        assert s.DMSO.corr(s.Crizo) > 0.6

    def test_subtraction_removes_most_of_the_growth_component(self, wide):
        s = wide.dropna(subset=["DMSO", "Crizo"])
        assert abs(s.DMSO.corr(s.delta_Crizo)) < abs(s.DMSO.corr(s.Crizo))


class TestKnownBiology:
    """The atlas must recover the MET resistance variants seen in patients."""

    KNOWN = ["D1228N", "D1228H", "D1228V", "D1228Y",
             "Y1230C", "Y1230H", "Y1230D", "Y1230S", "G1163R"]

    def test_activation_loop_variants_rank_high_under_a_type_I_inhibitor(self, wide):
        d = wide[wide.background == "exon14del"].dropna(subset=["delta_Crizo"]).copy()
        d["rank"] = d.delta_Crizo.rank(ascending=False, method="min")
        hit = d[d.aa_change.isin(self.KNOWN)]
        assert len(hit) >= 8
        assert hit["rank"].median() < len(d) * 0.05

    def test_type_II_inhibitor_is_less_affected_by_activation_loop_variants(self, wide):
        """D1228/Y1230 confer less resistance to a type II inhibitor; if this
        inverts, the inhibitor columns have probably been mislabelled."""
        d = wide[wide.background == "WT"].dropna(
            subset=["delta_Crizo", "delta_A458"]).copy()
        hit = d[d.aa_change.str.startswith(("D1228", "Y1230"))]
        assert len(hit) >= 4
        assert hit.delta_Crizo.median() > hit.delta_A458.median()
