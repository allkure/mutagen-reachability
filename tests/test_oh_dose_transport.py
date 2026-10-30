"""Regression tests for the Oh supplementary extraction and dose analysis.

These lock down the facts the manuscript's §2.6 claims rest on, so that a
change to the extraction or the aggregation cannot silently move a reported
number.  They run against the committed derived files, not against a live
recomputation, which is the point: if a rerun changes them, the test fails
and the number in the paper must be revisited.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "oh_2026"
DERIVED = ROOT / "data" / "derived"


@pytest.fixture(scope="module")
def dose() -> dict:
    p = DERIVED / "oh_dose_transport.json"
    if not p.exists():
        pytest.skip("run scripts/oh_dose_transport.py first")
    return json.loads(p.read_text())


class TestExtraction:
    """The extracted tables must match the published workbook's shape."""

    @pytest.mark.parametrize("fname,n_rows,n_missense,res_lo,res_hi", [
        ("oh_s2_main.csv", 3208, 2173, 1058, 1388),
        ("oh_s3_lorlatinib_doses.csv", 415, 278, 1172, 1215),
        ("oh_s5_three_drug_doses.csv", 420, 281, 1172, 1215),
        ("oh_s9_two_doses.csv", 645, 427, 1151, 1215),
    ])
    def test_shape(self, fname, n_rows, n_missense, res_lo, res_hi):
        d = pd.read_csv(RAW / fname)
        assert len(d) == n_rows
        assert int((d.Type == "Missense").sum()) == n_missense
        r = d.residue.dropna().astype(int)
        assert (r.min(), r.max()) == (res_lo, res_hi)

    def test_s2_matches_the_table_used_elsewhere(self):
        """The new extraction must reproduce the S2 CSV the pipeline already uses."""
        old = pd.read_csv(RAW / "table_s2_resistance_scores.csv")
        new = pd.read_csv(RAW / "oh_s2_main.csv")
        m = (old.merge(new, left_on="aa_change", right_on="AA_Change")
                .drop_duplicates("aa_change"))
        assert len(m) > 2500
        for a, b in [("alec_score", "Alectinib_Resistance_Score"),
                     ("lor_score", "Lorlatinib_Resistance_Score"),
                     ("zot_score", "Zotizalkib_Resistance_Score"),
                     ("fitness_score", "Fitness_Score")]:
            assert (m[a] - m[b]).abs().max() == pytest.approx(0.0, abs=1e-9)

    def test_no_missense_score_column_is_empty(self):
        for fname in ("oh_s3_lorlatinib_doses.csv", "oh_s5_three_drug_doses.csv",
                      "oh_s9_two_doses.csv"):
            d = pd.read_csv(RAW / fname)
            mis = d[d.Type == "Missense"]
            for c in [c for c in d.columns if c.endswith("Resistance_Score")]:
                assert mis[c].notna().sum() > 200, f"{fname}:{c}"


class TestDoseFindings:
    """The three claims §2.6 makes from the internal dose series."""

    def test_ranking_degrades_across_the_dose_span(self, dose):
        lo, hi = dose["headline"]["rho_low_vs_high_range"]
        assert 0.25 < lo < 0.40
        assert 0.55 < hi < 0.70
        assert dose["headline"]["rho_low_vs_high_median"] < 0.45

    def test_adjacent_steps_agree_better_than_the_full_span(self, dose):
        assert (dose["headline"]["rho_adjacent_median"]
                > dose["headline"]["rho_low_vs_high_median"])

    def test_cross_drug_matched_dose_beats_same_drug_across_span(self, dose):
        """The control that keeps us from calling the dose effect purely biological."""
        assert (dose["headline"]["cross_drug_matched_median_rho"]
                > dose["headline"]["rho_low_vs_high_median"])

    def test_matched_dose_rank_beats_unmatched(self, dose):
        b = dose["controls"]["B_cross_drug"]
        assert (b["matched_dose_rank"]["median_rho"]
                > b["unmatched_dose_rank"]["median_rho"])

    def test_topk_overlap_is_far_below_complete(self, dose):
        lo, hi = dose["headline"]["topk_overlap_range"]
        assert 0 < lo <= hi < dose["top_k"]


class TestConfounds:
    """Controls that rule out the two cheap explanations."""

    def test_sub_library_is_not_range_restricted(self, dose):
        """Attenuation by restricted range would invalidate the dose claim."""
        d = dose["controls"]["D_range_restriction"]
        assert d["range_restricted"] is False
        for drug, v in d["by_drug"].items():
            assert v["sd_ratio_window_over_full"] >= 1.0, drug

    def test_dispersion_grows_with_dose_in_the_primary_series(self, dose):
        """Low dose is the compressed one, so the atlas condition is the least
        discriminating. True for every S5 and S9 series; S3 is the exception."""
        c = dose["controls"]["C_sd_monotone_increasing"]
        assert all(v for k, v in c.items() if k.startswith(("S5", "S9")))

    def test_atlas_matches_the_low_dose_column_for_alec_and_lor(self, dose):
        """Establishes which column is the condition used as sigma elsewhere."""
        a = dose["controls"]["A_atlas_vs_dose"]
        assert a["alectinib"]["best_match_dose"] == "20.0"
        assert a["lorlatinib"]["best_match_dose"] == "5.0"

    def test_atlas_agreement_with_its_own_repeat_is_moderate(self, dose):
        """S2 vs the matched S5 condition: two runs of nominally the same assay."""
        a = dose["controls"]["A_atlas_vs_dose"]
        best = [max(v["rho_by_dose"].values()) for v in a.values()]
        assert all(0.40 < b < 0.70 for b in best), best
