"""The manuscript's tables are hand-typed LaTeX. Check them against source.

Table 2 carries 96 numbers and Table 4 another 15. Nothing checked either
until the pre-submission pass, and the equivalent prose numbers in Results
3.4 turned out to be stale by then. These tests close that gap so a
recomputation that moves a cell fails here rather than in review.

Table 2's rows are matched on their labels, not their order: the manuscript
deliberately omits two supply-alone rows that benchmark_r3.csv contains, so a
positional comparison silently desynchronises and reports every later cell as
wrong.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper" / "paper.tex"
DERIVED = ROOT / "data" / "derived"


def _table_body(label: str) -> list[list[str]]:
    """Cells of the tabular carrying `label`, one list per data row."""
    tex = PAPER.read_text()
    block = tex[tex.index(rf"\label{{{label}}}"):]
    block = block[:block.index(r"\end{tabular")]
    rows = []
    for line in block.split("\n"):
        if "&" not in line or line.lstrip().startswith("%"):
            continue
        if "textbf{" in line and "midrule" not in line:
            continue
        cells = [c.strip().replace("$", "").replace("\\", "")
                 for c in line.replace(r"\\", "").split("&")]
        rows.append(cells)
    return rows


def _fmt(x) -> str:
    return "---" if pd.isna(x) else f"{x:.2f}"


def _benchmark_key(anchor: str, ranking: str) -> tuple[str, str]:
    a, r = anchor.lower(), ranking.lower()
    if "doi" in a:
        anc = "doi"
    elif "nuvalent" in a or "industry" in a:
        src = a if "nuvalent" in a else r
        anc = "ind_" + ("alec" if "alec" in src else
                        "lor" if "lor" in src else "zot")
    elif "sot" in a:
        anc = "awad_sot"
    elif "ada" in a:
        anc = "awad_ada"
    else:
        raise AssertionError(f"unrecognised anchor: {anchor!r}")
    if "framework" in r or "log s +" in r:
        rk = "fw"
    elif "alone" in r and "log s" in r:
        rk = "supply"
    else:
        rk = "sigma"
    return anc, rk


def test_table2_matches_benchmark_csv():
    d = pd.read_csv(DERIVED / "benchmark_r3.csv")
    source = {_benchmark_key(r.anchor, r.ranking): r for r in d.itertuples()}

    checked = 0
    current_anchor = ""
    for cells in _table_body("tab:benchmark"):
        if len(cells) < 9:
            continue
        if cells[0]:
            current_anchor = cells[0]
        key = _benchmark_key(current_anchor, cells[1])
        assert key in source, f"Table 2 row {key} has no row in benchmark_r3.csv"
        r = source[key]
        expected = [_fmt(r._5), _fmt(r._6), _fmt(r._7), _fmt(r.AUROC),
                    _fmt(r.perm_p_at_100), _fmt(r.min_detectable_gain)]
        for name, got, want in zip(
                ["top10", "top50", "top100", "AUROC", "p", "MDG"],
                cells[3:9], expected):
            assert got == want, (
                f"Table 2 {key} {name}: manuscript has {got!r}, "
                f"benchmark_r3.csv has {want!r}")
            checked += 1
    assert checked == 96, f"expected 96 checked cells, checked {checked}"


def test_table4_matches_transport_json():
    j = json.loads((DERIVED / "sigma_transport.json").read_text())["transport"]
    rows = {c[0].strip().lower(): c for c in _table_body("tab:transport")
            if len(c) >= 5 and c[0].strip()}
    assert set(rows) == set(j), f"Table 4 drugs {set(rows)} vs {set(j)}"

    for drug, cells in rows.items():
        t = j[drug]
        assert int(cells[1]) == t["n"], f"{drug}: n"
        assert cells[2] == f'{t["spearman_rho"]:.2f}', f"{drug}: rho"
        lo, hi = t["boot95"]
        got = cells[3].replace("-", "-").replace(" to ", "|")
        want = f"{lo:.2f}|{hi:.2f}".replace("-", "-")
        assert got.replace("−", "-") == want, (
            f"{drug}: bootstrap interval, manuscript {cells[3]!r} vs "
            f"[{lo:.2f}, {hi:.2f}]")
        assert cells[4] == f'{t["median_inter_study_spread_log10"]:.2f}', (
            f"{drug}: spread")


@pytest.mark.parametrize("label", ["tab:benchmark", "tab:transport",
                                   "tab:enu-null", "tab:current-practice"])
def test_table_is_present_and_non_empty(label):
    assert _table_body(label), f"{label} parsed to no rows"


def test_supplement_track_matches_derived():
    """The supplementary track is a copy. Copies go stale.

    paper/supplement/reachability_track.json is what the interactive figure
    loads, and it is a byte copy of the pipeline's output. It sat at the old
    1,205-variant KRAS universe for eleven days after the classifier was
    corrected, because nothing compared the two. This does.
    """
    derived = DERIVED / "reachability_track.json"
    shipped = ROOT / "paper" / "supplement" / "reachability_track.json"
    assert shipped.exists(), f"missing {shipped}"
    a, b = json.loads(derived.read_text()), json.loads(shipped.read_text())
    for gene in a["by_gene"]:
        assert gene in b["by_gene"], f"supplement is missing {gene}"
        ga, gb = a["by_gene"][gene], b["by_gene"][gene]
        for field in ("n_class_forced_variants", "n_residues",
                      "residue_min", "residue_max"):
            assert ga[field] == gb[field], (
                f"{gene}.{field}: derived {ga[field]}, "
                f"supplement {gb[field]} (re-copy the file)")
    assert derived.read_bytes() == shipped.read_bytes(), (
        "supplement copy differs from the derived file; "
        "cp data/derived/reachability_track.json paper/supplement/")
