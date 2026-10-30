#!/usr/bin/env python3
"""Prospective audit of the Tiedt 2011 ENU MET screen against the Estevam atlas.

This is the third pipeline, and the only one where the chemistry-biased screen
reports its hits at the nucleotide level.  That changes what can be measured:
for KRAS and ALK the substitution class of a hit has to be inferred from the
amino-acid change, whereas Tiedt's Supplementary Tables S1 and S2 give the
actual base change and the number of clones carrying it.  The ENU class
spectrum is therefore observed rather than assumed.

Four analyses
-------------
A  Observed ENU class composition.  Clone-weighted six-class distribution of
   the base changes Tiedt reports, per screen, against the germline ENU
   reference and against the model-implied null of section 2.5.

B  The Y1230 single-codon experiment.  Under NVP-BVU972, 75 of 85 clones carry
   a substitution at one codon.  Every single-nucleotide route out of that
   codon is enumerable, so this is a controlled test of whether the mutagen's
   class weighting decides which variants appear or only how often.

C  Reachability of the chemistry-free priorities.  For each screen's drug, the
   top-K variants by atlas sigma, classified by how a one-step mutagen could
   reach them.

D  Recovery.  Where Tiedt's reported hits rank in the chemistry-free atlas.

Inputs
------
data/reference/codons_met.yaml           verified against NM_000245.2
data/raw/met_estevam2024/met_dms_wide.csv
data/reference/spectra.yaml

Outputs
-------
data/derived/met_audit.csv
data/derived/met_audit.json
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.codon_tools import sub_class, translate  # noqa: E402
from src.spectrum import load_spectra  # noqa: E402
from scripts.build_met_codons import TIEDT_S1, TIEDT_S2  # noqa: E402

CODONS = ROOT / "data" / "reference" / "codons_met.yaml"
DMS = ROOT / "data" / "raw" / "met_estevam2024" / "met_dms_wide.csv"
OUT_CSV = ROOT / "data" / "derived" / "met_audit.csv"
OUT_JSON = ROOT / "data" / "derived" / "met_audit.json"

CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
BASES = "ACGT"
BACKGROUND = "exon14del"      # primary; WT reported as sensitivity
SCREENS = {"NVP-BVU972": ("delta_NVP", TIEDT_S1),
           "AMG 458": ("delta_A458", TIEDT_S2)}
TOP_K = (10, 50, 100)
BETA = 1.0


def one_step_routes(codon: str, target_aa: str) -> list[tuple[str, str]]:
    """(base change label, pyrimidine class) for each single-nt route to target_aa."""
    out = []
    for i in range(3):
        for b in BASES:
            if b == codon[i]:
                continue
            mut = codon[:i] + b + codon[i + 1:]
            if translate(mut) == target_aa:
                out.append((f"{codon[i].lower()}{i}{b.lower()}", sub_class(codon[i], b)))
    return out


def classify(codon: str, target_aa: str) -> tuple[str, str | None, int]:
    """reachability category, forced class if any, number of one-step routes."""
    routes = one_step_routes(codon, target_aa)
    if not routes:
        return "two-step", None, 0
    classes = {c for _, c in routes}
    if len(classes) == 1:
        return "class-forced", classes.pop(), len(routes)
    return "multi-class SNV", None, len(routes)


def main() -> None:
    codons = {int(k): v for k, v in
              yaml.safe_load(CODONS.read_text())["codons"].items()}
    enu = dict(load_spectra(ROOT / "data" / "reference" / "spectra.yaml")
               ["ENU_Arnold2012_nominal"].values)
    dms_all = pd.read_csv(DMS)

    res: dict = {"background_primary": BACKGROUND, "beta": BETA,
                 "germline_enu": enu, "screens": {}, "y1230": {}, "sensitivity": {}}
    rows = []

    # ---------------------------------------------------------------- A, D --
    for screen, (sigma_col, table) in SCREENS.items():
        obs = defaultdict(int)          # class -> clones
        obs_var = defaultdict(int)      # class -> distinct base changes
        for residue, aa_change, pos, ref, alt, clones in table:
            cls = sub_class(ref.upper(), alt.upper())
            obs[cls] += clones
            obs_var[cls] += 1
            rows.append(dict(screen=screen, aa_change=aa_change, residue=residue,
                             cds_pos=pos, ref=ref, alt=alt, observed_class=cls,
                             clones=clones))
        n_clones = sum(obs.values())
        n_alt = sum(obs_var.values())

        # chi-square of clone-weighted composition against the germline spectrum
        exp = np.array([enu[c] for c in CLASSES]) * n_clones
        got = np.array([obs[c] for c in CLASSES], dtype=float)
        chi = float(((got - exp) ** 2 / np.maximum(exp, 1e-12)).sum())
        p_chi = float(stats.chi2.sf(chi, df=len(CLASSES) - 1))

        d = dms_all[dms_all.background == BACKGROUND].dropna(subset=[sigma_col]).copy()
        d["rank"] = d[sigma_col].rank(ascending=False, method="min").astype(int)
        hits = {aa for _, aa, _, _, _, _ in table}
        hit_ranks = (d[d.aa_change.isin(hits)][["aa_change", "rank"]]
                     .sort_values("rank").to_dict("records"))

        res["screens"][screen] = {
            "n_clones": n_clones,
            "n_distinct_base_changes": n_alt,
            "n_distinct_aa_changes": len({aa for _, aa, _, _, _, _ in table}),
            "observed_class_clones": {c: obs[c] for c in CLASSES},
            "observed_class_fraction": {c: obs[c] / n_clones for c in CLASSES},
            "germline_expected_clones": {c: enu[c] * n_clones for c in CLASSES},
            "chi2_vs_germline": chi, "p_chi2": p_chi,
            "classes_observed": sorted(c for c in CLASSES if obs[c] > 0),
            "classes_absent": sorted(c for c in CLASSES if obs[c] == 0),
            "atlas_n": int(len(d)),
            "hit_ranks": hit_ranks,
            "hits_found_in_atlas": len(hit_ranks),
            "hits_total": len(hits),
            "median_hit_rank": (int(np.median([h["rank"] for h in hit_ranks]))
                                if hit_ranks else None),
        }

        # ------------------------------------------------------------- C ----
        reach = defaultdict(int)
        for r in d.nlargest(max(TOP_K), sigma_col).itertuples():
            cat, _, _ = classify(codons[int(r.residue)], r.mutation)
            reach[cat] += 1
        res["screens"][screen]["top100_reachability"] = dict(reach)

    # ------------------------------------------------------------------ B --
    wt_codon = codons[1230]
    routes_by_aa: dict[str, list[tuple[str, str]]] = {}
    for i in range(3):
        for b in BASES:
            if b == wt_codon[i]:
                continue
            mut = wt_codon[:i] + b + wt_codon[i + 1:]
            aa = translate(mut)
            if aa in ("*", "Y"):
                continue
            routes_by_aa.setdefault(aa, []).append(
                (f"{wt_codon[i].lower()}{3*1230-2+i}{b.lower()}", sub_class(wt_codon[i], b)))

    found = {aa: clones for residue, aa_change, _, _, _, clones in TIEDT_S1
             if residue == 1230 for aa in [aa_change[-1]]}
    d = dms_all[dms_all.background == BACKGROUND].dropna(subset=["delta_NVP"]).copy()
    d["rank"] = d.delta_NVP.rank(ascending=False, method="min").astype(int)
    sig = d.set_index("aa_change")["rank"].to_dict()

    per_aa = {}
    for aa, rts in sorted(routes_by_aa.items()):
        supply = sum(enu[c] for _, c in rts)      # per-route supply summed
        per_aa[f"Y1230{aa}"] = {
            "routes": [{"change": ch, "class": c} for ch, c in rts],
            "summed_route_supply": supply,
            "clones_observed": found.get(aa, 0),
            "recovered": aa in found,
            "atlas_rank": sig.get(f"Y1230{aa}"),
        }

    n_y = sum(found.values())
    exp_c = {aa: v["summed_route_supply"] for aa, v in per_aa.items()}
    tot = sum(exp_c.values())
    exp_counts = {aa: n_y * v / tot for aa, v in exp_c.items()}
    obs_counts = {aa: v["clones_observed"] for aa, v in per_aa.items()}
    chi = float(sum((obs_counts[a] - exp_counts[a]) ** 2 / max(exp_counts[a], 1e-12)
                    for a in per_aa))
    res["y1230"] = {
        "wt_codon": wt_codon,
        "n_single_nt_missense_routes": len(per_aa),
        "n_recovered": sum(1 for v in per_aa.values() if v["recovered"]),
        "n_clones_at_codon": n_y,
        "per_variant": per_aa,
        "classes_spanned": sorted({c for v in per_aa.values()
                                   for c in [r["class"] for r in v["routes"]]}),
        "supply_range_of_those_classes": [
            min(enu[c] for v in per_aa.values() for c in [r["class"] for r in v["routes"]]),
            max(enu[c] for v in per_aa.values() for c in [r["class"] for r in v["routes"]])],
        "expected_clones_under_supply": exp_counts,
        "observed_clones": obs_counts,
        "chi2_within_codon": chi,
        "p_chi2_within_codon": float(stats.chi2.sf(chi, df=len(per_aa) - 1)),
        "missing_variant_atlas_rank": {
            k: v["atlas_rank"] for k, v in per_aa.items() if not v["recovered"]},
    }

    # --------------------------------------------------------- sensitivity --
    for screen, (sigma_col, table) in SCREENS.items():
        d = dms_all[dms_all.background == "WT"].dropna(subset=[sigma_col]).copy()
        d["rank"] = d[sigma_col].rank(ascending=False, method="min").astype(int)
        hits = {aa for _, aa, _, _, _, _ in table}
        hr = d[d.aa_change.isin(hits)]["rank"]
        res["sensitivity"][f"{screen}_WT_background"] = {
            "n_found": int(len(hr)), "n_total": len(hits),
            "median_rank": int(hr.median()) if len(hr) else None,
            "atlas_n": int(len(d)),
        }

    # --------------------------------------------------- homogeneity test --
    # The decisive control available only here: the SAME mutagen, cells and
    # protocol were used for both screens, so if the class composition of a
    # hit list were set by chemistry the two would agree. They do not.
    obs = np.array([[res["screens"][s]["observed_class_clones"][c] for c in CLASSES]
                    for s in SCREENS])
    keep = obs.sum(axis=0) > 0
    chi2, p_hom, dof, _ = stats.chi2_contingency(obs[:, keep])
    res["two_screen_homogeneity"] = {
        "screens": list(SCREENS),
        "table": {s: res["screens"][s]["observed_class_clones"] for s in SCREENS},
        "classes_tested": [c for c, k in zip(CLASSES, keep) if k],
        "chi2": float(chi2), "dof": int(dof), "p": float(p_hom),
        "interpretation": ("Same ENU protocol, same Ba/F3 TPR-MET cells, same "
                           "laboratory; the two drugs give different class "
                           "compositions, so composition is set by selection."),
    }

    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    OUT_JSON.write_text(json.dumps(res, indent=2, default=float))

    # ------------------------------------------------------------- report --
    for s, v in res["screens"].items():
        print(f"=== {s} ===")
        print(f"  {v['n_clones']} clones, {v['n_distinct_base_changes']} base changes, "
              f"{v['n_distinct_aa_changes']} aa changes")
        print("  observed class composition (clone-weighted) vs germline ENU:")
        for c in CLASSES:
            print(f"    {c}  obs {v['observed_class_clones'][c]:3d} "
                  f"({100*v['observed_class_fraction'][c]:5.1f}%)   "
                  f"germline expects {v['germline_expected_clones'][c]:5.1f} "
                  f"({100*enu[c]:4.1f}%)")
        print(f"  chi2 = {v['chi2_vs_germline']:.1f}, p = {v['p_chi2']:.3g}")
        print(f"  classes absent: {v['classes_absent']}")
        print(f"  hits found in atlas: {v['hits_found_in_atlas']}/{v['hits_total']}, "
              f"median rank {v['median_hit_rank']} of {v['atlas_n']}")
        print(f"  top-100 chemistry-free reachability: {v['top100_reachability']}")
        print()

    h = res["two_screen_homogeneity"]
    print(f"=== same mutagen, two drugs ===")
    print(f"  chi2 = {h['chi2']:.1f}, dof = {h['dof']}, p = {h['p']:.4f} "
          f"over classes {h['classes_tested']}")
    print()

    y = res["y1230"]
    print(f"=== Y1230 single-codon experiment (codon {y['wt_codon']}) ===")
    print(f"  {y['n_recovered']}/{y['n_single_nt_missense_routes']} single-nucleotide "
          f"missense routes recovered, {y['n_clones_at_codon']} clones")
    print(f"  classes spanned {y['classes_spanned']}, ENU supply "
          f"{y['supply_range_of_those_classes'][0]:.3f}-{y['supply_range_of_those_classes'][1]:.3f} "
          f"({y['supply_range_of_those_classes'][1]/y['supply_range_of_those_classes'][0]:.1f}-fold)")
    for k, v in y["per_variant"].items():
        cls = "/".join(r["class"] for r in v["routes"])
        print(f"    {k}  class {cls:8s} supply {v['summed_route_supply']:.3f}  "
              f"clones {v['clones_observed']:2d}  expected {y['expected_clones_under_supply'][k]:5.1f}  "
              f"atlas rank {v['atlas_rank']}")
    print(f"  within-codon chi2 = {y['chi2_within_codon']:.1f}, "
          f"p = {y['p_chi2_within_codon']:.3f}")
    print(f"  not recovered: {y['missing_variant_atlas_rank']}")


if __name__ == "__main__":
    main()
