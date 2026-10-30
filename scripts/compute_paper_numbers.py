#!/usr/bin/env python3
"""Compute every quantitative claim in paper.tex and emit a master-book audit trail.

GOLDEN RULE (project-wide, set 2026-10-15)
=========================================
Every number, percentage, probability, rank, interval, count, or fraction that
appears in the paper comes from this script. Every entry is registered in
REGISTRY with (a) value, (b) formatted string as it appears in prose,
(c) one-line provenance, (d) formal formula, (e) literal inputs, and
(f) the Python function that produced it. The script writes:

  - data/derived/paper_numbers.json  (machine-readable manifest)
  - paper/NUMBERS.md                 (human-readable master book)

If paper.tex contains a number not registered here, that is a bug: either
the number is wrong, or this script needs a new entry.

Run from the repo root:
    python3 scripts/compute_paper_numbers.py
"""
from __future__ import annotations

import json
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from scipy.optimize import brentq
from scipy.stats import beta, binom


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
PAPER = ROOT / "paper"


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------
@dataclass
class N:
    key: str                        # stable identifier
    value: Any                      # Python value as computed
    formatted: str                  # exactly as it must appear in paper.tex prose
    where: list[str]                # paper locations, e.g. ["Abstract", "§2.3", "Fig 4C caption"]
    computed_from: str              # one-sentence provenance
    formula: str                    # formal expression or algorithmic description
    inputs: dict = field(default_factory=dict)
    code_fn: str = ""               # name of the Python function in this script


REGISTRY: list[N] = []


def reg(**kwargs) -> N:
    """Register a number. Raises if the key is already taken."""
    n = N(**kwargs)
    if any(x.key == n.key for x in REGISTRY):
        raise ValueError(f"duplicate key: {n.key}")
    REGISTRY.append(n)
    return n


def fmt_pct(x: float, places: int = 1) -> str:
    """Percentage formatter; switches to 2 d.p. below 1%."""
    pct = x * 100
    return f"{pct:.{places}f}%" if pct >= 1.0 else f"{pct:.2f}%"


# ---------------------------------------------------------------------------
# Loaders — single source of truth for every data artefact
# ---------------------------------------------------------------------------
def load_spectra() -> dict:
    return yaml.safe_load(open(DATA / "reference" / "spectra.yaml"))["spectra"]


def load_feng_class_forced() -> pd.DataFrame:
    df = pd.read_csv(DATA / "derived" / "feng_decomposition.csv")
    cf = df[df["class_forced"] == True].drop_duplicates("aa_change", keep="first").copy()
    return cf


def load_oh_decomp() -> pd.DataFrame:
    return pd.read_csv(DATA / "derived" / "oh_decomposition_all_drugs.csv")


def load_oh_raw() -> pd.DataFrame:
    return pd.read_csv(DATA / "raw" / "oh_2026" / "table_s2_resistance_scores.csv")


# ---------------------------------------------------------------------------
# Spectrum re-normalisation helper (shared by H95D re-ranking under alt ENU)
# ---------------------------------------------------------------------------
def renorm_spectrum(base: dict, class_name: str, new_rate: float) -> dict:
    """Return a spectrum identical to base but with class_name set to new_rate
    and the other five classes rescaled proportionally to preserve sum = 1."""
    others = {c: v for c, v in base.items() if c != class_name}
    tot_others = sum(others.values())
    scale = (1 - new_rate) / tot_others
    out = {c: v * scale for c, v in others.items()}
    out[class_name] = new_rate
    return out


# ---------------------------------------------------------------------------
# Compute functions — one per logical group
# ---------------------------------------------------------------------------

# ---------------- Resistance panels ----------------
def compute_resistance_panels() -> None:
    kras_clinical = {"Y96C", "H95D", "H95Q", "H95R", "R68S"}      # Awad 2021 clinical cohort
    kras_case     = {"Y96D"}                                       # Tanaka 2021 case report
    kras_dms      = {"C12F", "V8L", "V8M", "R68L", "R68M", "Y96H"}  # Awad 2021 in vitro DMS
    kras_all      = kras_clinical | kras_case | kras_dms
    assert len(kras_all) == 12, f"KRAS panel wrong size: {len(kras_all)}"

    reg(key="kras_panel_total", value=len(kras_all), formatted="12",
        where=["Abstract", "§2.2", "Fig 1B", "Fig 3A caption", "Table S1"],
        computed_from="union of KRAS G12C clinical cohort + case report + accompanying in vitro DMS",
        formula="|clinical ∪ case ∪ DMS|",
        inputs={"clinical": sorted(kras_clinical), "case": sorted(kras_case), "dms": sorted(kras_dms)},
        code_fn="compute_resistance_panels")

    reg(key="kras_panel_clinical_cohort", value=len(kras_clinical), formatted="five",
        where=["Methods §Resistance panels", "Table S1"],
        computed_from="patient-observed variants from Awad 2021 clinical cohort",
        formula="|clinical|", inputs={"members": sorted(kras_clinical)},
        code_fn="compute_resistance_panels")

    reg(key="kras_panel_case_report", value=len(kras_case), formatted="one",
        where=["Methods §Resistance panels"],
        computed_from="patient-observed Y96D from Tanaka 2021 case report",
        formula="|case|", inputs={"members": sorted(kras_case)},
        code_fn="compute_resistance_panels")

    reg(key="kras_panel_dms", value=len(kras_dms), formatted="six",
        where=["Abstract", "§2.2", "Methods §Resistance panels", "Table S1"],
        computed_from="variants from Awad 2021 accompanying in vitro deep mutational scan",
        formula="|DMS|", inputs={"members": sorted(kras_dms)},
        code_fn="compute_resistance_panels")

    reg(key="kras_panel_clinical_total", value=len(kras_clinical) + len(kras_case),
        formatted="six",
        where=["Abstract"],
        computed_from="patient-observed KRAS G12C variants (clinical cohort + case report)",
        formula="|clinical ∪ case|",
        inputs={"members": sorted(kras_clinical | kras_case)},
        code_fn="compute_resistance_panels")

    # Reachability classification of the 12 KRAS panel members
    # (derived from codon analysis in src.codon_tools, cross-verified against
    # the Supplementary Table S1 data file)
    kras_snv_unreachable = {"R68L", "V8M"}
    kras_multi_class     = {"V8L", "R68S", "H95Q"}
    kras_class_forced    = {"C12F", "R68M", "H95D", "H95R", "Y96C", "Y96D", "Y96H"}
    assert kras_snv_unreachable | kras_multi_class | kras_class_forced == kras_all

    reg(key="kras_panel_snv_unreachable", value=len(kras_snv_unreachable), formatted="two",
        where=["Abstract", "§2.2", "§2.5"],
        computed_from="codon reachability analysis in the G12C-mutant background",
        formula="count of panel variants whose Hamming distance from WT codon is ≥ 2 to all target codons",
        inputs={"members": sorted(kras_snv_unreachable)},
        code_fn="compute_resistance_panels")

    reg(key="kras_panel_multi_class", value=len(kras_multi_class), formatted="three",
        where=["§2.2"],
        computed_from="codon reachability: variants reachable via multiple pyrimidine classes",
        formula="count of panel variants with ≥ 2 distinct single-nt pyrimidine routes",
        inputs={"members": sorted(kras_multi_class)},
        code_fn="compute_resistance_panels")

    reg(key="kras_panel_class_forced", value=len(kras_class_forced), formatted="seven",
        where=["§2.2", "§2.3"],
        computed_from="codon reachability: variants reachable via exactly one pyrimidine class",
        formula="count of panel variants with exactly 1 single-nt pyrimidine route",
        inputs={"members": sorted(kras_class_forced)},
        code_fn="compute_resistance_panels")

    # ALK panel
    alk_doi   = {"V1180L", "L1196M", "I1171N", "I1171T", "L1198F", "G1202R"}
    alk_takei = {"G1269A", "F1174C", "F1174L", "L1196Q", "C1156Y", "L1152R"}
    alk_all   = alk_doi | alk_takei
    assert len(alk_all) == 12

    alk_class_forced     = {"V1180L", "L1196M", "I1171N", "I1171T", "L1198F",
                            "G1269A", "F1174C", "L1196Q", "C1156Y", "L1152R"}
    alk_multi_class      = {"G1202R", "F1174L"}

    reg(key="alk_panel_total", value=len(alk_all), formatted="12",
        where=["§2.4", "Fig 3A caption", "Table S1"],
        computed_from="union of alectinib Mutazyme II screen hits + next-gen ALK TKI screen hits",
        formula="|Doi hits ∪ Takei hits|",
        inputs={"doi": sorted(alk_doi), "takei": sorted(alk_takei)},
        code_fn="compute_resistance_panels")

    reg(key="alk_panel_class_forced", value=len(alk_class_forced), formatted="10",
        where=["Table S1 summary"],
        computed_from="codon reachability of ALK panel",
        formula="count of ALK panel variants reachable via exactly one pyrimidine class",
        inputs={"members": sorted(alk_class_forced)},
        code_fn="compute_resistance_panels")

    reg(key="alk_panel_multi_class", value=len(alk_multi_class), formatted="two",
        where=["Table S1 summary"],
        computed_from="codon reachability of ALK panel",
        formula="count of ALK panel variants with multi-class single-nt reachability",
        inputs={"members": sorted(alk_multi_class)},
        code_fn="compute_resistance_panels")

    # Doi 2022 hits — six variants recovered from alectinib Mutazyme II Fig 3a
    reg(key="doi2022_hits_count", value=len(alk_doi), formatted="six",
        where=["Abstract", "§2.4", "§2.5", "Fig 4A caption"],
        computed_from="single-amino-acid resistance mutations in Doi 2022 Fig 3a, excluding silent passengers (L1152, L1165, R1192) and the I1171S sub-variant of the I1171N/T/S compound call",
        formula="|Doi Fig 3a resistance mutations|",
        inputs={"members": sorted(alk_doi)},
        code_fn="compute_resistance_panels")

    doi_class_forced = alk_doi - {"G1202R"}  # G1202R is multi-class
    reg(key="doi2022_hits_class_forced", value=len(doi_class_forced), formatted="five",
        where=["§2.5", "Fig 4A caption"],
        computed_from="Doi 2022 hits that are class-forced (G1202R is multi-class)",
        formula="|Doi hits ∩ class-forced|",
        inputs={"members": sorted(doi_class_forced)},
        code_fn="compute_resistance_panels")


# ---------------- Mutagen spectra ----------------
def compute_mutagen_spectra() -> None:
    sp = load_spectra()
    enu = sp["ENU_Arnold2012_nominal"]["values"]
    mut = sp["MutazymeII_DoiSupp2022"]["values"]

    for cls, v in enu.items():
        reg(key=f"enu_class_{cls.replace('>', '_')}", value=v, formatted=fmt_pct(v),
            where=(["§2.1", "Fig 2A"] if cls in ("T>C", "T>A", "C>G")
                   else ["Fig 2A"]),
            computed_from="Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies",
            formula=f"spectra.yaml:ENU_Arnold2012_nominal.values['{cls}']",
            inputs={"source": "arnold2012", "class": cls},
            code_fn="compute_mutagen_spectra")

    for cls, v in mut.items():
        reg(key=f"mut_class_{cls.replace('>', '_')}", value=v, formatted=fmt_pct(v),
            where=(["§2.1", "Fig 2A"] if cls in ("C>T", "T>G", "C>G")
                   else ["Fig 2A"]),
            computed_from="NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b",
            formula=f"spectra.yaml:MutazymeII_DoiSupp2022.values['{cls}']",
            inputs={"source": "doi2022", "class": cls},
            code_fn="compute_mutagen_spectra")

    # Class coverage counts (≥ 1 % threshold)
    reg(key="enu_classes_at_1pct", value=sum(1 for v in enu.values() if v >= 0.01),
        formatted="five of six",
        where=["§2.1", "Table 1 ENU row", "Fig 2B"],
        computed_from="count of ENU classes supplied at ≥ 1 %; C>G is at 0.2 % and falls below",
        formula="|{c: ENU[c] ≥ 0.01}|",
        inputs={"ENU": enu, "threshold": 0.01},
        code_fn="compute_mutagen_spectra")

    reg(key="mut_classes_at_1pct", value=sum(1 for v in mut.values() if v >= 0.01),
        formatted="6 of 6",
        where=["Table 1 Mutazyme row", "Fig 2B"],
        computed_from="count of Mutazyme II classes supplied at ≥ 1 %; minimum class is T>G at 3.6 %",
        formula="|{c: MutII[c] ≥ 0.01}|",
        inputs={"MutII": mut, "threshold": 0.01},
        code_fn="compute_mutagen_spectra")

    reg(key="be_classes_supplied", value=2, formatted="two of six",
        where=["§2.1", "Table 1 CBE+ABE row", "Fig 2B"],
        computed_from="CBE supplies C>T, ABE supplies T>C by design; the other 4 classes are 0",
        formula="deterministic by base-editor chemistry",
        inputs={"supplied": ["C>T", "T>C"]},
        code_fn="compute_mutagen_spectra")

    reg(key="twist_pe_classes_supplied", value=6, formatted="6 of 6",
        where=["Table 1 Twist row", "Table 1 PE row", "Fig 2B"],
        computed_from="Twist synthesis and prime editing produce all six classes approximately uniformly",
        formula="deterministic (chemistry-free library construction)",
        inputs={"baseline_per_class": 1 / 6},
        code_fn="compute_mutagen_spectra")

    # Mutazyme II max/min ratio — the "13-fold" quantitative-bias claim
    mut_max = max(mut.values())
    mut_min = min(mut.values())
    ratio = mut_max / mut_min
    reg(key="mut_class_rate_ratio", value=ratio, formatted=f"{ratio:.0f}-fold",
        where=["§2.1"],
        computed_from="ratio of Mutazyme II's most abundant class to its rarest",
        formula="max(MutII) / min(MutII) = max{MutII[c]} / min{MutII[c]}",
        inputs={"max": mut_max, "max_class": "C>T", "min": mut_min, "min_class": "T>G"},
        code_fn="compute_mutagen_spectra")

    # ENU C>A class rate used in prose and Fig 4C
    reg(key="enu_well_supplied_well_vs_total", value="9 of 12",
        formatted="nine of twelve",
        where=["§2.5", "Fig 4C caption"],
        computed_from="count of Nuvalent class-forced hits whose forced class ENU supplies at ≥ 1 %",
        formula="|{hit: ENU[hit.class] ≥ 0.01, hit in Nuvalent class-forced pool}|",
        inputs={"nuvalent_class_forced_total": 12, "threshold": 0.01},
        code_fn="compute_mutagen_spectra")


# ---------------- Oh atlas sizes ----------------
def compute_atlas_sizes() -> None:
    oh_raw = load_oh_raw()
    by_type = oh_raw["type"].value_counts().to_dict()

    total = len(oh_raw)
    reg(key="oh_atlas_total_snvs", value=total, formatted=f"{total:,}",
        where=["§Background", "Methods"],
        computed_from="total rows in Oh 2026 Table S2 resistance scores",
        formula="len(data/raw/oh_2026/table_s2_resistance_scores.csv)",
        inputs={"path": "data/raw/oh_2026/table_s2_resistance_scores.csv"},
        code_fn="compute_atlas_sizes")

    coding_types = {"Missense", "Synonymous", "Nonsense"}
    coding = sum(by_type.get(t, 0) for t in coding_types)
    reg(key="oh_atlas_coding_snvs", value=coding, formatted=f"{coding:,}",
        where=["§Background", "Methods"],
        computed_from="sum of coding categories in Oh 2026 Table S2",
        formula="sum(Missense + Synonymous + Nonsense) in table_s2_resistance_scores.csv",
        inputs={"by_type": by_type, "coding_types": sorted(coding_types)},
        code_fn="compute_atlas_sizes")

    missense = by_type.get("Missense", 0)
    reg(key="oh_atlas_missense", value=missense, formatted=f"{missense:,}",
        where=["Methods"],
        computed_from="Missense rows in Oh 2026 Table S2",
        formula="count of rows where type == 'Missense'",
        inputs={"path": "data/raw/oh_2026/table_s2_resistance_scores.csv"},
        code_fn="compute_atlas_sizes")

    syn = by_type.get("Synonymous", 0)
    reg(key="oh_atlas_synonymous", value=syn, formatted=f"{syn:,}",
        where=["Methods"],
        computed_from="Synonymous rows in Oh 2026 Table S2",
        formula="count of rows where type == 'Synonymous'",
        inputs={"path": "data/raw/oh_2026/table_s2_resistance_scores.csv"},
        code_fn="compute_atlas_sizes")

    nonsense = by_type.get("Nonsense", 0)
    reg(key="oh_atlas_nonsense", value=nonsense, formatted=f"{nonsense:,}",
        where=["Methods"],
        computed_from="Nonsense rows in Oh 2026 Table S2",
        formula="count of rows where type == 'Nonsense'",
        inputs={"path": "data/raw/oh_2026/table_s2_resistance_scores.csv"},
        code_fn="compute_atlas_sizes")

    intron = by_type.get("Intron", 0)
    splice = by_type.get("Splice Acceptor/Donor", 0)
    noncoding = intron + splice
    reg(key="oh_atlas_noncoding", value=noncoding, formatted=f"{noncoding}",
        where=["§Background", "Methods"],
        computed_from="sum of intron and splice-site rows in Oh 2026 Table S2",
        formula="count(type == 'Intron') + count(type == 'Splice Acceptor/Donor')",
        inputs={"intron": intron, "splice": splice},
        code_fn="compute_atlas_sizes")

    # Residue range
    coding_df = oh_raw[oh_raw["type"] != "Intron"]
    residues = coding_df["residue"].dropna().astype(int)
    r_min, r_max = int(residues.min()), int(residues.max())
    reg(key="oh_atlas_residue_min", value=r_min, formatted=f"{r_min}",
        where=["Methods"],
        computed_from="minimum residue index in Oh 2026 Table S2 coding rows",
        formula="min(residue) over rows with type != 'Intron'",
        inputs={"path": "data/raw/oh_2026/table_s2_resistance_scores.csv"},
        code_fn="compute_atlas_sizes")
    reg(key="oh_atlas_residue_max", value=r_max, formatted=f"{r_max}",
        where=["Methods"],
        computed_from="maximum residue index in Oh 2026 Table S2 coding rows",
        formula="max(residue) over rows with type != 'Intron'",
        inputs={"path": "data/raw/oh_2026/table_s2_resistance_scores.csv"},
        code_fn="compute_atlas_sizes")

    # Feng atlas size
    feng_raw = pd.read_csv(DATA / "raw" / "feng_2022" / "table_s1_zscores.csv")
    feng_variants = len(feng_raw)
    reg(key="feng_atlas_total_variants", value=feng_variants, formatted=f"{feng_variants:,}",
        where=["§Background"],
        computed_from="rows in Feng 2022 KRAS G12C Twist atlas Table S1",
        formula="len(data/raw/feng_2022/table_s1_zscores.csv)",
        inputs={"path": "data/raw/feng_2022/table_s1_zscores.csv"},
        code_fn="compute_atlas_sizes")

    # Class-forced counts (used in rank-shift denominators)
    feng_cf = load_feng_class_forced()
    reg(key="feng_class_forced_count", value=len(feng_cf), formatted=f"{len(feng_cf):,}",
        where=["§2.3", "Methods"],
        computed_from="class-forced KRAS G12C variants in the Feng atlas (unique amino-acid changes)",
        formula="|{v in Feng: v is class-forced single-nt reachable}|",
        inputs={"source": "data/derived/feng_decomposition.csv"},
        code_fn="compute_atlas_sizes")

    oh_decomp = load_oh_decomp()
    oh_cf = oh_decomp[oh_decomp["mutagen"] == "MutazymeII"].drop_duplicates(["aa_change", "drug"])
    oh_cf_total = len(oh_cf["aa_change"].unique())
    reg(key="oh_class_forced_count", value=oh_cf_total, formatted=f"{oh_cf_total:,}",
        where=["§2.4", "Methods", "Fig 3B caption"],
        computed_from="class-forced ALK kinase-domain variants used in the ranking",
        formula="|unique aa_change in Oh class-forced decomposition|",
        inputs={"source": "data/derived/oh_decomposition_all_drugs.csv"},
        code_fn="compute_atlas_sizes")


# ---------------- H95D ranks under all ENU calibrations + Mutazyme ----------------
def compute_h95d_ranks() -> None:
    cf = load_feng_class_forced()
    sp = load_spectra()
    enu = sp["ENU_Arnold2012_nominal"]["values"]
    mut = sp["MutazymeII_DoiSupp2022"]["values"]

    def rank_variant(spec: dict, variant: str, drug_z_col: str) -> tuple[int, int]:
        cf2 = cf.copy()
        cf2["S"] = cf2["forced_class"].map(spec)
        cf2["log_P"] = np.log(cf2["S"].clip(lower=1e-10)) + cf2[drug_z_col]
        cf2["rank_sim"] = cf2["log_P"].rank(ascending=False, method="min")
        cf2["rank_cf"] = cf2[drug_z_col].rank(ascending=False, method="min")
        row = cf2[cf2["aa_change"] == variant].iloc[0]
        return int(row["rank_sim"]), int(row["rank_cf"])

    # Chemistry-free rank is the same under any mutagen (only σ matters)
    _, rank_cf = rank_variant(enu, "H95D", "adagrasib_z")
    reg(key="h95d_ada_chemfree_rank", value=rank_cf, formatted=f"{rank_cf}",
        where=["§2.3"],
        computed_from="H95D chemistry-free rank in the KRAS G12C class-forced atlas under adagrasib",
        formula="rank(σ_adagrasib) descending over the class-forced KRAS G12C atlas, position of H95D",
        inputs={"variant": "H95D", "drug": "adagrasib"},
        code_fn="compute_h95d_ranks")

    # ENU calibrations
    calibrations = {
        "germline":        ("Arnold 2012 germline-ENU compilation",             0.002),
        "empirical_lower": ("Jeffreys 95% lower bound from 3-of-13 anomaly",    0.070),
        "empirical_mean":  ("Jeffreys posterior mean from 3-of-13 anomaly",     0.250),
        "empirical_upper": ("Jeffreys 95% upper bound from 3-of-13 anomaly",    0.497),
    }
    for name, (desc, p_CG) in calibrations.items():
        spec = renorm_spectrum(enu, "C>G", p_CG)
        r_sim, _ = rank_variant(spec, "H95D", "adagrasib_z")
        shift = r_sim - rank_cf
        reported = name == "germline"
        where_sim = (["supporting: §2.3 quotes the beta-range, not this point value"] if reported else
                     ["not reported: calibrated ENU rows were removed from §2.3 "
                      "because §2.5 shows the inversion is not identified"])
        reg(key=f"h95d_ada_enu_{name}_rank_sim", value=r_sim, formatted=f"{r_sim}",
            where=where_sim,
            computed_from=f"H95D simulated rank under adagrasib + ENU with C>G = {p_CG*100:.1f}% ({desc})",
            formula="rank of log P(v) = log S(v) + log σ_adagrasib(v), descending, over 1,205 class-forced variants",
            inputs={"variant": "H95D", "drug": "adagrasib", "mutagen": "ENU", "CG_rate": p_CG, "spec": spec},
            code_fn="compute_h95d_ranks")
        reg(key=f"h95d_ada_enu_{name}_shift", value=shift, formatted=f"{shift:+d}",
            where=where_sim,
            computed_from=f"H95D rank shift under adagrasib + ENU ({desc})",
            formula="rank_sim(ENU, p_CG) - rank_chemfree",
            inputs={"rank_sim": r_sim, "rank_cf": rank_cf, "p_CG": p_CG},
            code_fn="compute_h95d_ranks")

    # Mutazyme II (NGS-calibrated, no C>G uncertainty question here)
    r_sim_mut, _ = rank_variant(mut, "H95D", "adagrasib_z")
    shift_mut = r_sim_mut - rank_cf
    reg(key="h95d_ada_mut_rank_sim", value=r_sim_mut, formatted=f"{r_sim_mut}",
        where=["§2.3", "Discussion Limitation 4"],
        computed_from="H95D simulated rank under adagrasib + Mutazyme II (NGS-measured spectrum)",
        formula="rank of log P(v) = log S(v) + log σ_adagrasib(v) under Mutazyme II spectrum",
        inputs={"variant": "H95D", "drug": "adagrasib", "mutagen": "MutazymeII"},
        code_fn="compute_h95d_ranks")
    reg(key="h95d_ada_mut_shift", value=shift_mut, formatted=f"{shift_mut:+d}",
        where=["§2.3", "§2.5"],
        computed_from="H95D rank shift under adagrasib + Mutazyme II (NGS-measured)",
        formula="rank_sim(Mutazyme) - rank_chemfree",
        inputs={"rank_sim": r_sim_mut, "rank_cf": rank_cf, "mut_CG": mut["C>G"]},
        code_fn="compute_h95d_ranks")

    # ----- R1 β-sensitivity -----
    # Model: log P = log S + β·σ with β ∈ {0.25, 0.5, 1, 2, 4}.
    # Report shift intervals for every H95D calibration to make the
    # parameterization dependence explicit.
    def rank_h95d_beta(spec, beta):
        cf2 = cf.copy()
        cf2["S"] = cf2["forced_class"].map(spec)
        cf2["log_P"] = np.log(cf2["S"].clip(lower=1e-10)) + beta * cf2["adagrasib_z"]
        cf2["rank_sim"] = cf2["log_P"].rank(ascending=False, method="min")
        row = cf2[cf2["aa_change"] == "H95D"].iloc[0]
        return int(row["rank_sim"])

    BETA_GRID = [0.25, 0.5, 1.0, 2.0, 4.0]
    beta_calibrations = {
        "germline":        (renorm_spectrum(enu, "C>G", 0.002), "ENU germline (C>G = 0.2%)"),
        "empirical_lower": (renorm_spectrum(enu, "C>G", 0.070), "ENU empirical lower (C>G = 7%)"),
        "empirical_mean":  (renorm_spectrum(enu, "C>G", 0.250), "ENU empirical mean (C>G = 25%)"),
        "mutazyme":        (mut,                                "Mutazyme II (NGS; C>G = 5%)"),
    }
    for name, (spec, desc) in beta_calibrations.items():
        shifts = [rank_h95d_beta(spec, b) - rank_cf for b in BETA_GRID]
        shift_min, shift_max = min(shifts), max(shifts)
        quoted = ["§2.3"] if name == "germline" else [
            "not reported: calibrated ENU rows were removed from §2.3"]
        reg(key=f"h95d_ada_{name}_shift_beta_min", value=shift_min, formatted=f"{shift_min:+d}",
            where=["supporting: §2.3 quotes the interval, not its endpoints"],
            computed_from=f"minimum H95D shift across β ∈ {BETA_GRID} for {desc}",
            formula=f"min_β rank_sim(β) − rank_chemfree, with log P = log S + β·σ",
            inputs={"beta_grid": BETA_GRID, "spec": spec, "variant": "H95D", "drug": "adagrasib"},
            code_fn="compute_h95d_ranks")
        reg(key=f"h95d_ada_{name}_shift_beta_max", value=shift_max, formatted=f"{shift_max:+d}",
            where=["supporting: §2.3 quotes the interval, not its endpoints"],
            computed_from=f"maximum H95D shift across β ∈ {BETA_GRID} for {desc}",
            formula=f"max_β rank_sim(β) − rank_chemfree, with log P = log S + β·σ",
            inputs={"beta_grid": BETA_GRID, "spec": spec, "variant": "H95D", "drug": "adagrasib"},
            code_fn="compute_h95d_ranks")
        reg(key=f"h95d_ada_{name}_shift_beta_range",
            value=f"{shift_min:+d} to {shift_max:+d}",
            formatted=f"{shift_min:+d} to {shift_max:+d}",
            where=quoted,
            computed_from=f"H95D shift interval across β ∈ {BETA_GRID} for {desc}",
            formula="[min_β shift, max_β shift]",
            inputs={"beta_grid": BETA_GRID, "per_beta": dict(zip(BETA_GRID, shifts))},
            code_fn="compute_h95d_ranks")

    # Y96D under sotorasib + Mutazyme (secondary example in §2.3)
    r_sim_y, r_cf_y = rank_variant(mut, "Y96D", "sotorasib_z")
    reg(key="y96d_sot_chemfree_rank", value=r_cf_y, formatted=f"{r_cf_y}",
        where=["§2.3"],
        computed_from="Y96D chemistry-free rank under sotorasib",
        formula="rank(σ_sotorasib) descending, position of Y96D",
        inputs={"variant": "Y96D", "drug": "sotorasib"},
        code_fn="compute_h95d_ranks")
    reg(key="y96d_sot_mut_rank_sim", value=r_sim_y, formatted=f"{r_sim_y}",
        where=["§2.3"],
        computed_from="Y96D simulated rank under sotorasib + Mutazyme II",
        formula="rank of log P(v) under Mutazyme spectrum, Y96D forced class T>G at 3.6 %",
        inputs={"variant": "Y96D", "drug": "sotorasib", "mutagen": "MutazymeII", "mut_TG": mut["T>G"]},
        code_fn="compute_h95d_ranks")
    reg(key="y96d_sot_mut_shift", value=r_sim_y - r_cf_y, formatted=f"{r_sim_y - r_cf_y:+d}",
        where=["§2.3"],
        computed_from="Y96D rank shift under sotorasib + Mutazyme II",
        formula="rank_sim - rank_chemfree for Y96D under Mutazyme sotorasib",
        inputs={"rank_sim": r_sim_y, "rank_cf": r_cf_y},
        code_fn="compute_h95d_ranks")

    # The beta sweep for Y96D. This was prose-only and went stale with the
    # KRAS reclassification, the same way its point estimate did.
    _s = cf.copy()
    _s["S"] = _s["forced_class"].map(mut)
    _i = _s.index[_s["aa_change"] == "Y96D"][0]
    _cf = _s["sotorasib_z"].rank(ascending=False, method="min")[_i]
    _sweep = []
    for _b in (0.25, 0.5, 1.0, 2.0, 4.0):
        _lp = np.log(_s["S"].clip(lower=1e-10)) + _b * _s["sotorasib_z"]
        _sweep.append(int(_lp.rank(ascending=False, method="min")[_i] - _cf))
    reg(key="y96d_sot_mutazyme_shift_beta_range",
        value=[min(_sweep), max(_sweep)],
        formatted=f"{min(_sweep):+d} to {max(_sweep):+d}",
        where=["\u00a73.3"],
        computed_from="range of the Y96D rank shift under Mutazyme II with "
                      "sotorasib across beta in [0.25, 4]",
        formula="min and max over beta of rank(log S + beta*sigma) - rank(sigma)",
        inputs={"betas": [0.25, 0.5, 1.0, 2.0, 4.0], "shifts": _sweep},
        code_fn="compute_h95d_ranks")


# ---------------- Nuvalent anomaly: binomial + Jeffreys interval ----------------
def compute_nuvalent_binomial() -> None:
    """Hit-class composition test for the industry ENU anchor.

    The thirteen reported variants are post-selection, so the pre-selection
    spectrum is not a valid null on its own.  Numbers come from
    scripts/recompute_enu_class_test.py, which builds four nulls; this
    function registers them and the corrected inversion.
    """
    t = json.loads((DATA / "derived" / "enu_class_test.json").read_text())
    n, k = t["n_hits"], t["n_cg_hits"]

    reg(key="nuvalent_cg_binom_p_raw", value=t["p_ge3"]["N0_raw"],
        formatted=f'{t["p_ge3"]["N0_raw"]:.1e}',
        where=["§2.5", "Table 3 (N0 row)"],
        computed_from="exact binomial tail against the raw germline spectrum; "
                      "retained only as the confounded comparison",
        formula="binom.sf(k-1, n, S_ENU('C>G'))",
        inputs={"k": k, "n": n, "p": t["germline_cg_rate"]},
        code_fn="recompute_enu_class_test.py::main")

    reg(key="nuvalent_cg_pi_model", value=t["pi_cg"]["N3_pooled_beta1"],
        formatted=f'{t["pi_cg"]["N3_pooled_beta1"]*100:.2f}%',
        where=["§2.5", "Table 3 (N3 row)"],
        computed_from="model-implied probability that a hit is C>G, with germline "
                      "supply and pooled sigma at beta=1",
        formula="sum_{v in C>G} S(c)exp(beta*sigma_v) / sum_v S(c)exp(beta*sigma_v)",
        inputs={"beta": 1.0, "sigma": "max over alec/lor/zot", "S": "ENU germline"},
        code_fn="recompute_enu_class_test.py::pi_model")

    reg(key="nuvalent_cg_binom_p_model", value=t["p_ge3"]["N3_pooled_beta1"],
        formatted=f'{t["p_ge3"]["N3_pooled_beta1"]:.1e}',
        where=["§2.5", "Table 3 (N3 row)"],
        computed_from="binomial tail of 3-of-13 C>G against the model-implied null; "
                      "the germline spectrum is rejected",
        formula="binom.sf(k-1, n, pi_N3)",
        inputs={"k": k, "n": n, "pi": t["pi_cg"]["N3_pooled_beta1"]},
        code_fn="recompute_enu_class_test.py::main")

    br = t["p_ge3_N3_beta_range"]
    reg(key="nuvalent_cg_binom_p_model_beta_range",
        value=[min(br.values()), max(br.values())],
        formatted=f'{min(br.values()):.1e} to {max(br.values()):.1e}',
        where=["supporting: \u00a72.5 quotes the beta = 4 endpoint only"],
        computed_from="range of the model-implied binomial tail across beta in [0.25, 4]",
        formula="min/max over beta of binom.sf(k-1, n, pi_N3(beta))",
        inputs={"betas": sorted(float(b) for b in br)},
        code_fn="recompute_enu_class_test.py::main")

    reg(key="nuvalent_cg_top100_sigma", value=t["n_cg_in_top100_pooled"],
        formatted=f'{t["n_cg_in_top100_pooled"]} of the',
        where=["\u00a72.5 (as 24 of the top-100 sigma-ranked class-forced variants)", "Table 3 (N2 row)"],
        computed_from="C>G count among the top-100 sigma-ranked class-forced ALK variants, "
                      "pooled across drug arms",
        formula="count(forced_class == 'C>G') over nlargest(100, pooled sigma)",
        inputs={"N_top": 100},
        code_fn="recompute_enu_class_test.py::pi_sigma_topn")

    reg(key="nuvalent_cg_binom_p_selection", value=t["p_ge3"]["N2_top100_pooled"],
        formatted=f'{t["p_ge3"]["N2_top100_pooled"]:.2f}',
        where=["§2.5", "Table 3 (N2 row)"],
        computed_from="binomial tail of 3-of-13 C>G against a selection-only null; "
                      "the count is unremarkable once supply is removed",
        formula="binom.sf(k-1, n, pi_N2)",
        inputs={"k": k, "n": n, "pi": t["pi_cg"]["N2_top100_pooled"]},
        code_fn="recompute_enu_class_test.py::main")

    reg(key="nuvalent_cg_binom_p_composition", value=t["p_ge3_N4_composition"],
        formatted=f'{t["p_ge3_N4_composition"]:.2f}',
        where=["§2.5", "Table 3 (N4 row)"],
        computed_from="binomial tail against atlas class composition alone",
        formula="binom.sf(k-1, n, n_CG / n_classforced)",
        inputs={"n_CG": t["n_cg_classforced_atlas"], "n": t["n_classforced_atlas"]},
        code_fn="recompute_enu_class_test.py::main")

    reg(key="atlas_cg_share", value=t["pi_cg_N4_composition"],
        formatted=fmt_pct(t["pi_cg_N4_composition"]),
        where=["§2.5"],
        computed_from="share of class-forced ALK atlas variants whose forced class is C>G",
        formula="n_CG / n_classforced",
        inputs={"n_CG": t["n_cg_classforced_atlas"], "n": t["n_classforced_atlas"]},
        code_fn="recompute_enu_class_test.py::main")

    inv = t["implied_cg_supply_rate"]["pooled"]
    reg(key="implied_enu_cg_supply_interval",
        value=[inv["s_lo"], inv["s_hi"]],
        formatted=f'{inv["s_lo"]*100:.1f}%--{inv["s_hi"]*100:.1f}%',
        where=["§2.5"],
        computed_from="Jeffreys interval on the hit proportion mapped back through the "
                      "model-implied null to an implied C>G supply rate; reported only to "
                      "show it is narrower than the uncorrected inversion, and NOT adopted, "
                      "because the inversion presumes supply controls hit composition",
        formula="solve pi_N3(C>G; s) = p for s, at p = Jeffreys lo/hi",
        inputs={"jeffreys": t["jeffreys_hit_proportion"], "beta": 1.0},
        code_fn="recompute_enu_class_test.py::solve_s")

    pd_ = t["p_detect_cg_germline"]
    reg(key="cg_detect_prob_germline", value=pd_["p_detect"],
        formatted=f'{pd_["p_detect"]:.2f}',
        where=["§2.1", "§2.5", "Methods library depth"],
        computed_from="detection probability of a C>G variant at the germline ENU rate "
                      "in a 1e7 Ba/F3 library; shows a rare class is still supplied",
        formula="1 - exp(-N*mu*S*r)",
        inputs={k2: pd_[k2] for k2 in ("N", "mu", "S", "routes")},
        code_fn="recompute_enu_class_test.py::main")

    reg(key="cg_expected_carriers_germline", value=pd_["expected_carriers"],
        formatted=f'{pd_["expected_carriers"]:.1f}',
        where=["§2.1", "§2.5"],
        computed_from="expected carriers of a given C>G variant at the germline ENU rate",
        formula="N*mu*S*r",
        inputs={k2: pd_[k2] for k2 in ("N", "mu", "S", "routes")},
        code_fn="recompute_enu_class_test.py::main")


# ---------------- Two-nucleotide reachability ----------------
def compute_two_hit() -> None:
    """Replaces the 'impossible at any library depth' claim with arithmetic."""
    t = json.loads((DATA / "derived" / "two_hit_reachability.json").read_text())
    df = pd.read_csv(DATA / "derived" / "two_hit_reachability.csv")

    one = df[(df.mu == t["mu_primary"]) & (df.rounds == 1)]
    reg(key="two_hit_carriers_one_round",
        value=[float(one.expected_carriers.min()), float(one.expected_carriers.max())],
        formatted=f"{one.expected_carriers.min():.3f}--{one.expected_carriers.max():.3f}",
        where=["§2.2"],
        computed_from="expected carriers of R68L or V8M after a single mutagenesis round "
                      "in a 1e7 library, across ENU and Mutazyme II",
        formula="sum over routes of N*mu^2*S(c1)*S(c2)*2",
        inputs={"N": t["N"], "mu": t["mu_primary"], "rounds": 1},
        code_fn="two_hit_reachability.py::main")

    reg(key="two_hit_presence_one_round",
        value=float(one.p_at_least_one.max()),
        formatted=fmt_pct(one.p_at_least_one.max()),
        where=["supporting: \u00a72.2 states this as a presence probability below 2%"],
        computed_from="largest presence probability for R68L or V8M after one round",
        formula="1 - exp(-lambda)",
        inputs={"N": t["N"], "mu": t["mu_primary"], "rounds": 1},
        code_fn="two_hit_reachability.py::main")

    twenty = df[(df.mu == t["mu_primary"]) & (df.rounds == 20)]
    reg(key="two_hit_carriers_twenty_rounds",
        value=[float(twenty.expected_carriers.min()), float(twenty.expected_carriers.max())],
        formatted=f"{twenty.expected_carriers.min():.1f}--{twenty.expected_carriers.max():.1f}",
        where=["§2.2"],
        computed_from="expected carriers of R68L or V8M after twenty rounds",
        formula="lambda_simultaneous * R(R+1)/2",
        inputs={"N": t["N"], "mu": t["mu_primary"], "rounds": 20},
        code_fn="two_hit_reachability.py::main")

    folds = [t[m][v]["fold_below_one_nt"]
             for m in ("ENU", "MutazymeII") for v in ("R68L", "V8M")]
    reg(key="two_hit_fold_suppression",
        value=[min(folds), max(folds)],
        formatted=f"{min(folds):.0f}- to {max(folds):.0f}-fold",
        where=["§2.2", "Methods multi-nucleotide reachability", "Discussion"],
        computed_from="suppression of R68L and V8M relative to a one-step variant of "
                      "median class supply in the same library, ENU and Mutazyme II",
        formula="lambda_1nt(median class) / lambda_2nt",
        inputs={"N": t["N"], "mu": t["mu_primary"], "rounds": 20},
        code_fn="two_hit_reachability.py::main")

    reg(key="r68l_two_nt_routes", value=t["routes"]["R68L"],
        formatted=str(t["routes"]["R68L"]),
        where=["§2.2"],
        computed_from="two-substitution routes from AGG to a leucine codon",
        formula="enumerate codons at Hamming distance 2 translating to Leu",
        inputs={"wt_codon": "AGG", "target_aa": "L"},
        code_fn="two_hit_reachability.py::two_nt_routes")

    reg(key="v8m_two_nt_routes", value=t["routes"]["V8M"],
        formatted=str(t["routes"]["V8M"]),
        where=["§2.2"],
        computed_from="two-substitution routes from GTA to the methionine codon",
        formula="enumerate codons at Hamming distance 2 translating to Met",
        inputs={"wt_codon": "GTA", "target_aa": "M"},
        code_fn="two_hit_reachability.py::two_nt_routes")


# ---------------- Permutation-test resolution ----------------
def compute_power() -> None:
    """What the class-label permutation tests at n = 3-11 can actually detect."""
    b = pd.read_csv(DATA / "derived" / "benchmark_r3.csv")
    f = b[b.ranking.str.startswith("Framework")]

    reg(key="perm_mdg_range",
        value=[float(f.min_detectable_gain.min()), float(f.min_detectable_gain.max())],
        formatted=f"{f.min_detectable_gain.min()*100:.0f}--{f.min_detectable_gain.max()*100:.0f} pp",
        where=["§2.5", "Table 2"],
        computed_from="97.5th percentile of the permutation null's absolute top-100 recall "
                      "gain, across all six anchor-proxy combinations",
        formula="quantile(|null gain|, 0.975) per anchor",
        inputs={"n_iter": 5000, "anchors": int(len(f))},
        code_fn="recompute_benchmark_r3.py::permutation_p")

    reg(key="perm_gain_quantum_range",
        value=[float(f.gain_quantum.min()), float(f.gain_quantum.max())],
        formatted=f"{f.gain_quantum.min()*100:.0f}--{f.gain_quantum.max()*100:.0f} pp",
        where=["§2.5"],
        computed_from="smallest observable recall change, 1/n_hits, across anchors",
        formula="1 / n_hits",
        inputs={"n_hits": sorted(set(int(x) for x in f.n_hits))},
        code_fn="recompute_benchmark_r3.py::permutation_p")

    reg(key="perm_observed_gain_range",
        value=[float(f.observed_gain.min()), float(f.observed_gain.max())],
        formatted=f"{f.observed_gain.min()*100:.0f} to {f.observed_gain.max()*100:.0f} pp",
        where=["supporting: \u00a72.5 states the bound, "
               "perm_max_abs_gain, rather than the signed range"],
        computed_from="observed framework-minus-sigma top-100 recall gains across anchors",
        formula="recall@100(logS + sigma) - recall@100(sigma)",
        inputs={"K": 100},
        code_fn="recompute_benchmark_r3.py::run_anchor")

    # This is the quantity the prose actually quotes: the largest absolute
    # change in top-100 recall that adding the chemistry term produces on any
    # anchor. Registering the signed range alone left the stated bound
    # unchecked.
    reg(key="perm_max_abs_gain",
        value=float(f.observed_gain.abs().max()),
        formatted=f"{f.observed_gain.abs().max()*100:.0f} pp",
        where=["\u00a72.5"],
        computed_from="largest absolute framework-minus-sigma top-100 recall "
                      "gain over all six anchor-proxy combinations",
        formula="max |recall@100(logS + sigma) - recall@100(sigma)|",
        inputs={"K": 100, "anchors": int(len(f))},
        code_fn="recompute_benchmark_r3.py::run_anchor")

    # The per-anchor p-values are quoted individually in the prose, so they are
    # registered individually. Quoting them only as a range let a stale value
    # (an ALK proxy at p = 0.78, left behind by the codon-table rebuild) sit in
    # the text while the range it belonged to still matched.
    PERM_ANCHOR_KEYS = {
        "Doi 2022 (ALK Mutazyme \u00d7 alec)": "doi_alec",
        "Nuvalent (ALK ENU \u00d7 alec proxy)": "nuvalent_alec",
        "Nuvalent (ALK ENU \u00d7 lor proxy)": "nuvalent_lor",
        "Nuvalent (ALK ENU \u00d7 zot proxy)": "nuvalent_zot",
        "Awad clinical (KRAS Mut \u00d7 sot)": "awad_sot",
        "Awad clinical (KRAS Mut \u00d7 ada)": "awad_ada",
    }
    seen = set()
    for r in f.itertuples():
        slug = PERM_ANCHOR_KEYS.get(r.anchor)
        assert slug is not None, f"unregistered benchmark anchor: {r.anchor!r}"
        seen.add(slug)
        reg(key=f"perm_p_{slug}", value=float(r.perm_p_at_100),
            formatted=f"{r.perm_p_at_100:.2f}",
            where=["\u00a72.5", "Table 2"],
            computed_from=f"class-label permutation on the {r.anchor} anchor",
            formula="(#{|null gain| >= |observed gain|} + 1) / (n_iter + 1)",
            inputs={"n_iter": 5000, "K": 100, "n_hits": int(r.n_hits)},
            code_fn="recompute_benchmark_r3.py::permutation_p")
    assert seen == set(PERM_ANCHOR_KEYS.values()), sorted(
        set(PERM_ANCHOR_KEYS.values()) - seen)

    reg(key="perm_n_hits_range",
        value=[int(f.n_hits.min()), int(f.n_hits.max())],
        formatted=f"{f.n_hits.min()} to {f.n_hits.max()}",
        where=["\u00a72.5"],
        computed_from="class-forced observed hits per anchor",
        formula="count of named anchor hits present in the class-forced universe",
        inputs={"anchors": int(len(f))},
        code_fn="recompute_benchmark_r3.py::run_anchor")


# ---------------- Sigma transport ----------------
def compute_transport() -> None:
    """External test of sigma transport against measured Ba/F3 potency."""
    t = json.loads((DATA / "derived" / "sigma_transport.json").read_text())
    for drug, v in t["transport"].items():
        if "spearman_rho" not in v:
            continue
        reg(key=f"transport_rho_{drug}", value=v["spearman_rho"],
            formatted=f'{v["spearman_rho"]:.2f}',
            where=["\u00a72.6", "Table 4"],
            computed_from=f"Spearman rank correlation between atlas sigma and measured "
                          f"Ba/F3 log10 fold resistance for {drug}, over {v['n']} shared variants",
            formula="spearmanr(sigma_atlas, median log10(IC50_variant / IC50_WT) within study)",
            inputs={"n": v["n"], "boot95": v["boot95"], "p": v["p"],
                    "median_inter_study_spread_log10": v["median_inter_study_spread_log10"]},
            code_fn="sigma_transport_validation.py::main")
    for drug, v in t["proxy"].items():
        if "spearman_rho" not in v:
            continue
        reg(key=f"proxy_rho_{drug}", value=v["spearman_rho"],
            formatted=f'{v["spearman_rho"]:.2f}',
            where=["\u00a72.6"],
            computed_from=f"Spearman correlation between {drug} atlas sigma and measured "
                          f"Ba/F3 neladalkib fold resistance; ranks the proxy choices",
            formula="spearmanr(sigma_atlas[drug], log10 FC under neladalkib)",
            inputs={"n": v["n"], "boot95": v["boot95"], "p": v["p"]},
            code_fn="sigma_transport_validation.py::main")
    for drug, ratio in t["atlas_dose_over_wt_ic50"].items():
        reg(key=f"atlas_dose_ratio_{drug}", value=ratio, formatted=f"{ratio:.1f}x",
            where=["\u00a72.6", "Table 4"],
            computed_from=f"atlas assay concentration for {drug} divided by the measured "
                          f"Ba/F3 wild-type IC50",
            formula="atlas_dose_nM / median(WT Ba/F3 IC50_nM)",
            inputs={"dose_nM": t["atlas_dose_nM"][drug],
                    "wt_ic50_nM": t["wt_ic50_nM_by_drug"][drug]},
            code_fn="sigma_transport_validation.py::main")
    reg(key="transport_n_studies", value=t["n_studies"], formatted=str(t["n_studies"]),
        where=["\u00a72.6", "Methods external test of sigma transport"],
        computed_from="independent published studies contributing Ba/F3 ALK potency measurements",
        formula="nunique(paper_id)", inputs={}, code_fn="sigma_transport_validation.py::load_potency")


# ---------------- Internal dose transport ----------------
def compute_dose_transport() -> None:
    """Oh atlas internal dose series: how much the sigma ranking moves."""
    t = json.loads((DATA / "derived" / "oh_dose_transport.json").read_text())
    h = t["headline"]
    reg(key="dose_rho_low_high_range", value=h["rho_low_vs_high_range"],
        formatted=f'{h["rho_low_vs_high_range"][0]:.2f} to {h["rho_low_vs_high_range"][1]:.2f}',
        where=["\u00a72.6"],
        computed_from="Spearman between a drug's lowest and highest assayed concentration, "
                      "across the seven multi-dose series in Oh Tables S3, S5 and S9",
        formula="spearmanr(score_lowest_dose, score_highest_dose) per series",
        inputs={"n_series": h["n_series"], "median": h["rho_low_vs_high_median"]},
        code_fn="oh_dose_transport.py::main")

    reg(key="dose_rho_low_high_median", value=h["rho_low_vs_high_median"],
        formatted=f'{h["rho_low_vs_high_median"]:.2f}', where=["\u00a72.6"],
        computed_from="median across the seven dose series of the lowest-to-highest "
                      "rank correlation",
        formula="median over series", inputs={"n_series": h["n_series"]},
        code_fn="oh_dose_transport.py::main")

    reg(key="dose_topk_overlap_range", value=h["topk_overlap_range"],
        formatted=f'between {h["topk_overlap_range"][0]} and {h["topk_overlap_range"][1]} of the top 25',
        where=["\u00a72.6"],
        computed_from="variants shared between the top 25 at a drug's lowest and highest dose",
        formula="|top25(lowest) intersect top25(highest)|",
        inputs={"top_k": t["top_k"]}, code_fn="oh_dose_transport.py::main")

    reg(key="dose_rho_adjacent_median", value=h["rho_adjacent_median"],
        formatted=f'{h["rho_adjacent_median"]:.2f}', where=["\u00a72.6"],
        computed_from="median rank correlation between adjacent concentration steps",
        formula="median over all adjacent dose pairs", inputs={},
        code_fn="oh_dose_transport.py::main")

    b = t["controls"]["B_cross_drug"]
    reg(key="dose_cross_drug_matched_rho", value=b["matched_dose_rank"]["median_rho"],
        formatted=f'{b["matched_dose_rank"]["median_rho"]:.2f}', where=["\u00a72.6"],
        computed_from="median rank correlation between DIFFERENT drugs assayed at matched "
                      "dose rank; exceeds the same-drug across-dose figure, indicating "
                      "shared batch structure",
        formula="median spearmanr over cross-drug pairs with equal dose rank",
        inputs={"n_pairs": b["matched_dose_rank"]["n_pairs"]},
        code_fn="oh_dose_transport.py::main")

    reg(key="dose_cross_drug_unmatched_rho", value=b["unmatched_dose_rank"]["median_rho"],
        formatted=f'{b["unmatched_dose_rank"]["median_rho"]:.2f}', where=["\u00a72.6"],
        computed_from="same, for cross-drug pairs at unmatched dose rank",
        formula="median spearmanr over cross-drug pairs with unequal dose rank",
        inputs={"n_pairs": b["unmatched_dose_rank"]["n_pairs"]},
        code_fn="oh_dose_transport.py::main")

    a = t["controls"]["A_atlas_vs_dose"]
    best = [max(v["rho_by_dose"].values()) for v in a.values()]
    reg(key="atlas_vs_own_repeat_rho", value=[min(best), max(best)],
        formatted=f"{min(best):.2f}--{max(best):.2f}", where=["\u00a72.6"],
        computed_from="agreement between the main single-dose atlas and its best-matching "
                      "sub-library condition: two runs of nominally the same assay",
        formula="max over doses of spearmanr(S2 score, S5 score)",
        inputs={d: v["best_match_dose"] for d, v in a.items()},
        code_fn="oh_dose_transport.py::main")

    d4 = t["controls"]["D_range_restriction"]
    ratios = [v["sd_ratio_window_over_full"] for v in d4["by_drug"].values()]
    reg(key="dose_window_sd_ratio", value=ratios,
        formatted=f"{min(ratios):.1f} to {max(ratios):.1f}", where=["\u00a72.6"],
        computed_from="sigma dispersion in the exon-23 sub-library window relative to the "
                      "full atlas; above 1 rules out attenuation by range restriction",
        formula="sd(window) / sd(full atlas), per drug",
        inputs={"window": t["window"], "n_window": d4["n_window_missense"],
                "n_full": d4["n_full_missense"]},
        code_fn="oh_dose_transport.py::main")


# ---------------- Doi anchor per-variant ranks ----------------
def compute_doi_anchor_ranks() -> None:
    """Per-variant simulated ranks for the primary retrospective anchor.

    These were stated in the prose but not registered, so the rebuild of the
    ALK codon table moved I1171N from 106 to 118 without the audit noticing.
    """
    d = pd.read_csv(DATA / "derived" / "oh_decomposition_all_drugs.csv")
    u = d[(d.drug == "alectinib") & (d.mutagen == "MutazymeII")]
    doi = ["L1196M", "L1198F", "I1171T", "V1180L", "I1171N"]
    sub = u[u.aa_change.isin(doi)].set_index("aa_change")
    for v in doi:
        r = int(sub.loc[v, "rank_simulated"])
        reg(key=f"doi_{v.lower()}_rank_sim", value=r, formatted=str(r),
            where=["\u00a72.5"],
            computed_from=f"{v} rank under the Mutazyme-simulated alectinib ranking",
            formula="rank of log S + sigma, descending, over the class-forced ALK universe",
            inputs={"variant": v, "mutagen": "MutazymeII", "drug": "alectinib"},
            code_fn="notebook/05_prospective_predictions.py")
    worst = int(sub.rank_simulated.max())
    bound = 10 * ((worst // 10) + 1)
    reg(key="doi_anchor_topN_bound", value=bound, formatted=f"top-{bound}",
        where=["Background", "\u00a72.4", "\u00a72.5"],
        computed_from="smallest round bound containing all five class-forced Doi hits "
                      "in the Mutazyme-simulated ranking",
        formula="ceil(max rank / 10) * 10",
        inputs={"max_rank": worst}, code_fn="compute_doi_anchor_ranks")


# ---------------- MET audit ----------------
def compute_met_audit() -> None:
    """Tiedt 2011 ENU MET screen audited against the Estevam atlas."""
    t = json.loads((DATA / "derived" / "met_audit.json").read_text())

    h = t["two_screen_homogeneity"]
    reg(key="met_two_screen_chi2", value=h["chi2"], formatted=f'{h["chi2"]:.1f}',
        where=["\u00a72.7", "Fig 5A"],
        computed_from="chi-square homogeneity of the six-class clone composition between "
                      "the NVP-BVU972 and AMG 458 screens, which used the same ENU "
                      "protocol in the same cells",
        formula="scipy.stats.chi2_contingency(2 x 6 clone counts)",
        inputs={"dof": h["dof"], "p": h["p"], "table": h["table"]},
        code_fn="met_audit.py::main")

    reg(key="met_two_screen_p", value=h["p"], formatted=f'{h["p"]:.3f}',
        where=["Abstract", "\u00a72.7", "Fig 5A"],
        computed_from="p-value of that homogeneity test; the two screens differ, so the "
                      "class composition of a hit list is set by the drug, not the mutagen",
        formula="scipy.stats.chi2_contingency", inputs={"dof": h["dof"]},
        code_fn="met_audit.py::main")

    y = t["y1230"]
    reg(key="met_y1230_recovered", value=y["n_recovered"],
        formatted=f'{y["n_recovered"]} of {y["n_single_nt_missense_routes"]}',
        where=["\u00a72.7", "Fig 5B"],
        computed_from="single-nucleotide missense routes out of the Y1230 wild-type codon "
                      "that the NVP-BVU972 screen recovered",
        formula="enumerate one-step routes from the codon; intersect with Tiedt Table S1",
        inputs={"wt_codon": y["wt_codon"], "n_clones": y["n_clones_at_codon"]},
        code_fn="met_audit.py::main")

    lo, hi = y["supply_range_of_those_classes"]
    reg(key="met_y1230_supply_fold", value=hi / lo, formatted=f"{hi/lo:.1f}-fold",
        where=["\u00a72.7", "Fig 5B"],
        computed_from="ratio of the largest to smallest germline ENU supply among the "
                      "substitution classes those Y1230 routes span",
        formula="max(S_ENU(c)) / min(S_ENU(c)) over the routes' classes",
        inputs={"classes": y["classes_spanned"], "lo": lo, "hi": hi},
        code_fn="met_audit.py::main")

    # Within the Y1230 codon the supply term does predict abundance, which is
    # the one place in the paper where the weighting half of the decomposition
    # is tested at nucleotide resolution rather than at the ranking level.
    _tg = [k for k, v in y["per_variant"].items()
           if all(r["class"] == "T>G" for r in v["routes"])]
    _tg_exp = sum(y["expected_clones_under_supply"][k] for k in _tg)
    _tg_obs = sum(y["per_variant"][k]["clones_observed"] for k in _tg)
    reg(key="met_y1230_tg_expected_clones", value=_tg_exp,
        formatted=f"{_tg_exp:.1f}",
        where=["§2.7"],
        computed_from="clones the two T>G routes out of Y1230 are predicted to "
                      "carry under germline ENU supply alone",
        formula="sum of expected_clones_under_supply over the T>G routes",
        inputs={"routes": sorted(_tg)},
        code_fn="met_audit.py::main")

    reg(key="met_y1230_tg_observed_clones", value=_tg_obs,
        formatted=str(_tg_obs),
        where=["§2.7"],
        computed_from="clones the two T>G routes out of Y1230 actually carried "
                      "in the Tiedt NVP-BVU972 screen",
        formula="sum of clones_observed over the T>G routes",
        inputs={"routes": sorted(_tg)},
        code_fn="met_audit.py::main")

    reg(key="met_y1230_missing", value=list(y["missing_variant_atlas_rank"]),
        formatted=", ".join(y["missing_variant_atlas_rank"]) or "none",
        where=["\u00a72.7"],
        computed_from="one-step Y1230 variants the screen did not recover, with their "
                      "rank in the chemistry-free atlas",
        formula="routes minus recovered", inputs=y["missing_variant_atlas_rank"],
        code_fn="met_audit.py::main")

    for screen, v in t["screens"].items():
        key = screen.split()[0].lower().replace("-", "_")
        reg(key=f"met_{key}_clones", value=v["n_clones"], formatted=str(v["n_clones"]),
            where=["\u00a72.7", "Fig 5A"],
            computed_from=f"resistant clones reported for the {screen} screen",
            formula="sum of clone counts in Tiedt Supplementary Table",
            inputs={"n_base_changes": v["n_distinct_base_changes"],
                    "n_aa_changes": v["n_distinct_aa_changes"]},
            code_fn="met_audit.py::main")

        reg(key=f"met_{key}_median_rank", value=v["median_hit_rank"],
            formatted=str(v["median_hit_rank"]),
            where=["\u00a72.7"],
            computed_from=f"median chemistry-free atlas rank of the {screen} hits",
            formula="median(rank by delta sigma) over reported hits",
            inputs={"atlas_n": v["atlas_n"], "found": v["hits_found_in_atlas"],
                    "total": v["hits_total"]},
            code_fn="met_audit.py::main")

        reg(key=f"met_{key}_twostep_top100", value=v["top100_reachability"].get("two-step", 0),
            formatted=str(v["top100_reachability"].get("two-step", 0)),
            where=["\u00a72.7", "Fig 5C"],
            computed_from=f"variants in the chemistry-free top 100 for {screen} that need "
                          f"two or more nucleotide changes and so lie outside one-step reach",
            formula="classify each top-100 variant by one-step routes from its codon",
            inputs=v["top100_reachability"], code_fn="met_audit.py::main")

    reg(key="met_atlas_size", value=t["screens"][list(t["screens"])[0]]["atlas_n"],
        formatted=f'{t["screens"][list(t["screens"])[0]]["atlas_n"]:,}',
        where=["\u00a72.7"],
        computed_from="MET kinase-domain missense variants scored in the Estevam atlas "
                      "on the primary background",
        formula="rows of the wide table for that background",
        inputs={"background": t["background_primary"]}, code_fn="met_audit.py::main")


# ---------------- Software environment ----------------
def compute_env() -> None:
    e = json.loads((DATA / "derived" / "env_versions.json").read_text())
    reg(key="software_versions", value=e,
        formatted=f"Python {e['python']} with " +
                  ", ".join(f"{k} {v}" for k, v in e["packages"].items()),
        where=["Methods software and reproducibility"],
        computed_from="emitted from the running interpreter, not recorded by hand",
        formula="platform.python_version() and each package __version__",
        inputs={}, code_fn="env_versions.py::main")


# ---------------- Retrospective validation set sizes ----------------
def compute_retrospective_screens() -> None:
    # Takei: 5 V2 partners
    takei_v2 = ["L1256F", "E1129V", "E1129K", "C1156Y", "T1151R"]
    takei_c_t_class_forced = ["L1256F", "E1129K", "C1156Y"]
    reg(key="takei_v2_count", value=len(takei_v2), formatted="five",
        where=["§2.5", "Discussion §D7", "Fig 4B caption", "Methods"],
        computed_from="V2 compound partners recovered at 500 nM neladalkib on I1171N fixed-V1 background",
        formula="|Takei V2 hits at 500 nM neladalkib|",
        inputs={"members": takei_v2}, code_fn="compute_retrospective_screens")
    reg(key="takei_v2_c_t_count", value=len(takei_c_t_class_forced), formatted="three",
        where=["Abstract", "§2.5", "Fig 4B caption"],
        computed_from="Takei V2 partners class-forced into C>T (Mutazyme II's dominant class)",
        formula="|{v in Takei V2: forced_class(v) == 'C>T'}|",
        inputs={"members": takei_c_t_class_forced}, code_fn="compute_retrospective_screens")

    # Nuvalent: 13 aggregate, 12 class-forced (R1192P outside atlas)
    nuvalent = ["I1171N", "I1171S", "I1171T", "F1174C", "F1174V", "V1180L",
                "T1151K", "T1151R", "E1129K", "R1192P", "L1196M", "G1269A", "G1202R"]
    reg(key="nuvalent_total", value=len(nuvalent), formatted="thirteen",
        where=["§2.5", "Fig 4C caption", "Methods"],
        computed_from="single-amino-acid variants reported in Nuvalent AACR 2025 poster aggregate across arms",
        formula="|Nuvalent single-amino-acid hits (aggregate)|",
        inputs={"members": nuvalent}, code_fn="compute_retrospective_screens")
    # R1192P is outside the Oh atlas (not in codons.yaml); G1202R is multi-class
    nuvalent_classified = [v for v in nuvalent if v not in ("R1192P", "G1202R")]
    reg(key="nuvalent_class_forced", value=len(nuvalent_classified), formatted="twelve",
        where=["§2.5", "Fig 4C caption"],
        computed_from="Nuvalent hits that are class-forced single-nt reachable (R1192P outside atlas, G1202R multi-class)",
        formula="|Nuvalent hits ∩ class-forced|",
        inputs={"members": nuvalent_classified}, code_fn="compute_retrospective_screens")

    # Three C>G hits
    nuvalent_cg = ["V1180L", "T1151R", "G1269A"]
    reg(key="nuvalent_cg_count", value=len(nuvalent_cg), formatted="three",
        where=["Abstract", "§2.5", "Fig 4C caption"],
        computed_from="Nuvalent hits class-forced into C>G (ENU 0.2 % class per germline spectrum)",
        formula="|{v in Nuvalent: forced_class(v) == 'C>G'}|",
        inputs={"members": nuvalent_cg}, code_fn="compute_retrospective_screens")

    # 9 of 12 at ≥ 1 %
    nuvalent_well = [v for v in nuvalent_classified if v not in nuvalent_cg]
    reg(key="nuvalent_well_supplied_count", value=len(nuvalent_well), formatted="nine",
        where=["§2.5", "Fig 4C caption"],
        computed_from="Nuvalent class-forced hits whose forced class ENU supplies at ≥ 1 %",
        formula="|Nuvalent class-forced| - |Nuvalent C>G|",
        inputs={"members": nuvalent_well}, code_fn="compute_retrospective_screens")


# ---------------- Top-K recovery curves (summary numbers for Abstract) ----------------
def compute_topk_recovery() -> None:
    """Summary stats quoted in abstract/§2.4.

    Pulls recall@K values AND 95% Wilson CIs from the authoritative
    topk_recovery.csv (produced by scripts/recompute_topk_with_ci.py; R9).
    """
    df = pd.read_csv(DATA / "derived" / "topk_recovery.csv")
    # NB: mutagen labels in the derived CSV are {"ENU", "MutII", "BE", "Twist"}
    SNV_MUTAGENS = ["ENU", "MutII"]

    def _per_row_dict(sub):
        return {f"{r.drug}/{r.mutagen}": float(r.recovery) for r in sub.itertuples()}

    def _fmt_ci(lo: float, hi: float) -> str:
        return f"[{int(round(lo*100))}\\%--{int(round(hi*100))}\\%]"

    # ALK two SNV mutagens at top-100 (what the abstract quotes)
    alk100 = df[(df["gene"] == "ALK") & (df["mutagen"].isin(SNV_MUTAGENS)) & (df["K"] == 100)]
    lo, hi = float(alk100["recovery"].min()), float(alk100["recovery"].max())
    reg(key="alk_snv_top100_recovery_min", value=lo, formatted=fmt_pct(lo, 0),
        where=["§2.4"],
        computed_from="min top-100 recovery across ALK drugs for ENU and Mutazyme II",
        formula="min_{drug, mutagen in {ENU, MutII}} recovery(K=100, drug, mutagen)",
        inputs={"source": "data/derived/topk_recovery.csv", "gene": "ALK", "K": 100,
                "per_row": _per_row_dict(alk100)},
        code_fn="compute_topk_recovery")
    reg(key="alk_snv_top100_recovery_max", value=hi, formatted=fmt_pct(hi, 0),
        where=["§2.4"],
        computed_from="max top-100 recovery across ALK drugs for ENU and Mutazyme II",
        formula="max_{drug, mutagen in {ENU, MutII}} recovery(K=100, drug, mutagen)",
        inputs={"source": "data/derived/topk_recovery.csv", "gene": "ALK", "K": 100,
                "per_row": _per_row_dict(alk100)},
        code_fn="compute_topk_recovery")
    # Per-drug ALK top-100 (used in §2.4)
    alk100_mut = alk100[alk100["mutagen"] == "MutII"].set_index("drug")
    alk100_enu = alk100[alk100["mutagen"] == "ENU"].set_index("drug")
    for drug in ("alectinib", "lorlatinib", "zotizalkib"):
        v_mut = float(alk100_mut.loc[drug, "recovery"])
        v_enu = float(alk100_enu.loc[drug, "recovery"])
        lo_mut, hi_mut = float(alk100_mut.loc[drug, "ci_lo"]), float(alk100_mut.loc[drug, "ci_hi"])
        lo_enu, hi_enu = float(alk100_enu.loc[drug, "ci_lo"]), float(alk100_enu.loc[drug, "ci_hi"])
        short = {"alectinib": "alec", "lorlatinib": "lor", "zotizalkib": "zot"}[drug]
        reg(key=f"alk_{short}_mut_top100", value=v_mut,
            formatted=f"{fmt_pct(v_mut,0)} {_fmt_ci(lo_mut, hi_mut)}",
            where=["supporting: §2.4 quotes the range across drugs and mutagens"],
            computed_from=f"Mutazyme-simulated top-100 recovery for ALK, {drug}",
            formula="recovery(K=100, drug, MutII)",
            inputs={"source": "data/derived/topk_recovery.csv"},
            code_fn="compute_topk_recovery")
        reg(key=f"alk_{short}_enu_top100", value=v_enu,
            formatted=f"{fmt_pct(v_enu,0)} {_fmt_ci(lo_enu, hi_enu)}",
            where=["supporting: §2.4 quotes the range across drugs and mutagens"],
            computed_from=f"ENU-simulated top-100 recovery for ALK, {drug}",
            formula="recovery(K=100, drug, ENU)",
            inputs={"source": "data/derived/topk_recovery.csv"},
            code_fn="compute_topk_recovery")

    # KRAS G12C
    kras100 = df[(df["gene"] == "KRAS G12C") & (df["mutagen"].isin(SNV_MUTAGENS)) & (df["K"] == 100)]
    lo, hi = float(kras100["recovery"].min()), float(kras100["recovery"].max())
    reg(key="kras_snv_top100_recovery_min", value=lo, formatted=fmt_pct(lo, 0),
        where=["§2.4"],
        computed_from="min top-100 recovery across KRAS G12C drugs for ENU and Mutazyme II",
        formula="min_{drug, mutagen in {ENU, MutII}} recovery(K=100, drug, mutagen)",
        inputs={"source": "data/derived/topk_recovery.csv", "gene": "KRAS G12C", "K": 100,
                "per_row": _per_row_dict(kras100)},
        code_fn="compute_topk_recovery")
    reg(key="kras_snv_top100_recovery_max", value=hi, formatted=fmt_pct(hi, 0),
        where=["§2.4"],
        computed_from="max top-100 recovery across KRAS G12C drugs for ENU and Mutazyme II",
        formula="max_{drug, mutagen in {ENU, MutII}} recovery(K=100, drug, mutagen)",
        inputs={"source": "data/derived/topk_recovery.csv", "gene": "KRAS G12C", "K": 100,
                "per_row": _per_row_dict(kras100)},
        code_fn="compute_topk_recovery")

    # Base editing — top-10 recovery for each ALK drug (NOT zero for every drug!)
    be_alk10 = df[(df["gene"] == "ALK") & (df["mutagen"] == "BE") & (df["K"] == 10)]
    be_alk10_per_drug = be_alk10.set_index("drug")["recovery"].to_dict()
    be_alk10_min = min(be_alk10_per_drug.values())
    be_alk10_max = max(be_alk10_per_drug.values())
    reg(key="be_alk_top10_min", value=be_alk10_min, formatted=fmt_pct(be_alk10_min, 0),
        where=["§2.4"],
        computed_from="min CBE+ABE top-10 recovery across ALK drugs (alectinib shows 0)",
        formula="min_{drug in {alec, lor, zot}} recovery(K=10, drug, BE)",
        inputs={"source": "data/derived/topk_recovery.csv", "per_drug": be_alk10_per_drug},
        code_fn="compute_topk_recovery")
    reg(key="be_alk_top10_max", value=be_alk10_max, formatted=fmt_pct(be_alk10_max, 0),
        where=["§2.4"],
        computed_from="max CBE+ABE top-10 recovery across ALK drugs",
        formula="max_{drug in {alec, lor, zot}} recovery(K=10, drug, BE)",
        inputs={"source": "data/derived/topk_recovery.csv", "per_drug": be_alk10_per_drug},
        code_fn="compute_topk_recovery")
    # Fraction of each ALK drug's chemistry-free top-10 OUTSIDE what CBE+ABE can produce
    be_outside_min = 1 - be_alk10_max
    be_outside_max = 1 - be_alk10_min
    reg(key="be_alk_top10_outside_min", value=be_outside_min, formatted=fmt_pct(be_outside_min, 0),
        where=["§2.4"],
        computed_from="1 - max(BE top-10 recovery across ALK drugs)",
        formula="1 - max_{drug} recovery(K=10, drug, BE)",
        inputs={"source": "data/derived/topk_recovery.csv"},
        code_fn="compute_topk_recovery")
    reg(key="be_alk_top10_outside_max", value=be_outside_max, formatted=fmt_pct(be_outside_max, 0),
        where=["§2.4"],
        computed_from="1 - min(BE top-10 recovery across ALK drugs)",
        formula="1 - min_{drug} recovery(K=10, drug, BE)",
        inputs={"source": "data/derived/topk_recovery.csv"},
        code_fn="compute_topk_recovery")

    # BE ceiling in the practical rank range (K ≤ 100). The curve rises
    # toward 1.0 only as K approaches the full atlas size (by construction,
    # since the top-N includes all variants).
    be_both = df[df["mutagen"] == "BE"]
    for gene_name in ("ALK", "KRAS G12C"):
        tag = "alk" if gene_name == "ALK" else "kras"
        sub = be_both[(be_both["gene"] == gene_name) & (be_both["K"] <= 100)]
        peak = float(sub["recovery"].max())
        reg(key=f"be_{tag}_peak_at_k_le_100",
            value=peak, formatted=fmt_pct(peak, 0),
            where=["§2.4"],
            computed_from=f"max CBE+ABE recovery over K ≤ 100 across all {gene_name} drugs",
            formula=f"max_{{K ≤ 100, drug}} recovery(K, drug, BE) for gene={gene_name}",
            inputs={"source": "data/derived/topk_recovery.csv", "gene": gene_name, "K_max": 100},
            code_fn="compute_topk_recovery")


# ---------------- Output ----------------
def dump_json(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {n.key: asdict(n) for n in REGISTRY}
    path.write_text(json.dumps(payload, indent=2, default=str))


def render_markdown(path: Path) -> None:
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    lines = [
        "# Paper numbers — master book",
        "",
        f"*Last computed: {now}*",
        f"*Source script: `scripts/compute_paper_numbers.py`*",
        f"*Total registered: {len(REGISTRY)} entries*",
        "",
        "## Golden rule",
        "",
        "Every quantitative claim in `paper.tex` comes from this manifest. "
        "If a number appears in the paper but is not registered below, that is a bug: "
        "either the number is wrong, or this script needs a new entry. "
        "Numbers from citations (years, DOIs, journal volume/issue) are not registered here; "
        "they live in `refs.bib`.",
        "",
        "---",
        "",
        "## Index by key",
        "",
    ]
    for n in sorted(REGISTRY, key=lambda x: x.key):
        lines.append(f"- `{n.key}` = **{n.formatted}** — {n.computed_from}")
    lines.extend(["", "---", "", "## Grouped by paper location", ""])

    sections: dict[str, list[N]] = {}
    for n in REGISTRY:
        for loc in n.where:
            sections.setdefault(loc, []).append(n)

    for loc in sorted(sections.keys()):
        lines.append(f"### {loc}")
        lines.append("")
        for n in sections[loc]:
            lines.append(f"**`{n.key}`** = {n.formatted}")
            lines.append("")
            lines.append(f"- Value: `{n.value}`")
            lines.append(f"- Formula: `{n.formula}`")
            lines.append(f"- Computed from: {n.computed_from}")
            if n.inputs:
                inputs_fmt = {k: (v if not isinstance(v, dict) else "{...}") for k, v in n.inputs.items()}
                lines.append(f"- Inputs: `{inputs_fmt}`")
            lines.append(f"- Code: `scripts/compute_paper_numbers.py::{n.code_fn}`")
            lines.append("")
        lines.append("")
    path.write_text("\n".join(lines))


def compute_reachability_track() -> None:
    """Figure 3C: per-residue reachability of the top-ranked substitution.

    Only base editing blocks anything here, because the inventory is scoped to
    class-forced substitutions and the two SNV chemistries supply all six
    classes at a non-zero rate. That scope is what the caption has to state.
    """
    trk = json.loads((DATA / "derived" / "reachability_track.json").read_text())
    # The drug each Figure 3C column plots. Must match C_GENES in figure_3.py,
    # otherwise the registered count and the drawn panel disagree.
    PANEL_DRUG = {"ALK": "alectinib", "KRAS G12C": "adagrasib",
                  "MET": "NVP-BVU972"}
    shares = []
    for gene, G in trk["by_gene"].items():
        tag = gene.split()[0].lower()
        n_res = G["n_residues"]
        reg(key=f"reach_residues_{tag}", value=n_res, formatted=str(n_res),
            where=["Fig 3C"],
            computed_from=f"{gene} residues carrying at least one class-forced "
                          f"variant in the chemistry-free atlas",
            formula="count of distinct residues in the class-forced set",
            inputs={"transcript": G["transcript"]},
            code_fn="build_reachability_track.py::build")

        reg(key=f"reach_class_forced_{tag}", value=G["n_class_forced_variants"],
            formatted=f'{G["n_class_forced_variants"]:,}',
            where=["Fig 3C"],
            computed_from=f"{gene} class-forced variants per drug arm",
            formula="single-nucleotide variants whose every route falls in one class",
            inputs={"transcript": G["transcript"]},
            code_fn="build_reachability_track.py::load_variants")

        for drug in G["drugs"]:
            rows = [r for r in G["track"]
                    if r["mutagen"] == "CBE+ABE" and r["drug"] == drug]
            out = sum(1 for r in rows if not r["reachable"])
            shares.append(out / len(rows))
            if drug == PANEL_DRUG[gene]:
                reg(key=f"reach_be_blocked_{tag}", value=out,
                    formatted=f"{out} of {n_res}",
                    where=["Fig 3C", "\u00a73.2"],
                    computed_from=f"{gene} residues whose top-sigma class-forced "
                                  f"substitution CBE+ABE cannot make, under {drug}",
                    formula="count of residues where S_CBE+ABE of the top-sigma "
                            "class-forced substitution is zero",
                    inputs={"drug": drug, "n_residues": n_res,
                            "note": "the arm Figure 3C plots"},
                    code_fn="build_reachability_track.py::build")

    lo, hi = min(shares), max(shares)
    reg(key="reach_be_blocked_share_range",
        value=[lo, hi], formatted=f"{lo*100:.0f}\\%--{hi*100:.0f}\\%",
        where=["Fig 3C", "\u00a73.2"],
        computed_from="range of the CBE+ABE blocked-residue share across all "
                      "three genes and every drug arm",
        formula="min and max over (gene, drug) of blocked residues / residues",
        inputs={"n_arms": len(shares)},
        code_fn="build_reachability_track.py::build")


def compute_per_drug_top100() -> None:
    """The per-drug top-100 recoveries quoted in Results 3.4.

    These were prose-only until the reverse audit found all six disagreeing
    with topk_recovery.csv after the ALK codon-table rebuild moved the
    universe from 1,259 to 1,427. Registering them puts each one under the
    same check as every other number in the manuscript.
    """
    d = pd.read_csv(DATA / "derived" / "topk_recovery.csv")
    t = d[(d["K"] == 100) & (d["mutagen"].isin(["MutII", "ENU"]))]
    for r in t.itertuples():
        gene = "alk" if r.gene == "ALK" else "kras"
        mut = "mutii" if r.mutagen == "MutII" else "enu"
        reg(key=f"top100_{gene}_{mut}_{r.drug}",
            value=float(r.recovery),
            formatted=f"{r.recovery*100:.0f}\\% [{r.ci_lo*100:.0f}\\%--{r.ci_hi*100:.0f}\\%]",
            # The KRAS arms are quoted only as the 45%--51% range, so they are
            # registered for provenance but not checked against the prose.
            where=(["\u00a73.4"] if r.gene == "ALK"
                   else ["not reported individually: §3.4 quotes the range"]),
            computed_from=f"fraction of the chemistry-free top-100 that the "
                          f"{r.mutagen}-simulated {r.drug} screen also ranks in "
                          f"its top 100, with a Wilson 95% interval",
            formula="|sim top-100 ∩ chemistry-free top-100| / 100",
            inputs={"gene": r.gene, "N": int(r.N), "K": 100},
            code_fn="recompute_topk_with_ci.py::compute_curve")


def compute_alk_panel_exemplars() -> None:
    """The individual ALK variants Results 3.4 walks through by name.

    These were prose-only, and the ALK codon-table rebuild that moved the
    universe from 1,259 to 1,427 left four of them stale: L1196M's score,
    L1198F's chemistry-free rank, and both ranks for C1156Y and L1152R.
    Registering them means the next rebuild cannot do that silently.

    Only the values the prose actually quotes are checked against it; the
    rest are registered for provenance.
    """
    d = pd.read_csv(DATA / "derived" / "oh_decomposition_all_drugs.csv")
    m = d[(d["mutagen"] == "MutazymeII") & (d["drug"] == "alectinib")]
    QUOTES_SIGMA = {"L1196M", "L1198F", "V1180L"}
    UNQUOTED = ["not reported: the variant is named but its score is not"]
    for v in ["L1196M", "L1198F", "V1180L", "C1156Y", "L1152R"]:
        r = m[m["aa_change"] == v]
        assert len(r), f"{v} missing from the ALK decomposition"
        r = r.iloc[0]
        tag = v.lower()
        reg(key=f"alk_{tag}_sigma", value=float(r["alec_score"]),
            formatted=f'{r["alec_score"]:+.1f}',
            where=["\u00a73.4"] if v in QUOTES_SIGMA else UNQUOTED,
            computed_from=f"Oh atlas alectinib Resistance score for {v}",
            formula="sigma = Resistance score, chemistry-free",
            inputs={"forced_class": r["forced_class"]},
            code_fn="oh_decomposition")

        for which, col in (("cf", "rank_chemfree"), ("sim", "rank_simulated")):
            reg(key=f"alk_{tag}_{which}_rank", value=int(r[col]),
                formatted=str(int(r[col])), where=["\u00a73.4"],
                computed_from=(f"{v} {'chemistry-free' if which == 'cf' else 'Mutazyme II-simulated'} "
                               f"rank among the 1,427 class-forced ALK variants, alectinib"),
                formula="rank by sigma" if which == "cf"
                        else "rank by log S + beta*sigma at beta = 1",
                inputs={"forced_class": r["forced_class"],
                        "shift": int(r["rank_shift"])},
                code_fn="oh_decomposition")

    l = m[m["aa_change"] == "L1152R"].iloc[0]
    reg(key="alk_largest_panel_shift", value=int(l["rank_shift"]),
        formatted=f'+{int(l["rank_shift"])}', where=["Discussion"],
        computed_from="largest predicted rank shift among the ALK panel "
                      "variants under Mutazyme II with alectinib",
        formula="max over panel variants of rank_simulated - rank_chemfree",
        inputs={"variant": "L1152R", "forced_class": l["forced_class"]},
        code_fn="oh_decomposition")


def compute_be_panel_reachability() -> None:
    """How many KRAS G12C panel variants a base-editor pair can produce.

    This was prose-only and wrong: the manuscript said eight of twelve were
    excluded. Enumerating the one-step routes gives seven blocked by class,
    two needing a second nucleotide change, and three a CBE+ABE pair can
    actually make. V8L was the one the hand count missed, because both of
    its routes (GTA->TTA, GTA->CTA) fall in C>A and C>G.
    """
    import re as _re
    import yaml as _yaml
    sys.path.insert(0, str(ROOT))
    from src.harmonized_cohort import classify_amino_acid_change

    codons = _yaml.safe_load(
        (DATA / "reference" / "codons_kras.yaml").read_text()
    )["kras_4a_nm_033360_4"]["codons"]
    BE = {"C>A": 0.0, "C>G": 0.0, "C>T": 0.5,
          "T>A": 0.0, "T>C": 0.5, "T>G": 0.0}

    panel = pd.read_csv(PAPER / "suppl_s1_data.csv")
    panel = panel[panel["gene"] == "KRAS G12C"]
    supplied, class_blocked, two_step = [], [], []
    for r in panel.itertuples():
        m = _re.match(r"([A-Z])(\d+)([A-Z])", r.variant)
        pos = int(m.group(2))
        wt = codons.get(pos) or codons.get(str(pos))
        wt = wt["codon"] if isinstance(wt, dict) else wt
        if r.variant == "C12F":
            wt = "TGT"                      # the G12C background this panel sits on
        res = classify_amino_acid_change(m.group(1), m.group(3), wt)
        if not res["reachable_1nt"]:
            two_step.append(r.variant); continue
        classes = res.get("classes") or ([res["forced_class"]]
                                         if res.get("forced_class") else [])
        if max((BE.get(c, 0.0) for c in classes), default=0.0) > 0:
            supplied.append(r.variant)
        else:
            class_blocked.append(r.variant)

    total = len(panel)
    assert len(supplied) + len(class_blocked) + len(two_step) == total
    reg(key="be_kras_panel_class_blocked", value=len(class_blocked),
        formatted=f"{len(class_blocked)} of the {total}",
        # The Discussion writes the digits; the Abstract spells the same
        # count out as "seven of the twelve", so only one form is checked.
        where=["Discussion"],
        computed_from="KRAS G12C panel variants reachable in one nucleotide "
                      "step whose every route falls in a class CBE+ABE cannot make",
        formula="max(S_CBE+ABE(c)) over one-step routes equals zero",
        inputs={"variants": sorted(class_blocked)},
        code_fn="classify_amino_acid_change")

    reg(key="be_kras_panel_supplied", value=len(supplied),
        formatted=str(len(supplied)), where=["Discussion"],
        computed_from="KRAS G12C panel variants a CBE+ABE pair can produce",
        formula="max(S_CBE+ABE(c)) over one-step routes exceeds zero",
        inputs={"variants": sorted(supplied)},
        code_fn="classify_amino_acid_change")


def compute_remaining_prose_numbers() -> None:
    """Values Results and Discussion quote that no other block registered.

    Found by the reverse audit. Each is computed from the derived tables
    rather than copied from the prose.
    """
    # R68L is two-step, so it is absent from the class-forced decomposition.
    # Its scores come from the atlas directly.
    atlas = pd.read_csv(DATA / "raw" / "feng_2022" / "table_s1_zscores.csv")
    atlas["aa_change"] = (atlas["wt_aa"] + atlas["position"].astype(str)
                          + atlas["mut_aa"])
    hit = atlas[atlas["aa_change"] == "R68L"]
    assert len(hit) == 1, f"expected one R68L row in the atlas, got {len(hit)}"
    r68l = hit.iloc[0]
    for drug, col in (("sotorasib", "sotorasib_z"), ("adagrasib", "adagrasib_z")):
        reg(key=f"kras_r68l_z_{drug[:3]}", value=float(r68l[col]),
            formatted=f'{r68l[col]:+.2f}', where=["\u00a73.2"],
            computed_from=f"chemistry-free Z-score of R68L under {drug}",
            formula="Feng atlas average Z",
            inputs={"variant": "R68L"}, code_fn="feng_decomposition")

    met = json.loads((DATA / "derived" / "met_audit.json").read_text())
    for screen, key in (("NVP-BVU972", "nvp"), ("AMG 458", "amg")):
        sens = met["sensitivity"][f"{screen}_WT_background"]
        reg(key=f"met_wt_median_rank_{key}", value=sens["median_rank"],
            formatted=f'{sens["median_rank"]:,}', where=["\u00a73.7"],
            computed_from=f"median atlas rank of the {screen} hits on the "
                          f"wild-type MET background",
            formula="median rank by delta sigma, WT background",
            inputs={"found": sens["n_found"], "total": sens["n_total"]},
            code_fn="met_audit.py::main")

    amg = met["screens"]["AMG 458"]["observed_class_fraction"]
    reg(key="met_amg_ta_share", value=amg["T>A"],
        formatted=f'{amg["T>A"]*100:.0f}\\%', where=["\u00a73.7", "Fig 5A"],
        computed_from="share of AMG 458 clones whose substitution is T>A",
        formula="clones in class / total clones",
        inputs={"n_clones": met["screens"]["AMG 458"]["n_clones"]},
        code_fn="met_audit.py::main")

    y = met["y1230"]
    for v, d in sorted(y["per_variant"].items()):
        reg(key=f"met_y1230_clones_{v.lower()}", value=d["clones_observed"],
            formatted=str(d["clones_observed"]), where=["\u00a73.7", "Fig 5B"],
            computed_from=f"clones carrying {v} in the Tiedt NVP-BVU972 screen",
            formula="transcribed clone count",
            inputs={"classes": [r["class"] for r in d["routes"]]},
            code_fn="met_audit.py::main")

    tr = json.loads((DATA / "derived" / "sigma_transport.json").read_text())
    lor = tr["transport"]["lorlatinib"]["median_inter_study_spread_log10"]
    reg(key="lorlatinib_inter_study_spread", value=float(lor),
        formatted=f"{lor:.2f}", where=["\u00a73.6", "Table 4"],
        computed_from="median inter-study range of measured lorlatinib "
                      "potency on the same variant",
        formula="median over variants of max - min log10 fold change",
        inputs={"n": tr["transport"]["lorlatinib"]["n"]},
        code_fn="sigma_transport_validation.py")


def main() -> int:
    compute_resistance_panels()
    compute_mutagen_spectra()
    compute_atlas_sizes()
    compute_h95d_ranks()
    compute_nuvalent_binomial()
    compute_two_hit()
    compute_power()
    compute_transport()
    compute_dose_transport()
    compute_doi_anchor_ranks()
    compute_met_audit()
    compute_env()
    compute_retrospective_screens()
    compute_topk_recovery()
    compute_reachability_track()
    compute_per_drug_top100()
    compute_alk_panel_exemplars()
    compute_be_panel_reachability()
    compute_remaining_prose_numbers()

    json_out = DATA / "derived" / "paper_numbers.json"
    md_out = PAPER / "NUMBERS.md"
    dump_json(json_out)
    render_markdown(md_out)

    print(f"Registered {len(REGISTRY)} numbers.")
    print(f"  → {json_out}")
    print(f"  → {md_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
