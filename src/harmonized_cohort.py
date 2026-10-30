"""Ingestion and cohort filtering for Michele Ronco's harmonised ALK+ NSCLC cohort.

Three input tables, keyed on `patient_key`:
    - patients.csv       (one row per patient, 612 patients, 572 analysis_core)
    - mutations.csv      (3,554 non-synonymous calls across 1,858 genes)
    - treatment_lines.csv (954 lines, 823 ALK-TKI, in 366 patients)

Schema documented in `data/raw/harmonized_cohort/README.md` (shipped by MR).

Public interface
----------------
- `HarmonizedCohort`: dataclass carrying the three tables and provenance.
- `load_harmonized_cohort()`: read the three CSVs from disk.
- `build_luad_never_smoker_signature_cohort()`: never-smoker ALK+ patients,
  with ≥1 non-ALK non-synonymous mutation (feeds the LUAD never-smoker
  6-class spectrum).
- `build_crizotinib_refractory_cohort()`: ALK-TKI line is crizotinib and
  only crizotinib before the biopsy; ALK KD mutation documented at or
  after crizotinib progression.
- `build_alectinib_refractory_cohort()`: alectinib in the TKI sequence
  (prior crizotinib allowed); ALK KD mutation documented at or after
  alectinib progression.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

__all__ = [
    "HarmonizedCohort",
    "load_harmonized_cohort",
    "build_luad_never_smoker_signature_cohort",
    "build_crizotinib_refractory_cohort",
    "build_alectinib_refractory_cohort",
    "classify_amino_acid_change",
    "class_forced_6class_counts_per_patient",
]

ALK_SECOND_GEN_TKIS: frozenset[str] = frozenset(
    {"alectinib", "brigatinib", "ceritinib", "lorlatinib", "ensartinib"}
)

ALK_KD_RESIDUE_RANGE: tuple[int, int] = (1116, 1392)


@dataclass
class HarmonizedCohort:
    patients: pd.DataFrame
    mutations: pd.DataFrame
    treatment_lines: pd.DataFrame
    source_dir: Path


def load_harmonized_cohort(cohort_dir: str | Path) -> HarmonizedCohort:
    cohort_dir = Path(cohort_dir)
    pats = pd.read_csv(cohort_dir / "patients.csv")
    muts = pd.read_csv(cohort_dir / "mutations.csv")
    tx = pd.read_csv(cohort_dir / "treatment_lines.csv")
    return HarmonizedCohort(
        patients=pats,
        mutations=muts,
        treatment_lines=tx,
        source_dir=cohort_dir,
    )


# ---------------------------------------------------------------------------
# Cohort filters
# ---------------------------------------------------------------------------


def build_luad_never_smoker_signature_cohort(
    cohort: HarmonizedCohort,
    *,
    require_fusion_evidence: bool = True,
    require_mutations: bool = True,
    exclude_ctdna_only_sources: bool = True,
) -> pd.DataFrame:
    """Return patients that feed the LUAD never-smoker 6-class spectrum.

    Filters:
      - `analysis_core == True`.
      - `smoking_status == "Never"` (strict; "Ever" and NaN both excluded).
      - If `require_fusion_evidence`: `alk_evidence in {ngs_fusion,
        clinical_alk_positive}`. Setting False allows alk_specific_tki_only
        (TKI-only inference), which is a looser cohort.
      - If `require_mutations`: the patient has ≥1 non-synonymous mutation
        in `mutations.csv` with `source_type in {ngs}` (so we have a real
        tumour mutational catalogue; ctdna_* is small-panel plasma and
        misses the mutagenic process we want to characterise).
    """
    p = cohort.patients
    mask = p["analysis_core"] == True
    mask &= p["smoking_status"].astype(str).str.strip().str.lower() == "never"
    if require_fusion_evidence:
        mask &= p["alk_evidence"].isin(["ngs_fusion", "clinical_alk_positive"])
    sub = p[mask].copy()

    if require_mutations:
        mut_source = cohort.mutations
        if exclude_ctdna_only_sources:
            mut_source = mut_source[mut_source["source_type"] == "ngs"]
        with_muts = set(mut_source["patient_key"].unique())
        sub = sub[sub["patient_key"].isin(with_muts)]
    return sub


def _normalize_drug_list(value) -> set[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return set()
    s = str(value).lower()
    for sep in (",", ";", ">", "|"):
        s = s.replace(sep, ",")
    return {t.strip() for t in s.split(",") if t.strip()}


def _progressed_on_tki(cohort: HarmonizedCohort, tki_name: str) -> set[str]:
    """Patients whose treatment-line record indicates they stopped the
    named TKI (any time_on_drug_event == 1), implying progression or
    discontinuation."""
    tki_lines = cohort.treatment_lines
    tki_lines = tki_lines[tki_lines["is_alk_tki"] == True]
    tki_lines = tki_lines[
        tki_lines["alk_tki_agents"]
        .astype(str)
        .str.lower()
        .str.contains(tki_name.lower(), regex=False, na=False)
    ]
    # Event = 1 → the drug stopped (progression or death). Event = 0 → censored.
    stopped = tki_lines[tki_lines["time_on_drug_event"].fillna(0) == 1]
    return set(stopped["patient_key"].unique())


def _prior_tki_set_per_patient(cohort: HarmonizedCohort) -> dict[str, set[str]]:
    """Return patient_key -> set of ALK-TKI drug names ever received (any line)."""
    tx = cohort.treatment_lines[cohort.treatment_lines["is_alk_tki"] == True]
    out: dict[str, set[str]] = {}
    for pk, agents in tx.groupby("patient_key")["alk_tki_agents"]:
        all_drugs: set[str] = set()
        for a in agents:
            all_drugs |= _normalize_drug_list(a)
        out[pk] = all_drugs
    return out


def build_crizotinib_refractory_cohort(cohort: HarmonizedCohort) -> pd.DataFrame:
    """T1 clinic cohort — patients whose only prior ALK TKI exposure is
    crizotinib AND who stopped crizotinib (time_on_drug_event = 1 in that
    line, implying progression)."""
    stopped_criz = _progressed_on_tki(cohort, "crizotinib")
    prior_tkis = _prior_tki_set_per_patient(cohort)
    pass_patients = {
        pk
        for pk in stopped_criz
        if "crizotinib" in prior_tkis.get(pk, set())
        and not (prior_tkis[pk] & ALK_SECOND_GEN_TKIS)
    }
    p = cohort.patients
    sub = p[(p["analysis_core"] == True) & p["patient_key"].isin(pass_patients)].copy()
    return sub


def build_alectinib_refractory_cohort(cohort: HarmonizedCohort) -> pd.DataFrame:
    """T3 clinic cohort — patients who stopped alectinib (progression),
    regardless of prior crizotinib exposure."""
    stopped_alec = _progressed_on_tki(cohort, "alectinib")
    p = cohort.patients
    sub = p[(p["analysis_core"] == True) & p["patient_key"].isin(stopped_alec)].copy()
    return sub


# ---------------------------------------------------------------------------
# Class-forced amino-acid-change inference → 6-class pyrimidine spectrum
# ---------------------------------------------------------------------------

from .codon_tools import GENETIC_CODE, enumerate_snv_routes  # noqa: E402

_AA_TO_CODONS: dict[str, list[str]] = {}
for _codon, _aa in GENETIC_CODE.items():
    _AA_TO_CODONS.setdefault(_aa, []).append(_codon)


def classify_amino_acid_change(
    wt_aa: str,
    mut_aa: str,
    wt_codon: str | None = None,
) -> dict:
    """Return class-forced status for an amino-acid change.

    If `wt_codon` is given, routes are enumerated only from that codon
    (deterministic answer; use for genes with CDS-locked WT codons such as
    ALK residues in `data/reference/codons.yaml`). If `wt_codon` is None,
    routes are enumerated over every synonymous WT codon encoding `wt_aa`
    and the result is class-forced only when every synonymous codon's routes
    share a single pyrimidine class — the same conservative rule used for the
    Bradeen and Azam derivations.

    Returns dict with keys:
        - reachable_1nt : bool
        - class_forced : bool
        - forced_class : str or None
        - classes : list[str] (all classes reachable; length 1 if forced)
        - n_routes_per_wt_codon : "1" / "2" / "VARIES a..b"
    """
    if wt_aa not in _AA_TO_CODONS or len(mut_aa) != 1 or mut_aa not in "ACDEFGHIKLMNPQRSTVWY*":
        return {
            "reachable_1nt": False,
            "class_forced": False,
            "forced_class": None,
            "classes": [],
            "n_routes_per_wt_codon": "0",
        }

    if wt_codon is not None:
        codons = [wt_codon]
    else:
        codons = _AA_TO_CODONS[wt_aa]

    classes: set[str] = set()
    n_routes_per_codon: list[int] = []
    for codon in codons:
        rs = enumerate_snv_routes(codon, mut_aa)
        if rs:
            n_routes_per_codon.append(len(rs))
            for r in rs:
                classes.add(r.sub_class)
    if not classes:
        return {
            "reachable_1nt": False,
            "class_forced": False,
            "forced_class": None,
            "classes": [],
            "n_routes_per_wt_codon": "0",
        }
    if min(n_routes_per_codon) == max(n_routes_per_codon):
        notes = str(n_routes_per_codon[0])
    else:
        notes = f"VARIES {min(n_routes_per_codon)}..{max(n_routes_per_codon)}"
    return {
        "reachable_1nt": True,
        "class_forced": len(classes) == 1,
        "forced_class": next(iter(classes)) if len(classes) == 1 else None,
        "classes": sorted(classes),
        "n_routes_per_wt_codon": notes,
    }


def load_alk_codons(
    codons_yaml_path: str | Path = "data/reference/codons_alk_full.yaml",
) -> dict[int, str]:
    """Return a {residue -> WT codon} dict for the ALK kinase domain.

    The default source is the full 331-codon table read directly off the Oh
    2026 Table S11 EML4-ALK variant-1 wild-type construct, which is the
    sequence the prime-editing atlas was built in; every codon in it has been
    checked to encode the wild-type residue the atlas itself reports at that
    position (scripts/build_alk_codons.py).

    Passing the older `data/reference/codons.yaml` returns only the
    hand-curated slice, which covers 23 residues. Callers that fall back to
    that file get `None` for most residues, and `classify_amino_acid_change`
    then applies its conservative rule of assigning a forced class only when
    every synonymous wild-type codon agrees. That rule is correct but
    strictly weaker: it leaves a class unassigned wherever the synonymous
    codons differ, even though the construct fixes which one is present.
    """
    import yaml

    path = Path(codons_yaml_path)
    with open(path) as f:
        raw = yaml.safe_load(f)
    if "codons" in raw and "alk_nm_004304_5" not in raw:
        return {int(r): c for r, c in raw["codons"].items()}
    alk = raw["alk_nm_004304_5"]["codons"]
    return {int(r): v["codon"] for r, v in alk.items() if "codon" in v}


def class_forced_6class_counts_per_patient(
    mutations: pd.DataFrame,
    *,
    exclude_gene: str | None = "ALK",
    variant_classes: tuple[str, ...] = ("Missense_Mutation",),
) -> pd.DataFrame:
    """Tally per-patient 6-class pyrimidine substitution counts, keeping only
    class-forced amino-acid changes (so no CDS lookup is needed per gene).

    Parameters
    ----------
    mutations : DataFrame
        Columns at least: patient_key, gene, protein_change, variant_class.
        protein_change in HGVSp_Short-ish form "p.<WT><pos><MUT>".
    exclude_gene : str
        Gene to drop (default ALK — avoid selection coloring the signature).
    variant_classes : tuple
        Variant classes to keep. Default: only missense (SNVs). Nonsense is
        also a SNV but is more likely drug-selection-correlated.

    Returns
    -------
    DataFrame indexed by patient_key with columns = 6 pyrimidine classes
    plus `total` and `n_class_forced`.
    """
    from .codon_tools import SUB_CLASSES

    m = mutations[mutations["variant_class"].isin(variant_classes)].copy()
    if exclude_gene is not None:
        m = m[m["gene"].astype(str).str.upper() != exclude_gene.upper()]

    parsed = m["protein_change"].astype(str).str.extract(
        r"p\.(?P<wt>[A-Z\*])(?P<res>\d+)(?P<mut>[A-Z\*])"
    )
    valid = parsed["mut"].notna() & parsed["wt"].notna()
    m = m[valid].copy()
    m["wt_aa"] = parsed.loc[valid, "wt"].values
    m["mut_aa"] = parsed.loc[valid, "mut"].values

    results = {}
    for pk in m["patient_key"].unique():
        results[pk] = {c: 0 for c in SUB_CLASSES}
        results[pk]["n_class_forced"] = 0
        results[pk]["n_total_parsed"] = 0

    for row in m.itertuples(index=False):
        pk = row.patient_key
        results[pk]["n_total_parsed"] += 1
        info = classify_amino_acid_change(row.wt_aa, row.mut_aa)
        if info["class_forced"]:
            results[pk][info["forced_class"]] += 1
            results[pk]["n_class_forced"] += 1

    df = pd.DataFrame.from_dict(results, orient="index")
    df.index.name = "patient_key"
    df["total"] = df[list(SUB_CLASSES)].sum(axis=1)
    return df
