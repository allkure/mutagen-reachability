"""GENIE + MSK-IMPACT ingestion and cohort filtering.

Public interface
----------------
- `load_genie_release()`: load the four core GENIE files from a release dir.
- `attach_clinical_annotation()`: left-join the extended clinical table
  (smoking status, prior TKI exposure, biopsy timing) onto GENIE samples.
- `build_alk_fusion_luad_never_smoker()`: cohort A (feeds LUAD never-smoker
  signature fit for T1 and T3).
- `build_crizotinib_refractory_cohort()`: cohort B (clinic side of T1).
- `build_alectinib_refractory_cohort()`: cohort C (clinic side of T3).
- `variant_counts_from_cohort()`: tally per-residue amino-acid variant
  counts over the ALK kinase domain for a given cohort.

Data provenance assumptions and expected schemas are documented in
`data/raw/genie_msk/README.md`. If the uploaded files diverge from these
schemas, these loaders will raise with a clear message naming the
missing column so the fix is a one-line config edit, not a code rewrite.

None of these functions downloads data; they operate on files already on
disk. The project cloud container's proxy allowlist blocks Synapse and
cBioPortal direct downloads, so MR provisions the input files.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import pandas as pd

__all__ = [
    "GenieRelease",
    "load_genie_release",
    "attach_clinical_annotation",
    "build_alk_fusion_luad_never_smoker",
    "build_crizotinib_refractory_cohort",
    "build_alectinib_refractory_cohort",
    "variant_counts_from_cohort",
    "NEVER_SMOKER_LABELS",
    "ALK_SECOND_GEN_TKIS",
    "ALK_KD_RESIDUE_RANGE",
]

# ---------------------------------------------------------------------------
# Project-wide constants used in cohort definitions.
# ---------------------------------------------------------------------------

#: Values of SMOKING_STATUS treated as "never-smoker" after case-fold.
NEVER_SMOKER_LABELS: frozenset[str] = frozenset(
    {"never", "lifetime non-smoker", "nonsmoker", "non-smoker", "never-smoker"}
)

#: Second- and third-generation ALK TKIs. A crizotinib-refractory sample is
#: defined as prior crizotinib AND NO prior exposure to these.
ALK_SECOND_GEN_TKIS: frozenset[str] = frozenset(
    {"alectinib", "brigatinib", "ceritinib", "lorlatinib", "ensartinib"}
)

#: ALK kinase-domain residue range (type Ia on NM_004304.5). From Doi 2022.
ALK_KD_RESIDUE_RANGE: tuple[int, int] = (1116, 1392)


# ---------------------------------------------------------------------------
# Container for a loaded GENIE release
# ---------------------------------------------------------------------------


@dataclass
class GenieRelease:
    """A loaded AACR GENIE pooled release.

    Attributes
    ----------
    release_dir : Path
    mutations : pandas.DataFrame
    patients : pandas.DataFrame
    samples : pandas.DataFrame
    fusions : pandas.DataFrame | None
    panels : pandas.DataFrame | None
    extended_clinical : pandas.DataFrame | None  (attached post-load)
    """

    release_dir: Path
    mutations: pd.DataFrame
    patients: pd.DataFrame
    samples: pd.DataFrame
    fusions: pd.DataFrame | None = None
    panels: pd.DataFrame | None = None
    extended_clinical: pd.DataFrame | None = None
    release_tag: str = "unknown"


# ---------------------------------------------------------------------------
# GENIE release loader
# ---------------------------------------------------------------------------


def _check_required_columns(df: pd.DataFrame, required: Iterable[str], source: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"{source}: missing required columns {missing}. "
            f"Available: {list(df.columns)[:12]}... "
            "See data/raw/genie_msk/README.md for the expected schema."
        )


def load_genie_release(
    release_dir: str | Path,
    *,
    release_tag: str = "unknown",
    include_fusions: bool = True,
    include_panels: bool = True,
) -> GenieRelease:
    """Load the four core GENIE files from `release_dir`.

    Expected files:
      - data_mutations_extended.txt
      - data_clinical_patient.txt
      - data_clinical_sample.txt
      - data_fusions.txt (or data_sv.txt in later releases)
      - genomic_information.txt (panel info)

    Any lines starting with "#" are treated as header comments and skipped.
    """
    release_dir = Path(release_dir)
    if not release_dir.exists():
        raise FileNotFoundError(f"release_dir does not exist: {release_dir}")

    mut_path = release_dir / "data_mutations_extended.txt"
    pat_path = release_dir / "data_clinical_patient.txt"
    smp_path = release_dir / "data_clinical_sample.txt"
    mutations = pd.read_csv(mut_path, sep="\t", comment="#", low_memory=False)
    patients = pd.read_csv(pat_path, sep="\t", comment="#", low_memory=False)
    samples = pd.read_csv(smp_path, sep="\t", comment="#", low_memory=False)

    _check_required_columns(
        mutations,
        [
            "Hugo_Symbol",
            "Chromosome",
            "Start_Position",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
            "Variant_Classification",
            "HGVSp_Short",
            "Tumor_Sample_Barcode",
        ],
        source=str(mut_path),
    )
    _check_required_columns(patients, ["PATIENT_ID"], source=str(pat_path))
    _check_required_columns(
        samples,
        ["SAMPLE_ID", "PATIENT_ID", "CANCER_TYPE_DETAILED", "SEQ_ASSAY_ID"],
        source=str(smp_path),
    )

    fusions = None
    if include_fusions:
        for name in ("data_fusions.txt", "data_sv.txt"):
            p = release_dir / name
            if p.exists():
                fusions = pd.read_csv(p, sep="\t", comment="#", low_memory=False)
                break

    panels = None
    if include_panels:
        p = release_dir / "genomic_information.txt"
        if p.exists():
            panels = pd.read_csv(p, sep="\t", comment="#", low_memory=False)

    return GenieRelease(
        release_dir=release_dir,
        mutations=mutations,
        patients=patients,
        samples=samples,
        fusions=fusions,
        panels=panels,
        release_tag=release_tag,
    )


# ---------------------------------------------------------------------------
# Attach the richer clinical annotation (smoking, TKI exposure)
# ---------------------------------------------------------------------------


def attach_clinical_annotation(
    release: GenieRelease, extended_clinical: pd.DataFrame
) -> GenieRelease:
    """Left-join `extended_clinical` (keyed on PATIENT_ID) onto the GENIE
    samples table. The input is expected to carry at least:
        PATIENT_ID,
        SMOKING_STATUS,
        PRIOR_ALK_TKI,
        PROGRESSED_ON_ALK_TKI,
        ALK_FUSION_PARTNER (optional; falls back to GENIE fusions table),
        BIOPSY_TIMING (optional).
    Columns absent in the input are left absent — the cohort builders
    emit descriptive errors naming the missing field when they can't run.
    """
    required = ["PATIENT_ID", "SMOKING_STATUS", "PRIOR_ALK_TKI", "PROGRESSED_ON_ALK_TKI"]
    _check_required_columns(extended_clinical, required, source="extended_clinical")
    release.extended_clinical = extended_clinical.copy()
    return release


# ---------------------------------------------------------------------------
# Cohort builders — pure dataframe operations, no network, no side effects
# ---------------------------------------------------------------------------


def _normalize_tki_list(value) -> set[str]:
    """Normalize a free-text comma-separated TKI list to a lowercase set."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return set()
    return {t.strip().lower() for t in str(value).split(",") if t.strip()}


def _patients_with_alk_fusion(release: GenieRelease) -> set[str]:
    """Return PATIENT_IDs with at least one ALK fusion event in the release.

    Uses `release.fusions` if available; falls back to
    `release.extended_clinical["ALK_FUSION_PARTNER"]` non-null if present.
    """
    patient_ids: set[str] = set()

    if release.fusions is not None:
        f = release.fusions
        # The GENIE fusion-table schema varies by release. Try several
        # column shapes, in order of frequency.
        if "Hugo_Symbol" in f.columns and "Tumor_Sample_Barcode" in f.columns:
            alk_rows = f[f["Hugo_Symbol"].astype(str).str.upper() == "ALK"]
            alk_samples = set(alk_rows["Tumor_Sample_Barcode"].astype(str))
            smp2pat = release.samples.set_index("SAMPLE_ID")["PATIENT_ID"].to_dict()
            patient_ids.update(
                smp2pat.get(s) for s in alk_samples if s in smp2pat
            )
        elif "Site1_Hugo_Symbol" in f.columns:
            # Newer SV schema; also look at site 2
            site_mask = (
                (f["Site1_Hugo_Symbol"].astype(str).str.upper() == "ALK")
                | (f.get("Site2_Hugo_Symbol", pd.Series(dtype=str)).astype(str).str.upper() == "ALK")
            )
            alk_samples = set(f.loc[site_mask, "Sample_Id"].astype(str))
            smp2pat = release.samples.set_index("SAMPLE_ID")["PATIENT_ID"].to_dict()
            patient_ids.update(
                smp2pat.get(s) for s in alk_samples if s in smp2pat
            )

    if (
        release.extended_clinical is not None
        and "ALK_FUSION_PARTNER" in release.extended_clinical.columns
    ):
        ec = release.extended_clinical
        with_fusion = ec[ec["ALK_FUSION_PARTNER"].notna() & (ec["ALK_FUSION_PARTNER"] != "")]
        patient_ids.update(with_fusion["PATIENT_ID"].astype(str))

    patient_ids.discard(None)
    return patient_ids


def _luad_sample_ids(release: GenieRelease) -> set[str]:
    s = release.samples
    luad_mask = s["CANCER_TYPE_DETAILED"].astype(str).str.strip().str.lower() == "lung adenocarcinoma"
    return set(s.loc[luad_mask, "SAMPLE_ID"].astype(str))


def build_alk_fusion_luad_never_smoker(release: GenieRelease) -> pd.DataFrame:
    """Return the samples that feed the LUAD never-smoker signature fit.

    Filter:
      - CANCER_TYPE_DETAILED == Lung Adenocarcinoma
      - Patient has an ALK fusion event
      - SMOKING_STATUS case-insensitive in NEVER_SMOKER_LABELS

    Returns one row per SAMPLE_ID with columns:
      SAMPLE_ID, PATIENT_ID, SEQ_ASSAY_ID, SMOKING_STATUS,
      CANCER_TYPE_DETAILED, ALK_FUSION_PARTNER (if available)
    """
    if release.extended_clinical is None:
        raise ValueError(
            "LUAD-never-smoker cohort requires extended_clinical with SMOKING_STATUS; "
            "call attach_clinical_annotation() first."
        )

    luad_samples = _luad_sample_ids(release)
    alk_patients = _patients_with_alk_fusion(release)
    ec = release.extended_clinical
    smoking_lower = ec["SMOKING_STATUS"].astype(str).str.strip().str.lower()
    never_smoker_patients = set(
        ec.loc[smoking_lower.isin(NEVER_SMOKER_LABELS), "PATIENT_ID"].astype(str)
    )

    s = release.samples
    mask = (
        s["SAMPLE_ID"].astype(str).isin(luad_samples)
        & s["PATIENT_ID"].astype(str).isin(alk_patients)
        & s["PATIENT_ID"].astype(str).isin(never_smoker_patients)
    )
    cols = [c for c in ["SAMPLE_ID", "PATIENT_ID", "SEQ_ASSAY_ID", "CANCER_TYPE_DETAILED"] if c in s.columns]
    cohort = s.loc[mask, cols].copy()

    # Join extended_clinical for smoking status and fusion partner.
    add = ec[[c for c in ["PATIENT_ID", "SMOKING_STATUS", "ALK_FUSION_PARTNER"] if c in ec.columns]]
    cohort = cohort.merge(add, on="PATIENT_ID", how="left")
    return cohort


def build_crizotinib_refractory_cohort(release: GenieRelease) -> pd.DataFrame:
    """T1 clinic cohort.

    Filter:
      - ALK fusion
      - CANCER_TYPE_DETAILED == Lung Adenocarcinoma
      - PRIOR_ALK_TKI includes "crizotinib" AND no other ALK TKI
      - PROGRESSED_ON_ALK_TKI includes "crizotinib"
      - BIOPSY_TIMING == "post-TKI" if available (preferred)
    """
    if release.extended_clinical is None:
        raise ValueError(
            "Crizotinib-refractory cohort requires extended_clinical with "
            "PRIOR_ALK_TKI and PROGRESSED_ON_ALK_TKI; call attach_clinical_annotation()."
        )

    alk_patients = _patients_with_alk_fusion(release)
    luad_samples = _luad_sample_ids(release)
    ec = release.extended_clinical.copy()
    ec["_prior_tki_set"] = ec["PRIOR_ALK_TKI"].map(_normalize_tki_list)
    ec["_progressed_tki_set"] = ec["PROGRESSED_ON_ALK_TKI"].map(_normalize_tki_list)
    ec["_is_crizotinib_only_prior"] = ec["_prior_tki_set"].map(
        lambda s: ("crizotinib" in s) and (not (s & ALK_SECOND_GEN_TKIS))
    )
    ec["_progressed_on_crizotinib"] = ec["_progressed_tki_set"].map(
        lambda s: "crizotinib" in s
    )
    ec_pass = ec[ec["_is_crizotinib_only_prior"] & ec["_progressed_on_crizotinib"]]
    pass_patients = set(ec_pass["PATIENT_ID"].astype(str))

    s = release.samples
    mask = (
        s["SAMPLE_ID"].astype(str).isin(luad_samples)
        & s["PATIENT_ID"].astype(str).isin(alk_patients)
        & s["PATIENT_ID"].astype(str).isin(pass_patients)
    )
    cols = [c for c in ["SAMPLE_ID", "PATIENT_ID", "SEQ_ASSAY_ID"] if c in s.columns]
    cohort = s.loc[mask, cols].copy()
    cohort = cohort.merge(
        ec[[c for c in ["PATIENT_ID", "SMOKING_STATUS", "ALK_FUSION_PARTNER", "BIOPSY_TIMING"]
             if c in ec.columns]],
        on="PATIENT_ID",
        how="left",
    )

    # Prefer post-TKI biopsies when BIOPSY_TIMING is available (otherwise
    # keep both so the sample size does not collapse; this choice is
    # pre-registered in §5 of the pre-registration).
    if "BIOPSY_TIMING" in cohort.columns:
        timing = cohort["BIOPSY_TIMING"].astype(str).str.strip().str.lower()
        post = cohort[timing == "post-tki"]
        if len(post) >= 20:
            cohort = post
    return cohort


def build_alectinib_refractory_cohort(release: GenieRelease) -> pd.DataFrame:
    """T3 clinic cohort.

    Filter:
      - ALK fusion
      - CANCER_TYPE_DETAILED == Lung Adenocarcinoma
      - PROGRESSED_ON_ALK_TKI includes "alectinib" (prior crizotinib allowed)
    """
    if release.extended_clinical is None:
        raise ValueError(
            "Alectinib-refractory cohort requires extended_clinical with "
            "PROGRESSED_ON_ALK_TKI; call attach_clinical_annotation()."
        )

    alk_patients = _patients_with_alk_fusion(release)
    luad_samples = _luad_sample_ids(release)
    ec = release.extended_clinical.copy()
    ec["_progressed_tki_set"] = ec["PROGRESSED_ON_ALK_TKI"].map(_normalize_tki_list)
    ec_pass = ec[ec["_progressed_tki_set"].map(lambda s: "alectinib" in s)]
    pass_patients = set(ec_pass["PATIENT_ID"].astype(str))

    s = release.samples
    mask = (
        s["SAMPLE_ID"].astype(str).isin(luad_samples)
        & s["PATIENT_ID"].astype(str).isin(alk_patients)
        & s["PATIENT_ID"].astype(str).isin(pass_patients)
    )
    cols = [c for c in ["SAMPLE_ID", "PATIENT_ID", "SEQ_ASSAY_ID"] if c in s.columns]
    cohort = s.loc[mask, cols].copy()
    cohort = cohort.merge(
        ec[[c for c in ["PATIENT_ID", "SMOKING_STATUS", "ALK_FUSION_PARTNER", "BIOPSY_TIMING"]
             if c in ec.columns]],
        on="PATIENT_ID",
        how="left",
    )
    return cohort


# ---------------------------------------------------------------------------
# Variant-count tally over the kinase domain
# ---------------------------------------------------------------------------


def variant_counts_from_cohort(
    release: GenieRelease,
    cohort: pd.DataFrame,
    *,
    gene: str = "ALK",
    residue_range: tuple[int, int] = ALK_KD_RESIDUE_RANGE,
) -> pd.DataFrame:
    """Tally per-amino-acid-variant counts for the kinase-domain residues
    reported in `release.mutations` for the samples listed in `cohort`.

    Returns a DataFrame with columns:
      variant_ia  (e.g., "L1196M")
      wt_aa
      residue
      mut_aa
      n_isolates  (count across cohort samples)
      p_cohort    (count / total)

    One mutation per sample per residue; duplicate SAMPLE_ID×HGVSp_Short rows
    in the mutations table are deduplicated. Missense-and-indel inclusive;
    callers should filter by Variant_Classification upstream if needed.
    """
    cohort_samples = set(cohort["SAMPLE_ID"].astype(str))
    m = release.mutations
    m = m[m["Hugo_Symbol"].astype(str).str.upper() == gene.upper()]
    m = m[m["Tumor_Sample_Barcode"].astype(str).isin(cohort_samples)]
    m = m[m["Variant_Classification"].astype(str).isin({
        "Missense_Mutation",
        "Nonsense_Mutation",
        "In_Frame_Del",
        "In_Frame_Ins",
    })]
    # Parse residue number from HGVSp_Short like "p.L1196M" or "p.F1174L".
    hgvsp = m["HGVSp_Short"].astype(str).str.strip()
    parsed = hgvsp.str.extract(r"p\.(?P<wt>[A-Z\*])(?P<res>\d+)(?P<mut>[A-Z\*]|=|del|ins|fs)")
    parsed["res"] = pd.to_numeric(parsed["res"], errors="coerce")
    lo, hi = residue_range
    kd_mask = parsed["res"].between(lo, hi, inclusive="both")
    m = m.loc[kd_mask].copy()
    m["wt_aa"] = parsed.loc[kd_mask, "wt"].values
    m["residue"] = parsed.loc[kd_mask, "res"].astype(int).values
    m["mut_aa"] = parsed.loc[kd_mask, "mut"].values
    m["variant_ia"] = m["wt_aa"] + m["residue"].astype(str) + m["mut_aa"]

    # Deduplicate at (sample, variant_ia) level.
    m = m.drop_duplicates(subset=["Tumor_Sample_Barcode", "variant_ia"])

    tally = (
        m.groupby(["variant_ia", "wt_aa", "residue", "mut_aa"], as_index=False)
        .size()
        .rename(columns={"size": "n_isolates"})
    )
    total = int(tally["n_isolates"].sum())
    tally["p_cohort"] = tally["n_isolates"] / total if total > 0 else 0.0
    tally = tally.sort_values("n_isolates", ascending=False).reset_index(drop=True)
    tally.attrs["cohort_samples"] = len(cohort_samples)
    tally.attrs["total_isolates"] = total
    return tally
