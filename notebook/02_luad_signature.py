# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
# ---

# %% [markdown]
# # LUAD never-smoker signature fit (T1 and T3 mutagen side)
#
# Procedure:
#
# 1. Load the AACR GENIE pooled release (data_mutations_extended + clinical).
# 2. Attach the extended clinical annotation (smoking status, prior TKI).
# 3. Build the LUAD-never-smoker-with-ALK-fusion cohort.
# 4. For each patient, assemble their 96-channel SBS vector over the panel-
#    covered coding region. Uses trinucleotide context from GRCh38 chr2/chr9.
# 5. Fit COSMIC v3.3 SBS signatures to each patient's SBS96 vector by NNLS.
# 6. Collapse the fitted mix (× COSMIC signatures) back to the patient's
#    reconstructed 96-channel then to the 6-class pyrimidine spectrum.
# 7. Compute the cohort-level 6-class spectrum with mutation-count weighted
#    mean over patients and 10,000-bootstrap CI.
# 8. Write to data/derived/signatures/luad_neversmoker_alkfusion.yaml and
#    data/derived/signatures/luad_neversmoker_alkfusion_bootstrap.npz.
#
# This script is parameterized so re-running under a different cohort filter
# (e.g., never-smoker without ALK-fusion restriction, as a sensitivity check)
# is a two-line change.
#
# **Data provenance locked in `data/raw/genie_msk/README.md`.** Running this
# script requires:
#   - GENIE release files at data/raw/genie_msk/<release-tag>/
#   - COSMIC v3.3 SBS matrix at data/raw/cosmic_v33/COSMIC_v3.3.1_SBS_GRCh38.txt
#   - Extended clinical annotation (TSV) at
#     data/raw/genie_msk/<release-tag>/clinical_extended.tsv
#   - GRCh38 FASTA (or pyfaidx-indexed) with chr2 and chr9 covered.

# %%
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

NOTEBOOK_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = NOTEBOOK_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.io import (
    attach_clinical_annotation,
    build_alk_fusion_luad_never_smoker,
    load_genie_release,
)
from src.signatures import (
    SBS96_CHANNELS,
    bootstrap_cohort_spectrum,
    collapse_sbs96_to_6class,
    fit_sbs96_patient,
    load_cosmic_v33_sbs,
)

# %%
# Parameters. Edit for sensitivity analyses.
GENIE_RELEASE_DIR = PROJECT_ROOT / "data" / "raw" / "genie_msk" / "v16.1-public"
EXTENDED_CLINICAL_PATH = GENIE_RELEASE_DIR / "clinical_extended.tsv"
COSMIC_PATH = PROJECT_ROOT / "data" / "raw" / "cosmic_v33" / "COSMIC_v3.3.1_SBS_GRCh38.txt"
GRCH38_FASTA = PROJECT_ROOT / "data" / "raw" / "grch38" / "GRCh38.primary_assembly.chr2_chr9.fa"
OUT_YAML = PROJECT_ROOT / "data" / "derived" / "signatures" / "luad_neversmoker_alkfusion.yaml"
OUT_BOOT = PROJECT_ROOT / "data" / "derived" / "signatures" / "luad_neversmoker_alkfusion_bootstrap.npz"

# Which COSMIC signatures to allow in the fit. Setting this to None uses
# the full COSMIC v3.3 catalogue (recommended for primary; the never-smoker
# profile emerges from the NNLS, not from asserting SBS1+SBS5+SBS40a a priori).
# A restricted sensitivity run can set this to a small list.
RESTRICT_SIGNATURES = None  # or ["SBS1", "SBS5", "SBS40a"]

# %%
# 1. Load inputs. Skipped if files are missing — scaffolds the pipeline.
if not GENIE_RELEASE_DIR.exists():
    print(f"GENIE release dir not yet provisioned at {GENIE_RELEASE_DIR}.")
    print("See data/raw/genie_msk/README.md. Exiting after scaffolding check.")
    sys.exit(0)

release = load_genie_release(GENIE_RELEASE_DIR, release_tag="v16.1-public")
ec = pd.read_csv(EXTENDED_CLINICAL_PATH, sep="\t")
release = attach_clinical_annotation(release, ec)

cosmic = load_cosmic_v33_sbs(COSMIC_PATH)
print(f"COSMIC v3.3: {cosmic.shape[1]} signatures across {cosmic.shape[0]} channels")

# %%
# 2. Build the cohort.
cohort = build_alk_fusion_luad_never_smoker(release)
print(f"LUAD never-smoker + ALK-fusion cohort: {len(cohort)} samples, "
      f"{cohort['PATIENT_ID'].nunique()} patients")

# %%
# 3. Build per-patient SBS96 counts.
#
# We filter mutations to:
#   - samples in the cohort
#   - SNVs (len(ref)=len(alt)=1, both ACGT)
#   - coding regions covered by the sample's panel (optional denominator
#     correction via release.panels — in V1 we use all SNVs genome-wide
#     in the panel, which is approximately coding).
#
# Trinucleotide context is annotated by looking up the 5' and 3' bases
# from the GRCh38 FASTA at Start_Position - 1 and Start_Position + 1
# (for a SNV on chromosome CHROMOSOME).
try:
    from pyfaidx import Fasta
except ImportError:
    print("pyfaidx not installed. Run: pip install pyfaidx")
    sys.exit(1)

fasta = Fasta(str(GRCH38_FASTA), sequence_always_upper=True)
mut = release.mutations
cohort_samples = set(cohort["SAMPLE_ID"].astype(str))
snv_mask = (
    mut["Tumor_Sample_Barcode"].astype(str).isin(cohort_samples)
    & (mut["Reference_Allele"].astype(str).str.len() == 1)
    & (mut["Tumor_Seq_Allele2"].astype(str).str.len() == 1)
    & mut["Reference_Allele"].isin(list("ACGT"))
    & mut["Tumor_Seq_Allele2"].isin(list("ACGT"))
)
snv = mut.loc[snv_mask, [
    "Chromosome", "Start_Position", "Reference_Allele",
    "Tumor_Seq_Allele2", "Tumor_Sample_Barcode"
]].copy()


def _channel_for_row(row) -> str | None:
    chrom = str(row["Chromosome"])
    pos = int(row["Start_Position"])
    ref = row["Reference_Allele"]
    alt = row["Tumor_Seq_Allele2"]
    chrom_key = chrom if chrom.startswith("chr") else f"chr{chrom}"
    if chrom_key not in fasta:
        if chrom in fasta:
            chrom_key = chrom
        else:
            return None
    try:
        five_p = str(fasta[chrom_key][pos - 2]).upper()
        three_p = str(fasta[chrom_key][pos]).upper()
    except Exception:
        return None
    if five_p not in "ACGT" or three_p not in "ACGT":
        return None
    # Pyrimidine convention: flip to the pyrimidine strand if ref is A or G.
    if ref in "AG":
        from src.codon_tools import _COMP as COMP  # type: ignore
        ref = COMP[ref]
        alt = COMP[alt]
        five_p, three_p = COMP[three_p], COMP[five_p]
    return f"{five_p}[{ref}>{alt}]{three_p}"


snv["channel"] = snv.apply(_channel_for_row, axis=1)
snv = snv.dropna(subset=["channel"])
patient_sbs96 = (
    snv.groupby(["Tumor_Sample_Barcode", "channel"]).size().unstack(fill_value=0).T
)
patient_sbs96 = patient_sbs96.reindex(SBS96_CHANNELS, fill_value=0)
print(f"Per-patient SBS96 catalogue: {patient_sbs96.shape[1]} patients × 96 channels")
print(f"Total mutations across cohort: {int(patient_sbs96.values.sum())}")

# %%
# 4. Fit COSMIC signatures per patient.
fits = {
    sid: fit_sbs96_patient(
        patient_sbs96[sid],
        cosmic,
        sample_id=sid,
        restrict_to_signatures=RESTRICT_SIGNATURES,
    )
    for sid in patient_sbs96.columns
}

# Reconstructed 96-channel vectors (NNLS-smoothed; recommended input to the
# 6-class collapse because it dampens per-patient sampling noise).
reconstructed = pd.DataFrame(
    {sid: fits[sid].reconstructed_sbs96 for sid in fits},
).reindex(SBS96_CHANNELS, fill_value=0.0)
print(f"Reconstructed per-patient SBS96 shape: {reconstructed.shape}")

# %%
# 5. Cohort-level 6-class spectrum via bootstrap over patients.
result = bootstrap_cohort_spectrum(
    reconstructed,
    n_bootstrap=10_000,
    weight_by_mutation_count=True,
)

point = result["point_estimate"]
ci = result["ci95"]
print("\nLUAD never-smoker ALK-fusion 6-class spectrum (COSMIC-reconstructed):")
for cls in point.index:
    print(f"  {cls}: {point[cls]:.4f}  [95% CI {ci.loc[cls, 'lo']:.4f}, {ci.loc[cls, 'hi']:.4f}]")
print(f"n patients: {result['n_patients']}")

# %%
# 6. Write to spectra.yaml-compatible YAML and save bootstrap samples.
OUT_YAML.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_YAML, "w") as f:
    yaml.safe_dump(
        {
            "LUAD_neversmoker_alkfusion_fitted": {
                "values": {cls: float(point[cls]) for cls in point.index},
                "ci95": {
                    cls: [float(ci.loc[cls, "lo"]), float(ci.loc[cls, "hi"])]
                    for cls in point.index
                },
                "source": (
                    f"NNLS fit of COSMIC v3.3 SBS signatures to AACR GENIE "
                    f"{release.release_tag} LUAD-never-smoker-with-ALK-fusion "
                    f"patient SBS96 catalogues; cohort-level mutation-count-"
                    f"weighted mean with 10,000 bootstrap over patients."
                ),
                "method": "signatures.fit_sbs96_patient + bootstrap_cohort_spectrum",
                "n_patients": int(result["n_patients"]),
                "restrict_signatures": RESTRICT_SIGNATURES or "none (full COSMIC v3.3)",
            }
        },
        f,
        sort_keys=False,
    )
np.savez(OUT_BOOT, bootstrap_samples=result["bootstrap_samples"])
print(f"\nWrote {OUT_YAML}")
print(f"Wrote {OUT_BOOT}")
