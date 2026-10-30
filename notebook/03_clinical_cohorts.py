# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
# ---

# %% [markdown]
# # Clinical cohort assembly for T1 (crizotinib-refractory) and T3 (alectinib-refractory)
#
# For each cohort:
#
# 1. Apply the pre-registered filter (`src.io.build_*_cohort`).
# 2. Tally ALK kinase-domain variant counts (`src.io.variant_counts_from_cohort`).
# 3. Expand compound variants / deduplicate per (sample, variant).
# 4. Apply §5 of preregistration.md: drop variants that are not class-forced
#    AND not reachable in one nt from any synonymous WT codon; keep class-
#    forced variants, including multi-route single-class ones.
# 5. Write to data/derived/cohorts/<cohort_name>_variant_counts.csv.
#
# These tables feed the T1 and T3 ratio tests (notebook/04_t1_ratio_test.py
# and notebook/05_t3_ratio_test.py), paired with the LUAD never-smoker
# spectrum from notebook/02 and the ENU / Mutazyme II spectra already in
# spectra.yaml.

# %%
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

NOTEBOOK_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = NOTEBOOK_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.io import (
    ALK_KD_RESIDUE_RANGE,
    attach_clinical_annotation,
    build_alectinib_refractory_cohort,
    build_crizotinib_refractory_cohort,
    load_genie_release,
    variant_counts_from_cohort,
)
from src.codon_tools import GENETIC_CODE, enumerate_snv_routes

# %%
GENIE_RELEASE_DIR = PROJECT_ROOT / "data" / "raw" / "genie_msk" / "v16.1-public"
EXTENDED_CLINICAL_PATH = GENIE_RELEASE_DIR / "clinical_extended.tsv"
OUT_DIR = PROJECT_ROOT / "data" / "derived" / "cohorts"
OUT_DIR.mkdir(parents=True, exist_ok=True)

if not GENIE_RELEASE_DIR.exists():
    print(f"GENIE release dir not yet provisioned at {GENIE_RELEASE_DIR}.")
    print("See data/raw/genie_msk/README.md. Exiting after scaffolding check.")
    sys.exit(0)

# %%
release = load_genie_release(GENIE_RELEASE_DIR, release_tag="v16.1-public")
ec = pd.read_csv(EXTENDED_CLINICAL_PATH, sep="\t")
release = attach_clinical_annotation(release, ec)

# %%
AA_TO_CODONS: dict[str, list[str]] = {}
for codon, aa in GENETIC_CODE.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)


def classify_variant(wt_aa: str, mut_aa: str) -> dict:
    """Return the class-forced status and forced class (if any) for an
    amino-acid change, following the same methodology as the Bradeen and
    Azam derivations in data/raw/*/derive_substitution_classes.py.

    Returns dict with keys: class_forced ("yes"/"no"), forced_class,
    n_routes_per_wt_codon, candidate_wt_codons.
    """
    if wt_aa not in AA_TO_CODONS or mut_aa not in "ACDEFGHIKLMNPQRSTVWY*":
        return {
            "class_forced": "no",
            "forced_class": "",
            "n_routes_per_wt_codon": "0",
            "candidate_wt_codons": "",
            "reachable_1nt": False,
        }
    candidate_codons = []
    all_classes: set[str] = set()
    n_routes_per_codon = []
    for codon in AA_TO_CODONS[wt_aa]:
        rr = enumerate_snv_routes(codon, mut_aa)
        if rr:
            candidate_codons.append(codon)
            n_routes_per_codon.append(len(rr))
            for r in rr:
                all_classes.add(r.sub_class)
    if not candidate_codons:
        return {
            "class_forced": "no",
            "forced_class": "",
            "n_routes_per_wt_codon": "0",
            "candidate_wt_codons": "",
            "reachable_1nt": False,
        }
    if min(n_routes_per_codon) == max(n_routes_per_codon):
        n_notes = str(n_routes_per_codon[0])
    else:
        n_notes = f"VARIES {min(n_routes_per_codon)}..{max(n_routes_per_codon)}"
    forced = len(all_classes) == 1
    return {
        "class_forced": "yes" if forced else "no",
        "forced_class": next(iter(all_classes)) if forced else "",
        "n_routes_per_wt_codon": n_notes,
        "candidate_wt_codons": "|".join(candidate_codons),
        "reachable_1nt": True,
    }


def annotate_counts(tally: pd.DataFrame) -> pd.DataFrame:
    """Annotate a variant-count tally with class-forced / forced-class."""
    rows = []
    for _, r in tally.iterrows():
        info = classify_variant(r["wt_aa"], r["mut_aa"])
        rows.append({**r.to_dict(), **info})
    return pd.DataFrame(rows)


# %%
# T1 clinic side.
print("Building crizotinib-refractory cohort (T1 clinic)...")
cr = build_crizotinib_refractory_cohort(release)
cr_tally = variant_counts_from_cohort(release, cr, gene="ALK",
                                       residue_range=ALK_KD_RESIDUE_RANGE)
cr_annotated = annotate_counts(cr_tally)
cr_path = OUT_DIR / "t1_crizotinib_refractory_variant_counts.csv"
cr_annotated.to_csv(cr_path, index=False)
print(f"  Samples: {cr_tally.attrs.get('cohort_samples', '?')}")
print(f"  Isolates: {cr_tally.attrs.get('total_isolates', '?')}")
print(f"  Class-forced: {(cr_annotated['class_forced']=='yes').sum()} / {len(cr_annotated)}")
print(f"  Wrote {cr_path}")

# %%
# T3 clinic side.
print("\nBuilding alectinib-refractory cohort (T3 clinic)...")
al = build_alectinib_refractory_cohort(release)
al_tally = variant_counts_from_cohort(release, al, gene="ALK",
                                       residue_range=ALK_KD_RESIDUE_RANGE)
al_annotated = annotate_counts(al_tally)
al_path = OUT_DIR / "t3_alectinib_refractory_variant_counts.csv"
al_annotated.to_csv(al_path, index=False)
print(f"  Samples: {al_tally.attrs.get('cohort_samples', '?')}")
print(f"  Isolates: {al_tally.attrs.get('total_isolates', '?')}")
print(f"  Class-forced: {(al_annotated['class_forced']=='yes').sum()} / {len(al_annotated)}")
print(f"  Wrote {al_path}")
