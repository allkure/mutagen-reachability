# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
# ---

# %% [markdown]
# # T2 preliminary supply/observation table
#
# Builds the first-pass per-variant table for the T2 ratio test:
#
#   P_ENU(v)   = fraction of Bradeen 2006 imatinib isolates (2+4 μM pooled) at variant v
#   P_XL1R(v)  = fraction of Azam 2003 isolates (5 μM STI-571, Table 1) at variant v
#   S_ENU(v)   = sum over 1-nt routes from the WT codon encoding v's residue to a
#                codon encoding v's mutant AA, of the ENU spectrum class frequencies
#   S_XL1R(v)  = same under the XL-1 Red spectrum
#
# For class-forced variants (18/20 Bradeen; 50/76 Azam Table 1 — see
# derive_substitution_classes.py in each), the supply weight equals the mutagen's
# frequency at the forced class times the number of accessible routes under the
# canonical WT codon.
#
# No statistical conclusions are drawn here. This script is Phase 1 scaffolding;
# the ratio test machinery (bootstrap, permutation null, Pearson) will run in
# `02_t2_ratio_test.py` once (a) Rasila 2009 XL-1 Red primary spectrum is
# incorporated and (b) the preregistration is signed off for execution.

# %%
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

NOTEBOOK_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = NOTEBOOK_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))

from src.codon_tools import enumerate_snv_routes  # noqa: E402
from src.spectrum import load_spectra  # noqa: E402

# %%
# 1. Load the mutagen spectra.
spectra = load_spectra(PROJECT_ROOT / "data" / "reference" / "spectra.yaml")
enu = spectra["ENU_Arnold2012_nominal"]
xl1r = spectra["XL1Red_Rasila2009"]
xl1r_sens = spectra["XL1Red_BadranLiu2015_nominal"]
print(f"ENU spectrum ({enu.name}):")
for cls, freq in enu.values.items():
    print(f"  {cls}: {freq:.3f}")
print(f"XL-1R spectrum ({xl1r.name}):")
for cls, freq in xl1r.values.items():
    print(f"  {cls}: {freq:.3f}")

# %%
# 2. Load the Bradeen 2006 imatinib arm counts and compute 2+4 μM pooled fractions.
bradeen = pd.read_csv(
    PROJECT_ROOT / "data" / "raw" / "bradeen_2006" / "table1_imatinib_counts.csv",
    comment="#",
)
bradeen["n_lowstringency"] = bradeen["counts_2uM"] + bradeen["counts_4uM"]
bradeen_total = bradeen["n_lowstringency"].sum()
bradeen["p_bradeen"] = bradeen["n_lowstringency"] / bradeen_total
print(f"Bradeen 2+4 μM pooled total: {bradeen_total} isolates across {len(bradeen)} variants")

# %%
# 3. Load the Azam 2003 counts (Table 1, kinase domain), compound-rows handled below.
azam = pd.read_csv(
    PROJECT_ROOT / "data" / "raw" / "azam_2003" / "tables_1_2_counts.csv",
    comment="#",
)
# T2 primary pairing restricts to Table 1 (kinase domain) at 5 μM.
azam_kd_5um = azam[(azam["table"] == 1) & (azam["counts_5uM"] > 0)].copy()
print(f"Azam Table 1, counts_5uM > 0: {len(azam_kd_5um)} rows (some are compound sub-variants)")

# %%
# 4. Build a shared-variant table using the class-forced derivation output.
bradeen_classes = pd.read_csv(
    PROJECT_ROOT / "data" / "raw" / "bradeen_2006" / "bradeen_imatinib_variant_routes.csv"
)
azam_classes = pd.read_csv(
    PROJECT_ROOT / "data" / "raw" / "azam_2003" / "azam_variant_routes.csv"
)

# %%
# 5. Expand Azam compound rows (e.g., "Y253H/C/F") into per-sub-variant counts.
# Rule (preregistration §5): equal split across reachable sub-variants only.
# Unreachable sub-variants (V289S, T315G, E499I — those needing >1 nt change
# from every synonymous WT codon) are DROPPED from the expansion, with their
# share redistributed to the remaining reachable sub-variants. This keeps
# the Azam row's total count intact for the T2 statistic.
#
# Rationale: assigning counts to a sub-variant that cannot arise under the
# single-nt-SNV supply model would introduce model-inconsistent noise into
# the ratio test. The all-dominant-isoform alternative is retained for a
# §7 secondary analysis.

azam_classes_set = set(azam_classes["variant_ia"].tolist())


def expand_azam_compound_row(row):
    """Return list of (sub_variant_ia, assigned_count_5uM) for one Azam row."""
    v = str(row["variant_ia"])
    n = row["counts_5uM"]
    if n == 0 or pd.isna(n):
        return []
    if "/" not in v:
        return [(v, n)]
    # Compound: e.g. "Y253H/C/F" -> ["Y253H", "Y253C", "Y253F"]
    stem = v[:-1].rstrip("/").rstrip("ACDEFGHIKLMNPQRSTVWY")  # strip trailing mut AA
    # Easier: find the position where mut AA letters begin. Numbers end the stem.
    # Variant format: <WT><num><Mut1/Mut2/...>. Find last digit.
    i = len(v) - 1
    while i >= 0 and not v[i].isdigit():
        i -= 1
    stem = v[: i + 1]  # includes last digit
    subs = v[i + 1 :].split("/")
    sub_variants = [stem + s for s in subs]
    # Keep only sub-variants that are reachable in one nt (present in azam_classes_set)
    reachable = [s for s in sub_variants if s in azam_classes_set]
    if not reachable:
        return []
    per_sub = n / len(reachable)
    return [(s, per_sub) for s in reachable]


expanded_rows = []
for _, row in azam[azam["table"] == 1].iterrows():
    for sub, count in expand_azam_compound_row(row):
        expanded_rows.append({"variant_ia": sub, "counts_5uM_expanded": count})
azam_expanded_5um = (
    pd.DataFrame(expanded_rows)
    .groupby("variant_ia", as_index=False)["counts_5uM_expanded"]
    .sum()
)
azam_expanded_total = azam_expanded_5um["counts_5uM_expanded"].sum()
azam_expanded_5um["p_azam"] = (
    azam_expanded_5um["counts_5uM_expanded"] / azam_expanded_total
)
print(
    f"Azam Table 1 after equal-split compound-row expansion at 5 μM: "
    f"{azam_expanded_total:.1f} isolates across {len(azam_expanded_5um)} sub-variants"
)

# %%
# 6. Merge class-forced status (Bradeen) with Bradeen counts, and compute
# ENU supply weight for class-forced variants.
bradeen_merged = (
    bradeen.rename(columns={"variant_ia": "variant"})
    .merge(
        bradeen_classes.rename(
            columns={
                "variant_ia": "variant",
                "forced_class": "forced_class_bradeen",
                "class_forced_by_aa_change": "class_forced_bradeen",
                "candidate_wt_codons": "candidate_codons",
                "n_routes_per_wt_codon": "n_routes_per_codon_bradeen",
            }
        ),
        on="variant",
        how="left",
    )
    .copy()
)
# Coerce the per-codon route count (stored as string: "1" or "VARIES ...").
# For class-forced variants it is always "1" by construction of our variants.
bradeen_merged["n_routes_per_codon_bradeen_int"] = pd.to_numeric(
    bradeen_merged["n_routes_per_codon_bradeen"], errors="coerce"
)


def enu_supply_from_class_forced(row):
    """For class-forced variants under any synonymous WT codon:
        S_ENU(v) = n_routes_per_wt_codon * p_ENU(forced_class).
    This is the correct supply weight: the sum over ONE WT codon's routes of
    the class frequency, where every route lies in the forced class.
    For class-ambiguous variants, this returns None.
    """
    if row["class_forced_bradeen"] != "yes":
        return None
    n = row["n_routes_per_codon_bradeen_int"]
    if pd.isna(n):
        return None
    return n * enu[row["forced_class_bradeen"]]


def xl1r_supply_from_class_forced(row, spectrum=xl1r):
    if row["class_forced_bradeen"] != "yes":
        return None
    n = row["n_routes_per_codon_bradeen_int"]
    if pd.isna(n):
        return None
    return n * spectrum[row["forced_class_bradeen"]]


bradeen_merged["S_ENU"] = bradeen_merged.apply(enu_supply_from_class_forced, axis=1)
bradeen_merged["S_XL1R_Rasila"] = bradeen_merged.apply(
    xl1r_supply_from_class_forced, axis=1
)
bradeen_merged["S_XL1R_BadranLiu"] = bradeen_merged.apply(
    lambda r: xl1r_supply_from_class_forced(r, spectrum=xl1r_sens), axis=1
)
# Alias for backward-compatible column name used in cols list below.
bradeen_merged["S_XL1R"] = bradeen_merged["S_XL1R_Rasila"]

# %%
# 7. Attach Azam fractions (from the equal-split expanded table).
bradeen_merged = bradeen_merged.merge(
    azam_expanded_5um.rename(columns={"variant_ia": "variant"}),
    on="variant",
    how="left",
)

# %%
# 8. Compact view
cols = [
    "variant",
    "n_lowstringency",
    "p_bradeen",
    "counts_5uM_expanded",
    "p_azam",
    "class_forced_bradeen",
    "forced_class_bradeen",
    "n_routes_per_codon_bradeen_int",
    "S_ENU",
    "S_XL1R_Rasila",
    "S_XL1R_BadranLiu",
]
view = bradeen_merged[cols].sort_values("n_lowstringency", ascending=False)
pd.options.display.max_rows = 50
pd.options.display.float_format = "{:.4f}".format
print(view.to_string(index=False))

# %%
# 9. Write the preliminary table for inspection and for downstream statistics.
out = PROJECT_ROOT / "data" / "derived" / "t2_preliminary_supply_table.csv"
out.parent.mkdir(parents=True, exist_ok=True)
view.to_csv(out, index=False)
print(f"\nWrote {out}")

# %%
# 10. Sanity note: for the primary ratio test, we take only variants that are
# (a) class-forced in Bradeen, (b) observed in both Bradeen (n>=1) and Azam
# single-row at 5μM (n>=1), (c) not among the Azam-unreachable set
# {V289S, T315G, E499I} (already excluded by (a) if Bradeen doesn't observe them).
#
# This script does NOT run the ratio test yet. Statistics will run only after
# Rasila 2009 primary XL-1R numbers are substituted in spectra.yaml and the
# preregistration execution gate (§6) is honored.
primary_set = view[
    (view["class_forced_bradeen"] == "yes")
    & (view["n_lowstringency"] >= 1)
    & (view["counts_5uM_expanded"] >= 1)
].copy()
print(
    f"\nVariants in T2 primary candidate set: {len(primary_set)} "
    f"(Bradeen n_lowstringency ≥ 1, Azam counts_5uM_expanded ≥ 1, class-forced)"
)
print(primary_set[["variant", "n_lowstringency", "counts_5uM_expanded", "forced_class_bradeen"]].to_string(index=False))
