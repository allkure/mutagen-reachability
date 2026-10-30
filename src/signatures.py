"""Mutational-signature fitting (COSMIC v3.3) and 6-class spectrum collapse.

Public interface
----------------
- `SBS96_CHANNELS`: canonical 96-channel ordering (pyrimidine-referenced,
  trinucleotide context) used throughout the project.
- `load_cosmic_v33_sbs()`: load the COSMIC v3.3 signature matrix from disk.
- `fit_sbs96_patient()`: solve a non-negative, sum-constrained decomposition
  of a patient's 96-channel SBS vector into a mixture of COSMIC signatures.
- `collapse_sbs96_to_6class()`: aggregate a 96-channel vector to the 6-class
  pyrimidine-referenced spectrum used by `src.supply_weights.supply_weight`.
- `bootstrap_cohort_spectrum()`: cohort-level 6-class spectrum via
  per-patient fit then weighted-mean + bootstrap over patients.

The fit uses `scipy.optimize.nnls` (standard non-negative least squares),
not an NMF — the signature set is fixed to COSMIC v3.3, so this is a
constrained *regression*, not a factorization.

Random seed = `src.stats.SEED` (= 20261013) for reproducibility of all
bootstrap operations.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from .codon_tools import SUB_CLASSES
from .stats import SEED

__all__ = [
    "SBS96_CHANNELS",
    "SBS6_FROM_SBS96",
    "load_cosmic_v33_sbs",
    "fit_sbs96_patient",
    "collapse_sbs96_to_6class",
    "bootstrap_cohort_spectrum",
    "SignatureFit",
]


# ---------------------------------------------------------------------------
# SBS96 channel definition (COSMIC convention).
#
# A channel is "N[X>Y]M" where X>Y is a pyrimidine-referenced substitution
# and N, M are the 5' and 3' flanking bases on the pyrimidine strand.
# Canonical ordering: X>Y outer (C>A, C>G, C>T, T>A, T>C, T>G), then N inner
# (A, C, G, T) and M inner (A, C, G, T), yielding 6 * 16 = 96 channels.
# ---------------------------------------------------------------------------


def _build_sbs96_channels() -> list[str]:
    channels: list[str] = []
    for sub in SUB_CLASSES:
        ref, _, alt = sub.partition(">")
        for five_p in "ACGT":
            for three_p in "ACGT":
                channels.append(f"{five_p}[{ref}>{alt}]{three_p}")
    return channels


SBS96_CHANNELS: tuple[str, ...] = tuple(_build_sbs96_channels())
assert len(SBS96_CHANNELS) == 96


# Mapping each 96-channel label to its 6-class pyrimidine substitution class.
SBS6_FROM_SBS96: dict[str, str] = {}
for ch in SBS96_CHANNELS:
    SBS6_FROM_SBS96[ch] = ch[2:5]  # e.g., "A[C>T]G" -> "C>T"


# ---------------------------------------------------------------------------
# COSMIC v3.3 signature matrix loader
# ---------------------------------------------------------------------------


def load_cosmic_v33_sbs(path: str | Path) -> pd.DataFrame:
    """Load the COSMIC v3.3 SBS signature matrix from the Sanger TSV.

    The expected file has one row per 96 channel and one column per
    signature (SBS1, SBS2, SBS3, ..., SBS96 and the "a/b" variants).
    The first column is the 96-channel label in the "A[C>T]G" convention.
    Values are per-signature frequencies summing to ~1 per column.

    Parameters
    ----------
    path : str | Path
        Local filesystem path to `COSMIC_v3.3.1_SBS_GRCh38.txt` (or compatible).

    Returns
    -------
    pandas.DataFrame
        Rows indexed by SBS96 channel, columns are signature IDs.
        Verified to have 96 rows matching the expected channel order.
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t")

    # First column holds the channel label under names like "Type" or
    # "Substitution Type" or "MutationType"; detect robustly.
    channel_col = df.columns[0]
    df = df.set_index(channel_col)

    # Normalise channel labels to the canonical "N[X>Y]M" shape.
    # (COSMIC uses "A[C>T]G" format already, but defensive just in case.)
    missing = [c for c in SBS96_CHANNELS if c not in df.index]
    if missing:
        raise ValueError(
            f"{path}: COSMIC signature file missing {len(missing)} expected "
            f"channels; first missing = {missing[:3]}"
        )

    df = df.loc[list(SBS96_CHANNELS)]
    return df


# ---------------------------------------------------------------------------
# Non-negative decomposition of a patient catalogue onto COSMIC signatures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SignatureFit:
    """Result of fitting a patient's SBS96 catalogue to a COSMIC signature set.

    Attributes
    ----------
    sample_id : str
    mix : pandas.Series
        Signature contributions (non-negative, summing to approximately the
        patient's total mutation count); index = signature IDs.
    reconstructed_sbs96 : pandas.Series
        The 96-channel vector reconstructed from `mix @ cosmic_signatures.T`.
    residual : float
        L2 residual between observed and reconstructed catalogues.
    total_mutations : int
    """

    sample_id: str
    mix: pd.Series
    reconstructed_sbs96: pd.Series
    residual: float
    total_mutations: int

    @property
    def normalized_mix(self) -> pd.Series:
        """Signature contributions as fractions summing to 1."""
        total = self.mix.sum()
        return self.mix / total if total > 0 else self.mix


def fit_sbs96_patient(
    sbs96_counts: pd.Series,
    cosmic_signatures: pd.DataFrame,
    *,
    sample_id: str = "<unknown>",
    restrict_to_signatures: Iterable[str] | None = None,
) -> SignatureFit:
    """Fit a patient's 96-channel SBS counts to a mixture of COSMIC signatures.

    Solves  min_{c >= 0}  || sig_matrix @ c - sbs96_counts ||_2.
    Signature matrix columns are non-negative and sum to ~1 per column,
    so the fitted `c` is in count units (not fractions).

    Parameters
    ----------
    sbs96_counts : pandas.Series
        One patient's SBS96 counts, indexed by `SBS96_CHANNELS`.
    cosmic_signatures : pandas.DataFrame
        Rows = SBS96 channels (same ordering), columns = signatures.
    sample_id : str
        Label carried through for provenance.
    restrict_to_signatures : iterable of signature IDs, optional
        If given, the fit uses only these signatures (useful for LUAD
        never-smoker where we restrict to SBS1, SBS5, SBS40a plus a smoking
        prior at zero — but see §9 risk notes before restricting).

    Returns
    -------
    SignatureFit
    """
    if list(sbs96_counts.index) != list(SBS96_CHANNELS):
        sbs96_counts = sbs96_counts.reindex(SBS96_CHANNELS).fillna(0.0)
    if list(cosmic_signatures.index) != list(SBS96_CHANNELS):
        cosmic_signatures = cosmic_signatures.reindex(SBS96_CHANNELS).fillna(0.0)

    if restrict_to_signatures is not None:
        restrict = [s for s in restrict_to_signatures if s in cosmic_signatures.columns]
        sig_matrix = cosmic_signatures[restrict]
    else:
        sig_matrix = cosmic_signatures

    y = sbs96_counts.values.astype(float)
    A = sig_matrix.values.astype(float)

    c, residual = nnls(A, y)
    mix = pd.Series(c, index=sig_matrix.columns, name="contribution")
    reconstructed = pd.Series(A @ c, index=SBS96_CHANNELS, name="reconstructed")

    return SignatureFit(
        sample_id=sample_id,
        mix=mix,
        reconstructed_sbs96=reconstructed,
        residual=float(residual),
        total_mutations=int(sbs96_counts.sum()),
    )


# ---------------------------------------------------------------------------
# Collapse SBS96 vectors to 6-class pyrimidine spectrum
# ---------------------------------------------------------------------------


def collapse_sbs96_to_6class(sbs96: pd.Series) -> pd.Series:
    """Aggregate a 96-channel SBS vector into the 6-class pyrimidine spectrum.

    Parameters
    ----------
    sbs96 : pandas.Series
        Index = SBS96 channel labels. Values are counts OR fractions; the
        resulting 6-class vector is normalized to sum to 1 iff input sums > 0.

    Returns
    -------
    pandas.Series
        Index = ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G") in that order,
        values summing to 1 (or all zero if input sums to zero).
    """
    if list(sbs96.index) != list(SBS96_CHANNELS):
        sbs96 = sbs96.reindex(SBS96_CHANNELS).fillna(0.0)

    out = {cls: 0.0 for cls in SUB_CLASSES}
    for ch, value in sbs96.items():
        out[SBS6_FROM_SBS96[ch]] += float(value)

    total = sum(out.values())
    if total > 0:
        out = {k: v / total for k, v in out.items()}
    return pd.Series(out, index=list(SUB_CLASSES))


# ---------------------------------------------------------------------------
# Cohort-level spectrum with bootstrap over patients
# ---------------------------------------------------------------------------


def bootstrap_cohort_spectrum(
    patient_sbs96: pd.DataFrame,
    *,
    n_bootstrap: int = 10_000,
    weight_by_mutation_count: bool = True,
    seed: int = SEED,
) -> dict[str, object]:
    """Compute cohort-level 6-class spectrum with patient-level bootstrap CIs.

    The cohort spectrum is the patient-weighted mean of per-patient 6-class
    spectra. Weights are either per-patient mutation counts (recommended:
    rare-tumours-dont-dominate) or equal (every patient same weight).

    Parameters
    ----------
    patient_sbs96 : pandas.DataFrame
        Columns = patient / sample IDs. Rows = SBS96 channels (counts).
    n_bootstrap : int
        Number of patient-level bootstrap resamples (default 10,000).
    weight_by_mutation_count : bool
        If True, weight patients by their total mutation count (so a hyper-
        mutator doesn't swamp the signal); if False, equal weights.
    seed : int
        Random seed (defaults to project SEED).

    Returns
    -------
    dict with:
      - "point_estimate" : pandas.Series indexed by the 6 classes.
      - "bootstrap_samples" : (n_bootstrap, 6) ndarray.
      - "ci95" : pandas.DataFrame rows = classes, cols = ["lo", "hi"].
      - "n_patients" : int.
    """
    rng = np.random.default_rng(seed)
    sample_ids = list(patient_sbs96.columns)
    n = len(sample_ids)
    if n == 0:
        raise ValueError("empty patient cohort")

    per_patient_6class = pd.DataFrame(
        {sid: collapse_sbs96_to_6class(patient_sbs96[sid]) for sid in sample_ids}
    ).T  # rows = patients, cols = classes
    weights = (
        patient_sbs96.sum(axis=0).values if weight_by_mutation_count else np.ones(n)
    )
    if weights.sum() == 0:
        raise ValueError("all-zero patient weights; cannot compute cohort spectrum")

    def _weighted_mean(idx):
        sub = per_patient_6class.iloc[idx]
        w = weights[idx]
        if w.sum() == 0:
            return np.zeros(6)
        return (sub.values * w[:, None]).sum(axis=0) / w.sum()

    point = _weighted_mean(np.arange(n))

    boot = np.empty((n_bootstrap, 6))
    for b in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        boot[b] = _weighted_mean(idx)

    lo = np.quantile(boot, 0.025, axis=0)
    hi = np.quantile(boot, 0.975, axis=0)
    ci = pd.DataFrame(
        {"lo": lo, "hi": hi}, index=list(SUB_CLASSES)
    )
    point_series = pd.Series(point, index=list(SUB_CLASSES), name="point_estimate")

    return {
        "point_estimate": point_series,
        "bootstrap_samples": boot,
        "ci95": ci,
        "n_patients": n,
    }
