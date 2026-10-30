"""Statistics — bootstrap, permutation, correlation, with pre-registered seed.

All random operations in this project consume numpy.random.Generator seeded with
20261013. Any statistic in the paper must be reproducible from this seed.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

__all__ = [
    "SEED",
    "make_rng",
    "log_ratio",
    "bootstrap_correlation",
    "permutation_label_shuffle",
    "CorrelationResult",
]

SEED: int = 20261013


def make_rng(seed: int | None = None) -> np.random.Generator:
    """Return a numpy Generator seeded with the project-wide seed by default."""
    return np.random.default_rng(seed if seed is not None else SEED)


def log_ratio(a: np.ndarray, b: np.ndarray, pseudo: float = 0.5) -> np.ndarray:
    """Log-ratio with pseudocount for zero-count robustness.

    Zero cells in either `a` or `b` are handled by adding a pseudocount of 0.5
    to BOTH numerator and denominator (this is the Haldane-Anscombe correction).

    Parameters
    ----------
    a, b : np.ndarray
        Count vectors of equal length.
    pseudo : float, default 0.5
        Pseudocount added to each cell.

    Returns
    -------
    np.ndarray
        log((a + pseudo) / (b + pseudo))
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return np.log((a + pseudo) / (b + pseudo))


@dataclass(frozen=True)
class CorrelationResult:
    """Point estimate, bootstrap CI, and permutation P-value for a correlation."""

    r: float
    ci95: tuple[float, float]
    permutation_p: float
    n_variants: int
    n_bootstrap: int
    n_permutations: int

    def __str__(self) -> str:
        lo, hi = self.ci95
        return (
            f"r = {self.r:.3f}  95% CI ({lo:.3f}, {hi:.3f})  "
            f"permutation P = {self.permutation_p:.4g}  "
            f"n = {self.n_variants} variants"
        )


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson correlation; returns 0 for degenerate (zero-variance) inputs."""
    x = x - x.mean()
    y = y - y.mean()
    denom = np.sqrt((x * x).sum()) * np.sqrt((y * y).sum())
    if denom == 0:
        return 0.0
    return float((x * y).sum() / denom)


def bootstrap_correlation(
    observed_log_ratio: np.ndarray,
    predicted_log_ratio: np.ndarray,
    n_bootstrap: int = 10_000,
    rng: np.random.Generator | None = None,
) -> tuple[float, tuple[float, float]]:
    """Pearson r with bootstrap 95% CI by variant resampling.

    Variants are resampled with replacement; the Pearson correlation of the
    (observed, predicted) log-ratios is recomputed on each bootstrap sample.

    Parameters
    ----------
    observed_log_ratio : np.ndarray
        Observed log(P_A / P_B) per variant.
    predicted_log_ratio : np.ndarray
        Predicted log(S_A / S_B) per variant.
    n_bootstrap : int
        Number of bootstrap resamples.
    rng : np.random.Generator, optional
        Random generator; defaults to make_rng() with the project seed.

    Returns
    -------
    (r_point, (ci_lo, ci_hi))
        Point estimate and the 2.5/97.5 percentile bounds.
    """
    if rng is None:
        rng = make_rng()
    x = np.asarray(observed_log_ratio, dtype=float)
    y = np.asarray(predicted_log_ratio, dtype=float)
    n = x.size
    if n != y.size:
        raise ValueError(f"length mismatch: {n} observed vs {y.size} predicted")
    r_point = _pearson(x, y)
    rs = np.empty(n_bootstrap)
    for b in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        rs[b] = _pearson(x[idx], y[idx])
    lo, hi = np.percentile(rs, [2.5, 97.5])
    return r_point, (float(lo), float(hi))


def permutation_label_shuffle(
    counts_A: np.ndarray,
    counts_B: np.ndarray,
    predicted_log_ratio: np.ndarray,
    n_permutations: int = 10_000,
    rng: np.random.Generator | None = None,
    pseudo: float = 0.5,
) -> tuple[float, float]:
    """Two-sided permutation P-value for the chemistry hypothesis.

    Null: screen composition is independent of mutagen. Under the null, the
    total count at each variant is fixed, but the SHARE in screen A vs B is
    exchangeable across variants.

    Implementation: for each variant v with total t_v = counts_A[v] + counts_B[v],
    sample a_v* from Binomial(t_v, s*), where s* is the overall fraction in A
    (pooled across all variants). Recompute log-ratio and correlation.

    Parameters
    ----------
    counts_A, counts_B : np.ndarray
        Observed clone counts per variant.
    predicted_log_ratio : np.ndarray
        log(S_A / S_B) per variant (held fixed under the null).
    n_permutations : int
        Number of null draws.
    rng : np.random.Generator, optional
    pseudo : float
        Pseudocount for log_ratio.

    Returns
    -------
    (r_observed, p_value)
    """
    if rng is None:
        rng = make_rng()
    a = np.asarray(counts_A, dtype=int)
    b = np.asarray(counts_B, dtype=int)
    pred = np.asarray(predicted_log_ratio, dtype=float)
    if a.shape != b.shape or a.shape != pred.shape:
        raise ValueError("shape mismatch between counts and predicted log-ratio")
    total_A = a.sum()
    total_B = b.sum()
    pooled_share_A = total_A / (total_A + total_B)
    totals = a + b
    r_obs = _pearson(log_ratio(a, b, pseudo=pseudo), pred)
    r_null = np.empty(n_permutations)
    for p in range(n_permutations):
        a_null = rng.binomial(totals, pooled_share_A)
        b_null = totals - a_null
        r_null[p] = _pearson(log_ratio(a_null, b_null, pseudo=pseudo), pred)
    # two-sided P: fraction of null r at least as extreme in |r|
    p_value = float(np.mean(np.abs(r_null) >= abs(r_obs)))
    # Add-1 smoothing to avoid P = 0 on finite permutations
    p_value = (int(np.sum(np.abs(r_null) >= abs(r_obs))) + 1) / (n_permutations + 1)
    return float(r_obs), p_value
