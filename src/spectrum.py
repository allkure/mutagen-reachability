"""Mutagen substitution spectrum handling.

A spectrum is a dict mapping the six pyrimidine-referenced substitution classes
(C>A, C>G, C>T, T>A, T>C, T>G) to probabilities summing to 1 (modulo rounding).

Spectra are loaded from a single reference YAML file (data/reference/spectra.yaml)
with explicit provenance for every value. Loading enforces validation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .codon_tools import SUB_CLASSES

__all__ = ["Spectrum", "load_spectra", "validate_spectrum"]


@dataclass(frozen=True)
class Spectrum:
    """A six-class pyrimidine-referenced substitution spectrum with provenance.

    Attributes
    ----------
    name : str
        Short label used in figures and tables, e.g. "ENU_Arnold2012".
    values : dict[str, float]
        Keys are the six substitution classes; values sum to ~1.
    source : str
        Full citation or dataset identifier for the empirical measurement.
    method : str
        How the spectrum was measured (NGS of the library, germline survey, etc.).
    uncertainty : str
        One-line note on sampling or measurement uncertainty.
    """

    name: str
    values: dict[str, float]
    source: str
    method: str
    uncertainty: str = ""
    notes: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        validate_spectrum(self.values, name=self.name)

    def at_class(self, cls: str) -> float:
        """Return the probability mass at a given substitution class."""
        return self.values[cls]

    def at_origin(self, origin: str) -> float:
        """Return the total mass at substitutions originating at A:T or C:G base pairs."""
        if origin == "A:T":
            return sum(self.values[c] for c in ("T>A", "T>C", "T>G"))
        if origin == "C:G":
            return sum(self.values[c] for c in ("C>A", "C>G", "C>T"))
        raise ValueError(f"origin must be 'A:T' or 'C:G', got {origin!r}")

    def __getitem__(self, key: str) -> float:  # dict-like convenience
        return self.values[key]


def validate_spectrum(values: dict[str, float], name: str = "<spectrum>") -> None:
    """Check that a spectrum covers all six classes with non-negative values summing to ~1.

    Raises
    ------
    ValueError
        If the spectrum fails validation. Error message names the spectrum.
    """
    missing = set(SUB_CLASSES) - set(values)
    if missing:
        raise ValueError(f"{name}: missing classes {sorted(missing)}")
    extra = set(values) - set(SUB_CLASSES)
    if extra:
        raise ValueError(f"{name}: unknown classes {sorted(extra)}")
    negatives = {c: v for c, v in values.items() if v < 0}
    if negatives:
        raise ValueError(f"{name}: negative values {negatives}")
    total = sum(values.values())
    if not (0.95 <= total <= 1.05):
        raise ValueError(
            f"{name}: values sum to {total:.4f}; expected approximately 1"
        )


def load_spectra(path: str | Path) -> dict[str, Spectrum]:
    """Load every spectrum defined in a YAML reference file.

    The YAML schema is:

        spectra:
          ENU_Arnold2012:
            values: {"C>A": ..., ...}
            source: "..."
            method: "..."
            uncertainty: "..."
            notes: {...}
          MutazymeII_DoiSupp2022:
            ...

    Parameters
    ----------
    path : str | Path
        Location of the YAML file.

    Returns
    -------
    dict[str, Spectrum]
        Mapping from spectrum name to its Spectrum object.
    """
    path = Path(path)
    with path.open() as f:
        raw = yaml.safe_load(f)
    if raw is None or "spectra" not in raw:
        raise ValueError(f"{path}: missing top-level 'spectra' key")
    spectra: dict[str, Spectrum] = {}
    for name, spec in raw["spectra"].items():
        spectra[name] = Spectrum(
            name=name,
            values={k: float(v) for k, v in spec["values"].items()},
            source=str(spec.get("source", "")),
            method=str(spec.get("method", "")),
            uncertainty=str(spec.get("uncertainty", "")),
            notes=dict(spec.get("notes", {})),
        )
    return spectra
