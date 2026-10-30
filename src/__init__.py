"""Project-level package: cross-mutagen supply framework for TKI resistance screens.

Modules:
    codon_tools    — SNV route enumeration on RefSeq codons.
    spectrum       — 6-class substitution spectra with provenance.
    supply_weights — per-variant S_d(v) computation.
    stats          — bootstrap correlation, permutation null, pre-registered seed.

All random operations in this project use seed = 20261013 (stats.SEED).
"""

from . import codon_tools, spectrum, stats, supply_weights  # noqa: F401
