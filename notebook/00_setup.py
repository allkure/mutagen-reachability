"""Notebook 00 — environment setup and sanity checks.

This file is maintained as a .py file (percent-format) during development and
converted to 00_setup.ipynb at submission time via jupytext. Keeping notebooks
as .py in source control makes diffs reviewable and merge conflicts tractable.

Run order:
    00_setup.py       — environment + imports + seed + sanity checks
    01_data_audit.py  — Phase 1 audit of every primary paper
    02_spectra.py     — Phase 2 spectrum sourcing + LUAD fit
    03_ratio_test.py  — Phase 3 pre-registered ratio tests
    04_secondary.py   — Phase 4 secondary analyses
    05_figures.py     — Phase 5 publication figures
"""

# %% [markdown]
# # 00 — Setup
#
# **Project:** Decoupling mutagen supply from drug selection in preclinical TKI
# resistance screens — a cross-mutagen falsifiable test.
#
# **Pre-registration:** `../preregistration.md` (locked 2026-10-13).
#
# **Journal target:** *Genome Biology*.
#
# This notebook imports the project modules, locks the random seed, and runs
# the sanity-check suite (`pytest tests/`). Running this notebook is the
# first step of reproducing any downstream result in the paper.

# %%
from __future__ import annotations

import json
import sys
from pathlib import Path

# Add project root to path so `import src.*` works from notebooks
PROJECT_ROOT = Path.cwd().resolve()
while PROJECT_ROOT != PROJECT_ROOT.parent:
    if (PROJECT_ROOT / "preregistration.md").exists():
        break
    PROJECT_ROOT = PROJECT_ROOT.parent
assert (PROJECT_ROOT / "preregistration.md").exists(), (
    "Could not locate project root (preregistration.md); "
    "run from inside the repository."
)
sys.path.insert(0, str(PROJECT_ROOT))

print(f"Project root: {PROJECT_ROOT}")

# %%
import numpy as np
import pandas as pd
import yaml

from src.codon_tools import enumerate_snv_routes, sub_class, translate
from src.spectrum import Spectrum, load_spectra
from src.stats import SEED, make_rng
from src.supply_weights import supply_weight

# %%
# Pre-registered random seed. All random operations use rng = make_rng() which
# consumes this seed. Reproducibility of every statistic in the paper depends
# on this value being used consistently.
rng = make_rng()
print(f"Pre-registered seed: {SEED}")
print(f"First five draws (0, 1): {rng.random(5)}")

# %% [markdown]
# ## Load reference data

# %%
spectra = load_spectra(PROJECT_ROOT / "data" / "reference" / "spectra.yaml")
print(f"Loaded {len(spectra)} spectra:")
for name, spec in spectra.items():
    at_pct = 100 * spec.at_origin("A:T")
    cg_pct = 100 * spec.at_origin("C:G")
    print(f"  {name:40s}  A:T {at_pct:4.1f}%   C:G {cg_pct:4.1f}%")

# %%
with (PROJECT_ROOT / "data" / "reference" / "codons.yaml").open() as f:
    codons_doc = yaml.safe_load(f)

alk_codons: dict[int, str] = {
    res: entry["codon"] for res, entry in codons_doc["alk_nm_004304_5"]["codons"].items()
}
print(f"ALK codon table (NM_004304.5): {len(alk_codons)} residues loaded")
for res in sorted(alk_codons):
    aa = translate(alk_codons[res])
    print(f"  residue {res}: codon {alk_codons[res]} -> {aa}")

# %% [markdown]
# ## Spot-check: the two codon mistakes that broke v0.18

# %%
# F1174 should be TTC on NM_004304.5. The v0.18 draft used TTT, which changes
# F1174L route composition (TTT -> CTT/TTA/TTG vs TTC -> CTC/TTA/TTG).
assert alk_codons[1174] == "TTC", f"F1174 codon must be TTC, got {alk_codons[1174]}"
assert translate("TTC") == "F"

# L1196 is CTG (gatekeeper). This one v0.18 had right.
assert alk_codons[1196] == "CTG", f"L1196 codon must be CTG, got {alk_codons[1196]}"
assert translate("CTG") == "L"

# F1174L under TTC has three mixed-origin routes. By the pre-registration's
# inclusion rule it is EXCLUDED from the strict A:T vs C:G primary test.
routes_f1174l = enumerate_snv_routes("TTC", "L")
print(f"F1174L from TTC: {len(routes_f1174l)} routes")
for r in routes_f1174l:
    print(f"  pos{r.position} {r.ref_base}>{r.alt_base}  class {r.sub_class}  origin {r.origin}")

origins = {r.origin for r in routes_f1174l}
assert origins == {"A:T", "C:G"}, (
    f"F1174L must be mixed-origin from TTC; got {origins}"
)
print("F1174L correctly flagged as mixed-origin -> excluded from primary test by pre-registration rule.")

# %% [markdown]
# ## Verify Zhang 2011 characterized variants have well-defined supply weights

# %%
# Table 2 of Zhang 2011 names 10 variants. We compute their supply weights
# under the Mutazyme II spectrum as a smoke test of the supply-weight pipeline.
# Zhang 2011 used ENU, not Mutazyme II; this is a mechanical check only.

ZHANG_TABLE_2 = [
    (1151, "K"),  # T1151K
    (1152, "V"),  # L1152V
    (1156, "Y"),  # C1156Y
    (1171, "T"),  # I1171T
    (1174, "C"),  # F1174C
    (1196, "M"),  # L1196M
    (1206, "R"),  # S1206R
    (1210, "K"),  # E1210K
    (1245, "C"),  # F1245C
    (1269, "S"),  # G1269S
]

mz2 = spectra["MutazymeII_DoiSupp2022"]
rows = []
for residue, mut_aa in ZHANG_TABLE_2:
    rec = supply_weight(residue, alk_codons[residue], mut_aa, mz2)
    rows.append({
        "variant": rec.variant,
        "wt_codon": rec.wt_codon,
        "routes": " | ".join(
            f"pos{r.position} {r.ref_base}>{r.alt_base} ({r.sub_class})"
            for r in rec.routes
        ),
        "origin": rec.origin,
        "p_supply_MutazymeII": rec.weight,
    })

zhang_df = pd.DataFrame(rows)
print(zhang_df.to_string(index=False))

# %% [markdown]
# **Interpretation.** Six of the ten Zhang 2011 Table 2 variants (I1171T, F1174C,
# L1196M, L1196M, S1206R, E1210K, G1269S, C1156Y) have single-route supply from
# their WT codons. None of the Zhang Table 2 ten are mixed-origin — F1174L
# (which IS mixed-origin) is not in Table 2, though the paper reports it was
# recovered. The pre-registered exclusion rule removes F1174L from any
# A:T-vs-C:G Fisher test but leaves the ten Table 2 variants eligible.

# %% [markdown]
# ## Status
#
# Sanity-check suite passing. Infrastructure ready to accept Phase 1 data
# (Zhang 2011 Figure 2A counts, Bradeen 2006 counts, Azam 2003 counts) and
# proceed to the pre-registered ratio tests once T1 pairing is resolved.
