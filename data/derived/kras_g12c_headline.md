# KRAS G12C case study — headline results (codon-locked, 2026-10-14)

## Framework
- **σ atlas**: Feng et al. 2022 PNAS, DOI 10.1073/pnas.2120512119.
  Twist Biosciences synthetic saturation mutagenesis library of KRAS-G12C.
  Chemistry-free: every single-AA substitution enters the pool at uniform
  baseline (vs oligo synthesis bias only). 3,567 variants × 2 drugs
  (sotorasib + adagrasib) scored in NCI-H358.
- **Decomposition**: log P̂_simulated(v) = log S_mutagen(v) + Z_feng(v).
  Feng's Z is already normalised to control (DMSO Day 0); no separate
  fitness term needed.
- **Codon table**: KRAS-4A NM_033360.4, residues 1-150 ClinVar-anchored at
  V8, G12/C12, S17, Q61, R68, H95, Y96; residues 151-189 (4A HVR) consensus
  (immaterial to Awad clinical panel). Position 12 flipped to TGT in the
  library (G12C baked in).

## Awad 2021 clinical panel — codon-locked categorisation

For each of 12 Awad 2021 clinical resistance variants, SNV reachability from
the KRAS G12C WT codon:

| category | variants | fraction |
|----------|----------|----------|
| **SNV-UNREACHABLE** (require 2+ nt changes from WT codon) | R68L, V8M | 2 / 12 |
| MULTI-CLASS (SNV-reachable in 1 nt via multiple pyrimidine classes) | V8L, R68S, H95Q | 3 / 12 |
| CLASS-FORCED (SNV-reachable in exactly one class) | C12F, Y96C, Y96D, Y96H, H95D, H95R, R68M | 7 / 12 |

## Headline — the two kinds of chemistry undercount

### 1. SNV-UNREACHABLE (structural, chemistry-INDEPENDENT undercount)

| variant | WT codon | Feng Z (sot) | Feng Z (ada) | SNV reachable? |
|---------|----------|--------------|--------------|----------------|
| **R68L** | AGG | **+3.47 (#1 sot)** | **+4.34 (#1 ada)** | NO (AGG→CTR/TTR needs 2 nt) |
| V8M | GTA | +1.02 | +1.10 | NO (GTA→ATG needs 2 nt) |

**R68L is the #1-ranked variant in Feng's atlas under both drugs**, yet it
cannot arise from the WT codon AGG by any single-nucleotide change. No ENU,
Mutazyme, EMS, or any other SNV mutagen screen of KRAS G12C can ever produce
it. The Boehringer AACR 2022 poster is consistent with this: their ENU-based
Ba/F3 screen named only Y96D and Y96S as hits; R68L is absent despite being
the top clinically-relevant variant by Feng's chemistry-free measurement.
R68L required the field's shift to Twist-synthesised and prime-editing
libraries to be discoverable in a preclinical context.

### 2. CLASS-FORCED (chemistry-dependent undercount)

| variant | drug | class | rank_chemfree | rank_MutII | rank_ENU | shift_MutII | shift_ENU |
|---------|------|-------|---------------|------------|----------|-------------|-----------|
| **H95D** | adagrasib | C>G | **12** | 90 | **1060** | +78 | **+1048** |
| Y96D | sotorasib | T>G | 9 | 223 | 107 | +214 | +98 |
| Y96D | adagrasib | T>G | 23 | 375 | 286 | +352 | +263 |
| R68M | sotorasib | C>A | 81 | 154 | 433 | +73 | +352 |
| R68M | adagrasib | C>A | 81 | 163 | 427 | +82 | +346 |

**H95D under adagrasib + ENU** is the canonical case: C>G is the rarest class
in ENU (0.2% baseline), so a clinically prominent resistance variant drops
from rank 12 to rank 1060 under ENU simulation. No ENU-based KRAS G12C
preclinical screen has reported H95D despite its clinical presence in
adagrasib-treated patients (Awad 2021 NEJM).

**Y96D** was described by Tanaka 2021 Cancer Discovery as "a novel,
previously unidentified tyrosine-to-aspartate mutation" — precisely because
prior ENU/Mutazyme screens never found it (rank 107-223 in simulation vs
rank 9 chemistry-free). Twist-based screens (Boehringer AACR 2022) did
subsequently catch Y96D, consistent with chemistry-free = find-it.

## Interpretation for the paper

KRAS G12C is a stronger case than ALK because the SNV-UNREACHABLE category
exists: it demonstrates that chemistry undercount is not only quantitative
(rare class → low ranking) but also **structural** (no SNV can reach the
variant at all). R68L is Feng's top-1 under BOTH drugs — the single most
important variant for predicting sotorasib/adagrasib resistance by the
chemistry-free atlas — and it is structurally invisible to any SNV mutagen.

This is why the field needed Twist and prime-editing to catch what ENU
screens were systematically missing, and why the clinical panel precedes
the preclinical panel for KRAS G12C.

## Caveats stated honestly
- 7 / 12 clinical variants are class-forced; 3 are multi-class (chemistry
  effect is averaged across classes); 2 are SNV-unreachable (chemistry-
  independent). The headline predictions use only the 7 class-forced ones
  (plus the 2 unreachable ones as a separate structural claim).
- Feng Z is a single-replicate-normalised readout, noisier than Oh 2026's
  two-term atlas. Spearman ρ(simulated, chemistry-free) across all 981
  class-forced variants: 0.52 (MutII), 0.42 (ENU). Lower than ALK (0.73-0.88)
  but still highly significant.
- Codon table positions 151-189 (4A HVR) are consensus not ClinVar-locked;
  this is immaterial to Awad clinical variants which are all in 1-96.
- The AGG codon used for R68 is confirmed by ClinVar p.Arg68Ser c.204G>C
  (which anchors position 3 = G). If a given patient's R68L arose via a
  non-G12C-library route (e.g., compound serial selection), that is a
  different mechanism from a mutagen-screen single-step mutation and is
  outside this framework's scope.
