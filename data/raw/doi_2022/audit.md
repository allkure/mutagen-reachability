# Doi 2022 Target Oncol — Phase 1 data audit

**Audit date:** 2026-10-13 22:25 Europe/Rome
**Primary paper PDF:** `/root/.claude/uploads/.../93a30537-art3A10.10072Fs11523-022-00919-5.pdf`
**Supplementary PDF:** `/root/.claude/uploads/.../6e629416-11523_2022_919_MOESM1_ESM.pdf`

## Citation

Doi Y, Tagaya H, Noge A, Semba K.
**Prediction of resistance mutations against upcoming anaplastic lymphoma kinase inhibitors.**
*Targeted Oncology.* 2022;17:695-707. DOI: 10.1007/s11523-022-00919-5.

## Mutagenesis protocol

- **Mutagen:** GeneMorph II Random Mutagenesis Kit (Agilent/Stratagene) with Mutazyme II DNA polymerase.
- **Target amplicon:** EML4-ALK v1 cDNA, ALK kinase-domain region, GenBank AB663645 pos 1850-2383 (534 bp). Covers major clinical ALK resistance residues.
- **GC content of amplicon:** 58% (paper explicitly notes this biases the observed Mutazyme II spectrum).
- **Mutation frequency delivered:** 0.48 mutations/kb (NGS-measured on the pre-selection library).
- **Pre-selection library size:** ~7 × 10⁵ colonies (recovered from E. coli transformation of PCR products ligated into pMXs vector arms).
- **Cell line:** Ba/F3 expressing EML4-ALK v1, infected with retroviral library at >10% infection efficiency (3 × 10⁷ Ba/F3 cells + library virus).
- **Alectinib concentration:** 50 nM for 1 week.
- **Readout:** PCR amplicon from genomic DNA of surviving cells → NGS of the amplicon → per-position mutant-allele frequency.

## Mutazyme II substitution spectrum

Doi 2022 Supplementary Fig 2b shows the Mutazyme II substitution spectrum
measured by NGS on the pre-selection library (n=3 replicates):

| Pyrimidine class | % of substitutions |
|---|---:|
| C>T (= G>A, C>T) | 47 |
| T>C (= A>G, T>C) | 18 |
| C>A (= G>T, C>A) | 18 |
| T>A (= A>T, T>A) | 10 |
| C>G (= G>C, C>G) | 5 |
| T>G (= A>C, T>G) | 4 |

This is stored as `MutazymeII_DoiSupp2022` in `data/reference/spectra.yaml`
and cross-validated against Rasila 2009 Table 1 Mutazyme II arm — see
`MutazymeII_Rasila2009` entry.

## Alectinib-resistance screen (T3 preclinical side)

Fig 3a bar chart (per-nucleotide-position mutant-allele frequency) identifies
the following resistance mutations:

- **I1171N/T/S** (compound at codon 1171) — ~0.20 total at position 3514
- **V1180L** — ~0.21 at position 3534 (single highest peak)
- **L1196M** — ~0.03 at position 3584
- **L1198F** — ~0.02 at position 3594
- **G1202R** — ~0.02 at position 3604

Silent passengers also visible: L1152 silent, L1165 silent, R1192 silent
(these survived alectinib selection but are NOT resistance mutations; they
are dropped from the T3 analysis).

### Known gaps

- Exact numeric per-variant counts are NOT published for Fig 3a — only the
  figure. Our `fig3a_alectinib_counts.csv` records visual bar heights at
  ±2 percentage points precision.
- I1171N/T/S compound label is NGS-position-resolution ambiguity; we
  equal-split across the 3 reachable sub-variants (same rule as Azam 2003
  compound rows).
- The "others" category in Fig 3a is not visible as a labeled fraction.

## T3 clinic pairing

T3 clinic side is the alectinib-refractory subset of MR's harmonized ALK+
NSCLC cohort (data/raw/harmonized_cohort/), with ALK KD mutations
documented as `alk_mut_after_line` on an alectinib treatment line.
23 patients, 28 isolates (preliminary count 2026-10-13 23:40 Europe/Rome).

## Status

- Alectinib variant list extracted from Fig 3a
- Visual counts recorded in `fig3a_alectinib_counts.csv`
- Mutazyme II spectrum confirmed vs Rasila 2009 (two independent sources agree)
- T3 preclinical input ready; T3 ratio-test can run now
