# Bradeen et al. 2006 — Phase 1 data audit record

**Audit date:** 2026-10-13 21:02 Europe/Rome
**Primary paper text:** `/root/.claude/uploads/fc032df5-1689-5549-9823-8f122cba433e/813a52ec-attachment.txt`
**Audit status:** Methods and Table 1 complete; **Table S3 (nucleotide substitutions) STILL REQUIRED** for proper nucleotide-level supply calculation.

## Citation

Bradeen HA, Eide CA, O'Hare T, Johnson KJ, Willis SG, Lee FY, Druker BJ, Deininger MW.
**Comparison of imatinib mesylate, dasatinib (BMS-354825), and nilotinib (AMN107) in an
N-ethyl-N-nitrosourea (ENU)–based mutagenesis screen: high efficacy of drug combinations.**
*Blood.* DOI: 10.1182/blood-2006-02-004580.

## Mutagenesis protocol

- **Cell line:** Ba/F3-p210BCR-ABL (chronic phase CML construct)
- **Mutagen:** N-ethyl-N-nitrosourea (ENU), Sigma, dissolved at 50 mg/mL in DMSO, stored at -80°C
- **ENU exposure:** 50 μg/mL added to Ba/F3-p210BCR-ABL at 5 × 10⁶ cells/mL for 12–24 hours
- **Post-ENU:** cells washed 3× with RPMI, replated in complete medium, expanded 1 week
- **Resistance screen:** 1 × 10⁵ cells/well in 96-well plates, 200 μL complete media
- **Imatinib concentrations:** 2, 4, 8, 16 μM
- **Observation period:** at least 4 weeks
- **Readout:** DNA from pelleted resistant cells → PCR with specific primers (B2A + ABL4065) →
  Sanger sequencing (commercial, Agencourt) → Mutation Surveyor software analysis
- **Detection sensitivity:** 10–15% mutant allele in a wild-type background

## Headline results

- **763 wells** with growth → sequenced
- **681 wells** harboured BCR-ABL KD mutations
  - 636 wells with single mutations (93.4%)
  - 45 wells with 2 different mutations (6.6%) — counted as independent clones from mixed wells
- **726 total KD mutant clones** across all three drugs
- Imatinib: **N = 372**; Nilotinib: N = 173; Dasatinib: N = 181

## Table 1 imatinib arm — the one we extracted

- Covered in `table1_imatinib_counts.csv`
- 20 distinct amino acid variants
- Dominated by:
  - **Y253H: 168 isolates (45.16%)** — TAC → CAC via pos 1 T>C transition
  - **T315I: 139 isolates (37.37%)** — ACT → ATT via pos 2 C>T transition
  - Together 82% of all imatinib-resistant isolates
- D_eff computed from the paper's counts: **2.87** (matches v0.18 Table 1 reported value of 2.9)

## Comparison to clinical mutations (Results, line 79)

Paper's own claim:
- 25 of 30 clinically-observed residues reproduced in screen (83% residue coverage)
- 16 of 30 identical amino acid substitutions (64%)
- Weighted by clinical frequency: 92% residue coverage, 82% precise-mutation coverage
- 3 novel mutations not described clinically: L248R, F311I, L348M (note L348M is in Table S3,
  not Table 1 — reinforces the need for Table S3)

**This claim is important for our paper's framing**: Bradeen et al. position the ENU screen as
clinically predictive. Our thesis argues this is partly coincidence — ENU's dominant T>C and
C>T transitions overlap with the dominant tumour-signature classes (SBS1 is C>T-rich,
SBS40a/SBS5 include T>C), so the two spectra agree on hotspots like Y253H (T>C) and T315I
(C>T) even without any shared selection logic. The T2 ratio test is designed to disentangle
this.

## Comparability to Azam 2003 for T2 ratio test

| Parameter | Bradeen 2006 ENU | Azam 2003 XL-1R | Comparability |
|---|---|---|---|
| Cell line | Ba/F3-p210BCR-ABL | BaF3 + pEYKBA (BCR-ABL in pEYK3.1) | Both Ba/F3; vectors differ |
| Numbering | type Ia (clinical standard) | type Ib (adds +19) | Converge after Ib→Ia conversion |
| Imatinib range | 2, 4, 8, 16 μM | 5, 10 μM | 5μM Azam ≈ 2–4μM Bradeen pooled |
| Mutagen | ENU 50 μg/mL × 12–24 h on cells | XL-1R mutator strain on plasmid | Different in vivo chemistry |
| Readout | Sanger KD amplicon | Sanger of recovered plasmid | Both Sanger |
| Imatinib total N | 372 | 198 (5 μM alone) / 266 (both) | Both well-powered |
| D_eff | 2.87 (very narrow) | 56.5 at 5 μM (very broad) | Expected under hypothesis |

**Pairing recommendation for T2:** Pair **Azam 5 μM counts** with **Bradeen 2 + 4 μM pooled**
(both low-stringency, broad-permissive-set regimes). This maximises variant coverage and
minimises selection-pressure-mismatch confounding.

## Known gaps for T2 completion

| Gap | Required for | How to resolve |
|---|---|---|
| Bradeen Table S3 (nucleotide substitutions) | Proper nt-level supply calculation for variants with multiple nt routes (e.g., F317C/I/L/V, L248R/V) | Supplementary PDF from Blood |
| XL-1R substitution spectrum | Ratio prediction S_ENU/S_XL1R | Stratagene product literature, Greener 1994, or inference from Azam's own 266-mutation output |
| ABL1 RefSeq CDS (NM_005157.6) | Codon lookups for 90+ residues | WebFetch or Michele |

## Risks and open questions

1. **The 2-mutations-per-well interpretation.** Bradeen counted 45 wells with 2 mutations as
   90 independent clones. If any of these are actually co-mutations in one clone (compound
   resistance variants), the count is inflated by up to 45 for imatinib arm. The paper argues
   against compound interpretation based on allele frequency, but this is a potential source
   of noise in the primary count.
2. **ENU dose difference (Bradeen 50 μg/mL vs Zhang 100 μg/mL).** Affects absolute mutation
   rate, should not affect relative spectrum ratios. Flag in sensitivity analysis.
3. **Vector / construct differences** (Bradeen p210 full-length BCR-ABL vs Azam pEYKBA).
   May affect mutation sampling at non-kinase-domain positions, but kinase-domain output
   should be comparable.
4. **No L348M in Table 1 but mentioned as a novel variant.** Table S3 presumably contains it;
   another reason to source S3.

## Status

- ✅ Confirmed Bradeen 2006 as the ENU + ABL1 + imatinib paper
- ✅ Protocol documented with exact ENU dose and duration
- ✅ Table 1 imatinib arm per-variant counts extracted (`table1_imatinib_counts.csv`)
- ✅ Total 372 imatinib isolates, 20 distinct variants, D_eff = 2.87 verified
- ❌ Table S3 (nucleotide substitutions) NOT YET IN HAND
- ❌ XL-1R substitution spectrum NOT YET SOURCED

T2 can proceed to ratio-test stage once (a) Table S3 arrives OR we accept inferring nucleotide
identities from the single-route constraint for most variants, and (b) XL-1R spectrum is sourced
or substituted with a nominal.
