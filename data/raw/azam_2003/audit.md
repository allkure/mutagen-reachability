# Azam et al. 2003 — Phase 1 data audit record

**Audit date:** 2026-10-13
**Primary paper text:** `/root/.claude/uploads/fc032df5-1689-5549-9823-8f122cba433e/98d81dc6-attachment.txt` (provided by MR)
**Audit status:** Partial — paper methods and headline numbers extracted; **Tables 1, 2, 3 with per-variant counts STILL REQUIRED** for the ratio test.

## Citation

Azam M, Latek RR, Daley GQ.
**Mechanisms of Autoinhibition and STI-571/Imatinib Resistance Revealed by Mutagenesis of BCR-ABL.**
*Cell* 2003;112(6):831–843.

## Mutagenesis protocol

- Vector: pEYKBA = BCR-ABL cDNA cloned into pEYK3.1 retroviral vector (Koh et al. 2002)
- Mutagen: **XL-1 Red** (Stratagene) — *E. coli* strain deficient in three major DNA repair pathways
- Procedure: 1 μg BCR-ABL plasmid transformed into XL-1 Red, selected on Zeomycin for 18 hr,
  plasmid DNA isolated and transfected into HEK 293T
- Retrovirus titer: 10⁵–10⁶ IU/mL
- Target cells: BaF3 (IL-3-dependent hematopoietic)
- Selection: 5–10 μM STI-571/imatinib in soft agar, 10 days
- Readout: single colonies plucked, expanded in liquid STI-571; genomic DNA extracted; integrated
  provirus recovered as plasmid OR specific regions amplified by PCR; Sanger sequencing of the
  recovered inserts

## Headline results (verified from paper text)

- **283 drug-resistant colonies** recovered (NOT 175 as v0.18 draft had; v0.18 is wrong)
- **112 distinct amino acid substitutions**
- **90 residues** affected
- **68% of positions were identified more than once** (so ~61 residues with ≥2 isolates, ~29 residues with exactly 1 isolate)
- 65 variants reproduced by site-directed mutagenesis; 59/65 (91%) verified drug-resistant;
  6/65 (9%) false-positive rate

## Numbering convention

Azam uses **type Ib** human c-Abl numbering. Modern BCR-ABL clinical literature uses **type Ia**.
Ib numbering is +19 relative to Ia:

| Clinical name (Ia) | Azam name (Ib) | Residue |
|---|---|---|
| Y253 | Y272 | p-loop |
| E255 | E274 | p-loop (region) |
| T315 | T334 | gatekeeper |
| F317 | F336 | drug contact |
| M351 | M370 | SH2-α3 contact |
| H396 | H415 | activation loop |
| F382 | F401 | activation loop (kinase-dead when mutated) |

Our analysis will use type Ia numbering throughout (modern clinical standard) with Ib given in a
column of the final supplementary table for cross-reference.

## Variants named in paper text (per-variant counts PENDING table extraction)

### Kinase domain (Table 1 of the paper)

- **P-loop region:** L267R/V (Ia L248R/V; novel), Y272 variants (Ia Y253; patient), E274K/V (Ia E255K/V; patient implied)
- **Drug contact / ATP cleft:** T334 variants incl. T334I (Ia T315I; gatekeeper), F336 (Ia F317), A288V, M309L/T, G340W, I379F
- **β3 sheet (novel region in VAVK):** multiple novel variants
- **α1 catalytic helix (novel):** E298K, E300K
- **α3 helix (SH2 contact):** V357G, V358A/G, M362T (Ia M343T), A363V, Q365H, M370T/I (Ia M351T/I; patient), E371K
- **Activation loop:** V398E/A (Ia V379; near-patient site), L403M (Ia L384M; novel), M407I (Ia M388I; novel),
  G417R (Ia G398R; novel), H415R/P (Ia H396R/P; patient)
- **C lobe distal:** E513, E518K/I, I521, E528
- **SH3 contacts in kinase:** K313 (Ia K294), Q319 (Ia Q300; Q319H characterised), L320, E335

### Mutations outside kinase domain (Table 2)

- **Cap region (novel):** K70Q, E71K, N72D/T, A75V, G76E, P77R, E79K
- **SH3 FVALYD motif:** V86G, A87E/V, Y89D/H/N, D90N
- **SH3 WVPS motif:** W129C, P131 variants
- **SH3-SH2 connector:** V138G, S140R
- **SH2 domain:** S167C, S176T, V224E (3 false-positives); S206F/P (true positive; EF loop)
- **SH2-CD linker:** K238E, T243A, P249 variants

### Mutations in STI-571 contact sites (Table 3)

- V275A, M337A (kinase-active but drug-resistance indistinguishable from WT; i.e., the paper used them as a negative control for the ratio of mutation-to-resistance)
- I332, F401, K290, E305 (kinase-dead when mutated; not in the resistance list)

## Known missing from this text extract

**The actual Tables 1, 2, and 3 with per-variant isolate counts are not in the text I received.**
The text extract shows "Open table in a new tab" placeholders and summarises Table contents in prose
but does not reproduce the +, ++, +++ isolate-count notation per variant.

For the T2 ratio test, I need the per-variant counts from Tables 1 & 2 (Table 3 is a method control).
Candidates for sourcing: the Cell 2003 article HTML (behind reCAPTCHA for my fetcher); the Cell
supplementary materials; or the author's lab website.

## Risks and open questions

1. **Count discrepancy with v0.18.** The v0.18 draft wrote "n=175" for Azam 2003. The paper
   reports 283 drug-resistant colonies. v0.18 is wrong by an unknown source of error. Our Table 1
   must use 283.
2. **Variants recovered more than once but at unknown multiplicity.** 68% of 90 residues = ~61
   residues had ≥2 isolates. Of the 112 distinct variants, the top-recovery variants carry most
   of the ratio-test signal — we need the specific isolate counts for at least Y272 (Y253),
   T334 (T315), E274K (E255K), M370T (M351T), H415R (H396R), and L267R/V to compare against
   Bradeen 2006 ENU equivalents.
3. **Co-isolate mutations.** Azam notes "some recovered viruses carried two or more mutations"
   and that false-positive "piggy-back" mutations exist. Isolate counts in the tables may double-count
   residues when a single plasmid carries two mutations. This should be checked.
4. **XL-1 Red substitution spectrum.** The spectrum of this mutator strain is not published in a
   single clean reference I know of. Candidates: Stratagene product literature; Greener & Callahan
   1994; a Daley-lab characterisation if one exists. If no spectrum, T2 is reported with Azam counts
   but no quantitative ratio-prediction — only the Bradeen-vs-Azam qualitative pattern can be shown.

## Status

- ✅ Confirmed Azam 2003 is the XL-1 Red + ABL + imatinib screen
- ✅ Protocol and conditions documented
- ✅ Headline numbers (283 colonies, 112 distinct substitutions, 90 residues) recorded
- ❌ Per-variant counts from Tables 1 & 2 NOT YET IN HAND
- ❌ XL-1 Red substitution spectrum NOT YET SOURCED

Decision gate: whether T2 proceeds depends on obtaining Table 1 and Table 2, and on finding or
accepting a nominal XL-1 Red spectrum. If Tables unobtainable, T2 drops to qualitative secondary.
