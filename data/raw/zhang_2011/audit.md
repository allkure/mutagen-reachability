# Zhang et al. 2011 — Phase 1 data audit record

**Audit date:** 2026-10-13
**Primary paper text:** `/root/.claude/uploads/fc032df5-1689-5549-9823-8f122cba433e/5347b627-attachment.txt` (provided by MR)
**Audit status:** Partial — Table 2 extracted; Figure 2A per-variant counts PENDING.

## Citation

Zhang S, Wang F, Keats J, Zhu X, Ning Y, Wardwell SD, Moran L, Mohemmad QK, Anjum R,
Wang Y, Narasimhan NI, Dalgarno D, Shakespeare WC, Miret JJ, Clackson T, Rivera VM.
**Crizotinib-Resistant Mutants of EML4-ALK Identified Through an Accelerated
Mutagenesis Screen.** *Chemical Biology & Drug Design.* First published 31 October 2011.
DOI: `10.1111/j.1747-0285.2011.01239.x`.

Author affiliation: ARIAD Pharmaceuticals.

## Mutagenesis protocol

- Cell line: Ba/F3 expressing native EML4-ALK variant 1 (Genbank BAF73611.1, pMSCV-Neo)
- Mutagen: **N-ethyl-N-nitrosourea (ENU)**, Sigma Aldrich
- ENU exposure: 100 μg/mL, **overnight**
- Post-ENU selection: 96-well plates with crizotinib at 250, 500, 720, 1000, 1440, 2000 nM
- Growth medium: standard, **without IL-3** (parental Ba/F3 require IL-3; EML4-ALK-driven cells do not)
- Selection duration: 5 weeks
- Sequencing: Taq DyeDeoxy Terminator Cycle Sequencing (Sanger) of ALK kinase region from
  wells with substantial outgrowth

## Headline results

- **Total mutations recovered: 422**
- **Distinct residues mutated: 16**
- Site count by crizotinib concentration:
  - 500 nM: 15 sites
  - 720 nM: 8 sites
  - 1000 nM: 6 sites
  - 1440 nM: 2 sites
  - 2000 nM: 0 (no outgrowth)
- Residues at the highest crizotinib concentrations (highly resistant):
  C1156, I1171, F1174, L1196, S1206, G1269

## Characterized variants (Table 2 of paper)

Ten variants selected as "most frequently identified, each at a different residue":

| Variant | Crizotinib IC50 (nM) ENU clone | p-ALK IC50 (nM) ENU clone |
|---|---|---|
| T1151K | 231 ± 106 | ND |
| L1152V | 237 ± 25 | ND |
| C1156Y | 489 ± 48 | 408 ± 67 |
| I1171T | 393 ± 134 | ND |
| F1174C | 479 ± 40 | 165 ± 112 |
| **L1196M** | **981 ± 113** | **1162 ± 209** |
| S1206R | 681 ± 176 | 356 ± 173 |
| E1210K | 318 ± 143 | ND |
| F1245C | 425 ± 100 | ND |
| G1269S | 953 ± 213 | 1366 ± 52 |

Native EML4-ALK IC50: 132 ± 45 nM (crizotinib viability).

## F1174 cluster (text, line 115)

> "the mutations of F1174 to **Cys, Val, Ile, and Leu** were observed in our screen"

Therefore F1174C, F1174V, F1174I, F1174L all recovered. F1174L is **mixed-origin from TTC**
(three routes: pos-1 T>C A:T, pos-3 C>A C:G, pos-3 C>G C:G) and by the pre-registration
inclusion rule must be excluded from the strict A:T vs C:G primary test.

## Codon L1196 (gatekeeper)

**The paper names ONLY L1196M at codon 1196.** It is characterized in Table 2 as the single
most resistant of the three top mutants (L1196M, S1206R, G1269S), with crizotinib IC50
981 nM in ENU clones and 1215 nM in reintroduced clones.

**L1196Q and L1196P are not mentioned anywhere in the paper.**

The v0.18 draft of the previous manuscript claimed "13 L1196Q + 5 L1196M" from "Yoda 2018."
This paper is not Yoda 2018, and it does not report L1196Q at any count. The source of the
"13 Q + 5 M" datum is **unknown as of this audit** and must be either located in a different
primary paper or dropped from any derived claim.

## Phase 1 gaps for this paper

Pending data for a per-variant ratio test:

1. **Figure 2A per-variant counts** across the 15 sites detected at 500 nM crizotinib. Table 2
   gives IC50s for 10 of them but does not give counts. Fig 2A caption describes "the
   relative percentage of the indicated EML4-ALK kinase domain mutant among recovered
   clones" at each crizotinib concentration — so Fig 2A has what we need.

2. **Supplementary data** (Figure S1, S2 referenced; also potential supplementary tables with
   full variant lists).

3. **Residue identities for the 6 additional sites** not in Table 2 (the paper names V1180,
   R1181, D1268, Y1278, E1241 in the Discussion as structural contacts; we need to confirm
   whether all of these are among the 16 or only some).

## Known differences from Takei 2022 / Doi 2022

If T1 is reframed to use Zhang 2011 (ENU, crizotinib) vs any Doi 2022 arm:

- Doi 2022 screens **alectinib, repotrectinib, ensartinib** — NOT crizotinib. **σ cancellation
  across the pair requires drug matching, which fails here.** T1 as pre-registered needs a
  Mutazyme II crizotinib screen that we have not located.
