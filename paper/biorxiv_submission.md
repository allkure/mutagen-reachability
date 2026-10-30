# bioRxiv submission package — metadata and checklist

**Prepared:** 2026-10-16
**Author:** Michele Ronco
**Status:** Ready to submit (verify DOIs below before upload).

---

## 1. Title

Chemical reachability limits mutagen-based TKI-resistance screens

## 2. Running title (≤ 60 characters)

Reachability limits in preclinical TKI-resistance screens

## 3. Author list

| # | Name | Affiliation | ORCID | Email | Role |
|---|---|---|---|---|---|
| 1 | Michele Ronco | European Commission Joint Research Centre, Via E. Fermi 2749, I-21027 Ispra (VA), Italy | [0000-0002-2160-2452](https://orcid.org/0000-0002-2160-2452) | michele.ronco@ec.europa.eu | Corresponding author; sole author |

**Statement:** The views expressed are the author's and do not represent the position of the European Commission.

## 4. Abstract (mirrors `paper.tex` — regenerate on each change via `make abstract-sync`, or paste from the paper's Background/Results/Conclusions blocks)

> **Background.** Preclinical Ba/F3 mutagenesis screens are a first line of evidence for predicting kinase-inhibitor resistance, and their hit lists shape resistance panels for KRAS G12C, ALK, and MET inhibitors. Every mutagen in use has a non-uniform substitution spectrum, and whether that chemistry limits which resistance variants a screen can find, or only reorders them, has not been audited against a chemistry-free baseline.
>
> **Methods.** We separate a screen's expected output into a supply term, from codon tables and each mutagen's six-class spectrum, and a selection term, from chemistry-free deep mutational scans of KRAS G12C, ALK, and MET: P(v) ∝ S_mutagen(v)·exp(βσ(v)). We test supply by enumerating what each chemistry cannot reach, and selection by asking whether supply improves the ranking of the hits four retrospective screens reported.
>
> **Results.** Chemistry acts as a reachability threshold, not a ranking weight. Base editors, which produce only two of the six substitution classes, cannot supply eight of the twelve KRAS G12C panel variants or most of each ALK drug's top-ranked variants. A second limit binds every single-nucleotide mutagen: roughly three-quarters of each MET inhibitor's top 100 priorities need two or more nucleotide changes. Past that threshold supply does not improve the ranking on any anchor, and mutagens with very different spectra recover ALK variants at indistinguishable rates. Two MET screens sharing one mutagen and one cell line gave different class compositions (P = 0.029), so composition follows the drug, though clone counts within a single codon did track supply.
>
> **Conclusions.** The deliverable is the reachability inventory: which resistance variants a mutagen-and-drug pair cannot produce, enumerable from codon tables and a published spectrum before the screen runs. Reweighting the output afterward does not recover them, so chemistry belongs in a screen's design, not in its interpretation.

## 5. Keywords (6)

1. Preclinical resistance screens
2. Mutational signatures
3. Deep mutational scanning
4. KRAS G12C inhibitors
5. ALK inhibitors
6. Ba/F3 mutagenesis

## 6. bioRxiv subject categories

- **Primary:** Bioinformatics
- **Secondary:** Cancer Biology

(bioRxiv allows one primary + one secondary; alternatives if either is unavailable: Pharmacology and Toxicology, Systems Biology.)

## 7. License

**CC BY 4.0** (Attribution). Rationale: strongest attribution-preserving license that allows reuse and derivative work; required by *Genome Medicine* / BMC open-access workflows and consistent with the preprint community's expectations.

## 8. Funding statement

No external funding. The research was conducted outside the author's European Commission Joint Research Centre programme assignments.

## 9. Competing-interests statement

The author declares no competing interests.

## 10. Data and code availability

- **Primary data (reused):**
  - Feng et al. 2022 KRAS G12C Twist Biosciences saturation library Z-scores: PNAS Supplementary Dataset S1 (DOI 10.1073/pnas.2120512119).
  - Oh et al. 2026 ALK kinase-domain prime-editing atlas: *Genome Biology* 27:72 (Supplementary Table S2 and Table S11 WT CDS). DOI 10.1186/s13059-026-03977-4.
  - Doi et al. 2022 ALK Mutazyme II with alectinib Ba/F3 screen: *Targeted Oncology* 17:695–707 (DOI 10.1007/s11523-022-00919-5), Fig. 3a.
  - Awad et al. 2021 KRAS G12C clinical cohort and accompanying *in vitro* DMS: *New England Journal of Medicine* 384:2382–2393 (DOI 10.1056/NEJMoa2105281).
  - Tanaka et al. 2021 KRAS G12C case report: *Cancer Discovery* 11:1913–1922 (DOI 10.1158/2159-8290.CD-21-0365).
- **Derived data and analysis code:** will be deposited on Zenodo with a reserved DOI upon acceptance and tagged release of the private GitHub repository (available to reviewers on request). Includes codon tables, mutagen spectra, derived rank tables, and the figure-generating scripts.
- **Software versions used:** Python 3.11, pandas 2.2, numpy 1.26, matplotlib 3.8, scipy 1.11. Exact versions are pinned in `environment.yml` within the code repository.

## 11. Figures (4 main + 2 in-manuscript tables + Supplementary Table S1)

| # | Caption (short) | File |
|---|---|---|
| Fig. 1 | KRAS G12C illustrative decomposition and resistance-panel reachability. | `figure_1.pdf` / `figure_1.png` |
| Fig. 2 | Mutagen pyrimidine-class spectra and class-coverage score. | `figure_2.pdf` / `figure_2.png` |
| Fig. 3 | Chemistry blind-spots matrix and top-K recovery with 95 % Wilson CI bands. | `figure_3.pdf` / `figure_3.png` |
| Fig. 4 | Three-panel retrospective validation (Doi 2022 / Takei 2025 / Nuvalent 2025). | `figure_4.pdf` / `figure_4.png` |
| Table 1 | Preclinical TKI-resistance screens published 2022–2026 and their mutagen chemistries. | *in-manuscript* |
| Table 2 | Benchmark: σ-alone, framework (log S + σ), and supply-alone on four retrospective anchors, with class-label permutation p and three proxies for the industry anchor. | *in-manuscript* |
| Suppl. Table S1 | Per-variant provenance and codon reachability for both resistance panels. | `suppl_table_s1.pdf` |

## 12. Suggested reviewers

See cover letter for the four preferred referees (Fraser UCSF, Bloom Fred Hutch, Nik-Zainal Cambridge, Lin MGH) and the conflicts-of-interest list.

## 13. Suggested handling editor

The translational-methods character of the paper — genomics-grounded methodology directed at preclinical screens that feed into clinical resistance panels and treatment-sequencing decisions — fits the editorial profile of *Genome Medicine*'s precision-oncology desk.

## 14. Journal target after preprint

- **Preprint:** bioRxiv (this submission).
- **Journal:** *Genome Medicine*, Research Article.
- **Rationale (short):** The paper's primary payload is a predictive, auditable inventory of clinically-plausible resistance variants that specific preclinical mutagen chemistries cannot produce — a direct bridge from genomic data (chemistry-free DMS atlases plus codon-level reachability) to the clinical resistance-panel and treatment-sequencing decisions *Genome Medicine* publishes. The forward-looking neladalkib and next-generation KRAS G12C predictions, and the signature-concordance tagging that distinguishes smoking-NSCLC from never-smoker LUAD contexts, are specifically the translational framing *Genome Medicine* profiles.

---

## bioRxiv submission checklist

- [x] Title and running title drafted
- [x] Abstract (≤ 300 words)
- [x] Author / affiliation / email
- [x] **ORCID** — 0000-0002-2160-2452
- [x] Keywords
- [x] Subject categories
- [x] License (CC BY 4.0)
- [x] Funding statement
- [x] Competing-interests statement
- [x] Data & code availability (reviewers-only GitHub during review → Zenodo on acceptance)
- [x] Main manuscript PDF (`paper.pdf`)
- [x] Figure PDFs (`figure_1.pdf`, `figure_2.pdf`, `figure_3.pdf`, `figure_4.pdf`)
- [x] Supplementary Table S1 PDF (`suppl_table_s1.pdf`)
- [x] Markdown version of manuscript (`paper.md`) — for personal archive; bioRxiv upload is PDF
- [x] Cover letter PDF (`cover_letter.pdf`) — not required by bioRxiv but prepared for journal submission
- [x] **Oh et al. 2026 DOI verified** — 10.1186/s13059-026-03977-4 (Genome Biology 27:72, published Feb 2026)
- [x] **Doi et al. 2022 DOI verified** — 10.1007/s11523-022-00919-5 (Targeted Oncology 17(6):695-707)
- [ ] Zenodo-DOI reservation (do this after bioRxiv upload so the Zenodo record cites the preprint DOI)

---

## Order of operations

1. Fill in ORCID in `paper.tex` and `biorxiv_submission.md`.
2. ~~Verify Oh 2026 final DOI~~ — DONE (verified against the published record).
3. Recompile `paper.tex`, `suppl_table_s1.tex`, `cover_letter.tex` → PDFs.
4. Upload to bioRxiv: paper.pdf as main manuscript; figure PDFs as separate figure files; suppl_table_s1.pdf as supplementary file; this `biorxiv_submission.md` is for your records.
5. After bioRxiv DOI is issued, update the manuscript's "Preprint DOI" line in the Methods and resubmit to bioRxiv as a revision (standard practice).
6. Submit to *Genome Medicine* with `cover_letter.pdf`, `paper.pdf`, figure PDFs, `suppl_table_s1.pdf`, and the bioRxiv preprint DOI listed in both the Methods and the cover letter.
