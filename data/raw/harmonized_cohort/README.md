# ALK+ NSCLC patient records — harmonized across cohorts

> **The three data files described here are not distributed with this
> repository.** `patients.csv`, `mutations.csv` and `treatment_lines.csv`
> are built in part from AACR Project GENIE, GENIE BPC, MSK-IMPACT and
> MSK-CHORD, whose data use agreements do not permit redistribution, and
> they carry demographics and survival for 612 individuals. Obtain the
> sources directly and rebuild them with the scripts named below.
>
> Nothing in this repository reads them. The manuscript's pipeline
> reproduces every derived table byte for byte without them, and
> `src/harmonized_cohort.py` is imported elsewhere only for its codon
> utilities. This file is kept as schema documentation.

Built by `scripts/clinical/build_patient_dataset.py` (run it after `build_master_cohort.py`, `extract_gainor_2016.py` and `extract_madsen_2020.py`). It starts
from the deduplicated cohort in `data/cohort/alk_master_patients.csv` and adds:
- the clinical columns shared across cohorts
- mutations in every gene
- every recorded systemic drug
- one treatment-line table pooling all sources

Three tables, linked on `patient_key`:

| file | unit | rows |
|---|---|---|
| `patients.csv` | patient | 612 (572 with `analysis_core = True`) |
| `treatment_lines.csv` | patient × treatment line / regimen | 954 (823 ALK-TKI) in 366 patients |
| `mutations.csv` | patient × gene × protein change, non-synonymous | 3,554 in 505 patients, 1,858 genes |

**Start with `analysis_core == True`.** It excludes 18 ALK structural variants that do not look like
kinase-retaining fusions, and 22 MSK research-cohort patients whose non-`P-` IDs cannot be linked to
MSK-IMPACT, so they may already be in the table.

## Sources

| institution | core patients | contributes |
|---|---|---|
| MSK | 328 | demographics, IMPACT mutations, OS; MSK-CHORD drug timelines with derived rwPFS; 28 are also in GENIE BPC |
| OrigiMed (China) | 112 | demographics, stage, panel mutations; **no outcomes or drugs** |
| MGH (Gainor 2016) | 46 | published per-patient table: TKI sequence, time on drug, ALK mutation at progression; **no demographics, OS or PFS** |
| Aarhus/Aalborg, Denmark (Madsen 2020) | 24 | published per-patient table: demographics, TNM stage, smoking, ALK-TKI sequence with time on each drug; serial **plasma ctDNA** (NGS + ddPCR). **No PFS; OS not pooled** (see below) |
| DFCI, VICC (GENIE BPC) | 29 | chart-curated demographics, stage, regimens, per-regimen PFS, OS |
| NCI, POG, OncoSG, TCGA | 33 | demographics, WES/WGS mutations, OS; drugs for a few |

Deduplication: MSK `P-` IDs are shared across all MSK studies and GENIE. TCGA barcodes are shared between
`luad_tcga_pub` and the PanCan Atlas: 3 patients were counted twice upstream, and here they are merged
into `tcga:<barcode>`. When sources disagree on a field, the most curated source wins, in this order:
GENIE BPC, then MSK-CHORD, then the newest MSK cut, then the rest. The fields where this choice matters
carry a `*_source` column.

## `patients.csv`

| column | notes |
|---|---|
| `patient_key` | `P-xxxxxxx` (MSK), `tcga:…`, `genie_bpc_nsclc:…`, `gainor_2016:MGH…`, `madsen_2020:PT…`, otherwise `<study>:<id>` |
| `institution`, `source_studies` | where the patient comes from; every study that contains them |
| `alk_evidence` | `ngs_fusion` · `alk_specific_tki_only` (on a non-crizotinib ALK TKI, no DNA fusion call) · `clinical_alk_positive` (Gainor, Madsen: IHC/FISH/tissue NGS) · `uncertain_sv` |
| `analysis_core`, `possible_msk_duplicate` | see above |
| `sex` | Female / Male |
| `age_at_dx_years` | age at diagnosis; only from columns documented as at-diagnosis |
| `age_at_seq_years` | age at sequencing / biopsy |
| `race` | White / Asian / Black / Other. Blank includes "unknown/refused". OrigiMed reports none |
| `ethnicity` | Hispanic / Non-Hispanic |
| `smoking_status`, `smoking_pack_years` | Never / Ever. **MSK-CHORD's value is an NLP prediction** (see `smoking_status_source`); GENIE BPC public releases no smoking data |
| `stage_at_dx` | I–IV at diagnosis, from the most curated source available |
| `stage_iv_ever` | MSK-CHORD reports only the *highest stage ever recorded*, which is not stage at diagnosis, so it is stored here instead. Also True when `stage_at_dx` is IV |
| `histology`, `oncotree_code` | most frequent across the patient's samples |
| `n_sequenced_samples`, `sample_types`, `gene_panels` | sequencing metadata |
| `tmb_nonsynonymous_median` | median over samples. **Panel-dependent: do not compare across panels** |
| `alk_fusion_partner`, `alk_fusion_class` | 5' partner; `canonical` / `probable` / `uncertain` |
| `alk_kd_mutations`, `alk_kd_mutation_any` | ALK kinase-domain (aa 1116–1383) missense mutations, from any sample or time point (see `mutations.csv`) |
| `mut_<GENE>` | **1** = mutated (non-synonymous), **0** = the patient's panel covers the gene and no call was made, **blank** = not assessable. Madsen patients only ever get 1 or blank: a negative *plasma* ctDNA result does not show the tumour is wild-type. Genes are those on every MSK-IMPACT version that are mutated in ≥3 % of the cohort: TP53, SETD2, KMT2D, CDKN2A, ARID1A, KMT2C, EGFR, FAT1, SMARCA4, CREBBP |
| `n_nonsyn_mutations` | non-ALK non-synonymous calls (panel-dependent, like TMB); blank for ctDNA-only patients |
| `alk_tki_sequence`, `n_alk_tki_lines`, `first_alk_tki`, `first_alk_tki_time_on_drug_months` | taken from the preferred treatment source (see `treatment_lines.csv`) |
| `systemic_agents_ever`, `received_chemo`, `received_immunotherapy` | every systemic agent in any source. MSK-CHORD timelines can include drugs given for a *prior, different* cancer |
| `os_months`, `os_event`, `os_anchor` | **the anchor differs by source**: `first_sequencing`, `diagnosis_genie_bpc` or `advanced_diagnosis`. Stratify or adjust by anchor; never pool these blindly |
| `pfs_first_tki_months/_event/_source` | imaging-based PFS on the first recorded ALK TKI (GENIE BPC curated, else MSK-CHORD derived) |

Completeness among the 572 core patients:

| field | complete |
|---|---|
| sex | 92 % |
| histology | 91 % |
| TP53 status | 88 % |
| age at diagnosis | 75 % |
| smoking | 75 % |
| OS | 66 % |
| ALK-TKI sequence | 62 % |
| race | 61 % |
| first-TKI PFS | 48 % |
| stage at diagnosis | 47 % |

## `treatment_lines.csv`

One row per line or regimen. Patients in both MSK-CHORD and GENIE BPC have their lines from **both**
sources. Use `preferred_source == True` to keep one copy, which is the BPC curated version for the
24 patients with lines in both.

| column | notes |
|---|---|
| `line_id`, `source`, `provenance` | `curated` (BPC chart review), `derived` (computed from MSK-CHORD timelines), `published_table` (Gainor, Madsen), `raw` |
| `line_number`, `regimen_number` | order within patient × source. BPC regimens include chemo / IO; MSK-CHORD and Madsen rows are **ALK-TKI lines only** (Madsen `regimen_number` = the overall line number given in Table S1) |
| `drugs`, `is_alk_tki`, `alk_tki_agents`, `alk_tki_line`, `prior_alk_tkis` | lower-case generic names, `;`-joined for combinations |
| `time_anchor`, `start_day`, `stop_day` | days from the anchor named in `time_anchor`. The anchors differ by source |
| `time_on_drug_months`, `time_on_drug_event` | duration of the line. For Gainor and Madsen this is the only outcome. Madsen's own "PFS" is defined as start → PD, death or treatment cessation, i.e. this same time-to-treatment-failure (event 0 = still on drug at the 22 Mar 2019 cut) |
| `pfs_months`, `pfs_event`, `pfs_definition` | **primary PFS, imaging-based**: BPC `pfs_i_g`; MSK-CHORD first NLP progression scan |
| `pfs_alt_*` | alternative definition (BPC imaging-or-med-onc; MSK-CHORD confirmed rwPFS) |
| `os_from_line_*` | OS from the start of the line (BPC only) |
| `best_response` | Gainor Table S3 |
| `alk_mut_before_line` / `alk_mut_after_line` | ALK mutations (aa 1058–1400) documented **before** the line started (+30 d grace) vs found only afterwards |

**For mutation → outcome questions, use `alk_mut_before_line`.** Most ALK mutations were found on biopsies
taken *after* the drug failed, so patient-level "has mutation" vs PFS comparisons are biased. Only
**26 ALK-TKI lines in 18 patients** have the mutation documented beforehand: 14 from MSK-CHORD, 7 from BPC,
3 from Madsen (PT3 C1156Y → alectinib, PT5 L1196M → lorlatinib, PT9 L1196M → alectinib; measured in plasma
before the drug started) and 2 from Gainor.

## `mutations.csv`

Non-synonymous calls only: missense, nonsense, frameshift, in-frame, splice site, start/stop loss. They
are pooled over every sample and study and collapsed to patient × gene × `protein_change`. Further columns:
- `max_vaf`
- `n_samples`, `samples`, `studies`
- `alk_kinase_domain`
- `source_type`: `ngs` (tissue/panel), `published_table` (Gainor re-biopsy results), `ctdna_ngs` / `ctdna_ddpcr`
  (Madsen plasma; the ddPCR row is only added where NGS missed the variant: PT9 G1202R). Plasma VAFs are
  not comparable to tissue VAFs.

A non-ALK mutation absent from this table does **not** prove the patient was wild-type: check panel
coverage (the `mut_<GENE>` flags in `patients.csv` already do this).

## Madsen 2020 (added from `cancers-733499-supplementary.xlsx`)

Madsen AT et al., *Cancers* 2020;12(4):947, doi:10.3390/cancers12040947 (PMC7226192, CC BY). 24 ALK+ patients
enrolled prospectively at Aarhus and Aalborg University Hospitals, Denmark, 2015–2018. They cannot overlap
any other patient here, since every other patient-level source is from the US, China, Singapore or Canada.
Parsed by `scripts/clinical/extract_madsen_2020.py`. The full per-study detail, including the S2 ctDNA time
course, S3 ddPCR series, metastatic sites, TNM, study entry point and vital status, is in
`data/clinical/madsen_2020_cancers/processed/`.

Things the source does not settle, and how they were handled:
- **Drug names are not in Table S1** ("ALK TKI 1 … 4 (months)"). They were assigned from the S2 time-point
  labels and Figure 1's colour-coded trajectories. Figure 1's row labels PT5↔PT15 and PT16↔PT24 are swapped
  relative to S1/S2. S1, S2 and the text agree with each other, so their numbering is used.
- **PT17 and PT20 are inconsistent in S1.** PT17 lists one TKI, but S2 has a "prior to alectinib" sample.
  PT20's single 15-month line does not match its 21-month total. Only the S1 line is kept, and its
  `time_on_drug_event` is left blank.
- S1's "total time on ALK TKI" differs from the sum of the per-drug months for PT5 (by 1.0 mo) and PT6
  (by 0.3 mo). The per-drug values are used.
- **OS is not in the pooled `os_months`.** S1 gives months only for the 6 deceased patients ("Alive"
  otherwise, with no follow-up time), and the anchor is not stated. Pooling it would add deaths without
  their censored counterparts. It is kept as `os_months_if_deceased` and `vital_status` in
  `madsen_patients.csv`.
- Stage was derived from the TNM string (AJCC 8th) only where the components determine it. 19 patients are
  stage IV, matching the paper's text.

## Not included

- `data/cohort/shaw2019_pseudo_ipd.csv`: synthetic patients reconstructed from KM curves. It must never
  be mixed with real patients.
- Shaw 2019 waterfall observations: variant-level, not unique patients.
- Copy-number alterations (e.g. CDKN2A deletion) and non-ALK fusions. The raw files are in each study
  folder (`data_cna.txt`, `data_sv.txt`) if needed.
- Non-ALK-TKI drug lines from MSK-CHORD appear only in `systemic_agents_ever`. Turning the per-agent
  timeline into regimens would need a regimen-grouping rule that has not been defined.
