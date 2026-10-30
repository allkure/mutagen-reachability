# Phase 2 clinical-data provenance — GENIE + MSK-IMPACT

**Status:** data access pending user action. Framework code assumes files
land here with the schemas named below.

## Required inputs

Phase 2 needs three cuts of the AACR Project GENIE pooled release
(which includes MSK-IMPACT as a contributing center) plus two helper files.

### 1. The GENIE release version and files

- **Target release:** AACR GENIE **v16.1-public** (January 2025 cut) OR the
  closest later version MR has access to. The release version is logged in
  `data/raw/genie_msk/release_metadata.yaml` (to be created at ingest).
- **Required files** (from the public Synapse release, syn7222066):
  - `data_mutations_extended.txt` — all somatic mutations, MAF-like, with
    `Hugo_Symbol`, `Chromosome`, `Start_Position`, `End_Position`,
    `Reference_Allele`, `Tumor_Seq_Allele2`, `Variant_Classification`,
    `HGVSp_Short`, `HGVSc`, `Tumor_Sample_Barcode`,
    `t_ref_count`, `t_alt_count`.
  - `data_clinical_patient.txt` — patient-level, with `PATIENT_ID`,
    `SEX`, `PRIMARY_RACE`, `ETHNICITY`, `BIRTH_YEAR`.
  - `data_clinical_sample.txt` — sample-level, with `SAMPLE_ID`,
    `PATIENT_ID`, `CANCER_TYPE`, `CANCER_TYPE_DETAILED`, `ONCOTREE_CODE`,
    `AGE_AT_SEQ_REPORT`, `SEQ_ASSAY_ID`, `SEQ_YEAR`.
  - `data_fusions.txt` (or `data_sv.txt` for later releases) — structural
    variants including fusions, with `Hugo_Symbol`, `Fusion`,
    `Tumor_Sample_Barcode`.
  - `genomic_information.txt` — assay-panel information, with
    `SEQ_ASSAY_ID` → gene list, needed for denominator corrections.

### 2. MSK-IMPACT clinical annotation enrichment

GENIE's public release does NOT include smoking history or TKI treatment
history. We need these for both the LUAD never-smoker signature cohort
(T1/T3 mutagen source) and the crizotinib- / alectinib-refractory cohorts
(T1/T3 clinic side). Two possible sources:

- **MSK cBioPortal** `msk_impact_2017` or newer study — includes
  `SMOKING_STATUS` and some treatment annotation (file
  `data_clinical_patient.txt` or `data_clinical_sample.txt`).
- **MSK supplementary tables** published alongside LUAD-specific papers
  (e.g., Jordan 2017 *Cancer Discovery*, Rizvi 2018, Shaw 2016). For
  ALK-fusion-positive LUAD specifically, the Gainor 2016 and Lin 2019
  cohort tables may be the cleanest source.

MR to confirm which snapshot. The code expects a merged `clinical.tsv`
with columns `PATIENT_ID`, `SMOKING_STATUS` (never / former / current /
unknown), `ALK_FUSION_PARTNER`, `PRIOR_ALK_TKI` (comma-separated),
`PROGRESSED_ON_ALK_TKI` (comma-separated), `BIOPSY_TIMING` (pretreatment
/ post-TKI).

### 3. COSMIC v3.3 reference signatures

- **Source:** https://cancer.sanger.ac.uk/signatures/sbs/ (COSMIC v3.3 GRCh38)
- **File:** `COSMIC_v3.3.1_SBS_GRCh38.txt` (~550 KB, tab-separated, 96 rows
  × 72 signatures). Save as `data/raw/cosmic_v33/COSMIC_v3.3.1_SBS_GRCh38.txt`.
- **License:** CC0 for the signature matrices (verify in download header).

### 4. GRCh38 reference FASTA (for trinucleotide-context annotation)

- Only chr9 (ABL1) and chr2 (ALK) kinase-domain regions needed, not whole
  genome.
- Can be built from pyfaidx on Ensembl Human GRCh38 primary-assembly
  (chromosome 2: 29192773-29921586 covers ALK; chromosome 9: 130713946-130887675
  covers ABL1). Alternatively, UCSC twoBit regions.

## Sample-ID subsets we will construct

Three cohorts, written to `data/derived/cohorts/`:

### LUAD never-smoker signature cohort (feeds T1 and T3 ratio tests)

- `CANCER_TYPE_DETAILED == "Lung Adenocarcinoma"`
- `SMOKING_STATUS in {"Never", "never", "Lifetime Non-smoker"}`
- `ALK_FUSION in {any} == True` (i.e., sample has an ALK fusion event)
- Any treatment history (we want the mutagenic signature of the
  LUAD-never-smoker tumour, which should be dominated by SBS1/SBS5/SBS40a,
  independent of TKI exposure)
- Expected N: 100-300 patients (based on Lin 2019 cohort size estimates).

### Crizotinib-refractory cohort (T1 clinic side)

- `CANCER_TYPE_DETAILED == "Lung Adenocarcinoma"`
- `ALK_FUSION == True`
- `PRIOR_ALK_TKI` list includes "crizotinib" and NO second-/third-gen ALK
  TKI (alectinib, brigatinib, ceritinib, lorlatinib, ensartinib)
- `PROGRESSED_ON_ALK_TKI` includes "crizotinib"
- Expected N: 30-80 patients (post-crizotinib ALK kinase-domain mutations
  are the classic Gainor 2016 cohort).

### Alectinib-refractory cohort (T3 clinic side)

- `CANCER_TYPE_DETAILED == "Lung Adenocarcinoma"`
- `ALK_FUSION == True`
- `PROGRESSED_ON_ALK_TKI` includes "alectinib"
- Prior crizotinib ALLOWED (post-alectinib-refractory is almost always
  prior-crizotinib-refractory as well)
- Expected N: 50-150 patients.

For each cohort, the variant counts are tallied over the ALK kinase domain
residues F1116-G1392 (NM_004304.5) and filtered through the pre-registered
§5 inclusion/exclusion rules (class-forced single-nt-reachable).

## Scope note for the signature fit

The LUAD never-smoker signature fit uses the WHOLE-TUMOUR somatic mutation
catalogue (restricted to protein-coding regions covered by the GENIE panel
intersection, to avoid panel-size confounders). This gives us a 6-class
pyrimidine spectrum per patient, which is then averaged (with per-patient
weights) across the cohort to yield the LUAD-never-smoker spectrum for use
in T1 and T3 supply-weight calculations. See `src/signatures.py` and
`notebook/02_luad_signature.py`.
