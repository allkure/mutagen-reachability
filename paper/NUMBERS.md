# Paper numbers — master book

*Last computed: 2026-10-30 19:02:21*
*Source script: `scripts/compute_paper_numbers.py`*
*Total registered: 208 entries*

## Golden rule

Every quantitative claim in `paper.tex` comes from this manifest. If a number appears in the paper but is not registered below, that is a bug: either the number is wrong, or this script needs a new entry. Numbers from citations (years, DOIs, journal volume/issue) are not registered here; they live in `refs.bib`.

---

## Index by key

- `alk_alec_enu_top100` = **63% [53\%--72\%]** — ENU-simulated top-100 recovery for ALK, alectinib
- `alk_alec_mut_top100` = **65% [55\%--74\%]** — Mutazyme-simulated top-100 recovery for ALK, alectinib
- `alk_c1156y_cf_rank` = **132** — C1156Y chemistry-free rank among the 1,427 class-forced ALK variants, alectinib
- `alk_c1156y_sigma` = **+1.1** — Oh atlas alectinib Resistance score for C1156Y
- `alk_c1156y_sim_rank` = **55** — C1156Y Mutazyme II-simulated rank among the 1,427 class-forced ALK variants, alectinib
- `alk_l1152r_cf_rank` = **225** — L1152R chemistry-free rank among the 1,427 class-forced ALK variants, alectinib
- `alk_l1152r_sigma` = **+0.7** — Oh atlas alectinib Resistance score for L1152R
- `alk_l1152r_sim_rank` = **822** — L1152R Mutazyme II-simulated rank among the 1,427 class-forced ALK variants, alectinib
- `alk_l1196m_cf_rank` = **1** — L1196M chemistry-free rank among the 1,427 class-forced ALK variants, alectinib
- `alk_l1196m_sigma` = **+7.7** — Oh atlas alectinib Resistance score for L1196M
- `alk_l1196m_sim_rank` = **1** — L1196M Mutazyme II-simulated rank among the 1,427 class-forced ALK variants, alectinib
- `alk_l1198f_cf_rank` = **18** — L1198F chemistry-free rank among the 1,427 class-forced ALK variants, alectinib
- `alk_l1198f_sigma` = **+4.6** — Oh atlas alectinib Resistance score for L1198F
- `alk_l1198f_sim_rank` = **3** — L1198F Mutazyme II-simulated rank among the 1,427 class-forced ALK variants, alectinib
- `alk_largest_panel_shift` = **+597** — largest predicted rank shift among the ALK panel variants under Mutazyme II with alectinib
- `alk_lor_enu_top100` = **65% [55\%--74\%]** — ENU-simulated top-100 recovery for ALK, lorlatinib
- `alk_lor_mut_top100` = **72% [63\%--80\%]** — Mutazyme-simulated top-100 recovery for ALK, lorlatinib
- `alk_panel_class_forced` = **10** — codon reachability of ALK panel
- `alk_panel_multi_class` = **two** — codon reachability of ALK panel
- `alk_panel_total` = **12** — union of alectinib Mutazyme II screen hits + next-gen ALK TKI screen hits
- `alk_snv_top100_recovery_max` = **72%** — max top-100 recovery across ALK drugs for ENU and Mutazyme II
- `alk_snv_top100_recovery_min` = **63%** — min top-100 recovery across ALK drugs for ENU and Mutazyme II
- `alk_v1180l_cf_rank` = **8** — V1180L chemistry-free rank among the 1,427 class-forced ALK variants, alectinib
- `alk_v1180l_sigma` = **+5.3** — Oh atlas alectinib Resistance score for V1180L
- `alk_v1180l_sim_rank` = **15** — V1180L Mutazyme II-simulated rank among the 1,427 class-forced ALK variants, alectinib
- `alk_zot_enu_top100` = **63% [53\%--72\%]** — ENU-simulated top-100 recovery for ALK, zotizalkib
- `alk_zot_mut_top100` = **65% [55\%--74\%]** — Mutazyme-simulated top-100 recovery for ALK, zotizalkib
- `atlas_cg_share` = **17.0%** — share of class-forced ALK atlas variants whose forced class is C>G
- `atlas_dose_ratio_alectinib` = **1.7x** — atlas assay concentration for alectinib divided by the measured Ba/F3 wild-type IC50
- `atlas_dose_ratio_lorlatinib` = **3.7x** — atlas assay concentration for lorlatinib divided by the measured Ba/F3 wild-type IC50
- `atlas_dose_ratio_zotizalkib` = **30.6x** — atlas assay concentration for zotizalkib divided by the measured Ba/F3 wild-type IC50
- `atlas_vs_own_repeat_rho` = **0.55--0.62** — agreement between the main single-dose atlas and its best-matching sub-library condition: two runs of nominally the same assay
- `be_alk_peak_at_k_le_100` = **38%** — max CBE+ABE recovery over K ≤ 100 across all ALK drugs
- `be_alk_top10_max` = **10%** — max CBE+ABE top-10 recovery across ALK drugs
- `be_alk_top10_min` = **0.00%** — min CBE+ABE top-10 recovery across ALK drugs (alectinib shows 0)
- `be_alk_top10_outside_max` = **100%** — 1 - min(BE top-10 recovery across ALK drugs)
- `be_alk_top10_outside_min` = **90%** — 1 - max(BE top-10 recovery across ALK drugs)
- `be_classes_supplied` = **two of six** — CBE supplies C>T, ABE supplies T>C by design; the other 4 classes are 0
- `be_kras_panel_class_blocked` = **7 of the 12** — KRAS G12C panel variants reachable in one nucleotide step whose every route falls in a class CBE+ABE cannot make
- `be_kras_panel_supplied` = **3** — KRAS G12C panel variants a CBE+ABE pair can produce
- `be_kras_peak_at_k_le_100` = **100%** — max CBE+ABE recovery over K ≤ 100 across all KRAS G12C drugs
- `cg_detect_prob_germline` = **0.86** — detection probability of a C>G variant at the germline ENU rate in a 1e7 Ba/F3 library; shows a rare class is still supplied
- `cg_expected_carriers_germline` = **2.0** — expected carriers of a given C>G variant at the germline ENU rate
- `doi2022_hits_class_forced` = **five** — Doi 2022 hits that are class-forced (G1202R is multi-class)
- `doi2022_hits_count` = **six** — single-amino-acid resistance mutations in Doi 2022 Fig 3a, excluding silent passengers (L1152, L1165, R1192) and the I1171S sub-variant of the I1171N/T/S compound call
- `doi_anchor_topN_bound` = **top-120** — smallest round bound containing all five class-forced Doi hits in the Mutazyme-simulated ranking
- `doi_i1171n_rank_sim` = **118** — I1171N rank under the Mutazyme-simulated alectinib ranking
- `doi_i1171t_rank_sim` = **10** — I1171T rank under the Mutazyme-simulated alectinib ranking
- `doi_l1196m_rank_sim` = **1** — L1196M rank under the Mutazyme-simulated alectinib ranking
- `doi_l1198f_rank_sim` = **3** — L1198F rank under the Mutazyme-simulated alectinib ranking
- `doi_v1180l_rank_sim` = **15** — V1180L rank under the Mutazyme-simulated alectinib ranking
- `dose_cross_drug_matched_rho` = **0.75** — median rank correlation between DIFFERENT drugs assayed at matched dose rank; exceeds the same-drug across-dose figure, indicating shared batch structure
- `dose_cross_drug_unmatched_rho` = **0.57** — same, for cross-drug pairs at unmatched dose rank
- `dose_rho_adjacent_median` = **0.64** — median rank correlation between adjacent concentration steps
- `dose_rho_low_high_median` = **0.39** — median across the seven dose series of the lowest-to-highest rank correlation
- `dose_rho_low_high_range` = **0.30 to 0.63** — Spearman between a drug's lowest and highest assayed concentration, across the seven multi-dose series in Oh Tables S3, S5 and S9
- `dose_topk_overlap_range` = **between 5 and 15 of the top 25** — variants shared between the top 25 at a drug's lowest and highest dose
- `dose_window_sd_ratio` = **1.1 to 1.5** — sigma dispersion in the exon-23 sub-library window relative to the full atlas; above 1 rules out attenuation by range restriction
- `enu_class_C_A` = **9.0%** — Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- `enu_class_C_G` = **0.20%** — Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- `enu_class_C_T` = **17.9%** — Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- `enu_class_T_A` = **31.4%** — Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- `enu_class_T_C` = **36.6%** — Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- `enu_class_T_G` = **5.0%** — Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- `enu_classes_at_1pct` = **five of six** — count of ENU classes supplied at ≥ 1 %; C>G is at 0.2 % and falls below
- `enu_well_supplied_well_vs_total` = **nine of twelve** — count of Nuvalent class-forced hits whose forced class ENU supplies at ≥ 1 %
- `feng_atlas_total_variants` = **3,567** — rows in Feng 2022 KRAS G12C Twist atlas Table S1
- `feng_class_forced_count` = **981** — class-forced KRAS G12C variants in the Feng atlas (unique amino-acid changes)
- `h95d_ada_chemfree_rank` = **9** — H95D chemistry-free rank in the KRAS G12C class-forced atlas under adagrasib
- `h95d_ada_empirical_lower_shift_beta_max` = **+518** — maximum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU empirical lower (C>G = 7%)
- `h95d_ada_empirical_lower_shift_beta_min` = **+2** — minimum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU empirical lower (C>G = 7%)
- `h95d_ada_empirical_lower_shift_beta_range` = **+2 to +518** — H95D shift interval across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU empirical lower (C>G = 7%)
- `h95d_ada_empirical_mean_shift_beta_max` = **-3** — maximum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU empirical mean (C>G = 25%)
- `h95d_ada_empirical_mean_shift_beta_min` = **-4** — minimum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU empirical mean (C>G = 25%)
- `h95d_ada_empirical_mean_shift_beta_range` = **-4 to -3** — H95D shift interval across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU empirical mean (C>G = 25%)
- `h95d_ada_enu_empirical_lower_rank_sim` = **26** — H95D simulated rank under adagrasib + ENU with C>G = 7.0% (Jeffreys 95% lower bound from 3-of-13 anomaly)
- `h95d_ada_enu_empirical_lower_shift` = **+17** — H95D rank shift under adagrasib + ENU (Jeffreys 95% lower bound from 3-of-13 anomaly)
- `h95d_ada_enu_empirical_mean_rank_sim` = **5** — H95D simulated rank under adagrasib + ENU with C>G = 25.0% (Jeffreys posterior mean from 3-of-13 anomaly)
- `h95d_ada_enu_empirical_mean_shift` = **-4** — H95D rank shift under adagrasib + ENU (Jeffreys posterior mean from 3-of-13 anomaly)
- `h95d_ada_enu_empirical_upper_rank_sim` = **2** — H95D simulated rank under adagrasib + ENU with C>G = 49.7% (Jeffreys 95% upper bound from 3-of-13 anomaly)
- `h95d_ada_enu_empirical_upper_shift` = **-7** — H95D rank shift under adagrasib + ENU (Jeffreys 95% upper bound from 3-of-13 anomaly)
- `h95d_ada_enu_germline_rank_sim` = **839** — H95D simulated rank under adagrasib + ENU with C>G = 0.2% (Arnold 2012 germline-ENU compilation)
- `h95d_ada_enu_germline_shift` = **+830** — H95D rank shift under adagrasib + ENU (Arnold 2012 germline-ENU compilation)
- `h95d_ada_germline_shift_beta_max` = **+848** — maximum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU germline (C>G = 0.2%)
- `h95d_ada_germline_shift_beta_min` = **+15** — minimum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU germline (C>G = 0.2%)
- `h95d_ada_germline_shift_beta_range` = **+15 to +848** — H95D shift interval across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU germline (C>G = 0.2%)
- `h95d_ada_mut_rank_sim` = **56** — H95D simulated rank under adagrasib + Mutazyme II (NGS-measured spectrum)
- `h95d_ada_mut_shift` = **+47** — H95D rank shift under adagrasib + Mutazyme II (NGS-measured)
- `h95d_ada_mutazyme_shift_beta_max` = **+606** — maximum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for Mutazyme II (NGS; C>G = 5%)
- `h95d_ada_mutazyme_shift_beta_min` = **+1** — minimum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for Mutazyme II (NGS; C>G = 5%)
- `h95d_ada_mutazyme_shift_beta_range` = **+1 to +606** — H95D shift interval across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for Mutazyme II (NGS; C>G = 5%)
- `implied_enu_cg_supply_interval` = **2.6%--26.2%** — Jeffreys interval on the hit proportion mapped back through the model-implied null to an implied C>G supply rate; reported only to show it is narrower than the uncorrected inversion, and NOT adopted, because the inversion presumes supply controls hit composition
- `kras_panel_case_report` = **one** — patient-observed Y96D from Tanaka 2021 case report
- `kras_panel_class_forced` = **seven** — codon reachability: variants reachable via exactly one pyrimidine class
- `kras_panel_clinical_cohort` = **five** — patient-observed variants from Awad 2021 clinical cohort
- `kras_panel_clinical_total` = **six** — patient-observed KRAS G12C variants (clinical cohort + case report)
- `kras_panel_dms` = **six** — variants from Awad 2021 accompanying in vitro deep mutational scan
- `kras_panel_multi_class` = **three** — codon reachability: variants reachable via multiple pyrimidine classes
- `kras_panel_snv_unreachable` = **two** — codon reachability analysis in the G12C-mutant background
- `kras_panel_total` = **12** — union of KRAS G12C clinical cohort + case report + accompanying in vitro DMS
- `kras_r68l_z_ada` = **+4.34** — chemistry-free Z-score of R68L under adagrasib
- `kras_r68l_z_sot` = **+3.47** — chemistry-free Z-score of R68L under sotorasib
- `kras_snv_top100_recovery_max` = **54%** — max top-100 recovery across KRAS G12C drugs for ENU and Mutazyme II
- `kras_snv_top100_recovery_min` = **40%** — min top-100 recovery across KRAS G12C drugs for ENU and Mutazyme II
- `lorlatinib_inter_study_spread` = **0.32** — median inter-study range of measured lorlatinib potency on the same variant
- `met_amg_clones` = **60** — resistant clones reported for the AMG 458 screen
- `met_amg_median_rank` = **229** — median chemistry-free atlas rank of the AMG 458 hits
- `met_amg_ta_share` = **33\%** — share of AMG 458 clones whose substitution is T>A
- `met_amg_twostep_top100` = **75** — variants in the chemistry-free top 100 for AMG 458 that need two or more nucleotide changes and so lie outside one-step reach
- `met_atlas_size` = **5,458** — MET kinase-domain missense variants scored in the Estevam atlas on the primary background
- `met_nvp_bvu972_clones` = **85** — resistant clones reported for the NVP-BVU972 screen
- `met_nvp_bvu972_median_rank` = **171** — median chemistry-free atlas rank of the NVP-BVU972 hits
- `met_nvp_bvu972_twostep_top100` = **72** — variants in the chemistry-free top 100 for NVP-BVU972 that need two or more nucleotide changes and so lie outside one-step reach
- `met_two_screen_chi2` = **12.4** — chi-square homogeneity of the six-class clone composition between the NVP-BVU972 and AMG 458 screens, which used the same ENU protocol in the same cells
- `met_two_screen_p` = **0.029** — p-value of that homogeneity test; the two screens differ, so the class composition of a hit list is set by the drug, not the mutagen
- `met_wt_median_rank_amg` = **2,472** — median atlas rank of the AMG 458 hits on the wild-type MET background
- `met_wt_median_rank_nvp` = **495** — median atlas rank of the NVP-BVU972 hits on the wild-type MET background
- `met_y1230_clones_y1230c` = **31** — clones carrying Y1230C in the Tiedt NVP-BVU972 screen
- `met_y1230_clones_y1230d` = **1** — clones carrying Y1230D in the Tiedt NVP-BVU972 screen
- `met_y1230_clones_y1230f` = **0** — clones carrying Y1230F in the Tiedt NVP-BVU972 screen
- `met_y1230_clones_y1230h` = **20** — clones carrying Y1230H in the Tiedt NVP-BVU972 screen
- `met_y1230_clones_y1230n` = **19** — clones carrying Y1230N in the Tiedt NVP-BVU972 screen
- `met_y1230_clones_y1230s` = **4** — clones carrying Y1230S in the Tiedt NVP-BVU972 screen
- `met_y1230_missing` = **Y1230F** — one-step Y1230 variants the screen did not recover, with their rank in the chemistry-free atlas
- `met_y1230_recovered` = **5 of 6** — single-nucleotide missense routes out of the Y1230 wild-type codon that the NVP-BVU972 screen recovered
- `met_y1230_supply_fold` = **7.3-fold** — ratio of the largest to smallest germline ENU supply among the substitution classes those Y1230 routes span
- `met_y1230_tg_expected_clones` = **5.1** — clones the two T>G routes out of Y1230 are predicted to carry under germline ENU supply alone
- `met_y1230_tg_observed_clones` = **5** — clones the two T>G routes out of Y1230 actually carried in the Tiedt NVP-BVU972 screen
- `mut_class_C_A` = **18.0%** — NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- `mut_class_C_G` = **5.0%** — NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- `mut_class_C_T` = **47.0%** — NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- `mut_class_T_A` = **10.5%** — NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- `mut_class_T_C` = **18.0%** — NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- `mut_class_T_G` = **3.6%** — NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- `mut_class_rate_ratio` = **13-fold** — ratio of Mutazyme II's most abundant class to its rarest
- `mut_classes_at_1pct` = **6 of 6** — count of Mutazyme II classes supplied at ≥ 1 %; minimum class is T>G at 3.6 %
- `nuvalent_cg_binom_p_composition` = **0.39** — binomial tail against atlas class composition alone
- `nuvalent_cg_binom_p_model` = **4.7e-05** — binomial tail of 3-of-13 C>G against the model-implied null; the germline spectrum is rejected
- `nuvalent_cg_binom_p_model_beta_range` = **4.1e-10 to 4.7e-05** — range of the model-implied binomial tail across beta in [0.25, 4]
- `nuvalent_cg_binom_p_raw` = **2.3e-06** — exact binomial tail against the raw germline spectrum; retained only as the confounded comparison
- `nuvalent_cg_binom_p_selection` = **0.64** — binomial tail of 3-of-13 C>G against a selection-only null; the count is unremarkable once supply is removed
- `nuvalent_cg_count` = **three** — Nuvalent hits class-forced into C>G (ENU 0.2 % class per germline spectrum)
- `nuvalent_cg_pi_model` = **0.56%** — model-implied probability that a hit is C>G, with germline supply and pooled sigma at beta=1
- `nuvalent_cg_top100_sigma` = **24 of the** — C>G count among the top-100 sigma-ranked class-forced ALK variants, pooled across drug arms
- `nuvalent_class_forced` = **twelve** — Nuvalent hits that are class-forced single-nt reachable (R1192P outside atlas, G1202R multi-class)
- `nuvalent_total` = **thirteen** — single-amino-acid variants reported in Nuvalent AACR 2025 poster aggregate across arms
- `nuvalent_well_supplied_count` = **nine** — Nuvalent class-forced hits whose forced class ENU supplies at ≥ 1 %
- `oh_atlas_coding_snvs` = **2,953** — sum of coding categories in Oh 2026 Table S2
- `oh_atlas_missense` = **2,173** — Missense rows in Oh 2026 Table S2
- `oh_atlas_noncoding` = **255** — sum of intron and splice-site rows in Oh 2026 Table S2
- `oh_atlas_nonsense` = **115** — Nonsense rows in Oh 2026 Table S2
- `oh_atlas_residue_max` = **1388** — maximum residue index in Oh 2026 Table S2 coding rows
- `oh_atlas_residue_min` = **1058** — minimum residue index in Oh 2026 Table S2 coding rows
- `oh_atlas_synonymous` = **665** — Synonymous rows in Oh 2026 Table S2
- `oh_atlas_total_snvs` = **3,208** — total rows in Oh 2026 Table S2 resistance scores
- `oh_class_forced_count` = **1,427** — class-forced ALK kinase-domain variants used in the ranking
- `perm_gain_quantum_range` = **9--33 pp** — smallest observable recall change, 1/n_hits, across anchors
- `perm_max_abs_gain` = **33 pp** — largest absolute framework-minus-sigma top-100 recall gain over all six anchor-proxy combinations
- `perm_mdg_range` = **9--50 pp** — 97.5th percentile of the permutation null's absolute top-100 recall gain, across all six anchor-proxy combinations
- `perm_n_hits_range` = **3 to 11** — class-forced observed hits per anchor
- `perm_observed_gain_range` = **-33 to 0 pp** — observed framework-minus-sigma top-100 recall gains across anchors
- `perm_p_awad_ada` = **0.60** — class-label permutation on the Awad clinical (KRAS Mut × ada) anchor
- `perm_p_awad_sot` = **0.30** — class-label permutation on the Awad clinical (KRAS Mut × sot) anchor
- `perm_p_doi_alec` = **1.00** — class-label permutation on the Doi 2022 (ALK Mutazyme × alec) anchor
- `perm_p_nuvalent_alec` = **1.00** — class-label permutation on the Nuvalent (ALK ENU × alec proxy) anchor
- `perm_p_nuvalent_lor` = **1.00** — class-label permutation on the Nuvalent (ALK ENU × lor proxy) anchor
- `perm_p_nuvalent_zot` = **1.00** — class-label permutation on the Nuvalent (ALK ENU × zot proxy) anchor
- `proxy_rho_alectinib` = **0.21** — Spearman correlation between alectinib atlas sigma and measured Ba/F3 neladalkib fold resistance; ranks the proxy choices
- `proxy_rho_lorlatinib` = **0.57** — Spearman correlation between lorlatinib atlas sigma and measured Ba/F3 neladalkib fold resistance; ranks the proxy choices
- `proxy_rho_zotizalkib` = **0.57** — Spearman correlation between zotizalkib atlas sigma and measured Ba/F3 neladalkib fold resistance; ranks the proxy choices
- `r68l_two_nt_routes` = **2** — two-substitution routes from AGG to a leucine codon
- `reach_be_blocked_alk` = **172 of 273** — ALK residues whose top-sigma class-forced substitution CBE+ABE cannot make, under alectinib
- `reach_be_blocked_kras` = **121 of 188** — KRAS G12C residues whose top-sigma class-forced substitution CBE+ABE cannot make, under adagrasib
- `reach_be_blocked_met` = **183 of 287** — MET residues whose top-sigma class-forced substitution CBE+ABE cannot make, under NVP-BVU972
- `reach_be_blocked_share_range` = **63\%--66\%** — range of the CBE+ABE blocked-residue share across all three genes and every drug arm
- `reach_class_forced_alk` = **1,427** — ALK class-forced variants per drug arm
- `reach_class_forced_kras` = **981** — KRAS G12C class-forced variants per drug arm
- `reach_class_forced_met` = **1,486** — MET class-forced variants per drug arm
- `reach_residues_alk` = **273** — ALK residues carrying at least one class-forced variant in the chemistry-free atlas
- `reach_residues_kras` = **188** — KRAS G12C residues carrying at least one class-forced variant in the chemistry-free atlas
- `reach_residues_met` = **287** — MET residues carrying at least one class-forced variant in the chemistry-free atlas
- `software_versions` = **Python 3.13.16 with pandas 3.0.5, numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.2, pyyaml 6.0.3** — emitted from the running interpreter, not recorded by hand
- `takei_v2_c_t_count` = **three** — Takei V2 partners class-forced into C>T (Mutazyme II's dominant class)
- `takei_v2_count` = **five** — V2 compound partners recovered at 500 nM neladalkib on I1171N fixed-V1 background
- `top100_alk_enu_alectinib` = **63\% [53\%--72\%]** — fraction of the chemistry-free top-100 that the ENU-simulated alectinib screen also ranks in its top 100, with a Wilson 95% interval
- `top100_alk_enu_lorlatinib` = **65\% [55\%--74\%]** — fraction of the chemistry-free top-100 that the ENU-simulated lorlatinib screen also ranks in its top 100, with a Wilson 95% interval
- `top100_alk_enu_zotizalkib` = **63\% [53\%--72\%]** — fraction of the chemistry-free top-100 that the ENU-simulated zotizalkib screen also ranks in its top 100, with a Wilson 95% interval
- `top100_alk_mutii_alectinib` = **65\% [55\%--74\%]** — fraction of the chemistry-free top-100 that the MutII-simulated alectinib screen also ranks in its top 100, with a Wilson 95% interval
- `top100_alk_mutii_lorlatinib` = **72\% [63\%--80\%]** — fraction of the chemistry-free top-100 that the MutII-simulated lorlatinib screen also ranks in its top 100, with a Wilson 95% interval
- `top100_alk_mutii_zotizalkib` = **65\% [55\%--74\%]** — fraction of the chemistry-free top-100 that the MutII-simulated zotizalkib screen also ranks in its top 100, with a Wilson 95% interval
- `top100_kras_enu_adagrasib` = **54\% [44\%--63\%]** — fraction of the chemistry-free top-100 that the ENU-simulated adagrasib screen also ranks in its top 100, with a Wilson 95% interval
- `top100_kras_enu_sotorasib` = **51\% [41\%--61\%]** — fraction of the chemistry-free top-100 that the ENU-simulated sotorasib screen also ranks in its top 100, with a Wilson 95% interval
- `top100_kras_mutii_adagrasib` = **40\% [31\%--50\%]** — fraction of the chemistry-free top-100 that the MutII-simulated adagrasib screen also ranks in its top 100, with a Wilson 95% interval
- `top100_kras_mutii_sotorasib` = **42\% [33\%--52\%]** — fraction of the chemistry-free top-100 that the MutII-simulated sotorasib screen also ranks in its top 100, with a Wilson 95% interval
- `transport_n_studies` = **7** — independent published studies contributing Ba/F3 ALK potency measurements
- `transport_rho_alectinib` = **0.51** — Spearman rank correlation between atlas sigma and measured Ba/F3 log10 fold resistance for alectinib, over 18 shared variants
- `transport_rho_lorlatinib` = **0.17** — Spearman rank correlation between atlas sigma and measured Ba/F3 log10 fold resistance for lorlatinib, over 18 shared variants
- `transport_rho_zotizalkib` = **0.82** — Spearman rank correlation between atlas sigma and measured Ba/F3 log10 fold resistance for zotizalkib, over 13 shared variants
- `twist_pe_classes_supplied` = **6 of 6** — Twist synthesis and prime editing produce all six classes approximately uniformly
- `two_hit_carriers_one_round` = **0.005--0.017** — expected carriers of R68L or V8M after a single mutagenesis round in a 1e7 library, across ENU and Mutazyme II
- `two_hit_carriers_twenty_rounds` = **1.1--3.6** — expected carriers of R68L or V8M after twenty rounds
- `two_hit_fold_suppression` = **51- to 169-fold** — suppression of R68L and V8M relative to a one-step variant of median class supply in the same library, ENU and Mutazyme II
- `two_hit_presence_one_round` = **1.7%** — largest presence probability for R68L or V8M after one round
- `v8m_two_nt_routes` = **1** — two-substitution routes from GTA to the methionine codon
- `y96d_sot_chemfree_rank` = **7** — Y96D chemistry-free rank under sotorasib
- `y96d_sot_mut_rank_sim` = **139** — Y96D simulated rank under sotorasib + Mutazyme II
- `y96d_sot_mut_shift` = **+132** — Y96D rank shift under sotorasib + Mutazyme II
- `y96d_sot_mutazyme_shift_beta_range` = **+5 to +643** — range of the Y96D rank shift under Mutazyme II with sotorasib across beta in [0.25, 4]

---

## Grouped by paper location

### Abstract

**`kras_panel_total`** = 12

- Value: `12`
- Formula: `|clinical ∪ case ∪ DMS|`
- Computed from: union of KRAS G12C clinical cohort + case report + accompanying in vitro DMS
- Inputs: `{'clinical': ['H95D', 'H95Q', 'H95R', 'R68S', 'Y96C'], 'case': ['Y96D'], 'dms': ['C12F', 'R68L', 'R68M', 'V8L', 'V8M', 'Y96H']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`kras_panel_dms`** = six

- Value: `6`
- Formula: `|DMS|`
- Computed from: variants from Awad 2021 accompanying in vitro deep mutational scan
- Inputs: `{'members': ['C12F', 'R68L', 'R68M', 'V8L', 'V8M', 'Y96H']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`kras_panel_clinical_total`** = six

- Value: `6`
- Formula: `|clinical ∪ case|`
- Computed from: patient-observed KRAS G12C variants (clinical cohort + case report)
- Inputs: `{'members': ['H95D', 'H95Q', 'H95R', 'R68S', 'Y96C', 'Y96D']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`kras_panel_snv_unreachable`** = two

- Value: `2`
- Formula: `count of panel variants whose Hamming distance from WT codon is ≥ 2 to all target codons`
- Computed from: codon reachability analysis in the G12C-mutant background
- Inputs: `{'members': ['R68L', 'V8M']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`doi2022_hits_count`** = six

- Value: `6`
- Formula: `|Doi Fig 3a resistance mutations|`
- Computed from: single-amino-acid resistance mutations in Doi 2022 Fig 3a, excluding silent passengers (L1152, L1165, R1192) and the I1171S sub-variant of the I1171N/T/S compound call
- Inputs: `{'members': ['G1202R', 'I1171N', 'I1171T', 'L1196M', 'L1198F', 'V1180L']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`met_two_screen_p`** = 0.029

- Value: `0.029463225369071463`
- Formula: `scipy.stats.chi2_contingency`
- Computed from: p-value of that homogeneity test; the two screens differ, so the class composition of a hit list is set by the drug, not the mutagen
- Inputs: `{'dof': 5}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`takei_v2_c_t_count`** = three

- Value: `3`
- Formula: `|{v in Takei V2: forced_class(v) == 'C>T'}|`
- Computed from: Takei V2 partners class-forced into C>T (Mutazyme II's dominant class)
- Inputs: `{'members': ['L1256F', 'E1129K', 'C1156Y']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`

**`nuvalent_cg_count`** = three

- Value: `3`
- Formula: `|{v in Nuvalent: forced_class(v) == 'C>G'}|`
- Computed from: Nuvalent hits class-forced into C>G (ENU 0.2 % class per germline spectrum)
- Inputs: `{'members': ['V1180L', 'T1151R', 'G1269A']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`


### Background

**`doi_anchor_topN_bound`** = top-120

- Value: `120`
- Formula: `ceil(max rank / 10) * 10`
- Computed from: smallest round bound containing all five class-forced Doi hits in the Mutazyme-simulated ranking
- Inputs: `{'max_rank': 118}`
- Code: `scripts/compute_paper_numbers.py::compute_doi_anchor_ranks`


### Discussion

**`two_hit_fold_suppression`** = 51- to 169-fold

- Value: `[50.65856129685916, 168.86187098953056]`
- Formula: `lambda_1nt(median class) / lambda_2nt`
- Computed from: suppression of R68L and V8M relative to a one-step variant of median class supply in the same library, ENU and Mutazyme II
- Inputs: `{'N': 10000000.0, 'mu': 0.0001, 'rounds': 20}`
- Code: `scripts/compute_paper_numbers.py::two_hit_reachability.py::main`

**`alk_largest_panel_shift`** = +597

- Value: `597`
- Formula: `max over panel variants of rank_simulated - rank_chemfree`
- Computed from: largest predicted rank shift among the ALK panel variants under Mutazyme II with alectinib
- Inputs: `{'variant': 'L1152R', 'forced_class': 'T>G'}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`

**`be_kras_panel_class_blocked`** = 7 of the 12

- Value: `7`
- Formula: `max(S_CBE+ABE(c)) over one-step routes equals zero`
- Computed from: KRAS G12C panel variants reachable in one nucleotide step whose every route falls in a class CBE+ABE cannot make
- Inputs: `{'variants': ['C12F', 'H95D', 'H95Q', 'R68M', 'R68S', 'V8L', 'Y96D']}`
- Code: `scripts/compute_paper_numbers.py::classify_amino_acid_change`

**`be_kras_panel_supplied`** = 3

- Value: `3`
- Formula: `max(S_CBE+ABE(c)) over one-step routes exceeds zero`
- Computed from: KRAS G12C panel variants a CBE+ABE pair can produce
- Inputs: `{'variants': ['H95R', 'Y96C', 'Y96H']}`
- Code: `scripts/compute_paper_numbers.py::classify_amino_acid_change`


### Discussion Limitation 4

**`h95d_ada_mut_rank_sim`** = 56

- Value: `56`
- Formula: `rank of log P(v) = log S(v) + log σ_adagrasib(v) under Mutazyme II spectrum`
- Computed from: H95D simulated rank under adagrasib + Mutazyme II (NGS-measured spectrum)
- Inputs: `{'variant': 'H95D', 'drug': 'adagrasib', 'mutagen': 'MutazymeII'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`


### Discussion §D7

**`takei_v2_count`** = five

- Value: `5`
- Formula: `|Takei V2 hits at 500 nM neladalkib|`
- Computed from: V2 compound partners recovered at 500 nM neladalkib on I1171N fixed-V1 background
- Inputs: `{'members': ['L1256F', 'E1129V', 'E1129K', 'C1156Y', 'T1151R']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`


### Fig 1B

**`kras_panel_total`** = 12

- Value: `12`
- Formula: `|clinical ∪ case ∪ DMS|`
- Computed from: union of KRAS G12C clinical cohort + case report + accompanying in vitro DMS
- Inputs: `{'clinical': ['H95D', 'H95Q', 'H95R', 'R68S', 'Y96C'], 'case': ['Y96D'], 'dms': ['C12F', 'R68L', 'R68M', 'V8L', 'V8M', 'Y96H']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`


### Fig 2A

**`enu_class_C_A`** = 9.0%

- Value: `0.09`
- Formula: `spectra.yaml:ENU_Arnold2012_nominal.values['C>A']`
- Computed from: Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- Inputs: `{'source': 'arnold2012', 'class': 'C>A'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`enu_class_C_G`** = 0.20%

- Value: `0.002`
- Formula: `spectra.yaml:ENU_Arnold2012_nominal.values['C>G']`
- Computed from: Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- Inputs: `{'source': 'arnold2012', 'class': 'C>G'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`enu_class_C_T`** = 17.9%

- Value: `0.179`
- Formula: `spectra.yaml:ENU_Arnold2012_nominal.values['C>T']`
- Computed from: Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- Inputs: `{'source': 'arnold2012', 'class': 'C>T'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`enu_class_T_A`** = 31.4%

- Value: `0.314`
- Formula: `spectra.yaml:ENU_Arnold2012_nominal.values['T>A']`
- Computed from: Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- Inputs: `{'source': 'arnold2012', 'class': 'T>A'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`enu_class_T_C`** = 36.6%

- Value: `0.366`
- Formula: `spectra.yaml:ENU_Arnold2012_nominal.values['T>C']`
- Computed from: Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- Inputs: `{'source': 'arnold2012', 'class': 'T>C'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`enu_class_T_G`** = 5.0%

- Value: `0.05`
- Formula: `spectra.yaml:ENU_Arnold2012_nominal.values['T>G']`
- Computed from: Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- Inputs: `{'source': 'arnold2012', 'class': 'T>G'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`mut_class_C_A`** = 18.0%

- Value: `0.18`
- Formula: `spectra.yaml:MutazymeII_DoiSupp2022.values['C>A']`
- Computed from: NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- Inputs: `{'source': 'doi2022', 'class': 'C>A'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`mut_class_C_G`** = 5.0%

- Value: `0.05`
- Formula: `spectra.yaml:MutazymeII_DoiSupp2022.values['C>G']`
- Computed from: NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- Inputs: `{'source': 'doi2022', 'class': 'C>G'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`mut_class_C_T`** = 47.0%

- Value: `0.47`
- Formula: `spectra.yaml:MutazymeII_DoiSupp2022.values['C>T']`
- Computed from: NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- Inputs: `{'source': 'doi2022', 'class': 'C>T'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`mut_class_T_A`** = 10.5%

- Value: `0.105`
- Formula: `spectra.yaml:MutazymeII_DoiSupp2022.values['T>A']`
- Computed from: NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- Inputs: `{'source': 'doi2022', 'class': 'T>A'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`mut_class_T_C`** = 18.0%

- Value: `0.18`
- Formula: `spectra.yaml:MutazymeII_DoiSupp2022.values['T>C']`
- Computed from: NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- Inputs: `{'source': 'doi2022', 'class': 'T>C'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`mut_class_T_G`** = 3.6%

- Value: `0.036`
- Formula: `spectra.yaml:MutazymeII_DoiSupp2022.values['T>G']`
- Computed from: NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- Inputs: `{'source': 'doi2022', 'class': 'T>G'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`


### Fig 2B

**`enu_classes_at_1pct`** = five of six

- Value: `5`
- Formula: `|{c: ENU[c] ≥ 0.01}|`
- Computed from: count of ENU classes supplied at ≥ 1 %; C>G is at 0.2 % and falls below
- Inputs: `{'ENU': '{...}', 'threshold': 0.01}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`mut_classes_at_1pct`** = 6 of 6

- Value: `6`
- Formula: `|{c: MutII[c] ≥ 0.01}|`
- Computed from: count of Mutazyme II classes supplied at ≥ 1 %; minimum class is T>G at 3.6 %
- Inputs: `{'MutII': '{...}', 'threshold': 0.01}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`be_classes_supplied`** = two of six

- Value: `2`
- Formula: `deterministic by base-editor chemistry`
- Computed from: CBE supplies C>T, ABE supplies T>C by design; the other 4 classes are 0
- Inputs: `{'supplied': ['C>T', 'T>C']}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`twist_pe_classes_supplied`** = 6 of 6

- Value: `6`
- Formula: `deterministic (chemistry-free library construction)`
- Computed from: Twist synthesis and prime editing produce all six classes approximately uniformly
- Inputs: `{'baseline_per_class': 0.16666666666666666}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`


### Fig 3A caption

**`kras_panel_total`** = 12

- Value: `12`
- Formula: `|clinical ∪ case ∪ DMS|`
- Computed from: union of KRAS G12C clinical cohort + case report + accompanying in vitro DMS
- Inputs: `{'clinical': ['H95D', 'H95Q', 'H95R', 'R68S', 'Y96C'], 'case': ['Y96D'], 'dms': ['C12F', 'R68L', 'R68M', 'V8L', 'V8M', 'Y96H']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`alk_panel_total`** = 12

- Value: `12`
- Formula: `|Doi hits ∪ Takei hits|`
- Computed from: union of alectinib Mutazyme II screen hits + next-gen ALK TKI screen hits
- Inputs: `{'doi': ['G1202R', 'I1171N', 'I1171T', 'L1196M', 'L1198F', 'V1180L'], 'takei': ['C1156Y', 'F1174C', 'F1174L', 'G1269A', 'L1152R', 'L1196Q']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`


### Fig 3B caption

**`oh_class_forced_count`** = 1,427

- Value: `1427`
- Formula: `|unique aa_change in Oh class-forced decomposition|`
- Computed from: class-forced ALK kinase-domain variants used in the ranking
- Inputs: `{'source': 'data/derived/oh_decomposition_all_drugs.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`


### Fig 3C

**`reach_residues_alk`** = 273

- Value: `273`
- Formula: `count of distinct residues in the class-forced set`
- Computed from: ALK residues carrying at least one class-forced variant in the chemistry-free atlas
- Inputs: `{'transcript': 'NM_004304.5'}`
- Code: `scripts/compute_paper_numbers.py::build_reachability_track.py::build`

**`reach_class_forced_alk`** = 1,427

- Value: `1427`
- Formula: `single-nucleotide variants whose every route falls in one class`
- Computed from: ALK class-forced variants per drug arm
- Inputs: `{'transcript': 'NM_004304.5'}`
- Code: `scripts/compute_paper_numbers.py::build_reachability_track.py::load_variants`

**`reach_be_blocked_alk`** = 172 of 273

- Value: `172`
- Formula: `count of residues where S_CBE+ABE of the top-sigma class-forced substitution is zero`
- Computed from: ALK residues whose top-sigma class-forced substitution CBE+ABE cannot make, under alectinib
- Inputs: `{'drug': 'alectinib', 'n_residues': 273, 'note': 'the arm Figure 3C plots'}`
- Code: `scripts/compute_paper_numbers.py::build_reachability_track.py::build`

**`reach_residues_kras`** = 188

- Value: `188`
- Formula: `count of distinct residues in the class-forced set`
- Computed from: KRAS G12C residues carrying at least one class-forced variant in the chemistry-free atlas
- Inputs: `{'transcript': 'NM_033360.4'}`
- Code: `scripts/compute_paper_numbers.py::build_reachability_track.py::build`

**`reach_class_forced_kras`** = 981

- Value: `981`
- Formula: `single-nucleotide variants whose every route falls in one class`
- Computed from: KRAS G12C class-forced variants per drug arm
- Inputs: `{'transcript': 'NM_033360.4'}`
- Code: `scripts/compute_paper_numbers.py::build_reachability_track.py::load_variants`

**`reach_be_blocked_kras`** = 121 of 188

- Value: `121`
- Formula: `count of residues where S_CBE+ABE of the top-sigma class-forced substitution is zero`
- Computed from: KRAS G12C residues whose top-sigma class-forced substitution CBE+ABE cannot make, under adagrasib
- Inputs: `{'drug': 'adagrasib', 'n_residues': 188, 'note': 'the arm Figure 3C plots'}`
- Code: `scripts/compute_paper_numbers.py::build_reachability_track.py::build`

**`reach_residues_met`** = 287

- Value: `287`
- Formula: `count of distinct residues in the class-forced set`
- Computed from: MET residues carrying at least one class-forced variant in the chemistry-free atlas
- Inputs: `{'transcript': 'NM_000245.2'}`
- Code: `scripts/compute_paper_numbers.py::build_reachability_track.py::build`

**`reach_class_forced_met`** = 1,486

- Value: `1486`
- Formula: `single-nucleotide variants whose every route falls in one class`
- Computed from: MET class-forced variants per drug arm
- Inputs: `{'transcript': 'NM_000245.2'}`
- Code: `scripts/compute_paper_numbers.py::build_reachability_track.py::load_variants`

**`reach_be_blocked_met`** = 183 of 287

- Value: `183`
- Formula: `count of residues where S_CBE+ABE of the top-sigma class-forced substitution is zero`
- Computed from: MET residues whose top-sigma class-forced substitution CBE+ABE cannot make, under NVP-BVU972
- Inputs: `{'drug': 'NVP-BVU972', 'n_residues': 287, 'note': 'the arm Figure 3C plots'}`
- Code: `scripts/compute_paper_numbers.py::build_reachability_track.py::build`

**`reach_be_blocked_share_range`** = 63\%--66\%

- Value: `[0.6276595744680851, 0.6585365853658537]`
- Formula: `min and max over (gene, drug) of blocked residues / residues`
- Computed from: range of the CBE+ABE blocked-residue share across all three genes and every drug arm
- Inputs: `{'n_arms': 7}`
- Code: `scripts/compute_paper_numbers.py::build_reachability_track.py::build`


### Fig 4A caption

**`doi2022_hits_count`** = six

- Value: `6`
- Formula: `|Doi Fig 3a resistance mutations|`
- Computed from: single-amino-acid resistance mutations in Doi 2022 Fig 3a, excluding silent passengers (L1152, L1165, R1192) and the I1171S sub-variant of the I1171N/T/S compound call
- Inputs: `{'members': ['G1202R', 'I1171N', 'I1171T', 'L1196M', 'L1198F', 'V1180L']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`doi2022_hits_class_forced`** = five

- Value: `5`
- Formula: `|Doi hits ∩ class-forced|`
- Computed from: Doi 2022 hits that are class-forced (G1202R is multi-class)
- Inputs: `{'members': ['I1171N', 'I1171T', 'L1196M', 'L1198F', 'V1180L']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`


### Fig 4B caption

**`takei_v2_count`** = five

- Value: `5`
- Formula: `|Takei V2 hits at 500 nM neladalkib|`
- Computed from: V2 compound partners recovered at 500 nM neladalkib on I1171N fixed-V1 background
- Inputs: `{'members': ['L1256F', 'E1129V', 'E1129K', 'C1156Y', 'T1151R']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`

**`takei_v2_c_t_count`** = three

- Value: `3`
- Formula: `|{v in Takei V2: forced_class(v) == 'C>T'}|`
- Computed from: Takei V2 partners class-forced into C>T (Mutazyme II's dominant class)
- Inputs: `{'members': ['L1256F', 'E1129K', 'C1156Y']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`


### Fig 4C caption

**`enu_well_supplied_well_vs_total`** = nine of twelve

- Value: `9 of 12`
- Formula: `|{hit: ENU[hit.class] ≥ 0.01, hit in Nuvalent class-forced pool}|`
- Computed from: count of Nuvalent class-forced hits whose forced class ENU supplies at ≥ 1 %
- Inputs: `{'nuvalent_class_forced_total': 12, 'threshold': 0.01}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`nuvalent_total`** = thirteen

- Value: `13`
- Formula: `|Nuvalent single-amino-acid hits (aggregate)|`
- Computed from: single-amino-acid variants reported in Nuvalent AACR 2025 poster aggregate across arms
- Inputs: `{'members': ['I1171N', 'I1171S', 'I1171T', 'F1174C', 'F1174V', 'V1180L', 'T1151K', 'T1151R', 'E1129K', 'R1192P', 'L1196M', 'G1269A', 'G1202R']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`

**`nuvalent_class_forced`** = twelve

- Value: `11`
- Formula: `|Nuvalent hits ∩ class-forced|`
- Computed from: Nuvalent hits that are class-forced single-nt reachable (R1192P outside atlas, G1202R multi-class)
- Inputs: `{'members': ['I1171N', 'I1171S', 'I1171T', 'F1174C', 'F1174V', 'V1180L', 'T1151K', 'T1151R', 'E1129K', 'L1196M', 'G1269A']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`

**`nuvalent_cg_count`** = three

- Value: `3`
- Formula: `|{v in Nuvalent: forced_class(v) == 'C>G'}|`
- Computed from: Nuvalent hits class-forced into C>G (ENU 0.2 % class per germline spectrum)
- Inputs: `{'members': ['V1180L', 'T1151R', 'G1269A']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`

**`nuvalent_well_supplied_count`** = nine

- Value: `8`
- Formula: `|Nuvalent class-forced| - |Nuvalent C>G|`
- Computed from: Nuvalent class-forced hits whose forced class ENU supplies at ≥ 1 %
- Inputs: `{'members': ['I1171N', 'I1171S', 'I1171T', 'F1174C', 'F1174V', 'T1151K', 'E1129K', 'L1196M']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`


### Fig 5A

**`met_two_screen_chi2`** = 12.4

- Value: `12.420123432255783`
- Formula: `scipy.stats.chi2_contingency(2 x 6 clone counts)`
- Computed from: chi-square homogeneity of the six-class clone composition between the NVP-BVU972 and AMG 458 screens, which used the same ENU protocol in the same cells
- Inputs: `{'dof': 5, 'p': 0.029463225369071463, 'table': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_two_screen_p`** = 0.029

- Value: `0.029463225369071463`
- Formula: `scipy.stats.chi2_contingency`
- Computed from: p-value of that homogeneity test; the two screens differ, so the class composition of a hit list is set by the drug, not the mutagen
- Inputs: `{'dof': 5}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_nvp_bvu972_clones`** = 85

- Value: `85`
- Formula: `sum of clone counts in Tiedt Supplementary Table`
- Computed from: resistant clones reported for the NVP-BVU972 screen
- Inputs: `{'n_base_changes': 10, 'n_aa_changes': 10}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_amg_clones`** = 60

- Value: `60`
- Formula: `sum of clone counts in Tiedt Supplementary Table`
- Computed from: resistant clones reported for the AMG 458 screen
- Inputs: `{'n_base_changes': 15, 'n_aa_changes': 13}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_amg_ta_share`** = 33\%

- Value: `0.3333333333333333`
- Formula: `clones in class / total clones`
- Computed from: share of AMG 458 clones whose substitution is T>A
- Inputs: `{'n_clones': 60}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`


### Fig 5B

**`met_y1230_recovered`** = 5 of 6

- Value: `5`
- Formula: `enumerate one-step routes from the codon; intersect with Tiedt Table S1`
- Computed from: single-nucleotide missense routes out of the Y1230 wild-type codon that the NVP-BVU972 screen recovered
- Inputs: `{'wt_codon': 'TAT', 'n_clones': 75}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_supply_fold`** = 7.3-fold

- Value: `7.319999999999999`
- Formula: `max(S_ENU(c)) / min(S_ENU(c)) over the routes' classes`
- Computed from: ratio of the largest to smallest germline ENU supply among the substitution classes those Y1230 routes span
- Inputs: `{'classes': ['T>A', 'T>C', 'T>G'], 'lo': 0.05, 'hi': 0.366}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_clones_y1230c`** = 31

- Value: `31`
- Formula: `transcribed clone count`
- Computed from: clones carrying Y1230C in the Tiedt NVP-BVU972 screen
- Inputs: `{'classes': ['T>C']}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_clones_y1230d`** = 1

- Value: `1`
- Formula: `transcribed clone count`
- Computed from: clones carrying Y1230D in the Tiedt NVP-BVU972 screen
- Inputs: `{'classes': ['T>G']}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_clones_y1230f`** = 0

- Value: `0`
- Formula: `transcribed clone count`
- Computed from: clones carrying Y1230F in the Tiedt NVP-BVU972 screen
- Inputs: `{'classes': ['T>A']}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_clones_y1230h`** = 20

- Value: `20`
- Formula: `transcribed clone count`
- Computed from: clones carrying Y1230H in the Tiedt NVP-BVU972 screen
- Inputs: `{'classes': ['T>C']}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_clones_y1230n`** = 19

- Value: `19`
- Formula: `transcribed clone count`
- Computed from: clones carrying Y1230N in the Tiedt NVP-BVU972 screen
- Inputs: `{'classes': ['T>A']}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_clones_y1230s`** = 4

- Value: `4`
- Formula: `transcribed clone count`
- Computed from: clones carrying Y1230S in the Tiedt NVP-BVU972 screen
- Inputs: `{'classes': ['T>G']}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`


### Fig 5C

**`met_nvp_bvu972_twostep_top100`** = 72

- Value: `72`
- Formula: `classify each top-100 variant by one-step routes from its codon`
- Computed from: variants in the chemistry-free top 100 for NVP-BVU972 that need two or more nucleotide changes and so lie outside one-step reach
- Inputs: `{'two-step': 72, 'class-forced': 26, 'multi-class SNV': 2}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_amg_twostep_top100`** = 75

- Value: `75`
- Formula: `classify each top-100 variant by one-step routes from its codon`
- Computed from: variants in the chemistry-free top 100 for AMG 458 that need two or more nucleotide changes and so lie outside one-step reach
- Inputs: `{'class-forced': 23, 'two-step': 75, 'multi-class SNV': 2}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`


### Methods

**`oh_atlas_total_snvs`** = 3,208

- Value: `3208`
- Formula: `len(data/raw/oh_2026/table_s2_resistance_scores.csv)`
- Computed from: total rows in Oh 2026 Table S2 resistance scores
- Inputs: `{'path': 'data/raw/oh_2026/table_s2_resistance_scores.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`oh_atlas_coding_snvs`** = 2,953

- Value: `2953`
- Formula: `sum(Missense + Synonymous + Nonsense) in table_s2_resistance_scores.csv`
- Computed from: sum of coding categories in Oh 2026 Table S2
- Inputs: `{'by_type': '{...}', 'coding_types': ['Missense', 'Nonsense', 'Synonymous']}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`oh_atlas_missense`** = 2,173

- Value: `2173`
- Formula: `count of rows where type == 'Missense'`
- Computed from: Missense rows in Oh 2026 Table S2
- Inputs: `{'path': 'data/raw/oh_2026/table_s2_resistance_scores.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`oh_atlas_synonymous`** = 665

- Value: `665`
- Formula: `count of rows where type == 'Synonymous'`
- Computed from: Synonymous rows in Oh 2026 Table S2
- Inputs: `{'path': 'data/raw/oh_2026/table_s2_resistance_scores.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`oh_atlas_nonsense`** = 115

- Value: `115`
- Formula: `count of rows where type == 'Nonsense'`
- Computed from: Nonsense rows in Oh 2026 Table S2
- Inputs: `{'path': 'data/raw/oh_2026/table_s2_resistance_scores.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`oh_atlas_noncoding`** = 255

- Value: `255`
- Formula: `count(type == 'Intron') + count(type == 'Splice Acceptor/Donor')`
- Computed from: sum of intron and splice-site rows in Oh 2026 Table S2
- Inputs: `{'intron': 149, 'splice': 106}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`oh_atlas_residue_min`** = 1058

- Value: `1058`
- Formula: `min(residue) over rows with type != 'Intron'`
- Computed from: minimum residue index in Oh 2026 Table S2 coding rows
- Inputs: `{'path': 'data/raw/oh_2026/table_s2_resistance_scores.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`oh_atlas_residue_max`** = 1388

- Value: `1388`
- Formula: `max(residue) over rows with type != 'Intron'`
- Computed from: maximum residue index in Oh 2026 Table S2 coding rows
- Inputs: `{'path': 'data/raw/oh_2026/table_s2_resistance_scores.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`feng_class_forced_count`** = 981

- Value: `981`
- Formula: `|{v in Feng: v is class-forced single-nt reachable}|`
- Computed from: class-forced KRAS G12C variants in the Feng atlas (unique amino-acid changes)
- Inputs: `{'source': 'data/derived/feng_decomposition.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`oh_class_forced_count`** = 1,427

- Value: `1427`
- Formula: `|unique aa_change in Oh class-forced decomposition|`
- Computed from: class-forced ALK kinase-domain variants used in the ranking
- Inputs: `{'source': 'data/derived/oh_decomposition_all_drugs.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`takei_v2_count`** = five

- Value: `5`
- Formula: `|Takei V2 hits at 500 nM neladalkib|`
- Computed from: V2 compound partners recovered at 500 nM neladalkib on I1171N fixed-V1 background
- Inputs: `{'members': ['L1256F', 'E1129V', 'E1129K', 'C1156Y', 'T1151R']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`

**`nuvalent_total`** = thirteen

- Value: `13`
- Formula: `|Nuvalent single-amino-acid hits (aggregate)|`
- Computed from: single-amino-acid variants reported in Nuvalent AACR 2025 poster aggregate across arms
- Inputs: `{'members': ['I1171N', 'I1171S', 'I1171T', 'F1174C', 'F1174V', 'V1180L', 'T1151K', 'T1151R', 'E1129K', 'R1192P', 'L1196M', 'G1269A', 'G1202R']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`


### Methods external test of sigma transport

**`transport_n_studies`** = 7

- Value: `7`
- Formula: `nunique(paper_id)`
- Computed from: independent published studies contributing Ba/F3 ALK potency measurements
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::load_potency`


### Methods library depth

**`cg_detect_prob_germline`** = 0.86

- Value: `0.8646647167633873`
- Formula: `1 - exp(-N*mu*S*r)`
- Computed from: detection probability of a C>G variant at the germline ENU rate in a 1e7 Ba/F3 library; shows a rare class is still supplied
- Inputs: `{'N': 10000000.0, 'mu': 0.0001, 'S': 0.002, 'routes': 1}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`


### Methods multi-nucleotide reachability

**`two_hit_fold_suppression`** = 51- to 169-fold

- Value: `[50.65856129685916, 168.86187098953056]`
- Formula: `lambda_1nt(median class) / lambda_2nt`
- Computed from: suppression of R68L and V8M relative to a one-step variant of median class supply in the same library, ENU and Mutazyme II
- Inputs: `{'N': 10000000.0, 'mu': 0.0001, 'rounds': 20}`
- Code: `scripts/compute_paper_numbers.py::two_hit_reachability.py::main`


### Methods software and reproducibility

**`software_versions`** = Python 3.13.16 with pandas 3.0.5, numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.2, pyyaml 6.0.3

- Value: `{'python': '3.13.16', 'implementation': 'CPython', 'packages': {'pandas': '3.0.5', 'numpy': '2.5.3', 'scipy': '1.18.1', 'matplotlib': '3.11.2', 'pyyaml': '6.0.3'}}`
- Formula: `platform.python_version() and each package __version__`
- Computed from: emitted from the running interpreter, not recorded by hand
- Code: `scripts/compute_paper_numbers.py::env_versions.py::main`


### Methods §Resistance panels

**`kras_panel_clinical_cohort`** = five

- Value: `5`
- Formula: `|clinical|`
- Computed from: patient-observed variants from Awad 2021 clinical cohort
- Inputs: `{'members': ['H95D', 'H95Q', 'H95R', 'R68S', 'Y96C']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`kras_panel_case_report`** = one

- Value: `1`
- Formula: `|case|`
- Computed from: patient-observed Y96D from Tanaka 2021 case report
- Inputs: `{'members': ['Y96D']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`kras_panel_dms`** = six

- Value: `6`
- Formula: `|DMS|`
- Computed from: variants from Awad 2021 accompanying in vitro deep mutational scan
- Inputs: `{'members': ['C12F', 'R68L', 'R68M', 'V8L', 'V8M', 'Y96H']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`


### Table 1 CBE+ABE row

**`be_classes_supplied`** = two of six

- Value: `2`
- Formula: `deterministic by base-editor chemistry`
- Computed from: CBE supplies C>T, ABE supplies T>C by design; the other 4 classes are 0
- Inputs: `{'supplied': ['C>T', 'T>C']}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`


### Table 1 ENU row

**`enu_classes_at_1pct`** = five of six

- Value: `5`
- Formula: `|{c: ENU[c] ≥ 0.01}|`
- Computed from: count of ENU classes supplied at ≥ 1 %; C>G is at 0.2 % and falls below
- Inputs: `{'ENU': '{...}', 'threshold': 0.01}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`


### Table 1 Mutazyme row

**`mut_classes_at_1pct`** = 6 of 6

- Value: `6`
- Formula: `|{c: MutII[c] ≥ 0.01}|`
- Computed from: count of Mutazyme II classes supplied at ≥ 1 %; minimum class is T>G at 3.6 %
- Inputs: `{'MutII': '{...}', 'threshold': 0.01}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`


### Table 1 PE row

**`twist_pe_classes_supplied`** = 6 of 6

- Value: `6`
- Formula: `deterministic (chemistry-free library construction)`
- Computed from: Twist synthesis and prime editing produce all six classes approximately uniformly
- Inputs: `{'baseline_per_class': 0.16666666666666666}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`


### Table 1 Twist row

**`twist_pe_classes_supplied`** = 6 of 6

- Value: `6`
- Formula: `deterministic (chemistry-free library construction)`
- Computed from: Twist synthesis and prime editing produce all six classes approximately uniformly
- Inputs: `{'baseline_per_class': 0.16666666666666666}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`


### Table 2

**`perm_mdg_range`** = 9--50 pp

- Value: `[0.0909090909090909, 0.5]`
- Formula: `quantile(|null gain|, 0.975) per anchor`
- Computed from: 97.5th percentile of the permutation null's absolute top-100 recall gain, across all six anchor-proxy combinations
- Inputs: `{'n_iter': 5000, 'anchors': 6}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`

**`perm_p_doi_alec`** = 1.00

- Value: `1.0`
- Formula: `(#{|null gain| >= |observed gain|} + 1) / (n_iter + 1)`
- Computed from: class-label permutation on the Doi 2022 (ALK Mutazyme × alec) anchor
- Inputs: `{'n_iter': 5000, 'K': 100, 'n_hits': 5}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`

**`perm_p_nuvalent_alec`** = 1.00

- Value: `1.0`
- Formula: `(#{|null gain| >= |observed gain|} + 1) / (n_iter + 1)`
- Computed from: class-label permutation on the Nuvalent (ALK ENU × alec proxy) anchor
- Inputs: `{'n_iter': 5000, 'K': 100, 'n_hits': 11}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`

**`perm_p_nuvalent_lor`** = 1.00

- Value: `1.0`
- Formula: `(#{|null gain| >= |observed gain|} + 1) / (n_iter + 1)`
- Computed from: class-label permutation on the Nuvalent (ALK ENU × lor proxy) anchor
- Inputs: `{'n_iter': 5000, 'K': 100, 'n_hits': 11}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`

**`perm_p_nuvalent_zot`** = 1.00

- Value: `1.0`
- Formula: `(#{|null gain| >= |observed gain|} + 1) / (n_iter + 1)`
- Computed from: class-label permutation on the Nuvalent (ALK ENU × zot proxy) anchor
- Inputs: `{'n_iter': 5000, 'K': 100, 'n_hits': 11}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`

**`perm_p_awad_sot`** = 0.30

- Value: `0.2999400119976004`
- Formula: `(#{|null gain| >= |observed gain|} + 1) / (n_iter + 1)`
- Computed from: class-label permutation on the Awad clinical (KRAS Mut × sot) anchor
- Inputs: `{'n_iter': 5000, 'K': 100, 'n_hits': 3}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`

**`perm_p_awad_ada`** = 0.60

- Value: `0.6012797440511898`
- Formula: `(#{|null gain| >= |observed gain|} + 1) / (n_iter + 1)`
- Computed from: class-label permutation on the Awad clinical (KRAS Mut × ada) anchor
- Inputs: `{'n_iter': 5000, 'K': 100, 'n_hits': 4}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`


### Table 3 (N0 row)

**`nuvalent_cg_binom_p_raw`** = 2.3e-06

- Value: `2.253926009047811e-06`
- Formula: `binom.sf(k-1, n, S_ENU('C>G'))`
- Computed from: exact binomial tail against the raw germline spectrum; retained only as the confounded comparison
- Inputs: `{'k': 3, 'n': 13, 'p': 0.002}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`


### Table 3 (N2 row)

**`nuvalent_cg_top100_sigma`** = 24 of the

- Value: `24`
- Formula: `count(forced_class == 'C>G') over nlargest(100, pooled sigma)`
- Computed from: C>G count among the top-100 sigma-ranked class-forced ALK variants, pooled across drug arms
- Inputs: `{'N_top': 100}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::pi_sigma_topn`

**`nuvalent_cg_binom_p_selection`** = 0.64

- Value: `0.636406728920353`
- Formula: `binom.sf(k-1, n, pi_N2)`
- Computed from: binomial tail of 3-of-13 C>G against a selection-only null; the count is unremarkable once supply is removed
- Inputs: `{'k': 3, 'n': 13, 'pi': 0.24}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`


### Table 3 (N3 row)

**`nuvalent_cg_pi_model`** = 0.56%

- Value: `0.00556039021953432`
- Formula: `sum_{v in C>G} S(c)exp(beta*sigma_v) / sum_v S(c)exp(beta*sigma_v)`
- Computed from: model-implied probability that a hit is C>G, with germline supply and pooled sigma at beta=1
- Inputs: `{'beta': 1.0, 'sigma': 'max over alec/lor/zot', 'S': 'ENU germline'}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::pi_model`

**`nuvalent_cg_binom_p_model`** = 4.7e-05

- Value: `4.715801648246043e-05`
- Formula: `binom.sf(k-1, n, pi_N3)`
- Computed from: binomial tail of 3-of-13 C>G against the model-implied null; the germline spectrum is rejected
- Inputs: `{'k': 3, 'n': 13, 'pi': 0.00556039021953432}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`


### Table 3 (N4 row)

**`nuvalent_cg_binom_p_composition`** = 0.39

- Value: `0.38585678830594816`
- Formula: `binom.sf(k-1, n, n_CG / n_classforced)`
- Computed from: binomial tail against atlas class composition alone
- Inputs: `{'n_CG': 243, 'n': 1427}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`


### Table 4

**`transport_rho_alectinib`** = 0.51

- Value: `0.5149638802889577`
- Formula: `spearmanr(sigma_atlas, median log10(IC50_variant / IC50_WT) within study)`
- Computed from: Spearman rank correlation between atlas sigma and measured Ba/F3 log10 fold resistance for alectinib, over 18 shared variants
- Inputs: `{'n': 18, 'boot95': [-0.04275959677304641, 0.855497382198953], 'p': 0.028748754882630345, 'median_inter_study_spread_log10': 0.12928740120303056}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`transport_rho_lorlatinib`** = 0.17

- Value: `0.16615067079463364`
- Formula: `spearmanr(sigma_atlas, median log10(IC50_variant / IC50_WT) within study)`
- Computed from: Spearman rank correlation between atlas sigma and measured Ba/F3 log10 fold resistance for lorlatinib, over 18 shared variants
- Inputs: `{'n': 18, 'boot95': [-0.42618532491530653, 0.7229182515739769], 'p': 0.5099504603438955, 'median_inter_study_spread_log10': 0.3219465077435653}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`transport_rho_zotizalkib`** = 0.82

- Value: `0.8241758241758241`
- Formula: `spearmanr(sigma_atlas, median log10(IC50_variant / IC50_WT) within study)`
- Computed from: Spearman rank correlation between atlas sigma and measured Ba/F3 log10 fold resistance for zotizalkib, over 13 shared variants
- Inputs: `{'n': 13, 'boot95': [0.4871720116618075, 0.9602272727272727], 'p': 0.0005302228375973589, 'median_inter_study_spread_log10': 0.0}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`atlas_dose_ratio_alectinib`** = 1.7x

- Value: `1.7241379310344829`
- Formula: `atlas_dose_nM / median(WT Ba/F3 IC50_nM)`
- Computed from: atlas assay concentration for alectinib divided by the measured Ba/F3 wild-type IC50
- Inputs: `{'dose_nM': 20, 'wt_ic50_nM': 11.6}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`atlas_dose_ratio_lorlatinib`** = 3.7x

- Value: `3.7037037037037033`
- Formula: `atlas_dose_nM / median(WT Ba/F3 IC50_nM)`
- Computed from: atlas assay concentration for lorlatinib divided by the measured Ba/F3 wild-type IC50
- Inputs: `{'dose_nM': 5, 'wt_ic50_nM': 1.35}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`atlas_dose_ratio_zotizalkib`** = 30.6x

- Value: `30.612244897959187`
- Formula: `atlas_dose_nM / median(WT Ba/F3 IC50_nM)`
- Computed from: atlas assay concentration for zotizalkib divided by the measured Ba/F3 wild-type IC50
- Inputs: `{'dose_nM': 300, 'wt_ic50_nM': 9.799999999999999}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`lorlatinib_inter_study_spread`** = 0.32

- Value: `0.3219465077435653`
- Formula: `median over variants of max - min log10 fold change`
- Computed from: median inter-study range of measured lorlatinib potency on the same variant
- Inputs: `{'n': 18}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py`


### Table S1

**`kras_panel_total`** = 12

- Value: `12`
- Formula: `|clinical ∪ case ∪ DMS|`
- Computed from: union of KRAS G12C clinical cohort + case report + accompanying in vitro DMS
- Inputs: `{'clinical': ['H95D', 'H95Q', 'H95R', 'R68S', 'Y96C'], 'case': ['Y96D'], 'dms': ['C12F', 'R68L', 'R68M', 'V8L', 'V8M', 'Y96H']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`kras_panel_clinical_cohort`** = five

- Value: `5`
- Formula: `|clinical|`
- Computed from: patient-observed variants from Awad 2021 clinical cohort
- Inputs: `{'members': ['H95D', 'H95Q', 'H95R', 'R68S', 'Y96C']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`kras_panel_dms`** = six

- Value: `6`
- Formula: `|DMS|`
- Computed from: variants from Awad 2021 accompanying in vitro deep mutational scan
- Inputs: `{'members': ['C12F', 'R68L', 'R68M', 'V8L', 'V8M', 'Y96H']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`alk_panel_total`** = 12

- Value: `12`
- Formula: `|Doi hits ∪ Takei hits|`
- Computed from: union of alectinib Mutazyme II screen hits + next-gen ALK TKI screen hits
- Inputs: `{'doi': ['G1202R', 'I1171N', 'I1171T', 'L1196M', 'L1198F', 'V1180L'], 'takei': ['C1156Y', 'F1174C', 'F1174L', 'G1269A', 'L1152R', 'L1196Q']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`


### Table S1 summary

**`alk_panel_class_forced`** = 10

- Value: `10`
- Formula: `count of ALK panel variants reachable via exactly one pyrimidine class`
- Computed from: codon reachability of ALK panel
- Inputs: `{'members': ['C1156Y', 'F1174C', 'G1269A', 'I1171N', 'I1171T', 'L1152R', 'L1196M', 'L1196Q', 'L1198F', 'V1180L']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`alk_panel_multi_class`** = two

- Value: `2`
- Formula: `count of ALK panel variants with multi-class single-nt reachability`
- Computed from: codon reachability of ALK panel
- Inputs: `{'members': ['F1174L', 'G1202R']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`


### not reported individually: §3.4 quotes the range

**`top100_kras_enu_adagrasib`** = 54\% [44\%--63\%]

- Value: `0.54`
- Formula: `|sim top-100 ∩ chemistry-free top-100| / 100`
- Computed from: fraction of the chemistry-free top-100 that the ENU-simulated adagrasib screen also ranks in its top 100, with a Wilson 95% interval
- Inputs: `{'gene': 'KRAS G12C', 'N': 981, 'K': 100}`
- Code: `scripts/compute_paper_numbers.py::recompute_topk_with_ci.py::compute_curve`

**`top100_kras_mutii_adagrasib`** = 40\% [31\%--50\%]

- Value: `0.4`
- Formula: `|sim top-100 ∩ chemistry-free top-100| / 100`
- Computed from: fraction of the chemistry-free top-100 that the MutII-simulated adagrasib screen also ranks in its top 100, with a Wilson 95% interval
- Inputs: `{'gene': 'KRAS G12C', 'N': 981, 'K': 100}`
- Code: `scripts/compute_paper_numbers.py::recompute_topk_with_ci.py::compute_curve`

**`top100_kras_enu_sotorasib`** = 51\% [41\%--61\%]

- Value: `0.51`
- Formula: `|sim top-100 ∩ chemistry-free top-100| / 100`
- Computed from: fraction of the chemistry-free top-100 that the ENU-simulated sotorasib screen also ranks in its top 100, with a Wilson 95% interval
- Inputs: `{'gene': 'KRAS G12C', 'N': 981, 'K': 100}`
- Code: `scripts/compute_paper_numbers.py::recompute_topk_with_ci.py::compute_curve`

**`top100_kras_mutii_sotorasib`** = 42\% [33\%--52\%]

- Value: `0.42`
- Formula: `|sim top-100 ∩ chemistry-free top-100| / 100`
- Computed from: fraction of the chemistry-free top-100 that the MutII-simulated sotorasib screen also ranks in its top 100, with a Wilson 95% interval
- Inputs: `{'gene': 'KRAS G12C', 'N': 981, 'K': 100}`
- Code: `scripts/compute_paper_numbers.py::recompute_topk_with_ci.py::compute_curve`


### not reported: calibrated ENU rows were removed from §2.3

**`h95d_ada_empirical_lower_shift_beta_range`** = +2 to +518

- Value: `+2 to +518`
- Formula: `[min_β shift, max_β shift]`
- Computed from: H95D shift interval across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU empirical lower (C>G = 7%)
- Inputs: `{'beta_grid': [0.25, 0.5, 1.0, 2.0, 4.0], 'per_beta': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_empirical_mean_shift_beta_range`** = -4 to -3

- Value: `-4 to -3`
- Formula: `[min_β shift, max_β shift]`
- Computed from: H95D shift interval across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU empirical mean (C>G = 25%)
- Inputs: `{'beta_grid': [0.25, 0.5, 1.0, 2.0, 4.0], 'per_beta': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_mutazyme_shift_beta_range`** = +1 to +606

- Value: `+1 to +606`
- Formula: `[min_β shift, max_β shift]`
- Computed from: H95D shift interval across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for Mutazyme II (NGS; C>G = 5%)
- Inputs: `{'beta_grid': [0.25, 0.5, 1.0, 2.0, 4.0], 'per_beta': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`


### not reported: calibrated ENU rows were removed from §2.3 because §2.5 shows the inversion is not identified

**`h95d_ada_enu_empirical_lower_rank_sim`** = 26

- Value: `26`
- Formula: `rank of log P(v) = log S(v) + log σ_adagrasib(v), descending, over 1,205 class-forced variants`
- Computed from: H95D simulated rank under adagrasib + ENU with C>G = 7.0% (Jeffreys 95% lower bound from 3-of-13 anomaly)
- Inputs: `{'variant': 'H95D', 'drug': 'adagrasib', 'mutagen': 'ENU', 'CG_rate': 0.07, 'spec': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_enu_empirical_lower_shift`** = +17

- Value: `17`
- Formula: `rank_sim(ENU, p_CG) - rank_chemfree`
- Computed from: H95D rank shift under adagrasib + ENU (Jeffreys 95% lower bound from 3-of-13 anomaly)
- Inputs: `{'rank_sim': 26, 'rank_cf': 9, 'p_CG': 0.07}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_enu_empirical_mean_rank_sim`** = 5

- Value: `5`
- Formula: `rank of log P(v) = log S(v) + log σ_adagrasib(v), descending, over 1,205 class-forced variants`
- Computed from: H95D simulated rank under adagrasib + ENU with C>G = 25.0% (Jeffreys posterior mean from 3-of-13 anomaly)
- Inputs: `{'variant': 'H95D', 'drug': 'adagrasib', 'mutagen': 'ENU', 'CG_rate': 0.25, 'spec': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_enu_empirical_mean_shift`** = -4

- Value: `-4`
- Formula: `rank_sim(ENU, p_CG) - rank_chemfree`
- Computed from: H95D rank shift under adagrasib + ENU (Jeffreys posterior mean from 3-of-13 anomaly)
- Inputs: `{'rank_sim': 5, 'rank_cf': 9, 'p_CG': 0.25}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_enu_empirical_upper_rank_sim`** = 2

- Value: `2`
- Formula: `rank of log P(v) = log S(v) + log σ_adagrasib(v), descending, over 1,205 class-forced variants`
- Computed from: H95D simulated rank under adagrasib + ENU with C>G = 49.7% (Jeffreys 95% upper bound from 3-of-13 anomaly)
- Inputs: `{'variant': 'H95D', 'drug': 'adagrasib', 'mutagen': 'ENU', 'CG_rate': 0.497, 'spec': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_enu_empirical_upper_shift`** = -7

- Value: `-7`
- Formula: `rank_sim(ENU, p_CG) - rank_chemfree`
- Computed from: H95D rank shift under adagrasib + ENU (Jeffreys 95% upper bound from 3-of-13 anomaly)
- Inputs: `{'rank_sim': 2, 'rank_cf': 9, 'p_CG': 0.497}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`


### not reported: the variant is named but its score is not

**`alk_c1156y_sigma`** = +1.1

- Value: `1.11412523715659`
- Formula: `sigma = Resistance score, chemistry-free`
- Computed from: Oh atlas alectinib Resistance score for C1156Y
- Inputs: `{'forced_class': 'C>T'}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`

**`alk_l1152r_sigma`** = +0.7

- Value: `0.704426250435087`
- Formula: `sigma = Resistance score, chemistry-free`
- Computed from: Oh atlas alectinib Resistance score for L1152R
- Inputs: `{'forced_class': 'T>G'}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`


### supporting: §2.2 states this as a presence probability below 2%

**`two_hit_presence_one_round`** = 1.7%

- Value: `0.0167776607241758`
- Formula: `1 - exp(-lambda)`
- Computed from: largest presence probability for R68L or V8M after one round
- Inputs: `{'N': 10000000.0, 'mu': 0.0001, 'rounds': 1}`
- Code: `scripts/compute_paper_numbers.py::two_hit_reachability.py::main`


### supporting: §2.3 quotes the beta-range, not this point value

**`h95d_ada_enu_germline_rank_sim`** = 839

- Value: `839`
- Formula: `rank of log P(v) = log S(v) + log σ_adagrasib(v), descending, over 1,205 class-forced variants`
- Computed from: H95D simulated rank under adagrasib + ENU with C>G = 0.2% (Arnold 2012 germline-ENU compilation)
- Inputs: `{'variant': 'H95D', 'drug': 'adagrasib', 'mutagen': 'ENU', 'CG_rate': 0.002, 'spec': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_enu_germline_shift`** = +830

- Value: `830`
- Formula: `rank_sim(ENU, p_CG) - rank_chemfree`
- Computed from: H95D rank shift under adagrasib + ENU (Arnold 2012 germline-ENU compilation)
- Inputs: `{'rank_sim': 839, 'rank_cf': 9, 'p_CG': 0.002}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`


### supporting: §2.3 quotes the interval, not its endpoints

**`h95d_ada_germline_shift_beta_min`** = +15

- Value: `15`
- Formula: `min_β rank_sim(β) − rank_chemfree, with log P = log S + β·σ`
- Computed from: minimum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU germline (C>G = 0.2%)
- Inputs: `{'beta_grid': [0.25, 0.5, 1.0, 2.0, 4.0], 'spec': '{...}', 'variant': 'H95D', 'drug': 'adagrasib'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_germline_shift_beta_max`** = +848

- Value: `848`
- Formula: `max_β rank_sim(β) − rank_chemfree, with log P = log S + β·σ`
- Computed from: maximum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU germline (C>G = 0.2%)
- Inputs: `{'beta_grid': [0.25, 0.5, 1.0, 2.0, 4.0], 'spec': '{...}', 'variant': 'H95D', 'drug': 'adagrasib'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_empirical_lower_shift_beta_min`** = +2

- Value: `2`
- Formula: `min_β rank_sim(β) − rank_chemfree, with log P = log S + β·σ`
- Computed from: minimum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU empirical lower (C>G = 7%)
- Inputs: `{'beta_grid': [0.25, 0.5, 1.0, 2.0, 4.0], 'spec': '{...}', 'variant': 'H95D', 'drug': 'adagrasib'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_empirical_lower_shift_beta_max`** = +518

- Value: `518`
- Formula: `max_β rank_sim(β) − rank_chemfree, with log P = log S + β·σ`
- Computed from: maximum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU empirical lower (C>G = 7%)
- Inputs: `{'beta_grid': [0.25, 0.5, 1.0, 2.0, 4.0], 'spec': '{...}', 'variant': 'H95D', 'drug': 'adagrasib'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_empirical_mean_shift_beta_min`** = -4

- Value: `-4`
- Formula: `min_β rank_sim(β) − rank_chemfree, with log P = log S + β·σ`
- Computed from: minimum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU empirical mean (C>G = 25%)
- Inputs: `{'beta_grid': [0.25, 0.5, 1.0, 2.0, 4.0], 'spec': '{...}', 'variant': 'H95D', 'drug': 'adagrasib'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_empirical_mean_shift_beta_max`** = -3

- Value: `-3`
- Formula: `max_β rank_sim(β) − rank_chemfree, with log P = log S + β·σ`
- Computed from: maximum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU empirical mean (C>G = 25%)
- Inputs: `{'beta_grid': [0.25, 0.5, 1.0, 2.0, 4.0], 'spec': '{...}', 'variant': 'H95D', 'drug': 'adagrasib'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_mutazyme_shift_beta_min`** = +1

- Value: `1`
- Formula: `min_β rank_sim(β) − rank_chemfree, with log P = log S + β·σ`
- Computed from: minimum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for Mutazyme II (NGS; C>G = 5%)
- Inputs: `{'beta_grid': [0.25, 0.5, 1.0, 2.0, 4.0], 'spec': '{...}', 'variant': 'H95D', 'drug': 'adagrasib'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_mutazyme_shift_beta_max`** = +606

- Value: `606`
- Formula: `max_β rank_sim(β) − rank_chemfree, with log P = log S + β·σ`
- Computed from: maximum H95D shift across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for Mutazyme II (NGS; C>G = 5%)
- Inputs: `{'beta_grid': [0.25, 0.5, 1.0, 2.0, 4.0], 'spec': '{...}', 'variant': 'H95D', 'drug': 'adagrasib'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`


### supporting: §2.4 quotes the range across drugs and mutagens

**`alk_alec_mut_top100`** = 65% [55\%--74\%]

- Value: `0.65`
- Formula: `recovery(K=100, drug, MutII)`
- Computed from: Mutazyme-simulated top-100 recovery for ALK, alectinib
- Inputs: `{'source': 'data/derived/topk_recovery.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`

**`alk_alec_enu_top100`** = 63% [53\%--72\%]

- Value: `0.63`
- Formula: `recovery(K=100, drug, ENU)`
- Computed from: ENU-simulated top-100 recovery for ALK, alectinib
- Inputs: `{'source': 'data/derived/topk_recovery.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`

**`alk_lor_mut_top100`** = 72% [63\%--80\%]

- Value: `0.72`
- Formula: `recovery(K=100, drug, MutII)`
- Computed from: Mutazyme-simulated top-100 recovery for ALK, lorlatinib
- Inputs: `{'source': 'data/derived/topk_recovery.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`

**`alk_lor_enu_top100`** = 65% [55\%--74\%]

- Value: `0.65`
- Formula: `recovery(K=100, drug, ENU)`
- Computed from: ENU-simulated top-100 recovery for ALK, lorlatinib
- Inputs: `{'source': 'data/derived/topk_recovery.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`

**`alk_zot_mut_top100`** = 65% [55\%--74\%]

- Value: `0.65`
- Formula: `recovery(K=100, drug, MutII)`
- Computed from: Mutazyme-simulated top-100 recovery for ALK, zotizalkib
- Inputs: `{'source': 'data/derived/topk_recovery.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`

**`alk_zot_enu_top100`** = 63% [53\%--72\%]

- Value: `0.63`
- Formula: `recovery(K=100, drug, ENU)`
- Computed from: ENU-simulated top-100 recovery for ALK, zotizalkib
- Inputs: `{'source': 'data/derived/topk_recovery.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`


### supporting: §2.5 quotes the beta = 4 endpoint only

**`nuvalent_cg_binom_p_model_beta_range`** = 4.1e-10 to 4.7e-05

- Value: `[4.1106573099923305e-10, 4.715801648246043e-05]`
- Formula: `min/max over beta of binom.sf(k-1, n, pi_N3(beta))`
- Computed from: range of the model-implied binomial tail across beta in [0.25, 4]
- Inputs: `{'betas': [0.25, 0.5, 1.0, 2.0, 4.0]}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`


### supporting: §2.5 states the bound, perm_max_abs_gain, rather than the signed range

**`perm_observed_gain_range`** = -33 to 0 pp

- Value: `[-0.3333333333333333, 0.0]`
- Formula: `recall@100(logS + sigma) - recall@100(sigma)`
- Computed from: observed framework-minus-sigma top-100 recall gains across anchors
- Inputs: `{'K': 100}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::run_anchor`


### §2.1

**`enu_class_C_G`** = 0.20%

- Value: `0.002`
- Formula: `spectra.yaml:ENU_Arnold2012_nominal.values['C>G']`
- Computed from: Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- Inputs: `{'source': 'arnold2012', 'class': 'C>G'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`enu_class_T_A`** = 31.4%

- Value: `0.314`
- Formula: `spectra.yaml:ENU_Arnold2012_nominal.values['T>A']`
- Computed from: Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- Inputs: `{'source': 'arnold2012', 'class': 'T>A'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`enu_class_T_C`** = 36.6%

- Value: `0.366`
- Formula: `spectra.yaml:ENU_Arnold2012_nominal.values['T>C']`
- Computed from: Arnold et al. 2012 germline-ENU compilation of DNA-adduct conversion frequencies
- Inputs: `{'source': 'arnold2012', 'class': 'T>C'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`mut_class_C_G`** = 5.0%

- Value: `0.05`
- Formula: `spectra.yaml:MutazymeII_DoiSupp2022.values['C>G']`
- Computed from: NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- Inputs: `{'source': 'doi2022', 'class': 'C>G'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`mut_class_C_T`** = 47.0%

- Value: `0.47`
- Formula: `spectra.yaml:MutazymeII_DoiSupp2022.values['C>T']`
- Computed from: NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- Inputs: `{'source': 'doi2022', 'class': 'C>T'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`mut_class_T_G`** = 3.6%

- Value: `0.036`
- Formula: `spectra.yaml:MutazymeII_DoiSupp2022.values['T>G']`
- Computed from: NGS of untreated Mutazyme II library in Doi 2022 Suppl Fig 2b
- Inputs: `{'source': 'doi2022', 'class': 'T>G'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`enu_classes_at_1pct`** = five of six

- Value: `5`
- Formula: `|{c: ENU[c] ≥ 0.01}|`
- Computed from: count of ENU classes supplied at ≥ 1 %; C>G is at 0.2 % and falls below
- Inputs: `{'ENU': '{...}', 'threshold': 0.01}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`be_classes_supplied`** = two of six

- Value: `2`
- Formula: `deterministic by base-editor chemistry`
- Computed from: CBE supplies C>T, ABE supplies T>C by design; the other 4 classes are 0
- Inputs: `{'supplied': ['C>T', 'T>C']}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`mut_class_rate_ratio`** = 13-fold

- Value: `13.055555555555555`
- Formula: `max(MutII) / min(MutII) = max{MutII[c]} / min{MutII[c]}`
- Computed from: ratio of Mutazyme II's most abundant class to its rarest
- Inputs: `{'max': 0.47, 'max_class': 'C>T', 'min': 0.036, 'min_class': 'T>G'}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`cg_detect_prob_germline`** = 0.86

- Value: `0.8646647167633873`
- Formula: `1 - exp(-N*mu*S*r)`
- Computed from: detection probability of a C>G variant at the germline ENU rate in a 1e7 Ba/F3 library; shows a rare class is still supplied
- Inputs: `{'N': 10000000.0, 'mu': 0.0001, 'S': 0.002, 'routes': 1}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`

**`cg_expected_carriers_germline`** = 2.0

- Value: `2.0`
- Formula: `N*mu*S*r`
- Computed from: expected carriers of a given C>G variant at the germline ENU rate
- Inputs: `{'N': 10000000.0, 'mu': 0.0001, 'S': 0.002, 'routes': 1}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`


### §2.2

**`kras_panel_total`** = 12

- Value: `12`
- Formula: `|clinical ∪ case ∪ DMS|`
- Computed from: union of KRAS G12C clinical cohort + case report + accompanying in vitro DMS
- Inputs: `{'clinical': ['H95D', 'H95Q', 'H95R', 'R68S', 'Y96C'], 'case': ['Y96D'], 'dms': ['C12F', 'R68L', 'R68M', 'V8L', 'V8M', 'Y96H']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`kras_panel_dms`** = six

- Value: `6`
- Formula: `|DMS|`
- Computed from: variants from Awad 2021 accompanying in vitro deep mutational scan
- Inputs: `{'members': ['C12F', 'R68L', 'R68M', 'V8L', 'V8M', 'Y96H']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`kras_panel_snv_unreachable`** = two

- Value: `2`
- Formula: `count of panel variants whose Hamming distance from WT codon is ≥ 2 to all target codons`
- Computed from: codon reachability analysis in the G12C-mutant background
- Inputs: `{'members': ['R68L', 'V8M']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`kras_panel_multi_class`** = three

- Value: `3`
- Formula: `count of panel variants with ≥ 2 distinct single-nt pyrimidine routes`
- Computed from: codon reachability: variants reachable via multiple pyrimidine classes
- Inputs: `{'members': ['H95Q', 'R68S', 'V8L']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`kras_panel_class_forced`** = seven

- Value: `7`
- Formula: `count of panel variants with exactly 1 single-nt pyrimidine route`
- Computed from: codon reachability: variants reachable via exactly one pyrimidine class
- Inputs: `{'members': ['C12F', 'H95D', 'H95R', 'R68M', 'Y96C', 'Y96D', 'Y96H']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`two_hit_carriers_one_round`** = 0.005--0.017

- Value: `[0.005076, 0.01692]`
- Formula: `sum over routes of N*mu^2*S(c1)*S(c2)*2`
- Computed from: expected carriers of R68L or V8M after a single mutagenesis round in a 1e7 library, across ENU and Mutazyme II
- Inputs: `{'N': 10000000.0, 'mu': 0.0001, 'rounds': 1}`
- Code: `scripts/compute_paper_numbers.py::two_hit_reachability.py::main`

**`two_hit_carriers_twenty_rounds`** = 1.1--3.6

- Value: `[1.06596, 3.5532000000000004]`
- Formula: `lambda_simultaneous * R(R+1)/2`
- Computed from: expected carriers of R68L or V8M after twenty rounds
- Inputs: `{'N': 10000000.0, 'mu': 0.0001, 'rounds': 20}`
- Code: `scripts/compute_paper_numbers.py::two_hit_reachability.py::main`

**`two_hit_fold_suppression`** = 51- to 169-fold

- Value: `[50.65856129685916, 168.86187098953056]`
- Formula: `lambda_1nt(median class) / lambda_2nt`
- Computed from: suppression of R68L and V8M relative to a one-step variant of median class supply in the same library, ENU and Mutazyme II
- Inputs: `{'N': 10000000.0, 'mu': 0.0001, 'rounds': 20}`
- Code: `scripts/compute_paper_numbers.py::two_hit_reachability.py::main`

**`r68l_two_nt_routes`** = 2

- Value: `2`
- Formula: `enumerate codons at Hamming distance 2 translating to Leu`
- Computed from: two-substitution routes from AGG to a leucine codon
- Inputs: `{'wt_codon': 'AGG', 'target_aa': 'L'}`
- Code: `scripts/compute_paper_numbers.py::two_hit_reachability.py::two_nt_routes`

**`v8m_two_nt_routes`** = 1

- Value: `1`
- Formula: `enumerate codons at Hamming distance 2 translating to Met`
- Computed from: two-substitution routes from GTA to the methionine codon
- Inputs: `{'wt_codon': 'GTA', 'target_aa': 'M'}`
- Code: `scripts/compute_paper_numbers.py::two_hit_reachability.py::two_nt_routes`


### §2.3

**`kras_panel_class_forced`** = seven

- Value: `7`
- Formula: `count of panel variants with exactly 1 single-nt pyrimidine route`
- Computed from: codon reachability: variants reachable via exactly one pyrimidine class
- Inputs: `{'members': ['C12F', 'H95D', 'H95R', 'R68M', 'Y96C', 'Y96D', 'Y96H']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`feng_class_forced_count`** = 981

- Value: `981`
- Formula: `|{v in Feng: v is class-forced single-nt reachable}|`
- Computed from: class-forced KRAS G12C variants in the Feng atlas (unique amino-acid changes)
- Inputs: `{'source': 'data/derived/feng_decomposition.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`h95d_ada_chemfree_rank`** = 9

- Value: `9`
- Formula: `rank(σ_adagrasib) descending over the class-forced KRAS G12C atlas, position of H95D`
- Computed from: H95D chemistry-free rank in the KRAS G12C class-forced atlas under adagrasib
- Inputs: `{'variant': 'H95D', 'drug': 'adagrasib'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_mut_rank_sim`** = 56

- Value: `56`
- Formula: `rank of log P(v) = log S(v) + log σ_adagrasib(v) under Mutazyme II spectrum`
- Computed from: H95D simulated rank under adagrasib + Mutazyme II (NGS-measured spectrum)
- Inputs: `{'variant': 'H95D', 'drug': 'adagrasib', 'mutagen': 'MutazymeII'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_mut_shift`** = +47

- Value: `47`
- Formula: `rank_sim(Mutazyme) - rank_chemfree`
- Computed from: H95D rank shift under adagrasib + Mutazyme II (NGS-measured)
- Inputs: `{'rank_sim': 56, 'rank_cf': 9, 'mut_CG': 0.05}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`h95d_ada_germline_shift_beta_range`** = +15 to +848

- Value: `+15 to +848`
- Formula: `[min_β shift, max_β shift]`
- Computed from: H95D shift interval across β ∈ [0.25, 0.5, 1.0, 2.0, 4.0] for ENU germline (C>G = 0.2%)
- Inputs: `{'beta_grid': [0.25, 0.5, 1.0, 2.0, 4.0], 'per_beta': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`y96d_sot_chemfree_rank`** = 7

- Value: `7`
- Formula: `rank(σ_sotorasib) descending, position of Y96D`
- Computed from: Y96D chemistry-free rank under sotorasib
- Inputs: `{'variant': 'Y96D', 'drug': 'sotorasib'}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`y96d_sot_mut_rank_sim`** = 139

- Value: `139`
- Formula: `rank of log P(v) under Mutazyme spectrum, Y96D forced class T>G at 3.6 %`
- Computed from: Y96D simulated rank under sotorasib + Mutazyme II
- Inputs: `{'variant': 'Y96D', 'drug': 'sotorasib', 'mutagen': 'MutazymeII', 'mut_TG': 0.036}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`y96d_sot_mut_shift`** = +132

- Value: `132`
- Formula: `rank_sim - rank_chemfree for Y96D under Mutazyme sotorasib`
- Computed from: Y96D rank shift under sotorasib + Mutazyme II
- Inputs: `{'rank_sim': 139, 'rank_cf': 7}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`


### §2.4

**`alk_panel_total`** = 12

- Value: `12`
- Formula: `|Doi hits ∪ Takei hits|`
- Computed from: union of alectinib Mutazyme II screen hits + next-gen ALK TKI screen hits
- Inputs: `{'doi': ['G1202R', 'I1171N', 'I1171T', 'L1196M', 'L1198F', 'V1180L'], 'takei': ['C1156Y', 'F1174C', 'F1174L', 'G1269A', 'L1152R', 'L1196Q']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`doi2022_hits_count`** = six

- Value: `6`
- Formula: `|Doi Fig 3a resistance mutations|`
- Computed from: single-amino-acid resistance mutations in Doi 2022 Fig 3a, excluding silent passengers (L1152, L1165, R1192) and the I1171S sub-variant of the I1171N/T/S compound call
- Inputs: `{'members': ['G1202R', 'I1171N', 'I1171T', 'L1196M', 'L1198F', 'V1180L']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`oh_class_forced_count`** = 1,427

- Value: `1427`
- Formula: `|unique aa_change in Oh class-forced decomposition|`
- Computed from: class-forced ALK kinase-domain variants used in the ranking
- Inputs: `{'source': 'data/derived/oh_decomposition_all_drugs.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`doi_anchor_topN_bound`** = top-120

- Value: `120`
- Formula: `ceil(max rank / 10) * 10`
- Computed from: smallest round bound containing all five class-forced Doi hits in the Mutazyme-simulated ranking
- Inputs: `{'max_rank': 118}`
- Code: `scripts/compute_paper_numbers.py::compute_doi_anchor_ranks`

**`alk_snv_top100_recovery_min`** = 63%

- Value: `0.63`
- Formula: `min_{drug, mutagen in {ENU, MutII}} recovery(K=100, drug, mutagen)`
- Computed from: min top-100 recovery across ALK drugs for ENU and Mutazyme II
- Inputs: `{'source': 'data/derived/topk_recovery.csv', 'gene': 'ALK', 'K': 100, 'per_row': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`

**`alk_snv_top100_recovery_max`** = 72%

- Value: `0.72`
- Formula: `max_{drug, mutagen in {ENU, MutII}} recovery(K=100, drug, mutagen)`
- Computed from: max top-100 recovery across ALK drugs for ENU and Mutazyme II
- Inputs: `{'source': 'data/derived/topk_recovery.csv', 'gene': 'ALK', 'K': 100, 'per_row': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`

**`kras_snv_top100_recovery_min`** = 40%

- Value: `0.4`
- Formula: `min_{drug, mutagen in {ENU, MutII}} recovery(K=100, drug, mutagen)`
- Computed from: min top-100 recovery across KRAS G12C drugs for ENU and Mutazyme II
- Inputs: `{'source': 'data/derived/topk_recovery.csv', 'gene': 'KRAS G12C', 'K': 100, 'per_row': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`

**`kras_snv_top100_recovery_max`** = 54%

- Value: `0.54`
- Formula: `max_{drug, mutagen in {ENU, MutII}} recovery(K=100, drug, mutagen)`
- Computed from: max top-100 recovery across KRAS G12C drugs for ENU and Mutazyme II
- Inputs: `{'source': 'data/derived/topk_recovery.csv', 'gene': 'KRAS G12C', 'K': 100, 'per_row': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`

**`be_alk_top10_min`** = 0.00%

- Value: `0.0`
- Formula: `min_{drug in {alec, lor, zot}} recovery(K=10, drug, BE)`
- Computed from: min CBE+ABE top-10 recovery across ALK drugs (alectinib shows 0)
- Inputs: `{'source': 'data/derived/topk_recovery.csv', 'per_drug': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`

**`be_alk_top10_max`** = 10%

- Value: `0.1`
- Formula: `max_{drug in {alec, lor, zot}} recovery(K=10, drug, BE)`
- Computed from: max CBE+ABE top-10 recovery across ALK drugs
- Inputs: `{'source': 'data/derived/topk_recovery.csv', 'per_drug': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`

**`be_alk_top10_outside_min`** = 90%

- Value: `0.9`
- Formula: `1 - max_{drug} recovery(K=10, drug, BE)`
- Computed from: 1 - max(BE top-10 recovery across ALK drugs)
- Inputs: `{'source': 'data/derived/topk_recovery.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`

**`be_alk_top10_outside_max`** = 100%

- Value: `1.0`
- Formula: `1 - min_{drug} recovery(K=10, drug, BE)`
- Computed from: 1 - min(BE top-10 recovery across ALK drugs)
- Inputs: `{'source': 'data/derived/topk_recovery.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`

**`be_alk_peak_at_k_le_100`** = 38%

- Value: `0.3829787234042553`
- Formula: `max_{K ≤ 100, drug} recovery(K, drug, BE) for gene=ALK`
- Computed from: max CBE+ABE recovery over K ≤ 100 across all ALK drugs
- Inputs: `{'source': 'data/derived/topk_recovery.csv', 'gene': 'ALK', 'K_max': 100}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`

**`be_kras_peak_at_k_le_100`** = 100%

- Value: `1.0`
- Formula: `max_{K ≤ 100, drug} recovery(K, drug, BE) for gene=KRAS G12C`
- Computed from: max CBE+ABE recovery over K ≤ 100 across all KRAS G12C drugs
- Inputs: `{'source': 'data/derived/topk_recovery.csv', 'gene': 'KRAS G12C', 'K_max': 100}`
- Code: `scripts/compute_paper_numbers.py::compute_topk_recovery`


### §2.5

**`kras_panel_snv_unreachable`** = two

- Value: `2`
- Formula: `count of panel variants whose Hamming distance from WT codon is ≥ 2 to all target codons`
- Computed from: codon reachability analysis in the G12C-mutant background
- Inputs: `{'members': ['R68L', 'V8M']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`doi2022_hits_count`** = six

- Value: `6`
- Formula: `|Doi Fig 3a resistance mutations|`
- Computed from: single-amino-acid resistance mutations in Doi 2022 Fig 3a, excluding silent passengers (L1152, L1165, R1192) and the I1171S sub-variant of the I1171N/T/S compound call
- Inputs: `{'members': ['G1202R', 'I1171N', 'I1171T', 'L1196M', 'L1198F', 'V1180L']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`doi2022_hits_class_forced`** = five

- Value: `5`
- Formula: `|Doi hits ∩ class-forced|`
- Computed from: Doi 2022 hits that are class-forced (G1202R is multi-class)
- Inputs: `{'members': ['I1171N', 'I1171T', 'L1196M', 'L1198F', 'V1180L']}`
- Code: `scripts/compute_paper_numbers.py::compute_resistance_panels`

**`enu_well_supplied_well_vs_total`** = nine of twelve

- Value: `9 of 12`
- Formula: `|{hit: ENU[hit.class] ≥ 0.01, hit in Nuvalent class-forced pool}|`
- Computed from: count of Nuvalent class-forced hits whose forced class ENU supplies at ≥ 1 %
- Inputs: `{'nuvalent_class_forced_total': 12, 'threshold': 0.01}`
- Code: `scripts/compute_paper_numbers.py::compute_mutagen_spectra`

**`h95d_ada_mut_shift`** = +47

- Value: `47`
- Formula: `rank_sim(Mutazyme) - rank_chemfree`
- Computed from: H95D rank shift under adagrasib + Mutazyme II (NGS-measured)
- Inputs: `{'rank_sim': 56, 'rank_cf': 9, 'mut_CG': 0.05}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`

**`nuvalent_cg_binom_p_raw`** = 2.3e-06

- Value: `2.253926009047811e-06`
- Formula: `binom.sf(k-1, n, S_ENU('C>G'))`
- Computed from: exact binomial tail against the raw germline spectrum; retained only as the confounded comparison
- Inputs: `{'k': 3, 'n': 13, 'p': 0.002}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`

**`nuvalent_cg_pi_model`** = 0.56%

- Value: `0.00556039021953432`
- Formula: `sum_{v in C>G} S(c)exp(beta*sigma_v) / sum_v S(c)exp(beta*sigma_v)`
- Computed from: model-implied probability that a hit is C>G, with germline supply and pooled sigma at beta=1
- Inputs: `{'beta': 1.0, 'sigma': 'max over alec/lor/zot', 'S': 'ENU germline'}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::pi_model`

**`nuvalent_cg_binom_p_model`** = 4.7e-05

- Value: `4.715801648246043e-05`
- Formula: `binom.sf(k-1, n, pi_N3)`
- Computed from: binomial tail of 3-of-13 C>G against the model-implied null; the germline spectrum is rejected
- Inputs: `{'k': 3, 'n': 13, 'pi': 0.00556039021953432}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`

**`nuvalent_cg_binom_p_selection`** = 0.64

- Value: `0.636406728920353`
- Formula: `binom.sf(k-1, n, pi_N2)`
- Computed from: binomial tail of 3-of-13 C>G against a selection-only null; the count is unremarkable once supply is removed
- Inputs: `{'k': 3, 'n': 13, 'pi': 0.24}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`

**`nuvalent_cg_binom_p_composition`** = 0.39

- Value: `0.38585678830594816`
- Formula: `binom.sf(k-1, n, n_CG / n_classforced)`
- Computed from: binomial tail against atlas class composition alone
- Inputs: `{'n_CG': 243, 'n': 1427}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`

**`atlas_cg_share`** = 17.0%

- Value: `0.1702873160476524`
- Formula: `n_CG / n_classforced`
- Computed from: share of class-forced ALK atlas variants whose forced class is C>G
- Inputs: `{'n_CG': 243, 'n': 1427}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`

**`implied_enu_cg_supply_interval`** = 2.6%--26.2%

- Value: `[0.026107666101764928, 0.2616117716315627]`
- Formula: `solve pi_N3(C>G; s) = p for s, at p = Jeffreys lo/hi`
- Computed from: Jeffreys interval on the hit proportion mapped back through the model-implied null to an implied C>G supply rate; reported only to show it is narrower than the uncorrected inversion, and NOT adopted, because the inversion presumes supply controls hit composition
- Inputs: `{'jeffreys': '{...}', 'beta': 1.0}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::solve_s`

**`cg_detect_prob_germline`** = 0.86

- Value: `0.8646647167633873`
- Formula: `1 - exp(-N*mu*S*r)`
- Computed from: detection probability of a C>G variant at the germline ENU rate in a 1e7 Ba/F3 library; shows a rare class is still supplied
- Inputs: `{'N': 10000000.0, 'mu': 0.0001, 'S': 0.002, 'routes': 1}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`

**`cg_expected_carriers_germline`** = 2.0

- Value: `2.0`
- Formula: `N*mu*S*r`
- Computed from: expected carriers of a given C>G variant at the germline ENU rate
- Inputs: `{'N': 10000000.0, 'mu': 0.0001, 'S': 0.002, 'routes': 1}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::main`

**`perm_mdg_range`** = 9--50 pp

- Value: `[0.0909090909090909, 0.5]`
- Formula: `quantile(|null gain|, 0.975) per anchor`
- Computed from: 97.5th percentile of the permutation null's absolute top-100 recall gain, across all six anchor-proxy combinations
- Inputs: `{'n_iter': 5000, 'anchors': 6}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`

**`perm_gain_quantum_range`** = 9--33 pp

- Value: `[0.0909090909090909, 0.3333333333333333]`
- Formula: `1 / n_hits`
- Computed from: smallest observable recall change, 1/n_hits, across anchors
- Inputs: `{'n_hits': [3, 4, 5, 11]}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`

**`perm_max_abs_gain`** = 33 pp

- Value: `0.3333333333333333`
- Formula: `max |recall@100(logS + sigma) - recall@100(sigma)|`
- Computed from: largest absolute framework-minus-sigma top-100 recall gain over all six anchor-proxy combinations
- Inputs: `{'K': 100, 'anchors': 6}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::run_anchor`

**`perm_p_doi_alec`** = 1.00

- Value: `1.0`
- Formula: `(#{|null gain| >= |observed gain|} + 1) / (n_iter + 1)`
- Computed from: class-label permutation on the Doi 2022 (ALK Mutazyme × alec) anchor
- Inputs: `{'n_iter': 5000, 'K': 100, 'n_hits': 5}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`

**`perm_p_nuvalent_alec`** = 1.00

- Value: `1.0`
- Formula: `(#{|null gain| >= |observed gain|} + 1) / (n_iter + 1)`
- Computed from: class-label permutation on the Nuvalent (ALK ENU × alec proxy) anchor
- Inputs: `{'n_iter': 5000, 'K': 100, 'n_hits': 11}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`

**`perm_p_nuvalent_lor`** = 1.00

- Value: `1.0`
- Formula: `(#{|null gain| >= |observed gain|} + 1) / (n_iter + 1)`
- Computed from: class-label permutation on the Nuvalent (ALK ENU × lor proxy) anchor
- Inputs: `{'n_iter': 5000, 'K': 100, 'n_hits': 11}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`

**`perm_p_nuvalent_zot`** = 1.00

- Value: `1.0`
- Formula: `(#{|null gain| >= |observed gain|} + 1) / (n_iter + 1)`
- Computed from: class-label permutation on the Nuvalent (ALK ENU × zot proxy) anchor
- Inputs: `{'n_iter': 5000, 'K': 100, 'n_hits': 11}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`

**`perm_p_awad_sot`** = 0.30

- Value: `0.2999400119976004`
- Formula: `(#{|null gain| >= |observed gain|} + 1) / (n_iter + 1)`
- Computed from: class-label permutation on the Awad clinical (KRAS Mut × sot) anchor
- Inputs: `{'n_iter': 5000, 'K': 100, 'n_hits': 3}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`

**`perm_p_awad_ada`** = 0.60

- Value: `0.6012797440511898`
- Formula: `(#{|null gain| >= |observed gain|} + 1) / (n_iter + 1)`
- Computed from: class-label permutation on the Awad clinical (KRAS Mut × ada) anchor
- Inputs: `{'n_iter': 5000, 'K': 100, 'n_hits': 4}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::permutation_p`

**`perm_n_hits_range`** = 3 to 11

- Value: `[3, 11]`
- Formula: `count of named anchor hits present in the class-forced universe`
- Computed from: class-forced observed hits per anchor
- Inputs: `{'anchors': 6}`
- Code: `scripts/compute_paper_numbers.py::recompute_benchmark_r3.py::run_anchor`

**`doi_l1196m_rank_sim`** = 1

- Value: `1`
- Formula: `rank of log S + sigma, descending, over the class-forced ALK universe`
- Computed from: L1196M rank under the Mutazyme-simulated alectinib ranking
- Inputs: `{'variant': 'L1196M', 'mutagen': 'MutazymeII', 'drug': 'alectinib'}`
- Code: `scripts/compute_paper_numbers.py::notebook/05_prospective_predictions.py`

**`doi_l1198f_rank_sim`** = 3

- Value: `3`
- Formula: `rank of log S + sigma, descending, over the class-forced ALK universe`
- Computed from: L1198F rank under the Mutazyme-simulated alectinib ranking
- Inputs: `{'variant': 'L1198F', 'mutagen': 'MutazymeII', 'drug': 'alectinib'}`
- Code: `scripts/compute_paper_numbers.py::notebook/05_prospective_predictions.py`

**`doi_i1171t_rank_sim`** = 10

- Value: `10`
- Formula: `rank of log S + sigma, descending, over the class-forced ALK universe`
- Computed from: I1171T rank under the Mutazyme-simulated alectinib ranking
- Inputs: `{'variant': 'I1171T', 'mutagen': 'MutazymeII', 'drug': 'alectinib'}`
- Code: `scripts/compute_paper_numbers.py::notebook/05_prospective_predictions.py`

**`doi_v1180l_rank_sim`** = 15

- Value: `15`
- Formula: `rank of log S + sigma, descending, over the class-forced ALK universe`
- Computed from: V1180L rank under the Mutazyme-simulated alectinib ranking
- Inputs: `{'variant': 'V1180L', 'mutagen': 'MutazymeII', 'drug': 'alectinib'}`
- Code: `scripts/compute_paper_numbers.py::notebook/05_prospective_predictions.py`

**`doi_i1171n_rank_sim`** = 118

- Value: `118`
- Formula: `rank of log S + sigma, descending, over the class-forced ALK universe`
- Computed from: I1171N rank under the Mutazyme-simulated alectinib ranking
- Inputs: `{'variant': 'I1171N', 'mutagen': 'MutazymeII', 'drug': 'alectinib'}`
- Code: `scripts/compute_paper_numbers.py::notebook/05_prospective_predictions.py`

**`doi_anchor_topN_bound`** = top-120

- Value: `120`
- Formula: `ceil(max rank / 10) * 10`
- Computed from: smallest round bound containing all five class-forced Doi hits in the Mutazyme-simulated ranking
- Inputs: `{'max_rank': 118}`
- Code: `scripts/compute_paper_numbers.py::compute_doi_anchor_ranks`

**`takei_v2_count`** = five

- Value: `5`
- Formula: `|Takei V2 hits at 500 nM neladalkib|`
- Computed from: V2 compound partners recovered at 500 nM neladalkib on I1171N fixed-V1 background
- Inputs: `{'members': ['L1256F', 'E1129V', 'E1129K', 'C1156Y', 'T1151R']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`

**`takei_v2_c_t_count`** = three

- Value: `3`
- Formula: `|{v in Takei V2: forced_class(v) == 'C>T'}|`
- Computed from: Takei V2 partners class-forced into C>T (Mutazyme II's dominant class)
- Inputs: `{'members': ['L1256F', 'E1129K', 'C1156Y']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`

**`nuvalent_total`** = thirteen

- Value: `13`
- Formula: `|Nuvalent single-amino-acid hits (aggregate)|`
- Computed from: single-amino-acid variants reported in Nuvalent AACR 2025 poster aggregate across arms
- Inputs: `{'members': ['I1171N', 'I1171S', 'I1171T', 'F1174C', 'F1174V', 'V1180L', 'T1151K', 'T1151R', 'E1129K', 'R1192P', 'L1196M', 'G1269A', 'G1202R']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`

**`nuvalent_class_forced`** = twelve

- Value: `11`
- Formula: `|Nuvalent hits ∩ class-forced|`
- Computed from: Nuvalent hits that are class-forced single-nt reachable (R1192P outside atlas, G1202R multi-class)
- Inputs: `{'members': ['I1171N', 'I1171S', 'I1171T', 'F1174C', 'F1174V', 'V1180L', 'T1151K', 'T1151R', 'E1129K', 'L1196M', 'G1269A']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`

**`nuvalent_cg_count`** = three

- Value: `3`
- Formula: `|{v in Nuvalent: forced_class(v) == 'C>G'}|`
- Computed from: Nuvalent hits class-forced into C>G (ENU 0.2 % class per germline spectrum)
- Inputs: `{'members': ['V1180L', 'T1151R', 'G1269A']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`

**`nuvalent_well_supplied_count`** = nine

- Value: `8`
- Formula: `|Nuvalent class-forced| - |Nuvalent C>G|`
- Computed from: Nuvalent class-forced hits whose forced class ENU supplies at ≥ 1 %
- Inputs: `{'members': ['I1171N', 'I1171S', 'I1171T', 'F1174C', 'F1174V', 'T1151K', 'E1129K', 'L1196M']}`
- Code: `scripts/compute_paper_numbers.py::compute_retrospective_screens`


### §2.5 (as 24 of the top-100 sigma-ranked class-forced variants)

**`nuvalent_cg_top100_sigma`** = 24 of the

- Value: `24`
- Formula: `count(forced_class == 'C>G') over nlargest(100, pooled sigma)`
- Computed from: C>G count among the top-100 sigma-ranked class-forced ALK variants, pooled across drug arms
- Inputs: `{'N_top': 100}`
- Code: `scripts/compute_paper_numbers.py::recompute_enu_class_test.py::pi_sigma_topn`


### §2.6

**`transport_rho_alectinib`** = 0.51

- Value: `0.5149638802889577`
- Formula: `spearmanr(sigma_atlas, median log10(IC50_variant / IC50_WT) within study)`
- Computed from: Spearman rank correlation between atlas sigma and measured Ba/F3 log10 fold resistance for alectinib, over 18 shared variants
- Inputs: `{'n': 18, 'boot95': [-0.04275959677304641, 0.855497382198953], 'p': 0.028748754882630345, 'median_inter_study_spread_log10': 0.12928740120303056}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`transport_rho_lorlatinib`** = 0.17

- Value: `0.16615067079463364`
- Formula: `spearmanr(sigma_atlas, median log10(IC50_variant / IC50_WT) within study)`
- Computed from: Spearman rank correlation between atlas sigma and measured Ba/F3 log10 fold resistance for lorlatinib, over 18 shared variants
- Inputs: `{'n': 18, 'boot95': [-0.42618532491530653, 0.7229182515739769], 'p': 0.5099504603438955, 'median_inter_study_spread_log10': 0.3219465077435653}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`transport_rho_zotizalkib`** = 0.82

- Value: `0.8241758241758241`
- Formula: `spearmanr(sigma_atlas, median log10(IC50_variant / IC50_WT) within study)`
- Computed from: Spearman rank correlation between atlas sigma and measured Ba/F3 log10 fold resistance for zotizalkib, over 13 shared variants
- Inputs: `{'n': 13, 'boot95': [0.4871720116618075, 0.9602272727272727], 'p': 0.0005302228375973589, 'median_inter_study_spread_log10': 0.0}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`proxy_rho_alectinib`** = 0.21

- Value: `0.2087912087912088`
- Formula: `spearmanr(sigma_atlas[drug], log10 FC under neladalkib)`
- Computed from: Spearman correlation between alectinib atlas sigma and measured Ba/F3 neladalkib fold resistance; ranks the proxy choices
- Inputs: `{'n': 13, 'boot95': [-0.43016759776536306, 0.726501937984496], 'p': 0.493621819356885}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`proxy_rho_lorlatinib`** = 0.57

- Value: `0.5714285714285714`
- Formula: `spearmanr(sigma_atlas[drug], log10 FC under neladalkib)`
- Computed from: Spearman correlation between lorlatinib atlas sigma and measured Ba/F3 neladalkib fold resistance; ranks the proxy choices
- Inputs: `{'n': 13, 'boot95': [0.0, 0.898883485978661], 'p': 0.0413422673717235}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`proxy_rho_zotizalkib`** = 0.57

- Value: `0.5659340659340659`
- Formula: `spearmanr(sigma_atlas[drug], log10 FC under neladalkib)`
- Computed from: Spearman correlation between zotizalkib atlas sigma and measured Ba/F3 neladalkib fold resistance; ranks the proxy choices
- Inputs: `{'n': 13, 'boot95': [-0.05148088531187118, 0.9150201661305036], 'p': 0.04379259223545486}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`atlas_dose_ratio_alectinib`** = 1.7x

- Value: `1.7241379310344829`
- Formula: `atlas_dose_nM / median(WT Ba/F3 IC50_nM)`
- Computed from: atlas assay concentration for alectinib divided by the measured Ba/F3 wild-type IC50
- Inputs: `{'dose_nM': 20, 'wt_ic50_nM': 11.6}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`atlas_dose_ratio_lorlatinib`** = 3.7x

- Value: `3.7037037037037033`
- Formula: `atlas_dose_nM / median(WT Ba/F3 IC50_nM)`
- Computed from: atlas assay concentration for lorlatinib divided by the measured Ba/F3 wild-type IC50
- Inputs: `{'dose_nM': 5, 'wt_ic50_nM': 1.35}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`atlas_dose_ratio_zotizalkib`** = 30.6x

- Value: `30.612244897959187`
- Formula: `atlas_dose_nM / median(WT Ba/F3 IC50_nM)`
- Computed from: atlas assay concentration for zotizalkib divided by the measured Ba/F3 wild-type IC50
- Inputs: `{'dose_nM': 300, 'wt_ic50_nM': 9.799999999999999}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::main`

**`transport_n_studies`** = 7

- Value: `7`
- Formula: `nunique(paper_id)`
- Computed from: independent published studies contributing Ba/F3 ALK potency measurements
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py::load_potency`

**`dose_rho_low_high_range`** = 0.30 to 0.63

- Value: `[0.3017885086161378, 0.6314935064935066]`
- Formula: `spearmanr(score_lowest_dose, score_highest_dose) per series`
- Computed from: Spearman between a drug's lowest and highest assayed concentration, across the seven multi-dose series in Oh Tables S3, S5 and S9
- Inputs: `{'n_series': 7, 'median': 0.38887414804516635}`
- Code: `scripts/compute_paper_numbers.py::oh_dose_transport.py::main`

**`dose_rho_low_high_median`** = 0.39

- Value: `0.38887414804516635`
- Formula: `median over series`
- Computed from: median across the seven dose series of the lowest-to-highest rank correlation
- Inputs: `{'n_series': 7}`
- Code: `scripts/compute_paper_numbers.py::oh_dose_transport.py::main`

**`dose_topk_overlap_range`** = between 5 and 15 of the top 25

- Value: `[5, 15]`
- Formula: `|top25(lowest) intersect top25(highest)|`
- Computed from: variants shared between the top 25 at a drug's lowest and highest dose
- Inputs: `{'top_k': 25}`
- Code: `scripts/compute_paper_numbers.py::oh_dose_transport.py::main`

**`dose_rho_adjacent_median`** = 0.64

- Value: `0.6433936370601897`
- Formula: `median over all adjacent dose pairs`
- Computed from: median rank correlation between adjacent concentration steps
- Code: `scripts/compute_paper_numbers.py::oh_dose_transport.py::main`

**`dose_cross_drug_matched_rho`** = 0.75

- Value: `0.7506117071334464`
- Formula: `median spearmanr over cross-drug pairs with equal dose rank`
- Computed from: median rank correlation between DIFFERENT drugs assayed at matched dose rank; exceeds the same-drug across-dose figure, indicating shared batch structure
- Inputs: `{'n_pairs': 9}`
- Code: `scripts/compute_paper_numbers.py::oh_dose_transport.py::main`

**`dose_cross_drug_unmatched_rho`** = 0.57

- Value: `0.5668906011254523`
- Formula: `median spearmanr over cross-drug pairs with unequal dose rank`
- Computed from: same, for cross-drug pairs at unmatched dose rank
- Inputs: `{'n_pairs': 18}`
- Code: `scripts/compute_paper_numbers.py::oh_dose_transport.py::main`

**`atlas_vs_own_repeat_rho`** = 0.55--0.62

- Value: `[0.5503357540460461, 0.6169574886628025]`
- Formula: `max over doses of spearmanr(S2 score, S5 score)`
- Computed from: agreement between the main single-dose atlas and its best-matching sub-library condition: two runs of nominally the same assay
- Inputs: `{'alectinib': '20.0', 'lorlatinib': '5.0', 'zotizalkib': '0.9'}`
- Code: `scripts/compute_paper_numbers.py::oh_dose_transport.py::main`

**`dose_window_sd_ratio`** = 1.1 to 1.5

- Value: `[1.5290723518547962, 1.5346922882554717, 1.0807968335316527]`
- Formula: `sd(window) / sd(full atlas), per drug`
- Computed from: sigma dispersion in the exon-23 sub-library window relative to the full atlas; above 1 rules out attenuation by range restriction
- Inputs: `{'window': [1172, 1215], 'n_window': 284, 'n_full': 2173}`
- Code: `scripts/compute_paper_numbers.py::oh_dose_transport.py::main`


### §2.7

**`met_two_screen_chi2`** = 12.4

- Value: `12.420123432255783`
- Formula: `scipy.stats.chi2_contingency(2 x 6 clone counts)`
- Computed from: chi-square homogeneity of the six-class clone composition between the NVP-BVU972 and AMG 458 screens, which used the same ENU protocol in the same cells
- Inputs: `{'dof': 5, 'p': 0.029463225369071463, 'table': '{...}'}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_two_screen_p`** = 0.029

- Value: `0.029463225369071463`
- Formula: `scipy.stats.chi2_contingency`
- Computed from: p-value of that homogeneity test; the two screens differ, so the class composition of a hit list is set by the drug, not the mutagen
- Inputs: `{'dof': 5}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_recovered`** = 5 of 6

- Value: `5`
- Formula: `enumerate one-step routes from the codon; intersect with Tiedt Table S1`
- Computed from: single-nucleotide missense routes out of the Y1230 wild-type codon that the NVP-BVU972 screen recovered
- Inputs: `{'wt_codon': 'TAT', 'n_clones': 75}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_supply_fold`** = 7.3-fold

- Value: `7.319999999999999`
- Formula: `max(S_ENU(c)) / min(S_ENU(c)) over the routes' classes`
- Computed from: ratio of the largest to smallest germline ENU supply among the substitution classes those Y1230 routes span
- Inputs: `{'classes': ['T>A', 'T>C', 'T>G'], 'lo': 0.05, 'hi': 0.366}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_tg_expected_clones`** = 5.1

- Value: `5.136986301369864`
- Formula: `sum of expected_clones_under_supply over the T>G routes`
- Computed from: clones the two T>G routes out of Y1230 are predicted to carry under germline ENU supply alone
- Inputs: `{'routes': ['Y1230D', 'Y1230S']}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_tg_observed_clones`** = 5

- Value: `5`
- Formula: `sum of clones_observed over the T>G routes`
- Computed from: clones the two T>G routes out of Y1230 actually carried in the Tiedt NVP-BVU972 screen
- Inputs: `{'routes': ['Y1230D', 'Y1230S']}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_missing`** = Y1230F

- Value: `['Y1230F']`
- Formula: `routes minus recovered`
- Computed from: one-step Y1230 variants the screen did not recover, with their rank in the chemistry-free atlas
- Inputs: `{'Y1230F': 54}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_nvp_bvu972_clones`** = 85

- Value: `85`
- Formula: `sum of clone counts in Tiedt Supplementary Table`
- Computed from: resistant clones reported for the NVP-BVU972 screen
- Inputs: `{'n_base_changes': 10, 'n_aa_changes': 10}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_nvp_bvu972_median_rank`** = 171

- Value: `171`
- Formula: `median(rank by delta sigma) over reported hits`
- Computed from: median chemistry-free atlas rank of the NVP-BVU972 hits
- Inputs: `{'atlas_n': 5458, 'found': 10, 'total': 10}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_nvp_bvu972_twostep_top100`** = 72

- Value: `72`
- Formula: `classify each top-100 variant by one-step routes from its codon`
- Computed from: variants in the chemistry-free top 100 for NVP-BVU972 that need two or more nucleotide changes and so lie outside one-step reach
- Inputs: `{'two-step': 72, 'class-forced': 26, 'multi-class SNV': 2}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_amg_clones`** = 60

- Value: `60`
- Formula: `sum of clone counts in Tiedt Supplementary Table`
- Computed from: resistant clones reported for the AMG 458 screen
- Inputs: `{'n_base_changes': 15, 'n_aa_changes': 13}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_amg_median_rank`** = 229

- Value: `229`
- Formula: `median(rank by delta sigma) over reported hits`
- Computed from: median chemistry-free atlas rank of the AMG 458 hits
- Inputs: `{'atlas_n': 5458, 'found': 13, 'total': 13}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_amg_twostep_top100`** = 75

- Value: `75`
- Formula: `classify each top-100 variant by one-step routes from its codon`
- Computed from: variants in the chemistry-free top 100 for AMG 458 that need two or more nucleotide changes and so lie outside one-step reach
- Inputs: `{'class-forced': 23, 'two-step': 75, 'multi-class SNV': 2}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_atlas_size`** = 5,458

- Value: `5458`
- Formula: `rows of the wide table for that background`
- Computed from: MET kinase-domain missense variants scored in the Estevam atlas on the primary background
- Inputs: `{'background': 'exon14del'}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`


### §3.2

**`reach_be_blocked_alk`** = 172 of 273

- Value: `172`
- Formula: `count of residues where S_CBE+ABE of the top-sigma class-forced substitution is zero`
- Computed from: ALK residues whose top-sigma class-forced substitution CBE+ABE cannot make, under alectinib
- Inputs: `{'drug': 'alectinib', 'n_residues': 273, 'note': 'the arm Figure 3C plots'}`
- Code: `scripts/compute_paper_numbers.py::build_reachability_track.py::build`

**`reach_be_blocked_kras`** = 121 of 188

- Value: `121`
- Formula: `count of residues where S_CBE+ABE of the top-sigma class-forced substitution is zero`
- Computed from: KRAS G12C residues whose top-sigma class-forced substitution CBE+ABE cannot make, under adagrasib
- Inputs: `{'drug': 'adagrasib', 'n_residues': 188, 'note': 'the arm Figure 3C plots'}`
- Code: `scripts/compute_paper_numbers.py::build_reachability_track.py::build`

**`reach_be_blocked_met`** = 183 of 287

- Value: `183`
- Formula: `count of residues where S_CBE+ABE of the top-sigma class-forced substitution is zero`
- Computed from: MET residues whose top-sigma class-forced substitution CBE+ABE cannot make, under NVP-BVU972
- Inputs: `{'drug': 'NVP-BVU972', 'n_residues': 287, 'note': 'the arm Figure 3C plots'}`
- Code: `scripts/compute_paper_numbers.py::build_reachability_track.py::build`

**`reach_be_blocked_share_range`** = 63\%--66\%

- Value: `[0.6276595744680851, 0.6585365853658537]`
- Formula: `min and max over (gene, drug) of blocked residues / residues`
- Computed from: range of the CBE+ABE blocked-residue share across all three genes and every drug arm
- Inputs: `{'n_arms': 7}`
- Code: `scripts/compute_paper_numbers.py::build_reachability_track.py::build`

**`kras_r68l_z_sot`** = +3.47

- Value: `3.47`
- Formula: `Feng atlas average Z`
- Computed from: chemistry-free Z-score of R68L under sotorasib
- Inputs: `{'variant': 'R68L'}`
- Code: `scripts/compute_paper_numbers.py::feng_decomposition`

**`kras_r68l_z_ada`** = +4.34

- Value: `4.34`
- Formula: `Feng atlas average Z`
- Computed from: chemistry-free Z-score of R68L under adagrasib
- Inputs: `{'variant': 'R68L'}`
- Code: `scripts/compute_paper_numbers.py::feng_decomposition`


### §3.3

**`y96d_sot_mutazyme_shift_beta_range`** = +5 to +643

- Value: `[5, 643]`
- Formula: `min and max over beta of rank(log S + beta*sigma) - rank(sigma)`
- Computed from: range of the Y96D rank shift under Mutazyme II with sotorasib across beta in [0.25, 4]
- Inputs: `{'betas': [0.25, 0.5, 1.0, 2.0, 4.0], 'shifts': [643, 508, 132, 19, 5]}`
- Code: `scripts/compute_paper_numbers.py::compute_h95d_ranks`


### §3.4

**`top100_alk_enu_alectinib`** = 63\% [53\%--72\%]

- Value: `0.63`
- Formula: `|sim top-100 ∩ chemistry-free top-100| / 100`
- Computed from: fraction of the chemistry-free top-100 that the ENU-simulated alectinib screen also ranks in its top 100, with a Wilson 95% interval
- Inputs: `{'gene': 'ALK', 'N': 1427, 'K': 100}`
- Code: `scripts/compute_paper_numbers.py::recompute_topk_with_ci.py::compute_curve`

**`top100_alk_mutii_alectinib`** = 65\% [55\%--74\%]

- Value: `0.65`
- Formula: `|sim top-100 ∩ chemistry-free top-100| / 100`
- Computed from: fraction of the chemistry-free top-100 that the MutII-simulated alectinib screen also ranks in its top 100, with a Wilson 95% interval
- Inputs: `{'gene': 'ALK', 'N': 1427, 'K': 100}`
- Code: `scripts/compute_paper_numbers.py::recompute_topk_with_ci.py::compute_curve`

**`top100_alk_enu_lorlatinib`** = 65\% [55\%--74\%]

- Value: `0.65`
- Formula: `|sim top-100 ∩ chemistry-free top-100| / 100`
- Computed from: fraction of the chemistry-free top-100 that the ENU-simulated lorlatinib screen also ranks in its top 100, with a Wilson 95% interval
- Inputs: `{'gene': 'ALK', 'N': 1427, 'K': 100}`
- Code: `scripts/compute_paper_numbers.py::recompute_topk_with_ci.py::compute_curve`

**`top100_alk_mutii_lorlatinib`** = 72\% [63\%--80\%]

- Value: `0.72`
- Formula: `|sim top-100 ∩ chemistry-free top-100| / 100`
- Computed from: fraction of the chemistry-free top-100 that the MutII-simulated lorlatinib screen also ranks in its top 100, with a Wilson 95% interval
- Inputs: `{'gene': 'ALK', 'N': 1427, 'K': 100}`
- Code: `scripts/compute_paper_numbers.py::recompute_topk_with_ci.py::compute_curve`

**`top100_alk_enu_zotizalkib`** = 63\% [53\%--72\%]

- Value: `0.63`
- Formula: `|sim top-100 ∩ chemistry-free top-100| / 100`
- Computed from: fraction of the chemistry-free top-100 that the ENU-simulated zotizalkib screen also ranks in its top 100, with a Wilson 95% interval
- Inputs: `{'gene': 'ALK', 'N': 1427, 'K': 100}`
- Code: `scripts/compute_paper_numbers.py::recompute_topk_with_ci.py::compute_curve`

**`top100_alk_mutii_zotizalkib`** = 65\% [55\%--74\%]

- Value: `0.65`
- Formula: `|sim top-100 ∩ chemistry-free top-100| / 100`
- Computed from: fraction of the chemistry-free top-100 that the MutII-simulated zotizalkib screen also ranks in its top 100, with a Wilson 95% interval
- Inputs: `{'gene': 'ALK', 'N': 1427, 'K': 100}`
- Code: `scripts/compute_paper_numbers.py::recompute_topk_with_ci.py::compute_curve`

**`alk_l1196m_sigma`** = +7.7

- Value: `7.74452093727412`
- Formula: `sigma = Resistance score, chemistry-free`
- Computed from: Oh atlas alectinib Resistance score for L1196M
- Inputs: `{'forced_class': 'C>A'}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`

**`alk_l1196m_cf_rank`** = 1

- Value: `1`
- Formula: `rank by sigma`
- Computed from: L1196M chemistry-free rank among the 1,427 class-forced ALK variants, alectinib
- Inputs: `{'forced_class': 'C>A', 'shift': 0}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`

**`alk_l1196m_sim_rank`** = 1

- Value: `1`
- Formula: `rank by log S + beta*sigma at beta = 1`
- Computed from: L1196M Mutazyme II-simulated rank among the 1,427 class-forced ALK variants, alectinib
- Inputs: `{'forced_class': 'C>A', 'shift': 0}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`

**`alk_l1198f_sigma`** = +4.6

- Value: `4.55175771640876`
- Formula: `sigma = Resistance score, chemistry-free`
- Computed from: Oh atlas alectinib Resistance score for L1198F
- Inputs: `{'forced_class': 'C>T'}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`

**`alk_l1198f_cf_rank`** = 18

- Value: `18`
- Formula: `rank by sigma`
- Computed from: L1198F chemistry-free rank among the 1,427 class-forced ALK variants, alectinib
- Inputs: `{'forced_class': 'C>T', 'shift': -15}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`

**`alk_l1198f_sim_rank`** = 3

- Value: `3`
- Formula: `rank by log S + beta*sigma at beta = 1`
- Computed from: L1198F Mutazyme II-simulated rank among the 1,427 class-forced ALK variants, alectinib
- Inputs: `{'forced_class': 'C>T', 'shift': -15}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`

**`alk_v1180l_sigma`** = +5.3

- Value: `5.26448801756145`
- Formula: `sigma = Resistance score, chemistry-free`
- Computed from: Oh atlas alectinib Resistance score for V1180L
- Inputs: `{'forced_class': 'C>G'}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`

**`alk_v1180l_cf_rank`** = 8

- Value: `8`
- Formula: `rank by sigma`
- Computed from: V1180L chemistry-free rank among the 1,427 class-forced ALK variants, alectinib
- Inputs: `{'forced_class': 'C>G', 'shift': 7}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`

**`alk_v1180l_sim_rank`** = 15

- Value: `15`
- Formula: `rank by log S + beta*sigma at beta = 1`
- Computed from: V1180L Mutazyme II-simulated rank among the 1,427 class-forced ALK variants, alectinib
- Inputs: `{'forced_class': 'C>G', 'shift': 7}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`

**`alk_c1156y_cf_rank`** = 132

- Value: `132`
- Formula: `rank by sigma`
- Computed from: C1156Y chemistry-free rank among the 1,427 class-forced ALK variants, alectinib
- Inputs: `{'forced_class': 'C>T', 'shift': -77}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`

**`alk_c1156y_sim_rank`** = 55

- Value: `55`
- Formula: `rank by log S + beta*sigma at beta = 1`
- Computed from: C1156Y Mutazyme II-simulated rank among the 1,427 class-forced ALK variants, alectinib
- Inputs: `{'forced_class': 'C>T', 'shift': -77}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`

**`alk_l1152r_cf_rank`** = 225

- Value: `225`
- Formula: `rank by sigma`
- Computed from: L1152R chemistry-free rank among the 1,427 class-forced ALK variants, alectinib
- Inputs: `{'forced_class': 'T>G', 'shift': 597}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`

**`alk_l1152r_sim_rank`** = 822

- Value: `822`
- Formula: `rank by log S + beta*sigma at beta = 1`
- Computed from: L1152R Mutazyme II-simulated rank among the 1,427 class-forced ALK variants, alectinib
- Inputs: `{'forced_class': 'T>G', 'shift': 597}`
- Code: `scripts/compute_paper_numbers.py::oh_decomposition`


### §3.6

**`lorlatinib_inter_study_spread`** = 0.32

- Value: `0.3219465077435653`
- Formula: `median over variants of max - min log10 fold change`
- Computed from: median inter-study range of measured lorlatinib potency on the same variant
- Inputs: `{'n': 18}`
- Code: `scripts/compute_paper_numbers.py::sigma_transport_validation.py`


### §3.7

**`met_wt_median_rank_nvp`** = 495

- Value: `495`
- Formula: `median rank by delta sigma, WT background`
- Computed from: median atlas rank of the NVP-BVU972 hits on the wild-type MET background
- Inputs: `{'found': 7, 'total': 10}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_wt_median_rank_amg`** = 2,472

- Value: `2472`
- Formula: `median rank by delta sigma, WT background`
- Computed from: median atlas rank of the AMG 458 hits on the wild-type MET background
- Inputs: `{'found': 11, 'total': 13}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_amg_ta_share`** = 33\%

- Value: `0.3333333333333333`
- Formula: `clones in class / total clones`
- Computed from: share of AMG 458 clones whose substitution is T>A
- Inputs: `{'n_clones': 60}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_clones_y1230c`** = 31

- Value: `31`
- Formula: `transcribed clone count`
- Computed from: clones carrying Y1230C in the Tiedt NVP-BVU972 screen
- Inputs: `{'classes': ['T>C']}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_clones_y1230d`** = 1

- Value: `1`
- Formula: `transcribed clone count`
- Computed from: clones carrying Y1230D in the Tiedt NVP-BVU972 screen
- Inputs: `{'classes': ['T>G']}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_clones_y1230f`** = 0

- Value: `0`
- Formula: `transcribed clone count`
- Computed from: clones carrying Y1230F in the Tiedt NVP-BVU972 screen
- Inputs: `{'classes': ['T>A']}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_clones_y1230h`** = 20

- Value: `20`
- Formula: `transcribed clone count`
- Computed from: clones carrying Y1230H in the Tiedt NVP-BVU972 screen
- Inputs: `{'classes': ['T>C']}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_clones_y1230n`** = 19

- Value: `19`
- Formula: `transcribed clone count`
- Computed from: clones carrying Y1230N in the Tiedt NVP-BVU972 screen
- Inputs: `{'classes': ['T>A']}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`

**`met_y1230_clones_y1230s`** = 4

- Value: `4`
- Formula: `transcribed clone count`
- Computed from: clones carrying Y1230S in the Tiedt NVP-BVU972 screen
- Inputs: `{'classes': ['T>G']}`
- Code: `scripts/compute_paper_numbers.py::met_audit.py::main`


### §Background

**`oh_atlas_total_snvs`** = 3,208

- Value: `3208`
- Formula: `len(data/raw/oh_2026/table_s2_resistance_scores.csv)`
- Computed from: total rows in Oh 2026 Table S2 resistance scores
- Inputs: `{'path': 'data/raw/oh_2026/table_s2_resistance_scores.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`oh_atlas_coding_snvs`** = 2,953

- Value: `2953`
- Formula: `sum(Missense + Synonymous + Nonsense) in table_s2_resistance_scores.csv`
- Computed from: sum of coding categories in Oh 2026 Table S2
- Inputs: `{'by_type': '{...}', 'coding_types': ['Missense', 'Nonsense', 'Synonymous']}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`oh_atlas_noncoding`** = 255

- Value: `255`
- Formula: `count(type == 'Intron') + count(type == 'Splice Acceptor/Donor')`
- Computed from: sum of intron and splice-site rows in Oh 2026 Table S2
- Inputs: `{'intron': 149, 'splice': 106}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

**`feng_atlas_total_variants`** = 3,567

- Value: `3567`
- Formula: `len(data/raw/feng_2022/table_s1_zscores.csv)`
- Computed from: rows in Feng 2022 KRAS G12C Twist atlas Table S1
- Inputs: `{'path': 'data/raw/feng_2022/table_s1_zscores.csv'}`
- Code: `scripts/compute_paper_numbers.py::compute_atlas_sizes`

