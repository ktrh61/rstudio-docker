# Supplementary-file assembly map

This directory contains submission-assembly copies made from frozen analysis outputs. The formatting scripts read existing RDS objects; they do not rerun quality control, normalization, label shuffles, gene-level tests, or gene-set tests.

| Manuscript item | File | Source or formatting script | Status |
|---|---|---|---|
| Figure S1 | `figure_s1_gene_level_ma.png` | `figures/fig_ma_gene_bm.R` / existing `output/figures/fig_ma_gene_bm.png` | Available |
| Figure S2 | `figure_s2_set_level_null_calibration.png` | `figures/fig_d6_calibration.R` | Available |
| Table S1 | `table_s1_cohort_composition.csv` | `tables/tab_cohort_composition.R` | Available |
| Table S2 | `table_s2_normalization_diagnostics.csv` | `tables/supp_tab_normalization_diagnostics.R` | Available |
| Table S3 | — | Full descriptive ORA table object was not retained | Assembly pending; caption and text references are prepared, but the table was not reconstructed from prose |
| Table S4 | `table_s4_reo_panel.csv` | Frozen `output/reo_panel.csv` from `scripts/510_select_reo_pairs.R` | Available |
| Table S5 | `table_s5_software_versions.csv` | `tables/supp_tab_package_versions.R` | Available |
| Table S6 | `table_s6_gene_set_summary.csv` | `tables/tab_set_level_summary.R` | Available |
| Table S7 | `table_s7_complete_null_calibration.csv` | `tables/tab_set_level_summary.R` | Available |
| Table S8 | `table_s8_between_stratum_concordance.csv` | `tables/supp_tab_concordance.R` | Available |
| Supplementary Data 1 | `supplementary_data_1_set_level_results.csv` | `tables/supp_data_420_full.R` | Available |

The legends are collected in `../supplementary_results_for_5000_word_version.md`. Before submission, Table S3 must be exported from the final ORA result object and checked against the reported family-by-list counts. The image and table files in this directory are not a journal-formatted combined supplement; they are the mapped source files for that assembly step.
