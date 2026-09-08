# Baseline 2026-08-07 (pre-Storey-migration reference)

Extracted by `tests/reference/extract_baseline_20260807.R` from the processed/
artifacts as they stood before the phase-3 refactor and phase-4 Storey migration
(reorg plan v2). A full uncommitted snapshot of the source RDS files sits in
`processed/_baseline_20260807/`; `checksums.csv` ties this directory to it.

Two kinds of files, with different contracts:

- **Case-ID invariants** (`cases_main_bm.csv`, `unit_samples.csv`,
  `reo_eval_cases.csv`, `reo_panel.csv`): the phase-3 refactor must reproduce
  these exactly. Any change is a defect, not an improvement, until the
  divergent case and cause are explained.
- **Statistical audit records** (`gene_level_*.csv`, `gsea_12cell.csv`,
  `reo_summary.csv`): the OLD values for the phase-4 amendment diff. The Storey
  migration intentionally changes them; agreement is NOT a pass criterion.

Known mixed state, recorded as-is: `thyr_normalized_counts.rds` was regenerated
without Cook's removal (2026-07-24), while `thyr_expression_test.rds` and
`thyr_enrichment_test.rds` still rest on the earlier with-Cook's normalization
(verified near-identical at the signal level; see repo history).
