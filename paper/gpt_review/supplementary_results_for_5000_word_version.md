# Supplementary Results for the revised 5,000-word editorial trial

> Prepared from the results and verified numerical records underlying `paper/draft_manuscript.md`. This file restores supporting results omitted during compression of the main manuscript; it does not introduce a new analysis or supersede the canonical records in `paper/`.

## Gene-set analysis and complete-null assessment

After applying the 15–500-gene size filter, 6,141–6,242 sets were tested per contrast across the four collections. No set met the within-collection Benjamini–Hochberg threshold q<0.10 in any of the 16 contrast-by-collection cells. The smallest adjusted value was q=0.114 in the BRAF-tumor radiation-curated cell (Table S6; Supplementary Data 1).

Across the held-out complete-null assessment, the proportion of replicates producing at least one discovery ranged from 0.01 to 0.18 across the 16 cells, and 13 cells were at or below the nominal 0.10 level. The BRAF-tumor radiation-curated cell had a rate of 0.10 (descriptive 95% Clopper–Pearson interval, 0.049–0.176), and the BRAF-normal radiation-curated cell had a rate of 0.12 (0.064–0.200). The only descriptive interval lying entirely above 0.10 was for BRAF-normal/Hallmark, which had a rate of 0.18 (0.110–0.270). Descriptively pooled across cells, 102 of 1,600 replicates produced at least one discovery (0.064). Because the 100 pseudo-observations within a cell shared one null pool, their indicators were dependent; the Clopper–Pearson intervals treat them as independent and may therefore be too narrow (Fig. S2; Table S7).

In the single positive-control scenario, 195 genes from HALLMARK_ADIPOGENESIS were multiplied by 1.15 in the nine High-AS BRAF-tumor samples. The planted set had normalized enrichment score 2.28, permutation p=0.0002, and q=0.0101, ranking first among 50 Hallmark sets. No other Hallmark set met q<0.10. This result documents recovery in one coherent-signal configuration and is not a general power analysis.

## Over-representation annotation of the RET-tumor list

The over-representation analysis was a descriptive annotation under a gene-sampling reference, distinct from the subject-label permutation analysis above. Across Hallmark, C2 canonical pathways, C5 GO Biological Process, and the radiation-curated family, respectively, the numbers of sets meeting the family-by-list q<0.10 rule were 0/0/0/0 for genes expressed more highly in High-AS cases, 12/105/205/7 for genes expressed less highly in High-AS cases, and 6/37/71/2 for the combined list (Table S3).

The lower-expression list was concentrated in a recurring, overlapping theme involving proliferation, cell cycle, and DNA repair. Illustrative counts were 46 of 199 genes for E2F targets, 41 of 198 for G2M checkpoint, and 68 of 322 for Reactome DNA repair. The radiation-curated set with the smallest raw hypergeometric p-value contained 46 of 126 genes in the lower-expression list, compared with 6.4 expected under the gene-sampling reference. The radiation-curated findings therefore reflected the same cell-cycle theme rather than an independent line of evidence. These annotations remained hypothesis-generating and did not alter the absence of a q<0.10 result in the label-permutation gene-set analysis.

## Sex-chromosome annotation of the RET-tumor list

Among the 1,765 RET-tumor genes meeting the prespecified Storey q<0.10 rule, 57 were X-linked: 36 were expressed more highly and 21 less highly in High-AS cases. No Y-linked gene met the rule. The 15,621 genes tested in the contrast included 572 X-linked and two Y-linked genes. These counts were descriptive; no enrichment test was performed, and they do not exclude sex-related expression structure.

## REO pair construction

The top-500 construction ranking contained 318 genes expressed more highly and 182 expressed less highly in High-AS RET fusion-positive PTCs. One higher-expression gene had zero TPM in at least one construction sample and therefore failed the finite-log2-TPM requirement, leaving 317 higher- and 182 lower-expression genes eligible for pairing. Of their 57,694 cross-direction pairs, 153 met all prespecified stability, separation, and reversal criteria. Across these qualifying pairs, the absolute shift in median within-sample log2-TPM difference ranged from 1.159 to 4.700, and the High-AS reversal rate ranged from 0.53 to 0.87. Greedy selection with the prespecified gene-reuse and redundancy restrictions produced the final 10-pair panel. The maximum reversal score among dose-zero construction cases was 2; consequently, scores greater than 2 were classified as positive (Table S4).

The construction-band classifications and intermediate-band results are reported in the main text. Their interpretation remains asymmetric: construction-band separation establishes that the selected panel represents its source contrast but is not an unbiased estimate of classification performance; application to Low-AS and Mid-AS cases did not constitute independent external validation.

## Ancillary REO quality and purity diagnostics

The ancillary PC-OD screen identified no outliers among the intermediate-band tumors (Low-AS, 0/17; Mid-AS, 0/19). Independently of this result, the screen was non-exclusionary, and all 36 tumors remained in the REO application.

Relative purity could be estimated for 15 Low-AS and 16 Mid-AS tumors; two Low-AS and three Mid-AS tumors lacked the matched normal assay required by the pair-based estimator and remained in the REO analysis without a purity score. In the pooled paired RET purity run, median relative-purity scores were 0.690, 0.704, 0.739, and 0.814 in the dose-zero, Low-AS, Mid-AS, and High-AS bands, respectively. This diagnostic refitted the relative scale across all paired RET bands and is therefore not numerically interchangeable with the main-cohort purity scale.

Among the intermediate bands, Spearman correlations between purity and REO score were 0.684 in Low-AS cases, 0.345 in Mid-AS cases, and 0.538 when pooled. The unadjusted band–score correlation was 0.142, the band–purity correlation was 0.036, and the partial Spearman coefficient after adjustment for purity rank was 0.146 (descriptive one-sided permutation p=0.2162).

## External gene-list cross-reference

Among the validated anchors, 19 of 20 source-list-by-contrast cells had no overlap. The sole overlap was cross-tissue: S100A10 from the Dom et al. normal-tissue list occurred in the RET-tumor q<0.10 list, with expression lower in High-AS cases, opposite to the direction in the original normal-tissue report.

Five of the 39 Ory shared-tissue signature genes tested in each contrast occurred in the RET-tumor q<0.10 list: ATP5MF, MRPL52, NTHL1, URM1, and USE1. None occurred in RET-normal, BRAF-tumor, or BRAF-normal. Three of 40 tested genes from the Ory normal-tissue signature occurred in RET-tumor—PXDN, S100A10, and TESC—but none occurred in either normal-tissue contrast or in BRAF-tumor; the RET-tumor overlaps were therefore cross-tissue. Three of 46 tested genes from the Ory tumor-tissue signature occurred in RET-tumor—P2RY1, PLK2, and EHD4—whereas none of the 48 tested genes occurred in BRAF-tumor. No member of the tumor-tissue signature occurred in either normal-tissue q<0.10 list. These results were descriptive membership counts rather than enrichment tests and did not constitute validation of the multivariate signatures.

## Supplementary figure and table legends

**Figure S1 | Gene-level MA plots by analysis contrast.** The A axis is mean log2 counts per million across samples in each contrast, and the M axis is the display-only log2 fold change for High-AS minus dose-zero cases (labeled “Sporadic” in the analysis code). Points are colored by the Storey q<0.10 result from the rank-based Brunner–Munzel analysis; the fold change was not used for gene selection.

**Figure S2 | Held-out complete-null calibration of the gene-set procedure.** Each row shows, for one contrast-by-collection cell, the proportion of 100 held-out label-shuffle pseudo-observations that produced at least one discovery at q<0.10. Error bars are descriptive 95% Clopper–Pearson intervals, and the dashed line marks the nominal 0.10 level. Because pseudo-observations within a cell shared the same null pool, the intervals do not account for their dependence.

**Table S1 | Cohort composition by molecular classification, AS band, and sample pairing.** Counts are shown for all 440 cases using both the detailed driver classification and the broader driver-group classification. For each AS band, the table gives total, paired, and unpaired counts; absent combinations are displayed as zero.

**Table S2 | Expression filtering and DEGES-MUREN normalization diagnostics by contrast.** The table reports group sizes, the number of protein-coding genes before and after `filterByExpr`, the number of DEGES iterations, the third-iteration screening π0 estimate and Jaccard index, and the range of final normalization factors.

**Table S3 | Descriptive over-representation annotation of the RET-tumor gene list.** Results are reported separately for genes expressed more highly in High-AS cases, genes expressed less highly in High-AS cases, and their union across Hallmark, C2 canonical pathways, C5 Gene Ontology Biological Process, and the radiation-curated family. Set size, observed and expected overlap, one-sided hypergeometric p-value, and family-by-list Benjamini–Hochberg q-value are included. This table is hypothesis-generating and does not report the subject-label permutation gene-set test.

**Table S4 | Ten-pair relative-expression-ordering panel.** For each selected pair, the table gives Ensembl identifiers and gene symbols for the higher- and lower-direction genes, the absolute construction-band shift in median within-sample log2-TPM difference, the High-AS reversal rate, and the 10th percentile of the absolute within-sample difference among dose-zero construction cases.

**Table S5 | Software versions.** Versions of R packages used in the pinned publication container are listed.

**Table S6 | Gene-set results summarized by contrast and collection.** The number of sets passing the 15–500-gene size filter and the minimum within-collection Benjamini–Hochberg q-value are shown for each of the 16 contrast-by-collection cells.

**Table S7 | Held-out complete-null calibration by contrast and collection.** For each cell, the table reports the number of tested sets, number of pseudo-observations, number and proportion with at least one discovery at q<0.10, descriptive 95% Clopper–Pearson limits, mean number of discoveries, and maximum number in any pseudo-observation.

**Table S8 | Between-stratum concordance of signed gene-level statistics.** Normal-tissue and tumor-tissue comparisons are shown separately with the paired contrasts, number of shared genes, Spearman correlation, central 95% label-shuffle interval, two-sided shuffle p-value, and number of shuffles.

**Supplementary Data 1 | Complete label-permutation gene-set results.** All tested sets are provided for each contrast and collection with set size, enrichment score, normalized enrichment score, permutation p-value, within-collection Benjamini–Hochberg q-value, redundancy annotation, and leading-edge genes.
