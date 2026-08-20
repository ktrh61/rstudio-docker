# Supplementary Results for the revised 5,000-word editorial trial

> Prepared from the results and verified numerical records underlying `paper/draft_manuscript.md`. This file restores supporting results omitted during compression of the main manuscript; it does not introduce a new analysis or supersede the canonical records in `paper/`.

## Gene-set analysis and complete-null assessment

After applying the 15–500-gene size filter, 6,141–6,242 sets were tested per contrast across the four collections. No set met the within-collection Benjamini–Hochberg threshold q<0.10 in any of the 16 contrast-by-collection cells. The smallest adjusted value was q=0.114 in the BRAF-tumor radiation-curated cell.

Across the held-out complete-null assessment, the proportion of replicates producing at least one discovery ranged from 0.01 to 0.18 across the 16 cells, and 13 cells were at or below the nominal 0.10 level. The BRAF-tumor radiation-curated cell had a rate of 0.10 (descriptive 95% Clopper–Pearson interval, 0.049–0.176), and the BRAF-normal radiation-curated cell had a rate of 0.12 (0.064–0.200). The only descriptive interval lying entirely above 0.10 was for BRAF-normal/Hallmark, which had a rate of 0.18 (0.110–0.270). Descriptively pooled across cells, 102 of 1,600 replicates produced at least one discovery (0.064). Because the 100 pseudo-observations within a cell shared one null pool, their indicators were dependent; the Clopper–Pearson intervals treat them as independent and may therefore be too narrow.

In the single positive-control scenario, 195 genes from HALLMARK_ADIPOGENESIS were multiplied by 1.15 in the nine High-AS BRAF-tumor samples. The planted set had normalized enrichment score 2.28, permutation p=0.0002, and q=0.0101, ranking first among 50 Hallmark sets. No other Hallmark set met q<0.10. This result documents recovery in one coherent-signal configuration and is not a general power analysis.

## Over-representation annotation of the RET-tumor list

The over-representation analysis was a descriptive annotation under a gene-sampling reference, distinct from the subject-label permutation analysis above. Across Hallmark, C2 canonical pathways, C5 GO Biological Process, and the radiation-curated family, respectively, the numbers of sets meeting the family-by-list q<0.10 rule were 0/0/0/0 for genes expressed more highly in High-AS cases, 12/105/205/7 for genes expressed less highly in High-AS cases, and 6/37/71/2 for the combined list.

The lower-expression list was concentrated in a single correlated theme involving proliferation, cell cycle, and DNA repair. Illustrative counts were 46 of 199 genes for E2F targets, 41 of 198 for G2M checkpoint, and 68 of 322 for Reactome DNA repair. The leading radiation-curated set contained 46 of 126 genes in the lower-expression list, compared with 6.4 expected under the gene-sampling reference. The radiation-curated findings therefore reflected the same cell-cycle theme rather than an independent line of evidence. These annotations remained hypothesis-generating and did not alter the absence of a q<0.10 result in the label-permutation gene-set analysis.

## Sex-chromosome annotation of the RET-tumor list

Among the 1,765 RET-tumor genes meeting the prespecified Storey q<0.10 rule, 57 were X-linked: 36 were expressed more highly and 21 less highly in High-AS cases. No Y-linked gene met the rule. The 15,621 genes tested in the contrast included 572 X-linked and two Y-linked genes. These counts were descriptive; no enrichment test was performed, and they do not exclude sex-related expression structure.

## REO pair construction

The 500-gene construction pool comprised 317 genes expressed more highly and 182 expressed less highly in High-AS RET fusion-positive PTCs. From 57,694 candidate pairs, 153 met all prespecified stability, separation, and reversal criteria. Across these qualifying pairs, the shift in median within-sample log2-TPM difference ranged from 1.159 to 4.700, and the High-AS reversal rate ranged from 0.53 to 0.87. Greedy selection with the prespecified gene-reuse and redundancy restrictions produced the final 10-pair panel. The maximum reversal score among dose-zero construction cases was 2; consequently, scores greater than 2 were classified as positive.

The construction-band classifications and intermediate-band results are reported in the main text. Their interpretation remains asymmetric: construction-band separation establishes that the selected panel represents its source contrast but is not an unbiased estimate of classification performance; application to Low-AS and Mid-AS cases did not constitute independent external validation.

## REO quality and purity diagnostics

The mirrored PC-OD screen flagged no intermediate-band tumor (0 of 17 Low-AS and 0 of 19 Mid-AS), so all 36 eligible tumors remained in the REO application. Relative purity could be estimated for 15 Low-AS and 16 Mid-AS tumors; two Low-AS and three Mid-AS tumors lacked the matched normal assay required by the pair-based estimator and remained in the REO analysis without a purity score.

In the pooled paired RET purity run, median relative-purity scores were 0.690, 0.704, 0.739, and 0.814 in the dose-zero, Low-AS, Mid-AS, and High-AS bands, respectively. This diagnostic refitted the relative scale across all paired RET bands and is therefore not numerically interchangeable with the main-cohort purity scale. Among the intermediate bands, Spearman correlations between purity and REO score were 0.684 in Low-AS cases, 0.345 in Mid-AS cases, and 0.538 when pooled. The unadjusted band–score correlation was 0.142, the band–purity correlation was 0.036, and the partial Spearman coefficient after adjustment for purity rank was 0.146 (descriptive one-sided permutation p=0.2162).

After splitting the intermediate-band cases at median purity, the one-sided Mid-AS-over-Low-AS comparison gave Pr(Low<Mid)=0.580 and p=0.3069 in the lower-purity stratum and Pr(Low<Mid)=0.532 and p=0.4197 in the higher-purity stratum. These were diagnostics, not additional confirmatory tests, and they did not establish an AS-band association independent of purity.

## External gene-list cross-reference

Nineteen of 20 source-list-by-contrast cells had no overlap, including every tissue-matched comparison. The sole overlap was cross-tissue: S100A10 from the Dom et al. normal-tissue list occurred in the RET-tumor q<0.10 list, with expression lower in High-AS cases, opposite to the direction in the original normal-tissue report. The result remained descriptive and did not modify any inferential claim.
