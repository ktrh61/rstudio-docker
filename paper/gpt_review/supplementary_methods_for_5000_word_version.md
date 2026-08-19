# Supplementary Methods for the revised 5,000-word editorial trial

> Prepared on 2026-08-19 from the Supplementary Methods and method details in paper/draft_manuscript.md. This file is an editorial reallocation. A separate descriptive audit of finalized outputs is identified explicitly below; it did not rerun the inferential pipeline.

## Data sources and expression matrix

Gene-level STAR count files for the REBC-THYR project were obtained from the NCI Genomic Data Commons and md5-verified against the download manifest. The clinical table was Data S1 of Morton et al. (2021), read with all columns retained and no value edits. Missing markers were converted to NA, and a column was typed as numeric only when every non-missing value parsed as numeric.

Files were mapped to cases and biospecimens through the GDC API. One count assay per sample was assembled into a single matrix. All-zero genes were removed, yielding 58,448 genes and 906 samples. Exon-union gene lengths were derived from the GDC GENCODE v36 reference annotation.

Library strandedness was selected per sample from the totals of the two STAR stranded columns. A library was called stranded when the smaller total was no more than half the larger, in which case the larger column was used; otherwise the unstranded column was used. All 906 libraries were classified as reverse stranded, with smaller-to-larger ratios of 0.056–0.110 against the 0.5 threshold. GDC TPM and FPKM columns were not used. For the REO analysis, TPM was recalculated from the selected count assay and exon-union lengths.

## Exposure metric: Assigned Share

Assigned Share associated with the expected value of the excess relative risk was calculated with NIH IREP version 5.7.3. The fixed inputs were electrons E>15 keV, acute exposure rate, organ dose entered as a constant in cSv=mGy/10, recorded sex, exposure year 1986, birth year calculated as 1986 minus age at exposure, and surgery year calculated as birth year plus age at surgery. Remaining settings were the program defaults: user-defined uncertainty distribution Lognormal(1,1), 10,000 iterations, and random-number seed 99. Dose and age inputs in the versioned AS file were checked case by case against the clinical table.

IREP guidance identifies iodine-131 as an example for the electrons E>15 keV route. This route has a radiation effectiveness factor of 1.0 relative to reference higher-energy photons (Land et al. 2003), with independent microdosimetric support for similar biological effectiveness of internal iodine-131 and external gamma rays per absorbed dose (Sato et al. 2014). This addresses radiation quality, not the transfer of a thyroid risk model fitted to external photon exposures. The IREP thyroid model is based on atomic-bomb survivors and populations irradiated as children for medical indications (Ron et al. 1995; Kocher et al. 2008).

The acute convention followed two previous IREP applications in post-Chernobyl papillary carcinoma (Zurnadzhy et al. 2022; Bogdanova et al. 2022). In IREP, exposure rate affects the calculation through the dose and dose-rate effectiveness factor. Under the acute setting it approaches 1 at and above an uncertain reference dose distributed log-uniformly from 0.03 to 0.2 Gy; the High band's minimum dose was 188 mGy. A chronic calculation would reduce the numeric AS scale. Because the reported analysis uses an ordinal partition rather than a calibrated probability, alternative settings concern the dependence of the scale and band definition on the convention; they do not establish that the model transfer is correct.

Dose was entered as a point value, and dosimetric or risk-model uncertainty was not propagated into band membership. Dose-zero cases formed the separate reference group. Exposed cases were classified as Low (0<AS<33.3), Mid (33.3≤AS<66.6), and High (AS≥66.6); no case lay on a boundary. The composition of all 440 cases by driver, band, and pairing status is in Table S1.

## Quality control and analysis cohorts

The outlier procedure operated separately within each observed group-by-tissue submatrix after filterByExpr. It used unnormalized log counts per million with prior count 2. The most extreme sample on the first principal component was removed iteratively using the two-sided PC-OD rule at alpha=0.05 until no further sample was rejected (Nakayama et al. 2024). PC-OD was a QC and cohort-definition step, not a gene-level hypothesis test, and its finalized flags were held fixed in downstream inference. It flagged no RET tumor or normal sample; its only flag among the four main target groups was one High-AS BRAF tumor. Thus, it did not alter the primary RET-tumor cohort.

Relative tumor purity used the maximum-one proportion estimator from contamDE-lm, implemented in house for the proportion step only, without its differential-expression model (Shen et al. 2016; Ji et al. 2020). Protein-coding, filterByExpr-reduced paired counts were MUREN-normalized before estimation. Both dose-zero and High-AS samples were estimated jointly within each driver cohort, producing a common relative scale. Purity was treated as a relative score within the jointly estimated cases, not an absolute cellular fraction.

Case pairing used the GDC merged-aliquot expression assays. Each included case was required to have one unique merged tumor and one unique merged normal sample. The main cohort additionally required RET or BRAF classification, membership in the dose-zero or High-AS group, both tissues passing the outlier screen, and relative purity at least 0.6.

The REO training set comprised the RET subset of the main cohort. The evaluation set comprised RET tumors from the Low-AS and Mid-AS bands, whether or not a matched normal was available. Training outlier and purity exclusions were not applied to this evaluation set; corresponding metrics were reported as ancillary diagnostics and did not authorize exclusions.

## Analysis contrasts

The four High-AS versus dose-zero contrasts were RET-tumor (`R_Tumor` in the analysis files), RET-normal (`R_Normal`), BRAF-tumor (`B_Tumor`), and BRAF-normal (`B_Normal`). These labels refer to contrasts, not sample groups. In each contrast, group x was dose-zero and group y was High-AS. The Brunner–Munzel relative effect was

p = P(X<Y) + 0.5 P(X=Y).

Thus, p>0.5 indicated higher expression in High-AS cases. The signed effect used in Figure 2 was 2p−1. The RET-tumor Higher Criticism test was the single primary contrast-level test. The other three contrasts addressed separate secondary questions. Their p-values were not combined into a study-level decision, so no across-contrast adjustment was applied; all contrasts were reported, and no study-wide FDR or across-contrast family-wise error claim was made.

## Candidate confounders

The difference in age at surgery was summarized in each driver stratum with the Hodges–Lehmann median difference and the Brunner–Munzel relative effect P(Sporadic<High). Each estimate received a 95% percentile bootstrap interval from 9,999 within-group resamples using seed 19450809. No p-value was calculated. Age at exposure was undefined for the dose-zero group and was not compared between groups.

Sex was reported by group because sex-chromosome genes can differ with group composition. All selected genes were annotated for chromosome membership against GENCODE v36. RET fusion-partner counts were reported by band; the BRAF stratum contained BRAF V600E by construction. These summaries disclose possible covariate structure but do not test or remove confounding.

The finalized analysis objects contained no sequencing-batch, processing-batch, centre, or collection-period field suitable for a group-balance calculation. This unavailable information is treated as an unassessed source of technical confounding. The added primary RET-tumor audit reports the available relative-purity, age, sex, and fusion-partner variables without a second hypothesis test.

## Normalization

For each contrast, protein-coding genes were reduced with filterByExpr and normalized with an in-house three-iteration implementation of DEGES-MUREN (Kadota et al. 2012; Sun et al. 2013; Feng and Li 2021). In this two-group use, filterByExpr uses the group sizes to set the required number of expressed samples but does not use the identities of the samples meeting that expression rule; for fixed group sizes, the retained set is therefore invariant to relabelling. The MUREN step used its single-parameter pairwise least-trimmed-squares form. At each iteration, a permutation Brunner–Munzel screen produced Storey q-values with lambda=0.5. The larger of the q<0.10 set and the top 5% of genes by raw p-value was removed from scaling-factor estimation, implementing the floorPDEG guard. Tested-gene counts, screen pi0 estimates, Jaccard convergence, and scaling-factor ranges are in Table S2.

The DEGES exclusion screen uses the observed group labels, whereas the reported gene-level allocations hold the resulting normalized matrix fixed. Accordingly, “exact” refers to exhaustive enumeration of labels on the finalized matrix, not to end-to-end randomization of QC and normalization from raw counts. PC-OD is not part of the gene test and, in the primary RET cohort, removed no sample. The remaining scope qualification concerns the group-aware DEGES normalization and the exchangeability assumption for the retained samples.

## Gene-level differential expression

For each gene, the studentized Brunner–Munzel statistic was evaluated for every allocation of nx labels among n samples on the finalized normalized matrix. Exhaustive enumeration removed random Monte Carlo error and the 1/(B+1) resolution floor. It did not, by itself, remove the exchangeability condition of the permutation null or turn upstream QC and normalization into part of the enumerated procedure. Ties were retained in the effect estimate and statistic.

Gene-level p-values were converted to Storey q-values with the plug-in estimate of pi0 at fixed lambda=0.5 and threshold q<0.10 (Storey 2002). The fixed value avoided adaptive tuning. This definition specifies how the ranked gene list was constructed; it is not presented as a proof of FDR control under arbitrary gene dependence. We therefore refer to “genes meeting the prespecified Storey q rule” and reserve the contrast-level existence statement for the Higher Criticism label-permutation result.

## Descriptive magnitude and covariate audit

The additional audit read the finalized RET-tumor DGEList, gene-level result table, cohort table, and clinical table. It did not rerun QC, normalization, the Brunner–Munzel tests, Storey adjustment, or gene-set analysis. For the existing q<0.10 list, it summarized |2p−1| and the absolute mean log2 fold change. The latter was calculated for display only from DEGES-normalized log2 CPM with prior count 1 as the High-AS-group mean minus the dose-zero-group mean.

For sample-level context, principal components were calculated from all 15,621 tested genes and all 27 RET-tumor samples. Log2 CPM values were centered gene by gene and were not variance-scaled; no gene was selected by q-value or group difference for this plot. The same PC1–PC2 coordinates were displayed by group and sex, relative tumor-purity score, age at surgery, and RET fusion partner. Spearman correlations of PC1 and PC2 with age and purity were reported without p-values. No adjusted gene model was fitted; the available variables were reported descriptively because combining AS components, expression-derived purity, and potential biological intermediates as nuisance terms would define a different conditional estimand.

For contrast-level diagnostics, 9,999 label shuffles were generated with seed 19860426 and stored. Applying the plug-in estimator to each shuffle gave a null reference distribution for pi0. These same shuffles supplied the omnibus and gene-set references.

## Contrast-level omnibus inference

Higher Criticism was the prespecified primary omnibus statistic, with scan range alpha0=0.1 (Donoho and Jin 2004). Its p-value was calculated against the contrast's own label-shuffle distribution, so the analytic independence and sparsity assumptions of the original Higher Criticism null were not used. The RET-tumor Higher Criticism test was the single primary contrast-level test; the corresponding statistics in the other contrasts answered separate secondary questions. Count statistics at fixed cutoffs and the maximum statistic were retained as descriptive rows. The full rejection curve R(alpha) was retained to show how the displayed gene count varied across q thresholds.

Because no binary contrast label was defined, the omnibus p-values, per-gene q-values, and rejection curves should be reported as continuous evidence. “Support” and “no support” wording should not imply an undeclared cutoff.

## Gene-set level inference

Genes were ranked by tie-averaged normal scores of the signed Brunner–Munzel statistic. The weighted running sum used gseaParam=1 and was evaluated only at tie-block boundaries. Automated tests verified equality to the standard GSEA statistic for tie-free input (Subramanian et al. 2005).

The null comprised the 9,999 saved label shuffles per contrast. For each set, a sign-conditional permutation p-value was calculated and adjusted by Benjamini–Hochberg within collection. pi0 was fixed at 1 for the set-level analysis. The four families were Hallmark; C2 canonical pathways restricted to Reactome, WikiPathways, KEGG MEDICUS, BioCarta, and PID; C5 GO Biological Process; and a radiation-curated subset of C2:CGP. The radiation curation rule was fixed before the reported set-level run. Sets outside 15–500 genes were excluded.

The spike-in multiplied the 195 genes present from HALLMARK_ADIPOGENESIS by 1.15 in the nine High-AS BRAF-tumor samples. This is a single coherent-signal check and not a general power analysis.

For over-representation analysis, RET-tumor genes meeting Storey q<0.10 were split into higher, lower, and combined lists. One-sided hypergeometric tests used all genes tested in the RET-tumor contrast as the universe and BH adjustment within family-by-list. In the terminology of Goeman and Bühlmann (2007), the statistic is competitive, but its gene-sampling p-value does not refer to subject-level label randomization.

## Complete-null assessment of the set-level procedure

For each contrast, labels were shuffled B+R times with seeds derived from 19450809, using B=9,999 as a shared reference pool and R=100 held-out pseudo-observations. Each pseudo-observation was normalized against the pool, assigned per-set permutation p-values, adjusted within collection, and thresholded at q<0.10.

For each contrast-by-collection cell, the proportion of held-out global-null replicates yielding at least one discovery estimated P(at least one false discovery) under the complete null, where FDR and family-wise error coincide. Exact binomial intervals and mean discovery counts were reported in Table S7. The 100 indicators shared one null pool and were positively correlated; ordinary binomial intervals are therefore somewhat narrow.

This assessment informed selection of per-set permutation p-values with within-collection BH before the reported real-data run. An earlier pooled tail-ratio procedure yielded pooled complete-null rates of 0.140 and 0.221 for two variants, with a worst cell of 0.44, and was rejected. The adopted procedure yielded 0.045 in its development assessment. In the final complete-null assessment, the pooled rate was 0.064, with one cell-specific rate of 0.18. These results characterize complete-null behaviour; they do not prove FDR control for every partial-alternative configuration.

## Between-stratum concordance

For each tissue, signed Brunner–Munzel profiles were intersected by gene and compared between RET and BRAF with Spearman correlation. Labels were shuffled independently within each contrast 9,999 times, using contrast-specific seeds derived from 19450809, and the shuffled profiles were correlated to form the null reference interval. This interval represents the case in which neither stratum carries label-aligned structure. An observed correlation outside the interval does not distinguish exposure from a covariate shared across strata.

## REO panel construction and evaluation

TPM was recalculated from the selected assay and exon-union lengths. The candidate pool was the 500 genes with greatest absolute distance of the Brunner–Munzel effect from 0.5. The q<0.10 list was not used to define the pool.

For a gene pair, r was the within-sample log2-TPM difference. The dead zone was |r|<log2(1.2). A pair qualified when its sign was stable in dose-zero samples, with at most one exception among samples outside the dead zone; the 10th percentile of |r| was at least log2(1.5); and the pair reversed in more than 50% but fewer than 100% of High samples. Qualified pairs were ranked by the shift in median r. Greedy selection prohibited gene reuse and excluded a pair whose r profile had Spearman correlation at least 0.75 with a retained pair. The target was 10 pairs.

The reversal score counted pairs outside the dead zone that opposed the dose-zero reference sign. The classification boundary was above the maximum score observed among dose-zero training samples. Neither the panel nor its boundary was changed after construction.

The only prespecified independent comparison was Mid over Low by a one-sided Monte Carlo Brunner–Munzel test, using seed 19860426. Scores for the dose-zero and High training bands were descriptive and were not estimates of out-of-sample performance.

## REO diagnostics

The PC-OD screen used for training was mirrored on Low and Mid tumors, but its output was a count only and could not exclude evaluation cases. Relative purity was estimated on a common scale by pooling paired RET cases across all bands in one contamDE-lm proportion run. Purity was available for 31 of the 36 evaluation tumors. The remaining two Low and three Mid cases had a tumor assay but no matched normal assay, which the pair-based purity estimator requires; they remained in the REO evaluation.

Band and score were compared conditional on purity with a partial Spearman correlation and a one-sided permutation reference. Mid and Low were also compared within two strata split at median purity. These analyses were diagnostics and did not establish independence from purity.

## External anchor cross-reference

The cross-reference counted k of n source-list genes present among each contrast's tested genes and q<0.10 list. Abend et al. (2013) contributed separate tumor and normal lists. Abend et al. (2012) contributed a tumor-minus-normal paired-difference list, which had no matching estimand in the present contrasts. Dom et al. (2012) contributed a normal-tissue list, and CLIP2 was included as a single-gene anchor. Because the lists were small and differed in tissue, platform, and contrast, no enrichment statistic was calculated.

## Software, seeds, and reproducibility

The container used the immutable Ubuntu noble-20260410 base and matching apt snapshot, R 4.5.3 compiled against reference BLAS, and R packages from the 2026-04-09 CRAN and Bioconductor 3.22 snapshots. The bit-level reproduction contract targeted linux/amd64.

In-house components included exact-enumeration and Monte Carlo Brunner–Munzel statistics, the Storey plug-in estimator, tie-block enrichment, DEGES-MUREN normalization, the contamDE-lm proportion score, the principal-component outlier procedure, and REO panel construction. External packages included edgeR, SummarizedExperiment, msigdbr, limma, GenomicDataCommons, rtracklayer, Rcpp, and MASS; full versions are in Table S5.

The canonical inference seed was 19860426. Diagnostics used documented seeds derived from 19450809. The full pipeline was executed on two x86-64 machines with 1,819 md5-identical raw data files. Both runs passed 415 tests without failures and produced identical primary artifacts.
