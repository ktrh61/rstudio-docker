# Supplementary Methods for the revised 5,000-word editorial trial

> Prepared on 2026-08-19 from the Supplementary Methods and method details in paper/draft_manuscript.md. This file is an editorial reallocation and does not supersede the canonical files in paper/.

## Data sources and expression matrix

We obtained gene-level STAR count files for the REBC-THYR project from the NCI Genomic Data Commons and verified them against the download manifest using MD5 checksums. Clinical data came from Data S1 of Morton et al. (2021), which we read with all columns retained and without editing any values. We converted missing markers to NA and typed a column as numeric only when every non-missing value parsed as numeric.

We mapped files to cases and biospecimens through the GDC API and assembled one count assay per sample into a single matrix. Removing all-zero genes yielded 58,448 genes and 906 samples. We derived exon-union gene lengths from the GDC GENCODE v36 reference annotation.

For each sample, we selected library strandedness from the totals of the two STAR stranded columns. A library was classified as stranded when the smaller total was no more than half the larger, in which case the larger column was used; otherwise the unstranded column was used. All 906 libraries were classified as reverse stranded, with smaller-to-larger ratios of 0.056–0.110 against the 0.5 threshold. GDC TPM and FPKM columns were not used. For the REO analysis, we recalculated TPM from the selected count assay and exon-union lengths.

## Exposure metric: Assigned Share

We calculated Assigned Share associated with the expected value of the excess relative risk using NIH IREP version 5.7.3. The fixed inputs were electrons E>15 keV, acute exposure rate, organ dose entered as a constant in cSv=mGy/10, recorded sex, exposure year 1986, birth year calculated as 1986 minus age at exposure, and surgery year calculated as birth year plus age at surgery. All other settings remained at the program defaults: user-defined uncertainty distribution Lognormal(1,1), 10,000 iterations, and random-number seed 99. We checked the dose and age inputs in the versioned AS file case by case against the clinical table.

IREP guidance identifies iodine-131 as an example for the electrons E>15 keV route. This route has a radiation effectiveness factor of 1.0 relative to reference higher-energy photons (Land et al. 2003), with independent microdosimetric support for similar biological effectiveness of internal iodine-131 and external gamma rays per absorbed dose (Sato et al. 2014). This addresses radiation quality, not the transfer of a thyroid risk model fitted to external photon exposures. The IREP thyroid model is based on atomic-bomb survivors and populations irradiated as children for medical indications (Ron et al. 1995; Kocher et al. 2008).

We retained the acute setting to preserve computational continuity with earlier IREP applications to post-Chernobyl papillary carcinoma from the same research programme (Zurnadzhy et al. 2022; Bogdanova et al. 2022). These studies document the provenance of the convention rather than independently validating acute coding for protracted internal iodine-131 exposure. In IREP, exposure rate affects the calculation through the dose and dose-rate effectiveness factor. Under the acute setting it approaches 1 at and above an uncertain reference dose distributed log-uniformly from 0.03 to 0.2 Gy; the High band's minimum dose was 188 mGy. A chronic calculation lowers the numeric AS scale. In an input-convention audit based only on the paired AS calculations and not on expression outcomes, chronic coding preserved the ordering of exposed cases, with no pairwise rank inversions. This audit establishes stability of the ordering, not the physical validity of acute coding or the invariance of fixed numerical cut-points across calculation scales.

Dose was entered as a point value, and dosimetric or risk-model uncertainty was not propagated into band membership. Dose-zero cases formed the separate reference group. Exposed cases were classified as Low (0<AS<33.3), Mid (33.3≤AS<66.6), and High (AS≥66.6); no case lay on a boundary. The cut-points were operational boundaries on the prespecified acute-calculation scale and were not interpreted as biological or aetiological probability thresholds. The composition of all 440 cases by driver, band, and pairing status is in Table S1.

## Quality control and analysis cohorts

The outlier procedure operated separately within each observed group-by-tissue submatrix after filterByExpr. It used unnormalized log counts per million with prior count 2. The most extreme sample on the first principal component was removed iteratively using the two-sided PC-OD rule at alpha=0.05 until no further sample was rejected (Nakayama et al. 2024). PC-OD was a QC and cohort-definition step, not a gene-level hypothesis test, and its finalized flags were held fixed in downstream inference. It flagged no tumor or normal sample in the RET fusion-positive stratum; its only flag among the four main target groups was one High-AS BRAF V600E-positive tumor. Thus, it did not alter the primary RET-tumor cohort.

We estimated relative tumor purity jointly for dose-zero and High-AS samples within each driver cohort using the maximum-one proportion estimator from contamDE-lm, implemented in house for the proportion step only and without its differential-expression model (Shen et al. 2016; Ji et al. 2020). Protein-coding, filterByExpr-reduced paired counts were MUREN-normalized before estimation. Joint estimation placed both groups on a common relative scale. We treated purity as a relative score within the jointly estimated cases, not as an absolute cellular fraction.

Case pairing used the GDC merged-aliquot expression assays. Each included case was required to have one unique merged tumor and one unique merged normal sample. The main cohort additionally required RET or BRAF classification, membership in the dose-zero or High-AS group, both tissues passing the outlier screen, and relative purity at least 0.6.

The REO construction set comprised the RET fusion-positive subset of the main cohort. The intermediate-band application set comprised RET fusion-positive PTCs from the Low-AS and Mid-AS bands, whether or not a matched normal was available. Construction-cohort outlier and purity exclusions were not applied to this set; corresponding metrics were reported as ancillary diagnostics and did not authorize exclusions.

## Analysis contrasts

The four High-AS versus dose-zero contrasts were RET-tumor (`R_Tumor` in the analysis files), RET-normal (`R_Normal`), BRAF-tumor (`B_Tumor`), and BRAF-normal (`B_Normal`). These labels refer to contrasts, not sample groups. In each contrast, group x was dose-zero and group y was High-AS. The Brunner–Munzel relative effect was

p = P(X<Y) + 0.5 P(X=Y).

Thus, p>0.5 indicated higher expression in High-AS cases. The signed effect used in Figure 2 was 2p−1. The RET-tumor Higher Criticism test was the single primary contrast-level test. The other three contrasts addressed separate secondary questions. Their p-values were not combined into a study-level decision, so no across-contrast adjustment was applied; all contrasts were reported, and no study-wide FDR or across-contrast family-wise error claim was made.

## Covariate disclosure and estimand

We intentionally incorporated age at exposure into AS as a modifier of radiation risk. Because age at exposure was undefined for the dose-zero group, we did not compare it between groups. The difference in age at surgery was summarized in each driver stratum with the Hodges–Lehmann median difference and the Brunner–Munzel relative effect P(Sporadic<High). Each estimate received a 95% percentile bootstrap interval from 9,999 within-group resamples using seed 19450809. No p-value was calculated; these summaries characterized group structure rather than testing confounding or demonstrating balance.

We reported sex by group because sex-chromosome genes can differ with group composition and annotated all selected genes for chromosome membership against GENCODE v36. RET fusion-partner counts were reported by band; the BRAF stratum contained BRAF V600E by construction. We fitted no common covariate-adjusted model: conditioning on the age inputs used to construct AS would change the marginal AS-band estimand, relative purity was expression-derived and used for eligibility, and fusion partner may be part of the biological pathway under study. Treating them collectively as nuisance terms would not, by itself, identify a radiation-specific effect.

The finalized analysis objects contained no sequencing-batch, processing-batch, centre, or collection-period field suitable for a group-balance calculation. This unavailable information is treated as an unassessed source of technical confounding. Available relative-purity, age, sex, and fusion-partner variables were reported descriptively without an additional hypothesis test.

## Normalization

For each contrast, we reduced protein-coding genes with filterByExpr and normalized their counts with an in-house three-iteration implementation of DEGES-MUREN (Kadota et al. 2012; Sun et al. 2013; Feng and Li 2021). In this two-group use, filterByExpr uses the group sizes to set the required number of expressed samples but does not use the identities of the samples meeting that expression rule; for fixed group sizes, the retained set is therefore invariant to relabelling. The MUREN step used its single-parameter pairwise least-trimmed-squares form. At each iteration, a permutation Brunner–Munzel screen produced Storey q-values with lambda=0.5. The larger of the q<0.10 set and the top 5% of genes by raw p-value was removed from scaling-factor estimation, implementing the floorPDEG guard. Tested-gene counts, screen pi0 estimates, Jaccard convergence, and scaling-factor ranges are in Table S2.

The DEGES exclusion screen uses the observed group labels, whereas the reported gene-level allocations hold the resulting normalized matrix fixed. Accordingly, “exact” refers to exhaustive enumeration of labels on the finalized matrix, not to end-to-end randomization of QC and normalization from raw counts. PC-OD is not part of the gene test and, in the primary RET cohort, removed no sample. The remaining scope qualification concerns the group-aware DEGES normalization and the exchangeability assumption for the retained samples.

## Gene-level differential expression

For each gene, we evaluated the studentized Brunner–Munzel statistic for every allocation of nx labels among n samples on the finalized normalized matrix. Exhaustive enumeration removed Monte Carlo error and the 1/(B+1) resolution floor. It did not, by itself, remove the exchangeability condition of the permutation null or turn upstream QC and normalization into part of the enumerated procedure. Ties were retained in the effect estimate and statistic.

Gene-level p-values were converted to Storey q-values with the plug-in estimate of pi0 at fixed lambda=0.5 and threshold q<0.10 (Storey 2002). The fixed value avoided adaptive tuning. This definition specifies how the ranked gene list was constructed; it is not presented as a proof of FDR control under arbitrary gene dependence. We therefore refer to “genes meeting the prespecified Storey q rule” and reserve the contrast-level existence statement for the Higher Criticism label-permutation result.

## Contrast-level omnibus inference

Higher Criticism was the prespecified primary omnibus statistic, with scan range alpha0=0.1 (Donoho and Jin 2004). Its p-value was calculated against the contrast's own label-shuffle distribution, so the analytic independence and sparsity assumptions of the original Higher Criticism null were not used. The RET-tumor Higher Criticism test was the single primary contrast-level test; the corresponding statistics in the other contrasts answered separate secondary questions. Count statistics at fixed cutoffs and the maximum statistic were retained as descriptive rows. The full rejection curve R(alpha) was retained to show how the displayed gene count varied across q thresholds.

For contrast-level diagnostics, 9,999 label shuffles were generated with seed 19860426 and stored. Applying the plug-in estimator to each shuffle gave a null reference distribution for pi0. These same shuffles supplied the omnibus and gene-set references.

Because we defined no binary contrast label, we report the omnibus p-values, per-gene q-values, and rejection curves as continuous evidence. We do not use “support” or “no support” wording to imply an undeclared cutoff.

## Gene-set level inference

We ranked genes by tie-averaged normal scores of the signed Brunner–Munzel statistic and evaluated the weighted running sum with gseaParam=1 only at tie-block boundaries. Automated tests verified equality to the standard GSEA statistic for tie-free input (Subramanian et al. 2005).

The null comprised the 9,999 saved label shuffles per contrast. For each set, a sign-conditional permutation p-value was calculated and adjusted by Benjamini–Hochberg within collection. pi0 was fixed at 1 for the set-level analysis. The four families were Hallmark; C2 canonical pathways restricted to Reactome, WikiPathways, KEGG MEDICUS, BioCarta, and PID; C5 GO Biological Process; and a radiation-curated subset of C2:CGP. The radiation curation rule was fixed before the reported set-level run. Sets outside 15–500 genes were excluded.

For the spike-in, we multiplied the 195 genes present from HALLMARK_ADIPOGENESIS by 1.15 in the nine High-AS BRAF-tumor samples. This was a single coherent-signal check, not a general power analysis.

For over-representation analysis, RET-tumor genes meeting Storey q<0.10 were split into higher, lower, and combined lists. One-sided hypergeometric tests used all genes tested in the RET-tumor contrast as the universe and BH adjustment within family-by-list. In the terminology of Goeman and Bühlmann (2007), the statistic is competitive, but its gene-sampling p-value does not refer to subject-level label randomization.

## Complete-null assessment of the set-level procedure

For each contrast, we shuffled labels B+R times with seeds derived from 19450809, using B=9,999 as a shared reference pool and R=100 held-out pseudo-observations. We normalized each pseudo-observation against the pool, calculated per-set permutation p-values, adjusted them within collection, and applied the q<0.10 threshold.

For each contrast-by-collection cell, the proportion of held-out global-null replicates yielding at least one discovery estimated P(at least one false discovery) under the complete null, where FDR and family-wise error coincide. Exact binomial intervals and mean discovery counts were reported in Table S7. The 100 indicators shared one null pool and were positively correlated; ordinary binomial intervals are therefore somewhat narrow.

This assessment informed selection of per-set permutation p-values with within-collection BH before the reported real-data run. An earlier pooled tail-ratio procedure yielded pooled complete-null rates of 0.140 and 0.221 for two variants, with a worst cell of 0.44, and we rejected it. The adopted procedure yielded 0.045 in its development assessment. In the final complete-null assessment, the pooled rate was 0.064, with one cell-specific rate of 0.18. These results characterize complete-null behaviour; they do not prove FDR control for every partial-alternative configuration.

## Between-stratum concordance

For each tissue, we intersected signed Brunner–Munzel profiles by gene and compared the RET and BRAF profiles using Spearman correlation. We shuffled labels independently within each contrast 9,999 times, using contrast-specific seeds derived from 19450809, and correlated the shuffled profiles to form the null reference interval. This interval represents the case in which neither stratum carries label-aligned structure. An observed correlation outside the interval does not distinguish exposure from a covariate shared across strata.

## REO panel construction and intermediate-band application

We recalculated TPM from the selected assay using exon-union lengths. The candidate pool comprised the 500 genes with the greatest absolute distance of the Brunner–Munzel effect from 0.5. We did not use the q<0.10 list to define the pool.

For a gene pair, r was the within-sample log2-TPM difference. The dead zone was |r|<log2(1.2). A pair qualified when its sign was stable in dose-zero samples, with at most one exception among samples outside the dead zone; the 10th percentile of |r| was at least log2(1.5); and the pair reversed in more than 50% but fewer than 100% of High samples. Qualified pairs were ranked by the shift in median r. Greedy selection prohibited gene reuse and excluded a pair whose r profile had Spearman correlation at least 0.75 with a retained pair. The target was 10 pairs.

The reversal score counted pairs outside the dead zone that opposed the dose-zero reference sign. The construction-derived classification rule defined scores greater than 2 as positive, placing the boundary above the maximum score observed among dose-zero construction samples. We did not change the panel or its classification rule after construction.

We defined binary and graded readouts for the panel. In the dose-zero and High-AS construction bands, the binary readout characterized separation between the groups used to construct the panel. Demonstrating this separation was a prerequisite for application elsewhere, but it was not an unbiased estimate of classification performance because the same cases informed pair selection and the dose-zero cases set the boundary. We then applied the classification rule unchanged to Low- and Mid-AS RET fusion-positive PTCs and treated the resulting classifications as descriptive application results.

We summarized graded scores across all four bands, with dose-zero and High-AS serving as construction anchors and the panel applied without refitting to Low-AS and Mid-AS. Within the two intermediate bands, the prespecified one-sided Monte Carlo Brunner–Munzel comparison of Mid-AS with Low-AS (seed 19860426) assessed whether scores tended to be higher in Mid-AS than in Low-AS. It did not test linearity, four-band monotonicity, or a dose–response form. Neither this comparison nor the four-band display constituted independent validation because the panel and score direction came from the construction bands.

## REO diagnostics

We applied the construction-cohort PC-OD screen to Low- and Mid-AS RET fusion-positive PTCs as a non-exclusion diagnostic, recording only the number flagged and retaining all intermediate-band cases. We estimated relative purity on a common scale by pooling paired RET fusion-positive cases across all bands in one contamDE-lm proportion run. Purity was available for 31 of the 36 intermediate-band tumors. The remaining two Low-AS and three Mid-AS cases had a tumor assay but no matched normal assay, which the pair-based purity estimator requires; they remained in the REO analysis.

We compared band and score conditional on purity using partial Spearman correlation with a one-sided permutation reference. We also compared Mid-AS and Low-AS within two strata split at median purity. These analyses were diagnostics and did not establish independence from purity.

## External gene-list comparison

We counted k of n source-list genes present among each contrast's tested genes and q<0.10 list. Abend et al. (2013) provided separate tumor and normal lists. Abend et al. (2012) provided a tumor-minus-normal paired-difference list, which had no matching estimand in the present contrasts. Dom et al. (2012) provided a normal-tissue list, and CLIP2 served as a single-gene reference. Because the lists were small and differed in tissue, platform, and contrast, no enrichment statistic was calculated.

## Software, seeds, and reproducibility

We built the container from the immutable Ubuntu noble-20260410 base and matching apt snapshot, compiled R 4.5.3 against reference BLAS, and installed R packages from the 2026-04-09 CRAN and Bioconductor 3.22 snapshots. The bit-level reproduction contract targeted linux/amd64.

In-house components included exact-enumeration and Monte Carlo Brunner–Munzel statistics, the Storey plug-in estimator, tie-block enrichment, DEGES-MUREN normalization, the contamDE-lm proportion score, the principal-component outlier procedure, and REO panel construction. External packages included edgeR, SummarizedExperiment, msigdbr, limma, GenomicDataCommons, rtracklayer, Rcpp, and MASS; full versions are in Table S5.

The canonical inference seed was 19860426. Diagnostics used documented seeds derived from 19450809. We independently executed the full pipeline on two x86-64 machines using 1,819 MD5-identical raw data files. Both runs passed 415 tests without failures and produced identical primary artifacts.
