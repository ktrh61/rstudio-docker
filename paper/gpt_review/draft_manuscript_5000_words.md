# Driver-conditioned transcriptomic differences across radiation-attributability bands in papillary thyroid carcinoma

> Revised editorial trial prepared on 2026-08-19 from paper/draft_manuscript.md and the canonical ledgers in paper/. It does not supersede those files. A narrowly scoped descriptive audit was added from the finalized analysis objects; it did not rerun QC, normalization, hypothesis testing, or multiple-testing adjustment and produced no alternative gene list. Before submission, add the exact GDC accession/data release, ethics statement, and declarations.

## Abstract

### Background

Radiation dose is associated with fusion rather than point-mutation drivers in papillary thyroid carcinoma after the Chernobyl accident. Because drivers strongly organize tumor expression, pooled analyses can confound exposure with driver composition. We tested expression differences by radiation-attributability band within driver strata.

### Methods

We analysed RNA-sequencing counts and clinical data from 440 REBC-THYR cases. Assigned Share (AS) from NIH IREP served only as an ordinal index. The main cohort comprised 63 paired cases in dose-zero or high-AS groups within RET-fusion and BRAF V600E strata. We used permutation-enumerated Brunner–Munzel statistics, Storey q-values, a pre-specified Higher Criticism omnibus statistic, and label-permutation gene-set tests. A 10-pair relative-expression-ordering panel trained on the RET extreme bands was evaluated without refitting in the intermediate bands.

### Results

In RET tumors, 1,765 of 15,621 genes met Storey q<0.10; omnibus p=0.0112. Their median absolute signed rank effect was 0.600 and median absolute display-only log2 fold change was 0.422. No gene met q<0.10 in RET normal tissue or BRAF tumors. BRAF normal tissue yielded one X-linked gene; omnibus p=0.0773. No gene set met q<0.10. The RET-tumor list was over-represented for proliferation, cell-cycle, and DNA-repair annotations under a gene-sampling reference. The independent REO Mid–Low comparison was uncertain (one-sided p=0.1127).

### Conclusions

An AS-band-associated expression difference was observed within RET-fusion tumors, but it does not establish a radiation-specific trace. Age and birth cohort, purity, fusion partner, sex, and unmeasured technical structure remain competing explanations. Normal-tissue and out-of-sample REO results were inconclusive.

## Introduction

Papillary thyroid carcinoma after the Chernobyl accident is a major human model of radiation carcinogenesis. Post-Chernobyl childhood tumors frequently carry RET fusions, with NCOA4-RET more common than in age-matched sporadic childhood tumors (Nikiforov et al. 1997; Rabes et al. 2000). Irradiation can induce RET rearrangements in thyroid cells (Ito et al. 1993; Caudill et al. 2005), participating loci are spatially proximate in thyroid nuclei (Nikiforova et al. 2000), and RET rearrangements are enriched among highly exposed atomic-bomb survivors (Hamatani et al. 2008). In the REBC-THYR cohort, increasing dose was associated with fusion rather than point-mutation drivers (Morton et al. 2021). Whether an expression difference remains after accounting for driver is unresolved.

This question requires driver conditioning. Tumor transcriptomes are strongly organized by driver, and driver prevalence varies with dose (Morton et al. 2021); a pooled exposed-versus-unexposed comparison therefore mixes exposure with driver composition. Previous studies reported dose-associated expression in tumor or normal tissue (Abend et al. 2012, 2013), CLIP2 as a tumor-tissue marker (Hess et al. 2011), and normal-tissue signatures in exposed cases (Dom et al. 2012; Ory et al. 2026). Most expression studies did not condition on driver, leaving open whether their signals persist within comparable driver backgrounds.

We hypothesized that comparing the extreme bands of a case-level radiation-attributability index might provide greater contrast than an unstratified dose model, without treating that index as an observed tumor aetiology. We used the NIH Interactive RadioEpidemiological Program Assigned Share (AS), derived from recorded thyroid dose and age, only to order exposed cases and define bands within driver strata.

Only the RET-fusion and BRAF V600E strata contained enough dose-zero and exposed cases for analysis. We examined four contrasts: high-AS versus dose-zero cases in tumor and normal tissue within each stratum. R_Tumor was the single primary contrast. The normal-tissue and BRAF contrasts addressed distinct secondary questions and were not combined into a study-wide significance claim. The interpretation framework was set before the finalized analysis produced the results reported here (Table 2). We also constructed a relative-expression-ordering (REO) panel in the RET extreme bands and evaluated it without refitting in the intermediate bands.

## Methods

### Data sources

Gene-level STAR counts for REBC-THYR were downloaded from the NCI Genomic Data Commons and verified against the download manifest. The clinical source was Data S1 of Morton et al. (2021), containing 440 cases. GDC files were mapped to cases and biospecimens through the GDC API and assembled into a matrix of 58,448 nonzero genes and 906 samples. Strandedness, gene-length derivation, file-selection rules, and the required accession and release identifiers are described in Supplementary Methods.

### Assigned Share

We calculated AS for exposed cases with the NIH IREP thyroid model, version 5.7.3, using each case's recorded thyroid dose and ages (Kocher et al. 2008). The input route followed IREP guidance for internal iodine-131 exposure, and the acute exposure-rate convention followed previous IREP applications to post-Chernobyl papillary carcinoma (Zurnadzhy et al. 2022; Bogdanova et al. 2022). The thyroid risk model was fitted to external photon exposures; applying it to internal iodine-131 is therefore a model-transfer assumption. Dose was entered as a point value, so dosimetric and model uncertainty did not propagate into band membership.

AS was used only as an ordinal design index. Dose-zero cases formed the reference group, labelled Sporadic in the analysis files. Exposed cases were classified as Low (0<AS<33.3), Mid (33.3≤AS<66.6), or High (AS≥66.6); no case lay on a boundary. Primary comparisons used High versus dose-zero group membership, not the numeric AS value, and AS was not interpreted as a calibrated probability that an individual tumor was radiogenic. The full input convention and the implications of alternative exposure-rate settings are in Supplementary Methods.

### Quality control and cohorts

Eligibility conditions were encoded as case-level flags before the reported inference. Samples were screened for outliers with a principal-component procedure before relative tumor purity was estimated. Purity was estimated from matched tumor-normal counts with the proportion component of contamDE-lm after MUREN normalization (Shen et al. 2016; Ji et al. 2020; Feng and Li 2021). Both exposure groups were estimated together within each driver cohort to place them on a common relative scale. Purity served as a within-cohort filter and diagnostic covariate, not as an absolute purity estimate.

The main cohort comprised RET- or BRAF-classified cases in the dose-zero or High group, with one paired tumor and normal sample, both tissues passing the outlier screen, and relative purity ≥0.6 (n=63). The RET subset was used to train the REO panel. RET tumors in the Low and Mid bands formed the evaluation set and were deliberately not excluded by the training outlier or purity rules (n=36). Their quality metrics were reported without exclusion authority.

### Analysis contrasts and interpretation

The four contrasts were R_Tumor, R_Normal, B_Tumor, and B_Normal. In each, group x was dose-zero and group y was High; the Brunner–Munzel relative effect was p=P(X<Y)+0.5P(X=Y), so values above 0.5 indicated higher expression in High. The R_Tumor Higher Criticism test was the single primary contrast-level test. R_Normal, B_Tumor, and B_Normal addressed separate secondary questions. Because their results were not combined to declare success of the study, no across-contrast adjustment was applied; all four results are reported, and no across-contrast family-wise or study-wide FDR claim is made.

The interpretation map is shown in Table 2. R_Tumor was the expected signal-bearing contrast. R_Normal addressed a possible normal-tissue signal; B_Normal supplied its cross-driver counterpart. B_Tumor was direction-agnostic and was not treated as a specificity control. The framework was set before the finalized analysis produced the results reported here.

### Candidate confounders

Age, sex, and RET fusion partner followed directly from the design and were reported rather than used to redefine eligibility. AS includes age at exposure, and the single 1986 exposure event induces age and birth-cohort structure across bands. We therefore report group distributions and Hodges–Lehmann differences in age at surgery with 95% bootstrap intervals, but do not treat them as tests of confounding. Sex composition was reported because sex-chromosome expression can differ mechanically between groups. RET fusion partners were reported because partner distributions have differed between radiation-associated and sporadic tumors (Nikiforov et al. 1997; Rabes et al. 2000).

These disclosures do not eliminate confounding. Residual tumor purity, age and birth cohort, sex, fusion partner, and unmeasured technical structure can align with the exposure label. The present estimates are therefore associations between AS-defined groups. No sequencing- or processing-batch field was available in the finalized analysis objects used for the additional audit; this absence is reported as a limitation rather than treated as evidence of balance.

### Normalization and gene-level inference

Protein-coding genes retained by filterByExpr were normalized separately for each contrast with an in-house implementation of three-iteration DEGES-MUREN (Kadota et al. 2012; Sun et al. 2013; Feng and Li 2021). Potential differentially expressed genes were removed from scaling-factor estimation with a permutation Brunner–Munzel screen at Storey q<0.10. Detailed filtering, the floorPDEG rule, convergence, and scaling factors are in Supplementary Methods and Table S2.

For each gene, we enumerated all C(n,nx) allocations of a studentized Brunner–Munzel statistic (two-sided) on the finalized normalized matrix. Enumeration was exhaustive and therefore introduced no Monte Carlo error. Its inferential interpretation remains conditional on exchangeability of the retained samples and does not make the preceding QC and normalization steps an end-to-end exact randomization procedure. Gene-level lists used the prespecified Storey q-value rule with a plug-in estimate of π0 at λ=0.5 and q<0.10 (Storey 2002). We describe the resulting genes as meeting this rule rather than claiming FDR control under arbitrary gene dependence.

Contrast-level evidence was summarized with a pre-specified Higher Criticism statistic scanning the lower 10% of the ordered p-values (Donoho and Jin 2004). Its p-value came from 9,999 label shuffles within the contrast (seed 19860426). Count and maximum statistics and the full rejection curve were retained as descriptive diagnostics. We report these p-values continuously and do not assign a binary label to each contrast.

### Descriptive magnitude and covariate audit

Without changing the existing R_Tumor gene list, we summarized the absolute signed rank effect, |2p−1|, and the absolute display-only mean log2 fold change among genes meeting the existing Storey rule. Mean log2 fold change was calculated from DEGES-normalized log2 counts per million with prior count 1 as the difference between the High and dose-zero group means; it was not used for gene selection. We also performed an unsupervised principal-component analysis of all 15,621 tested R_Tumor genes after gene-wise centering, without q-based feature selection or variance scaling, and overlaid group, sex, relative purity, age at surgery, and RET fusion partner. These calculations produced no new p-values and no second gene list.

### Gene-set inference

Genes were ranked by tie-averaged normal scores of the signed Brunner–Munzel statistic. An in-house weighted running-sum statistic was evaluated only at tie-block boundaries and equalled the standard GSEA statistic on tie-free input (Subramanian et al. 2005). The same 9,999 label shuffles used at gene level formed the set-level reference. Sign-conditional permutation p-values were adjusted by Benjamini–Hochberg within each collection, with q<0.10 and no cross-collection claim.

We tested Hallmark, selected C2 canonical-pathway subcollections (Reactome, WikiPathways, KEGG MEDICUS, BioCarta, and PID), C5 GO Biological Process, and a radiation-curated C2:CGP family; sets outside 15–500 genes were excluded. Complete-null operating characteristics were assessed with held-out label shuffles before the reported real-data run, and a one-set 1.15-fold spike-in provided a single-scenario positive-control check. These exercises did not establish FDR control under every partial-alternative configuration. Full procedures are in Supplementary Methods.

As a descriptive complement, the R_Tumor q<0.10 list was tested for one-sided hypergeometric over-representation in the same collections, using all tested R_Tumor genes as the universe and BH adjustment within family and direction. These gene-sampling p-values do not represent subject-level randomization; experiment-level interpretation rests on the label-permutation analysis (Goeman and Bühlmann 2007).

### Between-stratum concordance

For each tissue, signed gene-level Brunner–Munzel statistics from the RET and BRAF contrasts were compared over shared genes by Spearman correlation. Independently shuffled labels within each stratum supplied a reference interval. The normal-tissue comparison was hypothesis-bearing; the tumor comparison was a descriptive completion. A correlation outside the interval indicated shared label-aligned structure but could not distinguish exposure from a shared covariate.

### REO panel

The REO panel used within-sample gene ordering and therefore required no between-sample normalization (Geman et al. 2004; Wang et al. 2015; Guan et al. 2016). TPM was recalculated from the selected count assay and exon-union gene lengths. Candidate pairs were generated from the 500 RET-tumor genes with the largest Brunner–Munzel effect magnitude, rather than from the q<0.10 list. Pre-fixed stability, reversal, and redundancy rules selected 10 pairs without gene reuse. A sample's reversal score was the number of pairs, from 0 to 10, that exceeded a dead zone and reversed the dose-zero training sign. The training-derived classification boundary was a score greater than 2.

The fixed panel was applied without alteration to Low and Mid RET tumors. Because both extreme bands determined pair selection and the boundary, their apparent separation describes training fit and is not a performance estimate. The sole pre-specified out-of-sample comparison was a one-sided Monte Carlo Brunner–Munzel test of Mid over Low (seed 19860426). Outlier counts and purity analyses were ancillary diagnostics; they did not create exclusions or independent validation claims.

### External anchors and reproducibility

qRT-PCR-validated gene lists from Abend et al. (2012, 2013), Dom et al. (2012), and CLIP2 were cross-referenced against each contrast's q<0.10 list as descriptive membership counts, without an enrichment test.

The publication run used R 4.5.3 on Ubuntu 24.04 with reference BLAS/LAPACK 3.12.0, four workers, and fixed seeds in a date-pinned container. Two x86-64 machines used md5-identical raw inputs, passed 415 tests with no failures, and produced identical primary artifacts. Analysis code, versioned inputs, and the container recipe accompany the paper; package versions are listed in Table S5.

Generative AI assistants were used under author direction for code debugging and refactoring and for drafting and editing text. The authors specified and reviewed the scientific content, checked quantitative statements against versioned outputs and citations against the source ledger, made all scientific decisions, and accept responsibility for the manuscript.

## Results

### Cohort

Of 440 REBC-THYR cases, 63 met the main cohort definition: 12 R_Sporadic, 15 R_High, 27 B_Sporadic, and 9 B_High (Fig. 1; Table 1; Table S1). Most attrition resulted from driver classification and restriction to the extreme bands rather than pairing, outlier screening, or purity filtering. PC-OD flagged no RET tumor or normal sample and therefore did not alter the realized R_Tumor cohort; its only flag among the four main target groups was one B_High tumor. The REO evaluation set comprised 17 R_Low and 19 R_Mid tumors.

High cases were older at surgery than dose-zero cases by a Hodges–Lehmann estimate of 2.5 years in the RET stratum (95% CI −1.0 to 6.0) and 8.0 years in the BRAF stratum (3.0 to 12.0). Sex and fusion-partner distributions are shown in Table 1. These estimates describe group structure and are not tests that confounding is present or absent.

### Gene-level and contrast-level results

In R_Tumor, 1,765 of 15,621 tested genes met the prespecified Storey q<0.10 rule; 971 were higher and 794 lower in High. The Higher Criticism omnibus p-value was 0.0112 (Fig. 2; Table 3; Fig. S1). Among the 1,765 genes, the median absolute signed Brunner–Munzel effect was 0.600 (IQR 0.567–0.644), and the median absolute display-only mean log2 fold change was 0.422 (IQR 0.290–0.698). The corresponding median absolute log2 fold changes were 0.384 among genes higher in High and 0.573 among genes lower in High. The list contained 57 X-linked genes among 1,765, compared with 572 X-linked genes among 15,621 tested genes; no Y-linked gene was selected.

The unsupervised PCA used all tested R_Tumor genes. PC1 and PC2 explained 19.7% and 16.6% of total variance (Fig. S3). In the main-cohort analysis object, relative-purity medians were 0.783 in R_Sporadic and 0.822 in R_High; age-at-surgery medians were 20.5 and 23.0 years. Spearman correlations of PC1 with purity and age were −0.396 and −0.290, and those of PC2 were 0.226 and −0.042, respectively; these are descriptive overlays without p-values. Sex and RET-partner distributions are reported in Table 1 and shown on the same PCA coordinates.

R_Normal had no q<0.10 gene and HC p=0.3199. B_Tumor also had no q<0.10 gene and HC p=0.1815; it was not interpreted as a specificity control. B_Normal had one q<0.10 gene, the X-linked BHLHB9 (effect 0.967, q=0.013), with HC p=0.0773 and a descriptive maximum-statistic p=0.0125. Its π0 estimate was 0.727, compared with 0.955 and 0.943 in R_Normal and B_Tumor.

### Gene sets and list annotation

No set met q<0.10 in any of the 16 contrast-by-collection cells; the smallest adjusted value was 0.114 in B_Tumor's radiation-curated family (Table S6; Supplementary Data 1). In the complete-null assessment, at least one set was selected in 102 of 1,600 held-out contrast-by-collection replicates (0.064 pooled), while B_Normal/Hallmark had a cell-specific rate of 0.18 (95% CI 0.110–0.270). These quantities describe complete-null operation and are not a general FDR guarantee. In the single positive-control scenario, the planted set ranked first among 50 Hallmark sets (q=0.0101), and no other set met q<0.10 (Fig. S2; Table S7).

Under the separate gene-sampling reference, genes lower in R_High were over-represented in proliferation, cell-cycle, and DNA-repair annotations, including E2F targets (46/199), G2M checkpoint (41/198), and Reactome DNA repair (68/322). Radiation-curated flags reflected the same cell-cycle theme. No family met the over-representation threshold among genes higher in R_High (Table S3). These annotations were hypothesis-generating and did not change the permutation-based set-level result.

### Cross-stratum concordance

Normal-tissue contrast profiles had Spearman rho=0.376 across 15,459 genes, within the shuffled reference interval of −0.46 to 0.46 (two-sided p=0.1199). Because R_Normal itself showed no detectable structure, the data could not distinguish absent sharing from a shared signal below the detection limit.

The descriptive tumor comparison gave rho=0.459 across 15,560 genes, outside its shuffled interval of −0.39 to 0.39 (p=0.0197; Table S8). This was compatible with either shared label-aligned biology or shared covariate structure and was treated as hypothesis-generating.

### REO evaluation

Of 57,694 candidate pairs, 153 met the selection criteria and 10 formed the final panel (Table S4). In training, all 12 R_Sporadic cases scored below the classification boundary and 13 of 15 R_High cases scored above it; this separation was expected because both groups informed construction.

The fixed score had medians of 0, 1, 4, and 6 across Sporadic, Low, Mid, and High (Fig. 3). In the independent comparison, Mid exceeded Low with effect P(Low<Mid)=0.616 and one-sided p=0.1127. Thus the band ordering was descriptive and the Mid–Low difference remained uncertain.

No evaluation tumor was flagged by the mirrored outlier screen. The matched-pair purity estimator was available for 15 of 17 Low and 16 of 19 Mid cases; the other five cases had a tumor but no matched normal sample. Among the 31 tumors with purity scores, score and purity correlated at rho=0.538; band-score rho was 0.142 and was 0.146 after conditioning on purity (one-sided permutation p=0.2162). Band and purity correlated at rho=0.036. These ancillary diagnostics did not establish an AS-band effect independent of purity.

### External gene lists

Nineteen of 20 list-by-contrast cells had no overlap, including every tissue-matched comparison. The only overlap was cross-tissue: S100A10 from the Dom normal-tissue list occurred in R_Tumor with the direction opposite to the original report. The counts were descriptive and did not alter the primary interpretation.

## Discussion

The main observation was a broad expression difference between high-AS and dose-zero RET-fusion tumors: the primary omnibus p-value was 0.0112, and 1,765 genes met the analysis's prespecified Storey q<0.10 rule. The additional descriptive audit shows that the selected genes were not characterized only by small p-values: their median absolute rank effect was 0.600 and their median absolute display-only log2 fold change was 0.422. The result establishes label-aligned transcriptomic structure within the RET stratum under the reported analysis. It does not by itself identify the structure as a radiation-specific molecular trace. The comparison is observational, AS is a model-derived banding index rather than an observed aetiology, and measured or unmeasured covariates can align with its groups.

That distinction is important because several plausible covariates remain. Age and birth cohort are structurally linked to exposure status: the source cohort's dose-zero cases were born after the accident, whereas exposed cases were alive in 1986. Age at exposure also contributes to AS. The smaller age difference in the RET than BRAF stratum does not eliminate age as an explanation. Sex-chromosome membership did not dominate the R_Tumor list, and reported RET-partner counts showed no simple monotonic pattern, but these descriptive observations do not exclude confounding.

Tumor purity deserves particular caution, but the main-cohort imbalance was modest on the scale used for eligibility: medians were 0.822 in R_High and 0.783 in R_Sporadic. The unsupervised PCA likewise did not identify a single simple purity or age axis that accounted for both leading components, although the small cohort cannot exclude multivariable covariate structure. We did not fit a covariate-adjusted inferential model. Age contributes to the definition of AS; sex is a formal IREP input, although its numeric influence on the realized AS values was negligible; relative purity was an expression-derived eligibility variable; and RET fusion partner may be part of the biology linking exposure history to tumor expression. Entering these variables together as nuisance terms would therefore define a different conditional estimand rather than deconfound the prespecified marginal AS-band contrast. We instead disclose their observed distributions and limit the result to an AS-band-associated expression difference; it does not identify a radiation-specific effect independent of these variables. No sequencing- or processing-batch field was present in the finalized analysis objects, so technical confounding remains unassessed.

B_Normal illustrates why gene-level and contrast-level evidence should be separated. One X-linked gene met the gene-level threshold in the comparison with the largest sex imbalance, whereas the HC p-value was 0.0773. The BRAF groups also differed in age by an estimated 8 years. We therefore report the evidence without assigning a positive/negative label and regard confounding as the leading interpretation. The B_Normal/Hallmark global-null excess further limits set-level conclusions for that cell.

The permutation-based gene-set analysis selected no set. The over-representation analysis nevertheless showed that genes lower in R_High were concentrated in proliferation, cell-cycle, and DNA-repair annotations under a weaker gene-sampling reference. The two results concern different randomization schemes: the ORA describes the composition of an already selected list and ignores gene correlation, whereas the label-permutation analysis addresses subject labels. The annotation is therefore a biological lead, not independent confirmatory evidence.

The cross-stratum comparisons were also limited. Normal tissue provided no evidence that could distinguish absent sharing from a shared signal below the study's detection limit. The tumor profiles were more concordant than their shuffled reference, but this result was not the pre-specified normal-tissue test and is equally compatible with shared covariates, including purity. It should be treated as a hypothesis for a larger, covariate-balanced cohort.

The REO panel converted the RET extreme-band contrast into an individual score, but its extreme-band separation is training fit by construction. The only independent comparison was Mid versus Low; its effect was 0.616 and p=0.1127. The four-band medians are therefore a descriptive display, not evidence of classification performance or a dose response. Correlation between score and purity in the 31 paired evaluation cases adds uncertainty. Independent validation, preferably with fixed pair identities and complete purity and technical metadata, is required before the score can be considered a classifier or radiation-attributability marker.

This study has further limitations. The exposed BRAF group contained nine cases, restricting precision. The 1,765-gene list is the result of the specified Storey plug-in estimator; we do not claim that this proves FDR control under arbitrary transcriptomic dependence. PC-OD was a preceding QC step and removed no RET sample, so it did not affect the realized R_Tumor cohort. The exhaustive Brunner–Munzel enumeration is exact as an enumeration on the finalized matrix, but its inferential validity remains conditional on exchangeability; the group-aware DEGES normalization was not repeated within allocations, and the analysis should not be described as an end-to-end exact randomization of the raw-data pipeline. No wet-laboratory or independent cohort replication was available. The purity measure was relative within cohort, and IREP transferred an external-photon thyroid model to internal iodine-131 exposure. These limitations constrain both causal interpretation and generalizability.

In conclusion, high-AS and dose-zero RET-fusion tumors differed broadly in expression under the reported analysis, whereas normal-tissue, BRAF-tumor, gene-set, and independent REO evidence was limited or inconclusive. The result motivates further study of driver-conditioned radiation-associated expression, but it should presently be described as an AS-band association. Independent, covariate-informed validation is needed to determine whether it represents radiation-related biology.

## Main figure and table captions

**Figure 1 | Case flow.** Flow from 440 REBC-THYR cases to the 63-case main cohort and 36-case REO evaluation set, stratified by RET and BRAF driver.

**Figure 2 | Gene-level evidence by contrast.** Signed Brunner–Munzel effect, defined as 2p−1 from p=P(X<Y)+0.5P(X=Y), against −log10 permutation p. Genes meeting Storey q<0.10 are colored by direction.

**Figure 3 | REO reversal scores by AS band.** Scores from 0 to 10 for the fixed 10-pair panel. Dose-zero and High are training bands; Low and Mid are out-of-sample evaluation bands. The four-band pattern is descriptive.

**Supplementary Figure S3 | Descriptive R_Tumor PCA.** Unsupervised PCA of DEGES-normalized log2 CPM for all 15,621 tested genes, with identical sample coordinates displayed by AS group and sex, relative purity, age at surgery, and RET fusion partner. No q-based feature selection or hypothesis test was used.

**Table 1 | Case characteristics.** Group size, sex, age at surgery and exposure, relative purity, and RET fusion-partner composition, with age-difference estimates and confidence intervals.

**Table 2 | Interpretation map.** Standing and planned interpretation of the four contrasts, set before the finalized analysis produced the results reported here.

**Table 3 | Gene-level results.** Tested genes, π0, q<0.10 counts and directions, minimum permutation p, and Higher Criticism p by contrast.

## References

Use the complete, verified bibliography represented in paper/citation_ledger.md and paper/draft_manuscript.md, formatted to the selected journal. The current working reference list must be converted from ledger notes to a conventional reference section before submission.
