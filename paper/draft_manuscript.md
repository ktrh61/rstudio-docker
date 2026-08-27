# 原稿統合ドラフト(draft_manuscript — 論文へ収束する単一ファイル)

状態: draft(2026-08-14 統合。研究者指示: Methods/Results の対構造の検査可能化と同期機構の廃止のため、旧4ドラフト — draft_methods / draft_results / draft_methods_results_d6 / draft_ora_annotation — を本ファイルへ吸収し削除。内容の変更なし、移設+Discussion 骨組み追加のみ。履歴は git)。規約: 英語本文が正本。インライン【訳】は C2 ポートで段階的に廃止(対訳は内容確定後に共著者向け別ファイルとして再生成 — 2026-08-21 研究者方針)。Supplementary Methods/Results は paper/supplementary_material.md へ分離(2026-08-21、BJC 投稿構造)。数値は N-ID タグ、主張段落は C-ID タグ(検査対象は英語本文)。図表番号は 2026-08-16 に本置き済み(主: 図1–3・表1–3 / 補足: 図S1–S2・表S1–S8・データ1。BJC の Article 規定 = 主6品目以内に整合)。各節の由来・決定は末尾の「起草メモ(統合)」。

---

<!-- タイトルページ素材(BJC: title page 別建て・語数枠外) -->
**Title:** Driver-conditioned transcriptomic differences across radiation-attributability bands in papillary thyroid carcinoma

## Abstract

### Background

Radiation dose is associated with fusion rather than point-mutation drivers in post-Chernobyl papillary thyroid carcinoma, while driver composition structures expression. We asked whether expression differed between high estimated radiation-attribution cases and dose-zero cases within driver-defined groups.

### Methods

We calculated NIH IREP Assigned Share (AS) for ordinal banding only. We stratified by driver (RET fusion or BRAF V600E) and compared expression between High-AS and dose-zero cases in tumor and matched normal tissue; the RET fusion-positive tumor comparison was primary. We used permutation Brunner–Munzel tests, Higher Criticism, label-permutation gene-set inference, and a relative-expression-ordering (REO) panel constructed from extreme bands, applied unchanged to intermediate bands.

### Results

In RET fusion-positive PTCs, 1,765 of 15,621 genes met Storey q<0.10 <!-- N-16, N-15 -->, with omnibus p=0.0112 <!-- N-20 -->; other contrasts lacked concordant gene- and contrast-level evidence. No gene set met q<0.10 under label permutation, whereas over-representation analysis highlighted cell-cycle and DNA-repair annotations among genes expressed less highly in High-AS cases. REO score medians showed an ordered descriptive profile across bands (dose-zero/Low-AS/Mid-AS/High-AS, 0/1/4/6 <!-- N-41 -->), with the extremes as construction anchors.

### Conclusions

High-AS and dose-zero RET fusion-positive PTCs differed transcriptionally, compatible with the hypothesized within-driver AS-band association. This does not establish a radiation-specific signature, and the REO profile requires independent validation.

## Introduction

<!-- 2026-08-21 C2 ポート: gpt_review 試案 Intro を移植(C1-e 仮説スロット批准・凍結文言は新正準へ統合)。予告列の同期は C4。インライン対訳は廃止(f 保留方針) -->

Papillary thyroid carcinoma after the Chernobyl accident remains a major human model of radiation carcinogenesis. Its molecular hallmark was recognized early: childhood tumors arising after the accident frequently carry RET fusions, and their fusion-partner distribution, with NCOA4-RET predominating, differs from that of age-matched sporadic childhood tumors (Nikiforov et al. 1997; Rabes et al. 2000). The historical labels RET/PTC1 and RET/PTC3 correspond to CCDC6-RET and NCOA4-RET, respectively; we use the gene-partner forms throughout. Experimental irradiation can induce RET rearrangements in thyroid cells (Ito et al. 1993; Caudill et al. 2005), the participating loci lie close together in thyroid-cell nuclei (Nikiforova et al. 2000), and RET rearrangements have also been enriched among highly exposed atomic-bomb survivors (Hamatani et al. 2008). At cohort scale, the REBC-THYR study showed that increasing thyroid dose was associated with fusion rather than point-mutation drivers (Morton et al. 2021). Whether exposure is associated with tumor expression within a driver stratum remains unresolved.

Two features make that question difficult. First, tumor transcriptomes are strongly organized by driver, while driver prevalence itself varies with radiation dose (Morton et al. 2021). A pooled exposed-versus-unexposed analysis can therefore reproduce differences between fusion- and point-mutation-driven tumors rather than an exposure-associated difference within a comparable molecular background. Second, previous expression findings span different biological compartments and designs: dose-associated expression has been reported in tumor and normal tissue, including paired tumor-minus-normal analyses (Abend et al. 2012, 2013); CLIP2 was proposed as a tumor-tissue radiation marker (Hess et al. 2011); and exposed-case normal-tissue signatures have been reported without driver stratification (Dom et al. 2012; Ory et al. 2026). The CLIP2 study matched RET-rearrangement status as far as possible but focused on a copy-number marker and CLIP2 rather than genome-wide expression. Whether earlier expression signals persist within driver strata is therefore unclear.

The source cohort's continuous dose-expression analyses did not identify a corresponding expression signal (Morton et al. 2021), motivating a different but related question. We hypothesized that comparing extreme bands of a case-level radiation-attributability index might reveal structure diluted by an unstratified continuous-dose model. This reflects a mixture model in which dose changes the proportion of tumors arising through a radiation-associated initiating event without requiring expression in every tumor to vary monotonically with dose. The index was the Assigned Share (AS) from the NIH Interactive RadioEpidemiological Program (Kocher et al. 2008), calculated from recorded thyroid dose and age under a specified input convention. We used AS only to order exposed cases and define bands within driver strata. The bands were intended to enrich contrasting groups; they did not identify individual tumors as radiogenic or sporadic, and AS was not treated as an observed etiology or calibrated individual probability.

Only the RET fusion-positive and BRAF V600E-positive strata contained enough dose-zero and exposed cases for driver-stratified analysis. Within each stratum, we compared the High-AS band with dose-zero cases separately in tumor tissue and matched normal tissue; these four comparisons were prespecified. The RET fusion-positive tumor comparison was primary, reflecting both the prior expectation of radiation-related tumor initiation in this stratum and the cohort-level association of increasing dose with fusion rather than point-mutation drivers. BRAF therefore addressed a distinct within-driver question, not a negative or specificity control. The three secondary comparisons asked whether a corresponding difference occurred in normal tissue, whether expression differed in BRAF V600E-positive PTCs, and whether normal-tissue patterns were concordant across driver strata. Their results were not combined into a study-wide significance criterion. Table 1 gives the interpretation assigned to each possible pattern before the finalized analysis produced the results reported here. Separately, we constructed a relative-expression-ordering (REO) panel from the two RET extreme bands and applied it unchanged to Low-AS and Mid-AS RET fusion-positive PTCs.

## Methods

<!-- 2026-08-21 C2 ポート: gpt_review 試案 Methods を移植(C1-a AI開示復元・C1-d 目盛り要約復帰・chronic監査文は不採用)。写像検査・★アンカー同期は C4 -->
### Data sources

Gene-level STAR count files for REBC-THYR were downloaded from the National Cancer Institute (NCI) Genomic Data Commons (GDC) and verified against the download manifest. The clinical source was Data S1 of Morton et al. (2021), containing 440 cases <!-- N-11 -->. GDC files were mapped to cases and biospecimens through the GDC API and assembled into a matrix of 58,448 nonzero genes and 906 samples <!-- N-14 -->. Strandedness, gene-length derivation, file-selection, and manifest-verification procedures are in Supplementary Methods.

### Assigned Share

For each exposed case, we calculated the AS associated with the expected excess relative risk using the NIH IREP thyroid model, version 5.7.3, and the recorded thyroid dose and ages (Kocher et al. 2008). We used the IREP input route for internal iodine-131 exposure and a prespecified convention in which exposure rate was entered as acute. This setting defined the numerical AS scale used for band construction; it was not intended to reconstruct the temporal pattern of internal iodine-131 dose delivery. The thyroid risk model was fitted to external photon exposures; applying it to internal iodine-131 is therefore a model-transfer assumption. The exposure-rate setting enters the calculation only through the dose and dose-rate effectiveness factor, and the High-AS band's minimum dose (188 mGy <!-- N-84 -->) sits at the top of the model's uncertain reference-dose range <!-- N-83 -->, so a chronic calculation would compress the AS scale over the High band essentially by one common factor and could not materially reorder cases (Supplementary Methods). Dose was entered as a point value, so dosimetric and model uncertainty did not propagate into band membership.

AS was used only as a convention-dependent ordinal design index. Dose-zero cases formed the reference group, labeled Sporadic in the analysis files. Exposed cases were classified as Low (0<AS<33.3), Mid (33.3≤AS<66.6), or High (AS≥66.6); no case lay on a boundary <!-- N-69 -->. These were operational boundaries on the specified calculation scale. The four extreme-band contrasts used High-AS versus dose-zero group membership, not the numeric AS value, and AS was not interpreted as a calibrated probability that an individual tumor was radiogenic. The full input convention and its provenance are in Supplementary Methods.

### Quality control and cohorts

We represented the eligibility criteria as case-level flags and fixed them before the inferential analyses reported here. Principal-component-scores-based outlier detection (PC-OD) was applied separately within each observed group-by-tissue matrix after filtering genes with edgeR's filterByExpr, using unnormalized log counts per million with prior count 2 and iterative two-sided screening at α=0.05 (Chen et al. 2016; Chen et al. 2025; Nakayama et al. 2024). PC-OD preceded estimation of relative tumor purity so that an anomalous sample could not affect the purity fit. We estimated relative tumor purity jointly for dose-zero and High-AS cases within each driver cohort from matched tumor-normal counts, using the proportion component of contamDE-lm after MUREN normalization (Shen et al. 2016; Ji et al. 2020; Feng and Li 2021), thereby placing both groups on a common relative scale. Purity served as a within-cohort eligibility filter and diagnostic variable, not as an absolute cellular fraction.

The main cohort comprised RET fusion-positive or BRAF V600E-positive cases in the dose-zero or High-AS group, with one matched tumor–normal sample pair, both tissues passing the outlier screen, and relative purity ≥0.6 <!-- N-70 --> (n=63 <!-- N-08 -->). The RET fusion-positive subset supplied the REO construction bands. RET fusion-positive PTCs in the Low-AS and Mid-AS bands formed the intermediate-band application set and were deliberately not excluded by the construction-cohort outlier or purity rules (n=36 <!-- N-10 -->). Their quality metrics were reported without exclusion authority.

### Analysis contrasts and interpretation

We defined four contrasts comparing High-AS with dose-zero cases in tumor and matched normal tissue within each driver stratum. Their analysis-file labels were `R_Tumor` for tumor tissue from RET fusion-positive PTCs, `R_Normal` for their matched normal tissue, `B_Tumor` for tumor tissue from BRAF V600E-positive PTCs, and `B_Normal` for their matched normal tissue. These are contrast labels, not sample-group labels. Within each contrast, X denoted the dose-zero group and Y the High-AS group. The Brunner–Munzel relative effect (probabilistic index; Brunner and Munzel 2000), denoted θ, was θ=Pr(X<Y)+0.5Pr(X=Y); values of θ>0.5 indicated that expression tended to be higher in High-AS cases.

Higher Criticism (HC; Donoho and Jin 2004) was the single primary contrast-level test for RET-tumor. RET-normal assessed a corresponding normal-tissue difference, BRAF-normal supplied the other driver stratum for the prespecified normal-tissue comparison, and BRAF-tumor completed the driver-stratified tumor analysis. BRAF-tumor was direction-agnostic and was not treated as a specificity control. Because the secondary contrasts answered distinct questions and were not combined to declare overall study success, no across-contrast adjustment was applied; all four results are reported, and no across-contrast family-wise error rate or study-wide false discovery rate (FDR) is claimed.

Table 1 records the interpretation map assigned before the finalized analysis produced the results reported here. It designated RET-tumor as the expected signal-bearing contrast, specified the normal-tissue question, and required a confounding-first reading of potential BRAF-normal findings. Results outside this map are labeled hypothesis-generating.

### Covariate disclosure and estimand

For exposed cases, we incorporated age at exposure into AS as a radiation-risk modifier and did not additionally treat it as a nuisance covariate in the expression analysis. Because dose-zero cases were born after the accident whereas exposed cases were alive in 1986, the groups are structurally separated by birth cohort. We report group distributions and Hodges–Lehmann age-at-surgery differences with 95% bootstrap intervals <!-- N-63 --> to characterize this structure, not to test confounding or demonstrate balance. We also report sex composition because sex-chromosome expression can differ mechanically between groups, and RET fusion partners because partner distributions have differed between radiation-associated and sporadic tumors (Nikiforov et al. 1997; Rabes et al. 2000).

We fitted no common covariate-adjusted model because these variables lack a common nuisance interpretation. Conditioning on the age inputs to AS would replace the prespecified marginal AS-band contrast with an age-conditional estimand; relative purity was expression-derived and used for eligibility; and RET fusion partner may lie on the biological pathway connecting exposure history, tumor initiation, and expression. Their joint adjustment would answer a different question and would not identify a radiation-specific effect.

Relative tumor purity, age and birth cohort, sex, fusion partner, and unmeasured technical structure may also align with AS-band labels; this limits radiation-specific attribution, and the estimates describe associations between AS-defined groups. No sequencing- or processing-batch field was available in the analysis objects (Limitations).

### Normalization and gene-level inference

We retained protein-coding genes using filterByExpr and normalized their counts separately for each contrast using an in-house implementation of three-iteration DEGES-MUREN (Kadota et al. 2012; Sun et al. 2013; Feng and Li 2021). For scaling-factor estimation, we replaced the original model-based differential-expression screen with a permutation Brunner–Munzel screen at Storey q<0.10 (Storey 2002) and excluded the genes it flagged as potentially differentially expressed <!-- N-71 -->. Supplementary Methods and Table S2 report the detailed filtering, floorPDEG rule, convergence, and scaling factors <!-- N-15 -->.

For each gene, we enumerated all C(n,nx) allocations of a studentized Brunner–Munzel statistic (two-sided) on the finalized normalized matrix. Enumeration was exhaustive and therefore introduced no Monte Carlo error. Because the DEGES exclusion screen used the observed group labels, gene-level enumeration was conditional on the resulting normalized matrix and exchangeability of the retained samples. Gene-level lists used the prespecified Storey q-value rule with a plug-in estimate of the null proportion π0 at λ=0.5 and q<0.10 <!-- N-05 -->. We describe the resulting genes as meeting this rule rather than claiming FDR control under arbitrary gene dependence.

Contrast-level evidence was summarized with a prespecified Higher Criticism statistic scanning the lower 10% of the ordered p-values. We chose it to aggregate moderately small p-values rather than depend on one extreme gene or a fixed gene-count threshold. Its p-value came from 9,999 label shuffles within the contrast (seed 19860426) <!-- N-04 -->. Count and maximum statistics and the full rejection curve were retained as descriptive diagnostics. We report these p-values continuously and do not assign a binary label to each contrast. <!-- C-16 -->

### Gene-set inference

We ranked genes by tie-averaged normal scores of the signed Brunner–Munzel statistic and evaluated an in-house weighted running-sum statistic only at tie-block boundaries. On tie-free input, this statistic equaled the standard gene set enrichment analysis (GSEA) statistic (Subramanian et al. 2005) <!-- N-72 -->. We used the 9,999 saved contrast-specific label shuffles <!-- N-07 --> to form the set-level reference distribution and adjusted sign-conditional permutation p-values within each collection using the Benjamini–Hochberg procedure (BH; Benjamini and Hochberg 1995), applying q<0.10 without making a cross-collection claim.

We tested Hallmark, selected C2 canonical-pathway subcollections (Reactome, WikiPathways, KEGG MEDICUS, BioCarta, and PID), C5 GO Biological Process, and a radiation-curated C2:CGP family <!-- N-29, N-55 -->; sets outside 15–500 genes were excluded <!-- N-72 -->. <!-- C-06 -->
An earlier held-out complete-null assessment was used to select the set-level procedure before the reported real-data run, and the selected procedure's final operating characteristics were assessed with a 9,999-shuffle reference pool <!-- N-06 -->. A one-set 1.15-fold spike-in provided a single-scenario positive-control check <!-- N-30 -->. Full procedures are in Supplementary Methods.

<!-- C-14 -->
As a descriptive complement, the primary RET-tumor q<0.10 list was tested for one-sided hypergeometric over-representation in the same collections, using all genes tested in that contrast as the universe and Benjamini–Hochberg adjustment within each family-by-list combination <!-- N-59 -->. These gene-sampling p-values do not represent subject-level randomization; experiment-level interpretation rests on the label-permutation analysis (Goeman and Bühlmann 2007).

### Between-stratum concordance

For each tissue, we compared signed gene-level Brunner–Munzel statistics from the RET and BRAF contrasts across shared genes using Spearman correlation. We generated a reference interval by independently shuffling labels within each stratum. The normal-tissue comparison was hypothesis-bearing, whereas the tumor comparison was a descriptive completion. A correlation outside this interval indicated shared label-aligned structure but could not distinguish exposure from a shared covariate.

### REO panel

We based the REO panel on within-sample gene ordering, which does not require between-sample normalization (Geman et al. 2004; Wang et al. 2015; Guan et al. 2016). We recalculated transcripts per million (TPM) from the selected count assay using exon-union gene lengths. We generated candidate pairs from the 500 RET-tumor genes with the largest Brunner–Munzel effect magnitudes, rather than from the q<0.10 list, and applied prespecified stability, reversal, and redundancy rules to select 10 pairs without gene reuse <!-- N-37 -->. For each sample, the reversal score counted, from 0 to 10, the pairs whose expression difference lay outside the dead zone and had the sign opposite to the dose-zero construction reference. We set the classification cutoff at the maximum reversal score among dose-zero construction samples and classified scores above that cutoff as positive.

<!-- C-08 の手法側 -->
We defined binary and graded readouts for the panel. In the dose-zero and High-AS construction bands, the binary readout characterized separation between the groups used to construct the panel. Demonstrating this separation was a prerequisite for application elsewhere, but it was not an unbiased estimate of classification performance because the same cases informed pair selection and the dose-zero cases set the cutoff. We then applied the classification rule unchanged to Low-AS and Mid-AS RET fusion-positive PTCs and treated the resulting classifications as descriptive application results.

We summarized graded scores across all four bands, with dose-zero and High-AS serving as construction anchors and the panel applied without refitting to Low-AS and Mid-AS. Within the two intermediate bands, the prespecified one-sided Monte Carlo Brunner–Munzel comparison of Mid-AS with Low-AS (seed 19860426) <!-- N-75 --> assessed whether scores tended to be higher in Mid-AS than in Low-AS. It did not test linearity, four-band monotonicity, or a dose–response form. Neither this comparison nor the four-band display constituted independent validation because the panel and score direction came from the construction bands. <!-- C-09 --> We also applied PC-OD separately within the Low-AS and Mid-AS groups as a non-exclusion diagnostic. For intermediate-band tumors with matched normal tissue, we estimated relative purity on a common scale across all four RET groups and summarized the band–score association after rank-scale adjustment for purity using partial Spearman correlation and a descriptive one-sided permutation reference <!-- C-12 -->. We did not use either diagnostic to exclude cases or modify the panel.

### External gene-list comparison

We cross-referenced two classes of published radiation-associated gene lists against each contrast's q<0.10 list: qRT-PCR-validated lists from Abend et al. (2012, 2013) and Dom et al. (2012), with CLIP2 as a single-gene reference, and the shared-tissue, normal-tissue, and tumor-tissue multivariate transcriptomic signatures from Ory et al. (2026). <!-- C-13 --> We kept the two classes distinct and reported descriptive membership counts without performing an enrichment test <!-- N-53 -->.

### Reproducibility and AI use

We independently executed the full publication pipeline on two x86-64 machines using MD5-identical raw inputs <!-- N-52 --> and a date-pinned container with R 4.5.3 (R Core Team 2026) on Ubuntu 24.04 <!-- N-02 -->, reference BLAS/LAPACK 3.12.0 <!-- N-03 -->, four workers <!-- N-04 -->, and fixed seeds. Both runs passed 415 tests with no failures <!-- N-51 --> and produced identical primary artifacts. The analysis code, versioned inputs, and container recipe accompany the paper; Table S5 lists package versions.

The authors designed the study, specified the analysis and interpretation framework, and made the final scientific and editorial decisions. Under author supervision, Claude Code (Anthropic; Claude Fable 5) and Codex CLI (OpenAI; GPT-5.6 Sol) assisted with code development (principally debugging and refactoring), consistency checks, and drafting and editing manuscript text. The authors reviewed and approved all AI-assisted material retained in the reported analysis and manuscript and take full responsibility for the work; the scope and verification procedures are described in Supplementary Methods.

## Results

<!-- 2026-08-21 C2 ポート: gpt_review 試案 Results を移植(タグ移植込み。N-87/N-88 は新設 draft 行 — C3 で照合)。写像検査同期は C4 -->
### Cohort

<!-- C-10 --> Of 440 <!-- N-11 --> REBC-THYR cases, 63 <!-- N-08 --> met the main cohort definition (Fig. 1; Table 2; Table S1). The RET fusion-positive stratum contained 12 dose-zero cases (`R_Sporadic`) and 15 High-AS cases (`R_High`); the BRAF V600E-positive stratum contained 27 dose-zero cases (`B_Sporadic`) and nine High-AS cases (`B_High`) <!-- N-09 -->. Each contributed one tumor and one matched normal sample. Females comprised 10 of 12 dose-zero and 11 of 15 High-AS RET fusion-positive cases, compared with 23 of 27 and four of nine BRAF V600E-positive cases <!-- N-12, N-13 -->. RET fusion partners were CCDC6-RET in six and seven cases, NCOA4-RET in two and four, and other partners in four and four, respectively <!-- N-12 -->. Most attrition resulted from driver classification and restriction to the extreme bands rather than pairing, outlier screening, or purity filtering. PC-OD flagged no tumor or normal sample in the RET fusion-positive stratum and therefore did not alter the primary RET fusion-positive cohort. Its only flag among the four target groups was one High-AS BRAF V600E-positive tumor <!-- N-87 -->. The separate REO intermediate-band application set comprised 17 Low-AS and 19 Mid-AS RET fusion-positive PTCs <!-- N-10 -->.

High-AS cases were older at surgery than dose-zero cases by a Hodges–Lehmann estimate of 2.5 years in the RET stratum (95% confidence interval [CI], −1.0 to 6.0) and 8.0 years in the BRAF stratum (3.0 to 12.0) <!-- N-64, N-65, C-15 -->. Sex, fusion-partner, and relative purity distributions are summarized in Table 2 and its footnotes. These descriptive estimates neither establish nor exclude confounding.

### Gene-level and contrast-level results

In the primary RET-tumor contrast, 1,765 of 15,621 tested genes met the prespecified Storey q<0.10 rule <!-- N-16, N-15 -->; 971 were expressed more highly and 794 less highly in the High-AS group <!-- N-17, C-02 -->. The prespecified Higher Criticism omnibus p-value was 0.0112 <!-- N-20 -->, so the contrast-level evidence of label-aligned expression structure did not rely only on the size of the gene list (Fig. 2; Table 3; Fig. S1; Supplementary Data 2).

The three secondary contrasts gave different forms of limited evidence. RET-normal had no gene at q<0.10 and HC p=0.3199 <!-- N-16, N-20, C-03 -->. In BRAF-tumor (direction-agnostic, and not designed as a negative control; Methods), no gene met q<0.10 and the HC p-value was 0.1815 <!-- N-16, N-20 -->. These results established neither equivalence between the BRAF bands nor specificity of the RET-tumor finding <!-- C-04 -->. BRAF-normal had one gene at q<0.10, the X-linked BHLHB9 (relative effect 0.967, q=0.013) <!-- N-22 -->, while its contrast-level HC p-value was 0.0773 <!-- N-20 -->. A descriptive maximum-statistic p-value was 0.0125 <!-- N-23 -->. Its π0 estimate was 0.727, compared with 0.955 in RET-normal and 0.943 in BRAF-tumor <!-- N-19 -->. We report this divergence between gene-level and contrast-level evidence without assigning the contrast a binary label <!-- C-16 -->. Chromosome-level annotation of the reported gene lists is provided in Supplementary Results.

### Gene sets and list annotation

No set met q<0.10 in any of the 16 contrast-by-collection cells <!-- N-27, C-05 -->; the smallest adjusted value was 0.114 in the BRAF-tumor radiation-curated family <!-- N-28 --> (Table S6; Supplementary Data 1). In the complete-null assessment, at least one set was selected in 102 of 1,600 held-out contrast-by-collection replicates (descriptive pooled rate, 0.064) <!-- N-24, C-06 -->, while BRAF-normal/Hallmark had a cell-specific rate of 0.18 (descriptive 95% interval, 0.110–0.269 <!-- N-25 -->; Fig. S2; Table S7). These complete-null rates do not guarantee FDR control under partial alternatives. In the single positive-control scenario, the planted set ranked first among 50 Hallmark sets (q=0.0101) <!-- N-31 -->, and no other set met q<0.10 <!-- N-32 --> (Supplementary Results).

Under the separate gene-sampling reference, genes expressed less highly in High-AS RET fusion-positive PTCs were over-represented in proliferation, cell-cycle, and DNA-repair annotations, including E2F targets (46/199), G2M checkpoint (41/198) <!-- N-61 -->, and Reactome DNA repair (68/322) <!-- N-62 -->. Radiation-curated flags reflected the same cell-cycle theme rather than an independent theme <!-- N-60 -->. Among genes expressed more highly in High-AS RET fusion-positive PTCs, no set in any family met the over-representation threshold (Table S3) <!-- N-59, C-14 -->. These annotations were hypothesis-generating and did not change the null result of the label-permutation gene-set analysis.

### Between-stratum concordance

The RET-normal and BRAF-normal contrast profiles had Spearman rho=0.376 across 15,459 shared genes, within the shuffled reference interval of −0.46 to 0.46 (two-sided p=0.1199) <!-- N-34, C-07 -->. Because the RET-normal contrast itself showed no detectable structure, this prespecified comparison could not distinguish absence of a shared normal-tissue pattern from a shared pattern below the study's detection limit.

The corresponding tumor-profile comparison, included as a descriptive completion rather than a prespecified hypothesis test, gave rho=0.459 across 15,560 genes, outside its shuffled interval of −0.39 to 0.39 (p=0.0197; Table S8) <!-- N-33, C-17 -->. This pattern was compatible with either expression structure shared across driver backgrounds or covariates aligned with both exposure labels; it was therefore treated as hypothesis-generating.

### REO construction and intermediate-band application

<!-- C-08 --> Of 57,694 candidate pairs, 153 met the selection criteria <!-- N-35 --> and 10 formed the final panel <!-- N-37 --> (Table S4; Supplementary Results). The maximum reversal score among dose-zero construction cases was 2 <!-- N-38 -->. At the resulting cutoff, all 12 dose-zero RET fusion-positive cases scored at or below it and 13 of 15 High-AS RET fusion-positive cases scored above it <!-- N-38 -->. This separation was expected because both extreme bands informed pair selection and the dose-zero cases set the cutoff; it was not treated as an estimate of classification performance.

Applied unchanged, the construction-derived cutoff classified the Low-AS cases as nine negative and eight positive and the Mid-AS cases as eight negative and 11 positive <!-- N-42 -->. Graded-score medians were 0, 1, 4, and 6 across the dose-zero, Low-AS, Mid-AS, and High-AS bands (Fig. 3) <!-- N-41 -->. The prespecified Mid–Low comparison gave Pr(Low<Mid)=0.616 and one-sided p=0.1127 <!-- N-40 -->. This test left the directional difference between those bands uncertain; it did not test the shape of the four-band profile.

No intermediate-band tumor was flagged by the ancillary PC-OD screen <!-- N-43, C-09 -->. Relative purity could be estimated for 15 of 17 Low-AS and 16 of 19 Mid-AS cases <!-- N-44, N-10 -->; the other five cases had a tumor but no matched normal sample. Among the 31 tumors with purity scores, score and purity correlated at rho=0.538 <!-- N-45 -->; band–score rho was 0.142, and the partial Spearman coefficient after adjustment for purity rank was 0.146 (descriptive one-sided permutation p=0.2162) <!-- N-46, N-48, C-12 -->. Band and purity correlated at rho=0.036 <!-- N-48 -->. These ancillary diagnostics did not establish an AS-band effect independent of purity (Supplementary Results).

### External gene lists

Among the 20 validated-anchor-by-contrast cells, 19 had no overlap <!-- N-53 -->. The sole overlap was cross-tissue: S100A10 from the Dom normal-tissue list occurred in the RET-tumor list with the direction opposite to the original report <!-- N-53, N-54, C-13 -->. Among the multivariate Ory signatures, five of 39 tested shared-tissue genes and three of 46 tested tumor-tissue genes occurred in the RET-tumor list <!-- N-85 -->. The Ory normal-tissue signature had no overlap with either normal-tissue contrast, but three of its 40 tested genes occurred in RET-tumor. No Ory signature overlapped a BRAF q<0.10 list <!-- N-85 -->. These membership counts were descriptive and were not enrichment tests (Supplementary Results).

## Discussion

<!-- 2026-08-21 C2 ポート: gpt_review 試案 Discussion を移植(C1-c C-11 復元・凍結文再設置・タグ移植)。★アンカー同期は C4 -->
The principal finding was a broad expression difference between High-AS and dose-zero RET fusion-positive PTCs. By comparing cases within a driver-defined tumor class, the design prevented the overall RET-versus-BRAF composition from defining the primary contrast and directly tested whether expression structure associated with radiation-attributability bands remained after driver stratification. The prespecified Higher Criticism p-value of 0.0112 <!-- N-20 --> provided contrast-level evidence of label-aligned structure, while the Storey q<0.10 rule identified 1,765 genes for reporting <!-- N-16 -->. These results play complementary roles and constitute positive evidence for the working hypothesis within the prespecified design <!-- C-01 -->. They do not, by themselves, establish that the expression difference is specific to radiation. Throughout this section, readings drawn from the interpretation map (Methods, Table 1) were assigned to outcome patterns before the finalized analysis produced the results reported here; interpretations that go beyond the map are labeled as hypothesis-generating.

This driver-conditioned question differs from those addressed by earlier expression studies. Previous reports used tumor tissue, normal tissue, or paired tumor-minus-normal differences and varied in platform, exposure definition, and driver control (Abend et al. 2012, 2013; Hess et al. 2011; Dom et al. 2012; Ory et al. 2026). Cross-referencing showed limited correspondence: the validated anchors provided only one cross-tissue overlap, while the Ory shared-tissue and tumor signatures contained five and three genes, respectively, from the RET-tumor list <!-- N-85 -->. Given the breadth of that list, the multivariate selection of the Ory signatures in driver-unstratified exposure groups, and the absence of an enrichment test, these memberships provide context rather than independent support or validation, and do not establish a universal radiation signature <!-- C-13 -->.

The secondary contrasts define the scope of that finding rather than overturning it. RET-normal provided no detectable evidence for an accompanying expression difference in matched normal tissue, leaving the proposed long-term normal-tissue signal unresolved in this cohort <!-- C-03 -->. BRAF-tumor asked the parallel, direction-agnostic question within BRAF V600E-positive PTCs. Its limited evidence establishes neither equivalence between AS bands nor specificity of the RET-tumor result <!-- C-04 -->. In BRAF-normal, the single reported gene and the weaker contrast-level statistic diverged against a background of pronounced age and sex imbalance; this was also the only contrast with a calibration cell whose descriptive complete-null interval lay above the nominal level <!-- N-25 -->. This result was therefore treated as unresolved rather than as either confirmation or contradiction of the primary finding <!-- C-16 -->.

The functional annotation gives the RET-tumor result biological organization without supplying a separate proof. Under the gene-sampling reference used for over-representation analysis, genes expressed less highly in High-AS tumors were concentrated in proliferation, cell-cycle, and DNA-repair annotations. The radiation-curated flags belonged to the same overlapping cell-cycle theme and should not be counted as independent confirmations. No set crossed q<0.10 under the calibrated label-permutation analysis, which preserved subject-level grouping and gene dependence. Thus, the two analyses differ in the strength of inferential support they provide, but the over-representation result still identifies a coherent direction for mechanistic investigation <!-- C-14 -->.

The between-stratum comparisons add two further qualifications. Concordance between the normal-tissue profiles did not exceed its shuffled reference, but RET-normal itself contained no detectable structure; the absence of detectable sharing therefore cannot be distinguished from a shared pattern below the detection limit <!-- C-07 -->. By contrast, the RET and BRAF tumor profiles were more concordant than their shuffled reference. As a descriptive, hypothesis-generating observation, this is the pattern expected if part of the AS-band-associated expression structure extends across driver backgrounds. It is also compatible with covariates aligned with both tumor contrasts, so it motivates a cross-driver hypothesis rather than demonstrating an oncogene-independent radiation trace <!-- C-17 -->.

The REO analysis translated the group-level RET-tumor difference into panel-based per-case readouts. Separation of the construction bands showed that the selected pairs and cutoff represented their source contrast, which was necessary before applying them elsewhere but was not an estimate of predictive performance. Applying the panel unchanged to Low-AS and Mid-AS tumors extended the analysis beyond the extreme bands without refitting. The ordered median scores of 0, 1, 4, and 6 <!-- N-41 --> across dose-zero, Low-AS, Mid-AS, and High-AS groups were consistent with the hypothesized ordering. The specified Mid–Low comparison was also in that direction, although its uncertainty remained substantial (Pr(Low<Mid)=0.616; one-sided p=0.1127) <!-- N-40, C-08 -->. These observations neither establish a linear dose response nor validate a classifier, but they retain information that would be lost by interpreting the intermediate-band test alone as simply negative. The correlation between REO score and relative purity further indicates that the panel's biological basis remains to be resolved.

The main limitations concern attribution rather than the existence of the observed within-driver association. The design is also cross-sectional across driver-conditioned strata: an exposure-associated expression difference conditioned on driver cannot, in principle, distinguish a shared trace that feeds RET oncogenesis from an association between trace type and driver selection, and every claim in this paper is worded neutrally between the two <!-- C-11 -->. AS intentionally incorporates age at exposure as a radiation-risk modifier, while dose-zero and exposed cases also differ structurally in birth cohort and latency. The available data cannot separate those features from direct age- or cohort-related expression. Relative purity was expression-derived and used for eligibility, fusion partner may be part of the biological pathway, and sequencing- or processing-batch metadata were unavailable. We did not combine these variables in one adjusted model because they do not share a defensible interpretation as nuisance covariates: conditioning on inputs to AS would change the marginal estimand, whereas conditioning on purity or fusion partner could adjust for quantities affected by the process under study. Other constraints include the small High-AS BRAF group, transfer of the IREP external-photon thyroid model to internal iodine-131 exposure, reliance on the specified Storey estimator, and the exchangeability assumption underlying inference on the finalized normalized matrices. Exhaustive enumeration removed Monte Carlo error but did not constitute end-to-end randomization of normalization, and no independent-cohort or wet-laboratory replication was available. These limitations restrict causal and radiation-specific attribution; they do not convert the primary association into a null result.

In conclusion, High-AS and dose-zero RET fusion-positive PTCs differed broadly in expression after driver stratification. The coherent functional annotation, cross-driver tumor-profile concordance, and ordered REO profile provide complementary, although not independent, observations consistent with the working hypothesis that AS-associated expression structure remains detectable within a driver-defined tumor class. The present design cannot determine whether that structure is specific to radiation or is partly carried by age, birth cohort, purity, or technical factors. It nevertheless defines a documented, reproducible framework in which the biological and cross-driver implications of this association can be tested in future cohorts.

## Figure legends and table captions

<!-- 2026-08-21 C2 ポート: 試案キャプションへ更新。表1ブロック(旧・表2ブロック、2026-08-26 に引用順へ番号交換)(キャプション・行・規則)は正本の批准版を維持 — 位置づけ列・パターン規則は不可侵。対訳は廃止(f 保留方針) -->

**Figure 1 | Case flow.** From 440 <!-- N-11 --> REBC-THYR cases to the 63-case <!-- N-08 --> main cohort and 36-case <!-- N-10 --> REO intermediate-band application set, stratified by RET and BRAF driver.

**Figure 2 | Gene-level evidence by contrast.** Signed Brunner–Munzel effect, defined as 2θ−1, where θ=Pr(X<Y)+0.5Pr(X=Y), against −log10 permutation p. Genes meeting Storey q<0.10 are colored by direction.

**Figure 3 | REO reversal scores by AS band.** Scores from 0 to 10 for the 10-pair panel <!-- N-37 -->. Dose-zero and High-AS defined the panel; it was applied without refitting to Low-AS and Mid-AS. The four-band pattern is descriptive and is not a test of a linear dose response.

**Table 1 | Interpretation map for the four exposure contrasts, fixed before the finalized analysis produced the results reported here.** Each contrast carries a prespecified standing and its basis; the pattern rules fixed together with the map are given as table footnotes. No outcome pattern requires post-hoc labeling.

<!-- 表1(旧表2)本体+脚注規則(正本。Methods: Analysis contrasts から 2026-08-18 移設 — 投稿整形時に一方向切り出し、編集可能な第2コピーは作らない。規則の (a)-(f) 脚注化は研究者の書き直し待ち) -->
| Contrast | Prespecified standing | Basis |
| --- | --- | --- |
| R_Tumor | Hypothesis-bearing; signal expected | primary expectation of the working hypothesis (Introduction); the driver-conditioned gap left by Morton et al. (2021) |
| R_Normal | Hypothesis-bearing; signal possible | long-term molecular memory of exposure — dose-dependent expression within exposed normal tissue (Abend et al. 2013); exposed-versus-unexposed normal-tissue signatures (Ory et al. 2026) |
| B_Normal | Test cell for a shared trace | counterpart of the cross-stratum normal comparison; the shared signatures of Ory et al. (2026) were reported without driver stratification |
| B_Tumor | Direction-agnostic; not excludable | point-mutation drivers become less likely as dose rises in the source cohort (Morton et al. 2021) |

Pattern rules, fixed with the map: signal in R_Tumor is read as agreement with prior expectation; tumor signal confined to the BRAF stratum would be hypothesis-discordant and read with non-etiological explanations first; concordant signal in both normal contrasts is read as consistent with a shared glandular trace, discordant signal as suggesting a driver-linked trace type; normal-tissue signal confined to the BRAF stratum is read as confounding first; confined to the RET stratum, as unable to distinguish attenuated sharing from a driver-linked type; and an all-null outcome is reported as a bounded null. No outcome pattern requires post-hoc labeling.

**Table 2 | Case characteristics.** Group size, sex, age at surgery and age at exposure, and RET fusion-partner composition; footnotes give the age-difference estimates with confidence intervals and the group medians of relative tumor purity.

*Table 2 footnote.* Age-at-surgery differences compare High-AS with dose-zero cases. In the RET stratum, the Hodges–Lehmann difference was +2.5 years (95% percentile bootstrap CI, −1.0 to 6.0) and the Brunner–Munzel relative effect was θ=0.625 (0.400 to 0.828) <!-- N-64 -->. In the BRAF stratum, the corresponding estimates were +8.0 years (3.0 to 12.0) and θ=0.850 (0.681 to 0.973) <!-- N-65 -->. Here, θ=Pr(X<Y)+0.5Pr(X=Y), with X denoting dose-zero and Y High-AS. Intervals were obtained from 9,999 resamples drawn separately within each group (seed 19450809) <!-- N-63 -->; no p-values were calculated. Relative tumor purity is a within-cohort relative score (cohort maximum = 1) estimated separately per driver cohort and comparable within a stratum only; group medians were 0.783 (RET dose-zero), 0.822 (RET High-AS), 0.836 (BRAF dose-zero), and 0.922 (BRAF High-AS) <!-- N-88 -->.

**Table 3 | Gene-level results.** Tested genes, π0, q<0.10 counts and directions, minimum permutation p, and Higher Criticism p by contrast.

## 引用文献(整理 2026-08-14 — 仮書式、投稿書式の整合は後工程)

凡例: 全行とも書誌確認済み(PubMed または Crossref 照合、2026-08-14 完了。★は解消)。各行末は挿入位置。

### A. 本文が現に引くもの

- Morton LM, et al. Radiation-related genomic profile of papillary thyroid carcinoma after the Chernobyl accident. Science 2021;372:eabg2538. doi:10.1126/science.abg2538 — コホート出典(Methods: Data sources / Results §1)
- Storey JD. A direct approach to false discovery rates. J R Stat Soc B 2002;64:479–498. doi:10.1111/1467-9868.00346(Crossref 照合済み 2026-08-14)— λ=0.5 は同論文自身の計算で用いられた固定値(既定値ではない — 2026-08-16 原文照合で訂正。§9 の適応的選択は不使用。初出引用は Methods: Normalization の DEGES スクリーン、選定事由は Gene-level — 本文反映済み)
- Donoho D, Jin J. Higher criticism for detecting sparse heterogeneous mixtures. Ann Stat 2004;32:962–994. doi:10.1214/009053604000000265 — HC の選定事由(Methods: Gene-level)
- Wasserstein RL, Lazar NA. The ASA statement on p-values. Am Stat 2016;70:129–133. doi:10.1080/00031305.2016.1154108 — 非2値化(Methods: Gene-level)
- Greenland S, Senn SJ, Rothman KJ, Carlin JB, Poole C, Goodman SN, Altman DG. Statistical tests, P values, confidence intervals, and power: a guide to misinterpretations. Eur J Epidemiol 2016;31:337–350. doi:10.1007/s10654-016-0149-3 — 同上
- Amrhein V, Greenland S, McShane B. Scientists rise up against statistical significance. Nature 2019;567:305–307. doi:10.1038/d41586-019-00857-9 — 同上
- Goeman JJ, Bühlmann P. Analyzing gene expression data in terms of gene sets: methodological issues. Bioinformatics 2007;23:980–987. doi:10.1093/bioinformatics/btm051 — sampling-model 軸・GSEA/ORA の濃淡開示(Methods: Gene-set)
- Chatsirisupachai K, Lesluyes T, Paraoan L, Van Loo P, de Magalhães JP. An integrative analysis of the age-associated multi-omic landscape across cancers. Nat Commun 2021;12:2345. doi:10.1038/s41467-021-22560-y(PubMed 照合済み 2026-08-16、pmid 33879792)— TCGA 甲状腺癌の年齢対発現(fold change フィルタなし、BH 調整 p<0.05 で年齢関連遺伝子なし)(Discussion: C-14 段落)

### B. Methods の手法・実装引用(挿入位置は本文に反映済みまたは執筆時に確定)

**処理方針(2026-08-14 確定)**: (1) 手法の原典は、パッケージ使用か自作再実装かに関わらず初出箇所で引用する。(2) パッケージは実使用のもののみ名前を本文に出し、網羅は Supp の版表(docker/versions.tsv + session_info から機械生成)が担う。(3) 核心統計機構は自作実装であることを Methods(Software 節)で集中宣言 — コードが参照実体。各行の注記: 〔手法引用/自作実装〕= 論文は方法の典拠、実装はコード。〔使用パッケージ〕= 実際にそのパッケージを使用。

検定・多重性:

- Brunner E, Munzel U. The nonparametric Behrens-Fisher problem: asymptotic theory and a small-sample approximation. Biom J 2000;42:17–25. doi:10.1002/(SICI)1521-4036(200001)42:1<17::AID-BIMJ17>3.0.CO;2-U(Crossref 照合済み)— BM 検定の原典〔手法引用/自作実装 — exact/MC 置換枚挙。CRAN の漸近実装は不使用(棄却記録は手元台帳 T-08)〕(Methods: Gene-level — 本文反映済み)
- Benjamini Y, Hochberg Y. Controlling the false discovery rate. J R Stat Soc B 1995;57:289–300. doi:10.1111/j.2517-6161.1995.tb02031.x(Crossref 照合済み)— BH 調整(Methods: Gene-set)
- Hodges JL, Lehmann EL. Estimates of location based on rank tests. Ann Math Stat 1963;34:598–611. doi:10.1214/aoms/1177704172(Crossref 照合済み)— HL 推定量(Methods: Candidate confounders。慣用のため省略可)
- Nakayama Y, Yata K, Aoshima M. Test for high-dimensional outliers with principal component analysis. Jpn J Stat Data Sci 2024;7:739–766. doi:10.1007/s42081-024-00255-0(Crossref 照合済み 2026-08-15)— PC-OD の手法原典〔手法引用/自作実装 — lib/qc_pc_od.R〕(Methods: QC — 本文反映済み)

セット解析:

- Subramanian A, et al. Gene set enrichment analysis. PNAS 2005;102:15545–15550. doi:10.1073/pnas.0506580102 — 標準 GSEA 統計量・MSigDB〔手法引用/自作実装 — tie-block 拡張、tie-free 一致は自動テストで強制〕(Methods: Gene-set — 本文反映済み)
- Liberzon A, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst 2015;1:417–425. doi:10.1016/j.cels.2015.12.004 — Hallmark(Methods: Gene-set)

前処理・正規化・純度(いずれも Methods: Data sources / Normalization / QC):

- Dobin A, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 2013;29:15–21. doi:10.1093/bioinformatics/bts635 — STAR counts の来歴(GDC パイプライン言及の形でも可)
- Grossman RL, et al. Toward a shared vision for cancer genomic data. N Engl J Med 2016;375:1109–1112. doi:10.1056/NEJMp1607591 — GDC
- Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 2010;26:139–140. doi:10.1093/bioinformatics/btp616— filterByExpr / DGEList / CPM〔使用パッケージ〕
- Chen Y, Lun ATL, Smyth GK. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. F1000Res 2016;5:1438. doi:10.12688/f1000research.8987.2(Crossref 照合済み 2026-08-15)— filterByExpr の規則の出典(Methods: QC 初出 — 本文反映済み。以降の filterByExpr 出現は同規則)
- Kadota K, Nishiyama T, Shimizu K. A normalization strategy for comparing tag count data. Algorithms Mol Biol 2012;7:5. doi:10.1186/1748-7188-7-5 — DEGES
- Sun J, Nishiyama T, Shimizu K, Kadota K. TCC: an R package for comparing tag count data with robust normalization strategies. BMC Bioinformatics 2013;14:219. doi:10.1186/1471-2105-14-219 — iDEGES(3反復)の系譜〔手法引用/自作実装 — TCC パッケージ不使用、スクリーンは BM に置換(本文明記済み)〕
- Feng Y, Li LM. MUREN: a robust and multi-reference approach of RNA-seq transcript normalization. BMC Bioinformatics 2021;22:386. doi:10.1186/s12859-021-04288-0 — MUREN(※書誌照合で年・巻を訂正: 2020;21 → 2021;22)〔手法引用/自作実装〕(Methods: QC 初出で引用 — 本文反映済み)
- Shen Q, et al. contamDE: differential expression analysis of RNA-seq data for contaminated tumor samples. Bioinformatics 2016;32:705–712. doi:10.1093/bioinformatics/btv657 — contamDE 原法(比率スコアの親手法; 同論文の反復擬似尤度推定量は不使用)〔手法引用〕(本文反映済み)
- Ji Y, Yu C, Zhang H. contamDE-lm: linear model-based differential gene expression analysis using next-generation RNA-seq data from contaminated tumor samples. Bioinformatics 2020;36:2492–2499. doi:10.1093/bioinformatics/btaa006(PubMed 照合済み 2026-08-16、pmid 31917401)— 実装した比率推定量(最大1正規化)の出典〔手法引用/自作実装 — lib/purity_contamde.R〕(Methods: QC — 本文反映済み)

被曝指標・その他:

- Kocher DC, et al. Interactive RadioEpidemiological Program (IREP): a web-based tool for estimating probability of causation/assigned share. Health Phys 2008;95:119–147. doi:10.1097/01.HP.0000291191.49583.f7 — NIH IREP(Methods: Exposure metric。使用版 5.7.3。全文照合 2026-08-16: 放射性ヨウ素への言及ゼロ — 内部被曝の正当化は単独では担えない、下記 Land・Sato・IREP 文書と分担)
- Land CE, Gilbert E, Smith JM, et al. Report of the NCI-CDC Working Group to Revise the 1985 NIH Radioepidemiological Tables. NIH Publication No. 03-5387. Bethesda, MD: National Institutes of Health; 2003. — REF 表(electrons E>15keV = 単一値 1.0; 放射性核種β線は平均エネルギーで区分に割り付け — 脚注 c)・DDREF 規則(chronic 全域/acute は D_L = 0.03–0.2 Gy 対数一様以上で 1)・甲状腺モデル来歴(全文 PDF 逐語照合 2026-08-16 完了、md5 は N-83 — 再確認札は解除)(Methods: Exposure metric)
- NIH-IREP v5.7.3 オンライン文書「Guidance on Selection of Radiation Type — Internal Exposure」(radiationcalculators.cancer.gov/irep/ 配下、参照 2026-08-16)— ヨウ素131の electrons E>15keV 指定の典拠(2回名指し)。投稿書式での web 引用整形は後工程(Methods: Exposure metric)
- Sato T, Manabe K, Hamada N. Microdosimetric analysis confirms similar biological effectiveness of external exposure to gamma-rays and internal exposure to 137Cs, 134Cs, and 131I. PLoS One 2014;9:e99831. doi:10.1371/journal.pone.0099831(PubMed 照合済み 2026-08-16、pmid 24919099)— 内部ヨウ素131の RBE ≈ 1(保守的仮定下でも最大 1.04)の独立した物理的支持(Methods: Exposure metric)
- Ron E, Lubin JH, Shore RE, et al. Thyroid cancer after exposure to external radiation: a pooled analysis of seven studies. Radiat Res 1995;141:259–277.(PubMed 照合済み 2026-08-16、pmid 7871153、DOI なし)— IREP 甲状腺モデルの適合データ来歴(Methods: Exposure metric)
- Zurnadzhy L, et al. Clinicopathological implications of the BRAF V600E mutation in papillary thyroid carcinoma of Ukrainian patients exposed to the Chernobyl radiation in childhood: a study for 30 years after the accident. Front Med (Lausanne) 2022;9:882727. doi:10.3389/fmed.2022.882727(PubMed 照合済み 2026-08-16、pmid 35665338)— 同集団への先行 IREP 適用・acute 規約の来歴(**規約の先例としてのみ引用** — 共著者の先行研究のため等価性の根拠には使わない)(Methods: Exposure metric)
- Bogdanova T, et al. The relationship of the clinicopathological characteristics and treatment results of post-Chornobyl papillary thyroid microcarcinomas with the latency period and radiation exposure. Front Endocrinol (Lausanne) 2022;13:1078258. doi:10.3389/fendo.2022.1078258(PubMed 照合済み 2026-08-16、pmid 36589808)— 同上(本研究と同一の AS 出力行「expected value of ERR」を使用)(Methods: Exposure metric)
- Geman D, et al. Classifying gene expression profiles from pairwise mRNA comparisons. Stat Appl Genet Mol Biol 2004;3:Article19. doi:10.2202/1544-6115.1071 — 検体内順位比較(REO)の原型〔概念引用/パネル構築規則は自作〕(Methods: REO — 本文反映済み)
- Wang H, Sun Q, Zhao W, Qi L, et al. Individual-level analysis of differential expression of genes and pathways for personalized medicine. Bioinformatics 2015;31:62–68. doi:10.1093/bioinformatics/btu522(PubMed/Crossref 照合済み 2026-08-15、pmid 25165092)— 安定参照+逆転スキーム(RankComp)・正規化フリー性の系譜(REO の用語は同論文にない — 2026-08-16 照合で訂正、用語の典拠は Guan 2016)〔概念引用/パネル構築規則は自作〕(Methods: REO — 本文反映済み)
- Guan Q, Chen R, Yan H, Cai H, Guo Y, Li M, Li X, Tong M, Ao L, Li H, Hong G, Guo Z. Differential expression analysis for individual cancer samples based on robust within-sample relative gene expression orderings across multiple profiling platforms. Oncotarget 2016;7:68909–68920. doi:10.18632/oncotarget.11996(PubMed 照合済み 2026-08-16、pmid 27634898)— REO の用語の典拠(略号つき初出は Cai 2015, doi:10.18632/oncotarget.6260 — 引用は Guan の正準的定義のみ)(Methods: REO — 本文反映済み)
- Nikiforov YE, Rowland JM, Bove KE, Monforte-Munoz H, Fagin JA. Distinct pattern of ret oncogene rearrangements in morphological variants of radiation-induced and sporadic thyroid papillary carcinomas in children. Cancer Res 1997;57:1690–1694.(PubMed 照合済み 2026-08-15、pmid 9135009、DOI なし)— 放射線関連 PTC とパートナー特異性(RET/PTC3 優位)の原典・旧称の橋渡し(Methods: Candidate confounders — 本文反映済み)
- Rabes HM, Demidchik EP, Sidorow JD, et al. Pattern of radiation-induced RET and NTRK1 rearrangements in 191 post-Chernobyl papillary thyroid carcinomas: biological, phenotypic, and clinical implications. Clin Cancer Res 2000;6:1093–1103.(PubMed 照合済み 2026-08-15、pmid 10741739、DOI なし)— 同上・最大古典シリーズ(Methods: Candidate confounders / Intro — 本文反映済み)
- Ito T, et al. In vitro irradiation is able to cause RET oncogene rearrangement. Cancer Res 1993;53:2940–2943.(PubMed 照合済み 2026-08-15、pmid 8319199、DOI なし)— RET 再配列の実験的誘導の原典(Intro — 本文反映済み)
- Nikiforova MN, et al. Proximity of chromosomal loci that participate in radiation-induced rearrangements in human cells. Science 2000;290:138–141. doi:10.1126/science.290.5489.138(PubMed 照合済み 2026-08-15、pmid 11021799)— RET と融合相手の核内空間隣接(機構論)(Intro — 本文反映済み)
- Caudill CM, et al. Dose-dependent generation of RET/PTC in human thyroid cells after in vitro exposure to gamma-radiation. J Clin Endocrinol Metab 2005;90:2364–2369. doi:10.1210/jc.2004-1811(PubMed 照合済み 2026-08-15、pmid 15671095)— 誘導の線量依存性(Intro — 本文反映済み)
- Hamatani K, et al. RET/PTC rearrangements preferentially occurred in papillary thyroid cancer among atomic bomb survivors exposed to high radiation dose. Cancer Res 2008;68:7176–7182. doi:10.1158/0008-5472.CAN-08-0293(PubMed 照合済み 2026-08-15、pmid 18757433)— 独立曝露設定(原爆・外部γ線)でのヒト線量関連(Intro — 本文反映済み)
- R Core Team (2026). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/(citation()> 出力どおり 2026-08-16、rebc-r453:refblas 内 R 4.5.3 (2026-03-11) — 版数の書式は投稿誌指定に合わせ後工程)— (Methods: Software)
- その他の使用パッケージ(SummarizedExperiment・limma・GenomicDataCommons・rtracklayer・data.table・Matrix・Rcpp・MASS・parallel)は Supp の版表が網羅 —個別引用は原典論文を持つもののみ執筆時に判断〔使用パッケージ〕。事前収集済みの原典(PubMed 照合済み 2026-08-14): Bioconductor/SummarizedExperiment = Huber W, et al. Nat Methods 2015;12:115–121, doi:10.1038/nmeth.3252 / limma = Ritchie ME, et al. Nucleic Acids Res 2015;43:e47, doi:10.1093/nar/gkv007 / rtracklayer = Lawrence M, et al. Bioinformatics 2009;25:1841–1842, doi:10.1093/bioinformatics/btp328

### C. Intro / Discussion で必要になる見込み

外部アンカー(C-13。DOI は diagnostics/external_gene_anchors.csv の curation 由来):

- Abend M, et al. PLoS One 2012;7:e39103. doi:10.1371/journal.pone.0039103(PubMed 照合済み)— ペア差リスト(11 遺伝子)
- Abend M, et al. Iodine-131 dose-dependent gene expression. Br J Cancer 2013;109:2286–2294. doi:10.1038/bjc.2013.574(PubMed 照合済み)— 正常8/腫瘍6
- Dom G, et al. A gene expression signature distinguishes normal tissues of sporadic and radiation-induced papillary thyroid carcinomas. Br J Cancer 2012;107:994–1000. doi:10.1038/bjc.2012.302(PubMed 照合済み)— 正常7 = SERPINE1・DUSP1・TRIB1・S100A10・ANXA1・GNAL・RDH12(原文 Results と Fig. 2 legend で照合 2026-08-16、diagnostics/external_gene_anchors.csv と 1:1 一致 — 札解除。clade E の番号は SERPINE1 と legend に明記)
- Hess J, et al. Gain of chromosome band 7q11 in papillary thyroid carcinomas of young patients is associated with exposure to low-dose irradiation. PNAS 2011;108:9595–9600. doi:10.1073/pnas.1017137108(PubMed 照合済み)— CLIP2 (追加検証: doi:10.1038/onc.2014.311)

正常組織・対照・仮説:

- Ory C, et al. Post-Chornobyl thyroid papillary carcinomas display distinct past exposure and radiation-associated carcinogenesis molecular signatures at low and high thyroid doses. Sci Rep 2026;16. doi:10.1038/s41598-026-53030-4(PubMed 照合済み)— 正常組織対照・driver 非層別の隙間(Intro 予告 / Discussion C-07)
- Efanov AA, et al. Investigation of the relationship between radiation dose and gene mutations and fusions in post-Chernobyl thyroid cancer. J Natl Cancer Inst 2018;110:371–378. doi:10.1093/jnci/djx209(PubMed 照合済み)— driver 組成と線量の共変(Intro の層別根拠。Q-01)

---

## 起草メモ(統合 — 本文に載せない)

- **見出し階層の正規化(2026-08-14)**: H1 = 文書タイトルのみ、節 = H2、下位 = H3。
- **表記統一(2026-08-14)**: tumour → tumor(米式、混在19/47の多数派に統一 — 学術誌確定時に
  英式が必要なら機械置換で戻せる)。contamDE の大小文字も統一済み。「signature agreement」は
  本文不使用(Methods 節スキャフォールドの注記のみ)。
- **改行方式(2026-08-14、研究者選択)**: 本文散文は段落=1行(エディタの折返しで読む)。構造行(見出し・リスト・単独タグ行)のみ独立行。固定幅折返しは機能的根拠がないため廃止 — 以後の起草もこの方式で書く。
- **Methods への【訳】付与(2026-08-14、研究者指示)**: 全13小節に対訳を挿入(Results と同形式・N/C タグ併記。検査対象は英語本文)。
- **パッケージ・実装の処理方針(2026-08-14、研究者指示)**: 手法の原典は実装形態に関わらず初出で引用/パッケージは実使用のみ本文に名前、網羅は Supp 版表(Supp.Tab.4(仮)= docker/versions.tsv + session_info)/核心統計機構の自作実装は Software 節で集中宣言(BM は個別にも宣言 — CRAN 漸近実装の棄却は T-08、一次ログの所在を凍結前に特定)。本文5箇所+訳5箇所に反映済み。
- **引用文献の整理(2026-08-14、研究者指示)**: 3分類(A = 本文が現に引く / B = 手法・実装 / C = Intro・Disc 見込み)+挿入位置+★印(書誌の最終確認要 —凍結前に PubMed 照合)。外部アンカーの DOI は curation CSV から転記(確認済み扱い)。★付き実装系(MUREN・contamDE・TCC ほか)は書誌の細部に私の記憶依存があるため必ず照合すること。
- **統合の記録(2026-08-14、研究者指示)**: 旧4ドラフトを本ファイルへ吸収し削除。各段落の実体は本ファイルの1箇所のみとなり、正本/再掲の同期機構は廃止。以下のサブ節は旧ファイルの起草メモをそのまま移設したもの(歴史込み — 「正本はあちら」等の記述は統合前の状態を指す)。

### Methods 由来(旧 draft_methods.md)

- **未タグの数値ゼロを目標に、設定値・規則を N-66〜N-76 として台帳化した** (numbers_ledger セクション P。出典は全てコミット済みコードの行参照、計算なし)。
- D6 段落と ORA 段落は正本(draft_methods_results_d6.md / draft_ora_annotation.md)からの再掲。編集はあちらで行い、本稿へは機械的に同期する(検査トリガー1で照合)。
- 「Analysis contrasts」節は意図的に A/B 中立・機構言及なしの定義のみ。仮説の文言(radiogenic/sporadic の二値観を含む)は研究者領分のため Intro/Discussion 側に置く想定で、Methods には書いていない。
- 530 の記述は「graded, descriptive, no dose–response form assumed」で§0.5 第1回の水準制約に合わせた。「if the panel captures a radiation-attributable signal…」型の予測文(530 ヘッダにある)は Methods に載せず、載せるなら Intro の予告側。
- 純度プーリング(研究者決定 2026-08-13、同日改訂): script ヘッダ・Methods とも**理論構造で記述** — 相対尺度の共通化(群別 run は単一閾値の意味を壊す。群別で先に回して観察された事実)、共通参照仮定は driver 層別設計の前提に乗る、純度軸は曝露対比と別軸、役割は相対フィルタ+診断共変量のみ。設計時実測(0.93–0.99 / 0.99)は**二重統計のためライセンスに使わない** — run コミット(8eed384)の凍結ヘッダと N-76 行(不使用)にのみ保存。
- 手術年齢の共変量非投入の理由文(C-10 段落)は Q-03 の要旨の圧縮。年齢層別診断をしない理由の詳細は Discussion/limitation 側(Q-03/Q-13 — 判断点4は決着済み)に置く想定。
- **選定事由の記載方針(2026-08-14 確定)**: 分野の既定から外れる選択で査読者が「なぜ?」と聞くと予想される箇所に1文の事由を書く。慣行どおりの選択には書かない。
- **深度配分の原則(2026-08-14 確定、Mid>Low の往復から一般化)**: Results に載る全ての検定・read について「それが何を検定するのか」を Methods が1文で持ち、Results は結果+Methods 参照に留める。この原則で4点を追補: HC の役割文(対比水準の主張を担う — Q-16 の水準論を本文に可視化)/ ORA の濃淡開示一文(Q-15 決着の実装)/ 節題を「Between-arm concordance of the exposure contrast」へ(判断点2の用語決定; 2026-08-14 の用語統一 arm→group/stratum で「Between-stratum concordance」に改称)/ REO 診断 (iii) の動機文(純度駆動の可能性の分離)。この基準で4点を追記: BM 検定の採用(順序関係を原始とするコミットメント+完全枚挙)/ HC の選定(多数の弱い効果への感度)/ tie-block ES の理由(順序の注入禁止)/セットレベル π0=1 の非対称(依存下で plug-in が分散支配・反保守側に誤る)。n_perm 9,999 は慣行的値のため本文で事由を書かない(D4 の導出は手元)。
- λ=0.5 の選定事由を 410 節に追記(2026-08-14、研究者指示): D1 の批准済み導出の圧縮(π0=1 は仮説の検出チャネルを塞ぐ/固定 λ の保守性は帰無 p の周辺一様性のみで成立/適応的 λ は独立性仮定で適用外/平坦な対立密度下で低 λ が MSE 優位)。弱拡散前提は**選定事由(設計時の作業仮説)として書く**のであって形の主張ではない — Q-15 (2) (形に賭けない)と整合。「small per-gene effects attenuated by within-group case mixture」の文言はバイナリ観の反映だが、仮説の正式な文言化(Intro、研究者領分)と執筆時に整合させること。
- 非2値化の一文(C-16)を 410 節末尾に追加(2026-08-14、判断点1決着)。引用は ASA 声明系3本(DOI は claim_map C-16 の根拠列)。対比の陽性/非陽性ラベルは全対比で不使用 — C-01 は両水準結合(DEG かつ HC)、C-03/C-04 は「検出されなかった」の記述形。文言検査時に positive/negative 型のラベル語が紛れていないか確認する。
- 決定履歴語の掃引(2026-08-14、研究者指示): 「protocol amendment」「decision record」への言及を本文から削除(査読者は決定の行き来を考慮しない。再現性文は「code and versioned inputs sufficient to regenerate」に限定 — 公開対象の選別は研究者裁量、§0.5b 明確化と整合)。
- **予測マップ+解釈規則は論文内に提示が必要**(正本=論文 — 外部計画への参照で済ませない。非2値化の一文が「interpretation map, Methods」を指すため)。置き場所(Methods 内の小表 or Intro)は図表構成(判断点5)と同時に確定。(2026-08-15 決着: Methods 内小表〔表5(仮)、Analysis contrasts 末尾〕として掲載済み。Intro はフルスクラッチ起草済み — 仮説スロット1箇所が研究者記入待ち)
- IREP の入力規約は 130 ヘッダの写し(N-68)。データ提供側の線量そのものの来歴は REBC-THYR 原論文への引用で受ける(引用は Intro/Methods 冒頭、書誌は研究者)。
- Table/Supplementary 番号は全て(仮)。図表構成の確定(判断点5)後に振り直す。
- リント自己検査(10項目): (1) 対比横断 FDR 主張なし(「no study-wide FDR is claimed」と明文化)/ (2) A/B 中立(機構文言なし)/ (3) 線量依存の仮定なし(530 は descriptive・no dose–response form)/ (4) q<0.25 不使用 / (5) 特異性確認の文言なし / (6) 探索セル由来なし(ORA は hypothesis-generating と明示)/ (7) 年齢の p>0.05 型議論なし(推定+CI のみ、p 値なしを明記)/ (8) 全数値 N タグ済み(新規設定は N-66〜N-76)/ (9) WY-FWER 不使用 / (10) 主張水準は手続き記述のみで超えない。

### Results 由来(旧 draft_results.md)

- 英語版が正本、【訳】は研究者確認用の対訳(タグは両方に付けたが検査対象は英語版)。
- 再掲2段落の正本: C-06 = draft_methods_results_d6.md、C-14 = draft_ora_annotation.md。編集は正本側で行い本稿へ同期(トリガー1で照合)。
- **C-14 正本の数値取り違えを修正**(2026-08-14): radiation 例示が down の k(46/126)に combined の期待値(14.2)を付けていた → down の期待 6.4 に修正(N-60)。正本側も修正済み。
- 丸めの規約: rho は3桁表示(+0.376 = N-34 の +0.3756、帯 [−0.46, +0.46] = 同[−0.4615, +0.4580])。他の数値は台帳の桁のまま。凍結時に表示桁を台帳へ追記予定。
- B_Normal(C-16)は評語なしの連続量記述で、読み(規則5条件適用)は Discussion 送り。「他の2つの静かな対比より低い π0」は N-19 の数値の並置で、比較検定ではない。
- tumor 側 concordance(C-17)は本文ポインタ+Supp. Tab. 3(仮)+Discussion 送り(判断点2決着どおり)。
- ORA の可視性は現状「短報1段落」— 共著者交渉で拡張可(C-14 の可視性無制約)。拡張時は正本(draft_ora_annotation)側を先に改稿。
- **フロー提示の方針(2026-08-14 研究者承諾)**: 駆動遺伝子別の途中経過(N-08 の RET 73→…→27 / BRAF 175→…→36)は本文でなく**コホートフロー図(両層並記)**が担う。本文は合算チェーン+最大削減2段の半文+フロー図参照まで。図の最終形は作図時(判断点5)に確定。
- **Mid>Low 検定の提示(2026-08-14 研究者指摘への対応)**: 検定の問い(訓練帯は循環・訓練外は Low/Mid のみ・方向は設計が事前指定 → 過適合でないことの out-of-sample 検査)は**Methods の深度**として Methods の REO 評価段落に記載。Results は近隣と同密度の 1句+Methods 参照のみ(深度アンマッチの回避)。検定自体は維持(事前指定 read の非掲載は選択的報告になるため。役割は「勾配の不確かさの開示」で、感度解析ではない)。
- **REO ペアの生物学(2026-08-14 研究者決定)**: 本文・Results ではペア名を挙げず表参照のみ(器械としての提示)。共著者は生物学的解釈を好む傾向 — 求めがあれば Discussion/共著者ラウンドで扱う余地は残すが、その際は P9/P10 が再正規化で入れ替わった経緯(選定順位の端の揺れ)を必ず添えること。
- 図表番号は全て(仮)。判断点5(図表構成)確定後に振り直し+図表台帳と照合。
- リント自己検査(10項目): (1) 対比横断 FDR 主張なし / (2) A/B 中立(機構文言なし)/ (3) 線量依存の仮定なし(C-08 に明示の否定文)/ (4) q<0.25 不使用 / (5) B_Tumor の「特異性確認」文言なし(明示の否定文)/ (6) 探索セル産なし(ORA は hypothesis- generating 明示、C-17 は Disc 送り)/ (7) 年齢は推定+CI のみ(検定なしを明示)/ (8) 全数値 N タグ済み / (9) WY 不使用 / (10) 主張水準は C 行の水準列どおり(主張 = C-01 のみ、他は記述的観察・仮説生成)。

### D6 二層(旧 draft_methods_results_d6.md)

- 二層の分割方針: 本文には「較正済み・pooled 0.064・超過1セル開示・spike-in 回収」の 4事実のみ(査読者がサプリを開く前に見えるべき情報)。機構・CI 注意・選定経緯はサプリ。
- 「0.045」は**採用時測定**(2026-08-08、開発 B=999、実データ前)で、本番凍結値は 16セル表+pooled 0.064(N-56)。前者は選定の時系列防御、後者が論文の較正値。
- WY-FWER 0.112 は N-57 に参考として台帳化したが本文・サプリとも不使用(B.2)。
- Supplementary Results 最終文の受け先: **確定(2026-08-14、判断点1決着)** — Discussion の B_Normal 段落(C-16 の規則5条件適用段落)の一文で受ける。B_Normal の扱いは2値ラベルなしの連続量記述+規則5条件適用(claim_map C-16、Q-16)。
- リント自己検査: 対比横断 FDR 主張なし / A/B 文言なし / 線量依存なし / q<0.25 不使用 /全数値 N タグ済み / C-06 タグ済み / WY 本文不使用。

### ORA 注釈(旧 draft_ora_annotation.md)

- **数値修正(2026-08-14、Results 起草時の照合で検出)**: Results 文の radiation 例示が down リストの k(46/126)に combined の期待値(14.2)を付けていた取り違えを修正 —正しくは down リストの k=46/126・期待 6.4(N-60。combined は k=50/126・期待 14.2)。
- **2026-08-13 来歴訂正・根拠差し替え(研究者承諾 — claim_map C-14 改訂メモ参照)**: ORA は Thyroid 原プロトコルの設計解析の復元であり、GSEA 主・ORA 副の濃淡の根拠は「事前固定」ではなく G&B 2007 の sampling-model 軸+本文での一文開示。開示文の例(Methods 併置文の直後を想定): "The two procedures' p-values refer to different randomness: only the label-permutation null refers to the experiment actually performed, and experiment-level claims rest on it; the over-representation q-values describe the discovered list against a gene-sampling reference."可視性は無制約 — 本稿の「短報2文」は**下限**であり、Results 段落規模・図・Abstract 言及への拡張は共著者交渉で可(水準ラベル携行のみ非交渉)。転落基準(ORA q を実験レベルの単独主張に使わない)は Q-15 (3b)。濃淡開示の一文は **draft_methods の C-14 段落に配置済み**(2026-08-14。文言は上の例文と同一、G&B 引用付き)。

- Methods 併置用の1文(420 Methods の末尾を想定): "As a descriptive complement, the discovered R_Tumor list was annotated by one-sided hypergeometric over-representation against the identical set universe (up / down / combined lists; universe = the contrast's tested genes; BH within family × list), reported in full as hypothesis-generating annotation (Supplementary Table SY)." <!-- C-14, N-59 -->
- ガードレール3点(C-14 に記録済み): 単一テーマの明示 / 年齢・潜伏期の候補説明の名指し / 420 の同順位・閾値未達の併記。
- 「420 の同順位」の具体: 較正済み検定でも H ファミリー最上位は同じセット群(SPERMATOGENESIS p 0.0036 q 0.179、KRAS_SIGNALING_UP・E2F_TARGETS・G2M_CHECKPOINT が NES 負の上位 — thyr_enrichment_test.rds; min q は N-28)。本文で個別 NES を引く場合は N 行の追加が必要(現状は N-28 の min q のみ台帳化)。
- 年齢・潜伏期の受け: **判断点4決着により確定(2026-08-14 反映)** — 年齢層別診断は実施せず(Q-03・撤回済み3案は復活させない)、開示された推定(C-15、N-64)で受ける(旧: 文献側 Vriens 2011 = Q-13 (i) — 2026-08-16 の引用監査で撤去)。Discussion の [ ] は埋め済み(N-64 タグ)。潜伏期は年齢と不可分のため候補説明としての名指しのまま。
- Supp Table SY = Supp.Tab.2(仮)(図表台帳)。D6 の Supp Table SX と番号整合は図表構成確定時に振り直し。
- リント自己検査: 仮説生成の水準明示あり / A/B 中立(機構主張なし)/ 線量依存なし / q<0.25 不使用 / 対比横断 FDR なし / 全数値 N タグ(NES 個別引用は保留)/ C-14 タグ済み。
