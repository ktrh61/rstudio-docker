# 原稿統合ドラフト(draft_manuscript — 論文へ収束する単一ファイル)

状態: draft(2026-08-14 統合。研究者指示: Methods/Results の対構造の検査可能化と同期機構の廃止のため、旧4ドラフト — draft_methods / draft_results / draft_methods_results_d6 / draft_ora_annotation — を本ファイルへ吸収し削除。内容の変更なし、移設+Discussion 骨組み追加のみ。履歴は git)。規約: 英語本文が正本。インライン【訳】は C2 ポートで段階的に廃止(対訳は内容確定後に共著者向け別ファイルとして再生成 — 2026-08-21 研究者方針)。Supplementary Methods/Results は paper/supplementary_material.md へ分離(2026-08-21、BJC 投稿構造)。数値は N-ID タグ、主張段落は C-ID タグ(検査対象は英語本文)。図表番号は 2026-08-16 に本置き済み(主: 図1–3・表1–3 / 補足: 図S1–S2・表S1–S8・データ1。BJC の Article 規定 = 主6品目以内に整合)。各節の由来・決定は末尾の「起草メモ(統合)」。

---

<!-- タイトルページ素材(BJC: title page 別建て・語数枠外) -->
**Title:** Driver-conditioned transcriptomic differences across radiation-attributability bands in papillary thyroid carcinoma

## Abstract

### Background

Although fusion drivers become more frequent with dose in post-Chornobyl papillary thyroid carcinoma (PTC), gene expression has not been consistently associated with exposure. We therefore examined this association within driver strata.

### Methods

In the REBC-THYR cohort, cases were stratified by driver, and the highest band of the NIH IREP Assigned Share (AS; High-AS) was compared with dose-zero cases within each stratum in tumor and matched normal tissue. The primary comparison was in RET fusion-positive tumors, where radiation-related tumor initiation was expected. Gene-level differences were tested with permutation Brunner–Munzel tests and summarized per comparison by Higher Criticism, and a relative-expression-ordering (REO) panel constructed from the extreme bands of the RET fusion-positive stratum was applied unchanged to the intermediate bands.

### Results

In the primary comparison, 1,765 of 15,621 tested genes met Storey q<0.10 <!-- N-16, N-15 --> (Higher Criticism p=0.0112 <!-- N-20 -->), whereas the secondary comparisons yielded at most one gene. Median REO scores in the intermediate bands lay between the dose-zero and High-AS values <!-- N-41 -->.

### Conclusions

The primary comparison provides positive evidence for the working hypothesis that AS-associated expression structure remains detectable within a driver-defined tumor class. The design offers a reproducible framework for testing the biological and cross-driver implications of this association in future cohorts.

## Background

<!-- 2026-08-21 C2 ポート: gpt_review 試案 Intro を移植(C1-e 仮説スロット批准・凍結文言は新正準へ統合)。予告列の同期は C4。インライン対訳は廃止(f 保留方針) -->

After the Chornobyl accident, people exposed to iodine-131 in childhood carried a dose-dependent excess risk of thyroid cancer that persisted for two decades (Brenner et al. 2011), and post-Chornobyl papillary thyroid carcinoma (PTC) remains a major human model of radiation carcinogenesis. Molecular studies have characterized these tumors at two levels, the driver and the transcriptome.

At the driver level, RET rearrangements are the best characterized. They were frequent in post-accident and in sporadic childhood PTC alike, but the fusion-partner distribution differed: NCOA4::RET (historically RET/PTC3) predominated in post-accident tumors and CCDC6::RET (RET/PTC1) in age-matched sporadic tumors (Nikiforov et al. 1997), a predominance that was most marked in the most heavily contaminated regions and shifted toward CCDC6::RET with longer latency (Rabes et al. 2000). The evidence for a radiation origin extends to both fusion partners: RET rearrangements were induced experimentally in human thyroid cells, with CCDC6::RET the more frequent of the two (Caudill et al. 2005), RET lies close to each partner locus in thyroid-cell nuclei (Nikiforova et al. 2000; Gandhi et al. 2006), and RET rearrangements were enriched in tumors from highly exposed atomic-bomb survivors (Hamatani et al. 2008). Evidence of the same kind exists for ETV6::NTRK3 rearrangements, which were identified in post-Chornobyl PTC (Ricarte-Filho et al. 2013; Leeman-Neill et al. 2014), were more frequent there than in sporadic PTC after adjustment for age and sex, and were induced in human thyroid cells in vitro (Leeman-Neill et al. 2014). At cohort scale, the REBC-THYR study of post-Chornobyl PTC with individual thyroid dose estimates found that fusion drivers as a class became more frequent as dose increased, whereas point-mutation drivers became less common (Morton et al. 2021).

At the transcriptome level, a gene expression signature distinguished radiation-induced from sporadic thyroid tumors after external radiotherapy (Ory et al. 2011), and in post-Chornobyl material, comparisons of exposed with sporadic tumors (Detours et al. 2005, 2007; Ugolin et al. 2011), dose-related expression within exposed cases (Abend et al. 2012, 2013), and copy-number-linked and normal-tissue expression studies (Hess et al. 2011; Dom et al. 2012; Ory et al. 2026) have each reported candidate differences, and the genes reported by the tumor comparisons showed little consistency across studies (Dom et al. 2012; Abend et al. 2013). At cohort scale, the REBC-THYR study found expression shaped by the driver rather than by radiation dose (Morton et al. 2021). Unlike the driver level, the transcriptome level has therefore yielded no consistent exposure-associated difference.

Two findings from these studies have implications for how an exposure-associated expression difference can be detected: PTC transcriptomes are strongly organized by driver, in sporadic PTC (Cancer Genome Atlas Research Network 2014) as in the REBC-THYR cohort (Morton et al. 2021), and driver composition itself shifts with dose. A comparison of exposed with unexposed cases pooled across drivers therefore confounds exposure with composition and can reproduce the difference between fusion- and point-mutation-driven tumors. None of the studies above, including the cohort-scale analyses, conditioned on driver, except the CLIP2 study (Hess et al. 2011), which matched RET-rearrangement status as far as possible but addressed a single gene.

We therefore examined, in the REBC-THYR cohort (Morton et al. 2021), whether exposure is associated with tumor expression within a driver stratum. Dose raises the risk of PTC, but a tumor in an exposed person need not be radiogenic (Kocher et al. 2008): exposed cases are a mixture of radiogenic and sporadic tumors in which dose changes the radiogenic proportion. Given that cohort-scale expression was shaped by the driver rather than by dose, a continuous-dose model averages over the two kinds of tumor and dilutes any exposure-associated difference, whereas the dose-zero group and the highest-attributability group differ most in radiogenic proportion. We hypothesized that comparing the highest band of a case-level radiation-attributability index with dose-zero cases might reveal expression structure that a continuous-dose model dilutes. The index was the Assigned Share (AS) from the NIH Interactive RadioEpidemiological Program (IREP; Kocher et al. 2008), calculated from recorded thyroid dose and age; AS served only to order exposed cases into bands, which, together with the dose-zero cases, form four ordered categories, and the cases of one band within one driver stratum form a group (Methods). Because AS is calculated without reference to the driver and fusion drivers became more frequent with dose (Morton et al. 2021), the High-AS group within the RET stratum may contain a larger radiogenic proportion than its AS range implies.

Only the RET fusion-positive and BRAF V600E-positive strata supported the comparison; the ETV6::NTRK3 stratum did not (Methods). We prespecified four comparisons of the High-AS group with the dose-zero group, in tumor and in matched normal tissue within each stratum, with the RET fusion-positive tumor comparison as primary because radiation-related tumor initiation was expected in that stratum and fusion drivers became more frequent with dose. The BRAF V600E-positive and normal-tissue comparisons were secondary. Table 1 gives the interpretation assigned to each possible pattern before the finalized analysis produced the results reported here. Separately, we constructed a relative-expression-ordering (REO) panel from the RET dose-zero and High-AS groups and applied it unchanged to Low-AS and Mid-AS RET fusion-positive PTCs.

## Methods

<!-- 2026-08-21 C2 ポート: gpt_review 試案 Methods を移植(C1-a AI開示復元・C1-d 目盛り要約復帰・chronic監査文は不採用)。写像検査・★アンカー同期は C4 -->
### Data sources

Gene-level STAR count files for REBC-THYR were downloaded from the National Cancer Institute Genomic Data Commons (Heath et al. 2021). The clinical source was Data S1 of Morton et al. (2021), containing 440 cases <!-- N-11 -->. Files were mapped to cases and biospecimens through the Genomic Data Commons API and assembled into a matrix of 58,448 genes with nonzero counts and 906 samples <!-- N-14 -->. Strandedness, gene-length derivation, file-selection, and manifest-verification procedures are described in Supplementary Methods.

### Assigned Share

For each exposed case, the NIH IREP thyroid model, version 5.7.3, was used to calculate the AS associated with the expected excess relative risk from the recorded thyroid dose and ages (Kocher et al. 2008). AS was calculated with the IREP input route for internal iodine-131 exposure, with exposure rate entered as acute and dose as a point value. The full input convention and its provenance are given in Supplementary Methods. Dose-zero cases formed the reference group, labeled Sporadic in the analysis files, and exposed cases were classified as Low (0<AS<33.3), Mid (33.3≤AS<66.6), or High (AS≥66.6), with no case lying on a boundary <!-- N-69 -->. Throughout this paper, band denotes one of the four AS-defined categories (dose-zero, Low-AS, Mid-AS, and High-AS) and group denotes the cases of one band within one driver stratum. Each contrast compared the dose-zero group with the High-AS group of a stratum on the basis of group membership alone. We used AS only to order exposed cases and to define the bands, and we derived no case-level causal assignment from it.

### Quality control and cohorts

Principal-component-scores-based outlier detection (PC-OD; Nakayama et al. 2024) was applied separately within each observed group-by-tissue matrix after filtering genes with edgeR's filterByExpr (Chen et al. 2016; Chen et al. 2025), with settings given in Supplementary Methods. After PC-OD, relative tumor purity was estimated jointly for dose-zero and High-AS cases within each driver cohort from matched tumor–normal counts with contamDE-lm (Shen et al. 2016; Ji et al. 2020) after MUREN normalization (Feng and Li 2021), placing both groups on a common relative scale. Purity served as a within-cohort eligibility filter and diagnostic variable.

The main cohort comprised cases that were RET fusion-positive or BRAF V600E-positive (without another candidate driver mutation; Supplementary Methods), were in the dose-zero or High-AS group, had one matched tumor–normal sample pair with both tissues passing the outlier screen, and had relative purity ≥0.6 <!-- N-70 --> (n=63 <!-- N-08 -->). Its RET fusion-positive subset supplied the REO construction groups, whereas RET fusion-positive PTCs in the Low-AS and Mid-AS bands formed the REO application set (n=36 <!-- N-10 -->). We did not apply the outlier and purity exclusions to this set and reported the corresponding metrics descriptively. Because the ETV6::NTRK3 stratum had only seven dose-zero cases (with 25 exposed cases) <!-- N-103 -->, we did not calculate AS for it and did not analyze it.

### Analysis contrasts

In this cross-sectional analysis, we defined four contrasts comparing High-AS with dose-zero cases in tumor and matched normal tissue within each driver stratum: RET-tumor, RET-normal, BRAF-tumor, and BRAF-normal (labeled R_Tumor, R_Normal, B_Tumor, and B_Normal in the tables and analysis files). Within each contrast, X denotes the dose-zero group and Y the High-AS group, and the Brunner–Munzel relative effect (probabilistic index; Brunner and Munzel 2000) is θ=Pr(X<Y)+0.5Pr(X=Y). Values of θ>0.5 indicate that expression tends to be higher in High-AS cases. Higher Criticism (HC; Donoho and Jin 2004) was the single primary contrast-level test for RET-tumor. We analyzed the three secondary contrasts in the same way, reported them separately, and applied no across-contrast adjustment. The interpretation assigned to each outcome pattern (Table 1) was fixed before the finalized analysis produced the results reported here.

### Covariates

For exposed cases, age at exposure entered the IREP model as a radiation-risk modifier. We report group distributions of age at surgery, with Hodges–Lehmann differences and 95% bootstrap intervals <!-- N-63 -->, together with sex composition and RET fusion-partner composition. We fitted no common covariate-adjusted model, for the estimand reasons given in Supplementary Methods. No sequencing- or processing-batch field was available in the analysis objects.

### Normalization and gene-level inference

Counts were restricted to protein-coding genes, filtered with filterByExpr, and normalized separately for each contrast using an in-house implementation of three-iteration DEGES-MUREN (Kadota et al. 2012; Sun et al. 2013; Feng and Li 2021). For scaling-factor estimation, we replaced the original model-based differential-expression screen with a permutation Brunner–Munzel screen at Storey q<0.10 (Storey 2002) and excluded the genes it flagged as potentially differentially expressed <!-- N-71 -->. Supplementary Methods and Table S1 report the detailed filtering, floorPDEG rule, convergence, and scaling factors <!-- N-15 -->.

For each gene, a studentized Brunner–Munzel statistic (two-sided; Neubert and Brunner 2007) was evaluated for all C(n,nx) label allocations on the finalized normalized matrix. Gene-level lists used the Storey q-value rule at q<0.10 <!-- N-05 -->. The estimator settings are given in Supplementary Methods. Contrast-level evidence was summarized with a Higher Criticism statistic scanning the smallest 10% of the ordered p-values, with its p-value obtained from 9,999 label shuffles within the contrast <!-- N-04, N-05 -->. We report these p-values as continuous evidence and do not assign a binary label to each contrast. <!-- C-16 -->

### Gene-set inference

Genes were ranked by tie-averaged normal scores of the signed Brunner–Munzel statistic, and an in-house weighted running-sum statistic was evaluated only at tie-block boundaries. On tie-free input, this statistic equals the standard gene set enrichment analysis statistic (Subramanian et al. 2005) <!-- N-72 -->. The 9,999 saved contrast-specific label shuffles <!-- N-07 --> formed the set-level reference distribution, and sign-conditional permutation p-values were adjusted within each collection using the Benjamini–Hochberg procedure at q<0.10 (Benjamini and Hochberg 1995).

The collections tested were Hallmark, selected C2 canonical-pathway subcollections (Reactome, WikiPathways, KEGG MEDICUS, BioCarta, and PID), C5 Gene Ontology Biological Process, and a radiation-curated C2:CGP family <!-- N-29, N-55 -->. Sets outside 15–500 genes were excluded <!-- N-72 -->. <!-- C-06 -->
The set-level procedure was selected on a held-out complete-null assessment before the finalized real-data run. Its operating characteristics were assessed with a 9,999-shuffle reference pool <!-- N-06 --> and examined in a single positive-control scenario with a one-set 1.15-fold spike-in <!-- N-30 --> (Supplementary Methods).

<!-- C-14 -->
As a descriptive complement, the primary RET-tumor q<0.10 list was tested for one-sided hypergeometric over-representation in the same collections, with Benjamini–Hochberg adjustment within each family-by-list combination <!-- N-59 -->, with the gene universe specified in Supplementary Methods. Because these gene-sampling p-values do not represent subject-level randomization, experiment-level interpretation rests on the label-permutation analysis (Goeman and Bühlmann 2007).

### Concordance between the RET and BRAF contrasts within each tissue

For each tissue, signed gene-level Brunner–Munzel statistics from the RET and BRAF contrasts were compared across shared genes using Spearman correlation, with a reference interval generated by independently shuffling labels within each stratum. The normal-tissue comparison addressed the question specified in Table 1 of whether a normal-tissue change is shared across driver strata, whereas the tumor comparison was included for completeness.

### REO panel

The REO panel was based on within-sample gene ordering, which does not require between-sample normalization (Geman et al. 2004; Wang et al. 2015; Guan et al. 2016), using transcripts per million recalculated as described in Supplementary Methods. We generated candidate pairs from the 500 RET-tumor genes with the largest Brunner–Munzel effect magnitudes, rather than from the q<0.10 list. The selection rules required a consistent dose-zero sign outside a dead zone of |log2 ratio|<log2(1.2) with at most one exception, a 10th-percentile magnitude of at least log2(1.5), and a High-AS reversal rate above 50% but below 100%. Pairs that met these rules were ranked by the shift in median within-sample difference, and 10 pairs were selected without gene reuse, excluding pairs whose Spearman correlation with a retained pair had an absolute value ≥0.75 <!-- N-37 --> <!-- N-74, N-75 -->. For each sample, the reversal score counted, from 0 to 10, the pairs whose expression difference lay outside the dead zone and had the sign opposite to the dose-zero construction reference. The classification cutoff was set at the maximum reversal score among dose-zero construction samples, and scores above that cutoff were classified as positive.

<!-- C-08 の手法側 -->
In the dose-zero and High-AS construction groups, the binary readout characterized separation between the groups used to construct the panel. We then applied the classification rule unchanged to Low-AS and Mid-AS RET fusion-positive PTCs and summarized graded scores across all four bands, with the panel applied without refitting. Within the two intermediate bands, a one-sided Monte Carlo Brunner–Munzel comparison of Mid-AS with Low-AS <!-- N-75 --> assessed whether scores tended to be higher in Mid-AS than in Low-AS. <!-- C-09 --> We also computed PC-OD and relative-purity diagnostics for the intermediate bands (Supplementary Methods) <!-- C-12 --> but did not use them to exclude cases or modify the panel.

### External gene-list comparison

Two classes of published radiation-associated gene lists were cross-referenced against each contrast's q<0.10 list. The first class comprised qRT-PCR-validated lists from Abend et al. (2012, 2013) and Dom et al. (2012), with CLIP2 as a single-gene reference. The second comprised multivariate transcriptomic signatures: the shared-tissue, normal-tissue, and tumor-tissue signatures from post-Chornobyl cases (Ory et al. 2026), the 106-gene tumor signature that discriminated post-Chornobyl from sporadic PTCs in a reanalysis of an earlier series (Ugolin et al. 2011), and the 322-probe tumor signature that discriminated radiation-induced from sporadic thyroid tumors after external radiotherapy (Ory et al. 2011). <!-- C-13 --> We kept the two classes distinct and reported membership counts without an enrichment test <!-- N-53 -->.

### Reproducibility and AI use

Analyses were run in a date-pinned container with R 4.5.3 (R Core Team 2026) on Ubuntu 24.04 <!-- N-02 --> and reference BLAS/LAPACK 3.12.0 <!-- N-03 -->. Table S2 lists package versions. The analysis code, versioned inputs, and container recipe accompany the paper. Verification of the in-house implementations against reference implementations and an independent replication of the full pipeline on a second machine are described in Supplementary Methods.

The authors designed the study, specified the analysis and interpretation framework, and made the final scientific and editorial decisions. Under author supervision, Claude Code (Anthropic; Claude Fable 5 and 5.1) and Codex CLI (OpenAI; GPT-5.6 Sol and GPT-6 Astra) assisted with code development (principally debugging and refactoring), consistency checks, and drafting and editing manuscript text. The authors reviewed and approved all AI-assisted material retained in the reported analysis and manuscript and take full responsibility for the work. The scope and verification procedures are described in Supplementary Methods.

## Results

<!-- 2026-08-21 C2 ポート: gpt_review 試案 Results を移植(タグ移植込み。N-87/N-88 は新設 draft 行 — C3 で照合)。写像検査同期は C4 -->
### Cohort

<!-- C-10 --> Of the 440 <!-- N-11 --> REBC-THYR cases, 63 <!-- N-08 --> met the main cohort definition (Fig. 1; Table 2): the RET fusion-positive stratum contained 12 dose-zero and 15 High-AS cases (R_Sporadic and R_High), the BRAF V600E-positive stratum 27 and nine (B_Sporadic and B_High) <!-- N-09 -->. Ten of 12 dose-zero and 11 of 15 High-AS RET fusion-positive cases were female, as were 23 of 27 and four of nine BRAF V600E-positive cases <!-- N-12, N-13 -->. Among dose-zero and High-AS RET fusion-positive cases, the fusion was CCDC6::RET in six and seven, NCOA4::RET in two and four, and with another partner in four and four, respectively <!-- N-12 -->. Most attrition occurred at driver classification and at restriction to the extreme bands (Fig. 1). PC-OD flagged no tumor or normal sample in the RET fusion-positive stratum. Among the four target groups it flagged one sample, a High-AS BRAF V600E-positive tumor, and that case was excluded from the main cohort <!-- N-87, N-08 -->. The separate REO application set comprised 17 Low-AS and 19 Mid-AS RET fusion-positive PTCs <!-- N-10 -->.

High-AS cases were older at surgery than dose-zero cases by a Hodges–Lehmann estimate of 2.5 years in the RET stratum (95% confidence interval, −1.0 to 6.0) and 8.0 years in the BRAF stratum (3.0 to 12.0) <!-- N-64, N-65, C-15 -->. Sex, fusion-partner, and relative purity distributions are summarized in Table 2 and its footnotes.

### Gene-level and contrast-level results

In the primary RET-tumor contrast, 1,765 of 15,621 tested genes met the Storey q<0.10 rule <!-- N-16, N-15 -->. Of these, 971 were expressed at higher and 794 at lower levels in the High-AS group <!-- N-17, C-02 -->. The Higher Criticism omnibus p-value was 0.0112 <!-- N-20 --> (Fig. 2; Table 3; Fig. S1; Supplementary Data 1).

In the three secondary contrasts, RET-normal had no gene at q<0.10 and an HC p-value of 0.3199 <!-- N-16, N-20, C-03 -->, and BRAF-tumor likewise had no gene at q<0.10, with an HC p-value of 0.1815 <!-- N-16, N-20, C-04 -->. <!-- C-16 --> BRAF-normal had one gene at q<0.10, the X-linked BHLHB9 (relative effect 0.967, q=0.013) <!-- N-22 -->, and its contrast-level HC p-value was 0.0773 <!-- N-20 --> and the descriptive maximum-statistic p-value 0.0125 <!-- N-23 -->. Its π0 estimate of 0.727 was lower than the 0.955 in RET-normal and the 0.943 in BRAF-tumor <!-- N-19 -->. Of the 1,765 RET-tumor genes, 57 were X-linked and none Y-linked, compared with 572 X-linked and two Y-linked genes among the 15,621 tested (Supplementary Results) <!-- N-82 -->.

### Gene sets and list annotation

No set met q<0.10 in any of the 16 contrast-by-collection cells <!-- N-27, C-05 -->. The smallest adjusted value was 0.114 in the BRAF-tumor radiation-curated cell <!-- N-28 --> (Table S3; Supplementary Data 2). In the complete-null assessment, at least one set was selected in 102 of 1,600 held-out contrast-by-collection pseudo-observations (descriptive pooled rate, 0.064) <!-- N-24, C-06 -->, while BRAF-normal/Hallmark had a cell-specific rate of 0.18 (descriptive 95% interval, 0.110–0.269 <!-- N-25 -->; Fig. S2; Table S4). In the single positive-control scenario, the spike-in set ranked first among 50 Hallmark sets (q=0.0101) <!-- N-31 -->, and no other set met q<0.10 <!-- N-32 --> (Supplementary Results).

Under the separate gene-sampling reference, genes expressed at lower levels in High-AS RET fusion-positive PTCs were over-represented in proliferation, cell-cycle, and DNA-repair annotations, including E2F targets (46/199), G2M checkpoint (41/198) <!-- N-61 -->, and Reactome DNA repair (68/322) <!-- N-62 -->. In the radiation-curated family, seven sets met the rule, led by two Zhou cell-cycle sets of the ionizing-radiation response (46/126 and 30/84) <!-- N-60 -->. Among genes expressed at higher levels in High-AS RET fusion-positive PTCs, no set in any family met the over-representation threshold (Table S5) <!-- N-59, C-14 -->.

### Concordance between the RET and BRAF contrasts within each tissue

The RET-normal and BRAF-normal contrast profiles had Spearman rho=0.376 across 15,459 shared genes, within the shuffled reference interval of −0.46 to 0.46 (two-sided p=0.1199) <!-- N-34, C-07 -->.

The corresponding tumor-profile comparison yielded rho=0.459 across 15,560 genes, outside its shuffled interval of −0.39 to 0.39 (p=0.0197; Table S6) <!-- N-33, C-17 -->.

### REO construction and intermediate-band application

<!-- C-08 --> Of 57,694 candidate pairs, 153 met the selection rules <!-- N-35 --> and 10 formed the final panel <!-- N-37 --> (Table S7; Supplementary Results). The maximum reversal score among dose-zero construction samples was 2, and 13 of 15 High-AS RET fusion-positive cases scored above this cutoff <!-- N-38 -->.

Applied unchanged, the construction-derived cutoff classified the Low-AS cases as nine negative and eight positive and the Mid-AS cases as eight negative and 11 positive <!-- N-42 -->. Graded-score medians were 0, 1, 4, and 6 across the dose-zero, Low-AS, Mid-AS, and High-AS bands (Fig. 3) <!-- N-41 -->. The Mid–Low comparison yielded a relative effect θ=0.616 (Mid-AS relative to Low-AS) and one-sided p=0.1127 <!-- N-40 -->.

No intermediate-band tumor was flagged by the ancillary PC-OD screen <!-- N-43, C-09 -->. Relative purity could be estimated for 15 of 17 Low-AS and 16 of 19 Mid-AS cases <!-- N-44, N-10 -->. The other five cases had a tumor but no matched normal sample. Among the 31 tumors with purity scores, score and purity correlated at rho=0.538 <!-- N-45 -->. Band–score rho was 0.142, and the partial Spearman coefficient after adjustment for purity rank was 0.146 (descriptive one-sided permutation p=0.2162) <!-- N-46, N-48, C-12 -->. Band and purity correlated at rho=0.036 (Supplementary Results) <!-- N-48 -->.

### External gene lists

Among the qRT-PCR-validated lists, the Abend et al. (2013) tumor-tissue list and CLIP2, the two tissue-matched references for RET-tumor, had no member in the RET-tumor list, and neither had the Abend et al. (2012) list. Of the 20 cells formed by the five lists and the four contrasts, 18 were empty <!-- N-53 -->. The two overlaps were cross-tissue: S100A10 from the Dom et al. (2012) normal-tissue list and UBA3 from the Abend et al. (2013) normal-tissue list were present in the RET-tumor list, both with expression lower in High-AS cases <!-- N-54, N-104, C-13 -->. S100A10 had been reported higher in exposed normal tissue. UBA3 had been reported with a non-monotone categorical dose association in normal tissue whose shape the authors did not characterize, so its direction could not be compared <!-- N-104 -->. Among the multivariate signatures of Ory et al. (2026), five of 39 tested shared-tissue genes and three of 46 tested tumor-tissue genes were present in the RET-tumor list <!-- N-85 -->. Their normal-tissue signature had no overlap with either normal-tissue contrast, but three of its 40 tested genes were present in RET-tumor. Of the 90 tested genes of the post-Chornobyl signature of Ugolin et al. (2011), 10 were present in the RET-tumor list, four of them in the reported direction <!-- N-101, N-102 -->. Of the 265 tested genes of the post-radiotherapy signature (Ory et al. 2011), 23 were present in the RET-tumor list, six of them in the direction reported for radiation-induced relative to sporadic tumors <!-- N-99, N-100 -->. No signature from any of the three studies overlapped a BRAF q<0.10 list (Supplementary Results) <!-- N-85, N-99, N-101 -->.

## Discussion

<!-- 2026-08-21 C2 ポート: gpt_review 試案 Discussion を移植(C1-c C-11 復元・凍結文再設置・タグ移植)。★アンカー同期は C4 -->
In this study, High-AS and dose-zero RET fusion-positive PTCs differed extensively in expression. By comparing cases within a driver stratum, we prevented the overall composition of RET and BRAF drivers from defining the primary contrast and directly tested whether expression structure associated with the AS bands remained after driver stratification. The Higher Criticism p-value of 0.0112 <!-- N-20 -->, which summarizes an excess of moderately small p-values across genes, provided contrast-level evidence of structure aligned with the band labels, while the Storey q<0.10 rule identified 1,765 genes for reporting <!-- N-16 -->. Together, these results provide positive evidence for the working hypothesis within the design fixed in Table 1 <!-- C-01 -->. However, it remains unclear whether the difference is specific to radiation or reflects age, birth cohort, or other factors aligned with the AS bands. Throughout the Discussion, interpretations drawn from Table 1 (Methods) were assigned to outcome patterns before the finalized analysis produced the results reported here; interpretations beyond those in Table 1 are presented as hypothesis-generating.

The three secondary contrasts address the scope of this difference. In RET-normal, there was no detectable evidence for an accompanying expression difference in matched normal tissue, which leaves the proposed long-term normal-tissue signal unresolved in this cohort <!-- C-03 -->. In BRAF-tumor, for which Table 1 specified no direction, the null result at these group sizes establishes neither equivalence between the dose-zero and High-AS groups nor specificity of the RET-tumor result <!-- C-04 -->. In BRAF-normal, the single gene at q<0.10 and the weaker contrast-level statistic diverged against a background of pronounced age and sex imbalance, and this was the only contrast with a calibration cell whose descriptive complete-null interval lay above the nominal level <!-- N-25 -->. This result therefore remained unresolved and has no bearing on the primary finding <!-- C-16 -->.

Annotation of the q<0.10 gene list describes the composition of the list already identified and links the RET-tumor result, at the hypothesis-generating level, to proliferation and DNA repair. Under the gene-sampling reference used for over-representation analysis, genes expressed at lower levels in High-AS tumors were concentrated in proliferation, cell-cycle, and DNA-repair annotations, and the radiation-curated sets meeting the rule belonged to the same theme and so did not provide independent confirmation. Under the calibrated label-permutation analysis, no set met q<0.10. This divergence is the pattern expected for a co-expressed program (Goeman and Bühlmann 2007): gene sampling treats genes as independent units, so a proliferation program within a 1,765-gene list yields small p-values, whereas label permutation preserves gene dependence and asks whether the set-level shift exceeds what relabeling of the 27 cases produces. The over-representation result therefore points to this theme as a direction for further mechanistic studies rather than as experiment-level evidence <!-- C-14 -->.

Concordance between the RET and BRAF contrasts within each tissue asks whether the AS-band-associated profile is shared across driver strata. In normal tissue, concordance between the RET-normal and BRAF-normal profiles did not exceed its shuffled reference, but the RET-normal contrast itself had no gene at q<0.10 and no contrast-level evidence, so the absence of detectable sharing cannot be distinguished from a shared pattern below the detection limit <!-- C-07 -->. In tumor tissue, the RET-tumor and BRAF-tumor profiles were more concordant than their shuffled reference. As a descriptive, hypothesis-generating observation, this tumor concordance is the pattern expected if part of the AS-band-associated expression structure in RET fusion-positive tumors extends to BRAF V600E-positive tumors, and it is equally compatible with covariates aligned with the AS bands in both strata, such as tumor purity or age. The tumor concordance therefore motivates the hypothesis that a radiation-related expression change exists independently of the driver oncogene <!-- C-17 -->.

The REO panel gave each case a score derived from the group-level RET-tumor contrast. Two construction choices followed from the mixture interpretation set out in the Background. Because the AS bands were designed to enrich tumors arising through a radiation-associated initiating event without isolating them, we required reversal in more than half but not all High-AS construction cases for pair selection. Because quality-based exclusion within the intermediate bands could have altered their composition, we applied the panel to every eligible case in those bands without the construction-set outlier and purity filters. The graded median scores across the four bands were in the order that the mixture interpretation anticipates, and the specified Mid–Low comparison, with substantial uncertainty, was in the same direction <!-- C-08 -->. For the same reason given in the Background, each band within the RET stratum may contain a larger radiogenic proportion than its AS range implies. This interpretation is hypothesis-generating, and the ordered scores are consistent with it. They establish no linear dose response and validate no classifier, but they retain information that treating the intermediate-band test in isolation as negative would discard. The correlation between REO score and relative purity further indicates that the panel's biological basis remains to be resolved. Because the score depends only on the within-sample ordering of 10 gene pairs <!-- N-37 -->, it can be measured by targeted assays without genome-wide profiling (Geman et al. 2004), so validation in independent samples does not depend on reproducing the present sequencing platform.

Overlap with published gene lists was limited, and because the q<0.10 lists of the three secondary contrasts were empty or contained a single gene, the comparison was in effect with the RET-tumor list <!-- N-16 -->. The tissue-matched qRT-PCR references for RET-tumor had no member in it, S100A10 was opposite in direction to its normal-tissue report, and the UBA3 association reported by Abend et al. (2013) was not monotonic in dose and so gave no direction for comparison <!-- N-53, N-54, N-104 -->. The two multivariate tumor signatures differ in origin: the post-Chornobyl signature of Ugolin et al. (2011) came from a reanalysis of an earlier post-Chornobyl series, whereas the post-radiotherapy signature of Ory et al. (2011) was derived from tumors arising after external irradiation at high dose rates and was later applied in part to a post-Chornobyl series (Ory et al. 2013). Neither overlapped the RET-tumor list beyond what its size would produce, and the directions agreed in a minority of genes <!-- N-99, N-100, N-101, N-102 -->. Two features of the earlier designs may contribute. Both signatures came from contrasts of exposed with sporadic tumors without driver stratification, and in the learning set of Ory et al. (2011) the groups differed in driver composition, so part of each signature may reflect composition rather than exposure. In the series reanalyzed by Ugolin et al. (2011), the groups also differed in age and ethnicity: the exposed cases were young Ukrainian patients and the comparison cases older French patients. Detours et al. (2007) noted that confounding by age or ethnicity was not formally ruled out, and Ugolin et al. (2011) balanced age and sex, but not origin, within their learning set. In the present cohort the age difference between exposed and unexposed cases is in the opposite direction, because the dose-zero cases were born after the accident and were the younger group (Morton et al. 2021), so a gene whose expression varies with age or birth cohort could show opposite signs in the two studies. Dilution of exposed groups by sporadic tumors, which the mixture interpretation anticipates, would weaken a shared signal but would not reverse its sign, so dilution does not account for the disagreements. The published values of both signatures are mean log-expression differences computed in the validation space of the classification method rather than fold changes on the original scale, which further limits a direction comparison. Without an enrichment test, these memberships provide context only and do not establish independent support, validation, or a universal radiation signature <!-- C-13 -->.

The limitations concern the design, the exposure index, and the sample size, and they mainly affect attribution of the observed within-driver association. First, conditioning on the driver removes differences between drivers, but it does not show whether the expression structure found within the RET stratum is specific to RET fusion-positive tumors or shared by exposed tumors with other drivers. The BRAF-tumor contrast showed no detectable structure and the cross-driver concordance was hypothesis-generating, so this study cannot distinguish the two interpretations, and every claim in this paper holds under either <!-- C-11 -->. Resolving this question depends on other exposed cohorts with driver and dose data, because tumors arising decades after the exposure cannot be generated in an experimental system. Second, AS orders cases by estimated attributability and does not identify individual tumors as radiogenic. Its scale rests on a transfer of the IREP external-photon thyroid model to internal iodine-131 exposure under an acute exposure-rate convention that did not reconstruct the time course of iodine-131 dose, and the AS bands carry covariate structure that the data cannot remove. Age at exposure enters AS by design, and dose-zero and exposed cases differ structurally in birth cohort and latency, so the available data cannot separate these features from direct age- or cohort-related expression. Relative purity, fusion partner, and the unavailable batch metadata may likewise align with the band labels. We fitted no model adjusting jointly for these variables, because age enters AS by construction and conditioning on it would replace the marginal band contrast with an age-conditional estimand (Supplementary Methods). Third, group sizes were small (9–27 cases per group) <!-- N-09, N-10 -->. The exact and permutation-based procedures remained valid at these sizes, but power was limited, so the null results in the secondary contrasts remain compatible with differences that the study could not detect. The REO panel was selected from 27 construction cases <!-- N-10 -->. Although we did not tune the panel, a panel selected from so few cases may be specific to this sample. No independent-cohort or wet-laboratory replication was available, and validation in cohorts that lack exposed cases could test only the panel's specificity.

In conclusion, High-AS and dose-zero RET fusion-positive PTCs differed extensively in expression after driver stratification. The functional annotation, cross-driver tumor-profile concordance, and ordered REO profile provide additional observations consistent with the working hypothesis that AS-associated expression structure remains detectable within a driver-defined tumor class. All three were derived from the same RET-tumor contrast and therefore add no evidence independent of it. The present design cannot determine whether this expression structure is specific to radiation or is partly attributable to age, birth cohort, purity, or technical factors. It nevertheless provides a reproducible framework in which the biological and cross-driver implications of this association can be tested in future cohorts.

## Figure legends and table captions

<!-- 2026-08-21 C2 ポート: 試案キャプションへ更新。表1ブロック(旧・表2ブロック、2026-08-26 に引用順へ番号交換)(キャプション・行・規則)は正本の批准版を維持 — 位置づけ列・パターン規則は不可侵。対訳は廃止(f 保留方針) -->

**Figure 1 | Case flow.** From 440 <!-- N-11 --> REBC-THYR cases to the 63-case <!-- N-08 --> main cohort and 36-case <!-- N-10 --> REO application set, stratified by driver (RET fusion or BRAF V600E).

**Figure 2 | Gene-level evidence by contrast.** Signed Brunner–Munzel effect, defined as 2θ−1, where θ=Pr(X<Y)+0.5Pr(X=Y), against −log10 permutation p. Genes meeting Storey q<0.10 are colored by direction.

**Figure 3 | REO reversal scores by AS band.** Scores from 0 to 10 for the 10-pair panel <!-- N-37 -->. The panel was constructed from the dose-zero and High-AS groups and applied without refitting to the Low-AS and Mid-AS bands.

**Table 1 | Prespecified interpretation of the four exposure contrasts, fixed before the finalized analysis produced the results reported here.** Each contrast has a role and its basis; the interpretation rules fixed together with the table are given as footnotes.

<!-- 表1(旧表2)本体+脚注規則(正本。Methods: Analysis contrasts から 2026-08-18 移設 — 投稿整形時に一方向切り出し、編集可能な第2コピーは作らない。規則の (a)-(f) 脚注化は研究者の書き直し待ち) -->
| Contrast | Role | Basis |
| --- | --- | --- |
| R_Tumor | Tests the working hypothesis; signal expected | primary expectation of the working hypothesis (Background); the driver-conditioned gap left by Morton et al. (2021) |
| R_Normal | Tests the working hypothesis; signal possible | long-term molecular memory of exposure — dose-dependent expression within exposed normal tissue (Abend et al. 2013); signatures distinguishing exposed from unexposed normal tissue (Ory et al. 2026) |
| B_Normal | Tests whether an exposure-related normal-tissue change is shared across driver strata | counterpart of R_Normal; the shared signatures of Ory et al. (2026) were reported without driver stratification |
| B_Tumor | No direction specified; a signal cannot be ruled out | point-mutation drivers become less likely as dose rises in the source cohort (Morton et al. 2021) |

Interpretation rules, fixed with the table: a signal in R_Tumor is interpreted as agreement with the prior expectation; a tumor signal confined to the BRAF stratum would contradict the hypothesis and is interpreted with non-etiological explanations first; concordant signals in both normal contrasts are interpreted as consistent with an exposure-related change shared across the gland, and discordant signals as suggesting that the type of change depends on the driver; a normal-tissue signal confined to the BRAF stratum is interpreted first as confounding, and one confined to the RET stratum as unable to distinguish attenuated sharing from a driver-dependent type; and an all-null outcome is reported as a null result bounded by the study's detection limits. These rules cover every possible outcome pattern.

**Table 2 | Case characteristics.** Group size with the size of the stratum-by-band pool from which each group was drawn (cases with a tumor–normal pair in parentheses) <!-- N-95 -->, sex, age at surgery and age at exposure, and RET fusion-partner composition; footnotes give the age-difference estimates with confidence intervals and the group medians of relative tumor purity.

*Table 2 footnote.* Age-at-surgery differences compare High-AS with dose-zero cases. In the RET stratum, the Hodges–Lehmann difference was +2.5 years (95% percentile bootstrap confidence interval, −1.0 to 6.0) and the Brunner–Munzel relative effect was θ=0.625 (0.400 to 0.828) <!-- N-64 -->. In the BRAF stratum, the corresponding estimates were +8.0 years (3.0 to 12.0) and θ=0.850 (0.681 to 0.973) <!-- N-65 -->. Here, θ=Pr(X<Y)+0.5Pr(X=Y), with X denoting dose-zero and Y High-AS. Intervals were obtained from 9,999 resamples drawn separately within each group (seed 19450809) <!-- N-63 -->; no p-values were calculated. Relative tumor purity is a within-cohort relative score (cohort maximum = 1) estimated separately per driver cohort and comparable within a stratum only; group medians were 0.783 (RET dose-zero), 0.822 (RET High-AS), 0.836 (BRAF dose-zero), and 0.922 (BRAF High-AS) <!-- N-88 -->.

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

遺伝子セットの収録元(Supplementary Methods: Gene-set inference、2026-09-07 追加。MSigDB・msigdbr・Hallmark の既存引用を保持):

- Milacic M, Beavers D, Conley P, Gong C, Gillespie M, Griss J, et al. Nucleic Acids Res 2024;52:D672-D678. doi:10.1093/nar/gkad1025 — The Reactome Pathway Knowledgebase 2024 (Supplementary Methods: Gene-set inference、収録元の原典引用)。
- Agrawal A, Balcı H, Hanspers K, Coort SL, Martens M, Slenter DN, et al. Nucleic Acids Res 2024;52:D679-D689. doi:10.1093/nar/gkad960 — WikiPathways 2024: next generation pathway database (Supplementary Methods: Gene-set inference、収録元の原典引用)。
- Kanehisa M, Furumichi M, Sato Y, Matsuura Y, Ishiguro-Watanabe M. Nucleic Acids Res 2025;53:D672-D677. doi:10.1093/nar/gkae909 — KEGG: biological systems database as a model of the real world (Supplementary Methods: Gene-set inference、収録元の原典引用)。
- Nishimura D. Biotech Softw Internet Rep 2001;2:117-120. doi:10.1089/152791601750294344 — BioCarta (Supplementary Methods: Gene-set inference、収録元の原典引用)。
- Schaefer CF, Anthony K, Krupa S, Buchoff J, Day M, Hannay T, et al. Nucleic Acids Res 2009;37:D674-D679. doi:10.1093/nar/gkn653 — PID: the Pathway Interaction Database (Supplementary Methods: Gene-set inference、収録元の原典引用)。
- Ashburner M, Ball CA, Blake JA, Botstein D, Butler H, Cherry JM, et al. Nat Genet 2000;25:25-29. doi:10.1038/75556 — Gene ontology: tool for the unification of biology (Supplementary Methods: Gene-set inference、収録元の原典引用)。
- The Gene Ontology Consortium. Nucleic Acids Res 2026;54:D1779-D1792. doi:10.1093/nar/gkaf1292 — The Gene Ontology knowledgebase in 2026 (Supplementary Methods: Gene-set inference、収録元の原典引用)。

放射線関連 C2:CGP の原典(Supplementary Methods: Gene-set inference / Table S8、2026-09-07 追加。報告対象24セットの出典11報であり、本研究の独立した検証研究としては扱わない):

- Ghandhi SA, Yaghoubian B, Amundson SA. BMC Med Genomics 2008;1:63. doi:10.1186/1755-8794-1-63 — Global gene expression analyses of bystander and alpha particle irradiated normal human lung fibroblasts: synchronous and differential responses.
- Kyng KJ, May A, Stevnsner T, Becker KG, Kølvrå S, Bohr VA. Oncogene 2005;24:5026-5042. doi:10.1038/sj.onc.1208692 — Gene expression responses to DNA damage are altered in human aging and in Werner Syndrome.
- Macaeva E, Saeys Y, Tabury K, Janssen A, Michaux A, Benotmane MA, et al. Sci Rep 2016;6:19251. doi:10.1038/srep19251 — Radiation-induced alternative transcription and splicing events and their applicability to practical biodosimetry.
- Monnier Y, Farmer P, Bieler G, Imaizumi N, Sengstag T, Alghisi GC, et al. Cancer Res 2008;68:7323-7331. doi:10.1158/0008-5472.CAN-08-0841 — CYR61 and alphaVbeta5 integrin cooperate to promote invasion and metastasis of tumors growing in preirradiated stroma.
- Quintens R, Verreet T, Janssen A, Neefs M, Leysen L, Michaux A, et al. Biol Open 2015;4:331-344. doi:10.1242/bio.20149969 — Identification of novel radiation-induced p53-dependent transcripts extensively regulated during mouse brain development.
- Rashi-Elkeles S, Elkon R, Weizman N, Linhart C, Amariglio N, Sternberg G, et al. Oncogene 2006;25:1584-1592. doi:10.1038/sj.onc.1209189 — Parallel induction of ATM-dependent pro- and antiapoptotic signals in response to ionizing radiation in murine lymphoid tissue.
- Smirnov DA, Brady L, Halasa K, Morley M, Solomon S, Cheung VG. Genome Res 2012;22:332-339. doi:10.1101/gr.122044.111 — Genetic variation in radiation-induced cell death.
- Tsai MH, Chen X, Chandramouli GV, Chen Y, Yan H, Zhao S, et al. Oncogene 2006;25:622-632. doi:10.1038/sj.onc.1209082 — Transcriptional responses to ionizing radiation reveal that p53R2 protects against radiation-induced mutagenesis in human lymphoblastoid cells.
- Tsai MH, Cook JA, Chandramouli GV, DeGraff W, Yan H, Zhao S, et al. Cancer Res 2007;67:3845-3852. doi:10.1158/0008-5472.CAN-06-4250 — Gene expression profiling of breast, prostate, and glioma cells following single versus fractionated doses of radiation.
- Warters RL, Packard AT, Kramer GF, Gaffney DK, Moos PJ. Radiat Res 2009;172:82-95. doi:10.1667/RR1677.1 — Differential gene expression in primary human skin keratinocytes and fibroblasts in response to ionizing radiation.
- Zhou T, Chou J, Mullen TE, Elkon R, Zhou Y, Simpson DA, et al. Cell Cycle 2007;6:972-981. doi:10.4161/cc.6.8.4106 — Identification of primary transcriptional regulation of cell cycle-regulated genes upon DNA damage.

前処理・正規化・純度(いずれも Methods: Data sources / Normalization / QC):

- Law CW, Chen Y, Shi W, Smyth GK. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol 2014;15:R29. doi:10.1186/gb-2014-15-2-r29 — 純度推定に用いる limma-voom の個別手法引用(Supplementary Methods: Quality control and analysis cohorts)。
- Phipson B, Lee S, Majewski IJ, Alexander WS, Smyth GK. Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression. Ann Appl Stat 2016;10:946–963. doi:10.1214/16-AOAS920 — 同工程の頑健な経験ベイズ推定。既存の limma 全体の Ritchie et al. 2015 と、置換 p 値の Phipson and Smyth 2010 はそれぞれ別に保持する。
- Dobin A, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics 2013;29:15–21. doi:10.1093/bioinformatics/bts635 — STAR counts の来歴(GDC パイプライン言及の形でも可)
- Heath AP, Ferretti V, Agrawal S, et al. The NCI Genomic Data Commons. Nat Genet 2021;53:257–262. doi:10.1038/s41588-021-00791-5 — GDC の公式引用案内が指定するリソース論文(Methods: Data sources / Supplementary Methods: Data sources and expression matrix)。REBC-THYR の原研究 Morton et al. 2021 とアクセス用パッケージ Morgan and Davis 2025 の引用は別に保持する。
- Mudge JM, Carbonell-Sala S, Diekhans M, Martinez JG, Hunt T, Jungreis I, et al. GENCODE 2025: reference gene annotation for human and mouse. Nucleic Acids Res 2025;53:D966–D975. doi:10.1093/nar/gkae1078 — GENCODE 自体の公式代表引用(Supplementary Methods: Data sources and expression matrix)。使用データは GENCODE v36 のままであり、2025 年の注釈版へ更新したという意味ではない。
- Grossman RL, et al. Toward a shared vision for cancer genomic data. N Engl J Med 2016;375:1109–1112. doi:10.1056/NEJMp1607591 — GDC の旧候補(本文不使用。公式指定の Heath et al. 2021 を採用)
- Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 2010;26:139–140. doi:10.1093/bioinformatics/btp616— filterByExpr / DGEList / CPM〔使用パッケージ〕
- Chen Y, Lun ATL, Smyth GK. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. F1000Res 2016;5:1438. doi:10.12688/f1000research.8987.2(Crossref 照合済み 2026-08-15)— filterByExpr の規則の出典(Methods: QC 初出 — 本文反映済み。以降の filterByExpr 出現は同規則)
- Kadota K, Nishiyama T, Shimizu K. A normalization strategy for comparing tag count data. Algorithms Mol Biol 2012;7:5. doi:10.1186/1748-7188-7-5 — DEGES
- Sun J, Nishiyama T, Shimizu K, Kadota K. TCC: an R package for comparing tag count data with robust normalization strategies. BMC Bioinformatics 2013;14:219. doi:10.1186/1471-2105-14-219 — iDEGES(3反復)の系譜〔手法引用/自作実装 — TCC パッケージ不使用、スクリーンは BM に置換(本文明記済み)〕
- Feng Y, Li LM. MUREN: a robust and multi-reference approach of RNA-seq transcript normalization. BMC Bioinformatics 2021;22:386. doi:10.1186/s12859-021-04288-0 — MUREN(※書誌照合で年・巻を訂正: 2020;21 → 2021;22)〔手法引用/自作実装〕(Methods: QC 初出で引用 — 本文反映済み)
- Shen Q, et al. contamDE: differential expression analysis of RNA-seq data for contaminated tumor samples. Bioinformatics 2016;32:705–712. doi:10.1093/bioinformatics/btv657 — contamDE 原法(比率スコアの親手法; 同論文の反復擬似尤度推定量は不使用)〔手法引用〕(本文反映済み)
- Ji Y, Yu C, Zhang H. contamDE-lm: linear model-based differential gene expression analysis using next-generation RNA-seq data from contaminated tumor samples. Bioinformatics 2020;36:2492–2499. doi:10.1093/bioinformatics/btaa006(PubMed 照合済み 2026-08-16、pmid 31917401)— 実装した比率推定量(最大1正規化)の出典〔手法引用/自作実装 — lib/purity_contamde.R〕(Methods: QC — 本文反映済み)

被曝指標・その他:

- Kocher DC, et al. Interactive RadioEpidemiological Program (IREP): a web-based tool for estimating probability of causation/assigned share. Health Phys 2008;95:119–147. doi:10.1097/01.HP.0000291191.49583.f7 — NIH IREP(Methods: Exposure metric。使用版 5.7.3。全文照合 2026-08-16: 放射性ヨウ素への言及ゼロ — 内部被曝の正当化は単独では担えない、下記 Land・Sato・IREP 文書と分担)
- Land CE, Gilbert E, Smith JM, et al. Report of the NCI-CDC Working Group to Revise the 1985 NIH Radioepidemiological Tables. NIH Publication No. 03-5387. Bethesda, MD: National Institutes of Health; 2003. — REF 表(electrons E>15keV = 単一値 1.0; 放射性核種β線は平均エネルギーで区分に割り付け — 脚注 c)・DDREF 規則(chronic 全域/acute は D_L = 0.03–0.2 Gy 対数一様以上で 1)・甲状腺モデル来歴(全文 PDF 逐語照合 2026-08-16 完了、md5 は N-83 — 再確認の保留は解除)(Methods: Exposure metric)
- National Cancer Institute (n.d.). NIH-IREP: Guidance on Selection of Radiation Type. <https://radiationcalculators.cancer.gov/irep/model/help/help_totaldose/> (accessed 7 September 2026) — ヨウ素131の electrons E>15 keV 指定の直接の典拠(Supplementary Methods: Exposure metric: Assigned Share)。2026-09-07 に Web 引用整形の保留を解除。ページの発行年は不明(n.d.)であり、使用プログラム版 5.7.3 とは区別する。
- Sato T, Manabe K, Hamada N. Microdosimetric analysis confirms similar biological effectiveness of external exposure to gamma-rays and internal exposure to 137Cs, 134Cs, and 131I. PLoS One 2014;9:e99831. doi:10.1371/journal.pone.0099831(PubMed 照合済み 2026-08-16、pmid 24919099)— 内部ヨウ素131の RBE ≈ 1(保守的仮定下でも最大 1.04)の独立した物理的支持(Methods: Exposure metric)
- Ron E, Lubin JH, Shore RE, et al. Thyroid cancer after exposure to external radiation: a pooled analysis of seven studies. Radiat Res 1995;141:259–277.(PubMed 照合済み 2026-08-16、pmid 7871153、DOI なし)— IREP 甲状腺モデルの適合データ来歴(Methods: Exposure metric)
- Zurnadzhy L, et al. Clinicopathological implications of the BRAF V600E mutation in papillary thyroid carcinoma of Ukrainian patients exposed to the Chernobyl radiation in childhood: a study for 30 years after the accident. Front Med (Lausanne) 2022;9:882727. doi:10.3389/fmed.2022.882727(PubMed 照合済み 2026-08-16、pmid 35665338)— 同集団への先行 IREP 適用・acute 規約の来歴(**規約の先例としてのみ引用** — 共著者の先行研究のため等価性の根拠には使わない)(Methods: Exposure metric)
- Bogdanova T, et al. The relationship of the clinicopathological characteristics and treatment results of post-Chornobyl papillary thyroid microcarcinomas with the latency period and radiation exposure. Front Endocrinol (Lausanne) 2022;13:1078258. doi:10.3389/fendo.2022.1078258(PubMed 照合済み 2026-08-16、pmid 36589808)— 同上(本研究と同一の AS 出力行「expected value of ERR」を使用)(Methods: Exposure metric)
- Geman D, et al. Classifying gene expression profiles from pairwise mRNA comparisons. Stat Appl Genet Mol Biol 2004;3:Article19. doi:10.2202/1544-6115.1071 — 検体内順位比較(REO)の原型〔概念引用/パネル構築規則は自作〕(Methods: REO — 本文反映済み)
- Wang H, Sun Q, Zhao W, Qi L, et al. Individual-level analysis of differential expression of genes and pathways for personalized medicine. Bioinformatics 2015;31:62–68. doi:10.1093/bioinformatics/btu522(PubMed/Crossref 照合済み 2026-08-15、pmid 25165092)— 安定参照+逆転スキーム(RankComp)・正規化フリー性の系譜(REO の用語は同論文にない — 2026-08-16 照合で訂正、用語の典拠は Guan 2016)〔概念引用/パネル構築規則は自作〕(Methods: REO — 本文反映済み)
- Guan Q, Chen R, Yan H, Cai H, Guo Y, Li M, Li X, Tong M, Ao L, Li H, Hong G, Guo Z. Differential expression analysis for individual cancer samples based on robust within-sample relative gene expression orderings across multiple profiling platforms. Oncotarget 2016;7:68909–68920. doi:10.18632/oncotarget.11996(PubMed 照合済み 2026-08-16、pmid 27634898)— REO の用語の典拠(略号つき初出は Cai 2015, doi:10.18632/oncotarget.6260 — 引用は Guan の正準的定義のみ)(Methods: REO — 本文反映済み)
- Nikiforov YE, Rowland JM, Bove KE, Monforte-Munoz H, Fagin JA. Distinct pattern of ret oncogene rearrangements in morphological variants of radiation-induced and sporadic thyroid papillary carcinomas in children. Cancer Res 1997;57:1690–1694.(PubMed 照合済み 2026-08-15、pmid 9135009、DOI なし)— 放射線関連 PTC とパートナー特異性(RET/PTC3 優位)の原典・旧称の橋渡し(Methods: Candidate confounders — 本文反映済み)
- Rabes HM, Demidchik EP, Sidorow JD, et al. Pattern of radiation-induced RET and NTRK1 rearrangements in 191 post-Chernobyl papillary thyroid carcinomas: biological, phenotypic, and clinical implications. Clin Cancer Res 2000;6:1093–1103.(PubMed 照合済み 2026-08-15、pmid 10741739、DOI なし)— 同上・最大古典シリーズ(Methods: Candidate confounders / Intro — 本文反映済み)
- Ito T, et al. In vitro irradiation is able to cause RET oncogene rearrangement. Cancer Res 1993;53:2940–2943.(PubMed 照合済み 2026-08-15、pmid 8319199、DOI なし)— RET 再配列の実験的誘導の原典(本文未使用 — 2026-10-15 に Background から除去、Caudill 2005 が単独で担う)
- Neubert K, Brunner E. A studentized permutation test for the non-parametric Behrens–Fisher problem. Comput Stat Data Anal 2007;51:5192–5204. doi:10.1016/j.csda.2006.05.024 — studentized 置換 BM 検定の手法帰属(Methods: Normalization and gene-level inference — 本文反映済み 2026-09-08。Supp では既出)
- Nikiforova MN, et al. Proximity of chromosomal loci that participate in radiation-induced rearrangements in human cells. Science 2000;290:138–141. doi:10.1126/science.290.5489.138(PubMed 照合済み 2026-08-15、pmid 11021799)— RET と融合相手の核内空間隣接(機構論)(Intro — 本文反映済み)
- Gandhi M, Medvedovic M, Stringer JR, Nikiforov YE. Interphase chromosome folding determines spatial proximity of genes participating in carcinogenic RET/PTC rearrangements. Oncogene 2006;25:2360–2366. doi:10.1038/sj.onc.1209268 (PMID 16331264)— RET と NCOA4・H4(= CCDC6)の核内近接性。Nikiforova et al. 2000 と併記し、NCOA4 を含む記述を支える(Background — 本文反映済み。2026-09-08、PubMed 書誌・要旨照合)。
- Caudill CM, et al. Dose-dependent generation of RET/PTC in human thyroid cells after in vitro exposure to gamma-radiation. J Clin Endocrinol Metab 2005;90:2364–2369. doi:10.1210/jc.2004-1811(PubMed 照合済み 2026-08-15、pmid 15671095)— 誘導の線量依存性(Intro — 本文反映済み)
- Hamatani K, et al. RET/PTC rearrangements preferentially occurred in papillary thyroid cancer among atomic bomb survivors exposed to high radiation dose. Cancer Res 2008;68:7176–7182. doi:10.1158/0008-5472.CAN-08-0293(PubMed 照合済み 2026-08-15、pmid 18757433)— 独立曝露設定(原爆・外部γ線)でのヒト線量関連(Intro — 本文反映済み)
- R Core Team (2026). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. <https://www.R-project.org/(citation()> 出力どおり 2026-08-16、rebc-r453:refblas 内 R 4.5.3 (2026-03-11) — 版数の書式は投稿誌指定に合わせ後工程)— (Methods: Software)
- その他の使用パッケージ(SummarizedExperiment・limma・GenomicDataCommons・rtracklayer・data.table・Matrix・Rcpp・MASS・parallel)は Supp の版表が網羅 —個別引用は原典論文を持つもののみ執筆時に判断〔使用パッケージ〕。事前収集済みの原典(PubMed 照合済み 2026-08-14): Bioconductor/SummarizedExperiment = Huber W, et al. Nat Methods 2015;12:115–121, doi:10.1038/nmeth.3252 / limma = Ritchie ME, et al. Nucleic Acids Res 2015;43:e47, doi:10.1093/nar/gkv007 / rtracklayer = Lawrence M, et al. Bioinformatics 2009;25:1841–1842, doi:10.1093/bioinformatics/btp328

### C. Intro / Discussion で必要になる見込み

外部アンカー(C-13。DOI は diagnostics/external_gene_anchors.csv の curation 由来):

- Abend M, et al. PLoS One 2012;7:e39103. doi:10.1371/journal.pone.0039103(PubMed 照合済み)— ペア差リスト(11 遺伝子)
- Abend M, et al. Iodine-131 dose-dependent gene expression. Br J Cancer 2013;109:2286–2294. doi:10.1038/bjc.2013.574(PubMed 照合済み)— 正常8/腫瘍6
- Dom G, et al. A gene expression signature distinguishes normal tissues of sporadic and radiation-induced papillary thyroid carcinomas. Br J Cancer 2012;107:994–1000. doi:10.1038/bjc.2012.302(PubMed 照合済み)— 正常7 = SERPINE1・DUSP1・TRIB1・S100A10・ANXA1・GNAL・RDH12(原文 Results と Fig. 2 legend で照合 2026-08-16、diagnostics/external_gene_anchors.csv と 1:1 一致 — 保留解除。clade E の番号は SERPINE1 と legend に明記)
- Hess J, et al. Gain of chromosome band 7q11 in papillary thyroid carcinomas of young patients is associated with exposure to low-dose irradiation. PNAS 2011;108:9595–9600. doi:10.1073/pnas.1017137108(PubMed 照合済み)— CLIP2 (追加検証: doi:10.1038/onc.2014.311)

正常組織・対照・仮説:

- Ory C, et al. Post-Chornobyl thyroid papillary carcinomas display distinct past exposure and radiation-associated carcinogenesis molecular signatures at low and high thyroid doses. Sci Rep 2026;16. doi:10.1038/s41598-026-53030-4(PubMed 照合済み)— 正常組織対照・driver 非層別の隙間(Intro 予告 / Discussion C-07)
- Ory C, et al. Gene expression signature discriminates sporadic from post-radiotherapy-induced thyroid tumors. Endocr Relat Cancer 2011;18:193–206. doi:10.1677/ERC-10-0205(PubMed 照合済み 2026-10-15)— Background P1 橋渡し文(トランスクリプトーム水準の根拠、甲状腺)
- Hadj-Hamou NS, et al. A transcriptome signature distinguished sporadic from postradiotherapy radiation-induced sarcomas. Carcinogenesis 2011;32:929–934. doi:10.1093/carcin/bgr064(PubMed 照合済み 2026-10-15)— Background P1 橋渡し文(他組織)(本文未使用 — 2026-10-15 に Background から除去、Ory 2011 が橋渡しを単独で担う)
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
