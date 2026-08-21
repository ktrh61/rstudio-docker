# 原稿統合ドラフト(draft_manuscript — 論文へ収束する単一ファイル)

状態: draft(2026-08-14 統合。研究者指示: Methods/Results の対構造の検査可能化と同期機構の廃止のため、旧4ドラフト — draft_methods / draft_results / draft_methods_results_d6 / draft_ora_annotation — を本ファイルへ吸収し削除。内容の変更なし、移設+Discussion 骨組み追加のみ。履歴は git)。規約: 英語本文が正本。インライン【訳】は C2 ポートで段階的に廃止(対訳は内容確定後に共著者向け別ファイルとして再生成 — 2026-08-21 研究者方針)。Supplementary Methods/Results は paper/supplementary_material.md へ分離(2026-08-21、BJC 投稿構造)。数値は N-ID タグ、主張段落は C-ID タグ(検査対象は英語本文)。図表番号は 2026-08-16 に本置き済み(主: 図1–3・表1–3 / 補足: 図S1–S2・表S1–S8・データ1。BJC の Article 規定 = 主6品目以内に整合)。各節の由来・決定は末尾の「起草メモ(統合)」。

---

<!-- タイトルページ素材(BJC: title page 別建て・語数枠外) -->
**Title:** Driver-conditioned transcriptomic differences across radiation-attributability bands in papillary thyroid carcinoma

## Abstract

### Background

Radiation dose is associated with fusion rather than point-mutation drivers in post-Chernobyl papillary thyroid carcinoma, while driver composition structures expression. We therefore asked whether expression differed between cases with high estimated radiation attribution and dose-zero cases within driver-defined groups.

### Methods

We calculated Assigned Share (AS) with NIH IREP for ordinal banding only. We stratified by tumor driver (RET fusion or BRAF V600E) and compared expression between High-AS and dose-zero cases in tumor and matched normal tissue; the RET fusion-positive tumor comparison was primary. We used permutation Brunner–Munzel tests, Higher Criticism, label-permutation gene-set inference, and a relative-expression-ordering (REO) panel constructed using extreme-band cases and applied unchanged to intermediate-band cases.

### Results

In RET fusion-positive PTCs, 1,765 of 15,621 genes met Storey q<0.10 <!-- N-16, N-15 -->, with omnibus p=0.0112 <!-- N-20 -->; other contrasts lacked concordant gene- and contrast-level evidence. No permutation gene set met q<0.10, whereas over-representation analysis highlighted cell-cycle and DNA-repair annotations among genes expressed less in High-AS cases. The REO score showed an ordered descriptive profile across bands (dose-zero/Low-AS/Mid-AS/High-AS medians, 0/1/4/6 <!-- N-41 -->), with the extremes as construction anchors.

### Conclusions

High-AS and dose-zero RET fusion-positive PTCs differed transcriptionally, a result compatible with the hypothesized within-driver AS-band association. This does not establish a radiation-specific signature, and the REO profile requires independent validation.

## Introduction

<!-- 2026-08-21 C2 ポート: gpt_review 試案 Intro を移植(C1-e 仮説スロット批准・凍結文言は新正準へ統合)。予告列の同期は C4。インライン対訳は廃止(f 保留方針) -->

Papillary thyroid carcinoma after the Chernobyl accident remains a major human model of radiation carcinogenesis. Its molecular hallmark was recognized early: childhood tumors arising after the accident frequently carry RET fusions, and their fusion-partner distribution, with NCOA4-RET predominating, differs from that of age-matched sporadic childhood tumors (Nikiforov et al. 1997; Rabes et al. 2000). The historical labels RET/PTC1 and RET/PTC3 correspond to CCDC6-RET and NCOA4-RET, respectively; we use the gene-partner forms throughout. Experimental irradiation can induce RET rearrangements in thyroid cells (Ito et al. 1993; Caudill et al. 2005), the participating loci lie close together in thyroid-cell nuclei (Nikiforova et al. 2000), and RET rearrangements have also been enriched among highly exposed atomic-bomb survivors (Hamatani et al. 2008). At cohort scale, REBC-THYR showed that increasing thyroid dose was associated with fusion rather than point-mutation drivers (Morton et al. 2021). Whether exposure is associated with tumor expression within a driver stratum remains unresolved.

Two features make that question difficult. First, tumor transcriptomes are strongly organized by driver, while driver prevalence itself varies with radiation dose (Morton et al. 2021). A pooled exposed-versus-unexposed analysis can therefore reproduce differences between fusion- and point-mutation-driven tumors rather than an exposure-associated difference within a comparable molecular background. Second, previous expression findings span different biological compartments and designs: dose-associated expression has been reported in tumor and normal tissue, including paired tumor-minus-normal analyses (Abend et al. 2012, 2013); CLIP2 was proposed as a tumor-tissue radiation marker (Hess et al. 2011); and exposed-case normal-tissue signatures have been reported without driver stratification (Dom et al. 2012; Ory et al. 2026). Most did not stratify expression comparisons by driver. The CLIP2 study matched RET-rearrangement status as far as possible but focused on a copy-number marker and CLIP2 rather than genome-wide expression. Whether earlier expression signals persist within driver strata is therefore unclear.

The source cohort's continuous dose-expression analyses did not identify a corresponding expression signal (Morton et al. 2021), motivating a different but related question. We hypothesized that comparing extreme bands of a case-level radiation-attributability index might reveal structure diluted by an unstratified continuous-dose model. This reflects a mixture model in which dose changes the proportion of tumors arising through a radiation-associated initiating event without requiring expression in every tumor to vary monotonically with dose. The index was the Assigned Share (AS) from the NIH Interactive RadioEpidemiological Program (Kocher et al. 2008), calculated from recorded thyroid dose and age under a specified input convention. We used AS only to order exposed cases and define bands within driver strata. The bands were intended to enrich contrasting groups; they did not identify individual tumors as radiogenic or sporadic, and AS was not treated as an observed etiology or calibrated individual probability.

Only the RET fusion-positive and BRAF V600E-positive strata contained enough dose-zero and exposed cases for driver-stratified analysis. Within each stratum, we compared the High-AS band with dose-zero cases separately in tumor tissue and matched normal tissue; these four comparisons were prespecified. The RET fusion-positive tumor comparison was primary, reflecting both the prior expectation of radiation-related tumor initiation in this stratum and the cohort-level association of increasing dose with fusion rather than point-mutation drivers. BRAF therefore addressed a distinct within-driver question, not a negative or specificity control. The three secondary comparisons asked whether a corresponding difference occurred in normal tissue, whether expression differed in BRAF V600E-positive PTCs, and whether normal-tissue patterns were concordant across driver strata. Their results were not combined into a study-wide significance criterion. Table 2 gives the interpretation assigned to each possible pattern before the finalized analysis produced the results reported here. Separately, we constructed a relative-expression-ordering (REO) panel from the two RET extreme bands and applied it unchanged to Low- and Mid-AS RET fusion-positive PTCs.

## Methods

<!-- 2026-08-21 C2 ポート: gpt_review 試案 Methods を移植(C1-a AI開示復元・C1-d 目盛り要約復帰・chronic監査文は不採用)。写像検査・★アンカー同期は C4 -->
### Data sources

Gene-level STAR count files for REBC-THYR were downloaded from the National Cancer Institute (NCI) Genomic Data Commons (GDC) and verified against the download manifest. The clinical source was Data S1 of Morton et al. (2021), containing 440 cases <!-- N-11 -->. GDC files were mapped to cases and biospecimens through the GDC API and assembled into a matrix of 58,448 nonzero genes and 906 samples <!-- N-14 -->. Strandedness, gene-length derivation, file-selection, query, and manifest-verification procedures are described in Supplementary Methods.

### Assigned Share

For each exposed case, we calculated the AS associated with the expected excess relative risk using the NIH IREP thyroid model, version 5.7.3, and the recorded thyroid dose and ages (Kocher et al. 2008). We used the IREP input route for internal iodine-131 exposure and a prespecified convention in which exposure rate was entered as acute. This setting defined the numerical AS scale used for band construction; it was not intended to reconstruct the temporal pattern of internal iodine-131 dose delivery. The thyroid risk model was fitted to external photon exposures; applying it to internal iodine-131 is therefore a model-transfer assumption. The exposure-rate setting enters the calculation only through the dose and dose-rate effectiveness factor, and the High band's minimum dose (188 mGy <!-- N-84 -->) sits at the top of the model's uncertain reference-dose range <!-- N-83 -->, so a chronic calculation would compress the AS scale over the High band essentially by one common factor and could not materially reorder cases (Supplementary Methods). Dose was entered as a point value, so dosimetric and model uncertainty did not propagate into band membership.

AS was used only as a convention-dependent ordinal design index. Dose-zero cases formed the reference group, labeled Sporadic in the analysis files. Exposed cases were classified as Low (0<AS<33.3), Mid (33.3≤AS<66.6), or High (AS≥66.6); no case lay on a boundary <!-- N-69 -->. These were operational boundaries on the specified calculation scale, not biological or etiological probability thresholds. The four extreme-band contrasts used High-AS versus dose-zero group membership, not the numeric AS value, and AS was not interpreted as a calibrated probability that an individual tumor was radiogenic. The full input convention, its provenance, and an audit of alternative exposure-rate coding are in Supplementary Methods.

### Quality control and cohorts

We represented the eligibility criteria as case-level flags and fixed them before the inferential analyses reported here. Principal-component-scores-based outlier detection (PC-OD) was applied separately within each observed group-by-tissue matrix after filtering genes with edgeR's filterByExpr, using unnormalized log counts per million with prior count 2 and iterative two-sided screening at α=0.05 (Chen et al. 2016; Chen et al. 2025; Nakayama et al. 2024). PC-OD preceded relative tumor-purity estimation so that an anomalous sample could not affect the purity fit. We estimated relative tumor purity jointly for dose-zero and High-AS cases within each driver cohort from matched tumor-normal counts, using the proportion component of contamDE-lm after MUREN normalization (Shen et al. 2016; Ji et al. 2020; Feng and Li 2021), thereby placing both groups on a common relative scale. Purity served as a within-cohort eligibility filter and diagnostic variable, not as an absolute cellular fraction.

The main cohort comprised RET fusion-positive or BRAF V600E-positive cases in the dose-zero or High-AS group, with one paired tumor and normal sample, both tissues passing the outlier screen, and relative purity ≥0.6 <!-- N-70 --> (n=63 <!-- N-08 -->). The RET fusion-positive subset supplied the REO construction bands. RET fusion-positive PTCs in the Low-AS and Mid-AS bands formed the intermediate-band application set and were deliberately not excluded by the construction-cohort outlier or purity rules (n=36 <!-- N-10 -->). Their quality metrics were reported without exclusion authority.

### Analysis contrasts and interpretation

We defined four contrasts comparing High-AS with dose-zero cases in tumor and matched normal tissue within each driver stratum. Their analysis-file labels were `R_Tumor` for tumor tissue from RET fusion-positive PTCs, `R_Normal` for their matched normal tissue, `B_Tumor` for tumor tissue from BRAF V600E-positive PTCs, and `B_Normal` for their matched normal tissue. These are contrast labels, not sample-group labels. Within each contrast, X denoted the dose-zero group and Y the High-AS group. The Brunner–Munzel relative effect (probabilistic index; Brunner and Munzel 2000), denoted θ, was θ=Pr(X<Y)+0.5Pr(X=Y); values of θ>0.5 indicated that expression tended to be higher in High-AS cases.

Higher Criticism (HC; Donoho and Jin 2004) was the single primary contrast-level test for RET-tumor. RET-normal assessed a corresponding normal-tissue difference, BRAF-normal supplied the other driver stratum for the prespecified normal-tissue comparison, and BRAF-tumor completed the driver-stratified tumor analysis. BRAF-tumor was direction-agnostic and was not treated as a specificity control. Because the secondary contrasts answered distinct questions and were not combined to declare overall study success, no across-contrast adjustment was applied. All four results are reported, with no claim of an across-contrast family-wise error rate or study-wide false discovery rate (FDR).

Table 2 records the interpretation map assigned before the finalized analysis produced the results reported here. It designated RET-tumor as the expected signal-bearing contrast, specified the normal-tissue question, and required a confounding-first reading of potential BRAF-normal findings. Results outside this map are labeled hypothesis-generating.

### Covariate disclosure and estimand

For exposed cases, we incorporated age at exposure into AS as a radiation-risk modifier and did not additionally treat it as a nuisance covariate in the expression analysis. Because dose-zero cases were born after the accident whereas exposed cases were alive in 1986, the groups are structurally separated by birth cohort. We report group distributions and Hodges–Lehmann age-at-surgery differences with 95% bootstrap intervals <!-- N-63 --> to characterize this structure, not to test confounding or demonstrate balance. We also report sex composition because sex-chromosome expression can differ mechanically between groups, and RET fusion partners because partner distributions have differed between radiation-associated and sporadic tumors (Nikiforov et al. 1997; Rabes et al. 2000).

We fitted no common covariate-adjusted model because these variables lack a common nuisance interpretation. Conditioning on the age inputs to AS would replace the prespecified marginal AS-band contrast with an age-conditional estimand; relative purity was expression-derived and used for eligibility; and RET fusion partner may lie on the biological pathway connecting exposure history, tumor initiation, and expression. Their joint adjustment would answer a different question and would not identify a radiation-specific effect.

Relative tumor purity, age and birth cohort, sex, fusion partner, and unmeasured technical structure may also align with AS-band labels. This does not make intended AS inputs errors in the metric; it limits radiation-specific attribution. The estimates therefore describe associations between AS-defined groups. No sequencing- or processing-batch field was available in the analysis objects; this absence is reported as a limitation rather than treated as evidence of balance.

### Normalization and gene-level inference

We retained protein-coding genes using filterByExpr and normalized their counts separately for each contrast using an in-house implementation of three-iteration DEGES-MUREN (Kadota et al. 2012; Sun et al. 2013; Feng and Li 2021). For scaling-factor estimation, we replaced the original model-based differential-expression screen with a permutation Brunner–Munzel screen at Storey q<0.10 (Storey 2002) and excluded the genes it flagged as potentially differentially expressed <!-- N-71 -->. Supplementary Methods and Table S2 report the detailed filtering, floorPDEG rule, convergence, and scaling factors <!-- N-15 -->.

For each gene, we enumerated all C(n,nx) allocations of a studentized Brunner–Munzel statistic (two-sided) on the finalized normalized matrix. Enumeration was exhaustive and therefore introduced no Monte Carlo error. Because the DEGES exclusion screen used the observed group labels, gene-level enumeration was conditional on the resulting normalized matrix and exchangeability of the retained samples. Gene-level lists used the prespecified Storey q-value rule with a plug-in estimate of the null proportion π0 at λ=0.5 and q<0.10 <!-- N-05 -->. We describe the resulting genes as meeting this rule rather than claiming FDR control under arbitrary gene dependence.

Contrast-level evidence was summarized with a prespecified Higher Criticism statistic scanning the lower 10% of the ordered p-values. We chose it to aggregate moderately small p-values rather than depend on one extreme gene or a fixed gene-count threshold. Its p-value came from 9,999 label shuffles within the contrast (seed 19860426) <!-- N-04 -->. Count and maximum statistics and the full rejection curve were retained as descriptive diagnostics. We report these p-values continuously and do not assign a binary label to each contrast. <!-- C-16 -->

### Gene-set inference

We ranked genes by tie-averaged normal scores of the signed Brunner–Munzel statistic and evaluated an in-house weighted running-sum statistic only at tie-block boundaries. On tie-free input, this statistic equalled the standard gene set enrichment analysis (GSEA) statistic (Subramanian et al. 2005) <!-- N-72 -->. We used the 9,999 saved contrast-specific label shuffles <!-- N-07 --> to form the set-level reference distribution and adjusted sign-conditional permutation p-values within each collection using the Benjamini–Hochberg procedure (BH; Benjamini and Hochberg 1995), applying q<0.10 without making a cross-collection claim.

We tested Hallmark, selected C2 canonical-pathway subcollections (Reactome, WikiPathways, KEGG MEDICUS, BioCarta, and PID), C5 GO Biological Process, and a radiation-curated C2:CGP family <!-- N-29, N-55 -->; sets outside 15–500 genes were excluded <!-- N-72 -->. <!-- C-06 -->
An earlier held-out complete-null assessment was used to select the set-level procedure before the reported real-data run, and the selected procedure's final operating characteristics were assessed with a 9,999-shuffle reference pool <!-- N-06 -->. A one-set 1.15-fold spike-in provided a single-scenario positive-control check <!-- N-30 -->. These exercises did not establish FDR control under every partial-alternative configuration. Full procedures are in Supplementary Methods.

<!-- C-14 -->
As a descriptive complement, the primary RET-tumor q<0.10 list was tested for one-sided hypergeometric over-representation in the same collections, using all genes tested in that contrast as the universe and BH adjustment within each family-by-list combination <!-- N-59 -->. These gene-sampling p-values do not represent subject-level randomization; experiment-level interpretation rests on the label-permutation analysis (Goeman and Bühlmann 2007).

### Between-stratum concordance

For each tissue, we compared signed gene-level Brunner–Munzel statistics from the RET and BRAF contrasts across shared genes using Spearman correlation. We generated a reference interval by independently shuffling labels within each stratum. The normal-tissue comparison was hypothesis-bearing, whereas the tumor comparison was a descriptive completion. A correlation outside this interval indicated shared label-aligned structure but could not distinguish exposure from a shared covariate.

### REO panel

We based the REO panel on within-sample gene ordering, which does not require between-sample normalization (Geman et al. 2004; Wang et al. 2015; Guan et al. 2016). We recalculated transcripts per million (TPM) from the selected count assay using exon-union gene lengths. We generated candidate pairs from the 500 RET-tumor genes with the largest Brunner–Munzel effect magnitudes, rather than from the q<0.10 list, and applied prespecified stability, reversal, and redundancy rules to select 10 pairs without gene reuse <!-- N-37 -->. For each sample, the reversal score counted, from 0 to 10, the pairs whose expression difference lay outside the dead zone and had the sign opposite to the dose-zero construction reference. We set the classification cutoff at the maximum reversal score among dose-zero construction samples and classified scores above that cutoff as positive.

<!-- C-08 の手法側 -->
We defined binary and graded readouts for the panel. In the dose-zero and High-AS construction bands, the binary readout characterized separation between the groups used to construct the panel. Demonstrating this separation was a prerequisite for application elsewhere, but it was not an unbiased estimate of classification performance because the same cases informed pair selection and the dose-zero cases set the cutoff. We then applied the classification rule unchanged to Low- and Mid-AS RET fusion-positive PTCs and treated the resulting classifications as descriptive application results.

We summarized graded scores across all four bands, with dose-zero and High-AS serving as construction anchors and the panel applied without refitting to Low-AS and Mid-AS. Within the two intermediate bands, the prespecified one-sided Monte Carlo Brunner–Munzel comparison of Mid-AS with Low-AS (seed 19860426) <!-- N-75 --> assessed whether scores tended to be higher in Mid-AS than in Low-AS. It did not test linearity, four-band monotonicity, or a dose–response form. Neither this comparison nor the four-band display constituted independent validation because the panel and score direction came from the construction bands. <!-- C-09 --> We also applied PC-OD separately within the Low-AS and Mid-AS groups as a non-exclusion diagnostic. For intermediate-band tumors with matched normal tissue, we estimated relative purity on a common scale across all four RET groups and summarized the band–score association after rank-scale adjustment for purity using partial Spearman correlation and a descriptive one-sided permutation reference <!-- C-12 -->. We did not use either diagnostic to exclude cases or modify the panel.

### External gene-list comparison

We cross-referenced two classes of published radiation-associated gene lists against each contrast's q<0.10 list: qRT-PCR-validated lists from Abend et al. (2012, 2013) and Dom et al. (2012), with CLIP2 as a single-gene reference, and the shared-tissue, normal-tissue, and tumor-tissue multivariate transcriptomic signatures from Ory et al. (2026). <!-- C-13 --> We kept the two classes distinct and reported descriptive membership counts without performing an enrichment test <!-- N-53 -->.

### Reproducibility and AI use

We independently executed the full publication pipeline on two x86-64 machines using MD5-identical raw inputs <!-- N-52 --> and a date-pinned container with R 4.5.3 (R Core Team 2026) on Ubuntu 24.04 <!-- N-02 -->, reference BLAS/LAPACK 3.12.0 <!-- N-03 -->, four workers <!-- N-04 -->, and fixed seeds. Both runs passed 415 tests with no failures <!-- N-51 --> and produced identical primary artifacts. The analysis code, versioned inputs, and container recipe accompany the paper; Table S5 lists package versions.

Pipeline design, the prototype scripts and all analytical and interpretive decisions are the authors'. Generative AI assistants (Claude Fable 5, Anthropic; GPT-5.6sol, OpenAI) were used for script debugging and refactoring, and for drafting and editing manuscript text under the authors' fine-grained direction. The writing assistance was directed at statistical precision: stating estimands, inferential scope and the qualifiers each claim carries exactly, and holding technical usage consistent across the manuscript. Drafted passages render content the authors specified and ratified, and were verified against versioned ledgers — every quantitative statement against the analysis outputs, every citation against its source document, and the wording of claim-bearing sentences against an author-ratified claim map — with author approval at each step. No scientific decision was delegated to these tools. Only open-access data were handled by them; all resulting code passes the verification described above (test suite, cross-machine reproduction, automated equivalence checks), and the authors take full responsibility for the entire content.

## Results

### 1. Cohort(C-10, C-15)

<!-- C-10 -->
Of the 440 cases of the REBC-THYR cohort (Morton et al. 2021) <!-- N-11 -->, the main analysis cohort comprised 63 paired, driver-stratified cases — 9 B_High, 27 B_Sporadic, 15 R_High and 12 R_Sporadic <!-- N-09 --> — reached through the pre-specified flow of driver classification, band eligibility, pairing, outlier screening and purity thresholding (Table 1; Supplementary Table S1). Most of the reduction reflects the pre-specified eligibility restrictions — driver classification and the restriction to the extreme bands — rather than technical losses to pairing, outlier screening or purity; the full flow, stratified by driver, is shown in the Fig. 1 <!-- N-08 -->. The REO evaluation set added 36 paired RET tumors of the intermediate bands (17 R_Low, 19 R_Mid
<!-- N-10 -->). In both driver strata the High group sat somewhat older at surgery than the Sporadic group <!-- N-12, N-13 -->.
<!-- C-15 -->
The three disclosed candidate confounders (Methods) are reported in Table 1: the per-group age distributions, with interval estimates of the between-group age difference as table footnotes (disclosed, not tested <!-- N-64, N-65 -->), and the per-group sex and fusion-partner compositions. Age at exposure exists only for exposed cases and admits no between-group comparison <!-- N-63 -->.

【訳】REBC-THYR コホート(Morton et al. 2021)の 440 症例 <!-- N-11 --> のうち、主解析コホートは、driver 分類・帯適格性・ペア有無・外れ値スクリーン・純度閾値という事前指定のフロー(表1; 補足表S1)を経た 63 のペア付き driver 層別症例 — B_High 9・B_Sporadic 27・R_High 15・R_Sporadic 12
<!-- N-09 --> — である。削減の大半は事前指定の適格性制限 — driver 分類と両端帯への制限 — によるもので、ペア・外れ値スクリーン・純度による技術的損失ではない。フローの全体は driver 層別で図1に示す <!-- N-08 -->。REO 評価セットとして中間帯のペア付き RET 腫瘍 36 例(R_Low 17・R_Mid 19 <!-- N-10 -->)を加えた。いずれの driver 層でも High 群は Sporadic 群よりやや高い手術時年齢に分布した <!-- N-12, N-13 -->。開示対象の候補交絡因子3つ(Methods)は表1に報告する: 群別の年齢分布(群間年齢差の区間推定は表脚注 — 開示であって検定ではない <!-- N-64, N-65 -->)、および群別の性構成と融合パートナー構成。被曝時年齢は被曝症例にのみ定義され、群間比較は成立しない <!-- N-63 -->。

### 2. Gene-level differential expression(C-01〜C-04, C-16, 410)

<!-- C-01 -->
In R_Tumor — the contrast where the pre-specified interpretation map expects the signal — 1,765 of 15,621 tested genes differed between the High and Sporadic groups at Storey q < 0.10 <!-- N-16, N-15 -->, and the pre-specified contrast-level omnibus supported the presence of signal (Higher Criticism p = 0.0112
<!-- N-20 -->) (Fig. 2, Supplementary Fig. S1, Table 3). <!-- C-02 --> The discovered genes ran in both directions: 971 higher and 794 lower in the exposed group <!-- N-17 -->. The sex-chromosome annotation of the discovered list (Methods) counted 57 X-linked genes among the 1,765 (36 higher, 21 lower in the exposed group) and no Y-linked gene; the contrast's tested genes include 572 X-linked and 2 Y-linked <!-- N-82 -->. <!-- C-03 --> In R_Normal no gene reached q < 0.10 and the
omnibus lent no support (HC p = 0.3199) <!-- N-16, N-20 -->. <!-- C-04 --> B_Tumor likewise yielded no discovery (0 genes at q < 0.10; HC p = 0.1815)
<!-- N-16, N-20 -->; under the pre-fixed reading this cell is direction-agnostic, and its quiet is not read as a specificity control. <!-- C-16 -->
In B_Normal the evidence is reported as it stands, without a binary label: one gene crossed the gene-level threshold (BHLHB9 — an X-linked gene <!-- N-80 --> — effect 0.967, q = 0.013
<!-- N-22 -->); the pre-specified primary omnibus gave HC p = 0.0773 <!-- N-20 --> while the descriptive max-statistic row reached p = 0.0125
<!-- N-23 -->, and the contrast's π0 estimate (0.727) sat below those of the other two quiet contrasts (0.955, 0.943) <!-- N-19 -->. The reading pre-assigned
to this pattern is taken up in Discussion.

【訳】事前指定の解釈マップがシグナルを期待する対比である R_Tumor では、検定対象 15,621 遺伝子のうち 1,765 が Storey q < 0.10 で High 群と Sporadic 群の間で発現差を示し <!-- N-16, N-15 -->、事前指定の対比レベル・オムニバスがシグナルの存在を支持した(Higher Criticism p = 0.0112 <!-- N-20 -->)(図2・補足図S1・表3)。発見遺伝子は双方向に分布した: 被曝群で高発現 971・低発現 794 <!-- N-17 -->。発見リストの性染色体注釈(Methods)では、1,765 遺伝子中 57 が X 連鎖(被曝群で高発現 36・低発現 21)、Y 連鎖は 0 だった。当該対比の検定遺伝子には X 連鎖 572・Y 連鎖 2 が含まれる <!-- N-82 -->。R_Normal では q < 0.10 の遺伝子はなく、オムニバスの支持もなかった(HC p = 0.3199)<!-- N-16, N-20 -->。B_Tumor も同様に発見なし(q < 0.10 は 0 遺伝子; HC p = 0.1815)<!-- N-16, N-20 --> — 事前固定の読みに従い、このセルは方向不可知であり、その静けさを特異性の対照とは読まない。B_Normal の証拠は2値ラベルなしにそのまま報告する: 1遺伝子が遺伝子レベル閾値を越え(BHLHB9 — X 連鎖遺伝子 <!-- N-80 --> — effect 0.967、q = 0.013 <!-- N-22 -->)、事前指定の主オムニバスは HC p = 0.0773
<!-- N-20 -->、記述的な max 統計量行は p = 0.0125 <!-- N-23 -->、この対比の π0 推定値(0.727)は他の2つの静かな対比(0.955・0.943)より低かった <!-- N-19 -->。このパターンに
事前割当された読みは Discussion で扱う。

### 3. Gene-set level(C-05, C-06, 420 + D6)

<!-- C-05 -->
At the gene-set level, the calibrated test declared no set at q_bh < 0.10 in any of the 16 contrast × collection cells <!-- N-27 --> (four MSigDB collections, Methods; 6,141–6,242 sets per contrast after filtering <!-- N-28 -->). The smallest adjusted value anywhere was q = 0.114, in the B_Tumor × radiation cell <!-- N-28 --> (complete listing in Supplementary Data 1; Supplementary Table S6).
<!-- C-06 -->
Under signal-free inputs — label permutations pushed through the identical procedure (Methods) — the set-level machinery produced at least one discovery in 102 of 1,600 replicates pooled across the 16 contrast × collection cells (0.064; nominal 0.10 <!-- N-56 -->), with a single disclosed excess (B_Normal/Hallmark: 0.18, 95% CI 0.110–0.270 <!-- N-25 -->). Under the synthetic positive control — one Hallmark set inflated 1.15-fold in one group (Methods) — the planted set was recovered at rank 1 of 50 (q = 0.0101 <!-- N-31 -->), with no other set at q < 0.10 <!-- N-32 --> (per-cell calibration in Supplementary Fig. S2 and Supplementary Table S7).

【訳】遺伝子セットレベルでは、較正済み検定は 16 の対比 × collection セルのいずれでも q_bh < 0.10 のセットを宣言しなかった <!-- N-27 -->(MSigDB 4 コレクション、Methods 参照; フィルタ後は対比あたり 6,141–6,242 セット <!-- N-28 -->)。全セルで最小の調整値は B_Tumor × radiation セルの q = 0.114 だった <!-- N-28 -->(全結果は補足データ1・補足表S6)。シグナルなしの入力 — 同一手続きに通したラベル置換(Methods)— の下では、セットレベル機構は 16 の対比 × collection セル合算 1,600 レプリケート中 102 で1つ以上の発見を生じ(0.064; 名目 0.10 <!-- N-56 -->)、開示済みの超過は1セルのみだった(B_Normal/Hallmark: 0.18、95% CI 0.110–0.270 <!-- N-25 -->)。合成陽性対照 — 1つの Hallmark セットを一方の群で 1.15 倍(Methods)— の下では、埋め込んだセットが 50 中 rank 1(q = 0.0101
<!-- N-31 -->)で回収され、それ以外に q < 0.10 のセットはなかった <!-- N-32 -->
(セル別較正は補足図S2・補足表S7)。

### 4. Composition of the discovered list(C-14, 430)

<!-- C-14 -->
In the descriptive annotation of the discovered list (hypothesis-generating; Methods), the 794 genes lower in the exposed group were strongly concentrated in proliferation, cell-cycle and DNA-repair programs (E2F targets 46/199, G2M checkpoint 41/198, Reactome DNA repair 68/322 <!-- N-61, N-62 -->), a single theme that extends into the radiation-curated family, whose leading flagged sets are themselves cell-cycle genes responding to irradiation (46 of 126 in the down list, expected 6.4 <!-- N-60 -->). The 971 genes higher in the exposed group showed no such concentration in any curated family <!-- N-59 --> (full results in Supplementary Table S3).

【訳】発見済みリストの記述的注釈(仮説生成; Methods)では、被曝群で低発現の 794 遺伝子は増殖・細胞周期・DNA 修復プログラムに強く集中し(E2F targets 46/199・G2M checkpoint 41/198・Reactome DNA repair 68/322 <!-- N-61, N-62 -->)、この単一テーマは放射線キュレーション・ファミリーにも及ぶ — そこでフラグが立った上位セット自体が照射に応答する細胞周期遺伝子である(down リストで 126 中 46、期待 6.4 <!-- N-60 -->)。被曝群で高発現の 971 遺伝子には、どのファミリーにもそのような集中はなかった
<!-- N-59 -->(全結果は補足表S3)。

### 5. Between-stratum concordance of the exposure contrast(C-07, C-17)

<!-- C-07 -->
The pre-specified between-stratum comparison in normal tissue — Spearman correlation, across strata, of the per-gene signed statistics of the exposure contrast — gave rho = +0.376 over 15,459 shared genes, inside its label-shuffle reference interval ([−0.46, +0.46]; two-sided p = 0.1199)
<!-- N-34 -->. With no within-contrast signal in R_Normal <!-- N-16 -->, the pre-fixed reading applies: whether the two normal-tissue contrasts share an exposure trace is not identifiable here. <!-- C-17 --> The symmetric
tumor-pair comparison, computed as a design completion, is reported in Supplementary Table S8 <!-- N-33 --> and taken up as hypothesis-generating in Discussion.

【訳】正常組織における事前指定の層間比較 — 曝露対比の遺伝子別符号付き統計量を層間で Spearman 相関 — は、共有 15,459 遺伝子で rho = +0.376 となり、ラベルシャッフルの参照区間([−0.46, +0.46]; 両側 p = 0.1199)の内側だった <!-- N-34 -->。R_Normal に対比内シグナルがない <!-- N-16 --> ため、事前固定の読みが適用される: 2つの正常組織対比が曝露痕跡を共有するか否かは、ここでは識別できない。対称補完として計算した腫瘍ペアの比較は補足表S8に報告し <!-- N-33 -->、Discussion で仮説生成として扱う。

### 6. REO grading(C-08, C-09, C-12, 510–530)

<!-- C-08 -->
Panel construction evaluated 57,694 candidate pairs from the training pool's 317 up- and 182 down-genes; 153 passed all criteria, with median r shifts of 1.159–4.700 and reversal rates 0.53–0.87 <!-- N-35, N-36 -->. The 10-pair REO panel separated its training groups as designed (all 12 R_Sporadic negative, 13 of 15 R_High positive; boundary at score > 2
<!-- N-38, N-37 -->; panel composition in Supplementary Table S4 <!-- N-39 -->). Applied unchanged to the intermediate bands, the per-case reversal score rose in band order — medians 0 / 1 / 4 / 6 for Sporadic / Low / Mid / High <!-- N-41 --> (Fig. 3) — with R_Low classified 9 negative / 8 positive
and R_Mid 8 / 11 <!-- N-42 -->. The one comparison available out of sample (Mid vs Low; direction pre-specified, Methods) gave one-sided Brunner–Munzel p = 0.1127, effect P(Low<Mid) = 0.616 <!-- N-40 -->. The graded profile is reported as a descriptive observation; no dose–response form is assumed or claimed.
<!-- C-09 -->
The outlier screen mirrored from training flagged no case in either evaluation band (0 of 17 and 0 of 19 <!-- N-43 -->); the screen carries no exclusion authority, and the evaluation was run once on all eligible cases.
<!-- C-12 -->
Within the evaluation bands the reversal score correlated with tumor purity (pooled Spearman +0.538 <!-- N-45 -->), while the band–score correlation was +0.142 and, conditioned on purity, +0.146 (partial Spearman; one-sided permutation p = 0.2162) <!-- N-48, N-46 -->; band and purity themselves were nearly uncorrelated (+0.036 <!-- N-48 -->). Purity-stratified comparisons are given in Supplementary <!-- N-47 -->.

【訳】パネル構築では、訓練プールの up 317・down 182 遺伝子から 57,694 候補ペアを評価し、153 が全基準を通過した(中央値 r シフト 1.159–4.700、逆転率 0.53–0.87)<!-- N-35, N-36 -->。10 ペアの REO パネルは訓練群を設計どおり分離した(R_Sporadic は 12 例全て negative、R_High は 15 例中 13 が positive; 境界は score > 2 <!-- N-38, N-37 -->;パネル構成は補足表S4<!-- N-39 -->)。パネルを変更せず中間帯に適用すると、症例別の逆転スコアは帯順に上昇した — 中央値は Sporadic / Low / Mid / High で 0 / 1 / 4 / 6
<!-- N-41 -->(図3)— 分類は R_Low が negative 9 / positive 8、R_Mid が 8 / 11 <!-- N-42 -->。out-of-sample で検査可能な唯一の比較(Mid 対 Low; 方向は事前指定、
Methods 参照)は片側 Brunner–Munzel p = 0.1127、効果 P(Low<Mid) = 0.616 だった
<!-- N-40 -->。この段階的プロファイルは記述的観察として報告し、線量反応の形は仮定も主張もしない。訓練側から鏡映した外れ値スクリーンはどちらの評価帯でも該当例を出さなかった(17 例中 0・19 例中 0 <!-- N-43 -->)— スクリーンは除外権限を持たず、評価は適格全例で一度だけ実施した。評価帯内では逆転スコアは腫瘍純度と相関し(pooled Spearman +0.538 <!-- N-45 -->)、帯–スコア相関は +0.142、純度を条件付けると +0.146(偏 Spearman; 片側置換 p = 0.2162)<!-- N-48, N-46 -->、帯と純度自体はほぼ無相関だった(+0.036 <!-- N-48 -->)。純度層別の比較は Supplementary に示す <!-- N-47 -->。

### 7. External anchor cross-reference(C-13)

<!-- C-13(事実側 — 読みは Discussion) -->
Cross-referencing the externally validated radiation-associated gene lists (Methods) against each contrast's discovered genes returned zero overlap in 19 of 20 cells, including every tissue-matched cell <!-- N-53 -->. The single non-zero cell was cross-tissue: S100A10, from the Dom normal-tissue list, appeared among the R_Tumor discoveries, with direction opposite to the original report <!-- N-54 -->.

【訳】外部で検証済みの放射線関連遺伝子リスト(Methods)を各対比の発見遺伝子と照合した結果、20 セル中 19 で重なりはゼロであり、組織対応のセルは全てゼロだった <!-- N-53 -->。唯一の非ゼロセルは組織対応外だった: Dom の正常組織リスト由来の S100A10 が R_Tumor の発見遺伝子に現れ、方向は原報告と逆であった <!-- N-54 -->。

## Discussion

<!-- 執筆状態: 2026-08-15 全面改稿(研究者指示: 学術論文の体裁へ)。規律句は §1 で一括宣言し、各節の反復を撤去。C 対応は不変(C-13 の事実側は Results §7 へ移送、C-17/C-07 は concordance 節に統合、年齢単独節は候補交絡因子節へ改組)。見出しと順序は仮 -->

### Principal findings(C-01, C-04, C-05)

<!-- C-01/C-04/C-05 の受け+解釈水準の一括宣言(事前固定の読みの範囲を明示 — C-17 は対象外なので「マップ由来の読み」に限定) -->

Conditioned on driver and tested within each contrast, the exposure comparison produced gene-level signal in exactly one place: RET tumors, where 1,765 genes crossed q < 0.10 and the contrast-level omnibus supported the presence of signal <!-- N-16, N-20 -->. We read this as an exposure-associated expression difference within RET tumors <!-- C-01 -->. Two design properties discipline that reading. Under the exact label-permutation null, structure uncorrelated with the exposure label — including any general noisiness or instability of the RET-fusion transcriptome — is controlled at the declared error rates by construction; what could mimic the finding is covariate structure correlated with the label, and the disclosed candidates — age, sex and fusion-partner composition — are taken up below. And the association names no mechanism: the wording stays neutral as to whether a shared radiation trace feeds RET oncogenesis or the type of initiating lesion associates with driver selection (Limitations). Throughout this section, readings drawn from the interpretation map (Methods, Table 2) were assigned to outcome patterns before the finalized analysis produced the results reported here; interpretations that go beyond the map are labelled as hypothesis-generating. The remaining cells were quiet at the gene level — no discovery in R_Normal or B_Tumor, and B_Tumor's quiet is not read as a specificity control <!-- C-04 --> — while B_Normal showed a two-level divergence taken up below. At the set level, the calibrated test declared nothing in any of the 16 cells <!-- N-27 -->: whatever organization the discovered list displays along curated set boundaries (below), it does not reach the declared error rate under the label-permutation null <!-- C-05 -->.

【訳】driver で条件付け、各対比内で検定した結果、曝露比較が遺伝子レベルのシグナルを生んだのはただ一箇所 — RET 腫瘍であり、1,765 遺伝子が q < 0.10 を越え、対比レベルのオムニバスがシグナルの存在を支持した <!-- N-16, N-20 -->。我々はこれを RET 腫瘍内の曝露関連発現差として読む <!-- C-01 -->。この読みは設計の2性質が規律する。exact ラベル置換帰無の下では、曝露ラベルと無相関な構造 — RET 融合腫瘍のトランスクリプトームの一般的なノイズや不安定性を含む — は構成的に宣言誤り率で制御される。所見を模倣し得るのはラベルと相関する共変量構造であり、開示済みの候補 — 年齢・性・融合パートナー構成 — は後述で扱う。そして関連は機構を名指ししない: 共有放射線痕跡が RET 発がん機序に接続するのか、起始病変のタイプが driver 選択と連関するのかについて、文言は中立を保つ(Limitations)。本節を通じて、解釈マップ(Methods、表2)に由来する読みは確定解析が本稿で報告する結果を生成する前に結果パターンへ割り当てられたものであり、マップを越える解釈には仮説生成のラベルを付す。残りのセルは遺伝子レベルで静かであり — R_Normal と B_Tumor に発見はなく、B_Tumor の静けさは特異性の対照とは読まない <!-- C-04 --> — B_Normal は後述する二水準の乖離を示した。セットレベルでは較正済み検定は 16 セルのいずれでも何も宣言しなかった <!-- N-27 -->: 発見済みリストがキュレート済みセット境界に沿って示す組織化(後述)がどうであれ、ラベル置換帰無の下で宣言済み誤り率には達しない <!-- C-05 -->。

### B_Normal: a two-level divergence read as confounding first(C-16)

<!-- C-16: 規則5の条件適用。D6 超過 0.18(N-25)の受けもここ。「q<0.10 陽性論」への受けは Q-16 -->

B_Normal is the one contrast where the two evidence levels diverge: a single discovered gene (BHLHB9 <!-- N-22 -->) against a primary omnibus lending no support (HC p = 0.0773 <!-- N-20 -->). No binary label is assigned (Methods); the pattern receives the reading the map fixed in advance — confounding first <!-- C-16 -->. The stratum's disclosed covariate asymmetries supply the first-line candidates: the one discovered gene is X-linked <!-- N-80 --> in the comparison carrying the largest sex-composition imbalance (23 of 27 versus 4 of 9 female <!-- N-13 -->), and the same stratum carries the largest between-group age difference (+8.0 years [3.0 to 12.0] <!-- N-65 -->). Neither candidate is tested; both are disclosed. One calibration caveat travels with this cell: the only held-out excess of the set-level calibration sits in B_Normal/Hallmark (0.18, 95% CI 0.110–0.270 <!-- N-25 -->), and any reading of that cell carries it.

【訳】B_Normal は二つの証拠水準が乖離する唯一の対比である: 発見遺伝子1件(BHLHB9 <!-- N-22 -->)に対し、主オムニバスは支持を与えない(HC p = 0.0773 <!-- N-20 -->)。2値ラベルは付与しない(Methods)。このパターンは、マップが事前に固定した読み — 交絡第一 — を受ける <!-- C-16 -->。第一候補はこの層の開示済み共変量非対称が与える: 唯一の発見遺伝子は X 連鎖であり <!-- N-80 -->、この対比は性構成の非対称が最大(女性 27 例中 23 対 9 例中 4 <!-- N-13 -->)、かつ群間年齢差も最大(+8.0 年 [3.0〜12.0] <!-- N-65 -->)である。どちらの候補も検定せず、どちらも開示する。このセルには較正上の注意が1つ付随する: セットレベル較正の唯一の held-out 超過が B_Normal/Hallmark にあり(0.18、95% CI 0.110–0.270 <!-- N-25 -->)、当該セルのどの読みもこれを携行する。

### Candidate confounders(C-15, Q-13)

<!-- C-15 の開示推定。Methods の Candidate confounders と対構造。2026-08-16 引用監査により4観察バウンドを解体(Vriens 2011 は振幅ゲート由来の負・年齢幅が原典に不在、Coclet 1989 は年齢効果を確立しないため両者撤去; B_Tumor の null を内部対照に使う形も規則4と衝突するため撤去)。残すのは開示量とその非対称のみ — Q-13 も同日改訂 -->

The design discloses three candidate confounders rather than adjusting for them (Methods); what follows is disclosure of what each could carry, not a demonstration that it does not. Age at surgery differs between groups — Hodges–Lehmann +2.5 years [−1.0 to 6.0] in the RET stratum and +8.0 years [3.0 to 12.0] in the BRAF stratum <!-- N-64, N-65 -->. Age is not entered as a covariate because age at exposure is a component of the assigned-share definition, so adjusting for it would remove part of the exposure metric itself (Methods); by the same construction the bands carry an age structure by design, and no comparison internal to this cohort can separate radiation attributability from the biology of young age at exposure. Two disclosed quantities qualify the age candidate without eliminating it. The interval estimates are asymmetric across the strata: in the stratum carrying the reported finding the difference is small and its interval crosses zero, whereas in the stratum the protocol reads as confounding first (B_Normal, above) it is larger and its interval does not <!-- N-64, N-65 -->. And within the exposed RET bands the graded score does not follow age at surgery — the median falls 30 / 25 / 23 years across Low / Mid / High <!-- N-12 --> while the reversal score rises 1 / 4 / 6 <!-- N-41 --> — although this speaks only to the narrow possibility that older cases score higher, since the banding is itself a function of age at exposure. Age, with latency inseparable from it, therefore remains a named and unexcluded candidate for the reported difference (Limitations). The other two candidates are disclosed in the same spirit: sex is named for B_Normal above — while in R_Tumor the discovered list's sex-chromosome membership is proportionate to the tested background (57 of 1,765 versus 572 of 15,621 X-linked; Results <!-- N-82 -->) — and the fusion-partner composition shows no gradient across the bands (Table 1 <!-- N-12 -->).

【訳】本設計は候補交絡因子3つを調整でなく開示する(Methods)。以下は各候補が担い得る範囲の開示であって、担っていないことの証明ではない。手術時年齢は群間で異なる — Hodges–Lehmann 中央値差は RET 層で +2.5 年 [−1.0〜6.0]、BRAF 層で +8.0 年 [3.0〜12.0] <!-- N-64, N-65 -->。年齢を共変量に入れないのは、被曝時年齢が assigned share の定義の構成要素であり、調整が曝露指標そのものの一部除去になるためである(Methods)。同じ構成により帯は設計上の年齢構造を持ち、本コホート内部のどの比較も、放射線起因性と被曝時年齢の生物学とを分離できない。開示された2つの量が、年齢候補を排除せずに限定する。区間推定は層間で非対称である: 報告された所見を担う層では差は小さく区間はゼロを跨ぐ一方、プロトコルが交絡第一と読む層(上の B_Normal)では差は大きく区間はゼロを跨がない <!-- N-64, N-65 -->。また被曝 RET 帯内で段階的スコアは手術時年齢に追随しない — 中央値は Low / Mid / High で 30 / 25 / 23 歳と下がり <!-- N-12 -->、逆転スコアは 1 / 4 / 6 と上がる <!-- N-41 --> — ただしこれは「高齢例ほど高スコア」という狭い可能性に答えるだけであり、帯化自体が被曝時年齢の関数だからである。したがって年齢は(不可分の潜伏期とともに)、報告された差に対する名指しされた排除不能の候補として残る(Limitations)。残る2候補も同じ趣旨で開示する: 性は上の B_Normal で名指し済みであり — 一方 R_Tumor の発見リストの性染色体所属は検定対象の背景と比例的である(発見 1,765 中 X 連鎖 57、検定対象 15,621 中 572; Results <!-- N-82 -->)— 融合パートナー構成は帯を跨ぐ勾配を示さない(表1<!-- N-12 -->)。

### Reading the over-representation annotation(C-14)

<!-- C-14(年齢の一般論は候補交絡因子節、ここは内容特異の受け: 純度の逆向き + PTC 腫瘍組織で加齢による増殖減弱は未報告 — Chatsirisupachai 2021、研究者示唆 2026-08-16) -->

The over-representation annotation of the R_Tumor list is read under a gene-sampling null that ignores co-expression and is anticonservative relative to the calibrated set-level inference; its flags are candidates supported only under that weaker null. Three observations bound its meaning. First, the flags across all four families express a single theme — relative attenuation of proliferation, cell-cycle and DNA-repair programs in the exposed group — rather than several independent lines of evidence, and the radiation-curated flags ride the same cell-cycle content. Second, the calibrated set-level test ranked the same Hallmark sets first without crossing its threshold <!-- N-28 -->: the two procedures disagree about evidence, not about ordering. Third, whether the attenuation reflects exposure-associated biology or the groups' age and latency structure cannot be decided here, though two observations bear on the alternatives. Tumor purity runs the wrong way for a contamination artifact: the High band is the purest (0.814 versus 0.690 in Sporadic <!-- N-44 -->), yet the proliferation programs read as attenuated rather than strengthened in the exposed group <!-- N-61 -->. And an attenuation of proliferation and cell-cycle programs with increasing age is not an established property of papillary thyroid tumor tissue: a genome-wide analysis of TCGA thyroid carcinoma against patient age, over a span of decades and without any fold-change filter, reported no age-associated gene among the protein-coding genes tested at a Benjamini–Hochberg adjusted p < 0.05 — a stricter level than the one used here (Chatsirisupachai et al. 2021). No claim of the primary analyses rests on this annotation.

【訳】R_Tumor リストの過剰代表注釈は、共発現を無視する gene-sampling 帰無の下で読まれ、較正済みセットレベル推論に対して反保守的である。そのフラグは、この弱い帰無の下でのみ支持される候補である。その意味の範囲を3つの観察が画定する。第一に、4ファミリー全てのフラグは複数の独立した証拠線ではなく単一のテーマ — 被曝群における増殖・細胞周期・DNA 修復プログラムの相対的減弱 — を表現しており、放射線キュレーションのフラグも同じ細胞周期の内容に乗っている。第二に、較正済みセットレベル検定は同じ Hallmark セット群を、閾値を越えないまま最上位に置いた <!-- N-28 -->: 二つの手続きが食い違うのは証拠についてであり、順位についてではない。第三に、この減弱が曝露関連の生物学を反映するのか、群の年齢・潜伏期構造を反映するのかはここでは決められないが、2つの観察が選択肢に関わる。腫瘍純度は混入アーティファクト説と逆向きに働く: High 帯は最も高純度だが(0.814 対 Sporadic 0.690 <!-- N-44 -->)、増殖プログラムは被曝群で強化でなく減弱として読まれる <!-- N-61 -->。そして、加齢に伴う増殖・細胞周期プログラムの減弱は乳頭状甲状腺腫瘍組織で確立した性質ではない: TCGA 甲状腺癌を患者年齢に対して数十年幅・fold change フィルタなしで解析したゲノムワイド研究は、検定した protein-coding 遺伝子のうち Benjamini–Hochberg 調整 p < 0.05(本研究より厳しい水準)で年齢関連遺伝子を報告していない(Chatsirisupachai et al. 2021)。主解析のどの主張もこの注釈に依存しない。

### Concordance across strata(C-17, C-07)

<!-- C-17(tumor: 仮説生成・純度筆頭の識別不能・B_Tumor×radiation q0.114 の断片束ね)+ C-07(normal: 仮説担持・識別不能の宣言)を統合。読みの水準が対になる -->

The two cross-stratum comparisons carry different weights and return different verdicts. The tumor pair, computed as a design completion, resembles each other beyond label-exchange chance: Spearman rho +0.459 over 15,560 shared genes, outside the label-shuffle reference interval ([−0.39, +0.39]; two-sided p = 0.0197 <!-- N-33 -->; Supplementary Table S8) <!-- C-17 -->. As hypothesis generation only: this is what an exposure-associated component visible across driver backgrounds — an oncogene-independent trace — would look like, and it is compatible with B_Tumor's within-contrast quiet, which requires only sub-threshold structure aligned with the R_Tumor contrast. It is equally compatible with shared covariate structure, purity foremost (the exposed groups sit purer at least on the RET side <!-- N-44 -->), and the design cannot separate the two; an adjacent fragment — the closest any calibrated set-level cell came to discovery, B_Tumor × radiation (q = 0.114 <!-- N-28 -->) — points the same way without crossing any threshold. This comparison does not test the shared-glandular-memory hypothesis, and no claim of this study rests on it. The normal-tissue pair, by contrast, is the comparison that carried the hypothesis: long-term molecular memory of exposure (Abend et al. 2013) and shared normal-tissue signatures reported without driver stratification (Ory et al. 2026) predict that the two driver-conditioned normal contrasts should resemble each other, and this comparison is the device that makes that prediction testable under driver conditioning. It returned indeterminacy: rho +0.376, inside its reference interval <!-- N-34 -->, with no within-contrast signal in R_Normal <!-- N-16 -->. A null cross-stratum correlation between contrasts that individually carry no detectable within-contrast signal cannot distinguish absent sharing from sharing too weak to detect; whether the two strata share a normal-tissue exposure trace is therefore not identifiable in this cohort <!-- C-07 -->, rather than shown to be absent; the comparison neither supports nor contradicts the shared-memory reports, and it remains the design under which a larger cohort could decide the question.

【訳】二つの層間比較は異なる重みを持ち、異なる判定を返した。腫瘍ペアは設計上の補完として計算されたもので、ラベル交換の偶然を超えて互いに似ている: 共有 15,560 遺伝子で Spearman rho +0.459、ラベルシャッフルの参照区間([−0.39, +0.39])の外側、両側 p = 0.0197 <!-- N-33 -->(補足表S8)<!-- C-17 -->。仮説生成としてのみ: これは driver 背景を越えて見える曝露関連成分 — oncogene 非依存の痕跡 — が示すはずの姿であり、B_Tumor の対比内の静けさとも矛盾しない(R_Tumor の対比と向きの揃った閾値下構造があれば足りる)。同時に、共有された共変量構造 — 筆頭は純度(少なくとも RET 側で被曝群は高純度 <!-- N-44 -->)— とも等しく整合し、本設計は両者を分離できない。隣接する断片 — 較正済みセットレベルの全セルで発見に最も近づいた B_Tumor × radiation(q = 0.114 <!-- N-28 -->)— も、どの閾値も越えないまま同じ方向を指す。この比較は共有腺記憶仮説を検定せず、本研究のどの主張もこれに依存しない。対照的に、正常組織ペアは仮説を担った比較である: 被曝の長期分子記憶の報告(Abend et al. 2013)と、driver 非層別で報告された正常組織の共有署名(Ory et al. 2026)は、driver で条件付けた二つの正常対比が互いに似ることを予測し、この比較はその予測を driver 条件付けの下で検定可能にする装置である。結果は不確定に終わった: rho +0.376 は参照区間の内側 <!-- N-34 -->、R_Normal に対比内シグナルなし <!-- N-16 -->。各対比が単独では検出可能な対比内シグナルを持たないとき、層間相関の帰無は「共有の不在」と「検出には弱すぎる共有」を区別できない。したがって二つの層が正常組織の曝露痕跡を共有するか否かは、本コホートでは不在と示されたのではなく識別できない <!-- C-07 -->。この比較は共有記憶の報告を支持も否定もしない — そして、より大きなコホートがこの問いを決められる設計として残る。

### External anchors(C-13)

<!-- C-13(読み側 — 事実は Results §7)。対称読みの事前固定 -->

The external anchor counts (Results) are read as description under the symmetric rule fixed in advance. Zero overlap with the qRT-PCR-validated cores is the norm rather than the exception among the prior studies themselves — Dom et al. (2012) tested four earlier tumor-derived signatures, as gene sets, against their own exposed-versus-non-exposed normal-tissue contrast and found none enriched — and the single cross-tissue overlap (S100A10, direction reversed <!-- N-54 --> — Dom's reported list being wholly upregulated in exposed normal tissue, any down-in-exposed overlap from it reverses by construction) is noted without weight. No claim of this study moves with these counts <!-- C-13 -->.

【訳】外部アンカーの員数(Results)は、事前に固定した対称規則の下で記述として読む。qRT-PCR 検証済みコアとの重なりゼロは、先行研究どうしの間でも例外でなく通例であり — Dom et al. (2012) 自身が、先行4署名(いずれも腫瘍由来)をセットとして自らの被曝対非被曝の正常組織対比にかけ、いずれの濃縮も検出していない — 唯一の組織対応外の重なり(S100A10、方向逆転 <!-- N-54 --> — Dom の報告リストは全て被曝正常組織で上昇方向なので、被曝群低発現側の重なりは構成上必ず逆転になる)は重みを付けずに記す。本研究のどの主張もこの員数では動かない <!-- C-13 -->。

### Limitations(C-11)

<!-- C-11: A/B 非識別性の明文化。ほか: ウェット追試不能・REO の中間データ依存、被曝時年齢の群間比較不能(N-63) -->

Three limitations are structural. First, the design is cross-sectional over driver-conditioned strata: an exposure-associated expression difference conditioned on driver cannot, in principle, distinguish a shared trace that feeds RET oncogenesis from an association between trace type and driver selection, and every claim in this paper is worded neutrally between the two <!-- C-11 -->. Second, two practical constraints cap claim strength: no wet-laboratory replication is available, and the REO panel is defined on intermediate data products rather than raw counts alone. Third, disclosed asymmetries qualify generalization: age at exposure exists only for exposed cases <!-- N-63 -->, the exposed BRAF group is small (9 cases <!-- N-09 -->), tumor purity is a within-cohort relative quantity, and the single held-out calibration excess (B_Normal/Hallmark <!-- N-25 -->) travels with any reading of that cell. Age belongs to this third class rather than to the first two: because the unexposed cases of the source cohort were all born after the accident (Morton et al. 2021), exposure status is tied to birth cohort by construction, and a contribution of age — and of the latency inseparable from it — to the reported difference cannot be excluded by this design. It is disclosed (Table 1; Discussion) rather than adjusted for, for the reason given in Methods. Assigned share itself is a model-transfer quantity — its thyroid risk model is fitted to external photon exposure and applied here to internal iodine-131 (Methods) — and it is used only as an ordinal banding index, never as a calibrated probability.

【訳】3つの限界は構造的である。第一に、本デザインは driver 条件付き層の上の横断研究である: driver で条件付けた曝露関連の発現差は、RET 発がん機序に接続する共有痕跡と、痕跡タイプと driver 選択の連関とを原理的に識別できず、本論文の全ての主張は両者の間で中立に文言化されている <!-- C-11 -->。第二に、二つの実務的制約が主張の強さの上限を画す: ウェット実験による追試は利用できず、REO パネルは生カウント単独でなく中間データ産物の上に定義されている。第三に、開示済みの非対称が一般化を限定する: 被曝時年齢は被曝症例にのみ存在し <!-- N-63 -->、BRAF の被曝群は小さく(9 例 <!-- N-09 -->)、腫瘍純度はコホート内の相対量であり、held-out 較正の唯一の超過(B_Normal/Hallmark <!-- N-25 -->)は当該セルのどの読みにも付随する。年齢は最初の2種でなくこの第三の種類に属する: ソースコホートの非被曝例は全員が事故後の出生であるため(Morton et al. 2021)、曝露の身分は構成上、出生コホートと結びついている。報告された差への年齢 — および不可分の潜伏期 — の寄与は、本デザインでは排除できない。Methods に述べた理由により、調整でなく開示する(表1; Discussion)。assigned share 自体も転送モデルに基づく量である — その甲状腺リスクモデルは外部光子被曝に当てはめられたもので、ここでは内部ヨウ素131に適用している(Methods)— そして順序つき帯化指標としてのみ用い、較正された確率としては用いない。

## Figure legends and table captions

<!-- 2026-08-16 起草(判断点5)。研究者の吟味・書き直し前提の叩き台。表1キャプションの開示規律札(検定なし・p なし・群内ブートストラップ CI)は要件 -->

**Fig. 1 | Case flow of the analysis cohorts.** From the 440 cases of the REBC-THYR cohort <!-- N-11 --> to the main analysis cohort (n = 63 <!-- N-08 -->) and the REO evaluation set (n = 36 <!-- N-10 -->), stratified by driver (RET / BRAF). Boxes give per-step counts with per-step exclusions annotated; the side branch carries the paired RET tumors of the Low and Mid assigned-share bands into the out-of-sample REO evaluation.

**Fig. 2 | Gene-level evidence for the exposure contrast, per analysis contrast.** Signed Brunner–Munzel effect (P(X<Y) − P(X>Y); positive = higher in the exposed group) against −log10 exact permutation p, one panel per contrast. Genes at Storey q < 0.10 are colored by direction (R_Tumor: 1,765; B_Normal: 1; none elsewhere <!-- N-16 -->). The discovery call is rank-based; no fold change enters it.

**Fig. 3 | REO reversal scores against assigned share in RET tumors.** Per-case reversal score (0–10) of the 10-pair panel <!-- N-37 -->; training bands (Sporadic, High) open, out-of-sample evaluation bands (Low, Mid) filled. Assigned share is undefined for the unexposed Sporadic cases, drawn as a separate strip (stacked points count cases). Band medians 0 / 1 / 4 / 6 <!-- N-41 -->; the graded profile is reported descriptively and no dose–response form is assumed.

**Table 1 | Case characteristics of the analysis groups.** Main analysis cohort (R_Sporadic, R_High, B_Sporadic, B_High) and REO evaluation set (R_Low, R_Mid): group size, sex, age at surgery and age at exposure (median [range]), and driver composition (fusion partner within the RET stratum) <!-- N-12, N-13 -->. Footnotes give the between-group age-at-surgery differences (High − Sporadic) as Hodges–Lehmann median differences and the rank-based effect P(Sporadic < High), each with a 95% percentile bootstrap confidence interval (within-group resampling) <!-- N-64, N-65 -->. These are disclosures of magnitude and uncertainty: no hypothesis test is performed and no p-value is computed. Age at exposure is undefined for the unexposed Sporadic groups <!-- N-63 -->.

**Table 2 | Interpretation map for the four exposure contrasts, fixed before the finalized analysis produced the results reported here.** Each contrast carries a pre-specified standing and its basis; the pattern rules fixed together with the map are given as table footnotes. No outcome pattern requires post-hoc labelling.

<!-- 表2 本体+脚注規則(正本。Methods: Analysis contrasts から 2026-08-18 移設 — 投稿整形時に一方向切り出し、編集可能な第2コピーは作らない。規則の (a)-(f) 脚注化は研究者の書き直し待ち) -->
| Contrast | Pre-specified standing | Basis |
| --- | --- | --- |
| R_Tumor | Hypothesis-bearing; signal expected | primary expectation of the working hypothesis (Introduction); the driver-conditioned gap left by Morton et al. (2021) |
| R_Normal | Hypothesis-bearing; signal possible | long-term molecular memory of exposure — dose-dependent expression within exposed normal tissue (Abend et al. 2013); exposed-versus-unexposed normal-tissue signatures (Ory et al. 2026) |
| B_Normal | Test cell for a shared trace | counterpart of the cross-stratum normal comparison; the shared signatures of Ory et al. (2026) were reported without driver stratification |
| B_Tumor | Direction-agnostic; not excludable | point-mutation drivers become less likely as dose rises in the source cohort (Morton et al. 2021) |

Pattern rules, fixed with the map: signal in R_Tumor is read as agreement with prior expectation; tumor signal confined to the BRAF stratum would be hypothesis-discordant and read with non-aetiological explanations first; concordant signal in both normal contrasts is read as consistent with a shared glandular trace, discordant signal as suggesting a driver-linked trace type; normal-tissue signal confined to the BRAF stratum is read as confounding first; confined to the RET stratum, as not identifiable between attenuated sharing and a driver-linked type; and an all-null outcome is reported as a bounded null. No outcome pattern requires post-hoc labelling.

【訳】(表2本体: R_Tumor = 仮説担持・シグナル期待〔作業仮説の一次期待(Introduction); Morton et al. 2021 が残した driver 条件付きの隙間〕/ R_Normal = 仮説担持・シグナルあり得る〔被曝の長期分子記憶 — 被曝例内の線量依存発現(Abend et al. 2013)と被曝対非被曝の正常組織署名(Ory et al. 2026)〕/ B_Normal = 共有痕跡のテストセル〔層間正常比較の相手方; Ory et al. 2026 の共有署名は driver 非層別で報告された〕/ B_Tumor = 方向不可知・排除不能〔ソースコホートでは線量が上がるほど点変異 driver は生じにくい(Morton et al. 2021)〕)。マップと同時に固定したパターン規則: R_Tumor のシグナルは事前期待との一致として読む。BRAF 層に限局した腫瘍シグナルは仮説不協和であり、非病因的説明を第一として読む。両正常対比の一致シグナルは共有腺痕跡と整合、不一致シグナルは driver 連関の痕跡タイプの示唆として読む。BRAF 層に限局した正常組織シグナルは交絡第一で、RET 層に限局した場合は減衰した共有と driver 連関タイプの間で識別不能として読む。全対比 null は bounded null として報告する。どの結果パターンも事後のラベリングを要しない。

**Table 3 | Gene-level exposure contrast, per analysis contrast.** Tested genes after filtering, Storey π0 (plug-in, λ = 0.5), discoveries at q < 0.10 with direction (higher / lower in the exposed group), minimum exact permutation p, and the pre-specified primary omnibus (Higher Criticism) p from the contrast's own label-permutation null <!-- N-15, N-16, N-17, N-18, N-19, N-20 -->.

【訳】図1 | 解析コホートの症例フロー。REBC-THYR コホートの 440 症例 <!-- N-11 --> から、主解析コホート(n = 63 <!-- N-08 -->)と REO 評価セット(n = 36 <!-- N-10 -->)まで、driver(RET / BRAF)で層別して示す。各箱は段階別の症例数で、段階別の除外は注記に示す。側枝は Low・Mid 帯のペア付き RET 腫瘍を out-of-sample の REO 評価へ運ぶ。/ 図2 | 解析対比別の遺伝子レベル証拠。符号付き Brunner–Munzel 効果(P(X<Y) − P(X>Y); 正 = 被曝群で高発現)を −log10(exact 置換 p)に対し、対比ごとに1パネルで示す。Storey q < 0.10 の遺伝子を方向別に着色(R_Tumor: 1,765・B_Normal: 1・他はなし <!-- N-16 -->)。発見の判定は順位ベースであり、fold change は入らない。/ 図3 | RET 腫瘍における REO 逆転スコアと assigned share。10 ペアパネル <!-- N-37 --> の症例別逆転スコア(0–10)。訓練帯(Sporadic・High)は白抜き、out-of-sample 評価帯(Low・Mid)は塗り。assigned share は非被曝の Sporadic には定義されず、専用ストリップに描く(積み上げ点は員数)。帯別中央値は 0 / 1 / 4 / 6 <!-- N-41 -->。段階的プロファイルは記述的に報告し、線量反応の形は仮定しない。/ 表1 | 解析群の症例特性。主解析コホート(R_Sporadic・R_High・B_Sporadic・B_High)と REO 評価セット(R_Low・R_Mid)について、群サイズ・性・手術時年齢・被曝時年齢(中央値[範囲])・driver 構成(RET 層は融合パートナー別)を示す <!-- N-12, N-13 -->。脚注に群間の手術時年齢差(High − Sporadic)を Hodges–Lehmann 中央値差と順位効果 P(Sporadic < High) で、それぞれ 95% percentile ブートストラップ信頼区間(群内復元抽出)つきで示す <!-- N-64, N-65 -->。これらは大きさと不確かさの開示であり、仮説検定は行わず p 値は計算しない。被曝時年齢は非被曝の Sporadic 群には定義されない <!-- N-63 -->。/ 表2 | 4つの曝露対比の解釈マップ(確定解析が本稿で報告する結果を生成する前に固定)。各対比は事前指定の位置づけとその根拠を持ち、マップと同時に固定したパターン規則を表脚注に示す。どの結果パターンも事後のラベリングを要しない。/ 表3 | 解析対比別の遺伝子レベル曝露対比。フィルタ後の検定遺伝子数、Storey π0(plug-in、λ = 0.5)、q < 0.10 の発見数と方向(被曝群で高発現/低発現)、最小 exact 置換 p、および当該対比自身のラベル置換帰無から得た事前指定主オムニバス(Higher Criticism)の p <!-- N-15, N-16, N-17, N-18, N-19, N-20 -->。

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
