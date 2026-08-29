# Cover letter(BJC 投稿用 — 起草 2026-08-29、研究者確認待ち)

状態: draft。BJC GTA の必須要素(BJC への投稿・未発表・他誌未投稿の言明/対応著者の所属と連絡先/競合利益の言明/重要性と読者層への適合の簡潔な説明)を満たす構成。英語が正本。数値は本文と同一(N-16・N-20・N-08・N-41)。
研究者確認事項: (i) 日付、(ii) プレプリントの有無(現状「未掲載」で起草 — 掲載する場合は該当文を差し替え)、(iii) 署名の形、(iv) BJC 既刊への言及段落(Detours 2005/2007・Dom 2012・Abend 2013 — PubMed で誌名・年を照合済み 2026-08-29)と設計の転用可能性の段落の強度。Thyroid への旧投稿は BJC 規定上の開示義務なし(内容も解析も別物のため言及しない)。

---

[Date]

The Editors  
British Journal of Cancer

Dear Editors,

Please consider the enclosed manuscript, "Driver-conditioned transcriptomic differences across radiation-attributability bands in papillary thyroid carcinoma", for publication in the *British Journal of Cancer* as an original research article.

Papillary thyroid carcinoma after the Chernobyl accident remains a major human model of radiation carcinogenesis, yet whether exposure is associated with tumor expression *within* a driver stratum has been unresolved: driver composition itself changes with dose, so pooled exposed-versus-unexposed comparisons can reproduce the difference between fusion-driven and point-mutation-driven tumors rather than an exposure-associated difference. Using the REBC-THYR cohort (440 cases; open-access RNA-seq from the NCI Genomic Data Commons), we compared cases with high estimated radiation attributability (NIH IREP Assigned Share) with dose-zero cases separately within the RET fusion-positive and BRAF V600E-positive strata, under an interpretation map fixed before the results were produced. This extreme-band design reflects a mixture model in which dose changes the proportion of tumors arising through a radiation-associated initiating event, without requiring expression in every tumor to vary monotonically with dose; the source cohort's continuous dose–expression analyses had not identified a corresponding signal. RET fusion-positive tumors from High-AS cases differed broadly in expression from dose-zero tumors (1,765 genes at Storey q<0.10; Higher Criticism p=0.0112), with no corresponding evidence in the other three contrasts. Under the gene-sampling reference used for over-representation analysis, genes expressed less highly in High-AS tumors were concentrated in proliferation, cell-cycle and DNA-repair annotations — a coherent direction for mechanistic investigation, although no set crossed q<0.10 under the calibrated label-permutation analysis — and the RET and BRAF tumor profiles were more concordant than their shuffled reference, a hypothesis-generating observation that motivates a cross-driver hypothesis rather than demonstrating an oncogene-independent radiation trace. A relative-expression-ordering panel built from the extreme bands and applied unchanged to the intermediate bands gave ordered median scores across the four attributability bands. The design supports an association between attributability band and expression within a driver stratum; it does not establish a radiation-specific signature, and the panel requires independent validation. All inference used exact enumeration or permutation calibration, and the complete pipeline, run independently on two machines in a date-pinned container, produced identical primary artifacts; the analysis code, versioned inputs and container recipe accompany the paper.

The manuscript continues a line of work that the *British Journal of Cancer* has carried for two decades. The question of whether post-Chernobyl thyroid cancers bear a radiation-specific expression signature was raised in this journal by Detours et al. (2005, 2007), and two of the studies against which we cross-reference our results — the normal-tissue signature of Dom et al. (2012) and the iodine-131 dose-dependent expression changes of Abend et al. (2013) — were also published here. What has changed since is the cohort-scale finding that driver composition itself varies with dose, which makes the driver-conditioned form of the question the necessary next step; this manuscript takes that step in the same cohort that established the finding.

Beyond thyroid cancer, the design — comparison of the extreme bands of a per-case attributability metric within driver-defined tumor classes, under an interpretation map fixed in advance, with exact enumeration and permutation calibration at group sizes of ten to thirty — is transferable to other settings in which an exposure metric is available and oncogenic drivers structure expression. We therefore expect the results and the framework to interest the journal's readership in cancer epidemiology, cancer genomics and biostatistics as well as in thyroid oncology.

This manuscript is submitted to the *British Journal of Cancer* only. It reports original research that has not been published previously and has not been posted as a preprint, and it is not under consideration for publication elsewhere. All authors have approved the manuscript and agree to its submission. The authors declare no competing interests, and the authors received no specific funding for this work. Generative AI tools assisted with code development and text editing under author supervision, as documented in the Methods.

Corresponding author: Kotaro Harakawa, Department of Molecular Oncology and Diagnostic Medicine, Atomic Bomb Disease Institute, Nagasaki University, 1-12-4 Sakamoto, Nagasaki 852-8523, Japan; e-mail ktrh61@gmail.com; ORCID 0009-0004-1086-8046.

Yours sincerely,

Kotaro Harakawa
on behalf of the authors (Vladimir A. Saenko, Norisato Mitsutake)
