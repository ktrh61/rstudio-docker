# 原稿統合ドラフト(draft_manuscript — 論文へ収束する単一ファイル)

状態: draft(2026-08-14 統合。研究者指示: Methods/Results の対構造の検査可能化と同期機構の廃止のため、旧4ドラフト — draft_methods / draft_results / draft_methods_results_d6 / draft_ora_annotation — を本ファイルへ吸収し削除。内容の変更なし、移設+Discussion 骨組み追加のみ。履歴は git)。規約: 英語本文が正本。Results の【訳】は研究者確認用の対訳(正本でない)。数値は N-ID タグ、主張段落は C-ID タグ(検査対象は英語本文)。図表番号は全て(仮)。各節の由来・決定は末尾の「起草メモ(統合)」。

---

## Introduction

<!-- 2026-08-15 フルスクラッチ起草(研究者指示)。仮説スロット1箇所は研究者記入待ち。予告列はここで確定(claim_map に記録)。マップ実体は Methods の小表(判断点5 の片割れを決着) -->

Papillary thyroid carcinoma after the Chernobyl accident remains the most extensively characterized human model of radiation carcinogenesis. Its molecular hallmark was established early: radiation-related tumors carry RET fusions at high frequency (Nikiforov et al. 1997; Rabes et al. 2000), and the genomic landscape of the REBC-THYR cohort confirmed the pattern at scale — fusion drivers, but not point mutations, associate with radiation dose (Morton et al. 2021). The driver event, however, is where the molecular account largely stops. Whether radiation exposure leaves a detectable trace in the expression program of the resulting tumor — beyond determining which driver initiates it — is unresolved.

Two features of the data make this question harder than it appears. First, tumor transcriptomes are organized by their driver: an exposed-versus-sporadic comparison that pools drivers confounds exposure with driver composition, because driver prevalence itself varies with dose (Morton et al. 2021). Any claim about an exposure trace must therefore be made within driver-conditioned comparisons. Second, candidate expression traces of exposure have been reported in both tissue compartments: dose-dependent gene expression in tumors and in corresponding normal tissue (Abend et al. 2012; 2013), the tumor-tissue radiation marker CLIP2 (Hess et al. 2011), and, recently, shared normal-tissue signatures in exposed cases reported without driver stratification (Ory et al. 2026). None of these observations, however, was made under driver conditioning; whether any exposure trace survives it has not been tested.

<!-- 仮説段落: 2026-08-15 研究者承認の叩き台(書き直し前提)。バイナリ観+Morton POC フック — 線量の役割は起始確率(Morton の線形性を包含)、否定するのは発現の連続線量応答のみ -->
These genomic findings suggest that dose is the wrong axis for an expression trace. Radiation acts by initiating: dose raises the probability that a tumor arises from a radiation-induced driver event, linearly in the observed range (Morton et al. 2021), so an individual tumor is either radiogenic or sporadic, and dose moves the composition of a case series, not the expression of any single tumor. Any exposure trace should therefore segregate with case status — and be diluted along a dose axis, which is how the source cohort's own dose regressions of expression came up empty. We instead work on the per-case scale: radiation attributability, the probability-of-causation quantity into which Morton et al. (2021) themselves proposed folding radiation-associated characteristics, computed as the NIH IREP Assigned Share (Methods) and used to band cases within each driver stratum.

The design follows from two constraints. Driver conditioning is mandatory — tumor expression is organized by the driver — and it is feasible only where a stratum supplies both exposed and unexposed cases in numbers that support inference; in this cohort, only the RET-fusion and BRAF V600E strata do. The same dose–fusion association that motivates the study ranks the strata: the exposed RET bands carry the highest expected radiogenic fraction, while BRAF V600E, with no dose association, provides the low-prior comparison. We therefore test four pre-specified contrasts — exposed versus sporadic, in tumor and normal tissue, within each stratum — under an interpretation map fixed before the reported results existed (Methods, Table 5(仮)): the primary expectation falls on RET tumors, the normal-tissue contrasts carry the exposure-memory hypothesis, and their cross-stratum comparison makes the shared-memory reports testable under driver conditioning. A rank-based REO panel adds a per-case readout expected to grade with assigned share across the intermediate bands, and the set-level machinery is calibrated before touching the real contrasts (Methods). This paper reports the outcome of that program.

【訳】チェルノブイリ事故後の乳頭状甲状腺癌は、ヒト放射線発がんの最も詳細に特徴づけられたモデルであり続けている。その分子的特徴は早期に確立された: 放射線関連腫瘍は高頻度に RET 融合を持ち(Nikiforov et al. 1997; Rabes et al. 2000)、REBC-THYR コホートのゲノムランドスケープはこのパターンを大規模に確認した — 融合 driver は線量と関連し、点変異は関連しない(Morton et al. 2021)。しかし分子的な説明は driver 事象でほぼ止まっている。放射線曝露が — どの driver で始まるかを決めることを越えて — 結果として生じた腫瘍の発現プログラムに検出可能な痕跡を残すのかは、未解決である。

データの2つの性質が、この問いを見かけより難しくしている。第一に、腫瘍トランスクリプトームは driver によって組織化される: driver をプールした被曝対散発比較は、driver の頻度自体が線量とともに変わるため(Morton et al. 2021)、曝露と driver 構成を交絡させる。したがって曝露痕跡についてのどの主張も、driver 条件付き比較の内側で立てる必要がある。第二に、曝露の発現痕跡の候補は両組織区画で報告されてきた: 腫瘍と対応正常組織における線量依存的遺伝子発現(Abend et al. 2012; 2013)、腫瘍組織の放射線マーカー CLIP2(Hess et al. 2011)、そして近年、driver 非層別で報告された被曝例の正常組織共有署名(Ory et al. 2026)。しかしこれらの観察はいずれも driver 条件付けの下で行われておらず、いずれかの曝露痕跡がそれに耐えるのかは検定されてこなかった。

これらのゲノム所見は、線量が発現痕跡の軸として誤りであることを示唆する。放射線は起始によって作用する: 線量は、腫瘍が放射線誘発の driver 事象から生じる確率を観察範囲では線形に引き上げる(Morton et al. 2021)。したがって個々の腫瘍は radiogenic か sporadic かのいずれかであり、線量が動かすのは症例系列の構成であって、個々の腫瘍の発現ではない。曝露痕跡があるとすれば、それは症例の身分と共に分離するはずであり — 線量軸に沿っては希釈される。ソースコホート自身の発現の線量回帰が空振りに終わったのは、まさにこの形である。我々は代わりに症例単位の尺度で作業する: 放射線起因性 — Morton et al. (2021) 自身が放射線関連特性の統合先として挙げた原因確率の量 — を NIH IREP の Assigned Share として計算し(Methods)、各 driver 層内の症例の帯化に用いる。

設計は2つの制約から従う。driver 条件付けは必須であり — 腫瘍発現は driver によって組織化される — それが可能なのは、層が被曝例と非被曝例の両方を推論に耐える数で供給する場合に限られる。このコホートでそれを満たすのは RET 融合と BRAF V600E の層だけである。研究を動機づけたのと同じ線量–融合連関が層の序列も与える: RET の被曝帯は最も高い radiogenic 割合を期待させ、線量関連を持たない BRAF V600E は低事前確率の比較側となる。我々はしたがって、事前指定の4対比 — 各層内の、腫瘍と正常組織での被曝対散発 — を、報告結果が存在する前に固定した解釈マップの下で検定する(Methods、表5(仮)): 一次期待は RET 腫瘍に落ち、正常組織対比が曝露記憶仮説を担い、その層間比較が共有記憶の報告を driver 条件付けの下で検定可能にする。順位ベースの REO パネルが、中間帯で assigned share と共に段階化すると期待される症例別の読み出しを加え、セットレベル機構は実対比に触れる前に較正される(Methods)。本論文はこのプログラムの帰結を報告する。

## Methods

### Data sources and expression matrix

<!-- 写像元: 010, 020, 030, 110, 120 -->

Gene-level RNA-seq counts (STAR - Counts, open access) for the REBC-THYR project were downloaded from the NCI Genomic Data Commons and md5-verified against the download manifest. Gene lengths (exon-union length per gene) were derived from the GDC GENCODE v36 reference annotation <!-- N-66 -->. The clinical table is Data S1 of the REBC-THYR publication (440 cases
<!-- N-11 -->), read without dropping columns or editing values; missing markers were coerced to NA and columns typed numerically only where every non-missing value parses. Expression files were mapped to cases and biospecimen samples through the GDC API, and a single count assay per sample was assembled into one expression matrix (58,448 genes × 906 samples after dropping all-zero genes <!-- N-14 -->). Library strandedness was decided per sample from the STAR
stranded-column totals (ratio rule: the library is called stranded when the smaller stranded total is at most half the larger, and the larger column is used; otherwise the unstranded column is used <!-- N-67 -->). In practice the call was unambiguous — all 906 libraries were stranded (reverse), with ratios of 0.056–0.110 against the 0.5 threshold <!-- N-77 -->. TPM and FPKM columns were discarded at loading.

【訳】REBC-THYR プロジェクトの遺伝子レベル RNA-seq カウント(STAR - Counts、open access)を NCI Genomic Data Commons からダウンロードし、ダウンロードマニフェストに対して md5 検証した。遺伝子長(遺伝子ごとの exon-union 長)は GDC の GENCODE v36 参照アノテーションから導出した <!-- N-66 -->。臨床表は REBC-THYR 原論文の Data S1 (440 症例 <!-- N-11 -->)で、列の削除も値の編集も行わずに読み込み、欠損マーカーは NA に統一し、非欠損の全値が数値として解釈できる列のみ数値型とした。発現ファイルは GDC API で症例・生体検体に対応づけ、検体あたり単一のカウントアッセイを1つの発現行列に組み上げた(全ゼロ遺伝子の除去後で 58,448 遺伝子 × 906 検体 <!-- N-14 -->)。ライブラリの strand は検体ごとに STAR の stranded 2列の合計から判定した(比率規則: 小さい方の合計が大きい方の半分以下なら stranded と判定して大きい列を採用、それ以外は unstranded 列を採用 <!-- N-67 -->)。実際には判定は一義的で、全 906 ライブラリが stranded(reverse)判定、比は 0.056–0.110 と閾値 0.5 から大きく離れていた <!-- N-77 -->。TPM・FPKM 列は読み込み時に破棄した。

### Exposure metric: assigned share

<!-- 写像元: 130, 140 -->

Radiation attributability per exposed case was quantified with the thyroid model of the National Institutes of Health Interactive RadioEpidemiological Program (NIH IREP, version 5.7.3; Kocher et al. 2008) as the Assigned Share associated with the expected value of the excess relative risk (AS), computed from each case's recorded thyroid dose and ages with a fixed input convention (electrons E > 15 keV; sex as recorded; dose in cSv = mGy/10; exposure year 1986; birth year = 1986 − age at exposure; surgery year = birth year + age at surgery), with all IREP assumptions and settings at their defaults (user-defined uncertainty distribution Lognormal (1,1); 10,000 iterations; random number seed 99) <!-- N-68 -->. The AS values are carried as a versioned input file with the analysis code; their dose and age inputs were verified case-by-case against the clinical table.

Cases were banded on AS with a single pre-fixed rule: unexposed cases (dose 0) form the Sporadic band; exposed cases fall in Low (0 < AS < 33.3), Mid (33.3 ≤ AS < 66.6), or High (AS ≥ 66.6); no case sits on a boundary
<!-- N-69 -->. The clinical composition of all 440 cases across driver × band × pairing is given in Table 1(仮) <!-- N-11 -->.

【訳】被曝症例ごとの放射線起因性は、NIH の Interactive RadioEpidemiological Program(NIH IREP、バージョン 5.7.3; Kocher et al. 2008)の甲状腺モデルにより、過剰相対リスク(ERR)の期待値に対応する Assigned Share(AS)として定量した。各症例の記録甲状腺線量と年齢から、固定した入力規約(電子 E > 15 keV; 性別は記録値のまま; 線量は cSv = mGy/10;被曝年 1986; 出生年 = 1986 − 被曝時年齢; 手術年 = 出生年 + 手術時年齢)で計算し、IREP の仮定・設定はすべて既定値のまま用いた(ユーザー定義不確かさ分布 Lognormal (1,1)、反復 10,000 回、乱数シード 99)
<!-- N-68 -->。AS 値は解析コードとともに版管理された入力ファイルとして持ち回り、線量・年齢の入力は臨床表と全例照合した。症例は単一の事前固定規則で AS 帯に割り付けた: 非被曝(線量 0)は Sporadic 帯、被曝症例は Low(0 < AS < 33.3)・Mid(33.3 ≤ AS < 66.6)・High(AS ≥ 66.6)。境界上の症例は存在しない <!-- N-69 -->。440 症例全体の driver × 帯 × ペア有無の臨床構成は
表1(仮)に示す <!-- N-11 -->。

### Quality control and analysis cohorts

<!-- 写像元: 210, 220, 230 -->

All eligibility conditions were finalized into explicit per-case flags before any inference, and every downstream script consumes the same flags. Sample outliers were screened before purity estimation, so that an anomalous sample cannot contaminate the purity fit, with PC-scores-based outlier detection (PC-OD; Nakayama et al. 2024; in-house implementation): within each group × tissue sub-matrix (filterByExpr-reduced, Chen et al. 2016; unnormalized log-CPM, since raw library sizes preserve the composition outliers normalization would mask), the most extreme sample on the first principal component was removed iteratively (two-sided, α = 0.05) until no sample was rejected <!-- N-78 -->. Relative tumor purity was then estimated per case with the contamDE method (Shen et al. 2016; in-house implementation) on MUREN-normalized (Feng and Li 2021) paired counts (protein-coding genes, filterByExpr), in one run per driver cohort with both groups (Sporadic and High) together: contamDE purities are relative within the set of pairs estimated jointly, so a pooled run is what gives both groups one common scale — reported as a maximum-one relative score within each cohort — on which a single threshold has a single meaning. The common tumor reference this assumes rests on the same premise as the driver stratification itself — tumor expression in both groups is dominated by the driver's biology — and the tumor-versus-normal contamination axis that purity measures is a different axis from the exposure contrast under test; the estimate serves only as a within-cohort relative filter and a diagnostic covariate.

Case pairing resolves, per case, the GDC merged-aliquot expression assays (sample identifiers suffixed _merged; one merged tumor and one merged normal, enforced unique) — a case lacking either is unpaired. The main analysis cohort comprises driver-classified (RET or BRAF) cases in the Sporadic or High band with a paired tumor/normal sample, both tissues outlier-clean, and pooled relative purity ≥ 0.6 <!-- N-70 --> (n = 63 <!-- N-08 -->); the realized case flow and group composition are reported in Results (cohort-flow figure(仮); Tables 1–2(仮)). The REO training set is the RET subset of this cohort, and the REO evaluation set is the paired RET tumors of the Low and Mid bands, deliberately left unfiltered by the outlier and purity conditions (its own QC is reported as diagnostics, below).

【訳】全ての適格条件は推論に先立って症例別の明示的フラグとして確定し、下流の全スクリプトが同一のフラグを読む。検体外れ値は純度推定より前に(異常検体が純度の当てはめを汚染しないよう)、PC スコアに基づく外れ値検出(PC-OD; Nakayama et al. 2024; 自作実装)でスクリーニングした: 群 × 組織の各部分行列(filterByExpr 縮約、Chen et al. 2016; 未正規化 log-CPM — 生のライブラリサイズは正規化が隠す組成外れ値を保存する)で、第1主成分上で最も極端な検体を反復的に除去し(両側、α = 0.05)、棄却が出なくなるまで続けた <!-- N-78 -->。相対腫瘍純度はその後、MUREN 正規化(Feng and Li 2021)したペア付きカウント(protein-coding、filterByExpr)に対する contamDE 法(Shen et al. 2016; 自作実装)で症例ごとに推定した — driver コホートごとに両群(Sporadic・High)を合わせた1回の実行である: contamDE の純度は同時に推定した集合内でのみ相対的なので、プールした実行こそが両群に単一の共通尺度(コホート内最大値 = 1 の相対スコアとして報告; 単一の閾値が単一の意味を持つ尺度)を与える。この仮定する共通腫瘍参照は driver 層別という設計自体の前提(両群の腫瘍発現は driver の生物学が支配する)に乗っており、純度が測る腫瘍対正常の混入軸は検定対象の曝露対比とは別の軸である。推定値の役割はコホート内の相対フィルタと診断用共変量のみである。症例のペア解決は GDC の merged-aliquot 発現アッセイ(検体 ID 接尾辞 _merged; 症例ごとに merged 腫瘍1・merged 正常1、一意性を強制)を用い、どちらかを欠く症例はペアなしとする。主解析コホートは、driver 分類済み(RET または BRAF)で Sporadic か High 帯、腫瘍/正常ペアを持ち、両組織とも外れ値クリーン、かつプール相対純度 ≥ 0.6
<!-- N-70 --> の症例からなる(n = 63 <!-- N-08 -->)。実現した症例フローと群構成は Results に報告する(コホートフロー図(仮); 表1–2(仮))。REO 訓練セットはこのコホートの RET 部分集合、REO 評価セットは Low・Mid 帯のペア付き RET 腫瘍で、意図的に外れ値・純度条件で濾過していない(その QC は後述の
診断として報告)。

### Analysis contrasts

<!-- 写像元: 310 header, lib/units.R -->

All inference is organized in four analysis contrasts — exposed versus sporadic within each driver stratum and tissue: R_Tumor, R_Normal (RET stratum), B_Tumor, B_Normal (BRAF stratum). Each contrast is tested within itself; no cross-contrast family-wise inference is performed, and no study-wide FDR is claimed. In every contrast, group x is Sporadic and group y is High, so effects > 0.5 indicate higher expression in the exposed group. Throughout, two-sample comparison uses the Brunner–Munzel statistic (Brunner and Munzel 2000): the protocol takes the order relation as its primitive — a commitment that recurs in the rank-based REO panel and the tie-invariant enrichment statistic below. Throughout, discovery thresholds sit at the 0.10 false-discovery level; the estimator behind each threshold is stated where that level of inference is introduced.

【訳】全ての推論は4つの解析対比 — 各 driver 層 × 組織内の被曝対散発対比: R_Tumor・R_Normal(RET 層)、B_Tumor・B_Normal(BRAF 層)— に組織化される。各対比は対比内で検定され、対比横断の family-wise 推論は行わず、研究全体の FDR は主張しない。全対比で群 x は Sporadic、群 y は High であり、効果 > 0.5 は被曝群での高発現を意味する。全体を通じて、2標本比較には Brunner–Munzel 統計量(Brunner and Munzel 2000)を用いる: 本プロトコルは順序関係を原始に取る — このコミットメントは後述の順位ベース REO パネルと tie 不変の濃縮統計量にも回帰する。全体を通じて、発見の閾値は誤発見水準 0.10 に置く。閾値の背後の推定量は、各推論水準の導入箇所で述べる。

<!-- 解釈マップの実体(計画v2 §0.6 批准マップの転記 — 判断点5 の置き場所決定: Methods 小表)。C-01/C-04/C-05/C-07/C-16 の読みの参照先 -->
The interpretation map fixed for the four contrasts is given in Table 5(仮); readings for every outcome pattern were assigned before the reported results existed.

| Contrast | Pre-specified standing | Basis |
|---|---|---|
| R_Tumor | Hypothesis-bearing; signal expected | primary expectation of the working hypothesis (Introduction); the driver-conditioned gap left by Morton et al. (2021) |
| R_Normal | Hypothesis-bearing; signal possible | long-term molecular memory of exposure (Abend et al. 2013; Ory et al. 2026) |
| B_Normal | Test cell for a shared trace | counterpart of the cross-stratum normal comparison; the shared signatures of Ory et al. (2026) were reported without driver stratification |
| B_Tumor | Direction-agnostic; not excludable | no dose association of BRAF V600E in the source cohort (Morton et al. 2021) |

Pattern rules, fixed with the map: signal in R_Tumor is read as agreement with prior expectation; tumor signal confined to the BRAF stratum would be hypothesis-discordant and read with non-aetiological explanations first; concordant signal in both normal contrasts is read as consistent with a shared glandular trace, discordant signal as suggesting a driver-linked trace type; normal-tissue signal confined to the BRAF stratum is read as confounding first; confined to the RET stratum, as not identifiable between attenuated sharing and a driver-linked type; and an all-null outcome is reported as a bounded null. No outcome pattern requires post-hoc labelling.

【訳】4対比について固定した解釈マップを表5(仮)に示す。全ての結果パターンへの読みは、報告結果が存在する前に割り当てられた。(表: R_Tumor = 仮説担持・シグナル期待〔作業仮説の一次期待(Introduction); Morton et al. 2021 が残した driver 条件付きの隙間〕/ R_Normal = 仮説担持・シグナルあり得る〔被曝の長期分子記憶の報告(Abend et al. 2013; Ory et al. 2026)〕/ B_Normal = 共有痕跡のテストセル〔層間正常比較の相手方; Ory et al. 2026 の共有署名は driver 非層別で報告された〕/ B_Tumor = 方向不可知・排除不能〔ソースコホートで BRAF V600E に線量関連なし(Morton et al. 2021)〕)。マップと同時に固定したパターン規則: R_Tumor のシグナルは事前期待との一致として読む。BRAF 層に限局した腫瘍シグナルは仮説不協和であり、非病因的説明を第一として読む。両正常対比の一致シグナルは共有腺痕跡と整合、不一致シグナルは driver 連関の痕跡タイプの示唆として読む。BRAF 層に限局した正常組織シグナルは交絡第一で、RET 層に限局した場合は減衰した共有と driver 連関タイプの間で識別不能として読む。全対比 null は bounded null として報告する。どの結果パターンも事後のラベリングを要しない。

### Candidate confounders

<!-- 写像元: diagnostics/age_arm_difference.R(年齢推定)、140/Table 2(組成)。性連鎖注釈は実行 Go 待ち -->

<!-- C-15 / C-10(年齢の共変量非投入) -->
Three candidate confounders of the exposure contrasts follow from the design and are disclosed rather than adjusted for; none is tested. The first is age: assigned share is a function of age at exposure, and with a single exposure event the banding itself induces age structure across the bands, so this candidate is guaranteed by construction. The between-group difference in age at surgery (High − Sporadic) is therefore estimated in each driver stratum by the Hodges–Lehmann median difference and the rank-based effect P(Sporadic < High) — the protocol's Brunner–Munzel effect estimator — each with a 95% percentile bootstrap confidence interval (within-group resampling, 9,999 replicates, seed 19450809 <!-- N-63 -->); this is a disclosure of magnitude and uncertainty, not a confounding test, and no p-value is computed. Age at exposure is structurally not comparable between groups (the Sporadic group is unexposed) <!-- N-63 -->. Age is not entered as a covariate, because age at exposure is a component of the AS definition and adjusting for it would remove part of the exposure metric itself; group age structure is reported descriptively (per-group median and range, Table 2(仮) <!-- N-12, N-13 -->).

The second is sex: the assigned-share values are numerically insensitive to the recorded sex input, so sex is not absorbed into the exposure metric, and sex-chromosome genes differ mechanically between groups of unequal sex composition; the per-group composition is therefore reported in full (Table 2(仮) <!-- N-12, N-13 -->), and discovered genes are annotated for sex-chromosome membership. The third is oncogene subtype: the literature associates radiation-related papillary carcinoma with particular RET fusion partners (NCOA4-RET, historically RET/PTC3; Nikiforov et al. 1997; Rabes et al. 2000), so the fusion-partner composition is reported per band (Table 2(仮) <!-- N-12 -->); the BRAF stratum carries a single mutation by construction.

【訳】曝露対比の候補交絡因子3つは設計から導かれ、調整でなく開示する。いずれも検定しない。第一は年齢である: assigned share は被曝時年齢の関数であり、曝露が単一事象であるため、帯化そのものが帯間に年齢構造を誘導する — この候補は構成により保証されている。したがって手術時年齢の群間差(High − Sporadic)を driver 層ごとに、Hodges–Lehmann 中央値差と順位ベースの効果 P(Sporadic < High)(本プロトコルの Brunner–Munzel 効果推定量)で推定し、それぞれに 95% percentile ブートストラップ信頼区間を付す(群内復元抽出、9,999 レプリケート、シード 19450809 <!-- N-63 -->)。これは大きさと不確かさの開示であって交絡の検定ではなく、p 値は計算しない。被曝時年齢は構造的に群間比較が成立しない(Sporadic 群は非被曝)<!-- N-63 -->。年齢は共変量に入れない — 被曝時年齢は AS の定義の構成要素であり、調整は曝露指標そのものの一部除去になるためである。群の年齢構造は記述的に報告する(群別中央値と範囲、表2(仮)<!-- N-12, N-13 -->)。第二は性である: assigned share の値は入力した性別に対して数値的に不感応であり、性は曝露指標に吸収されない。性構成が異なる群間では性染色体遺伝子が機械的に発現差を生むため、群別の性構成を全数報告し(表2(仮)<!-- N-12, N-13 -->)、発見遺伝子には性染色体所属の注釈を付す。第三は oncogene サブタイプである: 文献は放射線関連乳頭癌を特定の RET 融合パートナー(NCOA4-RET — 歴史的呼称 RET/PTC3; Nikiforov et al. 1997; Rabes et al. 2000)と関連づけるため、融合パートナー構成を帯別に報告する(表2(仮)<!-- N-12 -->)。BRAF 層は構成により単一変異である。
<!-- C-15 -->

### Normalization

<!-- 写像元: 310, lib/norm_deges.R -->

Each contrast's count matrix (protein-coding genes, filterByExpr; edgeR) was normalized with a DEGES scheme (Kadota et al. 2012; iterated form of Sun et al. 2013) reimplemented in-house: MUREN normalization (pairwise least-trimmed-squares regression mode) alternating with an exactly enumerated permutation Brunner–Munzel screen — substituted for the original scheme's model-based DEG screen — that removes potential DEGs from the scaling-factor estimation (Storey q < 0.10, Storey 2002; the same estimator as the gene-level inference below, with the TCC floorPDEG guard: the larger of the q-threshold set and the top 5% of genes by raw p is removed), iterated three times (iDEGES) <!-- N-71 -->. The per-contrast tested-gene counts, screen pi0 estimates, iteration convergence (Jaccard), and resulting scaling-factor ranges are reported in Supplementary(仮) <!-- N-15 -->.

【訳】各対比のカウント行列(protein-coding、filterByExpr; edgeR)は、自作再実装の DEGES 方式(Kadota et al. 2012; 反復形は Sun et al. 2013)で正規化した: MUREN 正規化(pairwise least-trimmed-squares 回帰モード)と、潜在的 DEG をスケーリング係数推定から除去する完全枚挙の置換 Brunner–Munzel スクリーン(原方式のモデルベース DEG スクリーンを置換; Storey q < 0.10、Storey 2002 — 後述の遺伝子レベル推論と同一の推定量 — に TCC の floorPDEG ガードを併用し、q 閾値集合と生 p 上位 5% の大きい方を除去)との交互反復を3回(iDEGES)<!-- N-71 -->。対比別の検定遺伝子数、スクリーンの pi0 推定値、反復収束(Jaccard)、スケーリング係数の範囲は Supplementary(仮)に報告する <!-- N-15 -->。

### Gene-level differential expression

<!-- 写像元: 410(遺伝子別検定・Storey・シャッフルセット定義・pi0 不確かさ) -->

Each gene was tested for an exposed-vs-sporadic difference with the exact permutation Brunner–Munzel test (two-sided; exact enumeration of all C(n, nx) allocations, so the gene-level p-values need no seed). The implementation, provided with the code, is an in-house permutation enumeration run here in its exact mode (its Monte Carlo mode enters only the REO evaluation, below) — no asymptotic approximation enters any p-value. The effect P(X<Y) reads directly as an exceedance probability, and at these group sizes the permutation distribution is exactly enumerable, so gene-level inference carries neither a count-model assumption nor Monte Carlo error. Gene-level inference is the Storey q-value on the exact p-values with the plug-in pi0 estimate at λ = 0.5, thresholded at the protocol-wide 0.10 level <!-- N-05 -->. The estimator was fixed a priori from the working hypothesis and the design rather than tuned: the hypothesis is a weak signal spread across many genes — small per-gene effects attenuated by within-group case mixture — so holding π0 at 1 would build the absence of exactly that signal into the correction; the conservativeness of the fixed-λ plug-in requires only marginal uniformity of the null p-values, which the exact test guarantees regardless of the dependence between genes, whereas adaptive choices of λ carry guarantees that assume independent or weakly dependent tests and do not apply to co-expressed genes sharing one label vector; and under a weak spread signal the alternative p-density is nearly flat, so raising λ buys little bias reduction at first-order variance cost —λ = 0.5 is the untuned Storey (2002) default. Contrast-level null distributions — for the pi0 uncertainty below, the omnibus test, and the set-level inference — use one shared set of 9,999 label shuffles drawn at seed 19860426 and saved with the results <!-- N-04, N-05 -->. pi0 is reported with permutation-calibrated uncertainty (the plug-in estimator applied to each saved label shuffle, giving the distribution of the estimate under the contrast's own null).

【訳】各遺伝子は exact 置換 Brunner–Munzel 検定(両側)で被曝対散発の差を検定した — 全 C(n, nx) 割当の完全枚挙のため遺伝子レベルの p 値にシードは不要である。実装はコードと共に提供する自作の置換枚挙であり、ここでは exact モードで実行する(モンテカルロモードは後述の REO 評価でのみ用いる)— どの p 値にも漸近近似は入らない。効果 P(X<Y) は超過確率として直接読め、この群サイズでは置換分布が完全枚挙できるため、遺伝子レベルの推論はカウントモデルの仮定もモンテカルロ誤差も持ち込まない。遺伝子レベルの推論は exact p 値に対する Storey q 値(plug-in pi0、λ = 0.5)で、プロトコル共通の 0.10 水準を閾値とする <!-- N-05 -->。推定量はチューニングでなく作業仮説と設計から a priori に固定した: 仮説は多数の遺伝子に広がる弱いシグナル(群内の症例混合で減衰した小さな遺伝子別効果)なので、π0 を 1 に置くことはまさにそのシグナルの不在を補正に組み込むことになる。固定 λ plug-in の保守性は帰無 p 値の周辺一様性のみを要件とし、これは exact 検定が遺伝子間の依存に関わらず保証する。一方、λ の適応的選択は独立ないし弱依存の検定を仮定する保証しか持たず、単一のラベルベクトルを共有する共発現遺伝子には適用できない。さらに弱く広がったシグナルの下では対立 p 密度はほぼ平坦で、λ を上げてもバイアス削減は僅かなのに分散コストは一次である — λ = 0.5 は無チューニングの Storey (2002) 既定値である。対比レベルの帰無分布 — 後述の pi0 不確かさ・オムニバス検定・セットレベル推論 — は、シード 19860426 で抽出し結果とともに保存した 9,999 回のラベルシャッフルの単一セットを共有する <!-- N-04, N-05 -->。pi0 は置換較正した不確かさ(保存済みの各ラベルシャッフルに plug-in 推定量を適用して得る、当該対比自身の帰無の下での推定量分布)を付して報告する。

### Contrast-level omnibus inference

<!-- 写像元: 410(オムニバス・R(α)); 非2値化は C-16 -->

The contrast-level inferential claim is carried by an omnibus permutation test: the question of whether a contrast carries any label-aligned signal is answered here, not by the size of its gene list. Higher Criticism (α0 = 0.1) is the pre-specified primary omnibus statistic, chosen for its sensitivity to many weak effects (Donoho & Jin 2004), with count- and max-type rows reported descriptively <!-- N-05 -->. The full rejection curve R(α) is retained so that results do not depend on displaying a single threshold count.

<!-- C-16 -->
No contrast-level binary significance label is assigned: contrast-level evidence is reported continuously — the per-gene q-values, the pre-specified omnibus p and the rejection curve — following the methodological guidance against dichotomizing evidence near a threshold (Wasserstein & Lazar 2016; Greenland et al. 2016; Amrhein et al. 2019), with the interpretation of each outcome pattern pre-assigned before the results were seen rather than left to post-hoc labeling (interpretation map, Methods).

【訳】対比水準の推論的主張はオムニバス置換検定が担う: 「対比にラベル整列シグナルがあるか」への回答はここで与えられ、遺伝子リストの大きさによってではない。Higher Criticism(α0 = 0.1)が事前指定の主オムニバス統計量であり(多数の弱い効果への感度で選定 — Donoho & Jin 2004)、count 型・max 型の行は記述的に報告する <!-- N-05 -->。棄却曲線 R(α) 全体を保持し、結果が単一閾値の件数表示に依存しないようにする。対比に2値の有意性ラベルは付与しない: 対比水準の証拠は連続量のまま — 遺伝子別 q 値、事前指定オムニバスの p、棄却曲線 — で報告する。これは閾値近傍の証拠の2値化を退ける方法論的勧告(Wasserstein & Lazar 2016; Greenland et al. 2016; Amrhein et al. 2019)に従うものであり、各結果パターンの解釈は結果を見る前に事前割当されている(解釈マップ、Methods)。<!-- C-16 -->

### Gene-set level inference

<!-- 写像元: 420, 430(ORA 注釈), lib/gsea_permutation.R, lib/gsea_collections.R -->

Set-level enrichment consumes the whole per-gene ranking (threshold-free; no DEG-list cut decides what is tested), ranked by tie-averaged normal scores of the signed Brunner–Munzel statistic. The enrichment score is the weighted running sum (gseaParam = 1) evaluated at tie-block boundaries only, so that no arbitrary tie-break injects an order the data do not contain; on tie-free input it equals the standard GSEA statistic (Subramanian et al. 2005), an equality enforced by an automated test in the analysis code <!-- N-72 -->. The implementation is in-house. The null reuses the identical label shuffles saved by the gene-level test (9,999 per contrast <!-- N-07 -->), so gene- and set-level results rest on one permutation null. Inference is the per-set, sign-conditional permutation p-value with Benjamini–Hochberg adjustment within each collection, q_bh < 0.10
<!-- N-72 -->; no cross-family claim is made. Unlike the gene level, π0 is held at its conservative bound of 1 here: tens to thousands of dependent p-values riding one collective drift mode leave the plug-in estimate variance-dominated, with errors pointing anticonservative exactly in the realizations where the procedure is most fragile. Four MSigDB collections were tested (msigdbr 26.1.0): Hallmark (50), C2:CP (3,910; Reactome, WikiPathways, KEGG MEDICUS, BioCarta and PID), C5:GO:BP (7,538), and a radiation-curated C2:CGP family (28) whose curation rule was fixed before this inference touched real data <!-- N-29, N-55 -->. Sets
outside the size window 15–500 were excluded <!-- N-72 -->.

<!-- C-06 -->
The complete gene-set inference was calibrated on held-out null replicates — label permutations pushed through the identical procedure — before being applied to the real contrasts, and the choice of the set-level FDR procedure itself was fixed by this calibration prior to any real-data run (Supplementary Methods) <!-- N-06 -->. A spike-in control (one Hallmark set × 1.15 in the 9 samples of one group) served as the sensitivity counterpart
<!-- N-30 -->. (サプリ層 = Supplementary Methods / Results 節) <!-- C-14 -->
As a descriptive complement, the discovered R_Tumor list was annotated by one-sided hypergeometric over-representation against the identical set universe (up / down / combined lists; universe = the contrast's tested genes; BH within family × list), reported in full as hypothesis-generating annotation (Supplementary Table SY) <!-- N-59 -->. The two procedures' p-values refer to different randomness: only the label-permutation null refers to the experiment actually performed, and experiment-level claims rest on it; the over-representation q-values describe the discovered list against a gene-sampling reference (Goeman & Bühlmann 2007).

【訳】セットレベルの濃縮は遺伝子別ランキング全体を入力とする(閾値フリー — 検定対象を DEG リストの切断が決めることはない)。ランキングは符号付き Brunner–Munzel 統計量の tie 平均正規スコア。濃縮スコアは tie ブロック境界でのみ評価する重み付き累積和(gseaParam = 1)で、データに存在しない順序を恣意的な tie-break で注入しないためであり、tie のない入力では標準 GSEA 統計量(Subramanian et al. 2005)に厳密に一致する — この一致は解析コード内の自動テストが強制する <!-- N-72 -->。実装は自作。帰無は遺伝子レベル検定が保存した同一のラベルシャッフルを再利用し(対比あたり 9,999 <!-- N-07 -->)、遺伝子レベルとセットレベルの結果は単一の置換帰無の上に立つ。推論はセット別・符号条件付きの置換 p 値に各コレクション内の Benjamini–Hochberg 調整、q_bh < 0.10 <!-- N-72 --> であり、family 横断の主張はしない。遺伝子レベルと異なりここでは π0 を保守的上限の 1 に置く: 単一の集団ドリフトモードに乗る数十〜数千の依存 p 値の下では plug-in 推定は分散支配となり、その誤差は手続きが最も脆い実現でちょうど反保守側を向く。MSigDB の4コレクションを検定した(msigdbr 26.1.0): Hallmark(50)、C2:CP(3,910; Reactome・WikiPathways・KEGG MEDICUS・BioCarta・PID の統合)、C5:GO:BP(7,538)、放射線キュレーションの C2:CGP ファミリー(28 — キュレーション規則は本推論が実データに触れる前に固定)<!-- N-29, N-55 -->。サイズ窓 15–500 の外のセットは除外した
<!-- N-72 -->。セットレベル推論一式は、実対比への適用前に held-out 帰無レプリケート(同一手続きに通したラベル置換)で較正し、セットレベル FDR 手続きの選択自体もこの較正により実データ実行前に固定した(Supplementary Methods)<!-- N-06 -->。spike-in 対照(1つの
Hallmark セット × 1.15 倍を一方の群の 9 検体に)が感度側のカウンターパートである
<!-- N-30 -->。<!-- C-06 -->
記述的な補完として、発見済み R_Tumor リストを同一のセット宇宙に対する片側超幾何の過剰代表で注釈した(up / down / 合算の3リスト; universe = 当該対比の検定遺伝子; family × list 内 BH)— 全結果を仮説生成の注釈として報告する(Supplementary Table SY)<!-- N-59 -->。2つの手続きの p 値は異なるランダム性を参照する: 実際に行った実験のランダム性を参照するのはラベル置換帰無のみであり、実験レベルの主張はそこに乗る。過剰代表の q 値は、gene-sampling の参照に対する発見済みリストの記述である(Goeman & Bühlmann 2007)。<!-- C-14 -->

### Between-stratum concordance of the exposure contrast

<!-- 写像元: diagnostics/signature_agreement.R(節題は判断点2決着に従い記述的言い回し。「signature agreement」は内部語として本文不使用) -->

For each tissue, the exposed-vs-sporadic contrast was summarized per gene by the signed Brunner–Munzel statistic separately in the RET and BRAF strata, and the two genome-wide profiles compared by Spearman correlation over shared genes — threshold-free and without any dose assumption. Each coefficient travels with a permutation-calibrated reference interval: labels were shuffled independently within each contrast (9,999 shuffles, per-contrast seeds from base 19450809 <!-- N-73 -->) and the shuffled profiles correlated, giving the null spread of rho when neither stratum carries label-aligned structure. The normal-tissue pair is the pre-specified hypothesis-bearing comparison; the tumor pair is a descriptive completion. An observed rho outside the interval indicates shared label-aligned structure; by itself it does not identify that structure as an exposure trace rather than a shared covariate.

【訳】各組織で、被曝対散発対比を RET 層・BRAF 層で別々に遺伝子別の符号付き Brunner–Munzel 統計量に要約し、2つのゲノムワイドプロファイルを共有遺伝子上の Spearman 相関で比較した — 閾値フリーで線量仮定なし。各係数には置換較正した参照区間が付く: 対比内で独立にラベルをシャッフルし(9,999 回、対比別シードは基底 19450809 から <!-- N-73 -->)、シャッフル済みプロファイル同士を相関させ、どちらの層もラベル整列構造を持たないときの rho の帰無分布を得る。正常組織ペアが事前指定の仮説担持比較であり、腫瘍ペアは記述的補完である。区間の外の rho は共有されたラベル整列構造を示すが、それ自体は、その構造が共有共変量ではなく曝露痕跡であることを識別しない。

### REO panel: construction and out-of-sample evaluation

<!-- 写像元: 510, 520, 530 -->

The REO panel provides a per-case, normalization-free readout of the training-band contrast (R_Sporadic versus R_High), built so that it can be examined out of sample in the intermediate bands; its construction and evaluation are as follows. Relative expression ordering (REO; Geman et al. 2004; Wang et al. 2015) works on within-sample gene ranks and is therefore free of between-sample normalization; expression was taken as TPM — recomputed from the selected count assay and the exon-union gene lengths (data acquisition, above), not the discarded GDC columns — because within-sample comparison across genes requires length normalization. Candidate pairs were generated in the RET tumor training groups from the top 500 genes by Brunner–Munzel effect magnitude |effect − 0.5| (a threshold-free pool; the q < 0.10 DEG set is deliberately not used). A pair qualifies when its within-sample log2-TPM difference r has a stable sign in Sporadic (dead zone |r| < log2(1.2); at most one exception among non-dead-zone samples; 10th percentile of |r| ≥ log2(1.5)) and reverses in more than 50% but not all of the High samples <!-- N-74 -->, ranked by the shift in median r. The final panel was selected greedily in rank order, excluding gene reuse and pairs whose per-sample r profiles correlate at Spearman ≥ 0.75 with a kept pair, to a target of 10 pairs
<!-- N-75, N-37 -->. A sample's reversal score is the number of panel pairs whose within-sample difference r clears the dead zone and runs opposite to the training-Sporadic reference sign (an integer 0–10); the classification boundary is training-derived (positive above the maximum reversal score among the training Sporadic samples), and the realized panel, boundary and construction yield are reported in Results and Table(仮) <!-- N-39 -->.

<!-- C-08 の手法側 -->
The finalized panel was then applied, untouched, to the intermediate-band RET tumors it was never trained on (R_Low, R_Mid <!-- N-10 -->) as a graded, out-of-sample check. Because both training bands entered panel construction, any separation within them is circular; Low and Mid are the only bands where the score can be examined out of sample, and the design fixes the expected direction in advance (the higher assigned-share band above the lower). The single pre-specified test is therefore the one-sided Brunner–Munzel comparison of Mid over Low (Monte Carlo, seed 19860426
<!-- N-75 -->). This evaluation does not alter the panel or its boundary, and the band-wise score profile is reported descriptively without assuming any dose–response form.

【訳】REO パネルは、訓練帯対比(R_Sporadic 対 R_High)の症例別・正規化フリーの読み出しを与える計器であり、中間帯で out-of-sample 検査ができるように構築する。構築と評価は以下の通りである。相対発現順序(REO; Geman et al. 2004; Wang et al. 2015)は検体内の遺伝子順位で動作し、検体間正規化から自由である。発現は TPM とした — 採用カウントアッセイと exon-union 遺伝子長から再計算したもので(前述のデータ取得節)、読み込み時に破棄した GDC 列ではない — 検体内での遺伝子間比較が長さ正規化を要求するためである。候補ペアは RET 腫瘍訓練群で、Brunner–Munzel 効果の大きさ |effect − 0.5|上位 500 遺伝子から生成した(閾値フリーのプール; q < 0.10 の DEG 集合は意図的に使わない)。ペアの適格条件: 検体内 log2-TPM 差 r が Sporadic で安定した符号を持ち(dead zone |r| < log2(1.2); 非 dead-zone 検体で例外は高々1; |r| の第10百分位 ≥ log2(1.5))、High 検体の 50% 超かつ全例未満で逆転すること <!-- N-74 -->。順位は
中央値 r のシフトで付けた。最終パネルは順位順の貪欲選定で、遺伝子の再使用と、採用済みペアと Spearman ≥ 0.75 で相関するペアを除外し、目標 10 ペアとした
<!-- N-75, N-37 -->。検体の逆転スコアは、パネルのうち検体内差 r が dead zone を越えかつ訓練 Sporadic の参照符号と逆向きのペア数である(0–10 の整数)。分類境界は訓練由来である(訓練 Sporadic 検体の最大逆転スコアを越えたら positive)。実現したパネル・境界・構築の歩留まりは Results と表(仮)に報告する <!-- N-39 -->。確定したパネルは、訓練に使っていない中間帯 RET 腫瘍(R_Low・R_Mid <!-- N-10 -->)へ変更なしで適用し、段階的な out-of-sample 検査とした。訓練2帯はパネル構築に入っているため帯内の分離は循環であり、スコアを out-of-sample で検査できる帯は Low と Mid のみで、設計が方向を事前に指定する(assigned share の高い帯が上)。したがって唯一の事前指定検定は Mid over Low の片側 Brunner–Munzel 比較である(モンテカルロ、シード 19860426 <!-- N-75 -->)。この評価はパネルにも境界にも変更を加えず、帯別スコアプロファイルは線量反応の形を仮定せず記述的に報告する。<!-- C-08 の手法側 -->

### REO evaluation diagnostics

<!-- 写像元: diagnostics/reo_lowmid_{outliers,purity,confound}.R -->

<!-- C-09 -->
Because the evaluation cohort is deliberately unfiltered, its QC is reported as diagnostics with no exclusion authority: (i) the training outlier screen (PC-OD) was mirrored on the R_Low/R_Mid tumors, and the evaluation was run once on all eligible cases — the screen reports counts only <!-- N-43 -->; (ii) tumor purity for the evaluation bands was estimated on the same common scale as training by pooling the whole RET cohort in one contamDE run
<!-- N-44 -->; <!-- C-12 --> (iii) because the reversal score correlates with purity within the evaluation bands, a graded band profile could in principle be purity-driven rather than band-driven — the two are separated by the partial Spearman correlation of band with score conditioning on purity (permutation reference) and by purity-stratified one-sided Brunner–Munzel comparisons of Mid over Low (two strata split at the median tumor purity) <!-- N-46, N-47 -->.

【訳】評価コホートは意図的に未濾過のため、その QC は除外権限を持たない診断として報告する <!-- C-09 -->: (i) 訓練側の外れ値スクリーン(PC-OD)を R_Low/R_Mid 腫瘍に鏡映し、評価は適格全例で一度だけ実行 — スクリーンは件数の報告のみ <!-- N-43 -->; (ii) 評価帯の腫瘍純度は、RET コホート全体を1回の contamDE にプールして訓練と同一の共通尺度で推定 <!-- N-44 -->; (iii) 逆転スコアは評価帯内で純度と相関するため、帯の段階的プロファイルは原理的には帯でなく純度駆動であり得る — 両者は、純度を条件付けた帯–スコアの偏 Spearman 相関(置換参照)と、純度層別(腫瘍純度の中央値で2層)の Mid over Low 片側 Brunner–Munzel 比較で分離する <!-- C-12, N-46, N-47 -->。

### External anchor cross-reference

<!-- 写像元: diagnostics/external_gene_anchors.R -->

<!-- C-13 -->
Externally validated radiation-associated gene lists (qRT-PCR-validated cores of Abend 2013, Abend 2012, Dom 2012, and CLIP2 <!-- N-53 -->) were cross-referenced against each contrast's q < 0.10 gene set as a descriptive membership count (k of n list genes among the contrast's tested genes), with no enrichment statistic: the lists are small, platforms and contrasts differ across sources, and the reading was fixed symmetrically in advance — any count is reported as description and no claim moves with the outcome.

【訳】外部で検証済みの放射線関連遺伝子リスト(Abend 2013・Abend 2012・Dom 2012 の qRT-PCR 検証コアと CLIP2 <!-- N-53 -->)を、各対比の q < 0.10 遺伝子集合と員数照合した(対比の検定遺伝子中のリスト遺伝子 n のうち k)— 濃縮統計量は計算しない:リストは小さく、プラットフォームと対比は出典間で異なり、読みは事前に対称に固定した — どの員数も記述として報告し、どの主張も結果で動かない。<!-- C-13 -->

### Software, seeds and reproducibility

<!-- 写像元: config.R, tests/, run 保全記録 -->

The publication run executed R 4.5.3 on Ubuntu 24.04 <!-- N-02 --> with the reference BLAS/LAPACK 3.12.0 <!-- N-03 -->, four workers and 9,999 label shuffles as the fixed reproduction contract <!-- N-04 -->, from a clean repository state <!-- N-01 -->. The environment is rebuilt from a fully date-pinned container build: an immutable Ubuntu base image (noble-20260410) with its matching apt snapshot, R 4.5.3 built from source against the reference BLAS, and all R packages installed from a dated repository snapshot (CRAN and Bioconductor 3.22; snapshot 2026-04-09) <!-- N-79 -->. The canonical inference seed is 19860426; diagnostics draw documented seeds from base 19450809 <!-- N-05, N-06 -->. The full pipeline was executed independently on two x86-64 machines: the 1,819 raw input files agree by md5 across the two machines <!-- N-52 --> and the pipeline's test suite passes on both (415 tests, 0 failures <!-- N-51 -->); primary artifacts were verified identical across machines. The reproduction contract, including this bit-level identity, is defined on x86-64 Linux — the platform the container build targets (linux/amd64) — and is not asserted for other architectures. The core statistical machinery — the exact/Monte-Carlo Brunner–Munzel test, the Storey plug-in estimator, the tie-block enrichment statistic with its permutation q-values, DEGES-MUREN normalization, the contamDE purity model, the principal-component outlier screen and the REO panel machinery — is implemented in-house in the analysis code rather than taken from packages, with method sources cited at first mention; external packages actually used (edgeR, SummarizedExperiment, msigdbr, limma, GenomicDataCommons, rtracklayer, Rcpp, MASS and others) are pinned in the container build and listed with versions in Supplementary Table(仮). Analysis code, versioned inputs and the container build recipe sufficient to regenerate the reported analyses accompany the paper. Pipeline design, the prototype scripts and all analytical decisions are the authors'. Generative AI coding assistants (Claude Fable 5, Anthropic; GPT-5.6sol, OpenAI) were used for script debugging and refactoring; only open-access data were handled by these tools, and all resulting code passes the verification described above (test suite, cross-machine reproduction, automated equivalence checks), for which the authors take full responsibility.

【訳】本番実行は R 4.5.3 / Ubuntu 24.04 <!-- N-02 -->、参照 BLAS/LAPACK 3.12.0
<!-- N-03 -->、ワーカー 4・ラベルシャッフル 9,999 を固定された再現契約とし<!-- N-04 -->、クリーンなリポジトリ状態から実行した <!-- N-01 -->。環境は全層を日付固定したコンテナビルドから再構築される: 不変タグの Ubuntu 基底イメージ(noble-20260410)と同日付の apt スナップショット、参照 BLAS に対してソースからビルドした R 4.5.3、日付スナップショット(CRAN + Bioconductor 3.22; 2026-04-09)から導入した全 R パッケージ <!-- N-79 -->。正準の推論
シードは 19860426、診断は基底 19450809 から文書化されたシードを引く <!-- N-05, N-06 -->。全パイプラインは2台の x86-64 マシンで独立に実行され、1,819 の生入力ファイルは二機間で md5 一致 <!-- N-52 -->、テストスイートは両機で通過(415 テスト、失敗 0
<!-- N-51 -->)、主要成果物は機械間で同一と検証された。この再現契約は、ビットレベル同一性を含めて x86-64 Linux(コンテナビルドが対象とするプラットフォーム = linux/amd64)上で定義されており、他のアーキテクチャについては主張しない。核心の統計機構 — exact/モンテカルロ Brunner–Munzel 検定、Storey plug-in 推定量、tie ブロック濃縮統計量とその置換 q 値、DEGES-MUREN 正規化、contamDE 純度モデル、主成分外れ値スクリーン、REO パネル機構 — はパッケージからでなく解析コード内の自作実装であり、手法の原典は初出箇所で引用する。実際に使用した外部パッケージ(edgeR・SummarizedExperiment・msigdbr・limma・GenomicDataCommons・rtracklayer・Rcpp・MASS ほか)はコンテナビルドで版固定され、Supplementary Table(仮)に版つきで一覧する。報告した解析の再生成に足る解析コード・版管理された入力・コンテナビルドレシピを論文に添付する。パイプラインの設計・原型スクリプト・全ての解析上の判断は著者らによる。スクリプトのデバッグとリファクタリングには生成 AI コーディング支援(Claude Fable 5, Anthropic; GPT-5.6sol, OpenAI)を用いた。これらのツールが扱ったのは open-access データのみであり、結果としての全コードは上述の検証(テストスイート・機械間再現・自動等価性検査)を通過している。責任は全て著者らが負う。

## Results

### 1. Cohort(C-10, C-15)

<!-- C-10 -->
Of the 440 cases of the REBC-THYR cohort (Morton et al. 2021) <!-- N-11 -->, the main analysis cohort comprised 63 paired, driver-stratified cases — 9 B_High, 27 B_Sporadic, 15 R_High and 12 R_Sporadic <!-- N-09 --> — reached through the pre-specified flow of driver classification, band eligibility, pairing, outlier screening and purity thresholding (Tables 1–2(仮)). Most of the reduction reflects the pre-specified eligibility restrictions — driver classification and the restriction to the extreme bands — rather than technical losses to pairing, outlier screening or purity; the full flow, stratified by driver, is shown in the cohort-flow figure(仮) <!-- N-08 -->. The REO evaluation set added 36 paired RET tumors of the intermediate bands (17 R_Low, 19 R_Mid
<!-- N-10 -->). In both driver strata the High group sat somewhat older at surgery than the Sporadic group <!-- N-12, N-13 -->.
<!-- C-15 -->
The three disclosed candidate confounders (Methods) are reported in Table 2(仮): the per-group age distributions, with interval estimates of the between-group age difference as table footnotes (disclosed, not tested <!-- N-64, N-65 -->), and the per-group sex and fusion-partner compositions. Age at exposure exists only for exposed cases and admits no between-group comparison <!-- N-63 -->.

【訳】REBC-THYR コホート(Morton et al. 2021)の 440 症例 <!-- N-11 --> のうち、主解析コホートは、driver 分類・帯適格性・ペア有無・外れ値スクリーン・純度閾値という事前指定のフロー(表1–2(仮))を経た 63 のペア付き driver 層別症例 — B_High 9・B_Sporadic 27・R_High 15・R_Sporadic 12
<!-- N-09 --> — である。削減の大半は事前指定の適格性制限 — driver 分類と両端帯への制限 — によるもので、ペア・外れ値スクリーン・純度による技術的損失ではない。フローの全体は driver 層別でコホートフロー図(仮)に示す <!-- N-08 -->。REO 評価セットとして中間帯のペア付き RET 腫瘍 36 例(R_Low 17・R_Mid 19 <!-- N-10 -->)を加えた。いずれの driver 層でも High 群は Sporadic 群よりやや高い手術時年齢に分布した <!-- N-12, N-13 -->。開示対象の候補交絡因子3つ(Methods)は表2(仮)に報告する: 群別の年齢分布(群間年齢差の区間推定は表脚注 — 開示であって検定ではない <!-- N-64, N-65 -->)、および群別の性構成と融合パートナー構成。被曝時年齢は被曝症例にのみ定義され、群間比較は成立しない <!-- N-63 -->。

### 2. Gene-level differential expression(C-01〜C-04, C-16, 410)

<!-- C-01 -->
In R_Tumor — the contrast where the pre-specified prediction map expects the signal — 1,765 of 15,621 tested genes differed between the High and Sporadic groups at Storey q < 0.10 <!-- N-16, N-15 -->, and the pre-specified contrast-level omnibus supported the presence of signal (Higher Criticism p = 0.0112
<!-- N-20 -->) (Figs. 1–2(仮), Table 3(仮)). <!-- C-02 --> The discovered genes ran in both directions: 971 higher and 794 lower in the exposed group <!-- N-17 -->. <!-- C-03 --> In R_Normal no gene reached q < 0.10 and the
omnibus lent no support (HC p = 0.3199) <!-- N-16, N-20 -->. <!-- C-04 --> B_Tumor likewise yielded no discovery (0 genes at q < 0.10; HC p = 0.1815)
<!-- N-16, N-20 -->; under the pre-fixed reading this cell is direction-agnostic, and its quiet is not read as a specificity control. <!-- C-16 -->
In B_Normal the evidence is reported as it stands, without a binary label: one gene crossed the gene-level threshold (BHLHB9 — an X-linked gene <!-- N-80 --> — effect 0.967, q = 0.013
<!-- N-22 -->); the pre-specified primary omnibus gave HC p = 0.0773 <!-- N-20 --> while the descriptive max-statistic row reached p = 0.0125
<!-- N-23 -->, and the contrast's π0 estimate (0.727) sat below those of the other two quiet contrasts (0.955, 0.943) <!-- N-19 -->. The reading pre-assigned
to this pattern is taken up in Discussion.

【訳】事前指定の予測マップがシグナルを期待する対比である R_Tumor では、検定対象 15,621 遺伝子のうち 1,765 が Storey q < 0.10 で High 群と Sporadic 群の間で発現差を示し <!-- N-16, N-15 -->、事前指定の対比レベル・オムニバスがシグナルの存在を支持した(Higher Criticism p = 0.0112 <!-- N-20 -->)(図1–2(仮)・表3(仮))。発見遺伝子は双方向に分布した: 被曝群で高発現 971・低発現 794 <!-- N-17 -->。R_Normal では q < 0.10 の遺伝子はなく、オムニバスの支持もなかった(HC p = 0.3199)<!-- N-16, N-20 -->。B_Tumor も同様に発見なし(q < 0.10 は 0 遺伝子; HC p = 0.1815)<!-- N-16, N-20 --> — 事前固定の読みに従い、このセルは方向不可知であり、その静けさを特異性の対照とは読まない。B_Normal の証拠は2値ラベルなしにそのまま報告する: 1遺伝子が遺伝子レベル閾値を越え(BHLHB9 — X 連鎖遺伝子 <!-- N-80 --> — effect 0.967、q = 0.013 <!-- N-22 -->)、事前指定の主オムニバスは HC p = 0.0773
<!-- N-20 -->、記述的な max 統計量行は p = 0.0125 <!-- N-23 -->、この対比の π0 推定値(0.727)は他の2つの静かな対比(0.955・0.943)より低かった <!-- N-19 -->。このパターンに
事前割当された読みは Discussion で扱う。

### 3. Gene-set level(C-05, C-06, 420 + D6)

<!-- C-05 -->
At the gene-set level, the calibrated test declared no set at q_bh < 0.10 in any of the 16 contrast × collection cells <!-- N-27 --> (four MSigDB collections, Methods; 6,141–6,242 sets per contrast after filtering <!-- N-28 -->). The smallest adjusted value anywhere was q = 0.114, in the B_Tumor × radiation cell <!-- N-28 --> (complete listing in Supp. Data 1(仮); Table 4(仮)).
<!-- C-06 -->
Under signal-free inputs — label permutations pushed through the identical procedure (Methods) — the set-level machinery produced at least one discovery in 102 of 1,600 replicates pooled across the 16 contrast × collection cells (0.064; nominal 0.10 <!-- N-56 -->), with a single disclosed excess (B_Normal/Hallmark: 0.18, 95% CI 0.110–0.270 <!-- N-25 -->). Under the synthetic positive control — one Hallmark set inflated 1.15-fold in one group (Methods) — the planted set was recovered at rank 1 of 50 (q = 0.0101 <!-- N-31 -->), with no other set at q < 0.10 <!-- N-32 --> (per-cell calibration in Supp. Tab. 1(仮)).

【訳】遺伝子セットレベルでは、較正済み検定は 16 の対比 × collection セルのいずれでも q_bh < 0.10 のセットを宣言しなかった <!-- N-27 -->(MSigDB 4 コレクション、Methods 参照; フィルタ後は対比あたり 6,141–6,242 セット <!-- N-28 -->)。全セルで最小の調整値は B_Tumor × radiation セルの q = 0.114 だった <!-- N-28 -->(全結果は Supp. Data 1(仮)・表4(仮))。シグナルなしの入力 — 同一手続きに通したラベル置換(Methods)— の下では、セットレベル機構は 16 の対比 × collection セル合算 1,600 レプリケート中 102 で1つ以上の発見を生じ(0.064; 名目 0.10 <!-- N-56 -->)、開示済みの超過は1セルのみだった(B_Normal/Hallmark: 0.18、95% CI 0.110–0.270 <!-- N-25 -->)。合成陽性対照 — 1つの Hallmark セットを一方の群で 1.15 倍(Methods)— の下では、埋め込んだセットが 50 中 rank 1(q = 0.0101
<!-- N-31 -->)で回収され、それ以外に q < 0.10 のセットはなかった <!-- N-32 -->
(セル別較正は Supp. Tab. 1(仮))。

### 4. Composition of the discovered list(C-14, 430)

<!-- C-14 -->
In the descriptive annotation of the discovered list (hypothesis-generating; Methods), the 794 genes lower in the exposed group were strongly concentrated in proliferation, cell-cycle and DNA-repair programs (E2F targets 46/199, G2M checkpoint 41/198, Reactome DNA repair 68/322 <!-- N-61, N-62 -->), a single theme that extends into the radiation-curated family, whose leading flagged sets are themselves cell-cycle genes responding to irradiation (46 of 126 in the down list, expected 6.4 <!-- N-60 -->). The 971 genes higher in the exposed group showed no such concentration in any curated family <!-- N-59 --> (full results in Supp. Tab. 2(仮)).

【訳】発見済みリストの記述的注釈(仮説生成; Methods)では、被曝群で低発現の 794 遺伝子は増殖・細胞周期・DNA 修復プログラムに強く集中し(E2F targets 46/199・G2M checkpoint 41/198・Reactome DNA repair 68/322 <!-- N-61, N-62 -->)、この単一テーマは放射線キュレーション・ファミリーにも及ぶ — そこでフラグが立った上位セット自体が照射に応答する細胞周期遺伝子である(down リストで 126 中 46、期待 6.4 <!-- N-60 -->)。被曝群で高発現の 971 遺伝子には、どのファミリーにもそのような集中はなかった
<!-- N-59 -->(全結果は Supp. Tab. 2(仮))。

### 5. Between-stratum concordance of the exposure contrast(C-07, C-17)

<!-- C-07 -->
The pre-specified between-stratum comparison in normal tissue — Spearman correlation, across strata, of the per-gene signed statistics of the exposure contrast — gave rho = +0.376 over 15,459 shared genes, inside its label-shuffle reference interval ([−0.46, +0.46]; two-sided p = 0.1199)
<!-- N-34 -->. With no within-contrast signal in R_Normal <!-- N-16 -->, the pre-fixed reading applies: whether the two normal-tissue contrasts share an exposure trace is not identifiable here. <!-- C-17 --> The symmetric
tumor-pair comparison, computed as a design completion, is reported in Supp. Tab. 3(仮) <!-- N-33 --> and taken up as hypothesis-generating in Discussion.

【訳】正常組織における事前指定の層間比較 — 曝露対比の遺伝子別符号付き統計量を層間で Spearman 相関 — は、共有 15,459 遺伝子で rho = +0.376 となり、ラベルシャッフルの参照区間([−0.46, +0.46]; 両側 p = 0.1199)の内側だった <!-- N-34 -->。R_Normal に対比内シグナルがない <!-- N-16 --> ため、事前固定の読みが適用される: 2つの正常組織対比が曝露痕跡を共有するか否かは、ここでは識別できない。対称補完として計算した腫瘍ペアの比較は Supp. Tab. 3(仮)に報告し <!-- N-33 -->、Discussion で仮説生成として扱う。

### 6. REO grading(C-08, C-09, C-12, 510–530)

<!-- C-08 -->
Panel construction evaluated 57,694 candidate pairs from the training pool's 317 up- and 182 down-genes; 153 passed all criteria, with median r shifts of 1.159–4.700 and reversal rates 0.53–0.87 <!-- N-35, N-36 -->. The 10-pair REO panel separated its training groups as designed (all 12 R_Sporadic negative, 13 of 15 R_High positive; boundary at score > 2
<!-- N-38, N-37 -->; panel composition in Table(仮) <!-- N-39 -->). Applied unchanged to the intermediate bands, the per-case reversal score rose in band order — medians 0 / 1 / 4 / 6 for Sporadic / Low / Mid / High <!-- N-41 --> (Fig. 3(仮)) — with R_Low classified 9 negative / 8 positive
and R_Mid 8 / 11 <!-- N-42 -->. The one comparison available out of sample (Mid vs Low; direction pre-specified, Methods) gave one-sided Brunner–Munzel p = 0.1127, effect P(Low<Mid) = 0.616 <!-- N-40 -->. The graded profile is reported as a descriptive observation; no dose–response form is assumed or claimed.
<!-- C-09 -->
The outlier screen mirrored from training flagged no case in either evaluation band (0 of 17 and 0 of 19 <!-- N-43 -->); the screen carries no exclusion authority, and the evaluation was run once on all eligible cases.
<!-- C-12 -->
Within the evaluation bands the reversal score correlated with tumor purity (pooled Spearman +0.538 <!-- N-45 -->), while the band–score correlation was +0.142 and, conditioned on purity, +0.146 (partial Spearman; one-sided permutation p = 0.2162) <!-- N-48, N-46 -->; band and purity themselves were nearly uncorrelated (+0.036 <!-- N-48 -->). Purity-stratified comparisons are given in Supplementary <!-- N-47 -->.

【訳】パネル構築では、訓練プールの up 317・down 182 遺伝子から 57,694 候補ペアを評価し、153 が全基準を通過した(中央値 r シフト 1.159–4.700、逆転率 0.53–0.87)<!-- N-35, N-36 -->。10 ペアの REO パネルは訓練群を設計どおり分離した(R_Sporadic は 12 例全て negative、R_High は 15 例中 13 が positive; 境界は score > 2 <!-- N-38, N-37 -->;パネル構成は表(仮)<!-- N-39 -->)。パネルを変更せず中間帯に適用すると、症例別の逆転スコアは帯順に上昇した — 中央値は Sporadic / Low / Mid / High で 0 / 1 / 4 / 6
<!-- N-41 -->(図3(仮))— 分類は R_Low が negative 9 / positive 8、R_Mid が 8 / 11 <!-- N-42 -->。out-of-sample で検査可能な唯一の比較(Mid 対 Low; 方向は事前指定、
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

Conditioned on driver and tested within each contrast, the exposure comparison produced gene-level signal in exactly one place: RET tumors, where 1,765 genes crossed q < 0.10 and the contrast-level omnibus supported the presence of signal <!-- N-16, N-20 -->. We read this as an exposure-associated expression difference within RET tumors <!-- C-01 -->. Two design properties discipline that reading. Under the exact label-permutation null, structure uncorrelated with the exposure label — including any general noisiness or instability of the RET-fusion transcriptome — is controlled at the declared error rates by construction; what could mimic the finding is covariate structure correlated with the label, and the disclosed candidates — age, sex and fusion-partner composition — are bounded below. And the association names no mechanism: the wording stays neutral as to whether a shared radiation trace feeds RET oncogenesis or the type of initiating lesion associates with driver selection (Limitations). Throughout this section, readings drawn from the interpretation map (Methods, Table 5(仮)) were assigned to outcome patterns before the reported results existed; interpretations that go beyond the map are labelled as hypothesis-generating. The remaining cells were quiet at the gene level — no discovery in R_Normal or B_Tumor, and B_Tumor's quiet is not read as a specificity control <!-- C-04 --> — while B_Normal showed a two-level divergence taken up below. At the set level, the calibrated test declared nothing in any of the 16 cells <!-- N-27 -->: whatever organization the discovered list displays along curated set boundaries (below), it does not reach the declared error rate under the label-permutation null <!-- C-05 -->.

【訳】driver で条件付け、各対比内で検定した結果、曝露比較が遺伝子レベルのシグナルを生んだのはただ一箇所 — RET 腫瘍であり、1,765 遺伝子が q < 0.10 を越え、対比レベルのオムニバスがシグナルの存在を支持した <!-- N-16, N-20 -->。我々はこれを RET 腫瘍内の曝露関連発現差として読む <!-- C-01 -->。この読みは設計の2性質が規律する。exact ラベル置換帰無の下では、曝露ラベルと無相関な構造 — RET 融合腫瘍のトランスクリプトームの一般的なノイズや不安定性を含む — は構成的に宣言誤り率で制御される。所見を模倣し得るのはラベルと相関する共変量構造であり、開示済みの候補 — 年齢・性・融合パートナー構成 — は後述で画定される。そして関連は機構を名指ししない: 共有放射線痕跡が RET 発がん機序に接続するのか、起始病変のタイプが driver 選択と連関するのかについて、文言は中立を保つ(Limitations)。本節を通じて、解釈マップ(Methods、表5(仮))に由来する読みは報告結果が存在する前に結果パターンへ割り当てられたものであり、マップを越える解釈には仮説生成のラベルを付す。残りのセルは遺伝子レベルで静かであり — R_Normal と B_Tumor に発見はなく、B_Tumor の静けさは特異性の対照とは読まない <!-- C-04 --> — B_Normal は後述する二水準の乖離を示した。セットレベルでは較正済み検定は 16 セルのいずれでも何も宣言しなかった <!-- N-27 -->: 発見済みリストがキュレート済みセット境界に沿って示す組織化(後述)がどうであれ、ラベル置換帰無の下で宣言済み誤り率には達しない <!-- C-05 -->。

### B_Normal: a two-level divergence read as confounding first(C-16)

<!-- C-16: 規則5の条件適用。D6 超過 0.18(N-25)の受けもここ。「q<0.10 陽性論」への受けは Q-16 -->

B_Normal is the one contrast where the two evidence levels diverge: a single discovered gene (BHLHB9 <!-- N-22 -->) against a primary omnibus lending no support (HC p = 0.0773 <!-- N-20 -->). No binary label is assigned (Methods); instead the pattern receives the reading the map assigned to it in advance — confounding first <!-- C-16 -->. A shared glandular exposure trace should appear at least as clearly on the RET side, whose exposed group sits higher on the exposure metric (median assigned share 86.7 versus 75.9 in the exposed BRAF group <!-- N-81 -->); a normal-tissue signal confined to the BRAF stratum therefore points toward that stratum's covariate asymmetries, and both named candidates point the same way. Age: the between-group difference is largest here (+8.0 years [3.0–12.0] <!-- N-65 -->), normal thyroid is the tissue with an established age effect (Coclet et al. 1989), and the same stratum's tumors, carrying the same age difference, produced nothing (0 genes at q < 0.10; all set-level cells quiet <!-- N-16, N-27 -->) — in this stratum, age moves normal tissue, not tumor. Sex: the same stratum carries the largest sex-composition asymmetry (23 of 27 versus 4 of 9 female <!-- N-13 -->), and the one discovered gene is X-linked <!-- N-80 -->. Neither candidate is tested; both are disclosed. One calibration caveat travels with this cell: the only held-out excess of the set-level calibration sits in B_Normal/Hallmark (0.18, 95% CI 0.110–0.270 <!-- N-25 -->), and any reading of that cell carries it.

【訳】B_Normal は二つの証拠水準が乖離する唯一の対比である: 発見遺伝子1件(BHLHB9 <!-- N-22 -->)に対し、主オムニバスは支持を与えない(HC p = 0.0773 <!-- N-20 -->)。2値ラベルは付与しない(Methods)。代わりにこのパターンは、マップが事前に割り当てた読み — 交絡第一 — を受ける <!-- C-16 -->。共有された腺の曝露痕跡ならば、被曝群が曝露指標上でより高い RET 側(assigned share 中央値 86.7 対 BRAF 被曝群 75.9 <!-- N-81 -->)に少なくとも同等に現れるはずであり、BRAF 層に限局した正常組織シグナルはこの層の共変量非対称を指す — そして名指しされた2候補は同じ方向を向く。年齢: 群間差はここで最大(+8.0 年 [3.0–12.0] <!-- N-65 -->)、正常甲状腺は年齢効果が確立している組織であり(Coclet et al. 1989)、同じ層の腫瘍は同じ年齢差を抱えながら何も生まなかった(q < 0.10 は 0 遺伝子、セットレベルも全セル静か <!-- N-16, N-27 -->)— この層では年齢は正常組織を動かし、腫瘍を動かさない。性: 同じ層は性構成の非対称も最大であり(女性 27 例中 23 対 9 例中 4 <!-- N-13 -->)、唯一の発見遺伝子は X 連鎖である <!-- N-80 -->。どちらの候補も検定せず、どちらも開示する。このセルには較正上の注意が1つ付随する: セットレベル較正の唯一の held-out 超過が B_Normal/Hallmark にあり(0.18、95% CI 0.110–0.270 <!-- N-25 -->)、当該セルのどの読みもこれを携行する。

### Candidate confounders(C-15, Q-13)

<!-- Q-13 の4観察(Vriens はここに一本化 — C-14 節の重複を解消)+ C-15 の開示推定。Methods の Candidate confounders と対構造 -->

The design discloses three candidate confounders rather than adjusting for them (Methods); the disclosed quantities bound what each can explain. Age at surgery differs between groups — Hodges–Lehmann +2.5 years [−1.0 to 6.0] in the RET stratum and +8.0 years [3.0 to 12.0] in the BRAF stratum <!-- N-64, N-65 --> — and four observations bound its reach. First, the established thyroid age effect concerns the follicular turnover of normal tissue (Coclet et al. 1989), whereas the discovered differences sit in tumors, where a published contrast spanning roughly 25 years of age — an order of magnitude wider than the RET-stratum difference here — produced no genome-wide differential expression and no set-level findings (Vriens et al. 2011). Second, within the exposed RET bands the age ordering runs against the score ordering: median age at surgery falls 30 / 25 / 23 years across Low / Mid / High <!-- N-12 --> while the reversal score rises 1 / 4 / 6 <!-- N-41 -->. Third, purity runs the wrong way for an artifact: the High band is the purest (0.814 versus 0.690 in Sporadic <!-- N-44 -->), yet proliferation programs read as attenuated, not strengthened, in the exposed group <!-- N-61 -->. Fourth, the internal control: the stratum with the largest age difference produced no tumor-side signal at all (B_Tumor <!-- N-16, N-27 -->), so attributing the R_Tumor differences — with an age gap roughly a third that size — to age contradicts the protocol's own BRAF stratum. Age, with latency inseparable from it, therefore remains named where it is live — the B_Normal reading above and the annotation below — and is not entered as a covariate because age at exposure is a component of the assigned-share definition (Methods). The other two candidates raise no comparable bound: sex is named for B_Normal above, and the fusion-partner composition shows no gradient across the bands (Table 2(仮) <!-- N-12 -->).

【訳】本設計は候補交絡因子3つを調整でなく開示する(Methods)。開示された量が、各候補の説明可能範囲を画定する。手術時年齢は群間で異なり — Hodges–Lehmann 中央値差は RET 層で +2.5 年 [−1.0〜6.0]、BRAF 層で +8.0 年 [3.0〜12.0] <!-- N-64, N-65 --> — その射程は4つの観察が画定する。第一に、甲状腺で確立している年齢効果は正常組織の濾胞細胞回転に関するものであり(Coclet et al. 1989)、発見された差は腫瘍にある — 腫瘍では、ここの RET 層の差より一桁広いおよそ 25 年幅の年齢対比でも、ゲノムワイドな発現差もセットレベル所見も生まれなかった(Vriens et al. 2011)。第二に、被曝 RET 帯内で年齢の順序はスコアの順序と逆行する: 手術時年齢の中央値は Low / Mid / High で 30 / 25 / 23 歳と下がり <!-- N-12 -->、逆転スコアは 1 / 4 / 6 と上がる <!-- N-41 -->。第三に、純度はアーティファクト説と逆向きに働く: High 帯は最も高純度だが(0.814 対 Sporadic 0.690 <!-- N-44 -->)、増殖プログラムは被曝群でむしろ減弱として読まれる <!-- N-61 -->。第四に、内部対照: 年齢差が最大の層は腫瘍側シグナルを一切生まなかった(B_Tumor <!-- N-16, N-27 -->)。約3分の1の年齢差しか持たない R_Tumor の差を年齢に帰属させることは、本プロトコル自身の BRAF 層と矛盾する。したがって年齢は(不可分の潜伏期とともに)生きている場所 — 上の B_Normal の読みと後述の注釈 — でのみ候補として名指しされ、被曝時年齢が assigned share の定義の構成要素であるため共変量には入れない(Methods)。残る2候補には同等の画定は必要ない: 性は上の B_Normal で名指し済みであり、融合パートナー構成は帯を跨ぐ勾配を示さない(表2(仮)<!-- N-12 -->)。

### Reading the over-representation annotation(C-14)

<!-- C-14(Vriens・年齢の詳細は候補交絡因子節へ寄せ、ここは参照) -->

The over-representation annotation of the R_Tumor list is read under a gene-sampling null that ignores co-expression and is anticonservative relative to the calibrated set-level inference; its flags are candidates supported only under that weaker null. Three observations bound its meaning. First, the flags across all four families express a single theme — relative attenuation of proliferation, cell-cycle and DNA-repair programs in the exposed group — rather than several independent lines of evidence, and the radiation-curated flags ride the same cell-cycle content. Second, the calibrated set-level test ranked the same Hallmark sets first without crossing its threshold <!-- N-28 -->: the two procedures disagree about evidence, not about ordering. Third, whether the attenuation reflects exposure-associated biology or the groups' age and latency structure cannot be decided here; the age bounds above apply unchanged. No claim of the primary analyses rests on this annotation.

【訳】R_Tumor リストの過剰代表注釈は、共発現を無視する gene-sampling 帰無の下で読まれ、較正済みセットレベル推論に対して反保守的である。そのフラグは、この弱い帰無の下でのみ支持される候補である。その意味の範囲を3つの観察が画定する。第一に、4ファミリー全てのフラグは複数の独立した証拠線ではなく単一のテーマ — 被曝群における増殖・細胞周期・DNA 修復プログラムの相対的減弱 — を表現しており、放射線キュレーションのフラグも同じ細胞周期の内容に乗っている。第二に、較正済みセットレベル検定は同じ Hallmark セット群を、閾値を越えないまま最上位に置いた <!-- N-28 -->: 二つの手続きが食い違うのは証拠についてであり、順位についてではない。第三に、この減弱が曝露関連の生物学を反映するのか、群の年齢・潜伏期構造を反映するのかはここでは決められない — 上の年齢の画定がそのまま適用される。主解析のどの主張もこの注釈に依存しない。

### Concordance across strata(C-17, C-07)

<!-- C-17(tumor: 仮説生成・純度筆頭の識別不能・B_Tumor×radiation q0.114 の断片束ね)+ C-07(normal: 仮説担持・識別不能の宣言)を統合。読みの水準が対になる -->

The two cross-stratum comparisons carry different weights and return different verdicts. The tumor pair, computed as a design completion, resembles each other beyond label-exchange chance: Spearman rho +0.459 over 15,560 shared genes, outside the label-shuffle reference interval ([−0.39, +0.39]; two-sided p = 0.0197 <!-- N-33 -->; Supp. Tab. 3(仮)) <!-- C-17 -->. As hypothesis generation only: this is what an exposure-associated component visible across driver backgrounds — an oncogene-independent trace — would look like, and it is compatible with B_Tumor's within-contrast quiet, which requires only sub-threshold structure aligned with the R_Tumor contrast. It is equally compatible with shared covariate structure, purity foremost (the exposed groups sit purer at least on the RET side <!-- N-44 -->), and the design cannot separate the two; an adjacent fragment — the closest any calibrated set-level cell came to discovery, B_Tumor × radiation (q = 0.114 <!-- N-28 -->) — points the same way without crossing any threshold. This comparison does not test the shared-glandular-memory hypothesis, and no claim of this study rests on it. The normal-tissue pair, by contrast, is the comparison that carried the hypothesis: long-term molecular memory of exposure (Abend et al. 2013) and shared normal-tissue signatures reported without driver stratification (Ory et al. 2026) predict that the two driver-conditioned normal contrasts should resemble each other, and this comparison is the device that makes that prediction testable under driver conditioning. It returned indeterminacy: rho +0.376, inside its reference interval <!-- N-34 -->, with no within-contrast signal in R_Normal <!-- N-16 -->. A null cross-stratum correlation between contrasts that individually carry no detectable within-contrast signal cannot distinguish absent sharing from sharing too weak to detect; whether the two strata share a normal-tissue exposure trace is therefore not identifiable in this cohort <!-- C-07 -->, rather than shown to be absent; the comparison neither supports nor contradicts the shared-memory reports, and it remains the design under which a larger cohort could decide the question.

【訳】二つの層間比較は異なる重みを持ち、異なる判定を返した。腫瘍ペアは設計上の補完として計算されたもので、ラベル交換の偶然を超えて互いに似ている: 共有 15,560 遺伝子で Spearman rho +0.459、ラベルシャッフルの参照区間([−0.39, +0.39])の外側、両側 p = 0.0197 <!-- N-33 -->(Supp. Tab. 3(仮))<!-- C-17 -->。仮説生成としてのみ: これは driver 背景を越えて見える曝露関連成分 — oncogene 非依存の痕跡 — が示すはずの姿であり、B_Tumor の対比内の静けさとも矛盾しない(R_Tumor の対比と向きの揃った閾値下構造があれば足りる)。同時に、共有された共変量構造 — 筆頭は純度(少なくとも RET 側で被曝群は高純度 <!-- N-44 -->)— とも等しく整合し、本設計は両者を分離できない。隣接する断片 — 較正済みセットレベルの全セルで発見に最も近づいた B_Tumor × radiation(q = 0.114 <!-- N-28 -->)— も、どの閾値も越えないまま同じ方向を指す。この比較は共有腺記憶仮説を検定せず、本研究のどの主張もこれに依存しない。対照的に、正常組織ペアは仮説を担った比較である: 被曝の長期分子記憶の報告(Abend et al. 2013)と、driver 非層別で報告された正常組織の共有署名(Ory et al. 2026)は、driver で条件付けた二つの正常対比が互いに似ることを予測し、この比較はその予測を driver 条件付けの下で検定可能にする装置である。結果は不確定に終わった: rho +0.376 は参照区間の内側 <!-- N-34 -->、R_Normal に対比内シグナルなし <!-- N-16 -->。各対比が単独では検出可能な対比内シグナルを持たないとき、層間相関の帰無は「共有の不在」と「検出には弱すぎる共有」を区別できない。したがって二つの層が正常組織の曝露痕跡を共有するか否かは、本コホートでは不在と示されたのではなく識別できない <!-- C-07 -->。この比較は共有記憶の報告を支持も否定もしない — そして、より大きなコホートがこの問いを決められる設計として残る。

### External anchors(C-13)

<!-- C-13(読み側 — 事実は Results §7)。対称読みの事前固定 -->

The external anchor counts (Results) are read as description under the symmetric rule fixed in advance. Zero overlap with the qRT-PCR-validated cores is the norm rather than the exception among the prior studies themselves — Dom et al. (2012) found none of the earlier signatures enriched in their own contrast — and the single cross-tissue overlap (S100A10, direction reversed <!-- N-54 -->) is noted without weight. No claim of this study moves with these counts <!-- C-13 -->.

【訳】外部アンカーの員数(Results)は、事前に固定した対称規則の下で記述として読む。qRT-PCR 検証済みコアとの重なりゼロは、先行研究どうしの間でも例外でなく通例であり — Dom et al. (2012) 自身が先行署名の非濃縮を報告している — 唯一の組織対応外の重なり(S100A10、方向逆転 <!-- N-54 -->)は重みを付けずに記す。本研究のどの主張もこの員数では動かない <!-- C-13 -->。

### Limitations(C-11)

<!-- C-11: A/B 非識別性の明文化。ほか: ウェット追試不能・REO の中間データ依存、被曝時年齢の群間比較不能(N-63) -->

Three limitations are structural. First, the design is cross-sectional over driver-conditioned strata: an exposure-associated expression difference conditioned on driver cannot, in principle, distinguish a shared trace that feeds RET oncogenesis from an association between trace type and driver selection, and every claim in this paper is worded neutrally between the two <!-- C-11 -->. Second, two practical constraints cap claim strength: no wet-laboratory replication is available, and the REO panel is defined on intermediate data products rather than raw counts alone. Third, disclosed asymmetries qualify generalization: age at exposure exists only for exposed cases <!-- N-63 -->, the exposed BRAF group is small (9 cases <!-- N-09 -->), tumor purity is a within-cohort relative quantity, and the single held-out calibration excess (B_Normal/Hallmark <!-- N-25 -->) travels with any reading of that cell.

【訳】3つの限界は構造的である。第一に、本デザインは driver 条件付き層の上の横断研究である: driver で条件付けた曝露関連の発現差は、RET 発がん機序に接続する共有痕跡と、痕跡タイプと driver 選択の連関とを原理的に識別できず、本論文の全ての主張は両者の間で中立に文言化されている <!-- C-11 -->。第二に、二つの実務的制約が主張の強さの上限を画す: ウェット実験による追試は利用できず、REO パネルは生カウント単独でなく中間データ産物の上に定義されている。第三に、開示済みの非対称が一般化を限定する: 被曝時年齢は被曝症例にのみ存在し <!-- N-63 -->、BRAF の被曝群は小さく(9 例 <!-- N-09 -->)、腫瘍純度はコホート内の相対量であり、held-out 較正の唯一の超過(B_Normal/Hallmark <!-- N-25 -->)は当該セルのどの読みにも付随する。

## Supplementary Methods / Results

### Supplementary Methods — Calibration of the gene-set level inference

<!-- C-06 -->
Because gene-set q-values inherit the dependence structure of the expression matrix — genes are co-expressed, and sets overlap — nominal false discovery rates at the set level cannot be taken on faith. We therefore measured the operating characteristics of the entire set-level procedure on held-out null replicates before it was applied to the real contrasts, and we report that measurement alongside the results it calibrates.

For each analysis contrast, exposure labels were shuffled B + R times with a seed independent of the one used for inference (base seed 19450809; B = 9999, R = 100 <!-- N-06 -->), and the block enrichment statistic was computed once for all shuffles. The first B shuffles form a shared null pool; each of the R remaining shuffles was then treated as a pseudo-observation and pushed through exactly the inference applied to the real data — normalized enrichment against the pool, per-set permutation p-values, Benjamini–Hochberg adjustment within each collection, and the pre-specified threshold q < 0.10
<!-- N-06 -->. Under label exchange the pseudo-observations are exchangeable with the pool, so for every contrast × collection cell the share of replicates yielding at least one discovery estimates P(≥1 false discovery) under the complete null, where FDR and family-wise error coincide. We report this share with an exact binomial confidence interval and the mean discovery count per replicate (Supplementary Table SX). Because the R replicates share one null pool, their discovery indicators are weakly positively correlated and the binomial interval is accordingly somewhat narrow; we record this rather than correct it. The calibration also fixed the choice of the set-level inference itself, before any real-data run. The tail-ratio FDR estimated from pooled normalized enrichment scores — the inference originally specified — was measured miscalibrated on these data (pooled P(≥1) 0.140, and 0.221 for a restandardized variant, against the nominal 0.10; worst cell 0.44 <!-- N-58 -->), and was replaced by the per-set permutation p with
within-collection Benjamini–Hochberg adjustment, which measured 0.045 in the same pre-real-data setting <!-- N-57 -->.

【訳】遺伝子セットの q 値は発現行列の依存構造 — 遺伝子は共発現し、セットは重複する — を継承するため、セットレベルの名目 FDR を額面どおりに信用することはできない。そこで我々は、実対比への適用に先立ち、セットレベル手続き全体の動作特性を held-out 帰無レプリケートで測定し、その測定を較正対象の結果と並べて報告する。
各解析対比で、推論に使うものと独立のシードにより曝露ラベルを B + R 回シャッフルし(基底シード 19450809; B = 9999、R = 100 <!-- N-06 -->)、block 濃縮統計量を全シャッフルについて一度に計算した。最初の B 本が共有帰無プールをなし、残る R 本の各々を擬似観測として、実データに適用するのと厳密に同一の推論 — プールに対する正規化濃縮、セット別置換 p 値、各コレクション内の Benjamini–Hochberg 調整、事前指定閾値 q < 0.10 <!-- N-06 --> — に通した。ラベル交換の下で擬似観測はプールと交換可能なので、各対比 × collection セルについて「発見を1つ以上生んだレプリケートの割合」は、FDR と family-wise 誤りが一致する完全帰無下の P(偽発見 ≥ 1) を推定する。この割合を正確二項信頼区間およびレプリケートあたり平均発見数とともに報告する(Supplementary Table SX)。R 本のレプリケートは1つの帰無プールを共有するため発見指標は弱い正の相関を持ち、二項区間はその分やや狭い — これは補正せず記録する。この較正は、実データ実行に先立ちセットレベル推論そのものの選択も固定した。当初指定されていた、pooled 正規化濃縮スコアからの tail-ratio FDR は、このデータ上で較正不良と実測され(pooled P(≥1) 0.140、再標準化変種で 0.221、名目 0.10 に対して; 最悪セル 0.44 <!-- N-58 -->)、コレクション内 Benjamini–Hochberg を伴うセット別置換 p に置き換えられた — 後者は同じ実データ前の設定で 0.045 と測定された <!-- N-57 -->。

### Supplementary Results — 較正表の読み(Supp Table SX に併記する注記)

<!-- C-06 -->
Across the 16 contrast × collection cells, the share of null replicates producing at least one discovery at q < 0.10 ranged from 0.01 to 0.18 <!-- N-24 -->, and 13 of 16 cells were at or below the nominal 0.10. Two further cells in the radiation-curated family straddled it (B_Tumor 0.10, CI 0.049–0.176; B_Normal 0.12, CI 0.064–0.200 <!-- N-26 -->). Together with the spike-in recovery (NES 2.28, p 0.0002 <!-- N-31 -->), the two controls bound the machinery from both sides: null inputs do not generate discoveries beyond the nominal level, and a planted coherent signal of modest size is detected. The one calibration excess (B_Normal/Hallmark) is disclosed and is taken into account wherever set-level results for that cell are read.

【訳】16 の対比 × collection セル全体で、q < 0.10 の発見を1つ以上生んだ帰無レプリケートの割合は 0.01 から 0.18 の範囲にあり <!-- N-24 -->、16 セル中 13 は名目 0.10 以下だった。radiation ファミリーのさらに2セルがこれを跨いだ(B_Tumor 0.10、CI 0.049–0.176; B_Normal 0.12、CI 0.064–0.200 <!-- N-26 -->)。spike-in の回収(NES 2.28、p 0.0002 <!-- N-31 -->)と合わせて、2つの対照は機構を両側から挟む: 帰無入力は名目水準を超える発見を生まず、植えられた中程度の協調シグナルは検出される。唯一の較正超過(B_Normal/Hallmark)は開示済みであり、当該セルのセットレベル結果を読むあらゆる場所で考慮される。

---

## 引用文献(整理 2026-08-14 — 仮書式、投稿書式の整合は後工程)

凡例: 全行とも書誌確認済み(PubMed または Crossref 照合、2026-08-14 完了。★は解消)。各行末は挿入位置。

### A. 本文が現に引くもの

- Morton LM, et al. Radiation-related genomic profile of papillary thyroid carcinoma after the Chernobyl accident. Science 2021;372:eabg2538. doi:10.1126/science.abg2538 — コホート出典(Methods: Data sources / Results §1)
- Storey JD. A direct approach to false discovery rates. J R Stat Soc B 2002;64:479–498. doi:10.1111/1467-9868.00346(Crossref 照合済み 2026-08-14)— λ=0.5 の無チューニング既定(初出引用は Methods: Normalization の DEGES スクリーン、選定事由は Gene-level — 本文反映済み)
- Donoho D, Jin J. Higher criticism for detecting sparse heterogeneous mixtures. Ann Stat 2004;32:962–994. doi:10.1214/009053604000000265 — HC の選定事由(Methods: Gene-level)
- Wasserstein RL, Lazar NA. The ASA statement on p-values. Am Stat 2016;70:129–133. doi:10.1080/00031305.2016.1154108 — 非2値化(Methods: Gene-level)
- Greenland S, Senn SJ, Rothman KJ, Carlin JB, Poole C, Goodman SN, Altman DG. Eur J Epidemiol 2016;31:337–350. doi:10.1007/s10654-016-0149-3 — 同上
- Amrhein V, Greenland S, McShane B. Scientists rise up against statistical significance. Nature 2019;567:305–307. doi:10.1038/d41586-019-00857-9 — 同上
- Goeman JJ, Bühlmann P. Analyzing gene expression data in terms of gene sets: methodological issues. Bioinformatics 2007;23:980–987. doi:10.1093/bioinformatics/btm051 — sampling-model 軸・GSEA/ORA の濃淡開示(Methods: Gene-set)
- Vriens MR, et al. Cancer 2011;117:259–267. doi:10.1002/cncr.25369 — 腫瘍側の年齢対比(Discussion: C-14 段落。全文確認済み 2026-08-12)

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
- Shen Q, et al. contamDE: differential expression analysis of RNA-seq data for contaminated tumor samples. Bioinformatics 2016;32:705–712. doi:10.1093/bioinformatics/btv657 — 純度推定〔手法引用/自作実装〕(本文反映済み)

被曝指標・その他:

- Kocher DC, et al. Interactive RadioEpidemiological Program (IREP): a web-based tool for estimating probability of causation/assigned share. Health Phys 2008;95:119–147. doi:10.1097/01.HP.0000291191.49583.f7 — NIH IREP(Methods: Exposure metric。使用版 5.7.3)
- Geman D, et al. Classifying gene expression profiles from pairwise mRNA comparisons. Stat Appl Genet Mol Biol 2004;3:Article19. doi:10.2202/1544-6115.1071 — 検体内順位比較(REO)の原型〔概念引用/パネル構築規則は自作〕(Methods: REO — 本文反映済み)
- Wang H, Sun Q, Zhao W, Qi L, et al. Individual-level analysis of differential expression of genes and pathways for personalized medicine. Bioinformatics 2015;31:62–68. doi:10.1093/bioinformatics/btu522(PubMed/Crossref 照合済み 2026-08-15、pmid 25165092)— REO の用語・安定参照+逆転スキーム(RankComp)・正規化フリー性の系譜〔概念引用/パネル構築規則は自作〕(Methods: REO — 本文反映済み)
- Nikiforov YE, Rowland JM, Bove KE, Monforte-Munoz H, Fagin JA. Distinct pattern of ret oncogene rearrangements in morphological variants of radiation-induced and sporadic thyroid papillary carcinomas in children. Cancer Res 1997;57:1690–1694.(PubMed 照合済み 2026-08-15、pmid 9135009、DOI なし)— 放射線関連 PTC とパートナー特異性(RET/PTC3 優位)の原典・旧称の橋渡し(Methods: Candidate confounders — 本文反映済み)
- Rabes HM, Demidchik EP, Sidorow JD, et al. Pattern of radiation-induced RET and NTRK1 rearrangements in 191 post-Chernobyl papillary thyroid carcinomas: biological, phenotypic, and clinical implications. Clin Cancer Res 2000;6:1093–1103.(PubMed 照合済み 2026-08-15、pmid 10741739、DOI なし)— 同上・最大古典シリーズ(Methods: Candidate confounders — 本文反映済み)
- R Core Team. R: A language and environment for statistical computing. v4.5.3, 2026(citation() 形式で執筆時に確定)— (Methods: Software)
- その他の使用パッケージ(SummarizedExperiment・limma・GenomicDataCommons・rtracklayer・data.table・Matrix・Rcpp・MASS・parallel)は Supp の版表が網羅 —個別引用は原典論文を持つもののみ執筆時に判断〔使用パッケージ〕。事前収集済みの原典(PubMed 照合済み 2026-08-14): Bioconductor/SummarizedExperiment = Huber W, et al. Nat Methods 2015;12:115–121, doi:10.1038/nmeth.3252 / limma = Ritchie ME, et al. Nucleic Acids Res 2015;43:e47, doi:10.1093/nar/gkv007 / rtracklayer = Lawrence M, et al. Bioinformatics 2009;25:1841–1842, doi:10.1093/bioinformatics/btp328

### C. Intro / Discussion で必要になる見込み

外部アンカー(C-13。DOI は diagnostics/external_gene_anchors.csv の curation 由来):

- Abend M, et al. PLoS One 2012;7:e39103. doi:10.1371/journal.pone.0039103(PubMed 照合済み)— ペア差リスト(11 遺伝子)
- Abend M, et al. Iodine-131 dose-dependent gene expression. Br J Cancer 2013;109:2286–2294. doi:10.1038/bjc.2013.574(PubMed 照合済み)— 正常8/腫瘍6
- Dom G, et al. A gene expression signature distinguishes normal tissues of sporadic and radiation-induced papillary thyroid carcinomas. Br J Cancer 2012;107:994–1000. doi:10.1038/bjc.2012.302(PubMed 照合済み)— 正常7 (SERPINE1 のメンバー番号は原文 PDF で凍結前に要確認 — 既存の札)
- Hess J, et al. Gain of chromosome band 7q11 in papillary thyroid carcinomas of young patients is associated with exposure to low-dose irradiation. PNAS 2011;108:9595–9600. doi:10.1073/pnas.1017137108(PubMed 照合済み)— CLIP2 (追加検証: doi:10.1038/onc.2014.311)

正常組織・対照・仮説:

- Ory C, et al. Post-Chornobyl thyroid papillary carcinomas display distinct past exposure and radiation-associated carcinogenesis molecular signatures at low and high thyroid doses. Sci Rep 2026;16. doi:10.1038/s41598-026-53030-4(PubMed 照合済み)— 正常組織対照・driver 非層別の隙間(Intro 予告 / Discussion C-07)
- Coclet J, et al. Cell population kinetics in dog and human adult thyroid. Clin Endocrinol (Oxf) 1989;31:655–665. doi:10.1111/j.1365-2265.1989.tb01290.x(PubMed 照合済み)— 正常甲状腺の細胞回転(Discussion: C-16 段落・年齢段落)
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
- 年齢・潜伏期の受け: **判断点4決着により確定(2026-08-14 反映)** — 年齢層別診断は実施せず(Q-03・撤回済み3案は復活させない)、開示された推定(C-15、N-64)+文献側(Vriens 2011、Q-13 (i))で受ける。Discussion の [ ] は埋め済み(N-64 タグ)。潜伏期は年齢と不可分のため候補説明としての名指しのまま。
- Supp Table SY = Supp.Tab.2(仮)(図表台帳)。D6 の Supp Table SX と番号整合は図表構成確定時に振り直し。
- リント自己検査: 仮説生成の水準明示あり / A/B 中立(機構主張なし)/ 線量依存なし / q<0.25 不使用 / 対比横断 FDR なし / 全数値 N タグ(NES 個別引用は保留)/ C-14 タグ済み。
