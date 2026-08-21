# Supplementary Material(Methods / Results — 投稿時の単一別ファイル)

状態: draft(2026-08-21 に draft_manuscript.md から分離 — BJC の SI 単一ファイル規定に対応する投稿構造化。規約は本文と共通: 英語が正本・数値は N-ID タグ・主張は C-ID タグ。C2 ポートで gpt_review 試案の文面に順次置換予定)。

---

## Supplementary Methods / Results

### Supplementary Methods — Data sources and expression matrix

Gene lengths (exon-union length per gene) were derived from the GDC GENCODE v36 reference annotation <!-- N-66 -->. In the clinical table, missing markers were coerced to NA and columns were typed numerically only where every non-missing value parses. Library strandedness was decided per sample from the STAR stranded-column totals (ratio rule: the library is called stranded when the smaller stranded total is at most half the larger, and the larger column is used; otherwise the unstranded column is used <!-- N-67 -->). In practice the call was unambiguous — all 906 libraries were stranded (reverse), with ratios of 0.056–0.110 against the 0.5 threshold <!-- N-77 -->. TPM and FPKM columns were discarded at loading.

【訳】遺伝子長(遺伝子ごとの exon-union 長)は GDC の GENCODE v36 参照アノテーションから導出した <!-- N-66 -->。臨床表では欠損マーカーを NA に統一し、非欠損の全値が数値として解釈できる列のみ数値型とした。ライブラリの strand は検体ごとに STAR の stranded 2列の合計から判定した(比率規則: 小さい方の合計が大きい方の半分以下なら stranded と判定して大きい列を採用、それ以外は unstranded 列を採用 <!-- N-67 -->)。実際には判定は一義的で、全 906 ライブラリが stranded(reverse)判定、比は 0.056–0.110 と閾値 0.5 から大きく離れていた <!-- N-77 -->。TPM・FPKM 列は読み込み時に破棄した。

### Supplementary Methods — Exposure metric: assigned share

The fixed input convention is: electrons E > 15 keV; exposure rate acute; organ dose entered as a constant, in cSv = mGy/10; sex as recorded; exposure year 1986; birth year = 1986 − age at exposure; surgery year = birth year + age at surgery; remaining IREP assumptions and settings at their defaults (user-defined uncertainty distribution Lognormal (1,1); 10,000 iterations; random number seed 99) <!-- N-68 -->. The AS values are carried as a versioned input file with the analysis code; their dose and age inputs were verified case-by-case against the clinical table. On the disclosed properties: IREP's radiation-type guidance names iodine-131 explicitly as an exemplar for the electrons E > 15 keV entry, and the radiation-effectiveness equivalence is independently supported by microdosimetric analysis of internal iodine-131 against external γ-rays (Sato et al. 2014). The thyroid risk model behind IREP is fitted to atomic-bomb survivors and children irradiated in medical treatments (Ron et al. 1995; Kocher et al. 2008). The acute setting fixes the dose and dose-rate effectiveness factor at 1 for doses at or above the model's uncertain reference dose (log-uniform on 0.03–0.2 Gy; Land et al. 2003); the High band's minimum dose sits at the top of that range. The Sporadic band requires no model; computed under the chronic setting, the same partition would simply carry lower numeric labels.

【訳】固定した入力規約は次のとおり: 電子 E > 15 keV; 曝露率 acute; 臓器線量は定数(constant)として cSv = mGy/10 で入力; 性別は記録値のまま; 被曝年 1986; 出生年 = 1986 − 被曝時年齢; 手術年 = 出生年 + 手術時年齢; 残りの IREP の仮定・設定は既定値(ユーザー定義不確かさ分布 Lognormal (1,1)、反復 10,000 回、乱数シード 99)<!-- N-68 -->。AS 値は解析コードとともに版管理された入力ファイルとして持ち回り、線量・年齢の入力は臨床表と全例照合した。開示3性質の細目: IREP の放射線タイプ選択ガイダンスは electrons E > 15 keV の項でヨウ素131を例示として明示的に名指しし、放射線有効性の等価性は内部ヨウ素131と外部γ線のマイクロドジメトリ解析が独立に支持している(Sato et al. 2014)。IREP の背後の甲状腺リスクモデルは原爆被爆者と医療照射を受けた小児(Ron et al. 1995; Kocher et al. 2008)に当てはめられたものである。acute 設定は、モデルの不確かな参照線量(0.03–0.2 Gy の対数一様; Land et al. 2003)以上の線量で線量・線量率有効性係数を 1 に固定する。High 帯の最小線量はこの範囲の上端に位置する。Sporadic 帯はモデルを要さない。chronic 設定で計算していれば、同じ分割がより低い数値ラベルを持つだけである。

### Supplementary Methods — Quality control and analysis cohorts

The outlier screen operates within each group × tissue sub-matrix (filterByExpr-reduced; unnormalized log-CPM, since raw library sizes preserve the composition outliers normalization would mask): the most extreme sample on the first principal component was removed iteratively (two-sided, α = 0.05) until no sample was rejected <!-- N-78 -->. The purity implementation is the maximum-one proportion estimator of contamDE-lm — an in-house implementation of the proportion step only, without the accompanying differential-expression model — on protein-coding, filterByExpr-reduced paired counts; purities are relative within the set of pairs estimated jointly and are reported as a maximum-one relative score within each cohort. The common tumor reference this assumes rests on the same premise as the driver stratification itself — tumor expression in both groups is dominated by the driver's biology — and the tumor-versus-normal contamination axis that purity measures is a different axis from the exposure contrast under test. Case pairing resolves, per case, the GDC merged-aliquot expression assays (sample identifiers suffixed _merged; one merged tumor and one merged normal, enforced unique) — a case lacking either is unpaired.

【訳】外れ値スクリーンは群 × 組織の各部分行列(filterByExpr 縮約; 未正規化 log-CPM — 生のライブラリサイズは正規化が隠す組成外れ値を保存する)で動作する: 第1主成分上で最も極端な検体を反復的に除去し(両側、α = 0.05)、棄却が出なくなるまで続けた <!-- N-78 -->。純度の実装は contamDE-lm の最大1正規化の比率推定量 — 比率ステップのみの自作実装で、付随する発現差モデルは用いない — を protein-coding・filterByExpr 縮約のペア付きカウントに適用したものである。純度は同時に推定した集合内でのみ相対的であり、コホート内最大値 = 1 の相対スコアとして報告する。この仮定する共通腫瘍参照は driver 層別という設計自体の前提(両群の腫瘍発現は driver の生物学が支配する)に乗っており、純度が測る腫瘍対正常の混入軸は検定対象の曝露対比とは別の軸である。症例のペア解決は GDC の merged-aliquot 発現アッセイ(検体 ID 接尾辞 _merged; 症例ごとに merged 腫瘍1・merged 正常1、一意性を強制)を用い、どちらかを欠く症例はペアなしとする。

### Supplementary Methods — Candidate confounders

The age difference is estimated by the Hodges–Lehmann median difference and the rank-based effect P(Sporadic < High) — the protocol's Brunner–Munzel effect estimator — each with a 95% percentile bootstrap confidence interval (within-group resampling, 9,999 replicates, seed 19450809). Age at exposure is structurally not comparable between groups (the Sporadic group is unexposed). The assigned-share values are numerically insensitive to the recorded sex input, so sex is not absorbed into the exposure metric.

【訳】年齢差は Hodges–Lehmann 中央値差と順位ベースの効果 P(Sporadic < High)(本プロトコルの Brunner–Munzel 効果推定量)で推定し、それぞれに 95% percentile ブートストラップ信頼区間を付す(群内復元抽出、9,999 レプリケート、シード 19450809)。被曝時年齢は構造的に群間比較が成立しない(Sporadic 群は非被曝)。assigned share の値は入力した性別に対して数値的に不感応であり、性は曝露指標に吸収されない。

### Supplementary Methods — Normalization

The MUREN step uses the single-parameter form (pairwise least-trimmed-squares regression). The Brunner–Munzel screen uses the same Storey estimator as the gene-level inference, with the TCC floorPDEG guard: at each iteration, the larger of the q-threshold set and the top 5% of genes by raw p is removed from the scaling-factor estimation.

【訳】MUREN ステップは単一パラメータ形(pairwise least-trimmed-squares 回帰)を用いる。Brunner–Munzel スクリーンは遺伝子レベル推論と同一の Storey 推定量を用い、TCC の floorPDEG ガードを併用する: 各反復で、q 閾値集合と生 p 上位 5% の大きい方をスケーリング係数推定から除去する。

### Supplementary Methods — Gene-level differential expression

The permutation enumeration is provided with the analysis code. The effect P(X<Y) reads directly as an exceedance probability, and at these group sizes the permutation distribution is exactly enumerable. The choice of λ = 0.5 rests on three considerations. First, the conservativeness of the fixed-λ plug-in requires only marginal uniformity of the null p-values, which the exact test guarantees regardless of the dependence between genes, whereas adaptive choices of λ carry guarantees that assume independent or weakly dependent tests and do not apply to co-expressed genes sharing one label vector. Second, under a weak spread signal the alternative p-density is nearly flat, so raising λ buys little bias reduction at first-order variance cost. Third, λ = 0.5 is the fixed value Storey (2002) used for all of that paper's own calculations; the adaptive selection of λ proposed in its Section 9 is a tuning device deliberately not used here.

【訳】置換枚挙の実装はコードと共に提供する。効果 P(X<Y) は超過確率として直接読め、この群サイズでは置換分布が完全枚挙できる。λ = 0.5 の選定は3つの考慮に立つ。第一に、固定 λ plug-in の保守性は帰無 p 値の周辺一様性のみを要件とし、これは exact 検定が遺伝子間の依存に関わらず保証する。一方、λ の適応的選択は独立ないし弱依存の検定を仮定する保証しか持たず、単一のラベルベクトルを共有する共発現遺伝子には適用できない。第二に、弱く広がったシグナルの下では対立 p 密度はほぼ平坦で、λ を上げてもバイアス削減は僅かなのに分散コストは一次である。第三に、λ = 0.5 は Storey (2002) が自身の計算全体で用いた固定値であり、同論文 §9 が提案する適応的 λ 選択はチューニング装置として意図的に用いない。

### Supplementary Methods — Contrast-level omnibus inference

Because the Higher Criticism p-value is taken from the contrast's own label-permutation null, neither the sparsity regime nor the independence assumptions of the original optimality theory are invoked; the original leaves the scan range α0 as a free constant.

【訳】Higher Criticism の p 値は当該対比自身のラベル置換帰無から得るため、原典の最適性理論の疎性領域も独立性仮定も援用しない。原典は走査範囲 α0 を自由定数として残している。

### Supplementary Methods — Gene-set level inference

The weighted running sum uses gseaParam = 1; its equality with the standard GSEA statistic on tie-free input is enforced by an automated test in the analysis code. The shared null comprises 9,999 label shuffles per contrast. The spike-in control multiplied one Hallmark set by 1.15 in the 9 samples of one group <!-- N-30 -->. The over-representation annotation tests the up, down and combined lists, with Benjamini–Hochberg adjustment within family × list. In the taxonomy of Goeman & Bühlmann (2007), the enrichment statistic used here is itself competitive — it compares set members against the rest of the ranking — but its p-value is computed under the subject-sampling label-permutation null.

【訳】重み付き累積和は gseaParam = 1 を用いる。tie のない入力での標準 GSEA 統計量との一致は、解析コード内の自動テストが強制する。共有帰無は対比あたり 9,999 回のラベルシャッフルからなる。spike-in 対照は、1つの Hallmark セットを一方の群の 9 検体で 1.15 倍した <!-- N-30 -->。過剰代表注釈は up・down・合算の3リストを検定し、family × list 内で Benjamini–Hochberg 調整する。Goeman & Bühlmann (2007) の分類では本研究の濃縮統計量自体は competitive — セット員をランキングの残りと比較する — だが、その p 値は subject-sampling のラベル置換帰無で計算している。

### Supplementary Methods — Calibration of the gene-set level inference

<!-- C-06 -->
Because gene-set q-values inherit the dependence structure of the expression matrix — genes are co-expressed, and sets overlap — nominal false discovery rates at the set level cannot be taken on faith. We therefore measured the operating characteristics of the entire set-level procedure on held-out null replicates before it was applied to the real contrasts, and we report that measurement alongside the results it calibrates.

For each analysis contrast, exposure labels were shuffled B + R times with a seed independent of the one used for inference (base seed 19450809; B = 9999, R = 100 <!-- N-06 -->), and the block enrichment statistic was computed once for all shuffles. The first B shuffles form a shared null pool; each of the R remaining shuffles was then treated as a pseudo-observation and pushed through exactly the inference applied to the real data — normalized enrichment against the pool, per-set permutation p-values, Benjamini–Hochberg adjustment within each collection, and the pre-specified threshold q < 0.10
<!-- N-06 -->. Under label exchange the pseudo-observations are exchangeable with the pool, so for every contrast × collection cell the share of replicates yielding at least one discovery estimates P(≥1 false discovery) under the complete null, where FDR and family-wise error coincide. We report this share with an exact binomial confidence interval and the mean discovery count per replicate (Supplementary Table S7). Because the R replicates share one null pool, their discovery indicators are weakly positively correlated and the binomial interval is accordingly somewhat narrow; we record this rather than correct it. The calibration also fixed the choice of the set-level inference itself, before any real-data run. The tail-ratio FDR estimated from pooled normalized enrichment scores — the inference originally specified — was measured miscalibrated on these data (pooled P(≥1) 0.140, and 0.221 for a restandardized variant, against the nominal 0.10; worst cell 0.44 <!-- N-58 -->), and was replaced by the per-set permutation p with
within-collection Benjamini–Hochberg adjustment, which measured 0.045 in the same pre-real-data setting <!-- N-57 -->.

【訳】遺伝子セットの q 値は発現行列の依存構造 — 遺伝子は共発現し、セットは重複する — を継承するため、セットレベルの名目 FDR を額面どおりに信用することはできない。そこで我々は、実対比への適用に先立ち、セットレベル手続き全体の動作特性を held-out 帰無レプリケートで測定し、その測定を較正対象の結果と並べて報告する。
各解析対比で、推論に使うものと独立のシードにより曝露ラベルを B + R 回シャッフルし(基底シード 19450809; B = 9999、R = 100 <!-- N-06 -->)、block 濃縮統計量を全シャッフルについて一度に計算した。最初の B 本が共有帰無プールをなし、残る R 本の各々を擬似観測として、実データに適用するのと厳密に同一の推論 — プールに対する正規化濃縮、セット別置換 p 値、各コレクション内の Benjamini–Hochberg 調整、事前指定閾値 q < 0.10 <!-- N-06 --> — に通した。ラベル交換の下で擬似観測はプールと交換可能なので、各対比 × collection セルについて「発見を1つ以上生んだレプリケートの割合」は、FDR と family-wise 誤りが一致する完全帰無下の P(偽発見 ≥ 1) を推定する。この割合を正確二項信頼区間およびレプリケートあたり平均発見数とともに報告する(Supplementary Table S7)。R 本のレプリケートは1つの帰無プールを共有するため発見指標は弱い正の相関を持ち、二項区間はその分やや狭い — これは補正せず記録する。この較正は、実データ実行に先立ちセットレベル推論そのものの選択も固定した。当初指定されていた、pooled 正規化濃縮スコアからの tail-ratio FDR は、このデータ上で較正不良と実測され(pooled P(≥1) 0.140、再標準化変種で 0.221、名目 0.10 に対して; 最悪セル 0.44 <!-- N-58 -->)、コレクション内 Benjamini–Hochberg を伴うセット別置換 p に置き換えられた — 後者は同じ実データ前の設定で 0.045 と測定された <!-- N-57 -->。

### Supplementary Methods — Between-stratum concordance of the exposure contrast

Labels were shuffled independently within each contrast (9,999 shuffles, per-contrast seeds from base 19450809 <!-- N-73 -->) and the shuffled profiles correlated, giving the null spread of rho when neither stratum carries label-aligned structure.

【訳】対比内で独立にラベルをシャッフルし(9,999 回、対比別シードは基底 19450809 から <!-- N-73 -->)、シャッフル済みプロファイル同士を相関させ、どちらの層もラベル整列構造を持たないときの rho の帰無分布を得る。

### Supplementary Methods — REO panel: construction and out-of-sample evaluation

TPM was recomputed from the selected count assay and the exon-union gene lengths (data acquisition, above), not the discarded GDC columns. The candidate pool is the top 500 genes by |effect − 0.5|. A pair qualifies when its within-sample log2-TPM difference r has a stable sign in Sporadic (dead zone |r| < log2(1.2); at most one exception among non-dead-zone samples; 10th percentile of |r| ≥ log2(1.5)) and reverses in more than 50% but not all of the High samples <!-- N-74 -->, ranked by the shift in median r. The final panel was selected greedily in rank order, excluding gene reuse and pairs whose per-sample r profiles correlate at Spearman ≥ 0.75 with a kept pair <!-- N-75 -->. The reversal score applies the same dead zone.

【訳】TPM は採用カウントアッセイと exon-union 遺伝子長から再計算したもので(前述のデータ取得節)、読み込み時に破棄した GDC 列ではない。候補プールは |effect − 0.5| 上位 500 遺伝子。ペアの適格条件: 検体内 log2-TPM 差 r が Sporadic で安定した符号を持ち(dead zone |r| < log2(1.2); 非 dead-zone 検体で例外は高々1; |r| の第10百分位 ≥ log2(1.5))、High 検体の 50% 超かつ全例未満で逆転すること <!-- N-74 -->。順位は中央値 r のシフトで付けた。最終パネルは順位順の貪欲選定で、遺伝子の再使用と、採用済みペアと Spearman ≥ 0.75 で相関するペアを除外した <!-- N-75 -->。逆転スコアも同じ dead zone を用いる。

### Supplementary Methods — REO evaluation diagnostics

(i) The training outlier screen (PC-OD) was mirrored on the R_Low/R_Mid tumors, and the evaluation was run once on all eligible cases — the screen reports counts only <!-- N-43 -->. (ii) Tumor purity for the evaluation bands was estimated on the same common scale as training by pooling the whole RET cohort in one contamDE-lm run <!-- N-44 -->. (iii) Band and purity are separated by the partial Spearman correlation of band with score conditioning on purity (permutation reference) and by purity-stratified one-sided Brunner–Munzel comparisons of Mid over Low (two strata split at the median tumor purity) <!-- N-46, N-47 -->.

【訳】(i) 訓練側の外れ値スクリーン(PC-OD)を R_Low/R_Mid 腫瘍に鏡映し、評価は適格全例で一度だけ実行 — スクリーンは件数の報告のみ <!-- N-43 -->。(ii) 評価帯の腫瘍純度は、RET コホート全体を1回の contamDE-lm にプールして訓練と同一の共通尺度で推定 <!-- N-44 -->。(iii) 帯と純度は、純度を条件付けた帯–スコアの偏 Spearman 相関(置換参照)と、純度層別(腫瘍純度の中央値で2層)の Mid over Low 片側 Brunner–Munzel 比較で分離する <!-- N-46, N-47 -->。

### Supplementary Methods — External anchor cross-reference

The membership count is k of n list genes among the contrast's tested genes. Abend 2013 contributes per-tissue lists; Abend 2012 contributes a list defined on the tumor-minus-normal paired difference — an estimand with no counterpart contrast here; Dom 2012 is a normal-tissue list.

【訳】員数照合は、対比の検定遺伝子中のリスト遺伝子 n のうち k を数える。Abend 2013 は組織別リスト、Abend 2012 は腫瘍−正常のペア差で定義されたリスト — 本研究に対応する対比を持たない推定対象 — であり、Dom 2012 は正常組織のリストである。

### Supplementary Methods — Software, seeds and reproducibility

The container build is fully date-pinned: an immutable Ubuntu base image (noble-20260410) with its matching apt snapshot, R 4.5.3 built from source against the reference BLAS, and all R packages installed from a dated repository snapshot (CRAN and Bioconductor 3.22; snapshot 2026-04-09) <!-- N-79 -->. The bit-level reproduction contract targets the container's platform (linux/amd64). The in-house statistical components are the exact/Monte-Carlo Brunner–Munzel test, the Storey plug-in estimator, the tie-block enrichment statistic with its permutation q-values, DEGES-MUREN normalization, the contamDE-lm purity score, the principal-component outlier screen and the REO panel machinery. External packages actually used include edgeR, SummarizedExperiment, msigdbr, limma, GenomicDataCommons, rtracklayer, Rcpp and MASS (full versioned list in Supplementary Table S5).

【訳】コンテナビルドは全層を日付固定する: 不変タグの Ubuntu 基底イメージ(noble-20260410)と同日付の apt スナップショット、参照 BLAS に対してソースからビルドした R 4.5.3、日付スナップショット(CRAN + Bioconductor 3.22; 2026-04-09)から導入した全 R パッケージ <!-- N-79 -->。ビットレベルの再現契約はコンテナの対象プラットフォーム(linux/amd64)を対象とする。自作の統計構成要素は、exact/モンテカルロ Brunner–Munzel 検定・Storey plug-in 推定量・tie ブロック濃縮統計量とその置換 q 値・DEGES-MUREN 正規化・contamDE-lm 純度スコア・主成分外れ値スクリーン・REO パネル機構である。実際に使用した外部パッケージは edgeR・SummarizedExperiment・msigdbr・limma・GenomicDataCommons・rtracklayer・Rcpp・MASS ほか(版つき全リストは補足表S5)。

### Supplementary Results — 較正表の読み(補足表S7 に併記する注記)

<!-- C-06 -->
Across the 16 contrast × collection cells, the share of null replicates producing at least one discovery at q < 0.10 ranged from 0.01 to 0.18 <!-- N-24 -->, and 13 of 16 cells were at or below the nominal 0.10. Two further cells in the radiation-curated family straddled it (B_Tumor 0.10, CI 0.049–0.176; B_Normal 0.12, CI 0.064–0.200 <!-- N-26 -->). Together with the spike-in recovery (NES 2.28, p 0.0002 <!-- N-31 -->), the two controls bound the machinery from both sides: null inputs do not generate discoveries beyond the nominal level, and a planted coherent signal of modest size is detected. The one calibration excess (B_Normal/Hallmark) is disclosed and is taken into account wherever set-level results for that cell are read.

【訳】16 の対比 × collection セル全体で、q < 0.10 の発見を1つ以上生んだ帰無レプリケートの割合は 0.01 から 0.18 の範囲にあり <!-- N-24 -->、16 セル中 13 は名目 0.10 以下だった。radiation ファミリーのさらに2セルがこれを跨いだ(B_Tumor 0.10、CI 0.049–0.176; B_Normal 0.12、CI 0.064–0.200 <!-- N-26 -->)。spike-in の回収(NES 2.28、p 0.0002 <!-- N-31 -->)と合わせて、2つの対照は機構を両側から挟む: 帰無入力は名目水準を超える発見を生まず、植えられた中程度の協調シグナルは検出される。唯一の較正超過(B_Normal/Hallmark)は開示済みであり、当該セルのセットレベル結果を読むあらゆる場所で考慮される。

---

