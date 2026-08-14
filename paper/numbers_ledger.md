# 数値台帳(様式 v1 — 2026-08-12 研究者批准)

運用規則:

- 論文(本文・図表・補足)に載る**すべての数値**はここに1行を持ち、本文はN-IDで辿れる。
- 追加計算は原則発生しない前提のため、台帳は一度作って**凍結**する(状態=verified で凍結)。
- 「出典」は機械照合可能な粒度(ファイルパス+ログの検索キー or rds のフィールド)で書く。
- 起草は詰め合わせセッション(handoff_prompt_summary.md)の成果物を流用してよい。
  照合(値↔出典の突き合わせ)は Claude Code が機械的に行い、状態を verified に上げる。
- **追補計算(執筆フェーズ、批准待ち)**: 執筆・査読対応で新しい数値が必要になった場合も、
  既存数値と同じ規律に服す(固定入力のスクリプトで産出し、N 行として出典つきで照合)。
  それ以外の手続きは設けない。決定の記録先は paper/ の台帳群であり、
  **計画v2 への amendment 追記は行わない**(計画v2 はリポジトリ外の作業用メモ相当 —
  新たな決定の記録先にしない。誤読を招く記述の明確化は可 — 2026-08-14 研究者整理。
  正規の記録と履歴は git、正本は論文)。
- **本台帳は足場(検査機構)である**(2026-08-14、3バケツ整理): 論文の全数値の出典照合の
  ための装置であり、採録確定まで併走した後に削除する(履歴は git に残る)。論文に載らないが
  リバイス応答に使う数値は、削除前に手元台帳(objections_ledger)側で出典を自蔵させること
  (例: N-76 は済み)。

## 列定義

| 列 | 内容 |
| --- | --- |
| N-ID | N-01, N-02, … |
| 値 | 数値そのもの(丸め方も確定させる) |
| 定義 | 何の数値か1行 |
| 出典 | ファイルパス + 検索キー(ログ行の文字列 / rds フィールド) |
| 使用箇所 | 図表番号・節(複数可) |
| 状態 | draft / verified(照合済み) |

## パス規約(出典列)

- `run/` = `~/rebc_run_20260809/`(正準保全ディレクトリ)
- `repo/` = `~/rstudio-docker/`(rebuild/r453。成果物は Xeon 産と md5 一致同期済み — 計画v2 B.14)
- rds のフィールド参照は summary_20260809.md §15 の検算コマンド(A/A2/B)で抽出可能。
  rds はホスト R では読めないため、照合は正準イメージ rebc-r453:refblas のコンテナ内で行う。
- 照合実施: 2026-08-12(Claude Code)。全行を下記出典の一次ファイルと直接突き合わせた。

## 台帳本体

### A. 実行環境・設定(Methods)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-01 | rebuild/r453 @ 8eed384(tree clean) | 本番 run のリポジトリ状態 | run/xeon_provenance/git_state.txt:1-2(行1=`8eed384cf45f…`、行2=dirty 数 `0`) | Methods | verified |
| N-02 | R 4.5.3 / Ubuntu 24.04.4 LTS | 実行環境(採用機 Xeon) | run/xeon_provenance/session_info.txt:1-3 | Methods | verified |
| N-03 | 参照 BLAS 3.12.0 / LAPACK 3.12.0 | 線形代数バックエンド | run/xeon_provenance/session_info.txt:6-7 キー「libblas.so.3.12.0」 | Methods | verified |
| N-04 | WORKERS = 4 / N_PERM = 9999 | 並列度・置換回数(再現契約のコミット値) | `git show 8eed384:config.R` 行21「WORKERS <- 4L」・行64「N_PERM <- 9999L」 | Methods | verified |
| N-05 | BM exact・two-sided・n_perm 9999・perm_seed 19860426・Storey q<0.10(plug-in pi0, λ=0.5)・主オムニバス hc(α0=0.1) | 410 検定設定一式 | repo/processed/thyr_expression_test.rds `$config`(§15 A。md5 `a6e083f593d5600d8f4110a74a4d6ede`、run/xeon_results/processed/ の同名と一致) | Methods | verified |
| N-06 | B=9999 共有帰無・R=100 反復・seed 19450809・q<0.10 | D6 帰無較正の設定 | run/xeon_final_20260811/logs/d6_calibration.log:1 キー「B = 9999 shared null」 | Methods | verified |
| N-07 | 410 の置換インデックスを再利用(9999/unit) | 420 の帰無生成方式 | run/xeon_final_20260811/logs/420_test_gene_sets.log:1 キー「Reusing 410's permutation index」(unit 別ハッシュは N-05 の rds `$config$perm_index_hash`) | Methods | verified |

### B. コホート

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-08 | 440 → 248 → 77 → 70 → 69 → 63(n_RET 73→73→34→31→31→27、n_BRAF 175→175→43→39→38→36) | コホートフロー6段(all_cases → driver_classified → band_sporadic_or_high → paired → pcod_clean → purity_pass) | run/xeon_results/logs/230_finalize_analysis_cohorts.log:1-8 キー「Cohort flow (main BM)」 | Methods / フロー図表 | verified |
| N-09 | B_High 9 / B_Sporadic 27 / R_High 15 / R_Sporadic 12 | main BM 63症例の群内訳 | 同 log:9-12 キー「include_main_bm by group」 | Methods / Results | verified |
| N-10 | training 27 / evaluation 36(R_Low 17・R_Mid 19) | REO コホート | 同 log:13 キー「REO training: 27」 | Methods / Results | verified |
| N-11 | 臨床 440 症例; AS 帯 non_exposed 81 / no_reference 146 / (0,33.3) 100 / [33.3,66.6) 82 / [66.6,100] 31; paired 392 / unpaired 48 | 臨床全数と AS バンド構成(整合性検査通過) | run/xeon_results/logs/tab_cohort_composition.log:1-4(Driver×AS帯×pair の全225行は同 log:9-233) | Tab(コホート構成) | verified |
| N-12 | R系 4群(R_Sporadic/R_Low/R_Mid/R_High): n 12/17/19/15、女/男 10/2・13/4・14/5・11/4、手術時年齢中央値[範囲] 20.5[14–27]/30[22–44]/25[17–31]/23[14–31]、被曝時年齢 NA/12[6–19]/3[0–13]/2[0–12]、CCDC6-RET 6/8/12/7、NCOA4-RET 2/5/3/4、RET-OTHER 4/4/4/4 | 症例特性(R系、群×帯) | run/xeon_results/processed/thyr_analysis_cohorts.rds + thyr_clinical.rds を case_submitter_id=REBC_ID で結合(§15 B で再現。2026-08-12 コンテナ内で再実行し一致) | Tab(症例特性) | verified |
| N-13 | B系 2群(B_Sporadic/B_High): n 27/9、女/男 23/4・4/5、手術時年齢中央値[範囲] 24[11–29]/29[19–39]、被曝時年齢 NA/3[0–13]、Driver は全例 BRAF.MutV600E | 症例特性(B系) | 同上(§15 B) | Tab(症例特性) | verified |

### C. 310 正規化

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-14 | 58448 genes × 906 samples | SE(発現行列)の寸法 | run/xeon_results/logs/310_normalize_counts.log:1 キー「SE: 58448 genes x 906 samples」 | Methods | verified |
| N-15 | ユニット別: filterByExpr 後 15621/15466/16002/15863、pi0_hat(iter3) 0.593/0.956/0.943/0.727、収束 Jaccard(iter3) 0.962/0.967/0.993/0.987、norm.factors 範囲 [0.9058,1.1095]/[0.9106,1.0665]/[0.8772,1.0947]/[0.8988,1.0873](順に R_Tumor/R_Normal/B_Tumor/B_Normal) | 310 正規化サマリ(4ユニット) | 同 log:9-32(各「Unit <名>」ブロックの iteration 3 行と norm.factors 行) | Methods / Supp | verified |

### D. 410 遺伝子別検定(主結果)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-16 | 1765 / 0 / 0 / 1(R_Tumor/R_Normal/B_Tumor/B_Normal) | DEG 数(Storey q<0.10) | repo/processed/thyr_expression_test.rds `units[[u]]$genes$q_storey < 0.10`(§15 A); ログ照合 run/xeon_results/logs/fig_gene_bm_evidence.log:2-5 | 主結果 / Fig / Abstract | verified |
| N-17 | R_Tumor: up 971 / down 794; B_Normal: up 1 / down 0(up = effect>0.5 = High 群で高発現) | DEG の方向内訳 | 同 rds `genes$effect`(§15 A2) | Results | verified |
| N-18 | 4.03e-06 / 4.08e-05 / 1.68e-05 / 1.13e-06 | ユニット別 min p_exact | 同 rds `min(genes$p_exact)`(§15 A); ログ照合 fig_gene_bm_evidence.log:2-5 | Results / Supp | verified |
| N-19 | 0.593 / 0.955 / 0.943 / 0.727 | Storey pi0(410、plug-in λ=0.5) | 同 rds `units[[u]]$pi0$estimate`(§15 A) | Results / Supp | verified |
| N-20 | 0.0112 / 0.3199 / 0.1815 / 0.0773 | 主オムニバス HC p(α0=0.1) | 同 rds `units[[u]]$omnibus` の hc 行(§15 A) | 主結果 / Abstract | verified |
| N-21 | オムニバス全表: 4ユニット × {count α=1e-2/1e-3/1e-4, max, hc} の observed / 帰無中央値 / p(R_Tumor count@1e-2: 1026/82/0.0128 … B_Normal hc: 36.31/−0.41/0.0773。全20行は出典参照) | オムニバス検定の完全な表 | 同 rds `units[[u]]$omnibus`(§15 A の print 出力) | Supp | verified |
| N-22 | ENSG00000198908.12(BHLHB9)、effect 0.967、q_storey 0.0130、rank 1 | B_Normal 唯一の DEG | 同 rds B_Normal `genes`(q_storey<0.10 の1行。effect 0.9670782、q 0.01298213); symbol は raw の STAR counts tsv 行「ENSG00000198908.12→BHLHB9(protein_coding)」 | Results | verified |
| N-23 | 0.0125 | B_Normal オムニバスの最小 p(max 統計量; observed 17.70、帰無中央値 6.99) | 同 rds B_Normal `$omnibus` の max 行(§15 A) | Results | verified |

### E. D6 帰無較正(GSEA 較正)

正準は Xeon 産(2026-08-11、B.14)。i9 の 8/10 再計算版(run/i9canon/d6_calibration.log)と中身一致。

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-24 | 16セル全表(4ユニット×4コレクション: m_sets、p_any、95%CI、mean/max 発見数。p_any 範囲 0.01–0.18) | D6 較正の完全な表 | run/xeon_final_20260811/logs/d6_calibration.log:6-39(rds は repo/diagnostics/output/gsea_null_calibration.rds) | Supp | verified |
| N-25 | 0.18(95%CI 0.110–0.270) | 唯一の名目 0.10 超過: B_Normal/H の p_any(CI 下限も 0.10 超) | 同 log:19 キー「13 B_Normal H 50 100 18 0.18 0.110311」・log:36(ci_hi 0.2695) | Results / Disc | verified |
| N-26 | B_Tumor/radiation 0.10(CI 0.049–0.176)、B_Normal/radiation 0.12(CI 0.064–0.200) | CI が 0.10 を跨ぐ境界セル | 同 log:18・22(ci_hi は log:35・39) | Supp / Disc | verified |

### F. 420 遺伝子セット検定

正準は Xeon 産(2026-08-11、B.14)。i9 の 8/10 再計算版(run/i9canon/420_test_gene_sets.log)と中身一致。

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-27 | 0(16/16 セル) | 420 本試験の発見数(q<0.10、全 unit × collection) | run/xeon_final_20260811/logs/420_test_gene_sets.log:5-8(各行「q<0.10 0」×4) | 主結果 / Abstract | verified |
| N-28 | min q 表: R_Tumor(6192セット) H 0.179/C2:CP 0.928/C5:GO:BP 0.936/rad 0.302; R_Normal(6141) 0.380/1.000/1.000/0.980; B_Tumor(6242) 0.520/0.892/0.262/0.114; B_Normal(6223) 0.284/0.283/0.467/0.137 | ユニット×コレクションの min q と総セット数 | 同 log:5-8 | Supp | verified |
| N-29 | H=50, C2:CP=3910, C5:GO:BP=7538, C2:CGP:radiation=28 | フィルタ前のセット数 | 同 log:2-3 キー「Gene sets before filtering」 | Methods | verified |
| N-55 | C2:CP 検定済み内訳(R_Tumor, 2262セット): Reactome 1095 / WikiPathways 633 / KEGG MEDICUS 217 / BioCarta 134 / PID 183(未分類 0)。msigdbr 26.1.0 | C2:CP の構成(=Reactome・KEGG を含む統合であることの数値根拠。他 unit はフィルタ差で僅差) | 定義: lib/gsea_collections.R `GSEA_COLLECTION_SPECS`(CP:REACTOME/CP:WIKIPATHWAYS/CP:KEGG_MEDICUS/CP:BIOCARTA/CP:PID); 内訳: repo/processed/thyr_enrichment_test.rds `units$R_Tumor$pathway` の接頭辞集計(2026-08-12 コンテナ内で照合) | Methods / Q-06 | verified |

### G. spike-in 対照

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-30 | HALLMARK_ADIPOGENESIS 195 遺伝子 × 1.15 倍、B_Tumor High 腕 9 検体 | spike-in の設計 | run/xeon_results/logs/spikein.log:1 キー「Spiked HALLMARK_ADIPOGENESIS (195 genes present) by 1.15x」 | Methods | verified |
| N-31 | NES 2.28、p 0.0002、q_bh 0.0101(H 50 セット中 rank 1)、q<0.10 で回収 TRUE | spike-in の回収 | 同 log:2 | Results / Supp | verified |
| N-32 | 0 | spike 以外の Hallmark セットで q<0.10(偽陽性の余剰なし) | 同 log:3 キー「besides the spike at q<0.10: 0」 | Results / Supp | verified |

### H. signature agreement(R系×B系の順位相関)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-33 | rho +0.4589(共有 15560 遺伝子)、shuffle 帰無 95% 帯 [−0.3914, +0.3930]、両側 p 0.0197 | tumor(R_Tumor × B_Tumor)の一致 | run/xeon_results/logs/signature_agreement.log:3-4 キー「tumor R_Tumor x B_Tumor」 | Results | verified |
| N-34 | rho +0.3756(共有 15459 遺伝子)、95% 帯 [−0.4615, +0.4580]、両側 p 0.1199 | normal(R_Normal × B_Normal)の一致 | 同 log:1-2 キー「normal R_Normal x B_Normal」 | Results | verified |

### I. REO パネル(510 → 520 → 530)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-35 | 候補プール up 317 / down 182、57694 ペア評価、全基準通過 153 ペア | 510 候補選定 | run/xeon_results/logs/510_select_reo_pairs.log:4-6 | Methods | verified |
| N-36 | median_diff [1.159, 4.700]、reversal_rate [0.53, 0.87] | 通過 153 ペアの範囲 | 同 log:7 | Methods / Supp | verified |
| N-37 | パネル 10 ペア、境界 = R0 基準で score > 2 を positive | 520 パネル確定 | run/xeon_results/logs/520_finalize_reo_panel.log:7 「Selected panel pairs: 10」・log:9 「Boundary (R0-based)」 | Methods / Results | verified |
| N-38 | 訓練分類: R0 12/12 negative、R1 13/15 positive(2 negative); score 範囲 R0 [0,2] / R1 [0,10] | 520 訓練成績 | 同 log:3-6(分類表)・log:8(score 範囲) | Results | verified |
| N-39 | P1 DBH/PROM1、P2 ZNF560/CSTA、P3 CA4/CTHRC1、P4 ADRA2B/FCGR2B、P5 DNASE1L2/LOX、P6 PNMA8B/FCER1A、P7 CASKIN1/CD1C、P8 GPR62/IL13RA2、P9 NPAS1/PBK、P10 ICAM4/S100A8(up/down、median_diff 降順) | パネル 10 ペアの構成 | run/xeon_results/output/reo_panel.csv 全10行(二機バイト一致; Ensembl ID・median_diff・reversal_rate も同 CSV) | Tab(パネル構成) | verified |
| N-40 | 片側 BM p = 0.1127(mc)、効果 P(Low<Mid) = 0.616 | 530 Read B(out-of-sample: Mid > Low 逆転スコア) | run/xeon_results/logs/530_evaluate_reo_panel.log:3 | 主結果 / Abstract | verified |
| N-41 | スコア中央値: R_Sporadic(train) 0 / R_Low 1 / R_Mid 4 / R_High(train) 6 | 帯別逆転スコア | 同 log:9-13(median 列) | Results / Fig | verified |
| N-42 | R_Low 9 negative / 8 positive、R_Mid 8 negative / 11 positive | 帯別分類(境界 score>2) | 同 log:14-17 | Results | verified |

### J. REO 診断3種

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-43 | 外れ値 R_Low 0 / R_Mid 0; 除去前後とも Low 中央値 1.0 (n=17) / Mid 4.0 (n=19)、BM p=0.1127、P(Low<Mid)=0.616 で不変 | (a) 外れ値診断 | run/xeon_results/logs/reo_lowmid_outliers.log:2-3(0 outlier)・log:9-11(Gradient robustness) | Supp | verified |
| N-44 | 共通尺度純度の中央値: R_Sporadic 0.690 (n=15) / R_Low 0.704 (n=15) / R_Mid 0.739 (n=16) / R_High 0.814 (n=16)(ペア付き RET コホート) | (b) 純度診断: 帯別純度 | run/xeon_results/logs/reo_lowmid_purity.log:10(コホート n)・log:14-18(median) | Supp | verified |
| N-45 | Spearman(purity, score): R_Low +0.684 (n=15) / R_Mid +0.345 (n=16) / pooled +0.538 (n=31) | (b) 純度×スコアの相関 | 同 log:20-23 | Supp / Disc | verified |
| N-46 | 偏相関 band-score \| purity = +0.146、置換 p(片側) = 0.2162 | (c) 交絡診断の主結果 | run/xeon_results/logs/reo_lowmid_confound.log:5 キー「PARTIAL band-score \| purity」 | Results / Supp | verified |
| N-47 | lo_purity: p=0.3069、P(Low<Mid)=0.580(Low n=8/Mid n=7); hi_purity: p=0.4197、P=0.532(n=7/9) | (c) 純度層別の Mid>Low | 同 log:9-10 | Supp | verified |
| N-48 | band-score +0.142 / band-purity +0.036 / score-purity +0.538(n = Low 15 / Mid 16) | (c) 順位相関(素) | 同 log:1・4 | Supp | verified |

### K. 図付随数値・検査

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-49 | \|M\| 中央値: R_Tumor 0.183 / R_Normal 0.108 / B_Tumor 0.105 / B_Normal 0.121 | MA プロットの効果量サマリ(q<0.10 数は N-16 と一致) | run/i9canon/fig_ma_gene_bm.log:2-5(正準図 = repo/output/figures/fig_ma_gene_bm.png 2026-08-10 01:40 版; 図は B.14 で二機 md5 一致) | Fig(MA)caption | verified |
| N-50 | AS 中央値: R_Sporadic 0.0 / R_Low 15.6 / R_Mid 55.5 / R_High 86.7(スコア中央値 0/1/4/6 は N-41 と一致) | fig_reo_grading の帯別 assigned share | run/xeon_results/logs/fig_reo_grading.log:2-5 | Fig(REO)caption | verified |
| N-51 | [FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 415](両機) | テストスイート結果 | run/xeon_results/logs/tests.log:48; run/i9canon/tests.log:28 | Methods(再現性) | verified |
| N-52 | raw 1819 ファイル md5 全一致(唯一の差分は Xeon 側のみの logs/*.parcel 1件で発現データ本体ではない) | 二機の入力同一性 | run/xeon_provenance/i9_raw_md5.txt(1819行) vs xeon_raw_md5.txt(1820行)、差分は raw_diff.txt:1-2 | Methods(再現性) | verified |

### L. 外部アンカー照合(追補計算 2026-08-12、claim_map C-13)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
|---|---|---|---|---|---|
| N-53 | 20 セル中 19 で k=0(層A 5 リスト × 4 unit。全 33 entries が現行注釈に 1:1 解決、全遺伝子が全 unit の検定対象に存在) | 外部放射線関連遺伝子リスト(Abend 2013 正常8/腫瘍6、Abend 2012 ペア差11、Dom 2012 正常7、CLIP2)と DEG 集合(q<0.10)の員数照合 | diagnostics/output/external_gene_anchors.log:2「All symbols resolved」・:5-24(Membership summary); rds は diagnostics/output/external_gene_anchors.rds `$summary` | Disc / Supp | verified |
| N-54 | S100A10: R_Tumor で effect 0.167(High 群で低発現)、p_exact 0.00208、q_storey 0.079、rank 233 | 唯一の非ゼロセル(dom2012_normal × R_Tumor = **組織対応外セル**: 正常組織由来リストを腫瘍ユニットで照合したもの)。Dom 2012 では被曝側正常組織で上方制御であり、組織側・方向のいずれも原報告と一致しない | 同 log:17「dom2012_normal R_Tumor [cross-tissue] k=1」・:102(per-gene detail); rds `$detail` | Disc / Supp | verified |

### M. D6 較正の派生値・採用時測定(執筆用、2026-08-12 追加)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
|---|---|---|---|---|---|
| N-56 | 0.064(102/1,600) | 本番 D6 の pooled global-null 超過率: 16セル合算で「q<0.10 の発見が1つ以上出た」帰無レプリケートの割合(名目上限 0.10)。13/16 セルが 0.10 以下 | N-24 の表(run/xeon_final_20260811/logs/d6_calibration.log:7-22 の n_any_discovery 列)の算術和 12+1+1+9+7+3+2+6+8+1+4+10+18+3+5+12=102、102/(16×100)=0.06375 | Methods / Results / Q-10 | verified |
| N-57 | per-set p + BH: pooled P(≥1) = 0.045(参考: WY-FWER 0.112) | **採用時測定**(2026-08-08、実データ走行前・開発 B=999): 採用された推論の held-out 較正値 | repo/diagnostics/output/gsea_null_calibration_alternatives_20260808.log:40 キー「Pooled P(>=1): BH 0.045」 | Methods(選定経緯) | verified |
| N-58 | pooled tail-ratio 0.140 / 再標準化 0.221(最悪セル 0.44 = B_Normal/radiation) | **棄却測定**(2026-08-08、実データ走行前): 旧 D2 の tail-ratio 系 FDR の較正破綻(名目 0.10) | repo/diagnostics/output/gsea_null_calibration_restd_20260808.log:40 キー「plain 0.140 ; restandardized 0.221」・:22(0.44) | Methods(選定経緯)/ Q-10, Q-12 | verified |

### N. R_Tumor DEG の ORA 注釈(claim_map C-14 — 水準は仮説生成。初回実行 2026-08-12 は追補計算 = diagnostics/、2026-08-13 に scripts/430 へ編入し正準再実行 — 旧版と table/config identical・log 1-466 行バイト同一を確認済み)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
|---|---|---|---|---|---|
| N-59 | up: 0/0/0/0、down: 12/105/205/7、combined: 6/37/71/2(順に H / C2:CP / C5:GO:BP / C2:CGP:radiation。フラグ = family×list 内 BH q_bh<0.10) | ORA フラグ数の全マトリクス。up リストはフラグ皆無、down リスト(High 群で低発現側)に集中。family 別セット数は 50/2262/3856/24 で 420/D6 と一致(整合性検査) | repo/output/430_annotate_deg_ora.log:8-19(Summary ブロック); 全表は repo/processed/thyr_deg_ora_annotation.rds `$table` | Supp.Tab.2 / Disc(仮説生成) | verified |
| N-60 | radiation down 7 セット: ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_24HR k=46/126(期待 6.4、p 1.09e-27)・同 6HR k=30/84(p 3.53e-18)・RASHI_RESPONSE_TO_IONIZING_RADIATION_6 q 0.0081・SMIRNOV_RESPONSE_TO_IR_6HR_DN q 0.0081・GHANDHI_BYSTANDER_IRRADIATION_UP q 0.0196・WARTERS_RESPONSE_TO_IR_SKIN q 0.0497・MACAEVA_PBMC_RESPONSE_TO_IR q 0.0872(combined では ZHOU 2 セットのみ) | radiation ファミリーのフラグ詳細(全て down = High 群で低発現側) | 同 log:22-30 | Supp.Tab.2 / Disc(仮説生成) | verified |
| N-61 | H down 12 セットの上位: E2F_TARGETS k=46/199(p 1.89e-18)・G2M_CHECKPOINT k=41/198(p 8.26e-15)・EMT・MYC_TARGETS_V1・MITOTIC_SPINDLE・KRAS_SIGNALING_UP ほか; H combined 6 セット(E2F・G2M・MITOTIC_SPINDLE・DNA_REPAIR・SPERMATOGENESIS・MYC_V1) | Hallmark のフラグ詳細(増殖・細胞周期プログラムが down 側に集中) | 同 log:449-466 | Supp.Tab.2 / Disc(仮説生成) | verified |
| N-62 | C2:CP combined 上位例: REACTOME_RESOLUTION_OF_D_LOOP_STRUCTURES k=17/35・CELL_CYCLE_CHECKPOINTS k=63/286・DNA_REPAIR k=68/322・HOMOLOGY_DIRECTED_REPAIR k=37/135・DNA_DOUBLE_STRAND_BREAK_REPAIR k=42/165(いずれも q≤0.0001) | canonical pathways のフラグ上位(DNA 修復・HR・チェックポイント系のテーマ) | 同 log:31-45(down 側含む全 105/37 セットは rds `$table`) | Supp.Tab.2 / Disc(仮説生成) | verified |

### O. 腕間年齢差の推定(追補計算 2026-08-13、claim_map C-15 — 開示であり交絡検定ではない)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
|---|---|---|---|---|---|
| N-63 | HL 中央値差(High−Sporadic、年)+ BM 効果 P(Sporadic<High)(410 と同一推定量 `.bm_effect`、検定は不使用・p 値なし)、腕内復元抽出 percentile ブートストラップ B=9999・seed 19450809・95% CI(2.5/97.5)。被曝時年齢は Sporadic 全 NA のため腕間推定なし(NA 員数のみ記録) | 腕間年齢差推定の設定一式 | diagnostics/output/age_arm_difference.log:1(設定行)・:11-15(AGE_EXPOSURE NA 員数); rds は同 output/age_arm_difference.rds `$config` | Methods / Tab.2 脚注 | verified |
| N-64 | R 系(Sporadic n=12 vs High n=15): HL +2.5 年 [−1.0, 6.0]、P(Sporadic<High) 0.625 [0.400, 0.828](CI はゼロ差/0.5 を跨ぐ) | 手術時年齢の腕間差(R 系)。腕別中央値[範囲]は N-12 と一致(整合性検査) | 同 log:3-5; rds `$summary` R 行(hl 2.5 [−1, 6]、effect 0.6250 [0.4000, 0.8278]) | Results(コホート記述)/ Tab.2 脚注 | verified |
| N-65 | B 系(Sporadic n=27 vs High n=9): HL +8.0 年 [3.0, 12.0]、P(Sporadic<High) 0.850 [0.681, 0.973](CI はゼロ差/0.5 を跨がない) | 手術時年齢の腕間差(B 系)。腕別中央値[範囲]は N-13 と一致(整合性検査) | 同 log:7-9; rds `$summary` B 行(hl 8 [3, 12]、effect 0.8498 [0.6811, 0.9733]) | Results(コホート記述)/ Tab.2 脚注 | verified |

### P. 設定値・規則の台帳化(Methods 起草用 2026-08-13 — 出典はコミット済みコード、計算なし)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
|---|---|---|---|---|---|
| N-66 | GENCODE v36(GDC 参照 GTF)の exon-union 遺伝子長; 発現は STAR - Counts(open access) | 参照アノテーションとデータ種別 | scripts/020 ヘッダ(GENCODE v36)・scripts/010 ヘッダ(STAR - Counts) | Methods | verified |
| N-67 | stranded 列合計の小/大比 ≤ 0.5 なら stranded(大きい列を採用)、それ以外は unstranded 列 | ライブラリ strand 判定規則(検体別) | scripts/120 ヘッダ(Dobin's general rule, ratio form) | Methods | verified |
| N-68 | NIOSH-IREP 甲状腺モデル「AS associated with the expected value of ERR」; 電子 E>15keV、線量 cSv=mGy/10、被曝年 1986、出生年=1986−被曝時年齢、手術年=出生年+手術時年齢。値は研究者計算の正準 CSV | IREP 入力規約 | scripts/130:1-15 ヘッダ; 来歴・入力同一性監査は計画v2 B.11 | Methods | verified |
| N-69 | dose 0 → Sporadic; 0<AS≤33.3 Low; 33.3<AS<66.6 Mid; AS≥66.6 High(境界例なし、2026-07-28 検証) | AS 帯規則 | config.R:36-43(AS_LOW_MAX 33.3・AS_HIGH_MIN 66.6) | Methods | verified |
| N-70 | 0.6 | pooled 共通尺度の相対純度閾値(main BM 採用条件) | config.R:45-46(PURITY_THRESHOLD) | Methods | verified |
| N-71 | iDEGES 3 反復; スクリーン = 置換 BM + BH q<0.10(DEGES_FDR); スケーリングは MUREN; 前処理 protein_coding → filterByExpr | 310 正規化の設定 | scripts/310:37(ITERATION 3L)・config.R:55-57(DEGES_FDR 0.10)・310 ヘッダ | Methods | verified |
| N-72 | ランキング = 符号付き BM 統計量の tie-averaged normal scores; ES = gseaParam=1 の block 評価(tie-free 入力で標準 GSEA と一致); 推論 = per-set 符号条件付き置換 p + family 内 BH、q_bh<0.10; size 窓 15–500 | 420 セットレベル推論の設定 | scripts/420 ヘッダ; lib/gsea_collections.R:27-28(15L/500L); config.R:51-53(FDR_CUT 0.10) | Methods | verified |
| N-73 | unit 内独立ラベルシャッフル 9,999 回、unit 別 seed(基底 19450809)で shuffle プロファイル同士を相関(410 の perm_index は意図的に不使用 — 等 n unit で同一 index となるため) | 署名一致の帰無参照帯の設定 | diagnostics/signature_agreement.R:57(AGREEMENT_SEED_BASE)・:84,97(N_PERM 使用)・ヘッダ | Methods | verified |
| N-74 | 候補プール = \|effect−0.5\| 上位 500; dead zone \|r\|<log2(1.2); R0 = 非 dead-zone 検体で符号一致(例外≤1)かつ q10(\|r\|)≥log2(1.5); R1 = 逆転 >50% かつ <100% | REO 候補ペア選定規則 | scripts/510:44(N_CANDIDATES 500L)・:48(log2(1.5))・:100,130,136(規則実装); config.R:48-49(DEAD_ZONE log2(1.2)) | Methods | verified |
| N-75 | 貪欲選定: 遺伝子再使用禁止 + 既採用ペアとの Spearman <0.75、目標 10 ペア; 530 の Mid>Low 片側 BM は mc・seed 19860426(正準シード) | REO パネル確定と評価検定の設定 | scripts/520:30-31(TARGET_PANEL_SIZE 10L・CORRELATION_THRESHOLD 0.75); scripts/530:78(seed = SEED); config.R:5-8 | Methods | verified |
| N-76 | 群内純度順位保存 0.93–0.99(腕別推定との比較)、High vs Sporadic 腫瘍プロファイル相関 0.99 | 220 の腕プーリング妥当性の**来歴実測値**(設計時測定。一次ログではない — 証拠階層は設計選択の来歴)。二重統計に当たるためライセンスには使わない(研究者決定 2026-08-13)— プーリングの記述は理論構造(相対尺度の共通化・driver 支配前提・軸の別・役割の限定)で立てる | `git show 8eed384:scripts/220_estimate_tumor_purity.R` ヘッダ行8-10(run 時点の凍結版。現行ヘッダは 2026-08-13 に理論ライセンスへ書き換え、実測値の記録はこの行と run コミットに保存) | **不使用**(査読応答の受けとして保存) | verified |

## 図表台帳(図・表は1行ずつ — キャプションも検査対象)

| 列 | 内容 |
| --- | --- |
| F-ID | Fig.1, Tab.1, Supp.Fig.1, … |
| 生成 | 生成スクリプト(figures/xxx.R 等)と入力 rds |
| キャプション内の主張/数値 | 対応する C-ID / N-ID(キャプションも本文と同じ検査を受ける) |
| 状態 | draft / verified |

番号(Fig.1 等)は仮置き。図表構成の確定時に振り直す。C-ID は主張マップ記入後に対応付ける。

| F-ID | 生成 | キャプション内の主張/数値 | 状態 |
| --- | --- | --- | --- |
| Fig.1(仮)遺伝子別 BM 証拠 | figures/fig_gene_bm_evidence.R ← processed/thyr_expression_test.rds + thyr_se_raw.rds。正準 repo/output/figures/fig_gene_bm_evidence.png(二機バイト一致) | C-未、N-16・N-18(付随ログ fig_gene_bm_evidence.log:2-5) | draft |
| Fig.2(仮)MA プロット | figures/fig_ma_gene_bm.R ← processed/thyr_normalized_counts.rds + thyr_expression_test.rds + thyr_se_raw.rds。正準 repo/output/figures/fig_ma_gene_bm.png(2026-08-10 01:40 版。B.14 で二機 md5 一致) | C-未、N-16・N-49 | draft |
| Fig.3(仮)REO グレーディング | figures/fig_reo_grading.R ← processed/thyr_reo_panel.rds + thyr_reo_evaluation.rds + thyr_se_raw.rds + thyr_case_assigned_share.rds。正準 repo/output/figures/fig_reo_grading.png(二機バイト一致) | C-未、N-41・N-50 | draft |
| フロー図(仮)コホートフロー(6段、**両 driver 層並記** — 2026-08-14 研究者承諾: 層別途中経過は本文でなくこの図が担う) | figures/fig_cohort_flow.R ← thyr_cohort_flow.rds(230 出力)。数値は N-08(合算+RET/BRAF 別)。**未実行 — 実行 Go 待ち** | C-10、N-08 | draft |
| Tab.1(仮)臨床コホート構成(Driver × AS 帯 × pair) | tables/tab_cohort_composition.R ← thyr_clinical.rds + thyr_case_assigned_share.rds + SE | C-未、N-11 | draft |
| Tab.2(仮)解析コホート症例特性(群×帯: n・性別・年齢・サブタイプ) | tables/tab_case_characteristics.R ← thyr_analysis_cohorts.rds + thyr_clinical.rds(§15 B と同経路。**未実行 — 実行 Go 待ち**。初回実行で N-09/N-12/N-13 と照合) | C-未、N-12・N-13 | draft |
| Tab.3(仮)410 主結果(ユニット × 検定遺伝子数・DEG・HC p) | tables/tab_gene_level_summary.R ← thyr_expression_test.rds(**未実行**。初回実行で N-15〜N-20 と照合) | C-未、N-15(遺伝子数)・N-16・N-19・N-20 | draft |
| Tab.4(仮)420/D6 サマリ(ユニット × コレクション) | tables/tab_set_level_summary.R ← thyr_enrichment_test.rds + diagnostics/output/gsea_null_calibration.rds(**未実行**。較正 rds のセル表は防御的に探索 — 初回実行で N-24/N-27/N-28 と照合) | C-未、N-24〜N-28 | draft |
| Supp.Tab.1(仮)D6 較正 16 セル表(m_sets・p_any・95%CI・mean/max 発見数) | diagnostics/gsea_null_calibration.R ← 310 正規化。正準値 = run/xeon_final_20260811/logs/d6_calibration.log(rds は diagnostics/output/gsea_null_calibration.rds) | C-06、N-24(注記に N-25, N-26) | draft |
| Supp.Data.1(仮)420 全セット結果の完全開示(全 unit × family の pathway・size・ES・NES・p・q_bh 全量) | tables/supp_data_420_full.R ← thyr_enrichment_test.rds(整形のみ・計算なし。**未実行**)。Q-08/Q-10 予備線「可能性の棚は開示で遺す」の実体 | C-05 | draft |
| Supp.Tab.2(仮)R_Tumor DEG の ORA 注釈(family × list の全結果) | scripts/430_annotate_deg_ora.R ← thyr_expression_test.rds(430 編入 2026-08-13。正準 = processed/thyr_deg_ora_annotation.rds + output/430_annotate_deg_ora.log。旧 diagnostics 版と同値確認済み) | C-14、N-59〜N-62 | draft |
| Supp.Tab.3(仮)腕間 concordance(normal / tumor 両ペア: rho・シャッフル帯・p の2行+識別不能注記) | tables/supp_tab_concordance.R ← diagnostics/output/signature_agreement.rds(**未実行** — rds のフィールド名は初回実行で確認・調整)。正準値 = run/xeon_results/logs/signature_agreement.log | C-07、C-17、N-33・N-34 | draft |
| Supp.Tab.4(仮)使用パッケージ版表(パッケージ名・版の一覧 — 2026-08-14 引用整理に伴い追加) | tables/supp_tab_package_versions.R ← docker/versions.tsv(**未実行**。実走時の実効版 = run/xeon_provenance/session_info.txt は組版時に併記)。機械生成のみ・計算なし | Methods(Software 節が参照) | draft |

(番号・採否は図表構成の確定後に更新)

## 改訂メモ

- 2026-08-12: 様式 v1 起草(Claude Code、批准待ち)
- 2026-08-12: 図表台帳を追加 — 図表の来歴とキャプションを一級の検査対象に昇格(批准待ち)
- 2026-08-12: 研究者批准(様式v1確定)
- 2026-08-12: N-01〜N-52 を記入(素材 = run/summary_20260809.md、出典は一次ファイルに付け替え)。
  全行を一次出典と直接照合し verified。rds 由来(N-05, N-12, N-13, N-16〜N-23)は
  rebc-r453:refblas コンテナ内で §15 コマンド相当を再実行して照合。410 rds は repo と
  xeon_results で md5 一致を再確認。D6/420 の出典は Xeon 産正準(xeon_final_20260811/、B.14)を指す。
  図表台帳は既存図3枚+表4点(仮)を下書き(番号仮置き・C-ID 未対応付け)。
- 2026-08-12: 追補計算の規則を起草(Claude Code、批准待ち)。新規数値も既存と同じ規律、
  記録先は paper/、計画v2 への amendment は終了。
- 2026-08-12: 追補計算の初適用(N-53, N-54): 外部アンカー照合(C-13 で読みを実行前に固定)。
  diagnostics/external_gene_anchors.{csv,R} → 同 output/ の log+rds(乱数なし・検定なし・
  正準イメージ内で実行、exit 0)。表示のみの修正2件(k>0 セルの詳細の自己完結化、
  セルラベルの用語分離 [pair-diff]/[tissue-matched]/[cross-tissue] — 「間接」の
  二重使用を解消)を挟んで再実行、算術は不変。
- 2026-08-12: 追補計算の第2適用(N-59〜N-62): R_Tumor DEG の ORA 注釈(C-14 で水準・読みを
  実行前に固定)。diagnostics/deg_ora_annotation.R → 同 output/ の log+rds(乱数なし・
  正準イメージ内・exit 0)。family 別セット数 50/2262/3856/24 は 420/D6 と一致(整合性検査)。
- 2026-08-13: 追補計算の第3適用(N-63〜N-65): 腕間年齢差の CI つき推定(C-15 で読みを
  実行前に固定。研究者 Go 2026-08-13)。diagnostics/age_arm_difference.R → 同 output/ の
  log+rds(seed 19450809・B=9999・正準イメージ内・exit 0)。腕別 n・中央値[範囲]は
  N-09/N-12/N-13 と一致(整合性検査)。照合は log↔rds `$summary` の突き合わせで verified。
- 2026-08-13: セクション P(N-66〜N-76)を Methods 起草に伴い追加 — 本文が引く設定値・
  規則・来歴実測(N-76)の台帳化。追補計算ではない(計算なし、出典は全てコミット済み
  コードの行参照)。照合はコード直読で verified(Claude Code、2026-08-13)。N-76 のみ
  一次ログでなくヘッダ記載の来歴値である旨を行内に明記。
- 2026-08-13: N-76 を不使用に変更(研究者決定): Methods は純度推定の手続きの事実のみ
  記述し(共通相対尺度・1 コホート 1 run・閾値 0.6 = N-70)、プーリングの設計弁明と
  検証実測は本文に載せない。来歴は 220 ヘッダに保存、行は査読応答の受けとして残す。
- 2026-08-13(同日改訂、研究者決定): プーリングの記述は script ヘッダ・Methods とも
  **理論構造で立てる**に変更 — (1) ContamDE 純度は同時推定集合内でのみ相対比較可能で、
  腕別 run は尺度を分断し単一閾値の意味を壊す(腕別で先に回して観察された事実)、
  (2) 共通参照仮定は driver 層別設計自体の前提に乗る、(3) 純度軸は曝露対比と別軸、
  (4) 役割はコホート内相対フィルタ+診断共変量のみ。実測値(0.93–0.99/0.99)は
  二重統計のためライセンスに使わず、run コミットの凍結ヘッダと N-76 行にのみ保存。
  現行 220 ヘッダは同日この線で書き換え(コメントのみ、挙動不変。run 状態は N-01 の
  ハッシュで凍結済み)。N-76 出典を `git show 8eed384:` 参照に付け替え。
- 2026-08-13: ORA の 430 化(研究者指示): diagnostics/deg_ora_annotation.R を
  scripts/430_annotate_deg_ora.R へ git mv(履歴保持)し、ヘッダを新根拠構造
  (Thyroid 来歴・G&B sampling-model 軸・Q-15 (3b) 転落基準)で改稿、出力を
  processed/thyr_deg_ora_annotation.rds へ変更して正準イメージで再実行(exit 0)。
  同値検査: 旧 diagnostics 版と table・config が identical、log は 1-466 行バイト同一
  (相違は最終行の保存先パスのみ)→ N-59〜N-62 の値・行番号引用は不変のまま出典パスのみ
  付け替え。旧成果物(diagnostics/output/)は追補計算の記録として残置。
  セクション N の表題に来歴を追記
- 2026-08-14: 3バケツ整理の実装(研究者批准): 本台帳を足場(検査機構)と位置づけ、
  採録確定まで併走→削除の運用を明記。非掲載数値の出典自蔵は手元台帳側の責務
- 2026-08-14: 図表整形スクリプト7本を作成(全て整形のみ・計算なし・**未実行**):
  figures/fig_cohort_flow.R、tables/tab_case_characteristics.R・tab_gene_level_summary.R・
  tab_set_level_summary.R・supp_data_420_full.R・supp_tab_concordance.R・
  supp_tab_package_versions.R。図表台帳の該当行を実体パスへ付け替え。初回実行時に
  各 N 行との照合を行う(スクリプト内に照合対象を明記)。実行は研究者 Go 待ち
