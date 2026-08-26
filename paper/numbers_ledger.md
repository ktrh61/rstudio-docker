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
| N-07 | 410 の置換インデックスを再利用(9999/対比) | 420 の帰無生成方式 | run/xeon_final_20260811/logs/420_test_gene_sets.log:1 キー「Reusing 410's permutation index」(対比別ハッシュは N-05 の rds `$config$perm_index_hash`) | Methods | verified |

### B. コホート

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-08 | 440 → 248 → 77 → 70 → 69 → 63(n_RET 73→73→34→31→31→27、n_BRAF 175→175→43→39→38→36) | コホートフロー6段(all_cases → driver_classified → band_sporadic_or_high → paired → pcod_clean → purity_pass) | run/xeon_results/logs/230_finalize_analysis_cohorts.log:1-8 キー「Cohort flow (main BM)」 | Results / フロー図表(Methods は n = 63 のみ) | verified |
| N-09 | B_High 9 / B_Sporadic 27 / R_High 15 / R_Sporadic 12 | main BM 63症例の群内訳 | 同 log:9-12 キー「include_main_bm by group」 | Results | verified |
| N-10 | training 27 / evaluation 36(R_Low 17・R_Mid 19) | REO コホート | 同 log:13 キー「REO training: 27」 | Methods / Results | verified |
| N-11 | 臨床 440 症例; AS 帯 non_exposed 81 / no_reference 146 / (0,33.3) 100 / [33.3,66.6) 82 / [66.6,100] 31; paired 392 / unpaired 48 | 臨床全数と AS バンド構成(整合性検査通過) | run/xeon_results/logs/tab_cohort_composition.log:1-4(Driver×AS帯×pair の全225行は同 log:9-233) | Tab(コホート構成) | verified |
| N-12 | R系 4群(R_Sporadic/R_Low/R_Mid/R_High): n 12/17/19/15、女/男 10/2・13/4・14/5・11/4、手術時年齢中央値[範囲] 20.5[14–27]/30[22–44]/25[17–31]/23[14–31]、被曝時年齢 NA/12[6–19]/3[0–13]/2[0–12]、CCDC6-RET 6/8/12/7、NCOA4-RET 2/5/3/4、RET-OTHER 4/4/4/4 | 症例特性(R系、群×帯) | run/xeon_results/processed/thyr_analysis_cohorts.rds + thyr_clinical.rds を case_submitter_id=REBC_ID で結合(§15 B で再現。2026-08-12 コンテナ内で再実行し一致) | Methods / Tab(症例特性) | verified |
| N-13 | B系 2群(B_Sporadic/B_High): n 27/9、女/男 23/4・4/5、手術時年齢中央値[範囲] 24[11–29]/29[19–39]、被曝時年齢 NA/3[0–13]、Driver は全例 BRAF.MutV600E | 症例特性(B系) | 同上(§15 B) | Methods / Tab(症例特性) | verified |

### C. 310 正規化

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-14 | 58448 genes × 906 samples | SE(発現行列)の寸法 | run/xeon_results/logs/310_normalize_counts.log:1 キー「SE: 58448 genes x 906 samples」 | Methods | verified |
| N-15 | 対比別: filterByExpr 後 15621/15466/16002/15863、pi0_hat(iter3) 0.593/0.956/0.943/0.727、収束 Jaccard(iter3) 0.962/0.967/0.993/0.987、norm.factors 範囲 [0.9058,1.1095]/[0.9106,1.0665]/[0.8772,1.0947]/[0.8988,1.0873](順に R_Tumor/R_Normal/B_Tumor/B_Normal) | 310 正規化サマリ(4対比) | 同 log:9-32(各「Unit <名>」ブロックの iteration 3 行と norm.factors 行) | Methods / Supp | verified |

### D. 410 遺伝子別検定(主結果)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-16 | 1765 / 0 / 0 / 1(R_Tumor/R_Normal/B_Tumor/B_Normal) | DEG 数(Storey q<0.10) | repo/processed/thyr_expression_test.rds `units[[u]]$genes$q_storey < 0.10`(§15 A); ログ照合 run/xeon_results/logs/fig_gene_bm_evidence.log:2-5 | 主結果 / Fig / Abstract | verified |
| N-17 | R_Tumor: up 971 / down 794; B_Normal: up 1 / down 0(up = effect>0.5 = High 群で高発現) | DEG の方向内訳 | 同 rds `genes$effect`(§15 A2) | Results | verified |
| N-18 | 4.03e-06 / 4.08e-05 / 1.68e-05 / 1.13e-06 | 対比別 min p_exact | 同 rds `min(genes$p_exact)`(§15 A); ログ照合 fig_gene_bm_evidence.log:2-5 | Results / Supp | verified |
| N-19 | 0.593 / 0.955 / 0.943 / 0.727 | Storey pi0(410、plug-in λ=0.5) | 同 rds `units[[u]]$pi0$estimate`(§15 A) | Results / Supp | verified |
| N-20 | 0.0112 / 0.3199 / 0.1815 / 0.0773 | 主オムニバス HC p(α0=0.1) | 同 rds `units[[u]]$omnibus` の hc 行(§15 A) | 主結果 / Abstract | verified |
| N-21 | オムニバス全表: 4対比 × {count α=1e-2/1e-3/1e-4, max, hc} の observed / 帰無中央値 / p(R_Tumor count@1e-2: 1026/82/0.0128 … B_Normal hc: 36.31/−0.41/0.0773。全20行は出典参照) | オムニバス検定の完全な表 | 同 rds `units[[u]]$omnibus`(§15 A の print 出力) | Supp | verified |
| N-22 | ENSG00000198908.12(BHLHB9)、effect 0.967、q_storey 0.0130、rank 1 | B_Normal 唯一の DEG | 同 rds B_Normal `genes`(q_storey<0.10 の1行。effect 0.9670782、q 0.01298213); symbol は raw の STAR counts tsv 行「ENSG00000198908.12→BHLHB9(protein_coding)」 | Results | verified |
| N-23 | 0.0125 | B_Normal オムニバスの最小 p(max 統計量; observed 17.70、帰無中央値 6.99) | 同 rds B_Normal `$omnibus` の max 行(§15 A) | Results | verified |

### E. D6 帰無較正(GSEA 較正)

正準は Xeon 産(2026-08-11、B.14)。i9 の 8/10 再計算版(run/i9canon/d6_calibration.log)と中身一致。

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-24 | 16セル全表(4対比×4コレクション: m_sets、p_any、95%CI、mean/max 発見数。p_any 範囲 0.01–0.18) | D6 較正の完全な表 | run/xeon_final_20260811/logs/d6_calibration.log:6-39(rds は repo/diagnostics/output/gsea_null_calibration.rds) | Supp | verified |
| N-25 | 0.18(95%CI 0.110–0.269 — rds 実値 ci_hi 0.26947709。旧記載 0.270 は転記時の丸め誤り、2026-08-25 訂正) | 唯一の名目 0.10 超過: B_Normal/H の p_any(CI 下限も 0.10 超) | 同 log:19 キー「13 B_Normal H 50 100 18 0.18 0.110311」・log:36(ci_hi 0.2695) | Results / Disc | verified |
| N-26 | B_Tumor/radiation 0.10(CI 0.049–0.176)、B_Normal/radiation 0.12(CI 0.064–0.200) | CI が 0.10 を跨ぐ境界セル | 同 log:18・22(ci_hi は log:35・39) | Supp / Disc | verified |

### F. 420 遺伝子セット検定

正準は Xeon 産(2026-08-11、B.14)。i9 の 8/10 再計算版(run/i9canon/420_test_gene_sets.log)と中身一致。

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-27 | 0(16/16 セル) | 420 本試験の発見数(q<0.10、全対比 × collection) | run/xeon_final_20260811/logs/420_test_gene_sets.log:5-8(各行「q<0.10 0」×4) | 主結果 / Abstract | verified |
| N-28 | min q 表: R_Tumor(6192セット) H 0.179/C2:CP 0.928/C5:GO:BP 0.936/rad 0.302; R_Normal(6141) 0.380/1.000/1.000/0.980; B_Tumor(6242) 0.520/0.892/0.262/0.114; B_Normal(6223) 0.284/0.283/0.467/0.137 | 対比×コレクションの min q と総セット数 | 同 log:5-8 | Supp | verified |
| N-29 | H=50, C2:CP=3910, C5:GO:BP=7538, C2:CGP:radiation=28 | フィルタ前のセット数 | 同 log:2-3 キー「Gene sets before filtering」 | Methods | verified |
| N-55 | C2:CP 検定済み内訳(R_Tumor, 2262セット): Reactome 1095 / WikiPathways 633 / KEGG MEDICUS 217 / BioCarta 134 / PID 183(未分類 0)。msigdbr 26.1.0 | C2:CP の構成(=Reactome・KEGG を含む統合であることの数値根拠。他対比はフィルタ差で僅差) | 定義: lib/gsea_collections.R `GSEA_COLLECTION_SPECS`(CP:REACTOME/CP:WIKIPATHWAYS/CP:KEGG_MEDICUS/CP:BIOCARTA/CP:PID); 内訳: repo/processed/thyr_enrichment_test.rds `units$R_Tumor$pathway` の接頭辞集計(2026-08-12 コンテナ内で照合) | Methods / Q-06 | verified |

### G. spike-in 対照

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
| --- | --- | --- | --- | --- | --- |
| N-30 | HALLMARK_ADIPOGENESIS 195 遺伝子 × 1.15 倍、B_Tumor High 群 9 検体 | spike-in の設計 | run/xeon_results/logs/spikein.log:1 キー「Spiked HALLMARK_ADIPOGENESIS (195 genes present) by 1.15x」 | Methods / Supp Methods | verified |
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
| N-35 | \|effect−0.5\| 上位500は up 318 / down 182、全構築標本で finite log2 TPM を要求後は up 317 / down 182、全 cross-direction 57694 ペア中153ペアが全基準通過 | 510 候補選定 | `processed/thyr_expression_test.rds` の上位500直接照合 + scripts/510 の finite-gene screen + `processed/thyr_reo_candidate_pairs.rds`(153行); 旧 run/xeon_results/logs/510_select_reo_pairs.log:4-6 | Results / Supp Results | verified |
| N-36 | median_diff [1.159, 4.700]、reversal_rate [0.53, 0.87] | 通過 153 ペアの範囲 | 同 log:7 | Supp Results | verified |
| N-37 | パネル 10 ペア、境界 = score > 2 を positive(閾値 = 訓練 Sporadic の最大逆転スコア; ログ表記 R0-based)| 520 パネル確定 | run/xeon_results/logs/520_finalize_reo_panel.log:7 「Selected panel pairs: 10」・log:9 「Boundary (R0-based)」 | Methods / Results | verified |
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
| N-50 | AS 中央値: R_Sporadic **未定義**(非被曝 — 2026-08-15 図改修で専用ストリップ表示に変更、旧表記 0.0 は y=0 配置の名残)/ R_Low 15.6 / R_Mid 55.5 / R_High 86.7(スコア中央値 0/1/4/6 は N-41 と一致) | fig_reo_grading の帯別 assigned share | run/xeon_results/logs/fig_reo_grading.log:2-5 | Fig(REO)caption | verified |
| N-51 | [FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 415](両機) | テストスイート結果 | run/xeon_results/logs/tests.log:48; run/i9canon/tests.log:28 | Methods(再現性) | verified |
| N-52 | raw 1819 ファイル md5 全一致(唯一の差分は Xeon 側のみの logs/*.parcel 1件で発現データ本体ではない) | 二機の入力同一性 | run/xeon_provenance/i9_raw_md5.txt(1819行) vs xeon_raw_md5.txt(1820行)、差分は raw_diff.txt:1-2 | Methods(再現性) | verified |

### L. 外部遺伝子リスト照合(claim_map C-13)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
|---|---|---|---|---|---|
| N-53 | 検証済みアンカー20セル中19で k=0(5リスト × 4対比。全33 entries が現行注釈に1:1解決し、全遺伝子が全対比の検定対象) | qRT-PCR検証済みリスト(Abend 2013 正常8/腫瘍6、Abend 2012 ペア差11、Dom 2012 正常7)および CLIP2 と DEG 集合(q<0.10)の員数照合 | `diagnostics/external_gene_anchors.{csv,R}`; `diagnostics/output/external_gene_anchors.rds` の `$summary` | Results / Disc / Supp | verified |
| N-54 | S100A10: R_Tumor で effect 0.167(High 群で低発現)、p_exact 0.00208、q_storey 0.079、rank 233 | 検証済みアンカーで唯一の非ゼロセル(dom2012_normal × R_Tumor = 組織横断)。Dom 2012では被曝側正常組織で上方制御であり、組織側・方向とも原報告と一致しない | `diagnostics/output/external_gene_anchors.rds` の `$summary`, `$detail` | Results / Disc / Supp | verified |
| N-85 | Ory 2026 多変量署名12セル中9で k=0。非ゼロは全て R_Tumor: 組織共有署名 5/39、正常組織署名 3/40(組織横断)、腫瘍署名 3/46(組織対応)。元リストは50/45/64遺伝子。SPEN-AS1・MCTS2・SNORD47はGENCODE v36に対応せず、他の分母差は対比別検定集合への残存差による | Ory 2026 Supplementary Tables S2/S4/S6 と各対比 DEG 集合(q<0.10)の記述的員数照合。濃縮検定ではない | `diagnostics/ory2026_gene_signatures.csv`; `diagnostics/external_gene_anchors.R`; `diagnostics/output/external_gene_anchors.rds` の `$resolution`, `$summary` | Methods / Results / Disc / Supp | verified |
| N-86 | R_Tumor との重なり: 組織共有 = ATP5MF, MRPL52, NTHL1, URM1, USE1; 正常組織 = PXDN, S100A10, TESC; 腫瘍 = P2RY1, PLK2, EHD4 | N-85の遺伝子別内訳。正常組織署名の3件は組織横断、腫瘍署名の3件は組織対応。Ory署名は正常組織対比およびBRAF対比では k=0 | `diagnostics/output/external_gene_anchors.rds` の `$detail` | Results / Disc / Supp | verified |
| N-87 | PC-OD flag(主4群): R_Sporadic/R_High/B_Sporadic は tumor・normal とも 0、B_High tumor 1(RET 層の主解析集合は不変) | 主コホート QC フラグの群別員数(Results 記述) | processed/thyr_case_outliers.rds(集計 = paper/gpt_review/additional_pcod_flag_counts.csv、追加監査 2026-08-19 実行記録つき)| Results 1 | verified(2026-08-21 rds 直読: B_High tumor 1・他は3群×2組織すべて0) |
| N-88 | 主解析の相対純度中央値(driver コホート別の相対尺度・層内でのみ比較可能): R_Sporadic 0.783 (n=12) / R_High 0.822 (n=15) / B_Sporadic 0.836 (n=27) / B_High 0.922 (n=9)(主コホート解析オブジェクト由来 — REO 全帯プール尺度 N-44 とは別尺度で互換しない) | 主コホート4群の純度開示(表2脚注) | processed/thyr_analysis_cohorts.rds(R 系集計 = paper/gpt_review/additional_covariate_summary.csv とも一致)| Table 2 footnote | verified(2026-08-21 rds 直読: 4群全て) |
| N-89 | GDC 入力の固定: open-access STAR-Counts 906 ファイル、照会 2026-07-23 と 2026-08-09 で byte-identical、manifest md5 7defb0c5574453474c67dfac8367a589 | GDC manifest の保全(Supp Methods) | paper/gpt_review/gdc_manifest_rebc_thyr_star_counts.tsv(コミット 095699b。md5 と 906 行は 2026-08-21 に実測照合済み)| Supp Methods | verified |
| N-91 | REO Mid-Low 検定の MC 規模: 完全枚挙 C(36,17)=8,597,496,600 のため 999,999 回モンテカルロ(シード 19860426、plus-one) | REO 検定の実装定数(Supp Methods) | lib/stat_brunnermunzel.R:411(B=999999L 既定)+ scripts/530:75-77(method="auto"・B 未指定 → 既定適用。auto は C(36,17) 枚挙不能で MC 選択)。組合せ数は算術確認済み | Supp Methods | verified(コード水準。実行ログ照合は run 記録に委ねる) |
| N-92 | 取得照会日(2026-07-23・2026-08-09)はいずれも GDC Data Release 45.0(2025-12-04 公開)の期間内。次リリース 46.0 は 2026-08-10 公開で、45.0 と 46.0 の間に中間リリースなし。取得時にリリース番号のログはなく、帰属はリリースノートの日付照合による(計算なしの文書事実) | GDC データ状態の版明示(Supp Methods) | https://docs.gdc.cancer.gov/Data/Release_Notes/Data_Release_Notes/(照合 2026-08-21)。版固定原典 = NCI-GDC/gdc-docs コミット 6d9a957fd50d5c349e5b94201874837de0a75939 の docs/Data/Release_Notes/Data_Release_Notes.md(45.0/46.0/44.0 の Release Date を逐語確認) | Supp Methods | verified(2026-08-21 原文照合) |

### M. D6 較正の派生値・採用時測定(執筆用、2026-08-12 追加)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
|---|---|---|---|---|---|
| N-56 | 0.064(102/1,600) | 本番 D6 の pooled global-null 超過率: 16セル合算で「q<0.10 の発見が1つ以上出た」帰無レプリケートの割合(名目上限 0.10)。13/16 セルが 0.10 以下 | N-24 の表(run/xeon_final_20260811/logs/d6_calibration.log:7-22 の n_any_discovery 列)の算術和 12+1+1+9+7+3+2+6+8+1+4+10+18+3+5+12=102、102/(16×100)=0.06375 | Methods / Results / Q-10 | verified |
| N-57 | per-set p + BH: pooled P(≥1) = 0.045(参考: WY-FWER 0.112) | **採用時測定**(2026-08-08、実データ走行前・開発 B=999): 採用された推論の held-out 較正値 | repo/diagnostics/output/gsea_null_calibration_alternatives_20260808.log:40 キー「Pooled P(>=1): BH 0.045」 | Supp Methods(選定経緯) | verified |
| N-58 | pooled tail-ratio 0.140 / 再標準化 0.221(最悪セル 0.44 = B_Normal/radiation) | **棄却測定**(2026-08-08、実データ走行前): 旧 D2 の tail-ratio 系 FDR の較正破綻(名目 0.10) | repo/diagnostics/output/gsea_null_calibration_restd_20260808.log:40 キー「plain 0.140 ; restandardized 0.221」・:22(0.44) | Methods(選定経緯)/ Q-10, Q-12 | verified |

### N. R_Tumor DEG の ORA 注釈(claim_map C-14 — 水準は仮説生成。初回実行 2026-08-12 は追補計算 = diagnostics/、2026-08-13 に scripts/430 へ編入し正準再実行 — 旧版と table/config identical・log 1-466 行バイト同一を確認済み)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
|---|---|---|---|---|---|
| N-59 | up: 0/0/0/0、down: 12/105/205/7、combined: 6/37/71/2(順に H / C2:CP / C5:GO:BP / C2:CGP:radiation。フラグ = family×list 内 BH q_bh<0.10) | ORA フラグ数の全マトリクス。up リストはフラグ皆無、down リスト(High 群で低発現側)に集中。family 別セット数は 50/2262/3856/24 で 420/D6 と一致(整合性検査) | repo/output/430_annotate_deg_ora.log:8-19(Summary ブロック); 全表は repo/processed/thyr_deg_ora_annotation.rds `$table` | Supp.Tab.2 / Disc(仮説生成) | verified |
| N-60 | radiation down 7 セット: ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_24HR k=46/126(期待 6.4、p 1.09e-27)・同 6HR k=30/84(p 3.53e-18)・RASHI_RESPONSE_TO_IONIZING_RADIATION_6 q 0.0081・SMIRNOV_RESPONSE_TO_IR_6HR_DN q 0.0081・GHANDHI_BYSTANDER_IRRADIATION_UP q 0.0196・WARTERS_RESPONSE_TO_IR_SKIN q 0.0497・MACAEVA_PBMC_RESPONSE_TO_IR q 0.0872(combined では ZHOU 2 セットのみ) | radiation ファミリーのフラグ詳細(全て down = High 群で低発現側) | 同 log:22-30 | Supp.Tab.2 / Disc(仮説生成) | verified |
| N-61 | H down 12 セットの上位: E2F_TARGETS k=46/199(p 1.89e-18)・G2M_CHECKPOINT k=41/198(p 8.26e-15)・EMT・MYC_TARGETS_V1・MITOTIC_SPINDLE・KRAS_SIGNALING_UP ほか; H combined 6 セット(E2F・G2M・MITOTIC_SPINDLE・DNA_REPAIR・SPERMATOGENESIS・MYC_V1) | Hallmark のフラグ詳細(増殖・細胞周期プログラムが down 側に集中) | 同 log:449-466 | Supp.Tab.2 / Disc(仮説生成) | verified |
| N-62 | C2:CP combined 上位例: REACTOME_RESOLUTION_OF_D_LOOP_STRUCTURES k=17/35・CELL_CYCLE_CHECKPOINTS k=63/286・DNA_REPAIR k=68/322・HOMOLOGY_DIRECTED_REPAIR k=37/135・DNA_DOUBLE_STRAND_BREAK_REPAIR k=42/165(いずれも q≤0.0001) | canonical pathways のフラグ上位(DNA 修復・HR・チェックポイント系のテーマ) | 同 log:31-45(down 側含む全 105/37 セットは rds `$table`) | Supp.Tab.2 / Disc(仮説生成) | verified |

### O. 群間年齢差の推定(追補計算 2026-08-13、claim_map C-15 — 開示であり交絡検定ではない)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
|---|---|---|---|---|---|
| N-63 | HL 中央値差(High−Sporadic、年)+ BM 効果 P(Sporadic<High)(410 と同一推定量 `.bm_effect`、検定は不使用・p 値なし)、群内復元抽出 percentile ブートストラップ B=9999・seed 19450809・95% CI(2.5/97.5)。被曝時年齢は Sporadic 全 NA のため群間推定なし(NA 員数のみ記録) | 群間年齢差推定の設定一式 | diagnostics/output/age_arm_difference.log:1(設定行)・:11-15(AGE_EXPOSURE NA 員数); rds は同 output/age_arm_difference.rds `$config` | Methods / Tab.2 脚注 | verified |
| N-64 | R 系(Sporadic n=12 vs High n=15): HL +2.5 年 [−1.0, 6.0]、P(Sporadic<High) 0.625 [0.400, 0.828](CI はゼロ差/0.5 を跨ぐ) | 手術時年齢の群間差(R 系)。群別中央値[範囲]は N-12 と一致(整合性検査) | 同 log:3-5; rds `$summary` R 行(hl 2.5 [−1, 6]、effect 0.6250 [0.4000, 0.8278]) | Results(コホート記述)/ Tab.2 脚注 | verified |
| N-65 | B 系(Sporadic n=27 vs High n=9): HL +8.0 年 [3.0, 12.0]、P(Sporadic<High) 0.850 [0.681, 0.973](CI はゼロ差/0.5 を跨がない) | 手術時年齢の群間差(B 系)。群別中央値[範囲]は N-13 と一致(整合性検査) | 同 log:7-9; rds `$summary` B 行(hl 8 [3, 12]、effect 0.8498 [0.6811, 0.9733]) | Results(コホート記述)/ Tab.2 脚注 | verified |

### P. 設定値・規則の台帳化(Methods 起草用 2026-08-13 — 出典はコミット済みコード、計算なし)

| N-ID | 値 | 定義 | 出典 | 使用箇所 | 状態 |
|---|---|---|---|---|---|
| N-66 | GENCODE v36(GDC 参照 GTF)の exon-union 遺伝子長; 発現は STAR - Counts(open access) | 参照アノテーションとデータ種別 | scripts/020 ヘッダ(GENCODE v36)・scripts/010 ヘッダ(STAR - Counts) | Supp Methods | verified |
| N-67 | stranded 列合計の小/大比 ≤ 0.5 なら stranded(大きい列を採用)、それ以外は unstranded 列 | ライブラリ strand 判定規則(検体別) | scripts/120 ヘッダ(Dobin's general rule, ratio form) | Supp Methods | verified |
| N-68 | NIH IREP(v5.7.3 引用)甲状腺モデル「AS associated with the expected value of ERR」; 電子 E>15keV、**曝露率 acute**(既定値なしの必須選択 — 研究者確認 2026-08-16。公表済み Chernobyl POC 論文2本も acute)、**臓器線量は定数(Constant (value))として入力**(既定値なしの必須選択 — 研究者確認 2026-08-16。それ以上の線量情報を持たないため択一)、性別は記録値のまま、線量 cSv=mGy/10、被曝年 1986、出生年=1986−被曝時年齢、手術年=出生年+手術時年齢; 残りの IREP 設定は既定値(ユーザー定義不確かさ分布 Lognormal(1,1)・反復 10,000・乱数シード 99)。値は研究者計算の正準 CSV | IREP 入力規約 | scripts/130 ヘッダ; 来歴・入力同一性監査は計画v2 B.11 | Supp Methods | verified |
| N-69 | dose 0 → Sporadic; 0<AS<33.3 Low; 33.3≤AS<66.6 Mid; AS≥66.6 High(境界例なし、2026-07-28 検証) | AS 帯規則 | config.R:36-43(AS_LOW_MAX 33.3・AS_HIGH_MIN 66.6)・lib/cohort_design.R:71-73(境界の帰属: 33.3 は Mid 側 — 2026-08-16 照合で本行の旧記載を訂正、本文は実装と一致していた) | Methods | verified |
| N-70 | 0.6 | pooled 共通尺度の相対純度閾値(main BM 採用条件) | config.R:45-46(PURITY_THRESHOLD) | Methods | verified |
| N-71 | iDEGES 3 反復; スクリーン = 置換 BM(exact)→ Storey q<0.10(plug-in λ=0.5、DEGES_FDR)+ floorPDEG 0.05(q 閾値集合と生 p 上位 5% の大きい方を採用); スケーリングは MUREN(lts); 前処理 protein_coding → filterByExpr | 310 正規化の設定(2026-08-14 訂正: 従前の「BH q<0.10」は誤記 — 実装は storey_q) | lib/norm_deges.R:137-153(storey_q・floorPDEG 採用規則)・scripts/310:37-46(ITERATION 3L・FLOOR_PDEG 0.05・MUREN_METHOD lts・BM_METHOD exact)・config.R:57(DEGES_FDR 0.10) | Methods | verified |
| N-72 | ランキング = 符号付き BM 統計量の tie-averaged normal scores; ES = gseaParam=1 の block 評価(tie-free 入力で標準 GSEA と一致); 推論 = per-set 符号条件付き置換 p + family 内 BH、q_bh<0.10; size 窓 15–500 | 420 セットレベル推論の設定 | scripts/420 ヘッダ; lib/gsea_collections.R:27-28(15L/500L); config.R:51-53(FDR_CUT 0.10) | Methods | verified |
| N-73 | 対比内独立ラベルシャッフル 9,999 回、対比別 seed(基底 19450809)で shuffle プロファイル同士を相関(410 の perm_index は意図的に不使用 — 等 n 対比で同一 index となるため) | 署名一致の帰無参照区間の設定 | diagnostics/signature_agreement.R:57(AGREEMENT_SEED_BASE)・:84,97(N_PERM 使用)・ヘッダ | Supp Methods | verified |
| N-74 | 候補順位 = \|effect−0.5\| 上位500、effect方向に分割後に全構築標本で finite log2 TPM の遺伝子だけを残し全 up×down ペアを評価; dead zone \|r\|<log2(1.2); R0 = 非 dead-zone 検体の多数符号(逆符号例外≤1)かつ q10(\|r\|)≥log2(1.5); R1 reversal rate = 逆符号かつ非 dead-zone の件数 / R1全例、>50%かつ<100%; 順位は群間median r差の絶対値 | REO 候補ペア選定規則 | scripts/510:44(N_CANDIDATES 500L)・:110-151(方向分割、finite-gene screen、全cross-directionペア、規則、順位); config.R:48-49(DEAD_ZONE log2(1.2)) | Supp Methods | verified |
| N-75 | 貪欲選定: 遺伝子再使用禁止 + 既採用ペアとの Spearman <0.75、目標 10 ペア; 530 の Mid>Low 片側 BM は mc・seed 19860426(正準シード) | REO パネル確定と評価検定の設定 | scripts/520:30-31(TARGET_PANEL_SIZE 10L・CORRELATION_THRESHOLD 0.75); scripts/530:78(seed = SEED); config.R:5-8 | Methods | verified |
| N-76 | 群内純度順位保存 0.93–0.99(群別推定との比較)、High vs Sporadic 腫瘍プロファイル相関 0.99 | 220 の両群プーリング妥当性の**来歴実測値**(設計時測定。一次ログではない — 証拠階層は設計選択の来歴)。二重統計に当たるためライセンスには使わない(研究者決定 2026-08-13)— プーリングの記述は理論構造(相対尺度の共通化・driver 支配前提・軸の別・役割の限定)で立てる | `git show 8eed384:scripts/220_estimate_tumor_purity.R` ヘッダ行8-10(run 時点の凍結版。現行ヘッダは 2026-08-13 に理論ライセンスへ書き換え、実測値の記録はこの行と run コミットに保存) | **不使用**(査読応答の受けとして保存) | verified |
| N-77 | 906/906 ライブラリが stranded_second(reverse)判定、比 0.056–0.110(生値 min 0.0558948828・max 0.1095053129、閾値 0.5) | strand 判定の実測結果 — 閾値の恣意性が実務上不活性であることの開示 | meta/strand_selection_20260722_073758.tsv の selected・ratio 列の全数集計(906 行) | Supp Methods | verified |
| N-78 | PC-OD(Nakayama, Yata & Aoshima 2024, JJSDS 7:739–766)= 第1主成分の検体ローディングへの反復 Grubbs 検定(Grubbs 1969; 両側、α = 0.05)、棄却毎に最大絶対ローディングの検体を除去して再計算、棄却なしで停止; 入力は群 × 組織の部分行列(filterByExpr 縮約・未正規化 log-CPM、prior.count 2) | 210 外れ値スクリーンの規則一式 | lib/qc_pc_od.R:1-39(α 既定 0.05・Grubbs 臨界値・反復規則)・scripts/210:69-81(run_pcod: filterByExpr・log-CPM 設定)・210 ヘッダ(8 部分行列) | Supp Methods | verified |
| N-79 | コンテナビルドの日付固定: 基底 ubuntu:noble-20260410(不変タグ)+同日 apt snapshot(20260410T000000Z)、R パッケージは P3M 日付スナップショット 2026-04-09(CRAN + Bioconductor 3.22)、R 4.5.3 は参照 BLAS でソースビルド、OpenBLAS 非導入は設計 | 環境再構築の固定機構(再現性節) | Dockerfile:5-16(ヘッダ)・:29(FROM)・:33-37(apt snapshot)・:82-88(P3M/Bioc 3.22) | Supp Methods | verified |
| N-80 | BHLHB9(ENSG00000198908.12)は chrX(102,720,688–102,753,540; GENCODE v36) | B_Normal 唯一の DEG の性染色体所属(性候補の開示) | raw/reference/gencode.v36.annotation.gtf の gene 行(2026-08-15 実測照合: chrX・gene_name BHLHB9) | Results / Disc | verified |
| N-81 | High 群の AS 中央値: R_High 86.7 [68.2–93.7] / B_High 75.9 [68.5–93.6](参考: 線量中央値 583 / 526 mGy) | **本文不使用**(2026-08-16: B_Normal 段落の軽量化により撤去 — 発見1遺伝子に対し重すぎる論証だったため。算出値は Q-16 の台帳回答用に保存)。旧用途 = C-16 の「R 側が曝露指標上で高い」論拠(追補計算 2026-08-15、研究者 Go) | processed/thyr_analysis_cohorts.rds(include_main_bm、driver × band=High の assigned_share/dose_mgy 中央値 — 正準イメージ rebc-r453:refblas 内で算出。R_High 86.7 は N-50 と一致) | Disc | verified |
| N-82 | R_Tumor 発見 1,765 遺伝子中 chrX 57(被曝群で高発現 36・低発現 21)・chrY 0; 検定対象 15,621 中 chrX 572・chrY 2 | R_Tumor 発見リストの性染色体所属の開示(Methods の注釈宣言の受け、性候補の開示 — 員数のみ、検定なし) | diagnostics/sex_chromosome_annotation.R(GENCODE v36 GTF gene 行 × thyr_expression_test.rds、正準イメージ rebc-r453:refblas 内で算出 2026-08-16。B_Normal 唯一 DEG=BHLHB9 chrX の N-80 一致を stopifnot で確認) | Results | verified |
| N-83 | IREP 開示定数: REF(electrons E>15keV)= 単一値 1.0(「Single-valued at 1.0 (assumed to be same as value for reference higher-energy photons)」— 区分は Chronic or acute 共通)。放射性核種β線は**平均エネルギー**で区分に割り付ける(脚注 c は E<15keV 側にあり「including average energies of beta particles emitted by radionuclides」— I-131 平均β ~192 keV は E>15keV 区分)。DDREF: chronic は全線量域に適用/acute は D ≥ D_L で 1(D_L = 0.03–0.2 Gy の対数一様、以下はロジスティックで chronic 値へ接続、D_L 時点 >0.99)。甲状腺・乳腺は DDREF=1 に大きい確率の離散分布(平均 1.6 は Kocher 2008 帰属)。※Table IV.H.1 脚注 b は「0.2 cGy」と表記し本文 IV.F の 0.2 Sv と単位不整合(原文誤植疑い — どちらの読みでも本文の論法は成立、採用は Gy 読み = 保守側) | IREP 内部被曝開示の定数(Methods 開示段落。本文の数値は 1.0・0.03–0.2 Gy) | Land et al. 2003(NIH Pub 03-5387)全文 PDF — dceg.cancer.gov/tools/analysis/ircp/nih-no-03-5387.pdf、md5 643161533fb7e58c9f9e6d7c2d91407c(CDC docket の wayback 版と同一)。逐語照合 2026-08-16(ASCII85+Flate 全展開、研究者指示の先行確認) | Methods | verified |
| N-84 | 帯別線量: Low n=100 min 16 中央値 41(<200mGy 100 = 全例)/ Mid n=82 min 48 中央値 112(<200mGy 73)/ High n=31 **min 188** 中央値 526(<200mGy 3); 主コホート High n=24 min 188(<200mGy 1) | AS 帯別の線量分布(IREP 開示段落の High 帯最小線量の出典。P(D_L>0.188)≈3.3% は対数一様 ln(0.2/0.188)/ln(0.2/0.03) — 手元導出、本文不使用) | processed/thyr_analysis_cohorts.rds(band × dose_mgy、rebc-r453:refblas 内で算出 2026-08-16、研究者 Go の確認計算) | Methods | verified |

## 図表台帳(図・表は1行ずつ — キャプションも検査対象)

| 列 | 内容 |
| --- | --- |
| F-ID | Fig.1, Tab.1, Supp.Fig.1, … |
| 生成 | 生成スクリプト(figures/xxx.R 等)と入力 rds |
| キャプション内の主張/数値 | 対応する C-ID / N-ID(キャプションも本文と同じ検査を受ける) |
| 状態 | draft / verified |

番号(Fig.1 等)は仮置き。図表構成の確定時に振り直す。C-ID は主張マップ記入後に対応付ける。

**注記(2026-08-21)**: 図表の正式番号・生成スクリプト・出力の正は figures/manifest.csv
(全17行が実体を指す完備状態)。本節の F-ID(表5(仮)等)は旧仮置きのまま — 本台帳は
足場につき行の振り直しは行わず、参照は manifest とキャプション節に従うこと。

| F-ID | 生成 | キャプション内の主張/数値 | 状態 |
| --- | --- | --- | --- |
| Fig.1(仮)遺伝子別 BM 証拠 | figures/fig_gene_bm_evidence.R ← processed/thyr_expression_test.rds + thyr_se_raw.rds。正準 repo/output/figures/fig_gene_bm_evidence.png(二機バイト一致) | C-未、N-16・N-18(付随ログ fig_gene_bm_evidence.log:2-5) | draft |
| Fig.2(仮)MA プロット | figures/fig_ma_gene_bm.R ← processed/thyr_normalized_counts.rds + thyr_expression_test.rds + thyr_se_raw.rds。正準 repo/output/figures/fig_ma_gene_bm.png(2026-08-10 01:40 版。B.14 で二機 md5 一致) | C-未、N-16・N-49 | draft |
| 表5(仮)事前指定の解釈マップ(4対比の位置づけ+パターン規則) | 実体 = 本文内表(「Figure legends and table captions」節・表1(旧表2)キャプション直下 — 2026-08-18 移設、Methods は参照文のみ)。出典 = 計画v2 §0.6 批准マップの転記(スクリプトなし・数値なし)。2026-08-15 掲載(判断点5 の置き場所決定) | C-01/C-04/C-05/C-07/C-16 の読みの参照先 | draft |
| Fig.3(仮)REO グレーディング | figures/fig_reo_grading.R ← processed/thyr_reo_panel.rds + thyr_reo_evaluation.rds + thyr_se_raw.rds + thyr_case_assigned_share.rds。正準 repo/output/figures/fig_reo_grading.png(二機バイト一致) | C-未、N-41・N-50 | draft |
| フロー図(仮)コホートフロー(6段、**両 driver 層並記** — 2026-08-14 研究者承諾: 層別途中経過は本文でなくこの図が担う) | figures/fig_cohort_flow.R ← thyr_cohort_flow.rds(230 出力)。数値は N-08(合算+RET/BRAF 別)。**実行済み 2026-08-15**(rebc-r453:refblas・N 照合一致) | C-10、N-08 | draft |
| Tab.1(仮)臨床コホート構成(Driver × AS 帯 × pair) | tables/tab_cohort_composition.R ← thyr_clinical.rds + thyr_case_assigned_share.rds + SE | C-未、N-11 | draft |
| Tab.2(仮)解析コホート症例特性(群×帯: n・性別・年齢・サブタイプ) | tables/tab_case_characteristics.R ← thyr_analysis_cohorts.rds + thyr_clinical.rds(§15 B と同経路。**実行済み 2026-08-15**(rebc-r453:refblas・N 照合一致)。初回実行で N-09/N-12/N-13 と照合) | C-未、N-12・N-13 | draft |
| Tab.3(仮)410 主結果(対比 × 検定遺伝子数・DEG・HC p) | tables/tab_gene_level_summary.R ← thyr_expression_test.rds(**実行済み 2026-08-15**。N-15〜N-20 照合一致) | C-未、N-15(遺伝子数)・N-16・N-19・N-20 | draft |
| Tab.4(仮)420/D6 サマリ(対比 × コレクション) | tables/tab_set_level_summary.R ← thyr_enrichment_test.rds + diagnostics/output/gsea_null_calibration.rds(**実行済み 2026-08-15**。N-24/N-27/N-28 と照合) | C-未、N-24〜N-28 | draft |
| Supp.Tab.1(仮)D6 較正 16 セル表(m_sets・p_any・95%CI・mean/max 発見数) | diagnostics/gsea_null_calibration.R ← 310 正規化。正準値 = run/xeon_final_20260811/logs/d6_calibration.log(rds は diagnostics/output/gsea_null_calibration.rds) | C-06、N-24(注記に N-25, N-26) | draft |
| Supp.Data.1(仮)420 全セット結果の完全開示(全対比 × family の pathway・size・ES・NES・p・q_bh 全量) | tables/supp_data_420_full.R ← thyr_enrichment_test.rds(整形のみ・計算なし。**実行済み 2026-08-15**)。Q-08/Q-10 予備線「可能性の棚は開示で遺す」の実体 | C-05 | draft |
| Supp.Tab.2(仮)R_Tumor DEG の ORA 注釈(family × list の全結果) | scripts/430_annotate_deg_ora.R ← thyr_expression_test.rds(430 編入 2026-08-13。正準 = processed/thyr_deg_ora_annotation.rds + output/430_annotate_deg_ora.log。旧 diagnostics 版と同値確認済み) | C-14、N-59〜N-62 | draft |
| Supp.Tab.3(仮)層間 concordance(normal / tumor 両ペア: rho・参照区間・p の2行+識別不能注記) | tables/supp_tab_concordance.R ← diagnostics/output/signature_agreement.rds(**実行済み 2026-08-15** — 初回実行でフィールド名を確認し写像確定)。正準値 = run/xeon_results/logs/signature_agreement.log | C-07、C-17、N-33・N-34 | draft |
| Supp.Tab.4(仮)使用パッケージ版表(パッケージ名・版の一覧 — 2026-08-14 引用整理に伴い追加) | tables/supp_tab_package_versions.R ← docker/versions.tsv(**実行済み 2026-08-15**。実走時の実効版 = run/xeon_provenance/session_info.txt は組版時に併記)。機械生成のみ・計算なし | Methods(Software 節が参照) | draft |

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
- 2026-08-13: 追補計算の第3適用(N-63〜N-65): 群間年齢差の CI つき推定(C-15 で読みを
  実行前に固定。研究者 Go 2026-08-13)。diagnostics/age_arm_difference.R → 同 output/ の
  log+rds(seed 19450809・B=9999・正準イメージ内・exit 0)。群別 n・中央値[範囲]は
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
  群別 run は尺度を分断し単一閾値の意味を壊す(群別で先に回して観察された事実)、
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
- 2026-08-14: 図表整形スクリプト7本を作成(全て整形のみ・計算なし・当時未実行 — 2026-08-15 に全て実行済み):
  figures/fig_cohort_flow.R、tables/tab_case_characteristics.R・tab_gene_level_summary.R・
  tab_set_level_summary.R・supp_data_420_full.R・supp_tab_concordance.R・
  supp_tab_package_versions.R。図表台帳の該当行を実体パスへ付け替え。初回実行時に
  各 N 行との照合を行う(スクリプト内に照合対象を明記)。実行は研究者 Go 待ち
- 2026-08-14: N-77 追加(strand 判定の全数一様性、meta/strand_selection tsv の集計を照合済み)—
  比率規則に引用文献がない代わりに閾値不感応の実測を Methods で開示(研究者 Go)
- 2026-08-14: Methods 監査(研究者指示)に伴う訂正・追補: N-71 を実装に一致(BH→Storey q、
  floorPDEG 0.05・MUREN lts・exact 宣言を追記)、N-37 の R0 表記を「訓練 Sporadic の最大
  逆転スコア」記述へ、N-78(PC-OD 規則一式)を新設。本文側は同監査で A2 件+抜け8 件+
  軽微2 件を反映(ペア解決 _merged 規則・contamDE 前処理と max-one 尺度・REO スコア定義・
  TPM 自前再計算・純度中央値2層ほか)
- 2026-08-14: 深度配分の折衷(研究者決定、BJC/ERC/npj PO 様式を踏まえ Methods=規則+最終 n、
  実現値=Results): 症例フロー数列・群内訳・REO セット員数・REO 構築歩留まり(N-35/N-36)を
  Methods から Results へ移送(N-08 は Methods に n=63 のみ残置、N-35/N-36 は Results §6 冒頭へ
  新設移載)。使用箇所列を N-08/N-09/N-35/N-36 で更新
- 2026-08-14: 用語統一第2弾(unit→対比/contrast)を本台帳にも適用(記述列の「ユニット/unit」。
  出典列の rds アクセサ `units[[u]]` 等は実体フィールド名のため不変更)
- 2026-08-15: 図表整形7本を初回実行(研究者 Go、正準イメージ rebc-r453:refblas、全 exit 0)。
  N 照合全一致(N-08/09/12/13/15〜20/24/25/27/28/33/34/56)。調整2件: fig_cohort_flow に
  cairo 指定+拡張(除外注記・REO 評価セット枝・最終箱4群分割 — いずれも既存 N 値の範囲)、
  supp_tab_concordance は実構造(pairs$normal/tumor)へ写像確定 — 参照区間 = rho_null 中央 95%
  が N-33/N-34 を厳密再現することを実測確認。出力は output/(gitignore、再生成可能)
- 2026-08-15: fig_reo_grading を改修(研究者 Go): ジッター全廃(整数スコア・実 AS に偽変動を
  与えない)、R_Sporadic は AS 未定義として専用ストリップ(積み上げ=員数)表示、未発火だった
  NA→66.6 代入のデッドコードを stopifnot に置換(発火実績ゼロを実測確認 — 15/15 とも実値)。
  N-50 の R_Sporadic を未定義表記に訂正
- 2026-08-15: N-81 追加(High 群 AS 中央値 R 86.7 / B 75.9 — C-16 論拠「leans higher-dose」の
  無数値状態を解消する追補計算、固定入力・研究者 Go)。C-07 に識別不能の論理根拠1文を追加
- 2026-08-16: N-82 追加(R_Tumor 発見リストの性染色体員数 — 保留だった性連鎖注釈の計算を実施、
  研究者 Go)。診断 = diagnostics/sex_chromosome_annotation.R(出力は diagnostics/output/、
  git 追跡外・再生成可能)。読み規則の事前固定は不要と裁定済み(記述的診断)
- 2026-08-16: 引用監査に伴う Discussion 改訂: 年齢バウンドの4観察を解体(Vriens 2011・Coclet 1989
  を撤去 — 前者は 2倍ゲート由来の負で本研究の効果クラスを検出不能かつ「25年幅」が原典に不在、
  後者は成人定常回転の測定で年齢効果を確立しない)。Chatsirisupachai 2021 を C-14 に内容特異の
  形で追加。N-81 は本文不使用へ。主張水準は「画定」から「開示・排除不能」へ低下
- 2026-08-16: IREP 開示の完成(研究者指示の Land 2003 先行確認): 全文 PDF を取得し逐語照合
  (N-83 verified 化 — 放射性核種β線の脚注帰属を訂正、DDREF 規則 D_L = 0.03–0.2 Gy 対数一様を
  確定、脚注 b の cGy/Gy 不整合を記録)。N-84 追加(帯別線量 — High min 188 mGy、研究者 Go の
  確認計算)。本文の chronic 反実仮想を目盛り論法(High 帯上の共通倍率圧縮)に差し替え、
  Q-17 を新設(chronic 論争への5層の受け)
- 2026-08-16: 引用忠実性監査(全引用の出典照合、敵対的検証つき)に伴う台帳訂正2件:
  N-69 の帯境界を実装(lib/cohort_design.R:71-73)に合わせ訂正(33.3 は Mid 側。本文は元から
  実装と一致、境界例なしのため結果不変)/ N-68 に曝露率 acute を追記(IREP の既定値なし
  必須入力 — 研究者確認。監査の指摘で記録欠落が判明)
- 2026-08-16: 判断点5を決着(研究者 Go、BJC Article 規定 = 逐語確認済み 5,000語・主6品目)。
  主 = 図1(フロー)・図2(遺伝子証拠)・図3(REO)・表1(症例特性)・表2(解釈マップ)・
  表3(遺伝子レベル)/補足 = 図S1(MA)・図S2(D6較正 — 新規作図、N-24/25/26 照合一致)・
  表S1〜S8・データ1。本文の(仮)17参照を全置換、SX/SY 解消、Intro のマップ参照は番号なし化。
  作表3本実行(特性・遺伝子・セット — 台帳と全一致)。TabS2(正規化)・TabS4(REO パネル)は
  作表スクリプト未作成の札(manifest に pending 記録)
- 2026-08-18: 投稿規定の原典確認(研究者指示 — Methods 分量配分の検討入力。研究者決定:
  ETJ は注力対象外とし、BJC→ERC→npj PO の三段構えで設計目標は BJC 5,000 語に単一化):
  (1) BJC 現行 GTA(nature.com/bjc/authors-and-referees/gta、2026-08-18 生 HTML 逐語)—
  Article = "5,000 words (excluding abstract, references and figure legends)"・Tables/Figures
  "Max of 6"・structured abstract max 200 words・references "Typically max 60"。2026-08-16
  記録(5,000語・主6品目)と一致し、保留札「現行版 GTA 再確認」のうち語数・品目数は解除
  (残札 = 図2の4パネル=2品目換算の編集部確認、投稿時)。SI 条項 = "The article must be
  complete and self-explanatory without the Supplementary Information"(Methods 分割の設計制約)。
  (2) ERC 著者ガイドライン(journals.bioscientifica.com/erc/pages/author-guidelines —
  Cloudflare により機械取得不可のため研究者がブラウザ取得、Claude Code が全文検分)—
  Research = 6,000 語("Be limited to 6,000 words for submissions")、算入基礎 = "Word count of
  the full article, excluding references and figure legends"(abstract は算入)、abstract は
  単一段落 250 語、図 10・文献 60("as recommended" — 推奨)、Vancouver 数値引用、Supp 許容
  (査読対象・掲載料あり)、二重匿名査読(2026-01-26〜)、語数をタイトルページに表示。
  本日セッション前半の「ERC 5,000 語」(検索スニペット由来)は誤りと確定し本行で訂正。
  (3) npj PO content types(nature.com/npjprecisiononcology/content-types、2026-08-18 生 HTML
  逐語)— 語数・頁数の制限なし("We do not impose strict limits on word count or page
  numbers")、Methods は本文必須("All descriptions of Methods should be provided in the main
  manuscript file and not in the Supplementary Information")、Discussion は小見出し・
  limitations・conclusions 節不可、abstract 150 語、文献 60(緩)、図凡例 350 語/図。
  帰結: BJC 5,000(abstract 別枠)に収めれば ERC は自動適合(abstract 250 込みでも 6,000 に
  余裕)。npj は語数不問だが Methods→Supp が規定違反方向のため、分割は本文へ機械的に
  再統合可能な二層設計とする。原典 PDF(研究者保存 2026-08-18、いずれも
  /mnt/c/Users/kotaro/OneDrive/論文関連（説明用資料含）/ 配下): BJC投稿規定.pdf
  (md5 9a358740e9db29c34e78d56da4266d54)、ERC投稿規定.pdf
  (md5 637978fc9c2ed403dd31f583c461aaba)。※PDF は保存操作の記録であり、逐語の根拠は
  各ページ本文(BJC/npj は取得 HTML、ERC は研究者提供の全文)。
  札(新規): 生成 AI 使用の開示 — ERC は開示文言必須("disclose in a statement whether or
  not they have used generative AI ... the name, version, model, and source")、BJC も GTA に
  "Use of an LLM should be properly documented in the Methods section"。開示は本文 Software
  節に記載済み(Claude Fable 5・GPT-5.6sol — 用途・データ範囲・検証・責任)で、残作業は
  投稿時の各誌書式への転記のみ(ERC = Declarations 文言、BJC = Methods 記載)。当初
  「暫定結果閲覧歴の開示深度」と併記したが成熟度が異なるため分離(次エントリで同件決着)
- 2026-08-18: 「暫定結果閲覧歴の開示深度」札を決着(研究者決定: **論文レベル非開示**)。
  対象: 計画v2 §0.6 決定記録(2026-08-07)が「旧パイプラインの暫定結果(R_Tumor 有意・
  B_Tumor 静か・B_Normal 異常高)を閲覧済みの状態でマップ固定」と記録し、紙面文言
  "before the reported results existed" を超える開示の深さがコレスポ判断として保留されて
  いた件(0815 指示書の未遂リスト由来)。根拠: (1) 閲覧対象は修正前実装の開発時出力 —
  MUREN の fixed-reference/self-reference 意味論・推定量が原理論から逸脱しており
  (56a7b33 "Align optimized MUREN with original theory"、2026-07-23 = マップ固定の2週間前)、
  固定時点で欠陥は既知。報告プロトコルの暫定結果に当たらない(事象の実態は、自作スクリプト
  土台からの積み上げ開発における通常反復での出力閲覧)。(2) 旧観察の最顕著特徴 B_Normal
  異常高は修正後に非残存(報告結果 1 DEG・HC p 0.0773 = N-16, N-23)。(3) マップは全パターン
  事前割当のため、閲覧が読みの割当を歪める経路が構造的にない。期待の根拠は文献列
  (Morton の隙間・Abend/Ory)で旧出力を引かない。紙面の4箇所(Intro 設計段落・Methods
  マップ文・Disc §1 一括宣言・表2キャプション)の文言は正確なスコープで真であり不変更。
  付随処置: Q-18 新設(誠実な査読受け)/計画v2 §0.6 に日付つき文脈追記(事実保持・
  性格づけ訂正)/文言凍結検査を We チェックリスト全域項目に追加(強化系への滑り禁止)。
  検査で発見: Methods の C-16 文に旧形 "pre-assigned before the results were seen" が
  1箇所残存(2026-08-15 統一の取りこぼし、「Contrast-level omnibus」節)— 正準形への
  揃えは★保護文言のため研究者判断待ち(チェックリストで監視)
- 2026-08-18: 表2実体(表+規則段落+対訳該当部)を Methods 本文からキャプション節の
  表2キャプション直下へ移設(研究者 Go — 内容不変の移動のみ、Methods は参照文と
  その対訳を保持)。e1bc7d6 の方針(正本 = 論文内・一方向切り出し・第2コピーなし)は
  位置の更新のみで不変。本文節から図表素材が構造的に分離され、語数カウンタの特例除外は
  不要化。**構造基準値の改訂: 対訳分割により【訳】36 → 37**(以後この基準で検査)。
  figures/manifest.csv Tab2 行・本台帳の図表台帳行を同期
- 2026-08-18: Methods/Results の様式決定(研究者決定 — 共著者の「Results に Methods 的記述を
  厚めに」様式の予備検討): 規範的内容(事前固定宣言・推定量仕様・手続き規則)は Methods 前置を
  維持し、Results 受けは各節冒頭の方向づけ1文(「To assess X, we ... (Methods)」型)まで —
  Results 書き直し時の検査項目とする。根拠(3誌の原典逐語): BJC = Methods→Results 順・
  "The Results section should briefly present the experimental data" / ERC = Methods→Results
  順・"The results should read as a narrative leading the reader through the experiments"
  (物語様式の**許容**であり要求ではない)/ npj PO = Methods 末尾のため Results の方向づけは
  様式上の必然。補強2点: 移送しても語数は浮かない(BJC/ERC とも Methods・Results 同一枠)/
  事前(規則)対実現(値)の境界は本論文の装置であり深度配分(2026-08-14)を維持。
  共著者説明札3点: (1) BJC の Results 定義逐語、(2) 語数利得ゼロ、(3) 事前固定マップ設計に
  おける Methods 前置の装置性(表2 の RR 先例札と同系)。ERC で厚みを求められた場合の譲歩線 =
  方向づけ1→2文の追記のみ(規範的内容は移送しない)
- 2026-08-18: BJC→ERC 転送チェックリスト(リジェクト時の再調査回避のため記録。転送時に照合):
  (1) Abstract を構造化200語 → 単一段落 ≤250語へ組み直し(唯一の内容作業)、(2) 引用を
  Vancouver 数値順へ(et al. の前に最低3著者)、(3) 二重匿名化 — タイトルページを
  "Separate file not for review" で別ファイル化・ファイル名匿名化・自己言及の三人称化・
  本文から所属/資金情報を除去、(4) タイトルページに語数表示(full article、references と
  figure legends 除き)+タイトル ≤85字の確認、(5) AI 開示を ERC 必須文言へ転記(name・
  version・model・source)、(6) 体裁 = ダブルスペース・連続行番号・本文 Word 形式・
  合計 40MB 上限、(7) 補足データは査読対象+掲載料あり(不掲載希望は採択前に編集部へ)。
  無作業で通る項目: 語数(BJC 設計 ≈5,250 < 6,000)・図表(≤10 推奨)・文献(≤60 推奨、
  現47)。原典 = ERC 著者ガイドライン PDF(2026-08-18 保存、md5 637978fc9c2ed403dd31f583c461aaba)
- 2026-08-18: Methods 配分表(paper/methods_allocation.md — 足場、フェーズ3後に削除)を
  研究者批准(3決定込み: 対訳 = Supp 随伴+本文新規/文体 = We 混合様式で直接起草・
  チェックリストをバッチ内消し込み/Supp に本文同名小節を同順で新設)。工程の正式名称 =
  「著者指揮下の AI 文章支援+研究者検収」。**AI 開示を文章面へ拡張する決定**(研究者):
  開示文言は用途(デバッグ+著者指揮下の起草・編集)+目的文(統計的正確性 — 推定対象・
  推論スコープ・限定句・用語一貫)+機構列挙(著者指定内容の文章化、台帳3系照合 =
  数値↔解析出力・引用↔原典・主張文言↔批准済み主張マップ、各段階の著者承認)+
  科学的判断の委譲否定 — 形容詞でなく機構の名指しで丸投げ自走と識別する設計(研究者要求)。
  ERC Declarations 用短縮形も同時確定。実装は B1(Software 小節の改稿)。
  既存 AI 開示札(2026-08-18 前段、転記のみ)の内容を本決定で上書き
- 2026-08-18: Methods 圧縮 B1 実施(研究者検収・Go): Data sources 186→96・Normalization
  131→98・Software 326→331(AI 開示拡張 +80 が環境圧縮 −81 を相殺)。Supp Methods に
  同名3小節を新設(移動は逐語)。タグ移動 N-66/67/77/79 → Supp Methods(使用箇所列同期済み)。
  本文 7,129→7,011。**構造基準値の改訂: H3 40・【訳】40**(以後この基準)。
  投影と吸収在庫は paper/methods_allocation.md の予算注記に記録
- 2026-08-18: Methods 圧縮 B2 実施(研究者検収・Go): QC 396→235・Confounders 283→231・
  Between-stratum 123→108・REO diagnostics 135→78(計 937→652 — ★保護文の密度により
  機械的分割の床)。Supp Methods に同名4小節を新設(本文と同順)、Chen et al. 2016 の
  初出引用を QC→Normalization へ移動。タグ移動 N-78/N-73/N-43/N-44/N-46/N-47 →
  Supp Methods(台帳同期は N-73/78 のみ要 — 他は既に Supp 表記)。完全性検査:
  元7小節の全文トークン被覆で欠落ゼロ、本文→Supp ポインタ 25 箇所全解決。
  本文 7,011→6,730。**構造基準値の改訂: H3 44・【訳】44**(以後この基準)。
  現路線の着地見込み ~5,850 → ★文言スリム化パス+最終パスで ~5,000 へ(計画は配分表)
- 2026-08-18: Methods 圧縮 B3・B4 実施(研究者 Go — 機械的分割フェーズ完走): Gene-level
  348→212・Omnibus 201→181・Anchor 100→74/AS 521→341・Gene-set 423→358・REO panel
  411→329。Supp Methods に同名6小節を追加(計14小節、本文と同順)。タグ移動 N-68・N-74・
  N-30 → Supp。完全性被覆検査: 実質の圧縮は π0 理由文の規模修飾「tens to thousands of」の
  1件のみ(復元可、読み通しで判定)、他は保存。Methods 3,764→2,856・本文 7,350→6,221。
  **構造基準値の改訂: H3 50・【訳】50**(以後この基準)。5,000 との残差 1,221 は
  読み通しの追加 Supp 送り+★スリム化パス+Disc/Results/Intro 最終パスが分担(配分表に記録)
- 2026-08-21: **閲覧歴記録の訂正**(研究者決定 — GPT レビュー評価と研究者評価の一致):
  2026-08-18 の決着エントリとその基礎の計画v2 §0.6 決定記録(2026-08-07、AI 起草)は、
  スクリプト開発中の通常の出力確認(共著者・AI 支援との共有を含む通常の協働)を
  「旧パイプライン暫定結果の体系的閲覧」として過大に性格づけていた。特定の4対比パターンを
  取得・検討したという列挙は事実記録として承認されず撤回(研究者証言)。過剰に保守的な
  記録は紙面を実態以上にデータドリブンに読ませる方向に働くため、正確性の要請として訂正する。
  処置: Q-18 全面改訂/claim_map・計画v2 §0.6 に同旨追記/紙面の凍結文言を
  "before the reported results existed" → **"before the finalized analysis produced the
  results reported here"** へ更改(4箇所+対訳同期+チェックリスト更改。精密化であり
  強化でない)。08-18 の MUREN 実装欠陥論拠は主装甲から補助へ降格。旧エントリは残置
  (時系列ログの追記訂正)
- 2026-08-21: C1 裁定(研究者承認 a〜e・f 保留): (a) AI 開示は批准長文(目的文+機構列挙+
  委譲否定)を維持・ポートで復元/(b) 多重性 family 構成を採用 — R_Tumor = 単一 primary、
  他3対比 = 個別の secondary questions(D5 予測マップ方式への追補: マップ・水準列は不変、
  統計的提示の family 定義として導入。表2 B_Tumor 根拠列の Morton 線量-融合根拠は復元)/
  (c) C-11(A/B 非識別性)を Limitations に復元/(d) Q-17 目盛り論法要約1文を本文復帰
  (acute 先例は Supp、Q-17 先回り欄を同時改訂)/(e) Intro 仮説スロットの混合比・バイナリ観
  文言(gpt_review 試案)を批准/(f) 対訳レジームは保留 — 議論余地の解消後に共著者・関係者
  向け別ファイルとして再生成する方針。ポート後のインライン対訳は維持せず、【訳】基準値検査は
  ポート完了時に廃止予定
- 2026-08-21: 遡及記録(08-19/20 の研究者+GPT 作業、改訂メモ未記載分の補完): REO 中間帯
  PC-OD のライブラリサイズを本線と整合化(3a5fbdc 系 — 結果不変: N-43 = 0/17・0/19、C-09 は
  ratified 化)/Ory 2026 多変量署名の照合を追加(853c5b6 系 — N-85/86 新設・C-13 全面改訂・
  Q-07 改訂)/引用台帳をソフトウェア引用込みで拡充(正準イメージ内 citation() 由来、
  Hess 2011 の逐語 verified 含む)/GDC manifest の保全(095699b → N-89)
- 2026-08-21: **C2 ポート完了(B0〜B6)**: 正本を投稿構造化 — 本文 = draft_manuscript.md、
  Supp = supplementary_material.md(単一別ファイル)— し、gpt_review 試案へ全面移植。
  裁定反映: AI 開示は批准長文を復元/C-11 を Limitations に復元/Q-17 目盛り要約を本文復帰/
  chronic 入力規約監査の記述は不採用(実行来歴なし — 理論主張は目盛り要約が担う)/
  表2ブロック(位置づけ列・パターン規則)は批准版を維持。タグ移植 = 本文 N 68種+C 17種、
  Supp N 54種。新 N 行 = N-87〜N-91(N-89 verified、N-87/88/90/91 draft・C3 照合待ち)。
  実測: Abstract 203・本文 5,086(Intro 616 / Methods 2,273 / Results 1,116 / Disc 1,081)・
  Supp 5,714。インライン対訳は全廃(f 方針)。**構造基準値の改訂: 本文 H2 8、【訳】検査は廃止**。
  残作業 = C3(N-87/88/90/91 の一次照合・表S3 実体)+ C4(横断検査・−86語・Abstract −3語・
  足場退役)
- 2026-08-21: **C3 完了**: N-87/N-88 をコンテナ(rebc-r453:refblas、read-only)の rds 直読で
  verified 化(PC-OD フラグ = B_High tumor 1 のみ/RET 主コホート純度中央値 0.783/0.822)。
  0.44 の出典探索の結果、既存 N-58 が同一事実(0.140/0.221/0.44 = B_Normal/radiation)を保持と判明 — 新設した N-90 は重複のため行を削除し Supp のタグを N-58 へ付替(N-90 は欠番)。N-91 は 530 のコード水準で verified(auto→MC・B=999999 既定)。Supp の D6 セル別 CI は N-26 へ付替、N-01(クリーンリポジトリ状態)と N-47(純度層別比較の結果文)を Supp に補完。**表S3 を作表**:
  processed/thyr_deg_ora_annotation.rds の $table(18,576 行)は全セット×3リストを保持して
  おり(GPT レビューの「未保存」評価は誤り)、table_s3_ora_annotation.csv を出力、
  family×list の q<0.10 員数は N-59 と全一致。gpt_review/supplementary_files が全品目完備に
- 2026-08-21: **図表出力の整合検査(全品目)**: S1 = 両分類とも合計 440 ✓/S2 = 群構成
  12/15・27/9 ✓/S4 = 10 行(表の shift 範囲 2.641–4.700 は最終10ペアのもの — 本文の
  1.159–4.700 は適格153ペア = N-36 で別物・無矛盾)/S5 = edgeR 4.8.2・limma 3.66.0・
  msigdbr 26.1.0 ✓/S6 = 16 セル・min q 0.114・セット数 6,141–6,242 ✓/S7 = Σ発見 102・
  B_Normal/H 0.18 [0.110–0.270] 等セル値一致 ✓/S8 = N-33/N-34 と完全一致 ✓/Data1 =
  24,798 行 = 4対比セット数の総和 ✓/図 S1/S2 は output/figures と md5 同一 ✓/S4 凍結
  ソース(output/reo_panel.csv)実在 ✓。**発見・修正1件**: 表1キャプションが relative
  purity を約束していたが凍結済み表体(tab_case_characteristics.csv)に純度列がない —
  キャプションから削除(純度中央値は本文 N-88 が担う)。表に純度列を足す再作表は
  研究者オプションとして残す。**表S3 をスクリプト化**(tables/supp_tab_ora_annotation.R —
  出力は組版コピーと diff 同一)、manifest の TabS3 行を実体へ更新
- 2026-08-21: **図内ラベルの用語同期(案1・研究者 Go)**: 図1・図2・図3 のスクリプトの表示
  ラベルを本文用語へ差し替え再生成(整形のみ・データ不変。ログ値は台帳一致を再確認 —
  63・15,621/1,765・36 = 17/19 等): R_Sporadic 等 → dose-zero / Low-AS / Mid-AS / High-AS、
  training/evaluation → construction/application、"in exposed" → "in High-AS"、
  図1 の枝は「REO application set」へ。**図2サブタイトルに残存していた旧凍結文言
  ("fixed before the reported results existed")を図から除去**(事前固定の主張は
  キャプション・本文が担う — 図内に凍結文を持たせない方針)。3図とも目視検証済み
- 2026-08-21: 表1の純度開示を脚注方式で完成(研究者 Go — 列追加は尺度混在のため不採用):
  B 系の主コホート純度中央値を rds 直読で算出(B_Sporadic 0.836 (n=27)・B_High 0.922 (n=9))し
  N-88 を4群へ拡張。表1脚注に「driver コホート別相対尺度・層内でのみ比較可能」の注記つきで
  4群中央値を記載、Results の純度文は脚注へ移動(本文 −16語)、キャプションに脚注参照を復帰
- 2026-08-21: **Supplementary Data 2(遺伝子レベル全結果)を新設**(研究者起案・Go):
  tables/supp_data_gene_level_full.R(整形のみ)で全検定遺伝子×4対比 = 62,952 行を一方向出力。
  照合 = 行数 N-15・q<0.10 員数 N-16・B_Normal 1件 = BHLHB9(N-22)と全一致。設計理由 =
  データ1(セットレベル全開示)の遺伝子レベル対応物 — 非2値化・R(α) 保持の方針と整合し、
  後続研究の照合(本研究が Abend/Dom/Ory に行った形)を可能にする。manifest Data2 行・
  Supp 凡例・Results §2 参照・組版コピーを追加
- 2026-08-21: Abstract を 200→194 語へ短縮(研究者 Go — BJC 計数器の揺らぎに対する余裕確保):
  4 修正 = IREP 文の語順統合("with NIH IREP" を前置修飾へ)/panel 並列の "and"→カンマ/
  REO 文の "medians" を括弧外主語へ移動/Background 対比を "high estimated radiation-attribution
  ... cases" の複合修飾へ。限定句差分検査 = 脱落なし(estimated・radiation・descriptive・
  unchanged・"ordinal banding only"・construction anchors 全維持)。N タグ(N-15/16/20/41)不変
- 2026-08-21: **AI 開示の二層化**(研究者裁定 — C1-a「批准長文を Methods 末尾に維持」を更改。
  契機 = 利害関係のない専門読者(専門医・学位)の読解反応「過剰な AI 代筆と読まれ得る」+
  規範・先行例の確認: Nature Portfolio AI 方針は段階制で本稿の用途(起草・編集)は開示必須層、
  ただし詳細度は著者裁量/短文ジャンルの範 = Elsevier 指定2文型。裁定原則 = 批准事項でも
  目的に対し客観的により合理的な方法と矛盾すればリトラクトする): Methods は3文(設計・判断の
  著者帰属/ツール名・用途・細粒度指揮/委譲否定+著者管理検証への Supp ポインタ+全責任)へ
  短縮し、目的文・機構列挙(台帳3系照合・各段階承認)・データ範囲・コード検証は Supp Methods
  新設小節 Scope and verification of AI use へ逐語移設(欠落ゼロ、"described above" は
  Software 小節名参照へ付替、Software 小節末尾の open-access 重複文は新小節へ統合)。
  批准内容の全要素は紙面に存続 — 変更は本文/Supp の配分のみ。ERC Declarations 短縮形
  (2026-08-18)は不変、submission_declarations の転記先メモを同期
- 2026-08-21: N-92 新設+Supp Methods(Data sources)へ GDC Data Release 45.0 の1文追記
  (研究者 Go — 身内指摘を受けた版明示。取得日だけでは後年の再照会が別リリースに当たるため):
  紙面は最小形(45.0 と公開日のみ)、検証装置(46.0 の日付・中間リリース不在・帰属根拠 =
  日付照合)は N-92 が保持。46.0 での再照会は不実施(結果がどの決定にも影響しないため —
  sensitivity-column 試金石)。私的保存(PDF 等)は不採用 — 原典が版固定で公開されている
  (gdc-docs コミット固定パーマリンク、出典列に記載)。凍結 manifest md5 との役割分担 =
  リリース番号は状態識別子、md5 は内容固定
- 2026-08-21: N-92 の Supp 文から公開日括弧 "(released 4 December 2025)" を削除(研究者指摘・Go):
  原稿の内部慣行はバージョン付き資源を識別子のみで引用(R 4.5.3・GENCODE v36 等。日付併記は
  日付自体が識別子の場合のみ = CRAN snapshot・noble タグ)であり、released 型注釈は本箇所が
  全文で唯一の異質パターンだった。公開日は 45.0 から一意に復元可能な客観情報で識別力ゼロ。
  日付は N-92 出典側が保持(台帳不変)
- 2026-08-23: Codex CLI(gpt-5.6-sol)による正本2本の文法・文体チェックを全件処理 —
  指摘26+補完1 → 適用21・維持/不採用6・保留1(C-03 同時)。コミット = 42803eb 5bef049
  e56c56e 3b7166b a55043d 6535c50 809c56d(詳細は各コミット)
- 2026-08-25: 言語整合レビュー2本(研究者作成・未追跡)を全件裁定・処理 — コミット =
  6086d69 ecef30b 5409e13 40516d0 81a6919 b7f2063 67ea2cb(詳細は各コミットとレビューファイル)
- 2026-08-25: **記録規約を批准**(研究者 — 過去2プロジェクトの amendment 文脈肥大による
  破綻を受けた構造対策): 改訂メモは1エントリ1〜2行・詳細は git・閉流は圧縮・末尾の
  「現況」節を維持し通常セッションはそこだけ読む。claim_map 運用規則も同期改訂。
  本日の2閉流を上記2エントリへ圧縮(情報損失なし — 全文は git 履歴)
- 2026-08-25: ratified 化同時方針を解消(軽量運用で順序制約消滅 — 研究者裁定)し、裁定済み
  3件を適用: C-03 matched normal 統一・C-11 across・C-07 absence of detectable sharing(+3語)。
  C-04 文のみ (iv-a)/(iv-b) 再裁定待ち(候補 = レビュー 16:49 版 §1-5)
- 2026-08-25: C-04 文を (iv-a) で確定・適用(+5語 — 二重並列の解体と遺伝子/セット水準の混線
  解消、無命名・"lay above"・N-25 維持)。"establishes that" 裁定は対象文言の消滅(C2 ポートで
  軟化済み、残る establish は否定ガード2箇所のみ)により解消、attribution/attributability の
  表記差は性質名 vs 値記述として現状維持で閉鎖(いずれも研究者裁定)。通し読み前の未決文言論点はゼロ
- 2026-08-25: 再走査レビュー(18:21版)の3件を適用: "one matched tumor–normal sample pair"
  (−1語)・HC文の主体明確化("so the contrast-level evidence … did not rely"、+1語)・
  "Relative purity could be estimated for"(−1語、Supp と統一)。2-1(Abstract Conclusions)は
  不採用(標準的絶対構文・余裕3→2・within-driver 限定句の削除が見合わない — 研究者裁定)
- 2026-08-25: 機械走査3層(誤字辞書・語彙照合・LanguageTool 6.6 ローカル)を実施 — 実質指摘は
  Supp の独立節等位コンマ1件のみ(適用、研究者 Go)。他18件は誤検出として根拠つき棄却

## 現況(通常セッションはこの節と直近の改訂メモのみ読めば足りる — 2026-08-25 更新)

- 語数: 本文 4,931(余裕 69)/Abstract 197(余裕 3)。機械検査(凍結4・N/C タグ・引用31/警告0)全通過。
- 開いている残件(文言の未決論点はゼロ — 残りは批准・記入・フェーズ2 のみ):
  1. 通し読み = C-01〜17 ratified 化(研究者。C-09 のみ済)
  2. 条件付き札: Abstract 余裕≥4 回復時に "cases with high estimated radiation attribution"(+2語)を適用
  3. submission_declarations の【記入】6箇所(研究者)
  4. フェーズ2 送り(書式仕様は BJC GTA ライブ照合 2026-08-25 で確定):
     References 投稿規定整形(Vancouver・6著者超は先頭6+et al.・上付き・PMID 併記整理)/
     組版記号(C(n,nx)・n_X)/Word 化 = 1.5 行間・全ページ+全行番号・図凡例は References 後の
     別ページ・Additional Information 7節を規定順/表1〜3 は個別ファイル+キャプション/
     図は tiff・eps・jpg・bmp で 300 dpi 以上/Supp は可能なら一体 PDF+Data 1/2 に各 50 語要約/
     タイトルページ(≤150字・非結論形 — 現115字適合・対応著者 ORCID)/カバーレター(重要性・
     BJC 宣言・未発表言明・対応著者情報・COI 文)
- ratified 行の変更 = 検査経由+改訂メモ1〜2行+研究者 accept(再批准を兼ねる)。
  別格 = 凍結4箇所・表1ブロック(旧表2 — 2026-08-26 引用順へ番号交換)・事前固定に関する記述(紙面の事実主張が根拠)
- 2026-08-25: BJC Guide to Authors をライブ再照合(公開ページ、保存 PDF は文字アウトライン化で
  抽出不能): 既知(5,000/200・品目≤6・文献≤60・Vancouver 6著者+et al.・LLM=Methods 文書化)を
  再確認+新規詳細 = タイトル≤150字・非結論形(現115字適合)/1.5行間+全頁・全行番号/
  表は個別ファイル+キャプション/図凡例は References 後の別ページ/図 tiff・eps・jpg・bmp
  300dpi以上/Supp 可能なら1ファイル+各ファイル50語要約/カバーレター必須要素(BJC宣言・
  独自性・対応著者・COI文)。declarations の節順は規定と一致確認済み
- 2026-08-25: 図5枚に 600dpi LZW TIFF 出力を追加・正準コンテナで再生成(@575f493 — 従来
  PNG は dpi160/200 で BJC 最低300未達だった。N 照合ログ全一致・内容不変)。再生成ログ照合で
  **N-25 の CI 上端の丸め誤りを発見・訂正**(0.270→0.269、rds 実値 0.26947709。本文・Supp・
  N-25 行の3箇所、他セル表示は標準丸めで整合 — 外側丸め仮説は実値で棄却済み)
- 2026-08-25: **Word 化パイプラインを構築・検証**(paper/make_docx.py、pandoc 3.6.4 静的
  バイナリはリポジトリ外): 凡例節の References 後方移動・引用上付き化・1.5 行間・全行番号・
  頁番号フッタを機械注入、docx 本文と md のトークン照合差 0 件・構造検査 4 項目 OK。
  閲覧用コピーは OneDrive 論文関連/word_check/(捨てコピー運用)
- 2026-08-25: References 機械検査(31 件): BJC 著者規則(≤6 全列挙/>6 先頭 6+et al.)に対し
  台帳書誌が「筆頭+et al.」形 17 件・4 名+et al. 1 件・11 名全列挙 1 件、Ory 2026 の
  記事番号欠落 1 件を検出 — PMID 起点の著者リスト補完(フェーズ2 前倒し候補)として保留
- 2026-08-25: docx 組版慣行を実装(研究者裁定): References 1件1段落化(段落融合バグ修正)/
  ヒト遺伝子記号イタリック 119 run(解析ラベル中の RET/BRAF 含む・PTC 単体は疾患略語で対象外)/
  BRAF V600E → 上付き形(BJC 2024–25 掲載例の現行慣行と一致、共著者嗜好とも整合 — 正本不変・
  派生層変換)/見出し均一 12pt(節=太字・小節=太字イタリック)・区切り線除去
- 2026-08-25: **References 著者リストを BJC 規則へ補完**(研究者 Go): 主17件を PubMed
  メタデータ照会(PMID 一括・行↔筆頭著者の一致を機械検証)で更新 — ≤6名は全列挙(4件)、
  >6名は先頭6+et al.(13件)、Ory 2026 に記事番号 16:53030 を補記(DOI・記事ページで確認)。
  引用変換 31・警告 0 不変、著者規則の機械再検査フラグ 0。出典 = PubMed(照会 2026-08-25)
- 2026-08-26: 通読指摘(研究者): Intro の REBC-THYR 初出に記述子を付与 — "the REBC-THYR
  study showed"(+2語 → 本文 4,933・余裕 67。裸のプロジェクト識別子の初出+行為主体化を解消。
  他 3 出現は名詞修飾用法で不変更)
- 2026-08-26: 通読指摘(研究者): 懸垂ハイフン形 "Low- and Mid-AS" 3箇所(本文2・Supp1)を
  原稿優勢の完全表記 "Low-AS and Mid-AS"(12回)へ統一(±0語)
- 2026-08-26: **表番号を引用順へ交換**(通読指摘・研究者 Go — BJC 規定「表は本文引用順」に対し
  旧 Table 2(解釈マップ)の初出が Intro で旧 Table 1(症例特性)より先だったため): マップ =
  Table 1・症例特性 = Table 2 へ全ラベル交換(本文8箇所・±0語)、キャプション節を番号順に
  並べ替え、manifest・N-88 使用箇所・ガード名(表2ブロック→表1ブロック)を同期。凍結断片は
  ホスト文の番号のみ変更で不変。図は元から引用順適合
- 2026-08-26: **共著者版の表提示形式を決定**(研究者決定 — 責任著者判断): レビュー段階は
  3 表とも実体を埋め込んだ一体 docx(論拠 = 版ズレ回避・コメントの単一ファイル集約・表が
  小規模で分離の利得なし・表1 正本の第2コピー禁止規則との整合)。投稿時は BJC 規定どおり
  3 表を個別ファイルへ分離(フェーズ2)。make_docx に凍結 CSV の一方向埋め込みを実装
  (列名は素のまま — 表整形はフェーズ2)
- 2026-08-26: 著者情報を declarations に記入(研究者提供、綴りは PubMed 出版名義で照合 —
  Vladimir A. Saenko・Norisato Mitsutake。所属現名称の最終確認のみ研究者持ち)し、make_docx に
  タイトル直下の著者ブロック描画を実装(一方向 — 投稿時はタイトルページ別建て)。
  【記入】残 = Acknowledgements(親族相談中)・貢献・倫理・COI・資金の5箇所
- 2026-08-26: declarations をほぼ完成(研究者裁定): 所属 = 2025-04 改称名 Molecular Oncology
  and Diagnostic Medicine を採用(共著者検収で再確認)/倫理・Consent・COI・資金 = 承認済み
  定型文/**Authors' contributions = 実態指定どおり**(driver 設定・AS 分割・IREP パラメータ =
  VAS/NM、他全て = KH — 共著者ラウンドで本人ら確認)。残【記入】= Acknowledgements(親族
  相談中)+公開リポジトリ URL/DOI 2箇所(公開時確定)
