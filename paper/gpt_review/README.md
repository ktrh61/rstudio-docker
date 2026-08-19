# GPT manuscript review deliverables

このディレクトリは GPT によるレビューと派生成果物の保全用です。`paper/` 直下の正本を置き換えるものではありません。

## Revised review and manuscript

- `review.md`: user feedback 1–9を受けた再評価。多重性、PC-OD、純度差の評価を訂正し、後に不適切と判明した共変量調整案も明示的に撤回。
- `draft_manuscript_5000_words.md`: 修正版の英語本文試案。BJCの計数対象に相当するIntroductionからDiscussionまでは約4,753語、structured Abstractは見出し込み約190語（Markdown記号を除く空白区切り）。Abstractは本文5,000語とは別枠の200語上限に合わせた。
- `supplementary_methods_for_5000_word_version.md`: 修正版 Supplemental Methods 試案、2,918語。
- `proposed_record_corrections.md`: 開発時の出力確認を過大評価した内部記録の訂正案。正本ファイルには未適用。

## Additional descriptive audit

- `additional_descriptive_audit.R`: 今回作成・実行したRスクリプト。
- `additional_analysis_execution.md`: コンテナ、コマンド、入力ハッシュ、解析境界、主要結果の記録。
- `additional_descriptive_audit.md`: 主要結果の短い要約。
- `additional_effect_size_summary.csv`: 全検定遺伝子と既存 q<0.10 リストの効果量要約。
- `additional_top_gene_effects.csv`: 既存 p 値順上位20遺伝子の順位効果と表示用log2 fold change。
- `additional_sample_pca_covariates.csv`: R_Tumor 27例のPC座標と利用可能な共変量。
- `additional_covariate_summary.csv`: R_Sporadic/R_High別の年齢、純度、性、融合パートナー、PC要約。
- `additional_pca_covariate_summary.csv`: PC1/PC2と年齢・純度の記述的Spearman相関（p値なし）。
- `additional_pca_variance.csv`: PCA分散説明率。
- `additional_adjusted_design_feasibility.csv`: 共変量入り候補モデル行列のrankを確認した監査記録。代数的な推定可能性だけを示し、科学的妥当性または調整効果の識別可能性を示さない。
- `additional_pcod_flag_counts.csv`: 既存PC-OD出力の群別flag件数。

バイナリである `additional_rtumor_pca.png` はこの保全ディレクトリには含めていません。PCA座標と要約はCSVに保存され、図は保存済みスクリプトから再生成できます。

語数は空白区切りの `wc`/`awk` で算出したため、投稿システムの計数とは多少異なる場合があります。引用文献は最終投稿書式へ展開していません。

追加監査は確定済みRDSを読み、記述量だけを出力しました。QC、正規化、遺伝子検定、Storey調整、gene-set検定は再実行しておらず、新しいp値または代替遺伝子リストも作成していません。元の監査実行時にはリポジトリ内ファイルを変更せず、今回この保全ディレクトリだけを追加しました。
