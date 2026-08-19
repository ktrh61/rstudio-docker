# 追加記述解析の実行記録

> **本番不採用:** 本記録中の選択後の効果量要約と `stats::prcomp` PCAは投稿原稿から撤去した。PCAはHDLSS専用法でもPC-ODでもなく、解析法として承認されていない。ファイルは実行履歴としてのみ保全し、本番解析の根拠に用いない。CrossDataMatrixによる置換解析は実行していない。

実行日: 2026-08-19（Asia/Tokyo）

## 目的と解析境界

既存の R_Tumor 推論結果へ、効果量の大きさと利用可能な共変量の記述を追加することを目的とした。QC、純度推定、コホート確定、正規化、Brunner–Munzel 検定、Storey q 値、gene-set 検定は再実行していない。既存 q<0.10 リストを変更せず、新しい p 値、信頼区間、調整 q 値、代替遺伝子リストを作成していない。

## 実行環境

- Image: `rebc-r453:refblas`
- Image ID: `sha256:ca93a538976ff8d076ff3b11d135b675f124225724fb918562fa9066ee2967e6`
- CPU limit: 7
- Repository mount: `/workspace` read-only
- Output mount: `/review` read-write（リポジトリ外の一時ディレクトリ）

保全後の配置から再実行する場合の等価なコマンド（リポジトリルートで実行）:

```bash
review_output_dir="$(mktemp -d)"
docker run --rm --cpus=7 \
  -v "$PWD":/workspace:ro \
  -v "$review_output_dir":/review:rw \
  -w /workspace \
  rebc-r453:refblas \
  Rscript /workspace/paper/gpt_review/additional_descriptive_audit.R /review
```

元の実行では同じイメージ、CPU上限、読み取り専用リポジトリマウントを用い、リポジトリ外の作業ディレクトリを `/review` にマウントした。

初回の作図は、R が X11 なしで構築されていたため `png()` が停止した。出力先以外の変更はなく、スクリプトで `type="cairo"` を明示して全処理を再実行し完了した。

## 入力 SHA-256

```text
76b33e98b84378d68985e8c18d49e2c31d9db5bb9bcc66f46b9d3de86b5dc459  processed/thyr_expression_test.rds
c6cd587a8f1ef7a1351c50681ca0315d19a500665a7fd287f61545c0d2db5f8b  processed/thyr_normalized_counts.rds
a5060c51683b2f363aac994563982e0ba23d8f00f6643657be2c022b8c690480  processed/thyr_analysis_cohorts.rds
1ebf60eb939d4e7853cd1a8426bf58dcfb3401d760721b7e06ba4213ef471621  processed/thyr_clinical.rds
e2a6db87454086ccbfaf8f2b1a01039f3af758dff91744707f3e30db23362786  processed/thyr_case_outliers.rds
```

## 出力の要約

- 既存 q<0.10 リスト: 1,765遺伝子（Highで高い971、低い794）。
- 絶対符号付き順位効果: 中央値0.600、IQR 0.567–0.644。
- 絶対表示用 mean log2 fold change: 中央値0.422、IQR 0.290–0.698。
- 全15,621検定遺伝子 PCA: PC1 19.7%、PC2 16.6%。
- 主解析用相対純度中央値: R_Sporadic 0.783、R_High 0.822。
- 共変量入り候補モデル行列: n=27、7列、rank 7、残差自由度20。これは代数的なrank確認の記録に限られ、科学的に妥当な調整モデルの実現可能性を示さない。遺伝子モデルは未実行。
- 既存 PC-OD flag: RET tumor 0、RET normal 0。4主対象群中のflagは B_High tumor 1件のみ。

各数値の完全な表は `additional_*.csv`、簡潔な要約は `additional_descriptive_audit.md` に保存した。PCAのバイナリ図は保全対象外とした。座標と要約は履歴用CSVに残るが、上記スクリプトとともに本番用途では再利用しない。
