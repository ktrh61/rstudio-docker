# decisions

operational decision（方針として運用者の裁量で決めた、supersede で覆り得る設計選択）の単一情報源。

> 各エントリの `参照: [records](records.md) 「…」` は詳細の**所在標識**であり、AI はその内容を保持しない。AI は records を自発参照せず、詳細が必要な場合は当該区間の提示を運用者に要請する（取得手順の本体は [rules](rules.md) B-4）。

## 運用規約

- 本ファイルは **operational decision**（方針として裁量で決めた設計選択）の単一情報源である。**research knowledge**（外部の証拠＝文献・データ・統計理論が決める研究知見）は本ファイルに置かず、[profile](profile.md) 区画2 に置く。
- **1決定=1エントリ。追記のみ。** 確定後のエントリ本文は書き換えない。
- 決定を覆す場合は、**新エントリを追加**し、旧エントリの状態行を `Superseded by D-0XX` に変更する（本文は保持）。これが唯一許される既存エントリの変更である。
- 状態は3値：
    - `Accepted`（確定）
    - `Pending-verification`（方針は決定済み、実行検証・確認で確定する。昇格条件を明記）
    - `Superseded by D-0XX`（後続決定に置換）
- **論拠は簡潔に（原則一行）。decisions を肥大させない。** records 側に置くべき詳細（検証記録・negative data・経緯）を持つ決定に限り、その詳細は [records](records.md) に置き、当該エントリから参照する。一行で自足する決定に records 詳細や参照を用意する必要はない。
- 履歴は Git と supersede 連鎖が担う。**本ファイルに changelog を設けない。**
- 一方向フロー：[worklog](worklog.md) の課題が解決 → decisions を発行（または Pending-verification を Accepted へ昇格）。実施が完了 → worklog の手順を一行化し、記録すべき詳細があれば [records](records.md) へ移す。

---

### D-001  R 4.5.3 採用

- 状態: Accepted
- 決定: R は 4.5.3 を使う（4.5.1 を踏襲せず、4.6 系にも上げない）。
- 論拠: production 定石（最新の一つ前のマイナー系列の最終安定パッチ）を採り、査読で最終パッチでない理由を問われない状態にするため。
- 参照: [records](records.md) 「R 4.5.3 選定の詳細」

### D-002  Bioconductor 3.22 採用

- 状態: Accepted
- 決定: Bioconductor 3.22 を採用し、`version="3.22"`（実装形 `R_BIOC_VERSION="3.22"`）を明示する。当初は現結果の再現を優先し 3.21 としていたが、現結果保持の圧力を外し production 定石を一貫適用する立場へ移行して 3.22 に改めた。
- 論拠: R 4.5.3 期に current だった世代で、4.5.3 期の修正を取り込んだ枯れた状態にある。
- 参照: [records](records.md) 「Bioconductor 3.22 の詳細論拠」

### D-003  パッケージ版固定＝P3M 日付スナップショット 2026-04-09

- 状態: Accepted
- 決定: Posit Public Package Manager の 2026-04-09 スナップショットで Bioconductor・CRAN 双方を同日固定する。
- 論拠: 世代指定だけでは固定されないパッチ版を決定論的に固定し、install 時刻による版ずれを排して再現性を担保するため。
- 参照: [records](records.md) 「P3M 2026-04-09 確定の全経緯」

### D-004  apt 層の日付固定＝snapshot.ubuntu.com 2026-04-09・Level 2

- 状態: Superseded by D-018
- 決定: 新規 apt install の取得元を `snapshot.ubuntu.com` 2026-04-09 に固定する（Level 2＝新規 install のみ固定、`apt upgrade` はしない。取得可能であることを前提とする）。
- 論拠: R パッケージ層（P3M 2026-04-09）と apt 層の再現性水準を揃えるため。

### D-005  ベースイメージ＝rocker/rstudio:4.5.3

- 状態: Superseded by D-019
- 決定: 器のベースを `rocker/rstudio:4.5.3`（noble, R 4.5.3）とする。当初は `bioconductor/bioconductor_docker:RELEASE_3_22` を想定していたが、同イメージの実 R は 4.5.2 で production 定石と不一致、かつ本プロジェクトは CRAN＋野良改変を含み純 Bioconductor でないため、これを置換した。
- 論拠: R 版を主軸に据える。検証結果・パッケージ版の詳細は records 参照。
- 参照: [records](records.md) 「rocker/rstudio:4.5.3 採用の検証詳細」

### D-006  ディレクトリのフラット化

- 状態: Accepted
- 決定: `analysis_v7/` 接頭辞を廃し、リポジトリ直下を作業ルートとする（scripts/ data/ output/ utils/ config/ logs/ を直下）。source は接頭辞なしで統一する。
- 論拠: `analysis_v7/` は成り行きで付いた不自然な接頭辞であり、公開時に構造をまたぐパスを避けるため。

### D-007  スクリプト内バージョン表記・手書き日付・相互参照の廃止

- 状態: Accepted
- 決定: スクリプト内の版表記・手書き日付・他スクリプトへの版言及を削除し、Git に委譲する。連番（01, 02…）は実行順序として保持する。
- 論拠: 二重管理回避、および版表記が生成 AI の意図しない繰り上げを誘発するリスクの排除。

### D-008  開発/公開分業

- 状態: Accepted
- 決定: 現リポジトリは開発用（作業場）として維持する。公開用は後工程の別成果物として別途立てる。
- 論拠: 作業ツリーの清潔さと履歴保存は独立の軸。備考: 開発側でもスクリプト内部品質（パス規則性・フォールバック排除・暗黙依存除去）は是正する。

### D-009  環境情報記録は専用スクリプトで自動出力・Git 追跡

- 状態: Accepted
- 決定: sessionInfo・主要パッケージ版を専用スクリプトで自動出力する（setup.R に混ぜない）。
- 論拠: separation of concerns。備考: 番号・命名・位置は未決（worklog A）。

### D-010  metadata .rds は Git 追跡外・来歴はテキストで別途追跡

- 状態: Accepted
- 決定: metadata の `.rds` は `.gitignore`、再現性メタ情報は人間可読テキストで自動出力し Git 追跡する。
- 論拠: バイナリは diff 不能で Git 管理の主目的（変更の可視化）を果たせない。

### D-011  manifest 固定名化・fail-fast

- 状態: Accepted
- 決定: manifest を固定名 `manifest_gene_counts.tsv`、metadata を固定名 `metadata_gene_counts.rds`（追跡外）で出力する。02 は固定名を直接読み、無ければ `stop()`（フォールバックなし）。retry manifest は `.gitignore`。
- 論拠: タイムスタンプ増殖の解消・fail-fast。備考: 実装は worklog 4-2 の 01/02 作り直しで行う。

### D-012  プロセス文書を寿命で5分割

- 状態: Accepted
- 決定: プロセス文書を decisions/worklog/records/rules/profile の5ファイルに、役柄でなく寿命で分割する。
- 論拠: 同一の生きた事実が複数ファイルに実体として存在する二重管理を構造的に防ぐ。

### D-013  全文書を docs/ 配下で Git 管理・GitHub 越しにプロジェクト同期

- 状態: Accepted
- 決定: プロセス文書を docs/ 配下で Git 管理し、GitHub 越しにプロジェクトへ同期する（同期一系統）。旧方針「プロセス文書は Git 追跡外・手動でプロジェクトナレッジ更新」を supersede する。
- 論拠: Git（スクリプト）とプロジェクトナレッジ（文書）で同期機構が二系統だった状態自体が二重管理の温床だった。

### D-014  ファイル名は接頭辞なし・版サフィックスなし

- 状態: Accepted
- 決定: docs/ 内のファイル名に `REBC-THYR` 接頭辞と版サフィックスを付けない。
- 論拠: `REBC-THYR` は GDC Portal 由来でプロジェクト固有でなく、docs/ 配下ではパスで所属が自明。版は Git と supersede が担う。

### D-015  決定履歴は supersede 連鎖で表現・文書内 changelog 廃止

- 状態: Accepted
- 決定: 決定の履歴は decisions の supersede 連鎖と Git で表現し、文書内の「改訂の要点」等の changelog を廃止する。
- 論拠: 日付後継ファイル方式・Git と文書内 changelog の履歴機構二重化を解消する。

### D-016  BLAS スレッド1固定

- 状態: Pending-verification（昇格条件: worklog 4-2 で結果への影響の最終評価。実装済みだが評価待ち）
- 決定: norm_improved.R・05・06・07 で BLAS スレッドを1に固定する（`blas_set_num_threads(1L)` / `OPENBLAS_NUM_THREADS=1`）。
- 論拠: マルチスレッド BLAS の非決定的丸め順序を排除する再現性措置。
- 参照: [records](records.md) 「5-4 BLAS/RcppArmadillo」

### D-017  contamDE の p 値調整を qvalue→BH 是正

- 状態: Pending-verification（昇格条件: worklog 4-2 で該当スクリプト是正・結果影響評価）
- 決定: 純度推定の p 値調整を qvalue から BH（`p.adjust(method="fdr")`）へ戻す。
- 論拠: オリジナル `contamDE.lm` 忠実性とプロジェクト方針の二重に正当。
- 参照: [records](records.md) 「5-2 contamDE」

### D-018  apt 層の日付固定＝snapshot.ubuntu.com 2026-04-10・Level 2

- 状態: Accepted
- 決定: 新規 apt install の取得元を `snapshot.ubuntu.com` 2026-04-10 に固定する（Level 2＝新規 install のみ固定、`apt upgrade` はしない）。ベースイメージ `ubuntu:noble-20260410`（[D-019](decisions.md#d-019-ベースイメージubuntunoble-20260410-r-453-自前ソースビルドrstudio-廃止)）の apt 層ビルド日と一致させる。[D-004](decisions.md#d-004-apt-層の日付固定snapshotubuntucom-2026-04-09level-2) を supersede。
- 論拠: ベースの apt 層（4/10 ビルド）と snapshot 固定日を一致させ版ずれを排す。P3M（4/09）と1日ずれるが OS が1日新しい＝安全側。
- 参照: [records](records.md) 「ubuntu ベース転換と 4-1b 成果」

### D-019  ベースイメージ＝ubuntu:noble-20260410 ＋ R 4.5.3 自前ソースビルド（RStudio 廃止）

- 状態: Accepted
- 決定: 器のベースを `rocker/rstudio:4.5.3` から `ubuntu:noble-20260410`（日付タグ＝不変。究極の固定は digest 指定）＋ R 4.5.3 の CRAN ソースからの自前ビルドへ転換し、RStudio Server を廃止する。[D-005](decisions.md#d-005-ベースイメージrockerrstudio453) を supersede。
- 論拠: `ubuntu:24.04`（ローリングタグ）は再ビルドで apt 層が動き再現性が崩れる。日付タグ `noble-YYYYMMDD` は特定日ビルドで不変。プロダクトは対話環境を前提としないため RStudio は不要。OMEN で R 4.5.3 ビルド・明示49パッケージ導入・主要版一致を実証。
- 参照: [records](records.md) 「ubuntu ベース転換と 4-1b 成果」

### D-020  開発環境の BLAS 構成＝OpenBLAS pthread（update-alternatives 方式）

- 状態: Accepted（本番の参照 BLAS 移行は別途・未決）
- 決定: 開発環境（OMEN）の BLAS/LAPACK を外部 OpenBLAS pthread に置く（`update-alternatives` で `libblas.so.3`/`liblapack.so.3` を OpenBLAS pthread、priority 100・auto）。R は `--with-blas --with-lapack` でこれを参照。旧 rocker 環境と同一構成。本番で参照 BLAS に揃えるかは別途判断する。
- 論拠: 現行スクリプトが OpenBLAS 前提（旧 Dockerfile で `libopenblas0-pthread` 使用を確認、[D-016](decisions.md#d-016-blas-スレッド1固定) のスレッド1固定もこの層に効く）。開発の便法を本番判断に混ぜない。
- 参照: [records](records.md) 「5-4. BLAS / RcppArmadillo」「ubuntu ベース転換と 4-1b 成果」

### D-021  MUREN improved の参照選択・遺伝子フィルタを論文準拠へ是正

- 状態: Accepted
- 決定: `utils/norm_improved.R` の `include_self` 既定を FALSE から TRUE にし、`refs="saturated"` を自己を含む全サンプル参照（`seq_len(n_exp)`）へ戻す。あわせて `utils/utils_improved.R` の `filter_gene_l` を複合条件（`rowMaxs>=trim` かつ `rowSums(reads>0)>=2`）から `rowMaxs(reads)>=trim` 単独へ戻す。
- 論拠: 論文を一次・オリジナル `norm.R`/`utils.R` を補助として照合した結果、是正前の自己除外・複合フィルタが論文・オリジナル双方から逸脱していたため回復。
- 参照: [records](records.md) 「5-1 MUREN 参照選択」（決着詳細は材料B で同節に追記予定）
