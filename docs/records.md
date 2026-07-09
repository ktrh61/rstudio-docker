# records

追記のみ・write-once。[decisions](decisions.md) の各エントリが参照する詳細論拠、実施記録、negative data、スクリプト乖離分析の詳細を収める。確定後は書き換えない。

---

## 詳細論拠

### R 4.5.3 選定の詳細（D-001）

- **production 定石（主軸・単独で十分に強い）**: 実用分類器（REO パネル）を主張する production 的発想では「最新の一つ前のマイナー系列の最終安定パッチ」が推奨される。4.6 系が最新である現在、4.5 系の最終パッチ＝4.5.3 が該当し、査読対応上「なぜ最終パッチでないか」の追及を受けにくい。この論拠は RcppArmadillo 廃止検討の帰趨に依存せず単独で成立する。中期的に R 4.7 リリース後（4.6.3 が出る時期）には説得力が減じ得るが、投稿はそれ以前の見込み。
- **Bioconductor 世代の固定（主軸）**: R を 4.6 系に上げると対応 Bioconductor が 3.23 以降になり解析パッケージ群（edgeR・limma 等）が入れ替わる。4.5 系に留めることで 3.22 を選択でき、世代を固定できる。RcppArmadillo/BLAS の帰趨と独立に成立する頑健な論拠。
- **補助論拠（降格）**: 4.6 系は C++ デフォルトの C++20 化・C++11/14 打ち切り・非API除去という破壊的変更を含み、自前 C++ 実装（`CDM_fast3_arma_enhanced.cpp`, `contamde_gls.cpp`）を持つ本プロジェクトには移行リスクがある。ただし RcppArmadillo 廃止を検討中であり、廃止した場合この論拠は失われる。BLAS は R バージョンと独立のため 4.6 固有の回避論拠にできない。よって補助に降格した。4.5.1→4.5.3 は同一マイナー系列内のパッチ差で、結果を左右する数値計算の変更は含まない。

### Bioconductor 3.22 の詳細論拠と 3.21→3.22 改訂経緯（D-002）

- **改訂経緯**: 当初は 3.21（現結果の再現を最優先）としていた。しかし MUREN 実装・contamDE 経路の見直しが必要な公算が高く、スクリプト書き直しと結果変化を受け入れる状況にある。これにより「R に production 定石・Bioconductor に結果再現」という基準の不一致を解消し、両者に production 定石を一貫適用する立場へ移行した。
- **論拠**: (1) production 定石の一貫適用。(2) R 4.5.3 期に current であったのは 3.22 であり、この期の修正を取り込んだより枯れた状態にある。3.21 は 4.5.0 期の世代で 4.5.3 期の修正を受けていない。(3) 3.21 で得た現結果（DEG 数等）は 3.22 再実装で変わり得るが、実装見直しに伴う書き直しを前提とする以上、取り立てての代償にならない。(4) `BiocManager::install()` をバージョン無指定で呼ぶと意図せぬ世代が入り得るため、必ず `version="3.22"` を明示し日付スナップショット（D-003）でパッチ版を固定する。

### 三層構成（日付スナップショット＋イメージ凍結＋テキスト記録）の枠組み解説（D-004 を割った先）

再現性は次の三層で担保する。renv は不要（役割が三層でカバーされ、層を増やさないため）。

- **(A) 日付スナップショット固定**: 世代指定（`version="3.22"`）だけではパッチ版は固定されない（3.22 ブランチ内で edgeR 等にバグ修正パッチが降り得る）。P3M の日付スナップショットで取得元を固定すると CRAN・Bioconductor 双方のその日のパッチ版が一意に決まり、いつ・どの環境で install しても同一版が入る。FDA の R 提出パイロットでも採用された標準的手法。→ 個別層は D-003。
- **(B) Docker イメージ凍結**: 版固定した状態を、パッケージ版・システムライブラリ（OpenBLAS）・R 本体まで丸ごと固定する。本番で構築したイメージを `docker save`/`load` で他環境へ移送し、install 再実行による版ずれを回避する（イメージは GitHub 同期外のため Dockerfile の Git 管理＋イメージ移送を併用）。→ 器は D-005。
- **(C) テキスト記録**: イメージ内部は外から見えない。査読者・共著者への提示、Methods 転記、イメージ喪失時の復元のため、パッケージ版を人間可読・Git 追跡可能なテキストで別途記録する。→ 環境記録は D-009、metadata の扱いは D-010。

### P3M 2026-04-09 確定の全経緯（D-003）

- **3.23 リリース日 = 2026-04-29（水）**（Bioconductor Release Schedule / Release Announcement で一次確認。3.23 は R 4.6 対応）。これが「3.22 が current でなくなった日」。
- **固定日 = 2026-04-09**: 公開版 P3M の Setup で Bioconductor 3.22 (R4.5) を選択・Snapshots=Yes とした際、カレンダーの選択可能範囲は 2025-04〜2026-04、最終選択可能日が 2026-04-09。3.22 スコープで到達できる最新スナップショットとして採用。
- **パッケージ別 History 最終日との区別**: edgeR の History は 2025-12-30（4.8.2 追加日）、limma は 2025-11-04 で止まって見えるが、これは各パッケージの最終更新日にすぎない。スナップショット自体は 2026-04-09 まで存在し、日付固定時に各パッケージの当日時点最新版が一意に揃う。
- **edgeR 版と世代の対応（実機で裏取り）**: 3.22 スコープで edgeR は 4.8.2 が最新（Other Versions に 4.8.1, 4.8.0）。付番規則（偶数マイナー=release, 奇数=devel）により 4.8=3.22 / 4.10=3.23 / 4.11=3.24。当初 4/30 スナップショットで観測された edgeR 4.11.0 は 3.24 devel の版で、`R_BIOC_VERSION="3.22"` スコープ指定により混入しない。
- **P3M Setup Instructions 実文面（4項目）**: Repository URL = `https://packagemanager.posit.co/bioconductor/2026-04-09`（エンコード suffix なし）。R startup file に設定:
    - `options(BioC_mirror = "https://packagemanager.posit.co/bioconductor/2026-04-09")`
    - `options(BIOCONDUCTOR_CONFIG_FILE = "https://packagemanager.posit.co/bioconductor/2026-04-09/config.yaml")`
    - `Sys.setenv("R_BIOC_VERSION" = "3.22")`（BiocManager の自動世代アップグレード抑止＝ D-002 の「version=3.22 明示」の実装形）
    - `options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2026-04-09"))`
- **CRAN も同日 2026-04-09 に固定**: P3M が 3.22 に対し CRAN スナップショット `cran/2026-04-09` を自動提示。Bioconductor の CRAN 依存（Rcpp, ggplot2 等）を同日固定し依存関係全体を一意固定する。2026-04-09 は木曜（営業日）でスナップショット存在の見込み。CRAN 側の有効性は install 実行時に最終検証（worklog B）。
- 参考: Bioconductor の日付スナップショット対応は Posit Package Manager 2024.08.0（2024年8月）導入・2025.12.0 で全面対応。「日付スナップショットは CRAN のみ」は 2024年8月以前の認識。

### rocker/rstudio:4.5.3 採用の検証詳細（D-005）

- **変更**: `bioconductor/bioconductor_docker:RELEASE_3_22` → `rocker/rstudio:4.5.3`（Ubuntu 24.04 noble、R 4.5.3）。
- **理由**: (1) `bioconductor_docker:RELEASE_3_22` の実 R は 4.5.2 で production 定石（D-001）と不一致。(2) 本プロジェクトは CRAN パッケージと手元書き換えの野良スクリプト（MUREN `norm_improved.R`、contamDE `gls.cpp`）を含み純 Bioconductor でない。(3) R バージョン選定の根拠はイメージ選定の根拠より査読で説明しやすい。よって R 版を主軸に据えた rocker を器とする。
- **検証結果**（rocker/rstudio:4.5.3 上、使い捨てコンテナ `docker run --rm -i` による多段確認）: P3M 2026-04-09 設定が機能し BiocManager が Bioconductor 3.22 を認識。主要パッケージ版: edgeR 4.8.2、limma 3.66.0、fgsea 1.36.2、RcppArmadillo 15.2.4.1、brunnermunzel 2.0、clusterProfiler 4.18.4、GenomicDataCommons 1.34.1。
- **システムライブラリ暫定リスト（8件）**: `libpng-dev`, `libcurl4-openssl-dev`, `libxml2-dev`, `libfontconfig1-dev`, `zlib1g-dev`, `libicu-dev`, `cmake`, `libcairo2-dev`。`cmake` は `fs` 2.0.1 が libuv 1.52.0 をバンドルし cmake で静的ビルドするため必要（`libuv1-dev` ではない。fs 公式 SystemRequirements で確認）。**検証2パッケージ（clusterProfiler / GenomicDataCommons）限定の暫定リストであり、最終確定は 4-2 のスクリプト単独の依存分析による。**
- **contamDE / MUREN は器の検証対象から除外**（GitHub install ではなく手元書き換えスクリプトを使用するため）。

### ubuntu ベース転換と 4-1b 成果（2026-07-09・OMEN）（D-019 / D-018 / D-020）

- **環境の性格**: rocker/rstudio ベース（[D-005](decisions.md#d-005-ベースイメージrockerrstudio453)）を廃し、ubuntu ベース＋ R 4.5.3 自前ソースビルドへ転換。RStudio Server も廃止（プロダクトは対話環境を前提としない）。OS(apt=4/10)・R パッケージ(P3M=4/09)・R 本体(4.5.3) の三層を固定。本記録は OMEN 開発環境での 4-1b 成果であり、本番（ThinkPad）はこの材料を元に別途構築する。

- **1. ベースイメージ**: `ubuntu:noble-20260410`、digest `sha256:c4a8d5503dfb2a3eb8ab5f807da5bc69a85730fb49b5cfca2330194ebcc41c7b`（Created 2026-04-10）。日付タグ（不変）を採用。タグ名だけで再現性が担保されるが、究極の固定は digest 指定。根拠: `ubuntu:24.04`（ローリングタグ）は再ビルドで apt 層が動き再現性が崩れる。ポイントリリース明示タグ（24.04.4 等）は存在しない。日付タグ `noble-YYYYMMDD` は特定日ビルドで不変。P3M（4/09）近傍として 4/10 を選択。

- **2. apt 層**: 取得元 `https://snapshot.ubuntu.com/ubuntu/20260410T000000Z`（HTTPS、deb822 形式 `ubuntu.sources`、Suites = noble/noble-updates/noble-backports/noble-security）。ベースの apt 層（4/10 ビルド）と snapshot 固定日を一致させ版ずれを回避（[D-018](decisions.md#d-018-apt-層の日付固定snapshotubuntucom-2026-04-10level-2)）。
    - **ca-certificates 導入の順序問題（鶏卵問題）**: 素の ubuntu には ca-certificates が無く、snapshot は HTTPS 強制のため初回だけ TLS 検証ができない。対処: apt の TLS 検証を一時無効化（`Acquire::https::Verify-Peer=false`）して snapshot 4/10 から `ca-certificates` を導入 → 以降は正規の TLS 検証で運用。GPG 署名検証（keyring）は生きているため真正性は担保。この方式なら `ca-certificates`/`openssl`/`libssl3t64` も 4/10 版で揃い、archive 経由（6〜7月版混入）を避けられる。
    - **導入 apt パッケージ（すべて snapshot 4/10）**:
        - 証明書（初回、Verify-Peer 無効化で導入）: `ca-certificates`（依存で `openssl`, `libssl3t64` が入る。すべて 4/10 版、`libssl3t64` の upgrade なし）。
        - R ビルド依存: `build-essential`, `gfortran`, `libreadline-dev`, `libpcre2-dev`, `zlib1g-dev`, `libbz2-dev`, `liblzma-dev`, `libcurl4-openssl-dev`, `libicu-dev`, `libncurses-dev`。
        - BLAS（開発環境）: `libopenblas0-pthread`, `libopenblas-dev`, `liblapack-dev`（[D-020](decisions.md#d-020-開発環境の-blas-構成openblas-pthreadupdate-alternatives-方式)）。`update-alternatives` で `libblas.so.3`/`liblapack.so.3` が OpenBLAS pthread（priority 100、auto）を指す。
        - 作図（cairo、X11 なし）: `libcairo2-dev`, `libfontconfig1-dev`, `libfreetype6-dev`, `libpng-dev`, `libjpeg-dev`, `libtiff-dev`。`--with-x=no` でビルドするが cairo 有効化で png/pdf/jpeg/tiff が出力可能（`capabilities()` で png/jpeg/tiff/cairo すべて TRUE を実地確認）。
        - 補助: `xz-utils`, `wget`。
        - R パッケージビルド用追加（因果確立分）: `libxml2-dev`（GenomicFeatures/rtracklayer 等の XML 依存）、`cmake`（`fs` が libuv を静的ビルドするのに必要、既知）、**`libssl-dev`（今回初判明。R の `openssl` パッケージビルドに必須。無いと openssl→httr→httr2→AnnotationDbi→biomaRt/rtracklayer/clusterProfiler/GenomicFeatures 等が連鎖失敗。rocker には同梱されていたため従来は顕在化せず）**。← 鶏検証の「システムライブラリ暫定リスト（8件）」への追補。

- **3. R 本体**: `R-4.5.3.tar.gz`（`https://cran.r-project.org/src/base/R-4/R-4.5.3.tar.gz`）。SHA-256 = `aa5c1ed4293c7271ac513d654670356ac0e8a6ad5e42be014365d11150b5b8f2`（R Core 公式アナウンス値と一致）。configure: `--prefix=/usr/local --enable-R-shlib --with-x=no --with-cairo --with-blas --with-lapack`。configure サマリ: External libraries に BLAS(OpenBLAS)/LAPACK(in blas)、Additional capabilities に PNG/JPEG/TIFF/NLS/cairo/ICU、shared R library 有効。ビルド: `make -j$(nproc)` → `make install`（real 2m4s、OMEN）。実行時確認: R version 4.5.3 (2026-03-11)、`La_library` が openblas-pthread（`libopenblasp-r0.3.26.so`）。無害な警告: texinfo/LaTeX/browser/PDF ビューア不在によるマニュアル生成不可（解析実行に無関係）。

- **4. Rprofile.site**（`/usr/local/lib/R/etc/Rprofile.site`）: `.libPaths(c("/opt/r-extra-lib", .libPaths()))` で永続ライブラリを先頭に。P3M（2026-04-09、[D-003](decisions.md#d-003-パッケージ版固定p3m-日付スナップショット-2026-04-09)）: `options(repos = c(CRAN="https://packagemanager.posit.co/cran/2026-04-09"))`、`options(BioC_mirror=".../bioconductor/2026-04-09")`、`options(BIOCONDUCTOR_CONFIG_FILE=".../bioconductor/2026-04-09/config.yaml")`、`Sys.setenv("R_BIOC_VERSION"="3.22")`。Bioc 3.22 と 4/09 の対応: 3.22 のパッケージパスは P3M で 4/09 が取得可能な最終日（4/10 で `packages/3.22/bioc` が 404、edgeR 4.8.2 を 4/09 で実測）。config.yaml は 4/09 時点で既に `release_version: 3.23` だが、`R_BIOC_VERSION="3.22"` 明示で 3.22 が正しく解決（`BiocManager::version()` = 3.22 を実測）。

- **5. ライブラリ永続化**: ホスト `~/r-libs-r453` をコンテナ `/opt/r-extra-lib` にマウント（開発環境。本番では再検討）。`.libPaths()` 先頭が `/opt/r-extra-lib`、書き込み権限あり（作業ユーザー UID 1000 ＝ ホスト UID 1000）。

- **6. R パッケージ（明示49）**: `BiocManager` 経由で一括導入（`update=FALSE, ask=FALSE`）。依存含め総216。
    - CRAN(32): R.utils, ROCR, Rcpp, RcppArmadillo, RhpcBLASctl, UpSetR, assertthat, brunnermunzel, circlize, data.table, doSNOW, dplyr, foreach, future, future.apply, ggplot2, glmnet, gridExtra, httr, igraph, iterators, jsonlite, matrixStats, openxlsx, pROC, pheatmap, randomForest, readxl, robustbase, statmod, stringr, tidyr。
    - Bioconductor(17): AnnotationDbi, ComplexHeatmap, DESeq2, GenomicDataCommons, GenomicFeatures, GenomicRanges, ReactomePA, SummarizedExperiment, clusterProfiler, edgeR, enrichplot, limma, msigdbr, org.Hs.eg.db, qvalue, rtracklayer, txdbmaker。
    - 版確認（Bioc 3.22 解決の裏付け）: edgeR 4.8.2 / limma 3.66.0 / clusterProfiler 4.18.4 / GenomicDataCommons 1.34.1（rocker 検証時と一致）。
    - 除外（パッケージ化しない）: MUREN, contamDE — 手元書き換えの野良スクリプト（`norm_improved.R`, `contamde_purity_functions.R`, `contamde_gls.cpp`）として `/workspace` から使用。

- **7. コンテナ構成（開発）**: 名前 `rebc-r453-dev`、常駐 `sleep infinity`（素の ubuntu は CMD 無しだと終了）、マウント `~/rstudio-docker`→`/workspace`・`~/r-libs-r453`→`/opt/r-extra-lib`、作業ユーザー `ubuntu`（UID 1000、noble 標準既存ユーザーを流用。名前は本番で調整可）、ポート公開なし（RStudio 廃止）。

- **未確定・持ち越し（4-2 で育てる）**: (1) R パッケージ用システムライブラリは因果確立分のみ先回り導入済み。明示49は完走したが、4-2 でスクリプトを実際に動かす際にさらに別の system-lib 不足（例: `contamde_gls.cpp` のコンパイル依存）が出る可能性は残る。(2) 暫定 Dockerfile 化はこの材料を元に別途。(3) `contamde_gls.cpp`（RcppArmadillo 経由 C++）のコンパイル・実行確認は未実施（次段）。

---

## 実施記録

### reconciliation 実施記録（2026-07-04）

- **系列反転の判明**: `git fetch` により、リモート `origin/fix/v7-protocol-refinement` が 2025-12-18 に force-push で `b18af57`（macOS `ktrh61-MacBookAir` 由来、12コミット：05 v7.12、REO panel v7.1、ORA廃止、BH変更、GSEA 一本化、スクリプト番号再編 07=deges/08=pca可視化/13=reo_evaluation、旧13/14 archive 化 等）へ更新済みと判明。ThinkPad ローカルの `834f7e5` は macOS 更新を pull し忘れた旧系列（共通祖先 `7ffb3f1` で分岐、ローカル固有6・リモート固有12）。プロジェクトナレッジの主要所見（1145 DEGs, BH, 10 遺伝子ペア REO, ORA廃止）は macOS 系列を反映しており **b18af57 が正**と確定。
- **正系列の最新 = `87a75a0`**（`b18af57` に macOS の `00_prepare_gene_lengths.R` クリーンアップ〔push 後 2025-12-23、冗長な spot check 削除、本番出力に影響なし〕を1コミット追加した状態）。以降の全同期の基準点。
- **三環境＋リモートの正系列統一**:
    - macOS: 00 をコミット・push し正系列確定、以後作業中断。00 以外の untracked（`bioc_snapshot.tar`, `core`, utils 配下のコピー群等）は残置（後続項目）。
    - ThinkPad (WSL2): `git reset --hard 87a75a0` でローカル `fix` を正系列に一致。`v7-data-expansion` はローカル・リモート削除、stash は drop（タグ保全済み）。untracked 遺物なし。完了。
    - OMEN: `git pull --ff-only` で `87a75a0` に一致（ff 1コミット、固有 0）。untracked の v6 遺物3ファイル（`11_feature_selection_final_v6.R`, `12_final_radiation_classifier.R`, `CopyOf11..._v6_optionC_cleaned.R`、本番未参照を検索確認）を手動削除。完了。
- **退避タグ2件**（ローカル・リモート双方へ push、単一マシン依存を排除）:
    - `archive/thinkpad-pre-sync-2026-07-04` → `834f7e5`（ThinkPad 旧系列先端、局所コミット 834f7e5・01a8bea を含む）。
    - `archive/wip-v7-data-expansion-2025-09-17` → `b2f532d`（9月WIP の作業内容、基点 `8fd62c5`。stash が指す宙吊りコミットに注釈付きタグを付与し到達可能参照として保全）。
- **`v7-data-expansion` 削除・stash drop**: ブランチをローカル（`git branch -d`、マージ済み安全形式）・リモート（`git push origin --delete`）双方から削除（正系列に `0` 固有コミットで内容損失なし）。stash エントリを `git stash drop`（中身はタグで保全済み、drop 後も `archive/wip-...^{}` = `b2f532d` を確認）。

### rebuild/r453 作成記録

- ThinkPad WSL2 で `git switch -c rebuild/r453 87a75a0` により作成（2026-07-05）。作成前に作業ツリー clean・HEAD=`87a75a0` を確認。作成後、現在ブランチ=`rebuild/r453`・作業ツリー clean・`fix/v7-protocol-refinement` と `main` の残存を確認。OMEN も 87a75a0 土台で `rebuild/r453` に一致（2026-07-09）。git 操作はユーザーが手元で実行。

### 本記録文書自体の Git 配置に関する旧判断

- 旧判断では、プロセス文書の Git 配置を当面行わず手元保持＋プロジェクトナレッジで保持するとしていた（既存プロセス文書は Git 追跡外と確認済み）。この旧判断は [D-013](decisions.md#d-013-全文書を-docs-配下で-git-管理github-越しにプロジェクト同期)（docs/ 配下で Git 管理・同期一系統）で上書きされた。

---

## negative data

- **snapshot.ubuntu.com の 503 障害**（apt 層＝[D-004](decisions.md#d-004-apt-層の日付固定snapshotubuntucom-2026-04-09level-2) の残点と対応）: 検証時、`https://snapshot.ubuntu.com/ubuntu/20260409T000000Z` の noble / noble-updates / noble-backports / noble-security の InRelease 取得がすべて **503 Service Unavailable**（IP 185.125.189.69）で失敗し、当該実行の system library install が `Unable to locate package` で落ちた。503 の一時性／継続性は未判定。稼働確認が済むまで apt 層の日付固定は確定扱いにしない。
- **snapshot.ubuntu.com 4/10 到達の確認（2026-07-09・OMEN、上記 503 の後続）**: apt 固定日を 2026-04-10 へ移した（[D-018](decisions.md#d-018-apt-層の日付固定snapshotubuntucom-2026-04-10level-2)）際、`https://snapshot.ubuntu.com/ubuntu/20260410T000000Z` から HTTPS で `ca-certificates` ほか R ビルド依存一式の取得に成功（ca-certificates は Verify-Peer 一時無効化で初回導入、詳細は「ubuntu ベース転換と 4-1b 成果」）。4/09 で観測した 503 の恒常性は依然未判定だが、4/10 は少なくとも当日到達可を実測。
- **Claude Science 評価（2026-07-07、現時点で不適）**: 自律実行エージェント設計が明示的承認の統制哲学と衝突、環境管理が凍結済み再現性アーキテクチャと衝突。文献レビュー・原稿ドラフト用途に限れば将来的な可能性はあるが、コアパイプラインとは分離。
- **既に方針として確定した negative data（簡潔に）**: 旧系列（`834f7e5`）は pull 忘れの旧系列でありタグ退避で保全のうえ不使用。macOS で得た数値結果データは正規化の小数点差のため不使用（ただし正系列コードの出所が macOS であることとは別問題で、再現実行は WSL2〔amd64〕基準）。ORA は廃止し GSEA へ一本化。既存 amd64 スナップショットの掘り起こしはせず P71 でクリーン再構築。

---

## スクリプト乖離分析詳細

worklog には課題見出しのみを置き、以下の乖離分析の詳細をここに置く。いずれも結果に影響し得るため、是正は AI制御ルール（影響度明示のうえ事前承認）に従い 4-2 の各スクリプト単独実行と一体で行う。現結果保持の圧力は外れており、結果変化を受け入れる前提。

### 5-1. MUREN 参照選択（`norm_improved.R`）

- オリジナル MUREN の `saturated`（既定＝全サンプルを参照、自己を含む `1:n_exp`）に対し、`norm_improved.R` は既定で自己参照を除外（`include_self=FALSE` → `setdiff(seq_len(n_exp), k)`）し、さらに中央値近傍へ参照を絞る `refs_cap`（既定 Inf）を新設。これは丸め差でなく、正規化係数推定に入る参照集合の定義変更。
- オリジナルのドキュメントは「all samples as references」と明記し、コードに自己除外案のコメントアウト痕があるため、自己除外はオリジナルの設計選択に反する可能性が高い。ただし MUREN 論文（NAR）の参照選択定義に照らして「違反」か「許容される別実装」かは未判定。オリジナルは非登録の野良スクリプトのため、スクリプト準拠と論文準拠のどちらを基準とするかにも議論の余地。
- 未決: MUREN 論文の参照選択・自己参照の定義確認。是正方針（オリジナル準拠に戻すか、変更を論文的に正当化するか）。BLAS 系統（5-4）とは別系統。contamDE の MUREN 置換（5-2(1)）の妥当性にも波及する上流課題。

### 5-2. contamDE（`contamde_purity_functions.R`）

- **(1) MUREN 正規化への置換**: オリジナルの `limma_voom` 内部 size factor を MUREN 係数に差し替え（`voom(normalize.method="none")`）。パイプライン全体で MUREN を正規化標準とする方針との一貫性から、適切な改変として論文中に明示する方針。ただし妥当性は 5-1 の帰趨に依存。
- **(2) qvalue の残置（BH 是正が正当）**: 純度推定の p 値調整に `qvalue(pi0.method="bootstrap")` を使用。オリジナル `contamDE.lm` はここで `p.adjust(method="fdr")`＝BH を使用しており、かつプロジェクト全体の方針も qvalue から BH へ戻した経緯がある。BH 是正はオリジナル忠実性・プロジェクト方針の二重に正当（[D-017](decisions.md#d-017-contamde-の-p-値調整を-qvaluebh-是正)）。この箇所は top1000 遺伝子ゲートの発動条件（有意遺伝子数）に効き、純度推定 w_hat に伝播するため優先度が高い。

### 5-3. gls.cpp（`contamde_gls.cpp`）

- オリジナルの t 検定ループ（逆行列経由 `solve(t(x_w) %*% W.g %*% x_w)`）を、Cholesky 分解＋三角ソルブ＋`inv_sympd` に置換した C++ 実装。実行者の任意実行。
- 過去の `contamde_purity_functions.R` は gls.cpp 不使用時にフォールバックし、gls 非実行時に結果が変化（カットオフ周辺の精度に影響）した。**現在はフォールバックを廃止し、gls.cpp 未実行ならエラーで停止。現在の結果は gls.cpp に準じたもの。**
- 解法選択（逆行列経由より数値的に安定）は妥当な方向で、数値線形代数の定石に沿う。gls.cpp は RcppArmadillo/C++（R 4.6 の C++20 化・非API除去の対象）と BLAS 依存（Cholesky・inv_sympd）の交差点にあり、R 4.5.3 選定（4.6 回避）を裏側から補強する。
- 未決: gls 採用が結果に与えた影響の最終評価。RcppArmadillo 廃止（5-4）を選ぶ場合、この C++ 実装の扱いも連動。

### 5-4. BLAS / RcppArmadillo

- **BLAS スレッド1固定**（[D-016](decisions.md#d-016-blas-スレッド1固定)）: `norm_improved.R`・05・06・07 で `blas_set_num_threads(1L)`／`OPENBLAS_NUM_THREADS=1` を設定。オリジナル MUREN は BLAS 無指定（環境委任＝R 標準作法）。スレッド1固定は、マルチスレッド BLAS の非決定的丸め順序を排除する再現性のための正当な措置（オリジナルが規定しない層の明示的固定）。実装済み。
- **BLAS 実装差**: OpenBLAS 化により結果は元実装と異なり結論にも影響が及ぶ（AMD64/ARM で結果が変わるのと同種の浮動小数点差。アルゴリズム改変ではない）。最終結論への影響は軽微との記憶（ユーザー）、見直しに前向き。
- **RcppArmadillo/OpenBLAS 廃止（未決）**: 廃止は数値結果に影響する独立した重い設計課題。廃止は問題を別 BLAS 実装へ移すだけになり得るため、廃止でなく「本番環境固定（イメージ凍結）＋結論の頑健性提示」が筋との整理。使用箇所は、RcppArmadillo が 05 の CDM に限定（sourceCpp 1 箇所、呼び出し 1 箇所）、OpenBLAS はパイプライン全体の数値基盤（05/06/07/norm_improved に散在）。廃止するか否かは未決（考えた末に廃止しない可能性も残す）。
