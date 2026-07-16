# build_xwalk.R
# GDC /files API を使って、downloads 配下の RNA tsv を sample_id/case_id にひも付ける。
# その後、既存の sample マスター(DT)と照合し、クロスウォークを保存する。
# 依存: data.table, httr, jsonlite

# ==== 0) セットアップ ==========================================================
if (!requireNamespace("data.table", quietly=TRUE)) stop("Need data.table")
if (!requireNamespace("httr",       quietly=TRUE)) stop("Need httr")
if (!requireNamespace("jsonlite",   quietly=TRUE)) stop("Need jsonlite")
library(data.table); library(httr); library(jsonlite)

options(stringsAsFactors = FALSE)

# ---- パラメータ（環境に合わせて必要なら修正） ----
dir_in            <- "./data/raw/downloads"
pattern_tsv       <- "rna_seq\\.augmented_star_gene_counts\\.tsv(\\.gz)?$"
sample_master_tsv <- "./data/raw/sample.tsv"      # 既存のマスター（DTの元）
dr_tag            <- "DR40"                       # 出力ファイル名に付与
out_dir           <- file.path("outputs","xwalk")
batch_size        <- 100                          # /files POST の1バッチ

# 似たハイフン等の正規化（ID照合用）
norm <- function(x){
  x <- trimws(tolower(as.character(x)))
  gsub("[\u2010-\u2015]", "-", x)  # ハイフン類をASCIIに
}

# ==== 1) downloads のファイル列挙 =============================================
files <- data.table(file_name = list.files(
  dir_in, pattern = pattern_tsv, full.names = FALSE
))
if (nrow(files) == 0L) stop("downloads に対象ファイルが見つかりません。pattern を確認してください。")

# バッチ分割
batches <- split(files$file_name, ceiling(seq_along(files$file_name) / batch_size))

# ==== 2) GDC /files API: file_name → sample_id / aliquot_id / case_id =========
pull_by_files_post <- function(file_names){
  endpoint <- "https://api.gdc.cancer.gov/files"
  payload <- list(
    filters = list(op="in", content=list(field="file_name", value=file_names)),
    fields  = paste(c(
      "file_id","file_name",
      "cases.case_id","cases.submitter_id",
      "cases.samples.sample_id","cases.samples.submitter_id",
      "cases.samples.portions.analytes.aliquots.aliquot_id"
    ), collapse=","),
    format  = "TSV",
    size    = 10000
  )
  res <- httr::POST(
    endpoint,
    body = jsonlite::toJSON(payload, auto_unbox = TRUE),
    httr::add_headers(`Content-Type` = "application/json"),
    encode = "raw"
  )
  httr::stop_for_status(res)
  txt <- httr::content(res, as="text", encoding="UTF-8")
  if (!nzchar(trimws(txt)) || grepl("^\\s*\\{", txt)) return(data.table())
  data.table::fread(text = txt, sep = "\t", colClasses = "character", fill = TRUE)
}

# 列名の配列展開ゆれを吸収して (file_name, sample_id, aliquot_id, case_id) に整形
tidy_gdc_files <- function(tab, requested_files) {
  if (!nrow(tab)) return(data.table(file_name=character(), sample_id=character(),
                                    aliquot_id=character(), case_id=character()))
  tab2 <- copy(tab)
  
  samp_cols <- grep("samples(\\.|\\d).*sample_id$",   names(tab2), value = TRUE)
  aliq_cols <- grep("aliquots(\\.|\\d).*aliquot_id$", names(tab2), value = TRUE)
  case_cols <- grep("(?:^|\\.)cases?(?:\\.|\\d)*case_id$", names(tab2), value = TRUE)
  
  if (length(samp_cols)) tab2[, sample_id  := do.call(fcoalesce, c(.SD, list(NA_character_))), .SDcols = samp_cols]
  if (length(aliq_cols)) tab2[, aliquot_id := do.call(fcoalesce, c(.SD, list(NA_character_))), .SDcols = aliq_cols]
  if (length(case_cols)) tab2[, case_id    := do.call(fcoalesce, c(.SD, list(NA_character_))), .SDcols = case_cols]
  
  out <- tab2[file_name %chin% requested_files, .(file_name, sample_id, aliquot_id, case_id)]
  unique(out[!is.na(file_name)])
}

# 実行
message(sprintf("Querying GDC /files in %d batches ...", length(batches)))
xlist <- lapply(seq_along(batches), function(i){
  v <- batches[[i]]
  tab <- pull_by_files_post(v)
  tidy_gdc_files(tab, v)
})
xwalk <- unique(rbindlist(xlist, use.names = TRUE, fill = TRUE))
setorder(xwalk, file_name)

# 期待どおりに解決できているか簡易チェック
stopifnot(nrow(xwalk) == nrow(files))

# ==== 3) sample マスター読み込み（DT相当） & 照合 ============================
DT <- data.table::fread(
  sample_master_tsv, sep = "\t",
  na.strings = c("", "NA", "--", "'--"),
  colClasses = "character"
)

# 必要列だけ（存在するものだけを採用）
keep_cols <- intersect(c("sample_id","sample_submitter_id","case_id","case_submitter_id","sample_type"), names(DT))
DT  <- DT[, ..keep_cols]
DT[, sample_id := norm(sample_id)]

# xwalk 側も正規化列を作成して照合フラグ
xwalk[, sample_id_norm := norm(sample_id)]
xwalk[, case_id        := ifelse(is.na(case_id), NA_character_, norm(case_id))]
setindex(DT, sample_id)
xwalk[, in_DT := !is.na(sample_id_norm) & (sample_id_norm %chin% DT$sample_id)]

# ==== 4) スモークテスト & duplicate チェック ================================
smoke <- xwalk[, .(
  total_files          = .N,
  unique_sample_ids    = uniqueN(sample_id_norm),
  files_without_id     = sum(is.na(sample_id_norm)),
  present_in_DT        = sum(in_DT)
)]
print(smoke)

dups <- xwalk[, .N, by = sample_id_norm][N > 1]
if (nrow(dups)) {
  message("Duplicate sample_id detected: showing top 10")
  print(dups[order(-N)][1:10])
} else {
  message("duplicate sample_id: 0")
}

# ==== 5) 代表ファイル（_merged 優先）に縮約（※今回一意なら xwalk=canonical） ====
xwalk[, is_merged := grepl("_merged\\.", file_name)]
canonical <- xwalk[order(-is_merged, file_name)][!duplicated(sample_id_norm)]
# 読み込み用のフルパス
canonical[, file_path := file.path(dir_in, file_name)]

# DTメタの付加（あれば）
canon_meta <- merge(
  canonical[, .(file_name, file_path, sample_id = sample_id_norm, case_id, is_merged)],
  unique(DT), by = "sample_id", all.x = TRUE, sort = FALSE
)

# ==== 6) 保存（DRタグ付きで凍結） ===========================================
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
f_xwalk      <- file.path(out_dir, sprintf("file_aliquot_sample_case_%s.tsv", dr_tag))
f_canonical  <- file.path(out_dir, sprintf("file2sample_canonical_%s.tsv",       dr_tag))
f_canon_meta <- file.path(out_dir, sprintf("file2sample_with_meta_%s.tsv",       dr_tag))
f_missing    <- file.path(out_dir, sprintf("file2sample_missing_in_DT_%s.tsv",   dr_tag))

fwrite(xwalk[, .(file_name, sample_id = sample_id_norm, aliquot_id, case_id, in_DT, is_merged)],
       f_xwalk, sep = "\t", quote = FALSE)
fwrite(canonical[, .(file_name, sample_id = sample_id_norm, case_id, is_merged, file_path)],
       f_canonical, sep = "\t", quote = FALSE)
fwrite(canon_meta, f_canon_meta, sep = "\t", quote = FALSE)

missing <- xwalk[!in_DT, .(file_name, sample_id = sample_id_norm, case_id)]
if (nrow(missing)) fwrite(missing, f_missing, sep = "\t", quote = FALSE)

message("Saved:")
message(" - ", f_xwalk)
message(" - ", f_canonical)
message(" - ", f_canon_meta)
if (nrow(missing)) message(" - ", f_missing)

# ==== 7) 終了メッセージ =======================================================
message("Done. present_in_DT = ", smoke$present_in_DT, " / total_files = ", smoke$total_files)
