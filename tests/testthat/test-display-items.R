.repo_root <- normalizePath(
  file.path(testthat::test_path(), "..", ".."),
  mustWork = TRUE
)

# Display items (figures, tables, supplementary data): one script per item in
# figures/ or tables/, named by its display ID, reading only declared inputs
# and writing the ID-named file under output/. The mapping is derived from the
# file names and the header contract (first line "# <file>  (<item>)", then
# "# Input :" and "# Output:" lines); there is no separate manifest to keep in
# step. Computation stages in scripts/ write processed/ only.

display_id <- function(file) {
  sub("^((figure|table)_(\\d+|s\\d+)|supplementary_data_\\d+)_.*\\.R$", "\\1",
      basename(file))
}

header_of <- function(file) {
  lines <- readLines(file, warn = FALSE)
  body_start <- grep("^source\\(\"setup.R\"\\)", lines)[1]
  header <- lines[seq_len(if (is.na(body_start)) length(lines) else body_start - 1)]
  header[grepl("^#", header)]
}

test_that("every display-item script is named by its display ID and is unique", {
  figs <- list.files(file.path(.repo_root, "figures"), "\\.R$", full.names = TRUE)
  tabs <- list.files(file.path(.repo_root, "tables"), "\\.R$", full.names = TRUE)
  expect_gt(length(figs), 0)
  expect_gt(length(tabs), 0)
  expect_true(all(grepl("^figure_(\\d+|s\\d+)_[a-z0-9_]+\\.R$", basename(figs))))
  expect_true(all(grepl("^(table_(\\d+|s\\d+)|supplementary_data_\\d+)_[a-z0-9_]+\\.R$",
                        basename(tabs))))
  ids <- display_id(c(figs, tabs))
  expect_false(anyDuplicated(ids) > 0)
})

test_that("each display-item header declares its inputs and its ID-named output", {
  files <- c(
    list.files(file.path(.repo_root, "figures"), "\\.R$", full.names = TRUE),
    list.files(file.path(.repo_root, "tables"), "\\.R$", full.names = TRUE)
  )
  for (f in files) {
    h <- header_of(f)
    id <- display_id(f)
    expect_true(startsWith(h[1], paste0("# ", basename(f), "  (")), info = f)
    expect_true(any(grepl("^# Input :", h)), info = f)
    output_line <- grep("^# Output:", h, value = TRUE)
    expect_length(output_line, 1)
    expected <- if (startsWith(id, "figure_")) {
      sprintf("output/figures/%s.png", id)
    } else {
      sprintf("output/tables/%s.csv", id)
    }
    expect_true(grepl(expected, output_line, fixed = TRUE), info = f)
    inputs <- unlist(regmatches(
      h, gregexpr("(processed|raw|docker)/[A-Za-z0-9_./-]+\\.(rds|csv|tsv)", h)
    ))
    expect_gt(length(inputs), 0)
    missing <- inputs[!file.exists(file.path(.repo_root, inputs))]
    expect_length(missing, 0)
  }
})

test_that("computation stages write processed/ only and render scripts compute nothing", {
  stages <- list.files(file.path(.repo_root, "scripts"), "\\.R$", full.names = TRUE)
  stage_code <- unlist(lapply(stages, function(f) {
    l <- readLines(f, warn = FALSE)
    l[!grepl("^\\s*#", l)]
  }))
  expect_false(any(grepl("paths\\$output", stage_code)))
  renders <- c(
    list.files(file.path(.repo_root, "figures"), "\\.R$", full.names = TRUE),
    list.files(file.path(.repo_root, "tables"), "\\.R$", full.names = TRUE)
  )
  render_code <- unlist(lapply(renders, function(f) {
    l <- readLines(f, warn = FALSE)
    l[!grepl("^\\s*#", l)]
  }))
  expect_false(any(grepl("saveRDS\\(", render_code)))
})
