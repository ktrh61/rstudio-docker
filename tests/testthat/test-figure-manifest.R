.repo_root <- normalizePath(
  file.path(testthat::test_path(), "..", ".."),
  mustWork = TRUE
)

# Manifest consistency (reorg plan v2 s2.6): every figures/ and tables/ script
# appears exactly once, non-empty outputs are unique, and every source object
# exists under processed/. Console-only scripts carry an empty output.

test_that("figure/table manifest is consistent with the scripts on disk", {
  manifest <- read.csv(
    file.path(.repo_root, "figures", "manifest.csv"),
    stringsAsFactors = FALSE
  )
  expect_true(all(
    c(
      "display_id", "script", "output", "status", "source_objects",
      "manuscript_location"
    ) %in% names(manifest)
  ))

  on_disk <- c(
    file.path("figures", list.files(
      file.path(.repo_root, "figures"), pattern = "\\.R$"
    )),
    file.path("tables", list.files(
      file.path(.repo_root, "tables"), pattern = "\\.R$"
    ))
  )
  expect_setequal(manifest$script, on_disk)
  expect_false(anyDuplicated(manifest$script) > 0)

  outputs <- manifest$output[!is.na(manifest$output) & manifest$output != ""]
  expect_false(anyDuplicated(outputs) > 0)

  sources <- unique(unlist(strsplit(manifest$source_objects, ";", fixed = TRUE)))
  missing <- sources[
    !file.exists(file.path(.repo_root, "processed", sources))
  ]
  expect_length(missing, 0L)
})
