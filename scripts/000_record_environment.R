# 000_record_environment.R
# Record the execution environment of this run and the R packages the code
# tree loads, so that the software table (Table S2) is rendered from
# processed/ like every other display item and reflects the environment that
# actually ran the pipeline.
# Input : the code tree (scripts/, lib/, figures/, tables/, tests/testthat,
#         config.R, setup.R) -- scanned for library(), requireNamespace() and
#         pkg:: calls; comments do not count as use
#         docker/versions.tsv (pinned versions of the qualified container;
#         used only as a run-time guard, not as the recorded value)
# Output: processed/thyr_environment.rds
#           r_version, platform, os, blas, lapack, run_time, git_commit,
#           packages (package, version, role), files_scanned
#
# The recorded versions are the live values (packageVersion()) of the session
# that runs this script. The run stops if any loaded package differs from the
# pinned version, or has no pinned version, so that a run outside the
# qualified environment cannot silently produce the reported artifacts.
# Roles (first match in this order): pipeline = scripts/, lib/, config.R,
# setup.R; render = figures/, tables/; tests = tests/testthat.

source("setup.R")

# --- Code tree to scan -----------------------------------------------------
role_files <- list(
  pipeline = c(
    list.files(file.path(paths$root, "scripts"), "\\.R$", full.names = TRUE),
    list.files(file.path(paths$root, "lib"), "\\.R$", full.names = TRUE),
    file.path(paths$root, c("config.R", "setup.R"))
  ),
  render = c(
    list.files(file.path(paths$root, "figures"), "\\.R$", full.names = TRUE),
    list.files(file.path(paths$root, "tables"), "\\.R$", full.names = TRUE)
  ),
  tests = list.files(file.path(paths$root, "tests", "testthat"), "\\.R$",
                     full.names = TRUE)
)

packages_used_in <- function(files) {
  txt <- unlist(lapply(files, readLines, warn = FALSE))
  txt <- txt[!grepl("^\\s*#", txt)]
  hits <- c(
    regmatches(txt, gregexpr("library\\(([A-Za-z][A-Za-z0-9.]*)\\)", txt)),
    regmatches(txt, gregexpr("requireNamespace\\(\"([A-Za-z][A-Za-z0-9.]*)\"", txt)),
    regmatches(txt, gregexpr("\\b([A-Za-z][A-Za-z0-9.]*)::", txt))
  )
  hits <- unlist(hits)
  unique(gsub("^library\\(|^requireNamespace\\(\"|\\)$|\"$|::$", "", hits))
}

base_pkgs <- c(rownames(utils::installed.packages(priority = "base")), "R")
rows <- do.call(rbind, lapply(names(role_files), function(r) {
  p <- setdiff(packages_used_in(role_files[[r]]), base_pkgs)
  if (length(p) == 0) return(NULL)
  data.frame(package = p, role = r, stringsAsFactors = FALSE)
}))
rows$role <- factor(rows$role, levels = names(role_files))
rows <- rows[order(rows$role, tolower(rows$package)), ]
rows <- rows[!duplicated(rows$package), ]
rows$role <- as.character(rows$role)

# --- Live versions and the pinned-version guard ----------------------------
rows$version <- vapply(rows$package, function(p) {
  as.character(utils::packageVersion(p))
}, character(1))

pin_lines <- readLines(file.path(paths$root, "docker", "versions.tsv"))
pin_lines <- pin_lines[!grepl("^#", pin_lines) & nzchar(pin_lines)]
pinned <- utils::read.delim(text = paste(pin_lines, collapse = "\n"),
                            stringsAsFactors = FALSE)
missing <- setdiff(rows$package, pinned$package)
if (length(missing) > 0) {
  stop("Loaded packages without a pinned version in docker/versions.tsv: ",
       paste(missing, collapse = ", "))
}
want <- gsub("-", ".", pinned$version[match(rows$package, pinned$package)],
             fixed = TRUE)
bad <- rows$package[rows$version != want]
if (length(bad) > 0) {
  stop("Package versions differ from docker/versions.tsv: ",
       paste(sprintf("%s (%s != %s)", bad, rows$version[rows$package %in% bad],
                     want[rows$package %in% bad]), collapse = ", "))
}

# --- Session facts ---------------------------------------------------------
os_release <- if (file.exists("/etc/os-release")) {
  l <- readLines("/etc/os-release")
  gsub("^PRETTY_NAME=|\"", "", grep("^PRETTY_NAME=", l, value = TRUE))
} else {
  paste(Sys.info()[c("sysname", "release")], collapse = " ")
}
# Commit of the code tree, read from .git directly (the container carries no
# git binary). NA when the tree is not a git checkout.
read_git_commit <- function(root) {
  head_file <- file.path(root, ".git", "HEAD")
  if (!file.exists(head_file)) return(NA_character_)
  head <- readLines(head_file, n = 1, warn = FALSE)
  if (grepl("^[0-9a-f]{40}$", head)) return(head)
  ref <- sub("^ref: ", "", head)
  ref_file <- file.path(root, ".git", ref)
  if (file.exists(ref_file)) return(readLines(ref_file, n = 1, warn = FALSE))
  packed <- file.path(root, ".git", "packed-refs")
  if (file.exists(packed)) {
    hit <- grep(paste0(" ", ref, "$"), readLines(packed, warn = FALSE), value = TRUE)
    if (length(hit) == 1) return(sub(" .*$", "", hit))
  }
  NA_character_
}
git_commit <- read_git_commit(paths$root)

environment_record <- list(
  r_version = R.version.string,
  platform = R.version$platform,
  os = os_release,
  blas = extSoftVersion()[["BLAS"]],
  lapack = La_library(),
  run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  git_commit = git_commit,
  packages = rows[, c("package", "version", "role")],
  files_scanned = sub(paste0("^", paths$root, "/"), "",
                      unlist(role_files, use.names = FALSE))
)
rownames(environment_record$packages) <- NULL

cat(environment_record$r_version, "on", environment_record$os, "\n")
cat("BLAS  :", environment_record$blas, "\nLAPACK:", environment_record$lapack, "\n")
cat("Commit:", environment_record$git_commit, "\n")
print(environment_record$packages, row.names = FALSE)
cat("packages:", nrow(environment_record$packages), "| pinned but not loaded:",
    paste(setdiff(pinned$package, rows$package), collapse = ", "), "\n")

out <- file.path(paths$processed, "thyr_environment.rds")
saveRDS(environment_record, out)
cat("Saved:", out, "\n")
