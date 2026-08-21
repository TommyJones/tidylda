#!/usr/bin/env Rscript
# build-corpora.R --- build and persist the two benchmark corpora. Run once.
#
#   Rscript warp-planning/benchmarks/build-corpora.R [--force]
#
# The medium corpus is a fixed sample drawn under a fixed seed and PERSISTED,
# never resampled. Roadmap 6.1 is emphatic about this: if a later session
# resamples, its numbers are not comparable to the stored baseline and Phase 1
# has to be redone. Hence the refusal to overwrite without --force.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidytext)
})

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "bench-lib.R"))

force <- "--force" %in% commandArgs(trailingOnly = TRUE)
data_dir <- file.path(here, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

repo <- normalizePath(file.path(here, "..", ".."))

guard <- function(path) {
  if (file.exists(path) && !force) {
    stop("Refusing to overwrite ", basename(path), ".\n",
         "  This corpus is the fixed basis of the stored baseline. Resampling it\n",
         "  invalidates every comparison made against that baseline.\n",
         "  Pass --force only if you intend to redo Phase 1.",
         call. = FALSE)
  }
}


# ---- small: nih_sample_dtm, ships with the package --------------------------

small_path <- file.path(data_dir, "small.rds")
guard(small_path)

small_env <- new.env()
load(file.path(repo, "data", "nih_sample_dtm.rda"), envir = small_env)
small <- prune_dtm(get("nih_sample_dtm", envir = small_env))

saveRDS(small, small_path)


# ---- medium: fixed 1,000-document sample of nih ------------------------------
#
# `nih` is 68,508 x 44 and .Rbuildignore'd, so it is available for benchmarking
# and never shipped. Tokenization follows tests/testthat/test-utils.R:7-11 --- the
# unnest_tokens call is the reusable part; that test casts with cast_dtm, we use
# cast_sparse to land directly on a dgCMatrix.
#
# Note: slice_sample(n = 1000) yields 997 rows in the DTM. Three sampled
# abstracts contribute no tokens and drop out at cast time.

medium_path <- file.path(data_dir, "nih-1000.rds")
guard(medium_path)

nih_env <- new.env()
load(file.path(repo, "data", "nih.rda"), envir = nih_env)

set.seed(8675309)
sampled <- dplyr::slice_sample(get("nih", envir = nih_env), n = 1000)

medium <- sampled[, c("APPLICATION_ID", "ABSTRACT_TEXT")] |>
  tidytext::unnest_tokens(input = "ABSTRACT_TEXT", output = "word") |>
  dplyr::count(APPLICATION_ID, word) |>
  tidytext::cast_sparse(APPLICATION_ID, word, n) |>
  prune_dtm()

saveRDS(medium, medium_path)


# ---- manifest ----------------------------------------------------------------

manifest <- dplyr::bind_rows(
  dplyr::bind_cols(tibble::tibble(corpus = "small"),  corpus_stats(small)),
  dplyr::bind_cols(tibble::tibble(corpus = "medium"), corpus_stats(medium))
) |>
  dplyr::mutate(
    sample_seed = c(NA_integer_, 8675309L),
    min_df      = 5L,
    built_at    = Sys.time(),
    git_head    = git_head(repo)
  )

saveRDS(manifest, file.path(data_dir, "corpora-manifest.rds"))

cat("\nCorpora built:\n\n")
print(as.data.frame(manifest[, c("corpus", "D", "V", "N", "nnz")]))
cat("\nWritten to", data_dir, "\n")
