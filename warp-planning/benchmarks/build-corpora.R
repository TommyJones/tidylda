#!/usr/bin/env Rscript
# build-corpora.R --- build and persist the benchmark corpora.
#
#   Rscript warp-planning/benchmarks/build-corpora.R [--force]
#
# The medium corpus is a fixed sample drawn under a fixed seed and PERSISTED,
# never resampled. Roadmap 6.1 is emphatic about this: if a later session
# resamples, its numbers are not comparable to the stored baseline and Phase 1
# has to be redone. Hence: anything already on disk is left alone unless --force.
#
# Re-running is safe and is how a new corpus gets added --- each one is built
# only if missing.

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

# Each corpus is built only if it is missing. A corpus already on disk is the
# fixed basis of the stored baseline, and resampling it invalidates every
# comparison made against that baseline --- so an existing file is left alone and
# reported, not overwritten.
#
# This was a hard error until Phase 7 added the `large` corpus. Erroring made the
# script un-rerunnable, which is wrong once there is more than one corpus and a
# later phase needs to add another. --force still rebuilds everything, and still
# means "I intend to redo Phase 1".
needs_build <- function(path) {
  if (file.exists(path) && !force) {
    message("  ", basename(path), " exists, leaving it alone")
    return(FALSE)
  }
  TRUE
}


# ---- small: nih_sample_dtm, ships with the package --------------------------

small_path <- file.path(data_dir, "small.rds")

if (needs_build(small_path)) {
  small_env <- new.env()
  load(file.path(repo, "data", "nih_sample_dtm.rda"), envir = small_env)
  saveRDS(prune_dtm(get("nih_sample_dtm", envir = small_env)), small_path)
}

small <- readRDS(small_path)


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

# Loaded unconditionally: `large` below needs it even when `medium` is skipped.
nih_env <- new.env()
load(file.path(repo, "data", "nih.rda"), envir = nih_env)

if (needs_build(medium_path)) {
  set.seed(8675309)
  sampled <- dplyr::slice_sample(get("nih", envir = nih_env), n = 1000)

  sampled[, c("APPLICATION_ID", "ABSTRACT_TEXT")] |>
    tidytext::unnest_tokens(input = "ABSTRACT_TEXT", output = "word") |>
    dplyr::count(APPLICATION_ID, word) |>
    tidytext::cast_sparse(APPLICATION_ID, word, n) |>
    prune_dtm() |>
    saveRDS(medium_path)
}

medium <- readRDS(medium_path)


# ---- large: all 68,508 nih abstracts, min_df = 2 ------------------------------
#
# Built for Phase 7, and used ONLY for profiling and for reporting speedup ---
# never for the gate. There is no Gibbs baseline at this scale and generating one
# is not affordable, so statistical parity stays on small/medium against
# baseline-5abaa96.rds, exactly as Phases 2 through 5 ran it.
#
# WHY min_df = 2 RATHER THAN 5. The phase attacks an O(V*K) term, so the corpus
# has to have a V large enough for that term to be visible against O(N). At
# min_df = 5 the full corpus gives V ~ 43k against N ~ 14M, so V*K only overtakes
# N at K ~ 350. At min_df = 2 it is V ~ 80k, and V*K/N reaches 2.8 at K = 500 and
# 5.5 at K = 1000 --- the regime the design notes describe as "would hurt".
#
# min_df = 2 drops only terms appearing in exactly one document, which is mostly
# typos and one-off proper nouns. min_df = 1 would push V to ~164k and reach the
# design notes' 20x case, but keeping hapax legomena is an unusual modeling
# choice and would flatter the optimization: every such term is a word with
# n_w = 1, which is the case the O(1) draw handles best.
#
# No sampling here, so no seed: this is the whole corpus.

large_path <- file.path(data_dir, "nih-full.rds")

if (needs_build(large_path)) {
  get("nih", envir = nih_env)[, c("APPLICATION_ID", "ABSTRACT_TEXT")] |>
    tidytext::unnest_tokens(input = "ABSTRACT_TEXT", output = "word") |>
    dplyr::count(APPLICATION_ID, word) |>
    tidytext::cast_sparse(APPLICATION_ID, word, n) |>
    prune_dtm(min_df = 2) |>
    saveRDS(large_path)
}

large <- readRDS(large_path)


# ---- manifest ----------------------------------------------------------------

manifest <- dplyr::bind_rows(
  dplyr::bind_cols(tibble::tibble(corpus = "small"),  corpus_stats(small)),
  dplyr::bind_cols(tibble::tibble(corpus = "medium"), corpus_stats(medium)),
  dplyr::bind_cols(tibble::tibble(corpus = "large"),  corpus_stats(large))
) |>
  dplyr::mutate(
    sample_seed = c(NA_integer_, 8675309L, NA_integer_),
    min_df      = c(5L, 5L, 2L),
    built_at    = Sys.time(),
    git_head    = git_head(repo)
  )

saveRDS(manifest, file.path(data_dir, "corpora-manifest.rds"))

cat("\nCorpora built:\n\n")
print(as.data.frame(manifest[, c("corpus", "D", "V", "N", "nnz")]))
cat("\nWritten to", data_dir, "\n")
