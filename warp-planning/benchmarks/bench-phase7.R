#!/usr/bin/env Rscript
# bench-phase7.R --- per-iteration cost, one installed build, across K
#
#   Rscript warp-planning/benchmarks/bench-phase7.R <libpath> <label> <corpus>
#
# Phase 7 replaced the per-word alias table build in the word pass with an O(1)
# mixture draw. This times one build of the package so that two builds can be
# compared: the question is not only whether large-K models get faster but
# whether small-K models get SLOWER, since the change trades O(K) setup per word
# type for a scattered read per token.
#
# BUILD WITH -O2, NOT devtools. pkgbuild::compile_dll() defaults to debug = TRUE
# and compiles at -O0, and devtools::load_all() inherits that. Timings taken that
# way describe unoptimized code and systematically overstate the cost of simple
# K-loops, which are exactly what this change removes and exactly what -O2
# vectorizes. Install with R CMD INSTALL into a private library and load from
# there, which is what a user actually runs.
#
# Per-iteration cost is a difference between two fits at different iteration
# counts, with a wide spread; see profile-vk.R's header for why both the
# differencing and the width are necessary.

args   <- commandArgs(trailingOnly = TRUE)
libp   <- args[1]
label  <- args[2]
corpus <- if (length(args) >= 3) args[3] else "large"

suppressPackageStartupMessages(library(tidylda, lib.loc = libp))

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
data_dir <- file.path(here, "data")
paths <- c(small = "small.rds", medium = "nih-1000.rds", large = "nih-full.rds")
dtm <- readRDS(file.path(data_dir, paths[[corpus]]))

V <- ncol(dtm)
iters <- c(2L, 62L)
ks <- if (corpus == "large") c(10L, 50L, 200L, 500L, 1000L) else c(10L, 50L, 200L)

secs <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

for (k in ks) {
  alpha <- rep(0.1, k)
  init <- tidylda:::initialize_topic_counts(
    dtm = dtm, k = k, alpha = alpha,
    eta = matrix(0.05, nrow = k, ncol = V), threads = 12L
  )

  run <- function(n) {
    t0 <- Sys.time()
    invisible(tidylda:::fit_lda_warp(
      dtm_in = dtm, Cd_start = init$Cd_start, alpha_in = alpha,
      eta_in = matrix(0.05, 1, 1), iterations = n, burnin = -1L,
      calc_likelihood = FALSE, Beta_in = init$beta_initial,
      freeze_topics = FALSE, threads = 12L, verbose = FALSE
    ))
    secs(t0)
  }

  lo <- run(iters[1]); hi <- run(iters[2])
  cat(sprintf("RESULT %s %s k=%d per_iter=%.4f init=%.1f\n",
              label, corpus, k, (hi - lo) / (iters[2] - iters[1]),
              lo - (hi - lo) / (iters[2] - iters[1]) * iters[1]))
}
