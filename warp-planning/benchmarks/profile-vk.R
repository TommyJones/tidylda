#!/usr/bin/env Rscript
# profile-vk.R --- is the O(V*K) term real, and how big is it?
#
#   Rscript warp-planning/benchmarks/profile-vk.R [corpus] [threads]
#
# Phase 7 attacks per-word-type work in the word pass: the sampler does O(K)
# work for every word type every iteration, regardless of how many tokens that
# word has. The roadmap is explicit that this phase needs a measurement before it
# needs an implementation, and this is that measurement.
#
# METHOD, AND THE TWO WAYS IT GOES WRONG.
#
# Per-iteration cost is obtained by DIFFERENCING two fits at different iteration
# counts, never by dividing one fit's total by its iterations. Initialization at
# this scale is minutes, and folding it into a per-iteration number is exactly
# the error that made Phase 5 report an efficiency collapse that turned out not
# to exist.
#
# But differencing has its own failure mode, and the first version of this script
# hit it: if both runs are dominated by a large fixed cost, their difference is a
# small number between two big ones and timer noise swamps it. That version
# timed tidylda() end to end at 3 and 8 iterations --- roughly 220 s apiece at
# K = 1000, differing by 3 s --- and produced a NEGATIVE slope at p = 0.39.
#
# Two fixes, both necessary:
#
#   1. Call fit_lda_warp() directly, so the R-side Dirichlet draws and all of
#      tidylda()'s post-processing sit outside the timed region.
#   2. Widen the iteration spread to 20, so the difference is tens of seconds
#      rather than single digits.
#
# The engine's own fused initialization (D16) is still inside both timed runs and
# still cancels in the difference. It is reported separately because at K = 1000
# it is ~107 s, which is worth knowing about on its own.

suppressPackageStartupMessages({
  library(dplyr)
  devtools::load_all(quiet = TRUE)
})

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "bench-lib.R"))

args    <- commandArgs(trailingOnly = TRUE)
corpus  <- if (length(args) >= 1) args[1] else "large"
threads <- if (length(args) >= 2) as.integer(args[2]) else 12L

data_dir <- file.path(here, "data")
paths <- c(small = "small.rds", medium = "nih-1000.rds", large = "nih-full.rds")
dtm <- readRDS(file.path(data_dir, paths[[corpus]]))

D <- nrow(dtm); V <- ncol(dtm); N <- sum(dtm)

cat(sprintf("corpus %s: D=%d V=%d N=%d   threads=%d\n\n", corpus, D, V, N, threads))

# SPREAD IS PRECISION. per_iter is a difference between two timed runs, so its
# error is roughly (noise in each run) / spread. Wall-clock noise here is 1-3 s
# --- turbo, scheduling, page faults on multi-GB allocations --- so a spread of
# 20 leaves +/- 0.1-0.3 s on a ~1.5 s quantity, which is 10-20%.
#
# That is not academic. At spread 20 the two accumulate terms measured NEGATIVE
# at K = 500 and K = 1000, i.e. adding work made the sampler faster. A spread of
# 60 puts the difference at ~100 s and the error near 4%.
iters <- c(2L, 62L)
ks    <- c(50L, 200L, 500L, 1000L, 2000L)

# BURN-IN IS A DIMENSION, NOT A DETAIL. Two of the five O(K)-per-work-item terms
# --- the Cd_sum and Cv_sum accumulations --- are gated on `averaging`, which is
# `burnin > -1` (warp_lda.cpp:323). Profiling with burnin = -1 silently omits
# them, which the first run of this script did.
#
# They are the expensive two. At K = 1000 on this corpus Cv_sum is V*K longs
# (649 MB) and Cd_sum is D*K longs (546 MB), both read-modify-written every
# post-burn-in iteration on top of reading Cd and Cv themselves. A real fit
# spends most of its iterations in that state, so the burnin = 1 row is the one
# that describes production, and the burnin = -1 row is the floor.
burnins <- c(-1L, 1L)

secs <- function(t0) as.numeric(difftime(Sys.time(), t0, units = "secs"))

measure_k <- function(k) {
  alpha <- rep(0.1, k)

  # Once per K, outside every timed region and shared across burn-in settings.
  t0 <- Sys.time()
  init <- tidylda:::initialize_topic_counts(
    dtm = dtm, k = k, alpha = alpha,
    eta = matrix(0.05, nrow = k, ncol = V), threads = threads
  )
  r_init <- secs(t0)

  run <- function(n_iter, burnin) {
    t0 <- Sys.time()
    invisible(fit_lda_warp(
      dtm_in = dtm, Cd_start = init$Cd_start, alpha_in = alpha,
      eta_in = matrix(0.05, 1, 1),          # D20 scalar path
      iterations = n_iter, burnin = burnin, calc_likelihood = FALSE,
      Beta_in = init$beta_initial, freeze_topics = FALSE,
      threads = threads, verbose = FALSE
    ))
    secs(t0)
  }

  dplyr::bind_rows(lapply(burnins, function(b) {
    t_lo <- run(iters[1], b)
    t_hi <- run(iters[2], b)

    per_iter    <- (t_hi - t_lo) / (iters[2] - iters[1])
    engine_init <- t_lo - per_iter * iters[1]

    cat(sprintf(
      "  K=%4d burnin=%2d  per-iteration %7.3f s | engine init %7.1f s | VK/N %5.2f\n",
      k, b, per_iter, engine_init, V * k / N
    ))

    tibble::tibble(k = k, burnin = b, per_iter = per_iter,
                   engine_init = engine_init, r_init = r_init,
                   VK = as.numeric(V) * k)
  }))
}

res <- dplyr::bind_rows(lapply(ks, measure_k))

# ---- is the VK term real? ----------------------------------------------------
#
# N is constant across the sweep, so the O(N) term is absorbed into the intercept
# and the question reduces to whether per-iteration cost rises with V*K.

fits <- lapply(burnins, function(b) stats::lm(per_iter ~ VK, data = dplyr::filter(res, burnin == b)))
names(fits) <- as.character(burnins)

for (b in burnins) {
  f <- fits[[as.character(b)]]; co <- summary(f)$coefficients
  cat(sprintf("\nburnin=%2d:  per_iter = %.4g + %.4g * VK   (R^2 %.4f, slope p %.3g)\n",
              b, co[1, 1], co[2, 1], summary(f)$r.squared, co[2, 4]))

  dplyr::filter(res, burnin == b) |>
    dplyr::transmute(
      k, per_iter = round(per_iter, 3),
      vk_secs  = round(co[2, 1] * VK, 3),
      vk_share = sprintf("%5.1f%%", 100 * co[2, 1] * VK / per_iter),
      ceiling  = sprintf("%.2fx", 1 / (1 - co[2, 1] * VK / per_iter))
    ) |>
    as.data.frame() |>
    print(row.names = FALSE)
}

cat("\ncost of the two accumulate terms alone (burnin=1 minus burnin=-1):\n")
res |>
  dplyr::select(k, burnin, per_iter) |>
  tidyr::pivot_wider(names_from = burnin, values_from = per_iter,
                     names_prefix = "b") |>
  dplyr::transmute(
    k,
    no_accum = round(`b-1`, 3), with_accum = round(b1, 3),
    accum_secs = round(b1 - `b-1`, 3),
    accum_share = sprintf("%5.1f%%", 100 * (b1 - `b-1`) / b1)
  ) |>
  as.data.frame() |>
  print(row.names = FALSE)

saveRDS(
  list(corpus = corpus, D = D, V = V, N = N, threads = threads,
       iters = iters, burnins = burnins, results = res, fits = fits,
       git_head = git_head(normalizePath(file.path(here, "..", "..")))),
  file.path(data_dir, sprintf("profile-vk-%s.rds", corpus))
)

cat(sprintf("\nwritten to profile-vk-%s.rds\n", corpus))
