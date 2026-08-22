#!/usr/bin/env Rscript
# compare.R --- the merge gate for Phases 2, 3 and 5
#
#   Rscript warp-planning/benchmarks/compare.R <baseline.rds> <new.rds>
#   Rscript warp-planning/benchmarks/compare.R --self-test <baseline.rds>
#
# Applies the one-sided non-inferiority test of roadmap D14 to each metric in
# each cell of the grid, and reports an overall PASS/FAIL.
#
# --self-test exists because Phase 1 has only one arm, so the gate would
# otherwise ship without ever having executed --- and Phase 2 is the worst place
# to discover it is wrong, since we would suspect the new sampler first. It
# synthesizes a second arm from the baseline at known offsets and asserts that
# the gate reaches the known-correct verdict.

suppressPackageStartupMessages({
  library(dplyr)
})

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "bench-lib.R"))

# ---- entry points ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if ("--self-test" %in% args) {

  # Validates the gate's DECISION LOGIC, and nothing else.
  #
  # It deliberately does not use the real baseline. Whether the gate reaches the
  # right verdict and whether a given cell has the power to reach any verdict are
  # separate questions, and running the logic check on real data conflates them:
  # an underpowered cell fails even at zero difference, which looks like a broken
  # gate but is a correctly-behaving one reporting that it cannot tell.
  #
  # So: synthetic arms with a spread small enough that every case is
  # comfortably powered (mdd well under the margin), and offsets chosen to sit
  # unambiguously inside or outside it. Power against the real baseline is
  # reported separately, by report(), from real data.

  set.seed(1)
  n <- 20
  base_mean <- 0.20
  noise_sd <- 0.002          # mdd ~ 0.0016, margin ~ 0.010: amply powered

  synth <- function(shift) {
    tidyr::expand_grid(corpus = "synthetic", k = 10L, seed = seq_len(n)) |>
      dplyr::mutate(
        r2             = stats::rnorm(n, base_mean * (1 + shift), noise_sd),
        coherence_mean = stats::rnorm(n, base_mean * (1 + shift), noise_sd)
      )
  }

  cases <- list(
    `offset  0%  (inside margin)`  = list(shift =  0.00, expect = TRUE),
    `offset -2%  (inside margin)`  = list(shift = -0.02, expect = TRUE),
    `offset -8%  (outside margin)` = list(shift = -0.08, expect = FALSE),
    `offset +8%  (better, allowed)`= list(shift =  0.08, expect = TRUE)
  )

  base <- synth(0)
  ok <- TRUE

  for (nm in names(cases)) {
    cat("\n===", nm, "--- expect", ifelse(cases[[nm]]$expect, "PASS", "FAIL"), "===\n\n")
    cmp <- compare_runs(base, synth(cases[[nm]]$shift))
    got <- report(cmp)

    if (any(cmp$mdd > cmp$margin)) {
      ok <- FALSE
      cat("\n  *** SELF-TEST BROKEN: synthetic case is underpowered ***\n")
    }
    if (!identical(got, cases[[nm]]$expect)) {
      ok <- FALSE
      cat("\n  *** SELF-TEST FAILURE: got", ifelse(got, "PASS", "FAIL"), "***\n")
    }
  }

  cat("\n", if (ok) "Self-test passed." else "SELF-TEST FAILED.", "\n")
  quit(status = if (ok) 0 else 1)
}

if (length(args) < 2) {
  stop("usage: compare.R <baseline.rds> <new.rds>   |   compare.R --self-test <baseline.rds>",
       call. = FALSE)
}

base <- readRDS(args[1])
new  <- readRDS(args[2])

cat(sprintf("baseline: %s (%s, %s)\n", args[1], base$engine, base$git_head))
cat(sprintf("new:      %s (%s, %s)\n\n", args[2], new$engine, new$git_head))

passed <- report(compare_runs(base$metrics, new$metrics))
cat("\nOverall:", if (passed) "PASS" else "FAIL", "\n")
quit(status = if (passed) 0 else 1)
