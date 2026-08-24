# bench-lib.R --- shared helpers for the Phase 1 statistical benchmark
#
# Sourcing this file has no side effects: it defines functions and constants and
# nothing else. The three driver scripts (build-corpora.R, run-benchmark.R,
# compare.R) source it and do the acting.
#
# See warp-planning/warplda-roadmap.md 6.1 for the specification this implements.

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tibble)
})


# ---- constants --------------------------------------------------------------

#' Fitting hyperparameters, pinned per engine and recorded in the results.
#'
#' Priors and metrics are identical across engines. **Iteration counts are not**,
#' because an iteration does not mean the same amount of work to both samplers:
#' collapsed Gibbs draws every token from its exact conditional, while warpLDA
#' takes one Metropolis-Hastings proposal per pass and rejects some fraction.
#' Matching iterations would compare mixing rates; matching *convergence* is what
#' D14 actually asks about.
#'
#' Measured: warp reaches Gibbs' 200-iteration quality at roughly 1200
#' iterations, within +/-5% on both metrics in all four cells (three within
#' +/-2%), with its log-likelihood plateauing at the same level. It is still
#' 1.5x to 7.2x faster in wall clock at that point, the advantage growing in K.
#'
#' `threads = 1` is load-bearing rather than incidental. Design notes 9(b)
#' records that the old sampler's batch loop is broken for `threads > 1` and was
#' deliberately left unfixed, since the rewrite deletes the code it lives in.
#' All parallelism in this harness is across fits, never inside one.
BENCH_HP <- list(
  gibbs = list(
    iterations      = 200,
    burnin          = 50,
    alpha           = 0.1,
    eta             = 0.05,
    calc_likelihood = TRUE,
    calc_r2         = TRUE,
    threads         = 1
  ),
  warp = list(
    iterations       = 1200,
    burnin           = 300,
    alpha            = 0.1,
    eta              = 0.05,
    calc_likelihood  = TRUE,
    calc_r2          = TRUE,
    threads          = 1,
    mh_steps         = 1,      # D18 default; reproduces the reference exactly
    likelihood_every = 10      # D11; a read-only diagnostic, not thinning
  )
)

#' Hyperparameters for one engine
bench_hp <- function(engine) {
  if (is.null(BENCH_HP[[engine]])) stop("No hyperparameters defined for engine '", engine, "'")
  BENCH_HP[[engine]]
}

#' The benchmark grid: two corpora, two topic counts, seeds calibrated per K.
#'
#' Seed counts are not uniform, and the asymmetry is measured rather than
#' guessed. From the Phase 1 Gibbs baseline, the minimum detectable difference
#' against a 5% margin:
#'
#'   metric      K=10                  K=50
#'   R^2         ratio 0.35 / 0.36     ratio 0.20 / 0.25    (n = 3 would do)
#'   coherence   ratio 1.53 / 1.88     ratio 0.74 / 0.77
#'
#' Coherence is the binding metric and it is underpowered only at K=10, where it
#' needs 48 (medium) and 71 (small) seeds. Between-seed spread in coherence
#' roughly doubles from K=50 to K=10 --- CV 0.046/0.048 against 0.096/0.117 ---
#' on both corpora, so this tracks topic count, not corpus size. Misspecified K
#' is the high-variance regime, and because the margin is *relative* and the
#' K=10 coherence mean is also lower, the margin shrinks at exactly the K where
#' the noise grows.
#'
#' 100 seeds at K=10 covers the 71 requirement with headroom, since n_needed is
#' itself estimated from a noisy sd. K=50 keeps 20; raising it would cost ~50
#' minutes of compute for a cell that already has 1.7x the power it needs.
BENCH_CORPORA <- c("small", "medium")
BENCH_KS      <- c(10L, 50L)
# K=50 raised from 20 to 100 after the Phase 4/4.5 gate run: small/K=50
# coherence came back at p = 0.168 with the point estimate at 74% of the margin
# and mdd (0.0063) close to it, which is the harness saying it cannot tell rather
# than that something is wrong. 6.1's prescription is to add seeds, not widen the
# margin, and these are the cheap cells -- 10 s and 104 s a fit.
BENCH_SEEDS   <- list("10" = 1:100, "50" = 1:100)

#' Seeds for a given K
seeds_for_k <- function(k) BENCH_SEEDS[[as.character(k)]]

#' Gate parameters (roadmap D14, as revised for one-sided non-inferiority).
BENCH_MARGIN_FRAC <- 0.05
BENCH_ALPHA       <- 0.05


# ---- corpus construction ----------------------------------------------------

#' Prune a document-term matrix to a fixed vocabulary
#'
#' Removes `tidytext::stop_words`, then drops terms appearing in fewer than
#' `min_df` documents. V --- which the whole O(VK) question depends on --- is
#' fixed by this and recorded in the manifest.
#'
#' `small` and `medium` use the default min_df = 5. `large` uses 2, which
#' roughly doubles V at the same N; see build-corpora.R for why. That makes
#' `large` not directly comparable to the other two, which is fine: it exists
#' to profile the VK term, never to run the gate.
#'
#' @param dtm a dgCMatrix with colnames
#' @param min_df minimum document frequency for a term to be retained
#' @return a dgCMatrix with a subset of the input's columns
prune_dtm <- function(dtm, min_df = 5) {
  stopifnot(!is.null(colnames(dtm)), min_df >= 1)

  keep <- !colnames(dtm) %in% tidytext::stop_words$word
  dtm <- dtm[, keep, drop = FALSE]

  doc_freq <- Matrix::colSums(dtm > 0)
  dtm[, doc_freq >= min_df, drop = FALSE]
}


#' Summarise a corpus in one row
#'
#' D documents, V terms, N total tokens, nnz nonzero cells. Recorded alongside
#' the results so a later session can confirm it is comparing like with like.
corpus_stats <- function(dtm) {
  tibble::tibble(
    D   = nrow(dtm),
    V   = ncol(dtm),
    N   = sum(dtm),
    nnz = length(dtm@x)
  )
}


# ---- the grid ---------------------------------------------------------------

#' Stable identifier for one fit, used as its result filename
#'
#' Seeds are zero-padded to three digits so that lexical order matches numeric
#' order --- 6.1 anticipates adding seeds if 20 proves underpowered.
fit_id <- function(corpus, k, seed) {
  sprintf("%s-k%02d-seed%03d", corpus, k, seed)
}


#' Enumerate the benchmark grid, one row per fit
#'
#' @param seeds `NULL` to use the per-K calibration in `BENCH_SEEDS`, or an
#'   explicit vector to apply the same seeds to every K (used by `--probe` and
#'   by `--seeds=`).
bench_grid <- function(corpora = BENCH_CORPORA,
                       ks = BENCH_KS,
                       seeds = NULL) {
  dplyr::bind_rows(lapply(ks, function(k) {
    tidyr::expand_grid(
      corpus = corpora,
      k = k,
      seed = if (is.null(seeds)) seeds_for_k(k) else seeds
    )
  })) |>
    dplyr::mutate(id = fit_id(.data$corpus, .data$k, .data$seed)) |>
    dplyr::arrange(.data$corpus, .data$k, .data$seed)
}


# ---- fitting ----------------------------------------------------------------

#' Pull the two benchmark metrics out of a fitted tidylda object
#'
#' Both metrics are computed by `tidylda()` itself: `model$r2` comes from
#' `calc_lda_r2()` (R/utils.R:770) and `model$summary$coherence` from
#' `calc_prob_coherence(m = 5)` (R/utils.R:497). We deliberately read them off
#' the model rather than recomputing, so the harness measures exactly what the
#' package reports.
#'
#' Both are wrapped in `try()` inside `new_tidylda()` (R/utils.R:731, :697), so a
#' failure arrives as a `try-error` object rather than as an error. This function
#' turns that into NA plus a status string rather than letting it propagate as a
#' nonsense number.
#'
#' Per-topic coherence can itself be NA (a top-5 term that never co-occurs). We
#' average with `na.rm = TRUE` but record how many topics were dropped, so a
#' quietly-degenerate cell is visible in the results rather than invisible.
extract_metrics <- function(model) {
  status <- character(0)

  r2 <- model$r2
  if (is.null(r2) || inherits(r2, "try-error") || !is.numeric(r2)) {
    r2 <- NA_real_
    status <- c(status, "r2_failed")
  }

  coh_mean <- NA_real_
  coh_n_na <- NA_integer_
  if (is.null(model$summary) || inherits(model$summary, "try-error")) {
    status <- c(status, "summary_failed")
  } else {
    coh <- model$summary$coherence
    coh_n_na <- sum(is.na(coh))
    if (coh_n_na == length(coh)) {
      status <- c(status, "coherence_all_na")
    } else {
      coh_mean <- mean(coh, na.rm = TRUE)
      if (coh_n_na > 0) status <- c(status, "coherence_partial_na")
    }
  }

  tibble::tibble(
    r2              = as.numeric(r2),
    coherence_mean  = coh_mean,
    coherence_n_na  = as.integer(coh_n_na),
    status          = if (length(status)) paste(status, collapse = ",") else "ok"
  )
}


#' Fit one model and return its metrics and likelihood curve
#'
#' Reproducibility is per fit, not per worker: `set.seed(seed)` is called here,
#' inside whatever process ends up running this fit. A fit's result is therefore
#' independent of which core ran it and of how many cores were used --- the
#' harness-level analogue of D12's per-work-item seeding.
#'
#' The likelihood curve is kept because 6.1 requires validating the iteration
#' setting by inspecting it. It is a few kilobytes per fit.
#'
#' Only the engine `tidylda()` currently dispatches to can be fitted here. As of
#' Phase 2 that is warpLDA; the Gibbs arm is the stored baseline from Phase 1
#' (`baseline-5abaa96.rds`) and is not re-runnable from this harness.
#'
#' @return a list with `metrics` (one row) and `log_likelihood` (a tibble)
fit_once <- function(dtm, k, seed, hp = bench_hp("warp")) {
  set.seed(seed)

  args <- list(
    data            = dtm,
    k               = k,
    iterations      = hp$iterations,
    burnin          = hp$burnin,
    alpha           = hp$alpha,
    eta             = hp$eta,
    calc_likelihood = hp$calc_likelihood,
    calc_r2         = hp$calc_r2,
    threads         = hp$threads,
    verbose         = FALSE
  )
  # Engine-specific extras travel through tidylda()'s dots.
  for (nm in c("mh_steps", "likelihood_every")) {
    if (!is.null(hp[[nm]])) args[[nm]] <- hp[[nm]]
  }

  elapsed <- system.time(model <- do.call(tidylda::tidylda, args))

  list(
    metrics = dplyr::bind_cols(
      extract_metrics(model),
      tibble::tibble(elapsed_sec = unname(elapsed[["elapsed"]]))
    ),
    log_likelihood = model$log_likelihood
  )
}


# ---- per-fit files ----------------------------------------------------------
#
# One file per fit, written the moment that fit completes. This is what makes a
# run crash-resumable and what lets the runner skip work already done. The
# schema is defined by `run_grid_row()` and consumed by `collect_results()`;
# both live here so the two cannot drift apart.

#' Run one row of the grid and persist it
#'
#' Returns the fit's id invisibly. Errors are caught and recorded rather than
#' thrown: one pathological cell should not take down a multi-hour grid, and a
#' failed fit is itself a result worth seeing.
run_grid_row <- function(row, dtm, out_dir, hp = bench_hp("warp")) {
  path <- file.path(out_dir, paste0(row$id, ".rds"))

  fit <- try(fit_once(dtm = dtm, k = row$k, seed = row$seed, hp = hp), silent = TRUE)

  if (inherits(fit, "try-error")) {
    result <- list(
      metrics = tibble::tibble(
        r2 = NA_real_, coherence_mean = NA_real_, coherence_n_na = NA_integer_,
        status = paste0("fit_error: ", conditionMessage(attr(fit, "condition"))),
        elapsed_sec = NA_real_
      ),
      log_likelihood = NULL
    )
  } else {
    result <- fit
  }

  result$metrics <- dplyr::bind_cols(
    tibble::tibble(
      id = row$id, corpus = row$corpus, k = row$k, seed = row$seed
    ),
    result$metrics
  )

  saveRDS(result, path)
  invisible(row$id)
}


#' Read every per-fit result file in a directory into one tibble
#'
#' Each fit writes its own file the moment it completes, which is what makes a
#' run crash-resumable. This reassembles them.
#'
#' @return a list with `metrics` (one row per fit) and `curves` (named list of
#'   likelihood tibbles)
collect_results <- function(dir) {
  files <- sort(list.files(dir, pattern = "\\.rds$", full.names = TRUE))
  if (length(files) == 0) {
    stop("No result files found in ", dir)
  }

  fits <- lapply(files, readRDS)

  list(
    metrics = dplyr::bind_rows(lapply(fits, `[[`, "metrics")),
    curves  = stats::setNames(
      lapply(fits, `[[`, "log_likelihood"),
      vapply(fits, function(f) f$metrics$id, character(1))
    )
  )
}


# ---- the gate ---------------------------------------------------------------

#' One-sided non-inferiority test of a new sampler against the baseline
#'
#' Roadmap D14, as revised. With margin delta = `margin_frac` * mean(baseline):
#'
#'   H0: mu_new - mu_base <= -delta      vs      H1: mu_new - mu_base > -delta
#'
#' Rejecting H0 says the new sampler is not worse than the baseline by more than
#' delta. PASS is a rejection at level `alpha`.
#'
#' Two things this deliberately is not:
#'
#'   * Not a paired test. Pairing across seeds is genuine only while both engines
#'     share an initialization path; at Phase 4 initialization moves into C++ and
#'     seed i stops producing a common starting state. An unpaired Welch test is
#'     valid in both regimes, so the gate does not silently change meaning
#'     mid-project.
#'   * Not a superiority test (H0: diff <= 0). A correct port produces a true
#'     difference near zero, which such a test would fail to reject roughly 95%
#'     of the time. The margin is what makes "no worse than" testable.
#'
#' The upper side is computed too, but only as a diagnostic: a new sampler that
#' scores implausibly better than the baseline is more likely to be one that is
#' not really sampling than a free win.
#'
#' @param x_new,x_base numeric vectors of per-seed metric values
#' @return a one-row tibble
noninf_test <- function(x_new, x_base,
                        margin_frac = BENCH_MARGIN_FRAC,
                        alpha = BENCH_ALPHA) {
  x_new  <- x_new[is.finite(x_new)]
  x_base <- x_base[is.finite(x_base)]
  stopifnot(length(x_new) >= 2, length(x_base) >= 2)

  base_mean <- mean(x_base)
  margin <- margin_frac * abs(base_mean)

  lower <- stats::t.test(x_new, x_base, alternative = "greater", mu = -margin)
  upper <- stats::t.test(x_new, x_base, alternative = "less",    mu =  margin)

  tibble::tibble(
    n_new        = length(x_new),
    n_base       = length(x_base),
    mean_new     = mean(x_new),
    mean_base    = base_mean,
    sd_new       = stats::sd(x_new),
    sd_base      = stats::sd(x_base),
    diff         = mean(x_new) - base_mean,
    margin       = margin,
    p_noninf     = lower$p.value,
    ci_lower     = unname(lower$conf.int[1]),
    pass         = lower$p.value < alpha,
    p_not_better = upper$p.value,
    mdd          = min_detectable_diff(x_base, length(x_new), alpha)
  )
}


METRICS <- c(r2 = "r2", coherence = "coherence_mean")


# ---- the gate ----------------------------------------------------------------

#' Apply the gate to every cell x metric of two runs
#'
#' @param base,new metrics tibbles as stored by run-benchmark.R
#' @return one row per cell per metric
compare_runs <- function(base, new,
                         margin_frac = BENCH_MARGIN_FRAC,
                         alpha = BENCH_ALPHA) {
  cells <- base |>
    dplyr::distinct(.data$corpus, .data$k) |>
    dplyr::arrange(.data$corpus, .data$k)

  rows <- list()
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, ]
    pick <- function(df, col) {
      df[[col]][df$corpus == cell$corpus & df$k == cell$k]
    }
    for (nm in names(METRICS)) {
      rows[[length(rows) + 1]] <- dplyr::bind_cols(
        tibble::tibble(corpus = cell$corpus, k = cell$k, metric = nm),
        noninf_test(
          x_new       = pick(new, METRICS[[nm]]),
          x_base      = pick(base, METRICS[[nm]]),
          margin_frac = margin_frac,
          alpha       = alpha
        )
      )
    }
  }
  dplyr::bind_rows(rows)
}


#' Print a comparison table and return the overall verdict
report <- function(cmp) {
  out <- cmp |>
    dplyr::transmute(
      cell     = paste0(.data$corpus, "/k", .data$k),
      metric   = .data$metric,
      base     = round(.data$mean_base, 4),
      new      = round(.data$mean_new, 4),
      diff     = round(.data$diff, 4),
      margin   = round(.data$margin, 4),
      mdd      = round(.data$mdd, 4),
      sd_ratio = round(.data$sd_new / .data$sd_base, 2),
      p        = signif(.data$p_noninf, 3),
      verdict  = ifelse(.data$pass, "PASS", "FAIL"),
      better   = ifelse(.data$p_not_better >= BENCH_ALPHA, "check", "")
    )
  print(as.data.frame(out), row.names = FALSE)

  # The gate compares means. A sampler that mixes worse can match on the mean
  # while being more seed-dependent, and that is most likely at misspecified K,
  # where between-seed spread is already largest (coherence CV roughly doubles
  # from K=50 to K=10 in the Gibbs baseline). An inflated ratio is not a failure
  # --- it is where to look first if anything downstream seems off.
  inflated <- cmp$sd_new / cmp$sd_base > 1.5
  if (any(inflated, na.rm = TRUE)) {
    cat("\nNote: sd_ratio > 1.5 in", sum(inflated, na.rm = TRUE),
        "cell(s) --- the new sampler matches on the mean but is markedly more\n",
        "  seed-dependent. Suspect mixing, not bias.\n")
  }

  underpowered <- cmp$mdd > cmp$margin
  if (any(underpowered)) {
    cat("\nNote: mdd > margin in",  sum(underpowered),
        "cell(s) --- the gate cannot resolve a",
        paste0(BENCH_MARGIN_FRAC * 100, "%"),
        "margin at this sample size.\n",
        "  Roadmap 6.1: add seeds rather than widen the margin.\n")
  }
  if (any(cmp$p_not_better >= BENCH_ALPHA)) {
    cat("\nNote: cells marked 'check' scored better than the baseline by more\n",
        "  than the margin. Not a failure, but a sampler that is not really\n",
        "  sampling can look like this. Worth a second look.\n")
  }

  all(cmp$pass)
}


#' Paired one-sided non-inferiority test
#'
#' Same hypothesis as `noninf_test()` --- H0: mu_new - mu_base <= -delta --- but
#' using per-seed paired differences.
#'
#' D14 chose the UNPAIRED test for the main grid, for a good reason: there, seed
#' i means something different to each engine, and the pairing that does exist in
#' Phases 2-3 evaporates at Phase 4 when initialization moves into the engine. A
#' test that silently changes meaning mid-project is worse than a weaker one.
#'
#' `tlda-compare.R` is a different experiment. It runs both engines from ONE set
#' of prepared inputs --- same tokens, same initial assignment, same prior --- so
#' the pairing is structural rather than incidental, and it stays valid however
#' initialization is later reorganized. That makes this the sharper instrument
#' where it applies: removing the seed-to-seed variation common to both arms
#' shrinks the standard deviation roughly fourfold, and with it the sample size
#' needed to resolve the margin.
#'
#' Report it alongside `noninf_test()`, not instead of it.
noninf_test_paired <- function(x_new, x_base,
                               margin_frac = BENCH_MARGIN_FRAC,
                               alpha = BENCH_ALPHA) {
  ok <- is.finite(x_new) & is.finite(x_base)
  x_new <- x_new[ok]; x_base <- x_base[ok]
  stopifnot(length(x_new) >= 2)

  base_mean <- mean(x_base)
  margin <- margin_frac * abs(base_mean)
  d <- x_new - x_base

  tt <- stats::t.test(d, alternative = "greater", mu = -margin)
  n <- length(d)
  df <- n - 1
  se <- stats::sd(d) / sqrt(n)

  tibble::tibble(
    n_pairs   = n,
    diff      = mean(d),
    sd_paired = stats::sd(d),
    margin    = margin,
    p_noninf  = tt$p.value,
    pass      = tt$p.value < alpha,
    mdd       = (stats::qt(1 - alpha, df) + stats::qt(0.8, df)) * se,
    p_any     = stats::t.test(d)$p.value
  )
}


#' Smallest true difference the gate could detect at this sample size
#'
#' Reported alongside the baseline so 6.1's power question --- "whether 20 seeds
#' resolves a 5% margin is empirical and Phase 1 will reveal it" --- has a
#' number attached. If `mdd` exceeds `margin`, the gate cannot resolve the
#' margin at this n, and 6.1 says to add seeds rather than widen the margin.
#'
#' Two-sample, one-sided, assuming the new arm's spread resembles the baseline's.
min_detectable_diff <- function(x_base, n_new, alpha = BENCH_ALPHA, power = 0.8) {
  n_base <- length(x_base)
  df <- n_base + n_new - 2
  se <- stats::sd(x_base) * sqrt(1 / n_base + 1 / n_new)
  (stats::qt(1 - alpha, df) + stats::qt(power, df)) * se
}


# ---- small utilities --------------------------------------------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

#' The commit the working tree is on, for recording in results
git_head <- function(repo = ".") {
  out <- try(
    system2("git", c("-C", repo, "rev-parse", "--short", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    silent = TRUE
  )
  if (inherits(out, "try-error") || length(out) == 0) NA_character_ else out[1]
}
