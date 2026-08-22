#!/usr/bin/env Rscript
# tlda-compare.R --- does the warpLDA engine match collapsed Gibbs under the
# tLDA MATRIX prior, which is the feature the package exists for?
#
#   Rscript warp-planning/benchmarks/tlda-compare.R [--workers=20] [--seeds=20]
#                                                   [--probe] [--collect-only]
#
# The main grid (run-benchmark.R) only ever exercises a scalar eta. This script
# covers the transfer case, and it can only be written NOW: it needs both
# engines, and Phase 6 deletes fit_lda_c().
#
# METHOD. Per cell and seed: split the corpus in half by document, fit a base
# model on the first half, build eta^(t) from it exactly as refit.tidylda.R:225
# does, then run BOTH engines on ONE set of prepared inputs. Sharing the
# initialization means the sampler is the only difference between the arms --
# same tokens, same starting assignment, same prior.

suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
  library(Matrix)
})

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "bench-lib.R"))
repo <- normalizePath(file.path(here, "..", ".."))

args <- commandArgs(trailingOnly = TRUE)
opt <- function(n, d) {
  h <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (length(h) == 0) d else sub(paste0("^--", n, "="), "", h[1])
}
n_workers    <- as.integer(opt("workers", 20))
n_seeds      <- as.integer(opt("seeds", 20))
# Per-K seed counts, calibrated from the 20-seed run: coherence at K=10 had
# mdd 1.8x (medium) and 3.2x (small) its margin, so 20 seeds cannot resolve it.
TLDA_SEEDS   <- list("10" = 1:200, "50" = 1:20)
probe        <- "--probe" %in% args
collect_only <- "--collect-only" %in% args

out_dir <- file.path(here, "results", "tlda")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

corpora <- list(
  small  = readRDS(file.path(here, "data", "small.rds")),
  medium = readRDS(file.path(here, "data", "nih-1000.rds"))
)

# This script reaches into internals (format_alpha, initialize_topic_counts,
# new_tidylda, fit_lda_c), so unlike run-benchmark.R it needs export_all.
suppressMessages(pkgload::load_all(repo, quiet = TRUE, export_all = TRUE))


#' Build the transfer scenario for one cell and seed
#'
#' Returns everything both engines need, already aligned, plus the held-out
#' half they are fitted on. `prior_weight = 1` reproduces refit()'s default
#' scaling: eta^(t) = w_star * beta^(t-1), with w_star the total posterior mass
#' behind each topic at t-1.
prepare_transfer <- function(dtm, k, seed) {
  set.seed(seed)
  half <- sample(nrow(dtm), floor(nrow(dtm) / 2))
  d1 <- dtm[half, , drop = FALSE]
  d2 <- dtm[-half, , drop = FALSE]

  # Keep both halves on the full vocabulary so no realignment is needed and the
  # comparison is purely about the sampler.
  base <- tidylda(
    data = d1, k = k, iterations = 200, burnin = 50,
    alpha = 0.1, eta = 0.05, calc_likelihood = FALSE, verbose = FALSE
  )

  # refit.tidylda.R:225 -- the tLDA recursion, eta^(t) = a * w_star * beta^(t-1)
  w_star <- rowSums(base$counts$Cv) + ncol(dtm) * 0.05
  eta_t  <- 1 * w_star * base$beta
  dimnames(eta_t) <- dimnames(base$beta)

  alpha <- format_alpha(base$alpha, k = k)
  eta   <- list(eta = eta_t, eta_class = "matrix")

  theta_initial <- predict(base, d2, method = "dot")

  set.seed(seed)
  counts <- initialize_topic_counts(
    dtm = d2, k = k, alpha = alpha$alpha, eta = eta$eta,
    beta_initial = base$beta, theta_initial = theta_initial,
    freeze_topics = FALSE, threads = 1
  )

  list(dtm = d2, k = k, alpha = alpha, eta = eta, counts = counts)
}


#' Score a fitted engine output the way a user would see it
score <- function(lda, prep, burnin) {
  model <- new_tidylda(
    lda = lda, dtm = prep$dtm, burnin = burnin, is_prediction = FALSE,
    alpha = prep$alpha, eta = prep$eta, optimize_alpha = FALSE,
    calc_r2 = TRUE, calc_likelihood = FALSE, call = NULL, threads = 1
  )
  m <- extract_metrics(model)
  m
}


run_tlda_row <- function(row, dtm) {
  path <- file.path(out_dir, paste0(row$id, ".rds"))
  res <- try({
    prep <- prepare_transfer(dtm, row$k, row$seed)
    cm <- list(Docs = prep$counts$Docs, Zd_in = prep$counts$Zd,
               Cd_in = prep$counts$Cd, Cv_in = prep$counts$Cv,
               Ck_in = prep$counts$Ck, alpha_in = prep$alpha$alpha,
               eta_in = prep$eta$eta)

    hg <- bench_hp("gibbs"); hw <- bench_hp("warp")

    set.seed(row$seed)
    tg <- system.time(g <- do.call(fit_lda_c, c(cm, list(
      iterations = hg$iterations, burnin = hg$burnin, optimize_alpha = FALSE,
      calc_likelihood = FALSE, Beta_in = prep$counts$Cv, freeze_topics = FALSE,
      threads = 1, verbose = FALSE))))[["elapsed"]]

    set.seed(row$seed)
    tw <- system.time(w <- do.call(fit_lda_warp, c(cm, list(
      iterations = hw$iterations, burnin = hw$burnin, calc_likelihood = FALSE,
      Beta_in = prep$counts$Cv, freeze_topics = FALSE,
      mh_steps = 1L, verbose = FALSE))))[["elapsed"]]

    bind_rows(
      bind_cols(tibble::tibble(engine = "gibbs"), score(g, prep, hg$burnin),
                tibble::tibble(elapsed_sec = tg)),
      bind_cols(tibble::tibble(engine = "warp"),  score(w, prep, hw$burnin),
                tibble::tibble(elapsed_sec = tw))
    )
  }, silent = TRUE)

  if (inherits(res, "try-error")) {
    res <- tibble::tibble(
      engine = c("gibbs", "warp"), r2 = NA_real_, coherence_mean = NA_real_,
      coherence_n_na = NA_integer_,
      status = paste0("error: ", conditionMessage(attr(res, "condition"))),
      elapsed_sec = NA_real_)
  }

  saveRDS(bind_cols(tibble::tibble(id = row$id, corpus = row$corpus,
                                   k = row$k, seed = row$seed), res), path)
  invisible(row$id)
}


grid <- if (probe) bench_grid(seeds = 1L) else
  dplyr::bind_rows(lapply(BENCH_KS, function(k)
    bench_grid(ks = k, seeds = TLDA_SEEDS[[as.character(k)]])))
done <- file.exists(file.path(out_dir, paste0(grid$id, ".rds")))
todo <- grid[!done, ]

cat(sprintf("tLDA grid=%d done=%d todo=%d workers=%d (each row runs BOTH engines)\n",
            nrow(grid), sum(done), nrow(todo), n_workers))

if (!collect_only && nrow(todo) > 0) {
  started <- Sys.time()
  invisible(parallel::mclapply(seq_len(nrow(todo)), function(i) {
    row <- todo[i, ]
    run_tlda_row(row, corpora[[row$corpus]])
  }, mc.cores = n_workers, mc.preschedule = FALSE))
  cat(sprintf("Finished in %.1f min\n",
              as.numeric(difftime(Sys.time(), started, units = "mins"))))
}


# ---- report ------------------------------------------------------------------

files <- sort(list.files(out_dir, pattern = "\\.rds$", full.names = TRUE))
if (length(files) == 0) stop("no results in ", out_dir)
res <- bind_rows(lapply(files, readRDS))

bad <- res |> filter(.data$status != "ok")
if (nrow(bad) > 0) { cat("\nnon-ok rows:\n"); print(as.data.frame(bad[, c("id", "engine", "status")])) }

cat("\nPer-cell means:\n\n")
print(as.data.frame(res |> group_by(.data$corpus, .data$k, .data$engine) |>
  summarise(n = n(), r2 = mean(.data$r2, na.rm = TRUE),
            coh = mean(.data$coherence_mean, na.rm = TRUE),
            sec = median(.data$elapsed_sec, na.rm = TRUE), .groups = "drop")))

if (!probe) {
  cat("\nGate (warp vs gibbs, matrix eta):\n\n")
  base <- res |> filter(.data$engine == "gibbs")
  new  <- res |> filter(.data$engine == "warp")
  cmp <- compare_runs(base, new)
  passed <- report(cmp)

  # Both engines here run from ONE shared initialization, so the pairing is
  # structural and the paired test is the sharper instrument -- see
  # noninf_test_paired() for why this is legitimate here and not in the main grid.
  cat("\nPaired gate (shared initialization, the sharper test):\n\n")
  cells <- base |> dplyr::distinct(.data$corpus, .data$k) |>
    dplyr::arrange(.data$corpus, .data$k)
  prs <- list()
  for (i in seq_len(nrow(cells))) {
    cl <- cells[i, ]
    for (nm in names(METRICS)) {
      col <- METRICS[[nm]]
      sel <- function(df) df[[col]][df$corpus == cl$corpus & df$k == cl$k][
        order(df$seed[df$corpus == cl$corpus & df$k == cl$k])]
      prs[[length(prs) + 1]] <- dplyr::bind_cols(
        tibble::tibble(cell = paste0(cl$corpus, "/k", cl$k), metric = nm),
        noninf_test_paired(sel(new), sel(base)))
    }
  }
  prs <- dplyr::bind_rows(prs)
  print(as.data.frame(prs |> dplyr::transmute(
    cell = .data$cell, metric = .data$metric, n = .data$n_pairs,
    diff = round(.data$diff, 5), margin = round(.data$margin, 5),
    mdd = round(.data$mdd, 5), p_noninf = signif(.data$p_noninf, 3),
    p_any = signif(.data$p_any, 3),
    verdict = ifelse(.data$pass, "PASS", "FAIL"))), row.names = FALSE)
  passed <- passed && all(prs$pass)
  cat("\nOverall:", if (passed) "PASS" else "FAIL", "\n")

  saveRDS(list(metrics = res, gate = cmp, hyperparams = list(
    gibbs = bench_hp("gibbs"), warp = bench_hp("warp")),
    paired = prs, git_head = git_head(repo), run_at = Sys.time()),
    file.path(here, "tlda-compare.rds"))
  cat("Assembled ->", file.path(here, "tlda-compare.rds"), "\n")
}
