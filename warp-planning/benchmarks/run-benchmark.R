#!/usr/bin/env Rscript
# run-benchmark.R --- run the benchmark grid for one sampler, persist the results
#
#   Rscript warp-planning/benchmarks/run-benchmark.R [options]
#
#     --engine=gibbs     which sampler; names the results subdirectory
#     --workers=20       parallel fits (never threads within a fit; see bench-lib)
#     --seeds=20         seeds per cell
#     --probe            seed 1 only (4 fits), report timings, do not assemble
#     --collect-only     skip fitting, just assemble what is already on disk
#     --out=FILE         assembled output; defaults per engine
#
# Resumability: each fit writes its own file the moment it completes, and fits
# whose file already exists are skipped. A crash at fit 70 of 80 costs one fit.

suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
})

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "bench-lib.R"))

repo <- normalizePath(file.path(here, "..", ".."))

# ---- arguments ---------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0) default else sub(paste0("^--", name, "="), "", hit[1])
}

engine       <- opt("engine", "gibbs")
n_workers    <- as.integer(opt("workers", 20))
n_seeds      <- opt("seeds", NA_character_)   # NA = per-K calibration
probe        <- "--probe" %in% args
collect_only <- "--collect-only" %in% args

# 5abaa96 is the last commit that touched the sampler; `warp` HEAD has moved on
# with docs-only commits. The name records the sampler state, and the actual HEAD
# is stored inside the file so there is no later ambiguity.
default_out <- c(gibbs = "baseline-5abaa96.rds")
out_file <- opt("out", if (engine %in% names(default_out)) {
  default_out[[engine]]
} else {
  paste0("run-", engine, ".rds")
})

out_dir <- file.path(here, "results", engine)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---- corpora and package -----------------------------------------------------

corpora <- list(
  small  = readRDS(file.path(here, "data", "small.rds")),
  medium = readRDS(file.path(here, "data", "nih-1000.rds"))
)

# Load (and compile) in the parent, before forking, so every worker inherits a
# built package rather than racing to build it.
suppressMessages(pkgload::load_all(repo, quiet = TRUE, export_all = FALSE))


# ---- the work ----------------------------------------------------------------

grid <- bench_grid(
  seeds = if (probe) 1L else if (!is.na(n_seeds)) seq_len(as.integer(n_seeds)) else NULL
)

done <- file.exists(file.path(out_dir, paste0(grid$id, ".rds")))
todo <- grid[!done, ]

cat(sprintf("engine=%s  grid=%d  done=%d  todo=%d  workers=%d\n",
            engine, nrow(grid), sum(done), nrow(todo), n_workers))

if (!collect_only && nrow(todo) > 0) {
  started <- Sys.time()

  # mc.preschedule = FALSE gives dynamic dispatch. Fit cost is wildly uneven
  # across the grid --- medium/K=50 is orders of magnitude more work than
  # small/K=10 --- so a static split would leave most cores idle at the end.
  invisible(parallel::mclapply(
    seq_len(nrow(todo)),
    function(i) {
      row <- todo[i, ]
      run_grid_row(row, dtm = corpora[[row$corpus]], out_dir = out_dir)
    },
    mc.cores = n_workers,
    mc.preschedule = FALSE
  ))

  cat(sprintf("Fitting finished in %.1f min\n",
              as.numeric(difftime(Sys.time(), started, units = "mins"))))
}


# ---- report / assemble -------------------------------------------------------

res <- collect_results(out_dir)

failed <- res$metrics |> dplyr::filter(.data$status != "ok")
if (nrow(failed) > 0) {
  cat("\nFits with a non-ok status:\n")
  print(as.data.frame(failed[, c("id", "status")]))
}

cat("\nPer-cell summary:\n\n")
summary_tbl <- res$metrics |>
  dplyr::group_by(.data$corpus, .data$k) |>
  # Summary names must not collide with the columns they read. Inside
  # summarise() the data mask updates as each expression is evaluated, so a
  # column named `r2` computed from `.data$r2` makes every later `.data$r2`
  # resolve to the length-1 result --- and sd() of one value is a silent NA.
  dplyr::summarise(
    n          = dplyr::n(),
    r2_mean    = mean(.data$r2, na.rm = TRUE),
    r2_sd      = stats::sd(.data$r2, na.rm = TRUE),
    coh_mean   = mean(.data$coherence_mean, na.rm = TRUE),
    coh_sd     = stats::sd(.data$coherence_mean, na.rm = TRUE),
    sec_median = stats::median(.data$elapsed_sec, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(summary_tbl))

if (probe) {
  total_sec <- sum(summary_tbl$sec_median * lengths(BENCH_SEEDS[as.character(summary_tbl$k)]))
  cat(sprintf(
    "\nProbe only; nothing assembled.\nProjected full grid: %.1f core-hours, ~%.1f wall-clock min at %d workers.\n",
    total_sec / 3600, total_sec / 60 / n_workers, n_workers
  ))
} else {
  saveRDS(
    list(
      engine      = engine,
      metrics     = res$metrics,
      curves      = res$curves,
      summary     = summary_tbl,
      corpora     = readRDS(file.path(here, "data", "corpora-manifest.rds")),
      hyperparams = BENCH_HP,
      git_head    = git_head(repo),
      run_at      = Sys.time(),
      session     = utils::sessionInfo()
    ),
    file.path(here, out_file)
  )
  cat("\nAssembled ->", file.path(here, out_file), "\n")
}
