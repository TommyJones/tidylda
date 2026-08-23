context("warpLDA engine")

# Tests for fit_lda_warp(), the Phase 2 sampler. See
# warp-planning/warplda-roadmap.md section 6.

dtm <- nih_sample_dtm[1:50, ]

# Build the inputs fit_lda_warp() expects, the same way tidylda_bridge() does.
# Since Phase 4 (D16) the engine takes the DTM directly and does its own
# initialization, so there is no lexicon to marshal.
setup <- function(k = 5, seed = 1) {
  v <- ncol(dtm)
  alpha <- format_alpha(0.1, k)
  eta <- format_eta(0.05, k, v)
  set.seed(seed)
  counts <- initialize_topic_counts(
    dtm = dtm, k = k, alpha = alpha$alpha, eta = eta$eta, threads = 1
  )
  list(k = k, v = v, alpha = alpha, eta = eta, counts = counts, n = sum(dtm))
}

run_warp <- function(s, seed = 42, ...) {
  set.seed(seed)
  args <- list(
    dtm_in = dtm, Cd_start = s$counts$Cd_start, alpha_in = s$alpha$alpha,
    eta_in = s$eta$eta, iterations = 30, burnin = 10,
    calc_likelihood = TRUE, Beta_in = s$counts$beta_initial, verbose = FALSE
  )
  do.call(fit_lda_warp, utils::modifyList(args, list(...)))
}

test_that("output contract matches fit_lda_c", {
  s <- setup()
  m <- run_warp(s)

  expect_setequal(
    names(m),
    c("Cd", "Cv", "Ck", "Cd_mean", "Cv_mean", "Cd_sum", "Cv_sum",
      "log_likelihood", "alpha", "eta")
  )

  # Cd is documents x topics; Cv is topics x words. new_tidylda() and
  # posterior.tidylda() both depend on these orientations (roadmap section 5).
  expect_equal(dim(m$Cd), c(nrow(dtm), s$k))
  expect_equal(dim(m$Cv), c(s$k, s$v))
  expect_equal(dim(m$Cd_mean), c(nrow(dtm), s$k))
  expect_equal(dim(m$Cv_mean), c(s$k, s$v))
  expect_equal(nrow(m$log_likelihood), 3)
})


test_that("token counts are conserved and the marginals agree", {
  s <- setup()
  m <- run_warp(s)

  # Every token is assigned to exactly one topic, so all three views of the
  # count structure must total N.
  expect_equal(sum(m$Cd), s$n)
  expect_equal(sum(m$Cv), s$n)
  expect_equal(sum(m$Ck), s$n)

  # Ck is the topic marginal of both matrices.
  expect_equal(as.numeric(m$Ck), unname(colSums(m$Cd)))
  expect_equal(as.numeric(m$Ck), unname(rowSums(m$Cv)))

  # Document lengths are fixed, so the posterior mean over topics recovers them.
  expect_equal(unname(rowSums(m$Cd_mean)), unname(Matrix::rowSums(dtm)))
  expect_equal(sum(m$Cv_mean), s$n)
})


test_that("Cd is current on return even without the likelihood", {
  # C_word is current after the word pass, but C_doc is not: the word pass
  # reassigns tokens after the doc pass last rebuilt it. The likelihood block
  # happens to refresh C_doc, so a stale Cd would hide whenever the likelihood
  # is on. new_tidylda() reads Cd directly when burnin == -1.
  s <- setup()
  set.seed(42)
  m <- run_warp(s, burnin = -1, calc_likelihood = FALSE)
  expect_equal(sum(m$Cd), s$n)
  expect_equal(as.numeric(m$Ck), unname(colSums(m$Cd)))
})


test_that("set.seed makes fits reproducible", {
  s <- setup()
  a <- run_warp(s, seed = 7)
  b <- run_warp(s, seed = 7)
  d <- run_warp(s, seed = 8)

  expect_identical(a$Cd, b$Cd)
  expect_identical(a$Cv, b$Cv)
  expect_false(identical(a$Cd, d$Cd))
})


test_that("mh_steps runs at values above the default", {
  s <- setup()
  for (steps in c(1, 2, 4)) {
    m <- run_warp(s, mh_steps = steps)
    expect_equal(sum(m$Cd), s$n)
    expect_equal(as.numeric(m$Ck), unname(colSums(m$Cd)))
  }
})


test_that("likelihood_every controls evaluation frequency, not the chain", {
  s <- setup()

  every_1  <- run_warp(s, likelihood_every = 1)
  every_10 <- run_warp(s, likelihood_every = 10)

  expect_equal(ncol(every_1$log_likelihood), 30)
  expect_equal(as.numeric(every_1$log_likelihood[1, ]), 0:29)

  # 0, 10, 20, and the always-recorded final iteration.
  expect_equal(as.numeric(every_10$log_likelihood[1, ]), c(0, 10, 20, 29))

  # The interval is a diagnostic only: the chain is untouched, so the same seed
  # must give the same posterior regardless of how often it was measured.
  expect_identical(every_1$Cd, every_10$Cd)
  expect_identical(every_1$Cv, every_10$Cv)
})


test_that("burnin averaging is off when burnin is -1", {
  s <- setup()
  set.seed(42)
  m <- run_warp(s, iterations = 10, burnin = -1, calc_likelihood = FALSE)
  expect_equal(length(m$Cd_mean), 0)
  expect_equal(length(m$Cv_mean), 0)
})


test_that("invalid arguments are rejected", {
  s <- setup()
  expect_error(run_warp(s, mh_steps = 0), "mh_steps")
  expect_error(
    run_warp(s, iterations = 10, burnin = 10, calc_likelihood = FALSE),
    "burnin"
  )
})


test_that("vector and matrix eta both fit", {
  vec <- tidylda(data = dtm, k = 4, iterations = 20, burnin = 10,
                 eta = rep(0.05, ncol(dtm)), calc_likelihood = FALSE,
                 verbose = FALSE)
  expect_length(vec$eta, ncol(dtm))

  mat <- tidylda(data = dtm, k = 4, iterations = 20, burnin = 10,
                 eta = matrix(0.05, nrow = 4, ncol = ncol(dtm)),
                 calc_likelihood = FALSE, verbose = FALSE)
  expect_true(inherits(mat$eta, "matrix"))
  expect_equal(dim(mat$eta), c(4L, ncol(dtm)))
})


test_that("a constant matrix eta matches the scalar it is made of", {
  # The sharpest available check that generalizing to a matrix prior did not
  # move the target: a K x V matrix filled with the scalar describes exactly the
  # same model, so the two must agree beyond sampling noise.
  skip_on_cran()

  r2 <- function(eta, seed) {
    set.seed(seed)
    tidylda(data = dtm, k = 5, iterations = 300, burnin = 75, alpha = 0.1,
            eta = eta, calc_likelihood = FALSE, calc_r2 = TRUE,
            verbose = FALSE)$r2
  }
  seeds <- 1:6
  scal <- vapply(seeds, function(s) r2(0.05, s), numeric(1))
  matx <- vapply(seeds, function(s) r2(matrix(0.05, 5, ncol(dtm)), s), numeric(1))

  expect_lt(abs(mean(matx) - mean(scal)) / mean(scal), 0.05)
})


test_that("frozen topics leave the word side untouched", {
  s <- setup(k = 5)
  # beta_initial is already a K x V matrix of P(token|topic); normalize so the
  # frozen proposal is an exact distribution.
  beta <- s$counts$beta_initial
  beta <- beta / rowSums(beta)

  set.seed(42)
  m <- run_warp(s, calc_likelihood = FALSE, Beta_in = beta, freeze_topics = TRUE)

  # C^v and C_k are not part of the frozen target, so they are not built.
  expect_equal(length(m$Cv), 0)
  expect_equal(length(m$Cv_mean), 0)

  # The document side still has to be a valid assignment of every token.
  expect_equal(sum(m$Cd), s$n)
  expect_equal(unname(rowSums(m$Cd_mean)), unname(Matrix::rowSums(dtm)))
})


test_that("predict() runs on the warp engine and returns a distribution", {
  set.seed(1)
  m <- tidylda(data = dtm, k = 5, iterations = 40, burnin = 10,
               calc_likelihood = FALSE, verbose = FALSE)

  p <- predict(m, nih_sample_dtm[51:70, ], method = "gibbs",
               iterations = 40, burnin = 10, verbose = FALSE)

  expect_equal(nrow(p), 20)
  expect_equal(ncol(p), 5)
  expect_equal(unname(rowSums(p)), rep(1, 20), tolerance = 1e-8)
})


test_that("theta reflects an asymmetric alpha", {
  # new_tidylda() used to compute theta as t(t(Cd + alpha)), a double transpose
  # that added a K-vector to a D x K matrix and so recycled alpha diagonally
  # across topics rather than adding alpha[k] to topic k. Invisible for a
  # symmetric alpha; wrong for every other one.
  set.seed(1)
  m <- tidylda(data = dtm, k = 4, iterations = 30, burnin = 10,
               alpha = c(50, 0.01, 0.01, 0.01), eta = 0.05,
               calc_likelihood = FALSE, verbose = FALSE)

  # Topic 1 carries 5000x the prior mass of the others, so it must be the
  # largest component of theta in every document.
  expect_true(all(apply(m$theta, 1, which.max) == 1))
  expect_gt(mean(m$theta[, 1]), 3 * max(colMeans(m$theta[, -1, drop = FALSE])))
})


test_that("the sampler targets the LDA posterior, not merely a stationary one", {
  # A wrong eta_bar or a mis-derived acceptance ratio still yields a valid MCMC
  # chain -- it just converges somewhere else. Comparing against the collapsed
  # Gibbs sampler on the same data, from the same initialization, run long
  # enough that both have settled, is what distinguishes those two cases.
  skip_on_cran()

  s <- setup(k = 5, seed = 3)

  # Gibbs still needs the lexicon; the warp engine builds its own.
  lex <- create_lexicon(Cd_in = s$counts$Cd_start,
                        Beta_in = s$counts$beta_initial,
                        dtm_in = dtm, alpha = s$alpha$alpha,
                        freeze_topics = FALSE)

  set.seed(11)
  g <- fit_lda_c(
    Docs = lex$Docs, Zd_in = lex$Zd, Cd_in = lex$Cd, Cv_in = lex$Cv,
    Ck_in = lex$Ck, alpha_in = s$alpha$alpha,
    eta_in = s$eta$eta, iterations = 400, burnin = 100, optimize_alpha = FALSE,
    calc_likelihood = TRUE, Beta_in = s$counts$beta_initial,
    freeze_topics = FALSE, threads = 1, verbose = FALSE
  )

  set.seed(11)
  w <- run_warp(s, seed = 11, iterations = 2000, burnin = 500,
                calc_likelihood = TRUE, likelihood_every = 1)

  g_ll <- g$log_likelihood[2, ncol(g$log_likelihood)]
  w_ll <- w$log_likelihood[2, ncol(w$log_likelihood)]

  # Both should land on the same plateau. The tolerance is loose because these
  # are two independent chains, but it is far tighter than the gap a wrong
  # target would produce -- porting the reference's beta_bar = K * eta instead
  # of V * eta moves this by orders of magnitude more than 2%.
  expect_lt(abs(w_ll - g_ll) / abs(g_ll), 0.02)
})


test_that("initialization is informed by the priors, not uniform", {
  # Phase 4 verified the fused initialization against create_lexicon() by exact
  # equality. Phase 4.5 replaced lsamp_one()'s per-token O(K log K) sort with a
  # constant-work draw, so that comparison no longer holds and the test was
  # retired deliberately rather than left to rot.
  #
  # What replaces it tests the property D8 actually cares about, and which no
  # future change to the sampling algorithm should break: each token's starting
  # topic is drawn from P(z) proportional to beta[k, v] * (Cd_start[d, k] +
  # alpha[k]), so the initial assignment tracks the priors it was given. A
  # uniform-random start -- warpLDA's own, which tidylda discards -- would not.
  d <- nih_sample_dtm[1:20, ]
  k <- 4
  v <- ncol(d)
  alpha <- rep(0.1, k)

  # Flat beta, so the document prior alone decides. Give document i almost all
  # of its mass on topic (i mod k).
  beta_flat <- matrix(1 / v, nrow = k, ncol = v)
  Cd_start <- matrix(1e-6, nrow = nrow(d), ncol = k)
  favoured <- (seq_len(nrow(d)) - 1) %% k + 1
  for (i in seq_len(nrow(d))) Cd_start[i, favoured[i]] <- 1000

  set.seed(4)
  m <- fit_lda_warp(dtm_in = d, Cd_start = Cd_start, alpha_in = alpha,
                    eta_in = matrix(0.05, nrow = k, ncol = v),
                    iterations = 0, burnin = -1, calc_likelihood = FALSE,
                    Beta_in = beta_flat, freeze_topics = FALSE, verbose = FALSE)

  # Every token still accounted for.
  expect_equal(sum(m$Cd), sum(d))
  expect_equal(sum(m$Cv), sum(d))
  expect_equal(as.numeric(m$Ck), unname(colSums(m$Cd)))

  # Each document should have landed overwhelmingly on the topic it was steered
  # toward. The prior ratio is 1000 : 1e-6, so anything less than near-total
  # concentration means the initialization is ignoring it.
  got <- apply(m$Cd, 1, which.max)
  expect_equal(got, favoured)
  expect_gt(mean(m$Cd[cbind(seq_len(nrow(d)), favoured)] / rowSums(m$Cd)), 0.99)
})


test_that("initialization is reproducible under set.seed", {
  d <- nih_sample_dtm[1:20, ]
  k <- 4
  al <- format_alpha(0.1, k)
  et <- format_eta(0.05, k, ncol(d))

  init <- function(seed) {
    set.seed(seed)
    p <- initialize_topic_counts(dtm = d, k = k, alpha = al$alpha,
                                 eta = et$eta, threads = 1)
    fit_lda_warp(dtm_in = d, Cd_start = p$Cd_start, alpha_in = al$alpha,
                 eta_in = et$eta, iterations = 0, burnin = -1,
                 calc_likelihood = FALSE, Beta_in = p$beta_initial,
                 freeze_topics = FALSE, verbose = FALSE)$Cd
  }
  expect_identical(init(21), init(21))
  expect_false(identical(init(21), init(22)))
})


test_that("results are identical at any thread count (D12)", {
  # This is Phase 5's exit criterion and the reason C_k is snapshot within a
  # pass rather than updated atomically. Atomics would remove the race but leave
  # what each work item READS dependent on interleaving, so the answer would
  # drift with thread count. Here every work item sees the same C_k however
  # chunks land on threads, and the deltas merge by integer addition, which is
  # associative and exact.
  skip_on_cran()

  d <- nih_sample_dtm[1:60, ]
  k <- 6
  al <- format_alpha(0.1, k)
  et <- format_eta(0.05, k, ncol(d))
  set.seed(3)
  p <- initialize_topic_counts(dtm = d, k = k, alpha = al$alpha, eta = et$eta,
                               threads = 1)

  run <- function(th) {
    set.seed(77)
    fit_lda_warp(dtm_in = d, Cd_start = p$Cd_start, alpha_in = al$alpha,
                 eta_in = et$eta, iterations = 40, burnin = 15,
                 calc_likelihood = TRUE, Beta_in = p$beta_initial,
                 freeze_topics = FALSE, likelihood_every = 5, mh_steps = 1L,
                 threads = as.integer(th), verbose = FALSE)
  }

  ref <- run(1)
  for (th in c(2, 4, 8)) {
    m <- run(th)
    expect_identical(ref$Cd, m$Cd, info = paste("threads =", th))
    expect_identical(ref$Cv, m$Cv, info = paste("threads =", th))
    expect_identical(ref$Ck, m$Ck, info = paste("threads =", th))
    expect_identical(ref$Cd_mean, m$Cd_mean, info = paste("threads =", th))
    expect_identical(ref$log_likelihood, m$log_likelihood,
                     info = paste("threads =", th))
  }

  # And the invariants must still hold when threaded.
  m <- run(8)
  expect_equal(sum(m$Cd), sum(d))
  expect_equal(sum(m$Cv), sum(d))
  expect_equal(as.numeric(m$Ck), unname(colSums(m$Cd)))
})


test_that("tidylda() accepts threads without the stale batching warning", {
  # The old warning -- fewer than 100 documents per thread gives "a poor fit" --
  # described the abandoned batched Gibbs, where the partition changed the model.
  # Under D12 it cannot, so the warning would now be false.
  d <- nih_sample_dtm[1:40, ]
  expect_no_warning(
    m <- tidylda(d, k = 4, iterations = 20, burnin = 5, calc_likelihood = FALSE,
                 threads = 4, verbose = FALSE)
  )
  expect_s3_class(m, "tidylda")
})


test_that("initialization is identical at any thread count", {
  # Phase 5.5 parallelized initialization by swapping R's RNG for D12's
  # per-work-item generator. R's RNG is main-thread-only (roadmap section 5), so
  # while initialization drew from it the O(N*K) setup could not be threaded --
  # and it had become the Amdahl bottleneck, exceeding parallel sampling by 6x
  # at K=200.
  skip_on_cran()

  d <- nih_sample_dtm[1:50, ]
  k <- 8
  al <- format_alpha(0.1, k)
  et <- format_eta(0.05, k, ncol(d))
  set.seed(5)
  p <- initialize_topic_counts(dtm = d, k = k, alpha = al$alpha, eta = et$eta,
                               threads = 1)

  init <- function(th) {
    set.seed(123)
    fit_lda_warp(dtm_in = d, Cd_start = p$Cd_start, alpha_in = al$alpha,
                 eta_in = et$eta, iterations = 0, burnin = -1,
                 calc_likelihood = FALSE, Beta_in = p$beta_initial,
                 freeze_topics = FALSE, threads = as.integer(th), verbose = FALSE)
  }

  ref <- init(1)
  for (th in c(2, 4, 8)) {
    m <- init(th)
    expect_identical(ref$Cd, m$Cd, info = paste("threads =", th))
    expect_identical(ref$Cv, m$Cv, info = paste("threads =", th))
    expect_identical(ref$Ck, m$Ck, info = paste("threads =", th))
  }
  expect_equal(sum(ref$Cd), sum(d))
  expect_equal(sum(ref$Cv), sum(d))
})


test_that("initialization draws from the distribution it is supposed to", {
  # Checked against ground truth rather than against the previous implementation.
  # The target is known in closed form --
  #   P(z = k | d, v) proportional to beta[k, v] * (Cd_start[d, k] + alpha[k])
  # -- so the expected document-topic counts can be computed directly and
  # compared with what the sampler produces. That would catch the old and new
  # code being wrong in the same way, which diffing them against each other
  # could not.
  skip_on_cran()

  d <- nih_sample_dtm[1:30, ]
  k <- 5
  al <- format_alpha(0.1, k)
  et <- format_eta(0.05, k, ncol(d))
  set.seed(5)
  p <- initialize_topic_counts(dtm = d, k = k, alpha = al$alpha, eta = et$eta,
                               threads = 1)

  # Cd_start reaches C++ as an IntegerMatrix, so it is truncated there.
  Cd0 <- floor(p$Cd_start)
  dm <- as.matrix(d)

  expected <- matrix(0, nrow(d), k)
  for (i in seq_len(nrow(d))) {
    nz <- which(dm[i, ] > 0)
    w <- outer(rep(1, length(nz)), Cd0[i, ] + al$alpha) *
      t(p$beta_initial[, nz, drop = FALSE])
    expected[i, ] <- colSums(dm[i, nz] * (w / rowSums(w)))
  }

  reps <- 200
  realized <- matrix(0, nrow(d), k)
  for (r in seq_len(reps)) {
    set.seed(2000 + r)
    realized <- realized + fit_lda_warp(
      dtm_in = d, Cd_start = p$Cd_start, alpha_in = al$alpha, eta_in = et$eta,
      iterations = 0, burnin = -1, calc_likelihood = FALSE,
      Beta_in = p$beta_initial, freeze_topics = FALSE, threads = 1L,
      verbose = FALSE)$Cd
  }
  realized <- realized / reps

  # Totals are exact -- every token is assigned exactly once.
  expect_equal(sum(realized), sum(d))
  expect_equal(sum(expected), sum(d), tolerance = 1e-8)

  # A Poisson standard error overstates the true multinomial variance, so these
  # z-scores are conservative; anything beyond 5 would be a real discrepancy.
  z <- (realized - expected) / sqrt(pmax(expected, 1e-9) / reps)
  expect_lt(max(abs(z)), 5)
  expect_gt(cor(as.vector(expected), as.vector(realized)), 0.999)
})
