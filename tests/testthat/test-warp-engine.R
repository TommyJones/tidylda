context("warpLDA engine")

# Tests for fit_lda_warp(), the Phase 2 sampler. See
# warp-planning/warplda-roadmap.md section 6.

dtm <- nih_sample_dtm[1:50, ]

# Build the inputs fit_lda_warp() expects, the same way tidylda_bridge() does.
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
    Docs = s$counts$Docs, Zd_in = s$counts$Zd, Cd_in = s$counts$Cd,
    Cv_in = s$counts$Cv, Ck_in = s$counts$Ck, alpha_in = s$alpha$alpha,
    eta_in = s$eta$eta, iterations = 30, burnin = 10,
    calc_likelihood = TRUE, Beta_in = s$counts$Cv, verbose = FALSE
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
  m <- fit_lda_warp(
    Docs = s$counts$Docs, Zd_in = s$counts$Zd, Cd_in = s$counts$Cd,
    Cv_in = s$counts$Cv, Ck_in = s$counts$Ck, alpha_in = s$alpha$alpha,
    eta_in = s$eta$eta, iterations = 30, burnin = -1,
    calc_likelihood = FALSE, Beta_in = s$counts$Cv, verbose = FALSE
  )
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
  m <- fit_lda_warp(
    Docs = s$counts$Docs, Zd_in = s$counts$Zd, Cd_in = s$counts$Cd,
    Cv_in = s$counts$Cv, Ck_in = s$counts$Ck, alpha_in = s$alpha$alpha,
    eta_in = s$eta$eta, iterations = 10, burnin = -1,
    calc_likelihood = FALSE, Beta_in = s$counts$Cv, verbose = FALSE
  )
  expect_equal(length(m$Cd_mean), 0)
  expect_equal(length(m$Cv_mean), 0)
})


test_that("invalid arguments are rejected", {
  s <- setup()
  expect_error(run_warp(s, mh_steps = 0), "mh_steps")
  expect_error(
    fit_lda_warp(
      Docs = s$counts$Docs, Zd_in = s$counts$Zd, Cd_in = s$counts$Cd,
      Cv_in = s$counts$Cv, Ck_in = s$counts$Ck, alpha_in = s$alpha$alpha,
      eta_in = s$eta$eta, iterations = 10, burnin = 10,
      calc_likelihood = FALSE, Beta_in = s$counts$Cv, verbose = FALSE
    ),
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
  beta <- s$counts$Cv + 0.05
  beta <- beta / rowSums(beta)

  set.seed(42)
  m <- fit_lda_warp(
    Docs = s$counts$Docs, Zd_in = s$counts$Zd, Cd_in = s$counts$Cd,
    Cv_in = s$counts$Cv, Ck_in = s$counts$Ck, alpha_in = s$alpha$alpha,
    eta_in = s$eta$eta, iterations = 30, burnin = 10,
    calc_likelihood = FALSE, Beta_in = beta, freeze_topics = TRUE,
    verbose = FALSE
  )

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
  common <- list(
    Docs = s$counts$Docs, Zd_in = s$counts$Zd, Cd_in = s$counts$Cd,
    Cv_in = s$counts$Cv, Ck_in = s$counts$Ck, alpha_in = s$alpha$alpha
  )

  set.seed(11)
  g <- do.call(fit_lda_c, c(common, list(
    eta_in = s$eta$eta, iterations = 400, burnin = 100, optimize_alpha = FALSE,
    calc_likelihood = TRUE, Beta_in = s$counts$Cv, freeze_topics = FALSE,
    threads = 1, verbose = FALSE
  )))

  set.seed(11)
  w <- do.call(fit_lda_warp, c(common, list(
    eta_in = s$eta$eta, iterations = 2000, burnin = 500,
    calc_likelihood = TRUE, Beta_in = s$counts$Cv, likelihood_every = 1,
    verbose = FALSE
  )))

  g_ll <- g$log_likelihood[2, ncol(g$log_likelihood)]
  w_ll <- w$log_likelihood[2, ncol(w$log_likelihood)]

  # Both should land on the same plateau. The tolerance is loose because these
  # are two independent chains, but it is far tighter than the gap a wrong
  # target would produce -- porting the reference's beta_bar = K * eta instead
  # of V * eta moves this by orders of magnitude more than 2%.
  expect_lt(abs(w_ll - g_ll) / abs(g_ll), 0.02)
})
