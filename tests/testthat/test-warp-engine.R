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
    eta_in = as.matrix(s$eta$eta), iterations = 30, burnin = 10,
    calc_likelihood = TRUE, Beta_in = s$counts$beta_initial, verbose = FALSE
  )
  do.call(fit_lda_warp, utils::modifyList(args, list(...)))
}

test_that("the output contract is what new_tidylda expects", {
  s <- setup()
  m <- run_warp(s)

  expect_setequal(
    names(m),
    c("Cd", "Cv", "Ck", "Cd_mean", "Cv_mean",
      "log_likelihood", "alpha", "eta")
  )

  # These orientations are a contract: new_tidylda() and posterior.tidylda()
  # both index them directly (roadmap section 5). They were inherited from
  # fit_lda_c(), which Phase 6 deleted, so this test is now the only thing
  # holding them fixed.
  # Both are <major> by topics since D17: Cd documents-by-topics, Cv
  # words-by-topics. The engine no longer transposes Cv on output; new_tidylda()
  # does it once when forming the public, topics-by-words `beta`.
  #
  # Only the pair new_tidylda() will read is materialized. run_warp() uses
  # burnin = 10, so that is the averaged pair; the raw pair comes back 0 x 0.
  expect_equal(dim(m$Cd_mean), c(nrow(dtm), s$k))
  expect_equal(dim(m$Cv_mean), c(s$v, s$k))
  expect_equal(dim(m$Cd), c(0L, 0L))
  expect_equal(dim(m$Cv), c(0L, 0L))

  # And the other way round without burn-in.
  m_raw <- run_warp(s, burnin = -1)
  expect_equal(dim(m_raw$Cd), c(nrow(dtm), s$k))
  expect_equal(dim(m_raw$Cv), c(s$v, s$k))
  expect_equal(dim(m_raw$Cd_mean), c(0L, 0L))
  expect_equal(dim(m_raw$Cv_mean), c(0L, 0L))

  expect_equal(nrow(m$log_likelihood), 3)
})


test_that("token counts are conserved and the marginals agree", {
  s <- setup()
  m     <- run_warp(s)                 # burnin = 10: averaged pair
  m_raw <- run_warp(s, burnin = -1)    # raw pair

  # Every token is assigned to exactly one topic, so every view of the count
  # structure must total N, averaged or not.
  expect_equal(sum(m_raw$Cd), s$n)
  expect_equal(sum(m_raw$Cv), s$n)
  expect_equal(sum(m_raw$Ck), s$n)
  expect_equal(sum(m$Ck), s$n)

  # Ck is the topic marginal of both matrices.
  expect_equal(as.numeric(m_raw$Ck), unname(colSums(m_raw$Cd)))
  expect_equal(as.numeric(m_raw$Ck), unname(colSums(m_raw$Cv)))

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
  # burnin = -1 so the raw counts are the pair that gets materialized. Comparing
  # Cd under the default burnin would compare two 0 x 0 matrices, which are
  # identical no matter what the sampler did.
  a <- run_warp(s, seed = 7, burnin = -1)
  b <- run_warp(s, seed = 7, burnin = -1)
  d <- run_warp(s, seed = 8, burnin = -1)

  expect_identical(a$Cd, b$Cd)
  expect_identical(a$Cv, b$Cv)
  expect_false(identical(a$Cd, d$Cd))
  expect_false(identical(a$Cv, d$Cv))
})


test_that("mh_steps runs at values above the default", {
  s <- setup()
  for (steps in c(1, 2, 4)) {
    m <- run_warp(s, mh_steps = steps, burnin = -1)
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

  # The document side still has to be a valid assignment of every token. This
  # run uses the default burnin = 10, so the averaged pair is the one that
  # exists; Cd itself is 0 x 0.
  expect_equal(sum(m$Cd_mean), s$n)
  expect_equal(unname(rowSums(m$Cd_mean)), unname(Matrix::rowSums(dtm)))
})


test_that("predict() runs on the warp engine and returns a distribution", {
  set.seed(1)
  m <- tidylda(data = dtm, k = 5, iterations = 40, burnin = 10,
               calc_likelihood = FALSE, verbose = FALSE)

  p <- predict(m, nih_sample_dtm[51:70, ], method = "mh",
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


# A textbook collapsed Gibbs sampler, written independently of the engine: no
# shared headers, no shared RNG, no shared data structures. Deliberately slow
# and obvious. This replaced a comparison against the package's own fit_lda_c(),
# which Phase 6 deleted -- and it is the stronger check, because a C++ reference
# sharing headers with the engine could agree with it for the wrong reason.
reference_gibbs <- function(dtm, k, alpha, eta, iterations, seed) {
  set.seed(seed)
  D <- nrow(dtm); V <- ncol(dtm)
  docs <- lapply(seq_len(D), function(d) rep(seq_len(V), times = dtm[d, ]))
  z <- lapply(docs, function(w) sample.int(k, length(w), replace = TRUE))

  Cd <- matrix(0, D, k); Cv <- matrix(0, k, V); Ck <- numeric(k)
  for (d in seq_len(D)) for (i in seq_along(docs[[d]])) {
    kk <- z[[d]][i]; v <- docs[[d]][i]
    Cd[d, kk] <- Cd[d, kk] + 1; Cv[kk, v] <- Cv[kk, v] + 1; Ck[kk] <- Ck[kk] + 1
  }
  for (it in seq_len(iterations)) {
    for (d in seq_len(D)) for (i in seq_along(docs[[d]])) {
      v <- docs[[d]][i]; old <- z[[d]][i]
      Cd[d, old] <- Cd[d, old] - 1
      Cv[old, v] <- Cv[old, v] - 1
      Ck[old]    <- Ck[old] - 1
      # The collapsed conditional. Note the normalizer is Ck + V*eta -- a sum
      # over the VOCABULARY, which is the term text2vec's reference gets wrong.
      p <- ((Cv[, v] + eta) / (Ck + V * eta)) * (Cd[d, ] + alpha)
      new <- sample.int(k, 1, prob = p)
      z[[d]][i] <- new
      Cd[d, new] <- Cd[d, new] + 1
      Cv[new, v] <- Cv[new, v] + 1
      Ck[new]    <- Ck[new] + 1
    }
  }
  list(Cd = Cd, Cv = Cv, Ck = Ck, docs = docs)
}

# Log probability of the corpus under the model. Invariant to topic relabelling,
# which beta and theta are not -- two samplers will land on different topic
# orderings, so they cannot be compared elementwise.
reference_loglik <- function(Cd, Cv, docs, alpha, eta) {
  theta <- Cd + rep(alpha, each = nrow(Cd)); theta <- theta / rowSums(theta)
  beta <- Cv + eta; beta <- beta / rowSums(beta)
  total <- 0
  for (d in seq_along(docs)) {
    for (v in docs[[d]]) total <- total + log(sum(theta[d, ] * beta[, v]))
  }
  total
}


test_that("the sampler targets the LDA posterior, not merely a stationary one", {
  # A wrong eta_bar or a mis-derived acceptance ratio still yields a valid MCMC
  # chain -- it just converges somewhere else. Comparing against an independent
  # collapsed Gibbs sampler, run long enough that both have settled, is what
  # distinguishes those two cases.
  skip_on_cran()

  d <- nih_sample_dtm[1:15, ]
  d <- d[, Matrix::colSums(d) > 0]
  d <- d[, order(Matrix::colSums(d), decreasing = TRUE)[1:40]]
  d <- d[Matrix::rowSums(d) > 0, ]
  k <- 3; a <- 0.1; e <- 0.05

  ref <- reference_gibbs(d, k, a, e, iterations = 300, seed = 11)
  ll_ref <- reference_loglik(ref$Cd, ref$Cv, ref$docs, a, e)

  set.seed(11)
  m <- tidylda(d, k = k, iterations = 3000, burnin = -1, alpha = a, eta = e,
               calc_likelihood = TRUE, verbose = FALSE)
  ll_engine <- tail(m$log_likelihood$log_likelihood, 1)

  # Loose, because these are two independent chains on a small corpus -- but far
  # tighter than a wrong target would produce. Porting text2vec's
  # beta_bar = K * eta instead of V * eta moves this by orders of magnitude more.
  expect_lt(abs(ll_engine - ll_ref) / abs(ll_ref), 0.02)
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
                 eta_in = as.matrix(et$eta), iterations = 0, burnin = -1,
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
                 eta_in = as.matrix(et$eta), iterations = 40, burnin = 15,
                 calc_likelihood = TRUE, Beta_in = p$beta_initial,
                 freeze_topics = FALSE, likelihood_every = 5, mh_steps = 1L,
                 threads = as.integer(th), verbose = FALSE)
  }

  # run() uses burnin = 15, so the averaged pair is what gets materialized and
  # Cd/Cv come back 0 x 0. Comparing those would pass no matter what the sampler
  # did, so compare the matrices that actually carry the chain's result.
  ref <- run(1)
  for (th in c(2, 4, 8)) {
    m <- run(th)
    expect_identical(ref$Cd_mean, m$Cd_mean, info = paste("threads =", th))
    expect_identical(ref$Cv_mean, m$Cv_mean, info = paste("threads =", th))
    expect_identical(ref$Ck, m$Ck, info = paste("threads =", th))
    expect_identical(ref$log_likelihood, m$log_likelihood,
                     info = paste("threads =", th))
  }

  # The raw counts are a stronger fingerprint still --- one integer per cell, no
  # averaging to blur a difference --- so check them too, on a run configured to
  # produce them.
  run_raw <- function(th) {
    set.seed(77)
    fit_lda_warp(dtm_in = d, Cd_start = p$Cd_start, alpha_in = al$alpha,
                 eta_in = as.matrix(et$eta), iterations = 40, burnin = -1,
                 calc_likelihood = TRUE, Beta_in = p$beta_initial,
                 freeze_topics = FALSE, likelihood_every = 5, mh_steps = 1L,
                 threads = as.integer(th), verbose = FALSE)
  }
  ref_raw <- run_raw(1)
  for (th in c(2, 8)) {
    m <- run_raw(th)
    expect_identical(ref_raw$Cd, m$Cd, info = paste("raw, threads =", th))
    expect_identical(ref_raw$Cv, m$Cv, info = paste("raw, threads =", th))
  }

  # And the invariants must still hold when threaded.
  m <- run_raw(8)
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
                 eta_in = as.matrix(et$eta), iterations = 0, burnin = -1,
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
      dtm_in = d, Cd_start = p$Cd_start, alpha_in = al$alpha, eta_in = as.matrix(et$eta),
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

# ---------------------------------------------------------------------------
# Collapsed joint log probability (Phase 6.6)
#
# The engine reports two quantities: a plug-in P(w | theta, beta) and the
# collapsed joint P(w, z | alpha, eta) with theta and beta integrated out. The
# joint replaced a plug-in Dirichlet density that had a sign error on both
# normalizers and, being a density, was unbounded above.
#
# Verified against an independent R implementation written from the algebra
# rather than from the C++, in the direct form with the full constants -- the
# engine accumulates the mathematically equivalent form in which zero counts
# cancel, so agreeing is a real check on that rearrangement.
# ---------------------------------------------------------------------------

collapsed_joint_reference <- function(Cd, Cv, Ck, alpha, eta) {
  D <- nrow(Cd)
  alpha_bar <- sum(alpha)
  eta_bar <- rowSums(eta)

  doc <- sum(lgamma(Cd + rep(alpha, each = D))) -
    sum(lgamma(rowSums(Cd) + alpha_bar)) +
    D * (lgamma(alpha_bar) - sum(lgamma(alpha)))

  wrd <- sum(lgamma(Cv + t(eta))) -
    sum(lgamma(Ck + eta_bar)) +
    sum(lgamma(eta_bar)) - sum(lgamma(eta))

  doc + wrd
}

fit_raw <- function(dtm, k, alpha, eta, scalar_eta, iterations = 25L) {
  init <- initialize_topic_counts(
    dtm = dtm, k = k, alpha = alpha, eta = eta, threads = 1L
  )
  fit_lda_warp(
    dtm_in = dtm,
    Cd_start = init$Cd_start,
    alpha_in = alpha,
    eta_in = if (scalar_eta) matrix(eta[1, 1], 1, 1) else eta,
    iterations = iterations,
    burnin = -1L,
    calc_likelihood = TRUE,
    Beta_in = init$beta_initial,
    freeze_topics = FALSE,
    threads = 1L,
    verbose = FALSE
  )
}

test_that("the collapsed joint matches an independent implementation", {
  dtm <- nih_sample_dtm[1:60, ]
  dtm <- dtm[, Matrix::colSums(dtm) > 0]
  V <- ncol(dtm)

  cases <- list(
    list(k = 5L, alpha = rep(0.1, 5), eta = matrix(0.05, 5, V), scalar = TRUE),
    list(k = 5L, alpha = c(.05, .1, .2, .3, .4), eta = matrix(0.05, 5, V), scalar = TRUE),
    list(k = 8L, alpha = rep(0.15, 8), eta = matrix(seq(0.01, 0.4, length.out = 8 * V), 8, V), scalar = FALSE)
  )

  for (cs in cases) {
    set.seed(42)
    raw <- fit_raw(dtm, cs$k, cs$alpha, cs$eta, cs$scalar)

    ll <- raw$log_likelihood
    expect_equal(nrow(ll), 3)

    got <- ll[3, ncol(ll)]
    want <- collapsed_joint_reference(
      raw$Cd, raw$Cv, as.numeric(raw$Ck), cs$alpha, cs$eta
    )

    # Tolerance is set by eta's float storage (D5), not by the rearrangement.
    expect_equal(got, want, tolerance = 1e-8)

    # A probability mass, not a density: it cannot be positive.
    expect_true(all(ll[3, ] < 0))

    # And the plug-in likelihood is negative too, for a different reason.
    expect_true(all(ll[2, ] < 0))
  }
})

test_that("both likelihood columns reach the model object", {
  dtm <- nih_sample_dtm[1:50, ]
  dtm <- dtm[, Matrix::colSums(dtm) > 0]

  set.seed(11)
  m <- tidylda(
    data = dtm, k = 4, iterations = 30, burnin = 10,
    calc_likelihood = TRUE, calc_r2 = FALSE, verbose = FALSE
  )

  expect_named(m$log_likelihood, c("iteration", "log_likelihood", "log_joint"))
  expect_true(all(m$log_likelihood$log_joint < 0))
  expect_true(all(is.finite(m$log_likelihood$log_joint)))

  # The two are different quantities and must not be silently the same column.
  expect_false(isTRUE(all.equal(
    m$log_likelihood$log_likelihood, m$log_likelihood$log_joint
  )))
})
