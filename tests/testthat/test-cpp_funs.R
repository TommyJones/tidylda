context("Test outputs of C++ functions directly")

# These checksums were written for fit_lda_c(), the collapsed Gibbs sampler
# Phase 6 removed. They are invariants of any correct LDA sampler rather than of
# that implementation, so they were re-pointed at the warpLDA engine rather than
# retired with it. Two matter more now than when they were written:
#
#   * the log-likelihood finiteness check is covered nowhere else;
#   * "alpha comes back unchanged" tested optimize_alpha's behaviour, and now
#     tests its removal under D7.

dtm <- nih_sample_dtm

k <- 10

alpha <- rep(0.1, k)

eta <- matrix(0.05, nrow = k, ncol = ncol(dtm))

priors <- initialize_topic_counts(
  dtm = dtm,
  k = k,
  alpha = alpha,
  eta = eta,
  threads = 1
)

set.seed(90210)

m <- fit_lda_warp(
  dtm_in = dtm,
  Cd_start = priors$Cd_start,
  alpha_in = alpha,
  eta_in = eta,
  iterations = 20,
  burnin = 10,
  calc_likelihood = TRUE,
  Beta_in = priors$beta_initial,
  freeze_topics = FALSE,
  likelihood_every = 1,
  mh_steps = 1L,
  threads = 1L,
  verbose = FALSE
)

sum_tokens <- sum(dtm)


test_that("checksums match expectation", {

  # Every token is assigned to exactly one topic, so all four views of the count
  # structure total N. Cd_mean and Cv_mean are post-burn-in averages, and an
  # average of quantities that each sum to N sums to N as well.
  expect_equal(sum(m$Cd), sum_tokens)

  expect_equal(sum(m$Cv), sum_tokens)

  expect_equal(sum(m$Cd_mean), sum_tokens)

  expect_equal(sum(m$Cv_mean), sum_tokens)

  # Ck is the topic marginal of both matrices.
  expect_equal(as.numeric(m$Ck), unname(colSums(m$Cd)))

  expect_equal(as.numeric(m$Ck), unname(colSums(m$Cv)))
})


test_that("the log likelihood is finite throughout", {

  # Underflow in beta or a zero-probability token would surface here as NaN or
  # -Inf long before it showed up in the fitted model.
  expect_equal(nrow(m$log_likelihood), 3)

  expect_true(all(is.finite(m$log_likelihood)))

  # Row 2 is the corpus log probability without the priors, and is negative.
  expect_true(all(m$log_likelihood[2, ] < 0))
})


test_that("alpha is returned unchanged", {

  # This asserted that optimize_alpha preserved the total. D7 removed
  # optimize_alpha entirely -- alpha is fixed for the whole run now -- so the
  # check is stronger: not merely the same sum, the same vector.
  expect_equal(as.numeric(m$alpha), alpha)
})
