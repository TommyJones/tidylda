context("Test outputs of C++ functions directly")

### Define some common objects ----
dtm <- nih_sample_dtm

k <- 10

alpha <- rep(0.1, k)

eta <- matrix(0.05, nrow = k, ncol = ncol(dtm))

# initialize_topic_counts() returns the starting priors; since Phase 4 (D16) the
# per-token lexicon is built inside the warpLDA engine rather than handed back
# to R. fit_lda_c() -- the collapsed Gibbs sampler these tests exercise, kept
# until Phase 6 -- still wants the old lexicon, so build it explicitly here.
priors <-
  initialize_topic_counts(
    dtm = dtm,
    k = k,
    alpha = alpha,
    eta = eta,
    threads = 1
  )

counts <- create_lexicon(
  Cd_in = priors$Cd_start,
  Beta_in = priors$beta_initial,
  dtm_in = dtm,
  alpha = alpha,
  freeze_topics = FALSE
)

m <- fit_lda_c(
  Docs = counts$Docs,
  Zd_in = counts$Zd,
  eta_in = eta,
  alpha_in = alpha,
  Cd_in = counts$Cd,
  Cv_in = counts$Cv,
  Ck_in = counts$Ck,
  Beta_in = counts$Cv, # ignored
  iterations = 20,
  burnin = 10,
  freeze_topics = FALSE,
  calc_likelihood = TRUE,
  optimize_alpha = TRUE,
  verbose = FALSE
)


test_that("checksums match expectation",{
  
  sum_tokens <- sum(dtm)
  
  expect_equal(sum(m$Cd), sum_tokens)
  
  expect_equal(sum(m$Cv), sum_tokens)
  
  expect_equal(sum(m$Cd_mean), sum_tokens)
  
  expect_equal(sum(m$Cv_mean), sum_tokens)
  
  
})


test_that("optimize_alpha doesn't break anything",{
  expect_equal(sum(m$alpha), sum(alpha))
  
  expect_true(sum(is.na(rowSums(m$log_likelihood))) == 0, "log likelihood check")
})