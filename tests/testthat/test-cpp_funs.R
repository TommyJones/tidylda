context("Test outputs of C++ functions directly")

### Define some common objects ----
dtm <- textmineR::nih_sample_dtm

k <- 10

alpha <- rep(0.1, k)

beta <- matrix(0.05, nrow = k, ncol = ncol(dtm))

counts <- 
  initialize_topic_counts(
    dtm = dtm, 
    k = 10,
    alpha = rep(0.1, 10), 
    beta = matrix(0.05, nrow = 10, ncol = ncol(dtm)),
    threads = 1
  )

m <- fit_lda_c(
  Docs = counts$Docs,
  Zd_in = counts$Zd,
  beta_in = beta,
  alpha_in = alpha,
  Cd_in = counts$Cd,
  Cv_in = counts$Cv,
  Ck_in = counts$Ck,
  Phi = counts$Cv, # ignored
  iterations = 200,
  burnin = 175,
  freeze_topics = FALSE,
  calc_likelihood = TRUE,
  optimize_alpha = TRUE
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