context("tidylda core tests")

### Define some common objects ----

dtm <- textmineR::nih_sample_dtm

d1 <- dtm[1:50, ]

### Tests for initial fitting of topic models ----

test_that("can fit lda models without error", {

  # if any of the below throw an error, you've got a problem...

  # scalar priors without optimizing alpha
  lda <- tidylda(
    dtm = d1,
    k = 4,
    iterations = 20, 
    burnin = 10,
    alpha = 0.1, 
    beta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE
  )
  
  # make sure r2 doesn't have a names element
  expect_null(names(lda$r2))

  expect_length(lda$alpha, 1)

  expect_length(lda$beta, 1)

  # while we're here... check dimensions and names of objects
  expect_s3_class(lda, "tidylda")

  expect_equal(sum(dim(lda$phi) == c(4, ncol(d1))), 2)

  expect_equal(sum(dim(lda$phi) == dim(lda$gamma)), 2)

  expect_equal(sum(dim(lda$theta) == c(nrow(d1), nrow(lda$phi))), 2)

  expect_setequal(colnames(lda$phi), colnames(d1))

  expect_setequal(rownames(lda$phi), colnames(lda$theta))

  expect_setequal(rownames(lda$theta), rownames(d1))

  # scalar priors optimizing alpha
  lda <- tidylda(
    dtm = d1,
    k = 4,
    iterations = 20, 
    burnin = 10,
    alpha = 0.1, 
    beta = 0.05,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = FALSE,
    return_data = FALSE
  )

  expect_length(lda$alpha, 4)

  # vector priors
  lda <- tidylda(
    dtm = d1,
    k = 4,
    iterations = 20, 
    burnin = 10,
    alpha = rep(0.1, 4), 
    beta = rep(0.05, ncol(d1)),
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = FALSE,
    return_data = FALSE
  )

  expect_length(lda$alpha, 4)

  expect_length(lda$beta, ncol(d1))

  # beta as matrix prior
  lda <- tidylda(
    dtm = d1,
    k = 4,
    iterations = 20, 
    burnin = 10,
    alpha = 0.1, 
    beta = matrix(0.05, nrow = 4, ncol = ncol(d1)),
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE
  )

  expect_true(inherits(lda$beta, "matrix"))
})

test_that("sparse priors for beta don't cause underflow failures", {
  m <- tidylda(
    dtm = textmineR::nih_sample_dtm,
    k = 10,
    iterations = 20,
    burnin = 15,
    alpha = 0.05,
    beta = 0.01,
    optimize_alpha = FALSE,
    calc_likelihood = TRUE,
    calc_r2 = FALSE,
    return_data = FALSE
  )
  
  expect_s3_class(m, "tidylda")
})

test_that("errors hit for malformed parameters", {

  # k = 1 is bad
  expect_error(
    tidylda(
      dtm = d1,
      k = 1,
      iterations = 20, burnin = 10,
      alpha = 0.1, beta = 0.05,
      optimize_alpha = TRUE,
      calc_likelihood = FALSE,
      calc_r2 = FALSE,
      return_data = FALSE
    ),
    regexp = "k must be 2 or greater"
  )



  # burnin >= iterations
  expect_error(tidylda(
    dtm = d1,
    k = 4,
    iterations = 20, burnin = 21,
    alpha = 0.1, beta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE
  ))

  # non-numeric k
  expect_error(tidylda(
    dtm = d1,
    k = "4",
    iterations = 20, burnin = 10,
    alpha = 0.1, beta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE
  ))

  # iterations not specified
  expect_error(tidylda(
    dtm = d1,
    k = 4
  ))

  # non-logical logicals
  expect_error(tidylda(
    dtm = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, beta = 0.05,
    optimize_alpha = "FALSE",
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE
  ))

  expect_error(tidylda(
    dtm = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, beta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = "FALSE",
    calc_r2 = FALSE,
    return_data = FALSE
  ))

  expect_error(tidylda(
    dtm = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, beta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = "FALSE",
    return_data = FALSE
  ))

  expect_error(tidylda(
    dtm = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, beta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = "FALSE"
  ))
  
  expect_error(tidylda(
    dtm = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, beta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE,
    threads = nrow(d1) + 1
  ), label = "threads > nrow(dtm)")
  
  expect_warning(tidylda(
    dtm = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, beta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE,
    threads = 2
  ), label = "nrow(dtm) / threads < 100")
})

test_that("parallelism works as expected", {
  suppressWarnings(
    lda <- tidylda(
      dtm = d1,
      k = 4,
      iterations = 20, burnin = 10,
      alpha = 0.1, beta = 0.05,
      optimize_alpha = FALSE,
      calc_likelihood = TRUE,
      calc_r2 = TRUE,
      return_data = FALSE,
      threads = 2
    )
  )

  
  expect_s3_class(lda, "tidylda")
  
})
