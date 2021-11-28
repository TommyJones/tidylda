context("tidylda core tests")

### Define some common objects ----

dtm <- nih_sample_dtm

d1 <- dtm[1:50, ]

### Tests for initial fitting of topic models ----

test_that("can fit lda models without error", {

  # if any of the below throw an error, you've got a problem...

  # scalar priors without optimizing alpha
  lda <- tidylda(
    data = d1,
    k = 4,
    iterations = 20, 
    burnin = 10,
    alpha = 0.1, 
    eta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE,
    verbose = FALSE
  )
  
  # make sure r2 is numeric since calc_r2 = TRUE
  expect_type(lda$r2, "double")
  
  # make sure r2 doesn't have a names element
  expect_null(names(lda$r2))

  # make sure that likelihood is correct since calc_likelihood = TRUE
  expect_s3_class(lda$log_likelihood, "tbl_df")
  
  expect_equal(ncol(lda$log_likelihood), 2)
  
  expect_equal(nrow(lda$log_likelihood), tail(lda$log_likelihood$iteration, 1) + 1)
  
  # while we're here... check dimensions and names of objects
  expect_s3_class(lda, "tidylda")

  expect_length(lda$alpha, 1)
  
  expect_length(lda$eta, 1)
  
  expect_equal(sum(dim(lda$beta) == c(4, ncol(d1))), 2)

  expect_equal(sum(dim(lda$beta) == dim(lda$lambda)), 2)

  expect_equal(sum(dim(lda$theta) == c(nrow(d1), nrow(lda$beta))), 2)

  expect_setequal(colnames(lda$beta), colnames(d1))

  expect_setequal(rownames(lda$beta), colnames(lda$theta))

  expect_setequal(rownames(lda$theta), rownames(d1))

  # scalar priors optimizing alpha
  lda <- tidylda(
    data = d1,
    k = 4,
    iterations = 20, 
    burnin = 10,
    alpha = 0.1, 
    eta = 0.05,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  )

  expect_length(lda$alpha, 4)

  # vector priors
  lda <- tidylda(
    data = d1,
    k = 4,
    iterations = 20, 
    burnin = 10,
    alpha = rep(0.1, 4), 
    eta = rep(0.05, ncol(d1)),
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  )

  expect_length(lda$alpha, 4)

  expect_length(lda$eta, ncol(d1))

  # eta as matrix prior
  lda <- tidylda(
    data = d1,
    k = 4,
    iterations = 20, 
    burnin = 10,
    alpha = 0.1, 
    eta = matrix(0.05, nrow = 4, ncol = ncol(d1)),
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  )

  expect_true(inherits(lda$eta, "matrix"))
})

test_that("sparse priors for eta don't cause underflow failures", {
  m <- tidylda(
    data = nih_sample_dtm,
    k = 10,
    iterations = 20,
    burnin = 15,
    alpha = 0.05,
    eta = 0.01,
    optimize_alpha = FALSE,
    calc_likelihood = TRUE,
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  )
  
  expect_s3_class(m, "tidylda")
})

test_that("errors hit for malformed parameters", {

  # k = 1 is bad
  expect_error(
    tidylda(
      data = d1,
      k = 1,
      iterations = 20, burnin = 10,
      alpha = 0.1, eta = 0.05,
      optimize_alpha = TRUE,
      calc_likelihood = FALSE,
      calc_r2 = FALSE,
      return_data = FALSE,
      verbose = FALSE
    ),
    regexp = "k must be 2 or greater"
  )

  # iterations not specified
  expect_error(
    tidylda(
      data = d1,
      k = 10,
      alpha = 0.1, eta = 0.05,
      optimize_alpha = TRUE,
      calc_likelihood = FALSE,
      calc_r2 = FALSE,
      return_data = FALSE,
      verbose = FALSE
    ),
    label = "iterations not specified"
  )
  

  # burnin >= iterations
  expect_error(tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 21,
    alpha = 0.1, eta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  ))

  # non-numeric k
  expect_error(tidylda(
    data = d1,
    k = "4",
    iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  ))

  # iterations not specified
  expect_error(tidylda(
    data = d1,
    k = 4
  ))

  # non-logical logicals
  expect_error(tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05,
    optimize_alpha = "FALSE",
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  ))

  expect_error(tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = "FALSE",
    calc_r2 = FALSE,
    return_data = FALSE,
    verbose = FALSE
  ))

  expect_error(tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = "FALSE",
    return_data = FALSE
  ))

  expect_error(tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = "FALSE"
  ))
  
  expect_error(tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE,
    threads = nrow(d1) + 1
  ), label = "threads > nrow(dtm)")
  
  expect_warning(tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1, eta = 0.05,
    optimize_alpha = FALSE,
    calc_likelihood = FALSE,
    calc_r2 = FALSE,
    return_data = FALSE,
    threads = 2,
    verbose = FALSE
  ), label = "nrow(dtm) / threads < 100")
  
  # data doesn't have column names
  d3 <- d1
  colnames(d3) <- NULL
  expect_error(
    tidylda(
      data = d3,
      k = 4,
      iterations = 20
    )
  )
  
})

# note as of this writing, not parallel, 
# but use of threads argument should not throw errors
test_that("parallelism works as expected", {
  suppressWarnings(
    lda <- tidylda(
      data = d1,
      k = 4,
      iterations = 20, burnin = 10,
      alpha = 0.1, eta = 0.05,
      optimize_alpha = FALSE,
      calc_likelihood = TRUE,
      calc_r2 = TRUE,
      return_data = FALSE,
      threads = 2,
      verbose = FALSE
    )
  )

  
  expect_s3_class(lda, "tidylda")
  
})
