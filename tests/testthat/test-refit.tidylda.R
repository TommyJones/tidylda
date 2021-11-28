context("tests of refit method for tidylda")

dtm <- nih_sample_dtm

d1 <- dtm[1:50, ]

d2 <- dtm[51:100, ]

# make sure we have different vocabulary for each data set
d1 <- d1[, Matrix::colSums(d1) > 0]

d2 <- d2[, Matrix::colSums(d2) > 0]

lda <- tidylda(
  data = d1,
  k = 4,
  iterations = 20, burnin = 10,
  alpha = 0.1, eta = 0.05,
  optimize_alpha = TRUE,
  calc_likelihood = TRUE,
  calc_r2 = TRUE,
  return_data = FALSE,
  verbose = FALSE
)



### Tests for refit.tidylda ----

test_that("can update models", {
  # continuation of the old model for another 20 iterations
  # matters because dtm lines up exactly with existing vocabulary etc.
  lda2 <- refit(
    object = lda,
    new_data = d1,
    iterations = 20,
    burnin = 10,
    prior_weight = NA,
    additional_k = 0,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE,
    verbose = FALSE
  )
  
  expect_named(lda2, names(lda))
  
  expect_equal(nrow(lda2$theta), nrow(d1))
  
  expect_equal(ncol(lda2$theta), ncol(lda$theta))
  
  expect_equal(sum(dim(lda2$beta) == dim(lda$beta)), 2)
  
  expect_equal(sum(dim(lda2$lambda) == dim(lda$lambda)), 2)
  
  
  # new data adding no extra topics no beta as prior
  lda2 <- refit(
    object = lda,
    new_data = d2,
    additional_k = 0,
    prior_weight = NA,
    iterations = 20,
    burnin = 10,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE,
    verbose = FALSE
  )
  
  expect_named(lda2, names(lda))
  
  expect_equal(nrow(lda2$theta), nrow(d2))
  
  expect_equal(ncol(lda2$theta), ncol(lda$theta))
  
  expect_equal(ncol(lda2$beta), length(union(colnames(d1), colnames(d2))))
  
  # 1 additonal topic and no beta as prior
  lda2 <- refit(
    object = lda,
    new_data = d2,
    additional_k = 1,
    prior_weight = NA,
    iterations = 20,
    burnin = 10,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE,
    verbose = FALSE
  )
  
  expect_named(lda2, names(lda))
  
  expect_equal(nrow(lda2$theta), nrow(d2))
  
  expect_equal(ncol(lda2$theta), ncol(lda$theta) + 1)
  
  expect_equal(ncol(lda2$beta), length(union(colnames(d1), colnames(d2))))
  
  # 3 additional topics and no beta as prior
  lda2 <- refit(
    object = lda,
    new_data = d2,
    additional_k = 3,
    prior_weight = NA,
    iterations = 20,
    burnin = 10,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE,
    verbose = FALSE
  )
  
  expect_named(lda2, names(lda))
  
  expect_equal(nrow(lda2$theta), nrow(d2))
  
  expect_equal(ncol(lda2$theta), ncol(lda$theta) + 3)
  
  expect_equal(ncol(lda2$beta), length(union(colnames(d1), colnames(d2))))
  
  # no additional topics and beta as prior
  lda2 <- refit(
    object = lda,
    new_data = d2,
    additional_k = 0,
    prior_weight = 1,
    iterations = 20,
    burnin = 10,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE,
    verbose = FALSE
  )
  
  expect_named(lda2, names(lda))
  
  expect_equal(nrow(lda2$theta), nrow(d2))
  
  expect_equal(ncol(lda2$theta), ncol(lda$theta))
  
  expect_equal(ncol(lda2$beta), length(union(colnames(d1), colnames(d2))))
  
  expect_true(inherits(lda2$eta, "matrix"))
  
  # 3 additonal topics and beta as prior
  lda2 <- refit(
    object = lda,
    new_data = d2,
    additional_k = 3,
    prior_weight = 1,
    iterations = 20,
    burnin = 10,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE,
    verbose = FALSE
  )
  
  expect_named(lda2, names(lda))
  
  expect_equal(nrow(lda2$theta), nrow(d2))
  
  expect_equal(ncol(lda2$theta), ncol(lda$theta) + 3)
  
  expect_equal(ncol(lda2$beta), length(union(colnames(d1), colnames(d2))))
  
  # update models with scalar eta
  
  # update models with matrix eta
  l1 <- tidylda(
    data = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1,
    eta = matrix(0.05, nrow = 4, ncol = ncol(d1)),
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE,
    verbose = FALSE
  )
  
  l2 <- refit(l1, d2, iterations = 20, verbose = FALSE)
  
  expect_equal(ncol(l2$eta), length(union(colnames(d1), colnames(d2))))
})

test_that("errors are thrown for malformed inputs to refit.tidylda", {
  
  # no vocabulary overlap between models
  nd <- rbind(numeric(10), numeric(10), numeric(10))
  
  colnames(nd) <- 1:10 # numbers means no vocab overlap
  
  lda2 <- refit(
    object = lda, 
    new_data = nd,
    iterations = 10,
    verbose = FALSE
  )
  
  expect_s3_class(lda2, "tidylda")
  
  # data doesn't have column names
  d3 <- d2
  colnames(d3) <- NULL
  expect_error(
    refit(
      object = lda,
      new_data = d3,
      iterations = 20,
      verbose = FALSE
    )
  )
  
  
  # iterations not specified
  expect_error(
    refit(
      object = lda,
      new_data = d2
    )
  )
  
  # burnin >= iterations
  expect_error(
    refit(
      object = lda,
      new_data = d2,
      additional_k = 3,
      beta_as_prior = TRUE,
      iterations = 20,
      burnin = 30,
      optimize_alpha = TRUE,
      calc_likelihood = TRUE,
      calc_r2 = TRUE,
      return_data = FALSE
    )
  )
  
  # additional_k is not numeric
  expect_error(
    refit(
      object = lda,
      new_data = d2,
      additional_k = "3",
      beta_as_prior = TRUE,
      iterations = 20,
      burnin = 10,
      optimize_alpha = TRUE,
      calc_likelihood = TRUE,
      calc_r2 = TRUE,
      return_data = FALSE
    )
  )
  
  # additional_k is less than zero
  expect_error(
    refit(
      object = lda,
      new_data = d2,
      additional_k = -3,
      beta_as_prior = TRUE,
      iterations = 20,
      burnin = 10,
      optimize_alpha = TRUE,
      calc_likelihood = TRUE,
      calc_r2 = TRUE,
      return_data = FALSE
    )
  )
  
  # iterations not specified
  expect_error(
    refit(
      object = lda,
      new_data = d2,
      additional_k = 3,
      beta_as_prior = TRUE
    )
  )
  
  # malformed prior weight
  expect_error(
    refit(
      object = lda,
      new_data = d2,
      additional_k = 3,
      prior_weight = TRUE,
      iterations = 20,
      burnin = 10,
      optimize_alpha = TRUE,
      calc_likelihood = TRUE,
      calc_r2 = TRUE,
      return_data = FALSE
    )
  )
  
  # logical things aren't logical
  expect_error(
    refit(
      object = lda,
      new_data = d2,
      additional_k = 3,
      prior_weight = 1,
      iterations = 20,
      burnin = 10,
      optimize_alpha = "TRUE",
      calc_likelihood = TRUE,
      calc_r2 = TRUE,
      return_data = FALSE
    )
  )
  
  expect_error(
    refit(
      object = lda,
      new_data = d2,
      additional_k = 3,
      beta_as_prior = TRUE,
      iterations = 20,
      burnin = 10,
      optimize_alpha = TRUE,
      calc_likelihood = "TRUE",
      calc_r2 = TRUE,
      return_data = FALSE
    )
  )
  expect_error(
    refit(
      object = lda,
      new_data = d2,
      additional_k = 3,
      beta_as_prior = TRUE,
      iterations = 20,
      burnin = 10,
      optimize_alpha = TRUE,
      calc_likelihood = TRUE,
      calc_r2 = "TRUE",
      return_data = FALSE
    )
  )
  
  expect_error(
    refit(
      object = lda,
      new_data = d2,
      additional_k = 3,
      beta_as_prior = TRUE,
      iterations = 20,
      burnin = 10,
      optimize_alpha = TRUE,
      calc_likelihood = TRUE,
      calc_r2 = TRUE,
      return_data = "FALSE"
    )
  )
})
