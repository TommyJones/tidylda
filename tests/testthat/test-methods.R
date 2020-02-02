context("tests of methods for tidylda")

dtm <- textmineR::nih_sample_dtm

d1 <- dtm[1:50, ]

d2 <- dtm[51:100, ]

# make sure we have different vocabulary for each data set
d1 <- d1[, Matrix::colSums(d1) > 0]

d2 <- d2[, Matrix::colSums(d2) > 0]

lda <- tidylda(
  dtm = d1,
  k = 4,
  iterations = 20, burnin = 10,
  alpha = 0.1, beta = 0.05,
  optimize_alpha = TRUE,
  calc_likelihood = TRUE,
  calc_r2 = TRUE,
  return_data = FALSE
)

### Tests for predictions ----

test_that("can make predictions without error", {
  # one row gibbs with burnin
  p <- predict(object = lda, new_data = d2[1, ], method = "gibbs", iterations = 20, burnin = 10)

  expect_equal(nrow(p), 1)

  expect_equal(ncol(p), ncol(lda$theta))

  expect_setequal(colnames(p), colnames(lda$theta))

  # multi-row gibbs with burnin
  p <- predict(object = lda, new_data = d2, method = "gibbs", iterations = 20, burnin = 10)

  expect_equal(nrow(p), nrow(d2))

  expect_equal(ncol(p), ncol(lda$theta))

  expect_setequal(colnames(p), colnames(lda$theta))

  # single row dot method
  p <- predict(object = lda, new_data = d2[1, ], method = "dot")

  expect_equal(nrow(p), 1)

  expect_equal(ncol(p), ncol(lda$theta))

  expect_setequal(colnames(p), colnames(lda$theta))

  # multi-row dot method
  p <- predict(object = lda, new_data = d2, method = "dot")

  expect_equal(nrow(p), nrow(d2))

  expect_equal(ncol(p), ncol(lda$theta))

  expect_setequal(colnames(p), colnames(lda$theta))
})

test_that("malformed args in predict throw errors", {

  # no iterations specified
  expect_error(
    predict(object = lda, new_data = d2, method = "gibbs")
  )

  # burnin >= iterations
  expect_error(
    predict(object = lda, new_data = d2, method = "gibbs", iterations = 5, burnin = 6)
  )

  # incorrect method
  expect_error(
    predict(object = lda, new_data = d2, method = "oops")
  )
})


### Tests for updates ----

test_that("can update models", {
  # continuation of the old model for another 20 iterations
  # matters because dtm lines up exactly with existing vocabulary etc.
  lda2 <- update(
    object = lda,
    dtm = d1,
    additional_k = 0,
    phi_as_prior = FALSE,
    iterations = 20,
    burnin = 10,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE
  )

  expect_named(lda2, names(lda))

  expect_equal(nrow(lda2$theta), nrow(d1))

  expect_equal(ncol(lda2$theta), ncol(lda$theta))

  expect_equal(sum(dim(lda2$phi) == dim(lda$phi)), 2)

  expect_equal(sum(dim(lda2$gamma) == dim(lda$gamma)), 2)


  # new data adding no extra topics no phi as prior
  lda2 <- update(
    object = lda,
    dtm = d2,
    additional_k = 0,
    phi_as_prior = FALSE,
    iterations = 20,
    burnin = 10,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE
  )

  expect_named(lda2, names(lda))

  expect_equal(nrow(lda2$theta), nrow(d2))

  expect_equal(ncol(lda2$theta), ncol(lda$theta))

  expect_equal(ncol(lda2$phi), length(union(colnames(d1), colnames(d2))))

  # 1 additonal topic and no phi as prior
  lda2 <- update(
    object = lda,
    dtm = d2,
    additional_k = 1,
    phi_as_prior = FALSE,
    iterations = 20,
    burnin = 10,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE
  )

  expect_named(lda2, names(lda))

  expect_equal(nrow(lda2$theta), nrow(d2))

  expect_equal(ncol(lda2$theta), ncol(lda$theta) + 1)

  expect_equal(ncol(lda2$phi), length(union(colnames(d1), colnames(d2))))

  # 3 additional topics and no phi as prior
  lda2 <- update(
    object = lda,
    dtm = d2,
    additional_k = 3,
    phi_as_prior = FALSE,
    iterations = 20,
    burnin = 10,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE
  )

  expect_named(lda2, names(lda))

  expect_equal(nrow(lda2$theta), nrow(d2))

  expect_equal(ncol(lda2$theta), ncol(lda$theta) + 3)

  expect_equal(ncol(lda2$phi), length(union(colnames(d1), colnames(d2))))

  # no additional topics and phi as prior
  lda2 <- update(
    object = lda,
    dtm = d2,
    additional_k = 0,
    phi_as_prior = TRUE,
    iterations = 20,
    burnin = 10,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE
  )

  expect_named(lda2, names(lda))

  expect_equal(nrow(lda2$theta), nrow(d2))

  expect_equal(ncol(lda2$theta), ncol(lda$theta))

  expect_equal(ncol(lda2$phi), length(union(colnames(d1), colnames(d2))))


  # 3 additonal topics and phi as prior
  lda2 <- update(
    object = lda,
    dtm = d2,
    additional_k = 3,
    phi_as_prior = TRUE,
    iterations = 20,
    burnin = 10,
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE
  )

  expect_named(lda2, names(lda))

  expect_equal(nrow(lda2$theta), nrow(d2))

  expect_equal(ncol(lda2$theta), ncol(lda$theta) + 3)

  expect_equal(ncol(lda2$phi), length(union(colnames(d1), colnames(d2))))

  # update models with scalar beta

  # update models with matrix beta
  l1 <- tidylda(
    dtm = d1,
    k = 4,
    iterations = 20, burnin = 10,
    alpha = 0.1,
    beta = matrix(0.05, nrow = 4, ncol = ncol(d1)),
    optimize_alpha = TRUE,
    calc_likelihood = TRUE,
    calc_r2 = TRUE,
    return_data = FALSE
  )

  l2 <- update(l1, d2, iterations = 20)

  expect_equal(ncol(l2$beta), length(union(colnames(d1), colnames(d2))))
})

test_that("errors are thrown for malformed inputs to update.tidylda", {

  # burnin >= iterations
  expect_error(
    update(
      object = lda,
      dtm = d2,
      additional_k = 3,
      phi_as_prior = TRUE,
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
    update(
      object = lda,
      dtm = d2,
      additional_k = "3",
      phi_as_prior = TRUE,
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
    update(
      object = lda,
      dtm = d2,
      additional_k = -3,
      phi_as_prior = TRUE,
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
    update(
      object = lda,
      dtm = d2,
      additional_k = 3,
      phi_as_prior = TRUE
    )
  )

  # logical things aren't logical
  expect_error(
    update(
      object = lda,
      dtm = d2,
      additional_k = 3,
      phi_as_prior = "TRUE",
      iterations = 20,
      burnin = 10,
      optimize_alpha = TRUE,
      calc_likelihood = TRUE,
      calc_r2 = TRUE,
      return_data = FALSE
    )
  )

  expect_error(
    update(
      object = lda,
      dtm = d2,
      additional_k = 3,
      phi_as_prior = TRUE,
      iterations = 20,
      burnin = 10,
      optimize_alpha = "TRUE",
      calc_likelihood = TRUE,
      calc_r2 = TRUE,
      return_data = FALSE
    )
  )

  expect_error(
    update(
      object = lda,
      dtm = d2,
      additional_k = 3,
      phi_as_prior = TRUE,
      iterations = 20,
      burnin = 10,
      optimize_alpha = TRUE,
      calc_likelihood = "TRUE",
      calc_r2 = TRUE,
      return_data = FALSE
    )
  )
  expect_error(
    update(
      object = lda,
      dtm = d2,
      additional_k = 3,
      phi_as_prior = TRUE,
      iterations = 20,
      burnin = 10,
      optimize_alpha = TRUE,
      calc_likelihood = TRUE,
      calc_r2 = "TRUE",
      return_data = FALSE
    )
  )

  expect_error(
    update(
      object = lda,
      dtm = d2,
      additional_k = 3,
      phi_as_prior = TRUE,
      iterations = 20,
      burnin = 10,
      optimize_alpha = TRUE,
      calc_likelihood = TRUE,
      calc_r2 = TRUE,
      return_data = "FALSE"
    )
  )
})

### Tests for the print method ----

test_that("print.tidylda behaves as expected", {

  # no error
  print(lda)

  # assignment creates a new object of class tidylda
  m <- print(lda)

  expect_true("tidylda" %in% class(m))

  expect_named(m, names(lda))

  # can modify digits
  m2 <- tidylda(dtm = d1, k = 5, iterations = 20, calc_r2 = TRUE)

  print(m2, digits = 2)
})

### tests for the glance method ----

test_that("glance.tidylda behaves nicely", {

  # well-formed call
  g <- glance(lda)

  expect_named(g, c(
    "num_topics", "num_documents", "num_tokens",
    "iterations", "burnin"
  ))

  expect_equal(g$num_topics, nrow(lda$phi))

  expect_equal(g$num_documents, nrow(lda$theta))

  expect_equal(g$num_tokens, ncol(lda$phi))

  expect_equal(g$iterations, lda$call$iterations)

  expect_equal(g$burnin, lda$call$burnin)

  # malformed call
  n <- new_tidylda(
    phi = matrix(0, nrow = 2, ncol = 2),
    theta = matrix(0, nrow = 2, ncol = 2),
    gamma = matrix(0, nrow = 2, ncol = 2),
    alpha = 0,
    beta = 0,
    summary = data.frame(
      topic = 0,
      prevalence = 0,
      coherence = 0,
      top_terms = "0",
      stringsAsFactors = FALSE
    ),
    call = "whee"
  )

  g <- glance(n)

  expect_named(g, c(
    "num_topics", "num_documents", "num_tokens",
    "iterations", "burnin"
  ))

  expect_equal(g$iterations, NA)

  expect_equal(g$burnin, NA)
})

test_that("glance works with updated models", {
  l2 <- update(lda, d2, iterations = 20)

  g <- glance(l2)

  expect_named(g, c(
    "num_topics", "num_documents", "num_tokens",
    "iterations", "burnin"
  ))

  expect_equal(g$num_topics, nrow(l2$phi))

  expect_equal(g$num_documents, nrow(l2$theta))

  expect_equal(g$num_tokens, ncol(l2$phi))

  expect_equal(g$iterations, l2$call$iterations)

  expect_equal(g$burnin, NA)
})

### test tidy methods ----

test_that("tidy.tidylda works as expected", {

  # tidy phi
  tidy_phi <- tidy(
    x = lda,
    matrix = "phi"
  )

  expect_named(tidy_phi, c("topic", "token", "phi"))

  expect_type(tidy_phi[[1]], "double")

  expect_type(tidy_phi[[2]], "character")

  expect_type(tidy_phi[[3]], "double")

  # log to tidy phi
  tidy_phi_log <- tidy(
    x = lda,
    matrix = "phi",
    log = TRUE
  )

  expect_named(tidy_phi_log, c("topic", "token", "log_phi"))

  expect_type(tidy_phi_log[[1]], "double")

  expect_type(tidy_phi_log[[2]], "character")

  expect_type(tidy_phi_log[[3]], "double")

  # tidy theta
  tidy_theta <- tidy(
    x = lda,
    matrix = "theta"
  )

  expect_named(tidy_theta, c("document", "topic", "theta"))

  expect_type(tidy_theta[[1]], "character")

  expect_type(tidy_theta[[2]], "double")

  expect_type(tidy_theta[[3]], "double")

  # tidy gamma
  tidy_gamma <- tidy(
    x = lda,
    matrix = "gamma"
  )

  expect_named(tidy_gamma, c("topic", "token", "gamma"))

  expect_type(tidy_gamma[[1]], "double")

  expect_type(tidy_gamma[[2]], "character")

  expect_type(tidy_gamma[[3]], "double")
})

test_that("tidy.tidylda throws errors for malformed inputs", {
  expect_error(
    tidy(
      x = lda,
      matrix = "WRONG"
    )
  )

  expect_error(
    tidy(
      x = lda,
      matrix = 1
    )
  )

  expect_error(
    tidy(
      x = lda,
      matrix = "phi",
      log = "WRONG"
    )
  )
})
