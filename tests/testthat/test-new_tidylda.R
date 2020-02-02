context("new_tidylda tests")

# global stuff
# set up objects with the correct formatting
phi <- matrix(0, nrow = 3, ncol = 5)

gamma <- phi

theta <- matrix(0, nrow = 7, ncol = 3)

alpha <- numeric(ncol(theta))

beta <- matrix(0, nrow = nrow(phi), ncol = ncol(phi))

summary <- tibble::as_tibble(data.frame(
  topic = 1,
  prevalence = 1,
  coherence = 1,
  top_terms = "whee",
  stringsAsFactors = FALSE
))

test_that("new_tidylda works as expected with expected inputs", {
  mylda <- new_tidylda(
    phi = phi,
    theta = theta,
    gamma = gamma,
    alpha = alpha,
    beta = beta,
    summary = summary,
    call = "blah blah blah"
  )

  expect_equal(class(mylda), "tidylda")

  expect_named(mylda, c("phi", "theta", "gamma", "alpha", "beta", "summary", "call"))
})

test_that("new_tidylda throws errors for wrong class inputs", {

  # class checks
  expect_error(new_tidylda(
    phi = 1,
    theta = theta,
    gamma = gamma,
    alpha = alpha,
    beta = beta,
    summary = summary,
    call = "blah"
  ))

  expect_error(new_tidylda(
    phi = phi,
    theta = 1,
    gamma = gamma,
    alpha = alpha,
    beta = beta,
    summary = summary,
    call = "blah"
  ))

  expect_error(new_tidylda(
    phi = phi,
    theta = theta,
    gamma = 1,
    alpha = alpha,
    beta = beta,
    summary = summary,
    call = "blah"
  ))

  expect_error(new_tidylda(
    phi = phi,
    theta = theta,
    gamma = gamma,
    alpha = "alpha",
    beta = beta,
    summary = summary,
    call = "blah"
  ))

  expect_error(new_tidylda(
    phi = phi,
    theta = theta,
    gamma = gamma,
    alpha = alpha,
    beta = "beta",
    summary = summary,
    call = "blah"
  ))
})

test_that("new_tidylda throws errors for wrong dimension inputs", {

  # compatibility of dimensionality
  expect_error(new_tidylda(
    phi = t(phi),
    theta = theta,
    gamma = gamma,
    alpha = alpha,
    beta = beta,
    summary = summary,
    call = "blah"
  ))

  expect_error(new_tidylda(
    phi = phi,
    theta = theta,
    gamma = t(gamma),
    alpha = alpha,
    beta = beta,
    summary = summary,
    call = "blah"
  ))

  expect_error(new_tidylda(
    phi = phi,
    theta = theta,
    gamma = gamma,
    alpha = 1:2,
    beta = beta,
    summary = summary,
    call = "blah"
  ))

  expect_error(new_tidylda(
    phi = phi,
    theta = theta,
    gamma = gamma,
    alpha = alpha,
    beta = beta[1:2, ],
    summary = summary,
    call = "blah"
  ))

  expect_error(new_tidylda(
    phi = phi,
    theta = theta,
    gamma = gamma,
    alpha = alpha,
    beta = 1:2,
    summary = summary,
    call = "blah"
  ))
})

test_that("new_tidylda throws errors for malformed summary", {

  # wrong class
  expect_error(new_tidylda(
    phi = phi,
    theta = theta,
    gamma = gamma,
    alpha = alpha,
    beta = beta,
    summary = "whee!",
    call = "blah"
  ))


  # wrong colnames, not enough columns etc
  expect_error(new_tidylda(
    phi = phi,
    theta = theta,
    gamma = gamma,
    alpha = alpha,
    beta = beta,
    summary = summary[, 1:2],
    call = "blah"
  ))
})
