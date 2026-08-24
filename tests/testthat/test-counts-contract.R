context("counts orientation and backward compatibility")

# D17 changed the exported counts to <major>-by-topics: Cd documents-by-topics,
# Cv words-by-topics, matching the orientation the engine holds them in. `beta`
# and `theta` are unchanged.
#
# These tests exist because a wrong orientation mostly fails quietly. The token
# labels still line up, the numbers are still probabilities, and only the
# topic-to-token association is scrambled -- which no dimension check would see.

d <- nih_sample_dtm[1:40, ]

set.seed(1)
m <- tidylda(d, k = 4, iterations = 30, burnin = 10,
             calc_likelihood = FALSE, verbose = FALSE)


test_that("counts are stored major-by-topics, beta and theta unchanged", {
  expect_equal(dim(m$counts$Cd), c(nrow(d), 4))
  expect_equal(dim(m$counts$Cv), c(ncol(d), 4))

  # The public matrices keep their historical orientation.
  expect_equal(dim(m$beta), c(4, ncol(d)))
  expect_equal(dim(m$theta), c(nrow(d), 4))

  # Marginals tie each matrix to the axis it is indexed by, which pins the
  # orientation rather than merely the shape.
  expect_equal(unname(rowSums(m$counts$Cd)), unname(Matrix::rowSums(d)))
  expect_equal(unname(rowSums(m$counts$Cv)), unname(Matrix::colSums(d)))

  # Totals are exact in both. Note the per-TOPIC marginals are not required to
  # agree between them: D10 accumulates Cd during the doc pass and Cv during the
  # word pass, each at the moment its own matrix is current, so they estimate
  # marginal expectations at different points in the cycle. Only the totals are
  # guaranteed, and asserting more than that fails against correct output.
  expect_equal(sum(m$counts$Cd), sum(d))
  expect_equal(sum(m$counts$Cv), sum(d))
})


test_that("posterior() draws associate tokens with the right topic", {
  # The check that matters: a transposed Cv would still give valid probability
  # vectors, just attached to the wrong topics.
  #
  # SEEDED, AND WITH ENOUGH DRAWS TO MEAN SOMETHING. This compares a Monte Carlo
  # mean against the exact one, so `times` sets its resolution. At the original
  # times = 10 the comparison was finer than the estimate: ranks 5 and 6 of beta
  # differ by about 11% here (0.0081 against 0.0072), which 10 draws reorder
  # easily. It passed only because the RNG happened to be in a favourable state
  # at this point in the file, and Phase 7's changed draw sequence exposed that.
  #
  # 200 draws resolve the gap, and the seed stops the result depending on
  # whatever ran before. The assertion is unchanged.
  set.seed(9001)
  p <- posterior(m, matrix = "beta", which = c(1, 3), times = 200)

  expect_setequal(unique(p$topic), c(1, 3))
  expect_equal(nrow(p), 2 * 200 * ncol(d))

  totals <- tapply(p$beta, list(p$topic, p$sample), sum)
  expect_equal(as.numeric(totals), rep(1, length(totals)))

  for (k in c(1, 3)) {
    from_post <- p[p$topic == k, ]
    from_post <- tapply(from_post$beta, from_post$token, mean)
    top_post <- names(sort(from_post, decreasing = TRUE))[1:5]
    top_beta <- names(sort(m$beta[k, ], decreasing = TRUE))[1:5]
    expect_setequal(top_post, top_beta)
  }
})


test_that("models saved with the old counts orientation still work", {
  # tidylda is on CRAN, so saved models carry the pre-D17 topics-by-words Cv.
  # refit() reads it on the tLDA critical path, where a wrong per-topic weight
  # would corrupt eta rather than error.
  old <- m
  old$counts$Cv <- t(m$counts$Cv)

  expect_equal(dim(old$counts$Cv), c(4, ncol(d)))
  expect_equal(dim(counts_cv(old)), dim(m$counts$Cv))
  expect_equal(unname(as.matrix(counts_cv(old))), unname(as.matrix(m$counts$Cv)))

  # And the new-format object is passed through untouched.
  expect_equal(unname(as.matrix(counts_cv(m))), unname(as.matrix(m$counts$Cv)))

  r_new <- refit(m, d, iterations = 15, prior_weight = 1, verbose = FALSE)
  r_old <- refit(old, d, iterations = 15, prior_weight = 1, verbose = FALSE)
  expect_equal(r_new$eta, r_old$eta)

  expect_equal(nrow(posterior(old, matrix = "beta", which = 2, times = 5)),
               5 * ncol(d))
})
