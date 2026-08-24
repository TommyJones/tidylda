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
  # WHAT THIS DOES AND DOES NOT GUARD. An earlier version of this comment said it
  # caught a transposed Cv. It cannot, and checking showed why: every read goes
  # through counts_cv(), which detects orientation and corrects it, so handing
  # posterior() a transposed Cv produces correct output rather than wrong output.
  # The shim and this test were written in the same phase and cancel each other.
  #
  # The real orientation guard is the marginals assertion in the test above ---
  # rowSums(Cv) against colSums(dtm) --- which fails on a K x V Cv because the
  # lengths no longer match. Verified.
  #
  # What this test is still worth: an end-to-end check that posterior() attaches
  # its draws to the topic they were asked for, across the whole path from the
  # stored counts through the Dirichlet draw to the returned tibble.
  #
  # TESTED BY CORRELATION, NOT BY RANKING. Two earlier versions of this compared
  # the top 5 tokens of the posterior mean against the top 5 of beta, and both
  # were fragile for the same reason: adjacent ranks are often nearly tied --
  # ranks 5 and 6 here differ by about 11% -- so the comparison is finer than
  # anything it is comparing. The first version failed once the RNG state at this
  # point in the file changed. Seeding it and raising the draw count fixed that
  # and still failed on macOS, because the FITTED MODEL differs by platform:
  # floating-point differences flip individual MH accept decisions and the chain
  # diverges from there. That is expected of MCMC and is not a defect --- D12
  # promises reproducibility across THREAD COUNTS, not across platforms.
  #
  # Correlation tests the association directly and is indifferent to near-ties.
  # A transposed Cv would drop same-topic correlation to cross-topic levels, and
  # the gap is enormous: about 0.9995 against 0.01-0.08. The argmax assertion
  # needs no threshold at all.
  set.seed(9001)
  times <- 100
  p <- posterior(m, matrix = "beta", which = c(1, 3), times = times)

  expect_setequal(unique(p$topic), c(1, 3))
  expect_equal(nrow(p), 2 * times * ncol(d))

  totals <- tapply(p$beta, list(p$topic, p$sample), sum)
  expect_equal(as.numeric(totals), rep(1, length(totals)))

  for (k in c(1, 3)) {
    from_post <- p[p$topic == k, ]
    from_post <- tapply(from_post$beta, from_post$token, mean)
    from_post <- as.numeric(from_post[colnames(m$beta)])

    # Correlation against every topic of beta, not just its own.
    cors <- vapply(
      seq_len(nrow(m$beta)),
      function(j) stats::cor(from_post, m$beta[j, ]),
      numeric(1)
    )

    # The draws must look most like the topic they came from. Threshold-free.
    expect_equal(which.max(cors), k)

    # And the resemblance must be strong, not merely strongest.
    expect_gt(cors[k], 0.9)
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
