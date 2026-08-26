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
  expect_equal(unname(Matrix::rowSums(m$counts$Cv)), unname(Matrix::colSums(d)))

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
  # Matrix::t() because a 0.1.0 Cv is a dgCMatrix; the pre-0.1.0 object this
  # simulates was a dense base matrix, so coerce back to match what those models
  # actually contain.
  old <- m
  old$counts$Cv <- as.matrix(Matrix::t(m$counts$Cv))

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


test_that("Cv is labelled with the model's vocabulary and topics", {
  expect_equal(rownames(m$counts$Cv), colnames(m$beta))
  expect_equal(colnames(m$counts$Cv), rownames(m$beta))
})


test_that("Cv is sparse and Cd is dense", {
  # Cv is sparse enough to be worth storing that way (7.7% nonzero without
  # burn-in, 23.2% with it); Cd is not (38.4% and 81.3%), and at the latter a
  # dgCMatrix would be LARGER than the dense form. Measured, not assumed.
  expect_s4_class(m$counts$Cv, "dgCMatrix")
  expect_true(is.matrix(m$counts$Cd))

  # Sparsity must not have cost anything: same values, same totals.
  expect_equal(sum(m$counts$Cv), sum(d))
  expect_equal(unname(Matrix::rowSums(m$counts$Cv)), unname(Matrix::colSums(d)))
})


test_that("counts_cv() gets the square case right in both directions", {
  # k == nrow(Cv) == ncol(Cv) is the one shape that cannot disambiguate itself,
  # so counts_cv() falls back on the labels. Before 0.1.0 gave Cv dimnames that
  # fallback could never succeed, and a correctly-oriented square Cv was
  # transposed into a wrong one.
  vocab <- paste0("w", 1:4)

  fake <- list(
    beta = matrix(0, nrow = 4, ncol = 4, dimnames = list(paste0("t", 1:4), vocab))
  )

  # Current layout: tokens by topics, labelled. Must come back untouched.
  # Sparse, as a 0.1.0 model actually stores it.
  current <- Matrix::Matrix(
    matrix(1:16, nrow = 4, dimnames = list(vocab, paste0("t", 1:4))),
    sparse = TRUE
  )
  fake$counts <- list(Cv = current)
  expect_equal(counts_cv(fake), current)

  # Pre-0.1.0 layout: topics by tokens, unlabelled. Must be transposed.
  old <- matrix(1:16, nrow = 4)
  fake$counts <- list(Cv = old)
  expect_equal(counts_cv(fake), t(old))
})


test_that("counts survive a vocabulary mismatch between model and data", {
  # Transfer learning means the model's vocabulary is the union of the base
  # model's and the new data's, so nrow(Cv) exceeds ncol(new_data). Nothing may
  # key off the data's vocabulary size.
  d1 <- nih_sample_dtm[1:30, ]
  d1 <- d1[, Matrix::colSums(d1) > 0]
  d2 <- nih_sample_dtm[31:60, ]
  d2 <- d2[, Matrix::colSums(d2) > 0]

  skip_if(ncol(d1) == ncol(d2), "corpora happen to share a vocabulary size")

  set.seed(1)
  base <- tidylda(d1, k = 5, iterations = 20, burnin = 5,
                  calc_r2 = FALSE, verbose = FALSE)
  moved <- refit(base, new_data = d2, iterations = 20, burnin = 5,
                 verbose = FALSE)

  # The model's vocabulary is the union, and larger than either input.
  expect_gt(nrow(moved$counts$Cv), ncol(d2))
  expect_equal(nrow(moved$counts$Cv), ncol(moved$beta))
  expect_equal(rownames(moved$counts$Cv), colnames(moved$beta))

  # The padded DTM must stay SPARSE. Until 0.1.0 refit() aligned vocabulary with
  # a dense matrix(0, ...) filler --- 5.4 GB on a 48,508-document corpus with
  # 15,000 model-only terms, enough to exhaust a 32 GB session before sampling
  # started. Nothing about the result's dimensions or values would reveal that,
  # since cbind() returns a dgCMatrix either way; only the transient allocation
  # differed. This asserts the property directly.
  expect_s4_class(pad_vocabulary(d2, setdiff(colnames(base$beta), colnames(d2))),
                  "dgCMatrix")

  # And the downstream consumers still work against it.
  expect_equal(nrow(counts_cv(moved)), ncol(moved$beta))
  p <- posterior(moved, matrix = "beta", which = 1, times = 5)
  expect_true(all(is.finite(p$beta)))
})


test_that("pad_vocabulary appends sparse zero columns and keeps names", {
  d <- nih_sample_dtm[1:10, 1:6]
  add <- c("zzz1", "zzz2", "zzz3")

  p <- pad_vocabulary(d, add)

  expect_s4_class(p, "dgCMatrix")
  expect_equal(ncol(p), ncol(d) + length(add))
  expect_equal(colnames(p), c(colnames(d), add))
  expect_equal(rownames(p), rownames(d))

  # The appended block is empty, and the original data is untouched.
  expect_equal(sum(p[, add]), 0)
  expect_equal(sum(p), sum(d))

  # Indexing the result by name is what the callers do next.
  expect_equal(colnames(p[, c("zzz2", colnames(d)[1])]), c("zzz2", colnames(d)[1]))

  # Nothing to add is a no-op rather than an error.
  expect_identical(pad_vocabulary(d, character(0)), d)
})


test_that("eta_row_sums matches the materialized form", {
  k <- 4
  Nv <- 500

  # Matrix prior: summed directly, so exactly rowSums().
  em <- list(eta = matrix(runif(k * Nv), nrow = k, ncol = Nv))
  expect_identical(eta_row_sums(em, k, Nv), rowSums(eta_matrix(em, k, Nv)))

  # Scalar prior: agrees with the materialized form, and returns one value per
  # topic. Not asserted bit-identical --- rowSums() accumulates in long double,
  # so the two diverge around 1e-16 once Nv passes a few thousand. The engine
  # stores eta as float (D5), which is nine orders of magnitude coarser, so the
  # difference cannot reach the sampler.
  es <- list(eta = 0.05)
  expect_length(eta_row_sums(es, k, Nv), k)
  expect_equal(eta_row_sums(es, k, Nv), rowSums(eta_matrix(es, k, Nv)))
})
