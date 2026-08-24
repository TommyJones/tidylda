// warp_init_sample.h --- drawing each token's initial topic
//
// Initialization samples one topic per token from
//
//     P(z = k) proportional to beta_hat[k][v] * (C^d_dk + alpha_k)
//
// held fixed within a document (D8, the informed initialization). The weights
// arrive as unnormalized LOG probabilities, one K-vector per distinct (d,v)
// pair, and are then drawn from once per repeat of that word.
//
// WHY THIS REPLACED lsamp_one(). That function (src/sample_int.h) drew a single
// variate by sorting the whole K-vector descending -- twice, once for the values
// and once for the index permutation -- then accumulating log-sum-exp across it
// with a log_add_exp per topic, each of which is an exp plus a log1p. Roughly
// three heap allocations and O(K log K) work per token, to produce one integer.
// On this package's medium benchmark corpus that was 3.90 s at K = 50: 4% of a
// 1200-iteration fit, but 20% of a 200-iteration one, and it grows with both N
// and K log K.
//
// The sort was not there for correctness. Accumulating in descending order is a
// numerical-stability device, and subtracting the maximum is the standard and
// stronger one -- it bounds every exponent at zero regardless of order, which is
// what matters when a transferred tLDA prior makes one topic dominate its
// column by many orders of magnitude.

#pragma once

#include <vector>
#include <cstddef>
#include <cmath>

namespace warp {

// Draw one category from unnormalized log weights.
//
// RNG-agnostic, following the convention sample.h established and warp_alias.h
// continues: the uniform variate is passed in rather than drawn, so this cannot
// reach for a generator that a caller did not intend.
//
//   lp       unnormalized log weights, length K
//   u        a U(0,1) variate
//   scratch  caller-owned buffer of length K, reused across calls
//
// One pass to find the maximum, one to exponentiate and accumulate, one to
// locate u. No allocation.
inline std::size_t sample_log_weights(const std::vector<double>& lp,
                                      double u,
                                      std::vector<double>& scratch) {
  const std::size_t K = lp.size();

  double max_lp = lp[0];
  for (std::size_t k = 1; k < K; k++) {
    if (lp[k] > max_lp) max_lp = lp[k];
  }

  // Every exponent is <= 0, so nothing overflows and the largest term is
  // exactly 1. Underflow of the small terms is harmless: they contribute
  // nothing to the total and cannot be selected.
  double total = 0.0;
  for (std::size_t k = 0; k < K; k++) {
    total += std::exp(lp[k] - max_lp);
    scratch[k] = total;
  }

  const double target = u * total;
  for (std::size_t k = 0; k < K; k++) {
    if (scratch[k] >= target) return k;
  }
  return K - 1;  // guards u arbitrarily close to 1 and rounding in the sum
}

}  // namespace warp
