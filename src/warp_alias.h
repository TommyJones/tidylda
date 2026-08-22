// warp_alias.h --- Vose alias table for O(1) categorical sampling
//
// warpLDA draws one categorical variate per token from a distribution that is
// fixed for the whole word (the word-proposal q_w) or for the whole run (the
// prior alpha, per D19). An alias table pays O(K) once to make every subsequent
// draw O(1): one indexed load and one comparison.
//
// RNG-AGNOSTIC BY CONSTRUCTION. Like sample_one() in sample.h, sample() takes
// its uniform variates as arguments rather than drawing them. This is what
// makes D12's "R's RNG only on the main thread" rule fall out structurally
// instead of requiring discipline at every call site -- the table has no
// generator to misuse.
//
// ALLOCATION. setup() reuses the object's buffers, so a single AliasTable
// hoisted out of the per-word loop costs zero allocations per word. The
// reference constructs a fresh AliasUrn and grows a fresh vector<double> with
// push_back for every word of every iteration (LDA.hpp:132-136), which is
// V * iterations heap allocations; design notes section 4.3 identifies that as
// its real inefficiency, ahead of the dense O(K) loop itself.

#pragma once

#include <vector>
#include <cstddef>

namespace warp {

class AliasTable {
public:
  // Build the table from unnormalized, non-negative weights.
  //
  // Vose's construction: scale the weights so they average 1, then repeatedly
  // pair a light entry with a heavy one until every entry is a two-outcome
  // mixture. O(K) time, O(K) space.
  void setup(const std::vector<double>& w) {
    const std::size_t k = w.size();
    prob_.resize(k);
    alias_.resize(k);
    small_.clear();
    large_.clear();

    double sum = 0.0;
    for (std::size_t i = 0; i < k; i++) sum += w[i];

    if (k == 0) return;

    // Defensive only: every caller passes C + eta with eta > 0, so the sum is
    // positive in practice. Falling back to uniform keeps a degenerate input
    // from producing NaN weights that would silently corrupt a chain.
    if (!(sum > 0.0)) {
      for (std::size_t i = 0; i < k; i++) { prob_[i] = 1.0; alias_[i] = i; }
      return;
    }

    const double scale = static_cast<double>(k) / sum;
    for (std::size_t i = 0; i < k; i++) {
      prob_[i] = w[i] * scale;                      // mean 1 across entries
      (prob_[i] < 1.0 ? small_ : large_).push_back(i);
    }

    while (!small_.empty() && !large_.empty()) {
      const std::size_t s = small_.back(); small_.pop_back();
      const std::size_t l = large_.back(); large_.pop_back();

      alias_[s] = l;
      prob_[l] = (prob_[l] + prob_[s]) - 1.0;       // donate s's deficit to l

      (prob_[l] < 1.0 ? small_ : large_).push_back(l);
    }

    // Anything still queued should be exactly 1 but may be a rounding step
    // away. Pinning it to 1 makes those entries deterministic rather than
    // very-nearly-deterministic.
    for (std::size_t i : large_) { prob_[i] = 1.0; alias_[i] = i; }
    for (std::size_t i : small_) { prob_[i] = 1.0; alias_[i] = i; }
    small_.clear();
    large_.clear();
  }

  // Draw one index. u1 picks the entry, u2 picks within it. Both U(0,1).
  std::size_t sample(double u1, double u2) const {
    const std::size_t k = prob_.size();
    std::size_t i = static_cast<std::size_t>(u1 * static_cast<double>(k));
    if (i >= k) i = k - 1;                          // guards u1 == 1 - epsilon
    return u2 < prob_[i] ? i : alias_[i];
  }

  std::size_t size() const { return prob_.size(); }

private:
  std::vector<double>      prob_;
  std::vector<std::size_t> alias_;
  std::vector<std::size_t> small_;   // scratch, retained to avoid reallocation
  std::vector<std::size_t> large_;
};

}  // namespace warp
