// warp_eta.h --- the tLDA matrix prior, laid out for the word pass
//
// Under tLDA each topic gets its own Dirichlet prior over words, so eta is a
// K x V matrix rather than a scalar. Two places in the sampler read it, both in
// the word pass and both needing only the column for the word being processed:
//
//   q_w construction   prob[k] = C_word[w][k] + eta(w,k)
//   pi_d               (C_word[w][k'] + eta(w,k')) / (C_word[w][k] + eta(w,k))
//
// So eta(w, .) must be contiguous in k, sitting beside the C_word[w] vector
// already in cache. Getting this wrong costs a cache miss per token; getting it
// right makes the matrix prior nearly free.
//
// THE LAYOUT IS FREE. R stores a K x V matrix column-major, so element (k,v)
// lives at k + v*K -- which means all K values for word v are ALREADY
// contiguous. D4 describes this as storing eta "column-major (V x K)", the
// transpose of tidylda's K x V, but that is a statement about matrix
// orientation, not about memory: no data is moved. This class copies straight
// through, downcasting to float.
//
// PRECISION (D5, D6). Stored as float, read as double. The K x V matrix is
// where the memory is, single precision halves it, and a prior does not need
// more -- at K=500, V=1e5 that is 400MB down to 200MB. at() returns double, so
// promotion happens on read and cannot be forgotten at a call site, keeping
// rounding out of the MH accept decision.
//
// eta_bar is a K-vector of DOUBLES, summed over the vocabulary from the
// double-precision input before the downcast: it is a sum of V terms and is
// added to C_k, which can exceed 2^24 where float stops representing integers
// exactly.

#pragma once

#include <vector>
#include <cstddef>

namespace warp {

class Eta {
public:
  // eta_kv points at a K x V column-major matrix of doubles (R's layout).
  Eta(const double* eta_kv, std::size_t k, std::size_t v)
      : k_(k), v_(v), e_(k * v), bar_(k, 0.0) {

    // Sum over the VOCABULARY for each topic, in double, before downcasting.
    // This is eta_bar[k] = sum_v eta[k][v] -- the normalizer of topic k's
    // Dirichlet over words. Note it is a sum over v, where alpha_bar is a sum
    // over k; conflating the two silently retargets the sampler.
    for (std::size_t w = 0; w < v_; w++) {
      const double* col = eta_kv + w * k_;
      for (std::size_t t = 0; t < k_; t++) {
        bar_[t] += col[t];
        e_[w * k_ + t] = static_cast<float>(col[t]);
      }
    }
  }

  // eta(w, k), promoted to double on read.
  double at(std::size_t w, std::size_t k) const { return e_[w * k_ + k]; }

  // Pointer to word w's contiguous K values, for the tight loops.
  const float* column(std::size_t w) const { return &e_[w * k_]; }

  double bar(std::size_t k) const { return bar_[k]; }
  const std::vector<double>& bar() const { return bar_; }

  std::size_t n_topics() const { return k_; }
  std::size_t n_words() const { return v_; }

private:
  std::size_t k_, v_;
  std::vector<float>  e_;     // K * V, word-major: e_[w*K + k]
  std::vector<double> bar_;   // K
};

}  // namespace warp
