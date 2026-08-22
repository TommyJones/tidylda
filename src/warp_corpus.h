// warp_corpus.h --- the dual-view token array warpLDA samples over
//
// warpLDA's efficiency comes from visiting the same tokens twice per iteration
// through two different orderings, so that each visit's working set is a single
// K-vector that stays in cache:
//
//   doc pass   tokens grouped by document -> working set is C_doc[d], a K-vector
//   word pass  tokens grouped by word     -> working set is C_word[w], a K-vector
//
// This class holds one token array plus CSR (by document) and CSC (by word)
// index structures over it. The doc pass walks the array contiguously; the word
// pass indirects through csc_token_.
//
// TOKEN LAYOUT. One flat uint16 array with stride (1 + mh_steps):
//
//   z_[i*stride]         old_z, the token's current topic
//   z_[i*stride + 1 + m] proposal m, pending resolution by the other pass
//
// At the default mh_steps = 1 that is 4 bytes per token with both fields on one
// cache line, and D18's "mh_steps * 2 bytes per token above the default" falls
// out of the stride. Interleaving rather than using parallel arrays matters most
// in the word pass, where access is scattered and a split layout would cost two
// cache misses per token instead of one.
//
// TOKENS ARE SORTED BY WORD WITHIN EACH DOCUMENT. Three things fall out:
// consecutive tokens in the doc pass often share a word; the CSC build is a
// counting sort; and the log-likelihood can iterate distinct (d,v) pairs by
// scanning runs, which is the O(nnz*K) rather than O(N*K) form that design
// notes section 7 requires. Reordering is safe -- the doc-proposal position
// branch draws a uniformly random token in the document, and the accept sweep
// is valid in any order.

#pragma once

#include <vector>
#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <numeric>

namespace warp {

typedef uint16_t topic_t;   // K < 65536; checked by the caller
typedef uint32_t word_t;
typedef uint32_t doc_t;
typedef std::size_t token_t;

class Corpus {
public:
  // Build from tidylda's create_lexicon() output.
  //
  //   docs[d]  word indices of document d, with repeats
  //   zd[d]    initial topic of each of those tokens (informed init, D8)
  //
  // Phase 4 (D16) replaces this constructor with a direct DTM walk, dropping
  // the 16-bytes-per-token vector<vector<size_t>> round trip through R.
  Corpus(const std::vector<std::vector<std::size_t>>& docs,
         const std::vector<std::vector<std::size_t>>& zd,
         std::size_t n_words,
         std::size_t mh_steps)
      : n_docs_(docs.size()),
        n_words_(n_words),
        mh_steps_(mh_steps),
        stride_(1 + mh_steps) {

    n_tokens_ = 0;
    for (const auto& d : docs) n_tokens_ += d.size();

    z_.assign(n_tokens_ * stride_, 0);
    word_of_.resize(n_tokens_);
    csr_.resize(n_docs_ + 1);

    // ---- CSR, sorting each document's tokens by word ----
    std::vector<std::size_t> order;
    token_t at = 0;
    csr_[0] = 0;
    for (doc_t d = 0; d < n_docs_; d++) {
      const std::size_t len = docs[d].size();

      order.resize(len);
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(),
                [&](std::size_t a, std::size_t b) {
                  return docs[d][a] < docs[d][b];
                });

      for (std::size_t j = 0; j < len; j++) {
        const std::size_t src = order[j];
        word_of_[at] = static_cast<word_t>(docs[d][src]);
        // Initialize every pending proposal to old_z. An acceptance ratio is a
        // valid MH step only if the proposal came from the matching proposal
        // distribution, and at iteration 1 no pass has run to produce one --- so
        // any other value here gets resolved as a proposal it never was, partly
        // discarding the informed initialization (D8). Equal values make the
        // first resolve a no-op via the skip in the accept loops.
        const topic_t k = static_cast<topic_t>(zd[d][src]);
        for (std::size_t m = 0; m < stride_; m++) z_[at * stride_ + m] = k;
        at++;
      }
      csr_[d + 1] = at;
    }

    // ---- CSC, by counting sort on word ----
    csc_.assign(n_words_ + 1, 0);
    for (token_t i = 0; i < n_tokens_; i++) csc_[word_of_[i] + 1]++;
    for (word_t w = 0; w < n_words_; w++) csc_[w + 1] += csc_[w];

    csc_token_.resize(n_tokens_);
    std::vector<token_t> fill(csc_.begin(), csc_.end() - 1);
    for (token_t i = 0; i < n_tokens_; i++) csc_token_[fill[word_of_[i]]++] = i;
  }

  // ---- token accessors ----
  topic_t& old_z(token_t i)             { return z_[i * stride_]; }
  topic_t  old_z(token_t i) const       { return z_[i * stride_]; }
  topic_t& proposal(token_t i, std::size_t m)       { return z_[i * stride_ + 1 + m]; }
  topic_t  proposal(token_t i, std::size_t m) const { return z_[i * stride_ + 1 + m]; }

  // ---- CSR view (doc pass): tokens of document d are [begin, end) ----
  token_t doc_begin(doc_t d) const { return csr_[d]; }
  token_t doc_end(doc_t d)   const { return csr_[d + 1]; }

  // ---- CSC view (word pass): word_token(i) for i in [begin, end) ----
  token_t word_begin(word_t w) const { return csc_[w]; }
  token_t word_end(word_t w)   const { return csc_[w + 1]; }
  token_t word_token(token_t i) const { return csc_token_[i]; }

  word_t word_of(token_t i) const { return word_of_[i]; }

  std::size_t n_docs()   const { return n_docs_; }
  std::size_t n_words()  const { return n_words_; }
  std::size_t n_tokens() const { return n_tokens_; }
  std::size_t mh_steps() const { return mh_steps_; }

private:
  std::size_t n_docs_, n_words_, n_tokens_, mh_steps_, stride_;

  std::vector<topic_t> z_;          // n_tokens * (1 + mh_steps)
  std::vector<word_t>  word_of_;    // n_tokens
  std::vector<token_t> csr_;        // n_docs + 1
  std::vector<token_t> csc_;        // n_words + 1
  std::vector<token_t> csc_token_;  // n_tokens
};

}  // namespace warp
