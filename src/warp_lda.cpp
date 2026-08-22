// warp_lda.cpp --- warpLDA sampler, single-threaded, scalar prior (Phase 2)
//
// Replaces fit_lda_c()'s O(N*K) collapsed Gibbs sweep with warpLDA's
// O(V*K + N) alternating Metropolis-Hastings scheme. See
// warp-planning/warplda-design-notes.md sections 2-4 for the derivation and
// warp-planning/warplda-roadmap.md section 6 for the phase plan.
//
// THE TWO PASSES, AND WHICH PROPOSAL EACH ONE HANDLES. Every token carries its
// current topic (old_z) plus one or more pending proposals. A pass DRAWS
// proposals of its own type and RESOLVES the ones the other pass left:
//
//   doc pass   working set C_doc[d]  | resolves WORD-proposals via pi_w
//                                    | draws    DOC-proposals  from q_d
//   word pass  working set C_word[w] | resolves DOC-proposals  via pi_d
//                                    | draws    WORD-proposals from q_w
//
// The cancellations that make this cheap:
//
//   pi_w = (C_doc[d][k'] + alpha[k']) / (C_doc[d][k] + alpha[k])
//          * (C_k + eta_bar) / (C_k' + eta_bar)
//
//        -- evaluated in the DOC pass, and NO v-indexed quantity survives. This
//           is what keeps the doc pass cache-resident, and it is a constraint to
//           preserve rather than a coincidence: it holds because q_w matches the
//           word factor of p exactly. Replacing q_w with a cheaper approximation
//           that ignores eta would resurrect eta_{kv} inside pi_w, where v varies
//           token to token, and defeat the whole design.
//
//   pi_d = (C_word[w][k'] + eta) / (C_word[w][k] + eta)
//          * (C_k + eta_bar) / (C_k' + eta_bar)
//
//        -- evaluated in the WORD pass, where C_word[w] and (from Phase 3) the
//           eta column for w are both already in cache.
//
// THE TWO PRIOR SUMS ARE OVER DIFFERENT INDICES. alpha_bar is a sum over topics,
// sum_k alpha_k. eta_bar is a K-vector summed over the VOCABULARY,
// eta_bar[k] = sum_v eta_{kv}. They are easy to conflate and the second is
// load-bearing: eta_bar appears only as an addition to C_k, so getting it wrong
// leaves a perfectly well-behaved MCMC chain that converges to a different
// posterior.
//
// FROZEN TOPICS (prediction). With beta held fixed the target becomes
// p(k) ~ beta_hat[k][v] * (C^d_dk + alpha_k), and both acceptance ratios get
// simpler rather than harder:
//
//   pi_w = (C^d_dk' + alpha_k') / (C^d_dk + alpha_k)      -- doc pass
//   pi_d = beta_hat[k'][v] / beta_hat[k][v]               -- word pass
//
// No C_k, no eta_bar, no C^v: those terms exist in training only because beta is
// being re-estimated. And because q_w ~ beta_hat[.][v] is fixed for the whole
// run, its alias tables are built once instead of per word per iteration, which
// is where the prediction path's speed actually comes from. All of this is
// selected by a runtime flag rather than a separate kernel -- measured at R's
// -O2 the branch is free (the loop is divides and scattered loads, and a
// loop-invariant branch predicts perfectly against that), and one kernel keeps
// the MH bookkeeping in a single copy that cannot drift out of sync.

#include <RcppArmadillo.h>
#include <progress.hpp>
#include <vector>
#include <cmath>

#include "warp_rng.h"
#include "warp_alias.h"
#include "warp_corpus.h"
#include "warp_eta.h"

using namespace Rcpp;
using namespace warp;

namespace {

// Draw a 64-bit master seed from R's stream, on the main thread, so that
// set.seed() governs the entire run (D12). Every per-work-item stream is
// derived from this value and never touches R's RNG again.
uint64_t draw_master_seed() {
  const uint64_t hi = static_cast<uint64_t>(R::unif_rand() * 4294967296.0);
  const uint64_t lo = static_cast<uint64_t>(R::unif_rand() * 4294967296.0);
  return (hi << 32) ^ lo;
}

}  // namespace


//' Fit an LDA model with the warpLDA sampler
//' @keywords internal
//' @description Metropolis-Hastings replacement for \code{fit_lda_c}. Phase 2:
//'   single-threaded, scalar \code{eta} only.
//' @param Docs List of vectors of word indices, one per document
//' @param Zd_in List of initial topic assignments matching \code{Docs}
//' @param Cd_in IntegerMatrix, documents by topics
//' @param Cv_in IntegerMatrix, topics by words
//' @param Ck_in Vector of token counts per topic
//' @param alpha_in Vector of prior parameters for topics over documents
//' @param eta_in NumericMatrix, topics by words. Prior for words over topics
//' @param iterations int number of sampling iterations
//' @param burnin int number of burn in iterations, -1 to disable averaging
//' @param calc_likelihood bool, calculate log likelihood?
//' @param likelihood_every int, evaluate the likelihood every n-th iteration
//' @param mh_steps int, Metropolis-Hastings proposals per token per pass
//' @param freeze_topics bool, hold topics fixed for prediction?
//' @param Beta_in NumericMatrix, topics by words. The fitted beta, used only
//'   when \code{freeze_topics = TRUE}
//' @param verbose bool, show a progress bar?
//' @return Returns a list with the same names as \code{fit_lda_c}.
// [[Rcpp::export]]
Rcpp::List fit_lda_warp(
    const std::vector<std::vector<std::size_t>>&  Docs,
    const std::vector<std::vector<std::size_t>>&  Zd_in,
    const IntegerMatrix&                          Cd_in,
    const IntegerMatrix&                          Cv_in,
    const std::vector<long>&                      Ck_in,
    const std::vector<double>&                    alpha_in,
    const NumericMatrix&                          eta_in,
    const std::size_t&                            iterations,
    const int&                                    burnin,
    const bool&                                   calc_likelihood,
    const NumericMatrix&                          Beta_in,
    const bool&                                   freeze_topics = false,
    const std::size_t&                            likelihood_every = 10,
    const std::size_t&                            mh_steps = 1,
    const bool&                                   verbose = true
) {

  const std::size_t K = Ck_in.size();
  const std::size_t D = Docs.size();
  const std::size_t V = Cv_in.ncol();

  if (K > 65535) stop("warpLDA engine supports at most 65535 topics");
  if (mh_steps < 1) stop("mh_steps must be at least 1");
  if (burnin >= 0 && static_cast<std::size_t>(burnin) >= iterations) {
    stop("burnin must be less than iterations");
  }

  Corpus corpus(Docs, Zd_in, V, mh_steps);
  const std::size_t N = corpus.n_tokens();

  // ---- priors ----
  const std::vector<double> alpha = alpha_in;
  double alpha_bar = 0.0;
  for (std::size_t k = 0; k < K; k++) alpha_bar += alpha[k];

  // eta and its per-topic vocabulary sums. Unused when topics are frozen, so
  // skip building them: at K=500, V=1e5 that is 200MB not allocated.
  const Eta eta = freeze_topics ? Eta(nullptr, 0, 0) : Eta(eta_in.begin(), K, V);

  // With topics frozen the word-proposal is q_w(k) ~ beta_hat[k][v], fixed for
  // the whole run. Same column-major reasoning as eta: R's K x V matrix already
  // stores each word's K values contiguously.
  const double* beta_hat = freeze_topics ? Beta_in.begin() : nullptr;

  // D19: alias table over alpha for the doc-proposal's prior branch. That
  // branch must draw proportional to alpha_k, which a uniform draw achieves only
  // when alpha is symmetric -- and tidylda permits a vector. Built once, since
  // D7 removed optimize_alpha and alpha is fixed for the run.
  AliasTable alpha_alias;
  alpha_alias.setup(alpha);

  // ---- counts ----
  // Cd is document-major, Cv is WORD-major internally so that the K counts for
  // one word are contiguous -- that contiguity is the word pass's whole point.
  // Output is transposed to the K x V the R side currently expects; D17 removes
  // that transpose in Phase 6.
  std::vector<int> Cd(D * K, 0);
  std::vector<int> Cv(freeze_topics ? 0 : V * K, 0);
  std::vector<long> Ck = Ck_in;

  const bool averaging = (burnin > -1);
  std::vector<long> Cd_sum(averaging ? D * K : 0, 0);
  std::vector<long> Cv_sum((averaging && !freeze_topics) ? V * K : 0, 0);

  // ---- likelihood scratch ----
  std::vector<double> beta_prob((calc_likelihood && !freeze_topics) ? V * K : 0);  // word-major
  std::vector<double> theta_prob(calc_likelihood ? K : 0);
  std::vector<double> beta_denom(calc_likelihood ? K : 0);
  std::vector<double> log_likelihood;

  double lgeta = 0.0, lgalpha = 0.0;
  if (calc_likelihood && !freeze_topics) {
    // Dirichlet log-normalizer, summed per topic: each topic now has its own
    // prior over words, so this no longer factors into one term times K.
    for (std::size_t k = 0; k < K; k++) {
      for (std::size_t v = 0; v < V; v++) lgeta += lgamma(eta.at(v, k));
      lgeta -= lgamma(eta.bar(k));
    }
    for (std::size_t k = 0; k < K; k++) lgalpha += lgamma(alpha[k]);
    lgalpha = (lgalpha - lgamma(alpha_bar)) * static_cast<double>(D);
  }

  std::vector<double> prob(K);  // hoisted out of the per-word loop
  AliasTable word_alias;        // reused for every word, so no per-word alloc

  // D9: with topics frozen, q_w(k) ~ beta_hat[k][v] does not change from one
  // iteration to the next, so its alias tables are built ONCE here rather than
  // rebuilt per word per iteration -- the dominant cost of the word pass, and
  // the real saving on the prediction path.
  //
  // Built only for words that actually occur. The prediction DTM is aligned to
  // the model's full vocabulary, but a word with no tokens is never sampled, so
  // the table count is bounded by distinct words present rather than by V.
  std::vector<AliasTable> frozen_alias;
  if (freeze_topics) {
    frozen_alias.resize(V);
    for (std::size_t w = 0; w < V; w++) {
      if (corpus.word_begin(w) == corpus.word_end(w)) continue;
      const double* bw = beta_hat + w * K;
      for (std::size_t k = 0; k < K; k++) prob[k] = bw[k];
      frozen_alias[w].setup(prob);
    }
  }

  const uint64_t master = draw_master_seed();
  Progress progress(iterations, verbose);

  // =========================================================================
  // Main loop
  // =========================================================================
  for (std::size_t t = 0; t < iterations; t++) {

    const bool accumulate = averaging && static_cast<int>(t) >= burnin;

    // -----------------------------------------------------------------------
    // DOC PASS: resolve word-proposals, draw doc-proposals
    // -----------------------------------------------------------------------
    for (std::size_t d = 0; d < D; d++) {
      Xoshiro256pp rng = work_item_rng(master, t, Pass::doc, d);

      const token_t begin = corpus.doc_begin(d);
      const token_t end   = corpus.doc_end(d);
      const std::size_t Ld = end - begin;
      if (Ld == 0) continue;

      int* Cd_d = &Cd[d * K];
      const double* eta_bar = freeze_topics ? nullptr : eta.bar().data();

      // 0. Rebuild C_doc[d] from old_z.
      for (std::size_t k = 0; k < K; k++) Cd_d[k] = 0;
      for (token_t i = begin; i < end; i++) Cd_d[corpus.old_z(i)]++;

      // 1. Resolve pending word-proposals via pi_w.
      for (token_t i = begin; i < end; i++) {
        for (std::size_t m = 0; m < mh_steps; m++) {
          const topic_t k  = corpus.old_z(i);
          const topic_t kp = corpus.proposal(i, m);
          if (kp == k) continue;

          // pi_w. The (C_k + eta_bar) factor exists only because beta is
          // being re-estimated; with topics frozen the word factor of p is the
          // fixed beta_hat, which cancels against q_w exactly, leaving just the
          // document ratio.
          const double accept =
              freeze_topics
                  ? ((Cd_d[kp] + alpha[kp]) / (Cd_d[k] + alpha[k]))
                  : ((Cd_d[kp] + alpha[kp]) / (Cd_d[k] + alpha[k])) *
                    ((Ck[k] + eta_bar[k]) / (Ck[kp] + eta_bar[kp]));

          if (rng.unif() < accept) {
            // C_doc must be maintained here, not just rebuilt in step 0: it is
            // accumulated into Cd_sum below and exported, so it has to keep
            // matching old_z through the whole accept loop. Dropping these two
            // lines biases the posterior mean without changing anything the
            // sampler itself reads.
            Cd_d[kp]++;
            Cd_d[k]--;

            if (!freeze_topics) {
              Ck[kp]++;
              Ck[k]--;
            }

            corpus.old_z(i) = kp;
          }
        }
      }

      // D10: accumulate now, while C_doc[d] is exactly the counts of old_z.
      if (accumulate) {
        long* s = &Cd_sum[d * K];
        for (std::size_t k = 0; k < K; k++) s[k] += Cd_d[k];
      }

      // 2. Draw doc-proposals from q_d proportional to C_doc[d][k] + alpha[k].
      // Drawn in O(1) without a table: with probability Ld/(Ld + alpha_bar)
      // copy a uniformly random token's topic (which is distributed exactly as
      // C_doc[d][k]/Ld), otherwise draw from the alpha alias table.
      const double p_position =
          static_cast<double>(Ld) / (static_cast<double>(Ld) + alpha_bar);

      for (token_t i = begin; i < end; i++) {
        for (std::size_t m = 0; m < mh_steps; m++) {
          if (rng.unif() < p_position) {
            corpus.proposal(i, m) = corpus.old_z(begin + rng.below(Ld));
          } else {
            corpus.proposal(i, m) =
                static_cast<topic_t>(alpha_alias.sample(rng.unif(), rng.unif()));
          }
        }
      }
    }

    // -----------------------------------------------------------------------
    // WORD PASS: resolve doc-proposals, draw word-proposals
    // -----------------------------------------------------------------------
    for (std::size_t w = 0; w < V; w++) {
      Xoshiro256pp rng = work_item_rng(master, t, Pass::word, w);

      const token_t begin = corpus.word_begin(w);
      const token_t end   = corpus.word_end(w);

      if (begin == end) continue;  // word absent from this corpus

      int* Cv_w = freeze_topics ? nullptr : &Cv[w * K];
      const double* beta_w = freeze_topics ? beta_hat + w * K : nullptr;

      // Hoist the row-invariant pointers. w is fixed for this whole block, so
      // recomputing w*K inside the token loop -- which runs once per token per
      // MH step -- is pure overhead against a matrix prior.
      //
      // eta_w points at FLOATS (D5 stores eta single-precision to halve the
      // largest allocation). Every read below is cast to double before entering
      // any arithmetic: `int + float` would otherwise evaluate the entire
      // acceptance ratio in single precision, which is exactly the drift D5
      // requires be kept out of the MH accept decision.
      const float*  eta_w   = freeze_topics ? nullptr : eta.column(w);
      const double* eta_bar = freeze_topics ? nullptr : eta.bar().data();

      // 0. Rebuild C_word[w] from old_z. Frozen topics never touch C^v.
      if (!freeze_topics) {
        for (std::size_t k = 0; k < K; k++) Cv_w[k] = 0;
        for (token_t i = begin; i < end; i++) Cv_w[corpus.old_z(corpus.word_token(i))]++;
      }

      // 1. Resolve pending doc-proposals via pi_d.
      for (token_t i = begin; i < end; i++) {
        const token_t tok = corpus.word_token(i);
        for (std::size_t m = 0; m < mh_steps; m++) {
          const topic_t k  = corpus.old_z(tok);
          const topic_t kp = corpus.proposal(tok, m);
          if (kp == k) continue;

          // pi_d. Frozen topics reduce to the ratio of fixed beta_hat entries:
          // the document factor of p is exactly q_d and cancels, and there is no
          // C_k normalizer because beta_hat is already a distribution.
          const double e_kp = freeze_topics ? 0.0 : static_cast<double>(eta_w[kp]);
          const double e_k  = freeze_topics ? 0.0 : static_cast<double>(eta_w[k]);

          const double accept =
              freeze_topics
                  ? (beta_w[kp] / beta_w[k])
                  : ((Cv_w[kp] + e_kp) / (Cv_w[k] + e_k)) *
                    ((Ck[k] + eta_bar[k]) / (Ck[kp] + eta_bar[kp]));

          if (rng.unif() < accept) {
            if (!freeze_topics) {
              Cv_w[kp]++;
              Cv_w[k]--;

              Ck[kp]++;
              Ck[k]--;
            }

            corpus.old_z(tok) = kp;
          }
        }
      }

      // D10: accumulate now, while C_word[w] is exactly the counts of old_z.
      if (accumulate && !freeze_topics) {
        long* s = &Cv_sum[w * K];
        for (std::size_t k = 0; k < K; k++) s[k] += Cv_w[k];
      }

      // 2. Draw word-proposals from q_w proportional to C_word[w][k] + eta.
      //
      // D15: this is the O(V*K) formulation -- a dense pass over all K topics
      // for every word type, every iteration. The O(N) alternative from the
      // warpLDA paper splits q_w into a sparse count part (C_word has at most
      // min(n_w, K) nonzeros) plus a dense prior part served by a shared table.
      // Under tLDA that dense part is eta's column for w, which is word-
      // specific, so recovering O(1) means precomputing V alias tables over
      // eta's columns. They would be built once, since eta is fixed, but cost
      // roughly 2x the size of eta in permanent memory -- about 400MB at
      // V = 1e5, K = 500, on top of eta itself. Accepted as-is for now; revisit
      // if profiling at high K shows this dominating.
      if (!freeze_topics) {
        for (std::size_t k = 0; k < K; k++)
          prob[k] = Cv_w[k] + static_cast<double>(eta_w[k]);
        word_alias.setup(prob);
      }
      const AliasTable& tbl = freeze_topics ? frozen_alias[w] : word_alias;

      for (token_t i = begin; i < end; i++) {
        const token_t tok = corpus.word_token(i);
        for (std::size_t m = 0; m < mh_steps; m++) {
          corpus.proposal(tok, m) =
              static_cast<topic_t>(tbl.sample(rng.unif(), rng.unif()));
        }
      }
    }

    // -----------------------------------------------------------------------
    // Log likelihood (D11)
    // -----------------------------------------------------------------------
    // Evaluated every likelihood_every-th iteration. THIS IS NOT THINNING: the
    // chain advances every iteration and every post-burnin iteration still
    // contributes to Cd_sum and Cv_sum above. Only this read-only diagnostic,
    // which feeds nothing the sampler uses, runs less often.
    if (calc_likelihood && !freeze_topics &&
        (t % likelihood_every == 0 || t == iterations - 1)) {

      // The count averages are marginal and need no synchronized joint sample,
      // but the likelihood does: theta comes from C_doc and beta from C_word at
      // the SAME point in the chain. C_word is current -- the word pass just
      // rebuilt and maintained every column. C_doc is stale, because the word
      // pass moved tokens after the doc pass last touched it. One O(N) sweep
      // with sequential access fixes that; it is roughly a sixth of an
      // iteration and is not the expensive part of this block.
      std::fill(Cd.begin(), Cd.end(), 0);
      for (std::size_t d = 0; d < D; d++) {
        for (token_t i = corpus.doc_begin(d); i < corpus.doc_end(d); i++) {
          Cd[d * K + corpus.old_z(i)]++;
        }
      }

      // Word-major outer loop. eta, Cv and beta_prob all store a word's K values
      // contiguously, so iterating topics inside words walks all three
      // sequentially; the transpose of this nest strides every one of them by K.
      for (std::size_t k = 0; k < K; k++) {
        beta_denom[k] = static_cast<double>(Ck[k]) + eta.bar(k);
      }

      double lp_eta = 0.0;
      for (std::size_t v = 0; v < V; v++) {
        const float* ev = eta.column(v);
        const int*   cv = &Cv[v * K];
        double*      bp = &beta_prob[v * K];
        for (std::size_t k = 0; k < K; k++) {
          const double e = static_cast<double>(ev[k]);
          const double p = (static_cast<double>(cv[k]) + e) / beta_denom[k];
          bp[k] = p;
          lp_eta += (e - 1.0) * std::log(p);
        }
      }
      lp_eta += lgeta;

      double lp_alpha = 0.0;
      double lpd = 0.0;

      for (std::size_t d = 0; d < D; d++) {
        const int* Cd_d = &Cd[d * K];

        double denom = 0.0;
        for (std::size_t k = 0; k < K; k++) denom += Cd_d[k] + alpha[k];
        for (std::size_t k = 0; k < K; k++) {
          theta_prob[k] = (Cd_d[k] + alpha[k]) / denom;
          lp_alpha += (alpha[k] - 1.0) * std::log(theta_prob[k]);
        }

        // Design notes section 7, mitigation 1: take the inner sum once per
        // distinct (d, v) pair rather than once per token occurrence, giving
        // O(nnz * K) instead of O(N * K). Tokens are word-sorted within a
        // document, so the distinct pairs are just the runs.
        token_t i = corpus.doc_begin(d);
        const token_t end = corpus.doc_end(d);
        while (i < end) {
          const word_t v = corpus.word_of(i);
          token_t j = i;
          while (j < end && corpus.word_of(j) == v) j++;
          const double count = static_cast<double>(j - i);

          const double* bp = &beta_prob[v * K];
          double p = 0.0;
          for (std::size_t k = 0; k < K; k++) p += theta_prob[k] * bp[k];

          lpd += count * std::log(p);
          i = j;
        }
      }
      lp_alpha += lgalpha;

      log_likelihood.push_back(static_cast<double>(t));
      log_likelihood.push_back(lpd);
      log_likelihood.push_back(lpd + lp_alpha + lp_eta);
    }

    if (verbose) progress.increment();
    Rcpp::checkUserInterrupt();
  }

  // =========================================================================
  // Assemble the return value, matching fit_lda_c's contract exactly
  // =========================================================================

  // C_word is current -- the word pass rebuilt and maintained every column and
  // nothing has moved since. C_doc is not: the word pass reassigned tokens after
  // the doc pass last touched it. The likelihood block above happens to leave a
  // fresh C_doc behind, but only when calc_likelihood is TRUE, so rebuild it
  // here unconditionally. new_tidylda() reads Cd directly whenever burnin == -1.
  std::fill(Cd.begin(), Cd.end(), 0);
  for (std::size_t d = 0; d < D; d++) {
    for (token_t i = corpus.doc_begin(d); i < corpus.doc_end(d); i++) {
      Cd[d * K + corpus.old_z(i)]++;
    }
  }

  // Cd is already D x K. Cv is word-major internally and must go out as K x V.
  IntegerMatrix Cd_out(D, K);
  for (std::size_t d = 0; d < D; d++)
    for (std::size_t k = 0; k < K; k++) Cd_out(d, k) = Cd[d * K + k];

  // Frozen topics never build C^v, so it goes out empty rather than as a
  // zero matrix that could be mistaken for a real count.
  // new_tidylda(is_prediction = TRUE) reads only Cd/Cd_mean and alpha.
  IntegerMatrix Cv_out(freeze_topics ? 0 : K, freeze_topics ? 0 : V);
  if (!freeze_topics)
    for (std::size_t v = 0; v < V; v++)
      for (std::size_t k = 0; k < K; k++) Cv_out(k, v) = Cv[v * K + k];

  NumericMatrix Cd_mean(averaging ? D : 0, averaging ? K : 0);
  const bool cv_out = averaging && !freeze_topics;
  NumericMatrix Cv_mean(cv_out ? K : 0, cv_out ? V : 0);
  IntegerMatrix Cd_sum_out(averaging ? D : 0, averaging ? K : 0);
  IntegerMatrix Cv_sum_out(cv_out ? K : 0, cv_out ? V : 0);

  if (averaging) {
    const double n_post = static_cast<double>(iterations - burnin);
    for (std::size_t d = 0; d < D; d++) {
      for (std::size_t k = 0; k < K; k++) {
        Cd_mean(d, k)    = static_cast<double>(Cd_sum[d * K + k]) / n_post;
        Cd_sum_out(d, k) = static_cast<int>(Cd_sum[d * K + k]);
      }
    }
    if (cv_out) {
      for (std::size_t v = 0; v < V; v++) {
        for (std::size_t k = 0; k < K; k++) {
          Cv_mean(k, v)    = static_cast<double>(Cv_sum[v * K + k]) / n_post;
          Cv_sum_out(k, v) = static_cast<int>(Cv_sum[v * K + k]);
        }
      }
    }
  }

  const std::size_t n_eval = log_likelihood.size() / 3;
  NumericMatrix ll_out(3, n_eval);
  for (std::size_t i = 0; i < n_eval; i++)
    for (std::size_t r = 0; r < 3; r++) ll_out(r, i) = log_likelihood[i * 3 + r];

  // eta is not modified by the sampler, so hand back the input rather than
  // spending a K x V copy on it.

  return Rcpp::List::create(
    _["Cd"]             = Cd_out,
    _["Cv"]             = Cv_out,
    _["Ck"]             = Ck,
    _["Cd_mean"]        = Cd_mean,
    _["Cv_mean"]        = Cv_mean,
    _["Cd_sum"]         = Cd_sum_out,
    _["Cv_sum"]         = Cv_sum_out,
    _["log_likelihood"] = ll_out,
    _["alpha"]          = alpha,
    _["eta"]            = eta_in
  );
}
