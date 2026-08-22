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
// A NOTE ON eta_bar, AND A BUG IN THE REFERENCE. eta_bar is
// sum_v eta[k][v] -- a sum over the VOCABULARY, so V*eta under a scalar prior.
// text2vec's LDA.hpp:65 sets beta_bar = n_topic*beta, i.e. K*beta, which is
// dimensionally wrong: alpha_bar is a sum over topics but beta_bar is a sum over
// words. On this package's medium benchmark corpus (V = 4443, eta = 0.05) the
// correct value is 222.15 against the reference's 2.5 at K = 50. Since eta_bar
// only ever appears added to C_k, that difference substantially changes how far
// the acceptance ratio is shrunk toward 1. Porting it faithfully would give a
// perfectly valid MCMC chain converging to the wrong posterior. tidylda's own
// collapsed Gibbs sampler has this right (lda_gibbs2.cpp:333).

#include <RcppArmadillo.h>
#include <progress.hpp>
#include <vector>
#include <cmath>

#include "warp_rng.h"
#include "warp_alias.h"
#include "warp_corpus.h"

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
//' @param eta_in Scalar prior parameter for words over topics
//' @param iterations int number of sampling iterations
//' @param burnin int number of burn in iterations, -1 to disable averaging
//' @param calc_likelihood bool, calculate log likelihood?
//' @param likelihood_every int, evaluate the likelihood every n-th iteration
//' @param mh_steps int, Metropolis-Hastings proposals per token per pass
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
    const double&                                 eta_in,
    const std::size_t&                            iterations,
    const int&                                    burnin,
    const bool&                                   calc_likelihood,
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
  const double eta = eta_in;
  double alpha_bar = 0.0;
  for (std::size_t k = 0; k < K; k++) alpha_bar += alpha[k];
  // Sum over the VOCABULARY. See the header note on the reference's beta_bar.
  const double eta_bar = static_cast<double>(V) * eta;

  // D19: alias table over alpha for the doc-proposal's prior branch. The
  // reference draws that branch uniformly (LDA.hpp:191), which is proportional
  // to alpha_k only when alpha is symmetric; tidylda permits a vector, so a
  // uniform draw would sample from the wrong prior while running perfectly.
  // Built once, since D7 removed optimize_alpha and alpha is now fixed.
  AliasTable alpha_alias;
  alpha_alias.setup(alpha);

  // ---- counts ----
  // Cd is document-major, Cv is WORD-major internally so that the K counts for
  // one word are contiguous -- that contiguity is the word pass's whole point.
  // Output is transposed to the K x V the R side currently expects; D17 removes
  // that transpose in Phase 6.
  std::vector<int> Cd(D * K, 0);
  std::vector<int> Cv(V * K, 0);
  std::vector<long> Ck = Ck_in;

  const bool averaging = (burnin > -1);
  std::vector<long> Cd_sum(averaging ? D * K : 0, 0);
  std::vector<long> Cv_sum(averaging ? V * K : 0, 0);

  // ---- likelihood scratch ----
  std::vector<double> beta_prob(calc_likelihood ? V * K : 0);  // word-major
  std::vector<double> theta_prob(calc_likelihood ? K : 0);
  std::vector<double> log_likelihood;

  double lgeta = 0.0, lgalpha = 0.0;
  if (calc_likelihood) {
    for (std::size_t v = 0; v < V; v++) lgeta += lgamma(eta);
    lgeta = (lgeta - lgamma(eta_bar)) * static_cast<double>(K);
    for (std::size_t k = 0; k < K; k++) lgalpha += lgamma(alpha[k]);
    lgalpha = (lgalpha - lgamma(alpha_bar)) * static_cast<double>(D);
  }

  std::vector<double> prob(K);  // hoisted out of the per-word loop
  AliasTable word_alias;        // reused for every word, so no per-word alloc

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

      // 0. Rebuild C_doc[d] from old_z.
      for (std::size_t k = 0; k < K; k++) Cd_d[k] = 0;
      for (token_t i = begin; i < end; i++) Cd_d[corpus.old_z(i)]++;

      // 1. Resolve pending word-proposals via pi_w.
      for (token_t i = begin; i < end; i++) {
        for (std::size_t m = 0; m < mh_steps; m++) {
          const topic_t k  = corpus.old_z(i);
          const topic_t kp = corpus.proposal(i, m);
          if (kp == k) continue;

          const double accept =
              ((Cd_d[kp] + alpha[kp]) / (Cd_d[k] + alpha[k])) *
              ((Ck[k] + eta_bar) / (Ck[kp] + eta_bar));

          if (rng.unif() < accept) {
            // The reference updates C_all here but NOT C_doc (LDA.hpp:168-178),
            // while its word pass does maintain C_word (:116-117). It never
            // reads C_doc again, so the omission is harmless there. We export
            // C_doc and accumulate it below, so without these two lines the
            // matrix would stop matching old_z the moment anything is accepted
            // -- a silent bias in the posterior mean, exactly the kind D10
            // assumes away.
            Cd_d[kp]++;
            Cd_d[k]--;

            Ck[kp]++;
            Ck[k]--;

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

      int* Cv_w = &Cv[w * K];

      // 0. Rebuild C_word[w] from old_z.
      for (std::size_t k = 0; k < K; k++) Cv_w[k] = 0;
      for (token_t i = begin; i < end; i++) Cv_w[corpus.old_z(corpus.word_token(i))]++;

      // 1. Resolve pending doc-proposals via pi_d.
      for (token_t i = begin; i < end; i++) {
        const token_t tok = corpus.word_token(i);
        for (std::size_t m = 0; m < mh_steps; m++) {
          const topic_t k  = corpus.old_z(tok);
          const topic_t kp = corpus.proposal(tok, m);
          if (kp == k) continue;

          const double accept =
              ((Cv_w[kp] + eta) / (Cv_w[k] + eta)) *
              ((Ck[k] + eta_bar) / (Ck[kp] + eta_bar));

          if (rng.unif() < accept) {
            Cv_w[kp]++;
            Cv_w[k]--;

            Ck[kp]++;
            Ck[k]--;

            corpus.old_z(tok) = kp;
          }
        }
      }

      // D10: accumulate now, while C_word[w] is exactly the counts of old_z.
      if (accumulate) {
        long* s = &Cv_sum[w * K];
        for (std::size_t k = 0; k < K; k++) s[k] += Cv_w[k];
      }

      if (end == begin) continue;  // word absent from the corpus

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
      for (std::size_t k = 0; k < K; k++) prob[k] = Cv_w[k] + eta;
      word_alias.setup(prob);

      for (token_t i = begin; i < end; i++) {
        const token_t tok = corpus.word_token(i);
        for (std::size_t m = 0; m < mh_steps; m++) {
          corpus.proposal(tok, m) =
              static_cast<topic_t>(word_alias.sample(rng.unif(), rng.unif()));
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
    if (calc_likelihood && (t % likelihood_every == 0 || t == iterations - 1)) {

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

      double lp_eta = 0.0;
      for (std::size_t k = 0; k < K; k++) {
        const double denom = static_cast<double>(Ck[k]) + eta_bar;
        for (std::size_t v = 0; v < V; v++) {
          const double p = (static_cast<double>(Cv[v * K + k]) + eta) / denom;
          beta_prob[v * K + k] = p;
          lp_eta += (eta - 1.0) * std::log(p);
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

  IntegerMatrix Cv_out(K, V);
  for (std::size_t v = 0; v < V; v++)
    for (std::size_t k = 0; k < K; k++) Cv_out(k, v) = Cv[v * K + k];

  NumericMatrix Cd_mean(averaging ? D : 0, averaging ? K : 0);
  NumericMatrix Cv_mean(averaging ? K : 0, averaging ? V : 0);
  IntegerMatrix Cd_sum_out(averaging ? D : 0, averaging ? K : 0);
  IntegerMatrix Cv_sum_out(averaging ? K : 0, averaging ? V : 0);

  if (averaging) {
    const double n_post = static_cast<double>(iterations - burnin);
    for (std::size_t d = 0; d < D; d++) {
      for (std::size_t k = 0; k < K; k++) {
        Cd_mean(d, k)    = static_cast<double>(Cd_sum[d * K + k]) / n_post;
        Cd_sum_out(d, k) = static_cast<int>(Cd_sum[d * K + k]);
      }
    }
    for (std::size_t v = 0; v < V; v++) {
      for (std::size_t k = 0; k < K; k++) {
        Cv_mean(k, v)    = static_cast<double>(Cv_sum[v * K + k]) / n_post;
        Cv_sum_out(k, v) = static_cast<int>(Cv_sum[v * K + k]);
      }
    }
  }

  const std::size_t n_eval = log_likelihood.size() / 3;
  NumericMatrix ll_out(3, n_eval);
  for (std::size_t i = 0; i < n_eval; i++)
    for (std::size_t r = 0; r < 3; r++) ll_out(r, i) = log_likelihood[i * 3 + r];

  // eta goes back out as the K x V matrix the R side still expects. D20 removes
  // this materialization in Phase 4; until then format_eta() has already built
  // one anyway, so this changes nothing about peak memory.
  NumericMatrix eta_out(K, V);
  std::fill(eta_out.begin(), eta_out.end(), eta);

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
    _["eta"]            = eta_out
  );
}
