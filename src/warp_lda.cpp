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
#include <RcppThread.h>
#include <memory>
#include <vector>
#include <cmath>

#include "warp_rng.h"
#include "warp_alias.h"
#include "warp_corpus.h"
#include "warp_eta.h"
#include "warp_init_sample.h"
#include "matrix_conversions.h"

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
//' @param dtm_in arma::sp_mat document term matrix, documents by words
//' @param Cd_start IntegerMatrix, documents by topics. Initial document-topic
//'   counts, \code{theta_initial * rowSums(dtm)} from the R side
//' @param alpha_in Vector of prior parameters for topics over documents
//' @param eta_in NumericMatrix, topics by words. Prior for words over topics
//' @param iterations int number of sampling iterations. Zero initializes and
//'   returns without sampling, which is how the initialization is inspected
//' @param burnin int number of burn in iterations, -1 to disable averaging
//' @param calc_likelihood bool, calculate log likelihood?
//' @param likelihood_every int, evaluate the likelihood every n-th iteration
//' @param mh_steps int, Metropolis-Hastings proposals per token per pass
//' @param freeze_topics bool, hold topics fixed for prediction?
//' @param Beta_in NumericMatrix, topics by words. The fitted beta, used only
//'   when \code{freeze_topics = TRUE}
//' @param threads int, number of worker threads. Results are identical at any
//'   thread count (D12), so this trades wall clock for cores and nothing else
//' @param verbose bool, show a progress bar?
//' @return Returns a list of counts and diagnostics. \code{Cd} and
//'   \code{Cd_mean} are documents by topics; \code{Cv} and \code{Cv_mean} are
//'   words by topics (D17). Only the pair the caller can use is materialized:
//'   \code{Cd}/\code{Cv} when \code{burnin} is -1, \code{Cd_mean}/
//'   \code{Cv_mean} otherwise. The other pair comes back 0 x 0.
// [[Rcpp::export]]
Rcpp::List fit_lda_warp(
    const arma::sp_mat&                           dtm_in,
    const IntegerMatrix&                          Cd_start,
    const std::vector<double>&                    alpha_in,
    const NumericMatrix&                          eta_in,
    const std::size_t&                            iterations,
    const int&                                    burnin,
    const bool&                                   calc_likelihood,
    const NumericMatrix&                          Beta_in,
    const bool&                                   freeze_topics = false,
    const std::size_t&                            likelihood_every = 10,
    const std::size_t&                            mh_steps = 1,
    const std::size_t&                            threads = 1,
    const bool&                                   verbose = true
) {

  const std::size_t K = alpha_in.size();
  const std::size_t D = dtm_in.n_rows;
  const std::size_t V = dtm_in.n_cols;

  if (K > 65535) stop("warpLDA engine supports at most 65535 topics");
  if (mh_steps < 1) stop("mh_steps must be at least 1");
  if (burnin >= 0 && static_cast<std::size_t>(burnin) >= iterations) {
    stop("burnin must be less than iterations");
  }

  // =========================================================================
  // Parallel scaffolding (Phases 5 and 5.5)
  // =========================================================================
  // Both the initialization and the two sampling passes partition cleanly on
  // their own index. Chunk count is sized for load balance only: more chunks
  // than threads lets the pool even out documents and words of very different
  // lengths, which matters because word frequencies are Zipfian.
  const std::size_t n_chunks =
      std::max<std::size_t>(1, std::min<std::size_t>(std::max(D, V), threads * 4));

  auto chunk_lo = [&](std::size_t total, std::size_t c) { return (total * c) / n_chunks; };
  auto chunk_hi = [&](std::size_t total, std::size_t c) { return (total * (c + 1)) / n_chunks; };

  std::unique_ptr<RcppThread::ThreadPool> pool;
  if (threads > 1) pool.reset(new RcppThread::ThreadPool(threads));

  auto run_parallel = [&](auto&& body) {
    if (pool) {
      pool->parallelFor(0, static_cast<int>(n_chunks), body);
      pool->wait();
    } else {
      for (std::size_t c = 0; c < n_chunks; c++) body(static_cast<int>(c));
    }
  };

  // Drawn on the main thread so set.seed() governs the whole run, and drawn
  // here because initialization now needs it too (Phase 5.5).
  const uint64_t master = draw_master_seed();

  // =========================================================================
  // Initialization (D16): build the token structure and sample each token's
  // starting topic, in one walk of the sparse DTM.
  //
  // This was create_lexicon(), which returned Docs and Zd to R as
  // vector<vector<size_t>> -- 16 bytes per token out and straight back in.
  // Nothing of the kind is materialized now.
  //
  // The initialization is INFORMED, not uniform (D8): each token is drawn from
  // P(z) proportional to beta_hat[k][v] * theta_hat[d][k], in log space. That is
  // what makes seeded and transfer-learned models work, and it is why warpLDA's
  // own uniform init is discarded wholesale.
  //
  // The draw order is load-bearing. qz is computed once per distinct (d,v) pair
  // and reused for every repeat of that word, C_doc is NOT updated within a
  // document, and tokens are visited documents-ascending then words-ascending.
  // Changing any of that changes which random numbers each token consumes.
  //
  // One uniform per token, drawn from the work item's own stream (D12) rather
  // than R's. That is what makes this parallel: R's RNG is main-thread-only
  // (section 5 invariant), so as long as initialization drew from it, it could
  // not be threaded -- and being O(N*K) and serial it capped total speedup near
  // 2x at K=200 however well the sampler scaled.
  //
  // Pass::init keeps this off the doc pass's iteration-0 stream. Sharing it
  // would drive a token's starting topic and its first proposal from the same
  // uniform.
  // =========================================================================
  const arma::sp_mat dtm = dtm_in.t();   // V x D: a column is now a document

  // Raw CSC arrays rather than iterators. arma's sparse iterators can trigger
  // lazy synchronization, which is not something to invoke from several threads
  // at once; after an explicit sync these are plain reads.
  dtm.sync();
  const arma::uword* col_ptr = dtm.col_ptrs;
  const arma::uword* row_idx = dtm.row_indices;
  const double*      val     = dtm.values;

  auto Cd0  = mat_to_vec(Cd_start, true);   // Cd0[d][k]
  auto Beta = mat_to_vec(Beta_in, true);    // Beta[k][v]

  double sum_alpha = 0.0;
  for (std::size_t k = 0; k < K; k++) sum_alpha += alpha_in[k];

  // Document lengths, then the CSR layout. O(nnz) and serial, but memory-bound
  // and roughly K times cheaper than the sampling it precedes. finalize()'s
  // counting sort below is the other serial piece; both are the next candidates
  // if initialization ever shows up in a profile again.
  std::vector<token_t> doc_len(D, 0);
  std::size_t N = 0;
  for (std::size_t d = 0; d < D; d++) {
    double nd = 0.0;
    for (arma::uword p = col_ptr[d]; p < col_ptr[d + 1]; p++) nd += val[p];
    doc_len[d] = static_cast<token_t>(nd);
    N += doc_len[d];
  }

  Corpus corpus(D, V, N, mh_steps);
  corpus.begin_build(doc_len);

  {
    // Per-chunk buffers, so nothing is allocated per token or shared.
    std::vector<std::vector<double>> qz_buf(n_chunks, std::vector<double>(K));
    std::vector<std::vector<double>> sc_buf(n_chunks, std::vector<double>(K));

    run_parallel([&](int c) {
      std::vector<double>& qz = qz_buf[c];
      std::vector<double>& scratch = sc_buf[c];

      for (std::size_t d = chunk_lo(D, c); d < chunk_hi(D, c); d++) {
        // Each document owns a disjoint slot range, so the fill needs no
        // coordination at all.
        Xoshiro256pp rng = work_item_rng(master, 0, Pass::init, d);
        token_t at = corpus.doc_begin(d);

        const double denom_term =
            std::log(static_cast<double>(doc_len[d]) + sum_alpha - 1.0);

        for (arma::uword pp = col_ptr[d]; pp < col_ptr[d + 1]; pp++) {
          const std::size_t v = row_idx[pp];
          const std::size_t count = static_cast<std::size_t>(val[pp]);
          if (count == 0) continue;

          for (std::size_t k = 0; k < K; k++) {
            qz[k] = std::log(Beta[k][v]) +
                    std::log(Cd0[d][k] + alpha_in[k]) - denom_term;
          }

          for (std::size_t i = 0; i < count; i++) {
            corpus.set_token(at++, static_cast<word_t>(v),
                             static_cast<topic_t>(
                                 sample_log_weights(qz, rng.unif(), scratch)));
          }
        }
      }
    });
    corpus.finalize();
  }

  // ---- priors ----
  const std::vector<double> alpha = alpha_in;
  double alpha_bar = 0.0;
  for (std::size_t k = 0; k < K; k++) alpha_bar += alpha[k];

  // eta and its per-topic vocabulary sums. Unused when topics are frozen, so
  // skip building them: at K=500, V=1e5 that is 200MB not allocated.
  // D20: a 1 x 1 eta_in means the user passed a scalar, and format_eta() did not
  // materialize K x V for it. Scalar mode holds one K-length array instead.
  const bool eta_is_scalar = (eta_in.nrow() == 1 && eta_in.ncol() == 1);
  const Eta eta =
      freeze_topics   ? Eta(nullptr, 0, 0)
      : eta_is_scalar ? Eta(static_cast<double>(eta_in(0, 0)), K, V)
                      : Eta(eta_in.begin(), K, V);

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

  // Derive all three count views from the topics just sampled.
  //
  // Ck is the one that MUST be built here: it is never rebuilt during a pass,
  // only maintained incrementally on every acceptance, so it has to start
  // exactly right. Cd and Cv are rebuilt by their own passes and would be
  // filled anyway -- deriving them now costs one O(N) sweep and makes the state
  // well defined at iterations = 0, which is what lets the initialization be
  // inspected and tested on its own.
  std::vector<long> Ck(K, 0);
  for (token_t i = 0; i < N; i++) Ck[corpus.old_z(i)]++;

  for (std::size_t d = 0; d < D; d++)
    for (token_t i = corpus.doc_begin(d); i < corpus.doc_end(d); i++)
      Cd[d * K + corpus.old_z(i)]++;

  if (!freeze_topics)
    for (token_t i = 0; i < N; i++)
      Cv[corpus.word_of(i) * K + corpus.old_z(i)]++;


  const bool averaging = (burnin > -1);
  std::vector<long> Cd_sum(averaging ? D * K : 0, 0);
  std::vector<long> Cv_sum((averaging && !freeze_topics) ? V * K : 0, 0);

  // C_k is READ-ONLY within a sampling pass. Each chunk accumulates its own
  // delta and the deltas are summed in afterwards. Three consequences, all
  // wanted:
  //
  //   * every work item sees the same C_k however chunks land on threads;
  //   * integer addition is associative and exact, so merge order is irrelevant;
  //   * a work item's delta contribution depends only on its own state and the
  //     snapshot, so even the CHUNK COUNT cannot change the result.
  //
  // Results therefore do not depend on `threads`, which is what D12 asks for.
  //
  // ATOMICS WOULD NOT DO. They remove the race but leave what each work item
  // READS dependent on interleaving, so results would drift with thread count --
  // exactly what D12 forbids. The requirement is determinism, not just safety.
  struct Scratch {
    std::vector<double> prob;      // q_w construction buffer
    AliasTable          alias;     // per-word proposal table
    std::vector<long>   ck_delta;  // this chunk's contribution to C_k
    void reset(std::size_t k) { prob.assign(k, 0.0); ck_delta.assign(k, 0); }
    void clear_delta() { std::fill(ck_delta.begin(), ck_delta.end(), 0); }
  };

  std::vector<Scratch> scratch(n_chunks);
  for (auto& sc : scratch) sc.reset(K);

  auto run_chunks = [&](auto&& body) {
    for (auto& sc : scratch) sc.clear_delta();
    run_parallel(body);
    for (const auto& sc : scratch)
      for (std::size_t k = 0; k < K; k++) Ck[k] += sc.ck_delta[k];
  };

  // ---- likelihood scratch ----
  std::vector<double> beta_prob((calc_likelihood && !freeze_topics) ? V * K : 0);  // word-major
  std::vector<double> theta_prob(calc_likelihood ? K : 0);
  std::vector<double> beta_denom(calc_likelihood ? K : 0);
  std::vector<double> log_likelihood;

  // lgamma of each prior entry, for the collapsed joint below. Precomputed
  // because the joint subtracts one of these per nonzero count, and under a
  // scalar prior (D20) there is exactly one distinct value to compute.
  std::vector<double> lg_alpha(calc_likelihood ? K : 0);
  std::vector<double> lg_eta((calc_likelihood && !freeze_topics)
                                 ? (eta.is_scalar() ? K : V * K)
                                 : 0);
  if (calc_likelihood) {
    for (std::size_t k = 0; k < K; k++) lg_alpha[k] = lgamma(alpha[k]);
    if (!freeze_topics) {
      if (eta.is_scalar()) {
        const float* e = eta.column(0);
        for (std::size_t k = 0; k < K; k++) lg_eta[k] = lgamma(static_cast<double>(e[k]));
      } else {
        for (std::size_t v = 0; v < V; v++) {
          const float* e = eta.column(v);
          for (std::size_t k = 0; k < K; k++)
            lg_eta[v * K + k] = lgamma(static_cast<double>(e[k]));
        }
      }
    }
  }

  std::vector<double> prob(K);  // used only for the one-time frozen build below

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

  Progress progress(iterations, verbose);

  // =========================================================================
  // Main loop
  // =========================================================================
  for (std::size_t t = 0; t < iterations; t++) {

    const bool accumulate = averaging && static_cast<int>(t) >= burnin;

    // -----------------------------------------------------------------------
    // DOC PASS: resolve word-proposals, draw doc-proposals
    // -----------------------------------------------------------------------
    run_chunks([&](int c) {
    Scratch& sc = scratch[c];
    long* ck_delta = sc.ck_delta.data();
    for (std::size_t d = chunk_lo(D, c); d < chunk_hi(D, c); d++) {
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
              ck_delta[kp]++;
              ck_delta[k]--;
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
    });

    // -----------------------------------------------------------------------
    // WORD PASS: resolve doc-proposals, draw word-proposals
    // -----------------------------------------------------------------------
    run_chunks([&](int c) {
    Scratch& sc = scratch[c];
    long* ck_delta = sc.ck_delta.data();
    for (std::size_t w = chunk_lo(V, c); w < chunk_hi(V, c); w++) {
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

              ck_delta[kp]++;
              ck_delta[k]--;
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
      // min(n_w, K) nonzeros) plus a dense prior part served by a precomputed
      // table, built once because eta is fixed.
      //
      // THAT ALTERNATIVE WAS BUILT AND MEASURED IN PHASE 7, AND IT IS SLOWER.
      // Do not try it again without reading roadmap 6.7 first.
      //
      // The mixture draw is what the doc pass uses (step 2 there), and porting
      // it here made the sampler 3-4x SLOWER at every K from 10 to 1000, on both
      // the medium and large benchmark corpora. The reason is the thing this
      // build gets right by accident: the table below is constructed immediately
      // before it is used, so it is in L1 for every one of the word's draws. The
      // mixture replaces that cache-hot lookup with old_z(word_token(begin +
      // random)) -- a scattered read into an N-element array, once per token.
      // At 14.7M tokens that is a miss per draw, and no amount of saved setup
      // pays for it. The new code was flat in K at ~1.56 s per iteration, which
      // is the signature of a cost that is pure memory latency.
      //
      // The premise was wrong too. Profiling at -O2 shows per-iteration cost is
      // FLAT in K (0.36-0.52 s from K=10 to K=1000 on the large corpus), so the
      // O(V*K) term this was meant to remove is not a bottleneck at all. An
      // earlier profile that found it growing had been built at -O0, where these
      // dense K-loops do not vectorize and look far more expensive than they are.
      //
      // D15 stands, now on evidence rather than estimate.
      if (!freeze_topics) {
        for (std::size_t k = 0; k < K; k++)
          sc.prob[k] = Cv_w[k] + static_cast<double>(eta_w[k]);
        sc.alias.setup(sc.prob);
      }
      const AliasTable& tbl = freeze_topics ? frozen_alias[w] : sc.alias;

      for (token_t i = begin; i < end; i++) {
        const token_t tok = corpus.word_token(i);
        for (std::size_t m = 0; m < mh_steps; m++) {
          corpus.proposal(tok, m) =
              static_cast<topic_t>(tbl.sample(rng.unif(), rng.unif()));
        }
      }
    }
    });

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

      // The word half of the collapsed joint (see the block comment below):
      //
      //   sum_k [ sum_v lgamma(Cv[v][k] + eta) - lgamma(Ck[k] + eta_bar[k]) ]
      //        + sum_k [ lgamma(eta_bar[k]) - sum_v lgamma(eta) ]
      //
      // accumulated in the form where the two sums over v cancel term by term
      // wherever the count is zero, since lgamma(0 + eta) - lgamma(eta) = 0.
      // Only nonzero counts reach lgamma, which is the expensive part.
      double joint_w = 0.0;
      const bool eta_scalar = eta.is_scalar();
      for (std::size_t v = 0; v < V; v++) {
        const float*  ev = eta.column(v);
        const double* lg = eta_scalar ? lg_eta.data() : &lg_eta[v * K];
        const int*    cv = &Cv[v * K];
        double*       bp = &beta_prob[v * K];
        for (std::size_t k = 0; k < K; k++) {
          const double e = static_cast<double>(ev[k]);
          const double c = static_cast<double>(cv[k]);
          bp[k] = (c + e) / beta_denom[k];
          if (cv[k] != 0) joint_w += lgamma(c + e) - lg[k];
        }
      }
      for (std::size_t k = 0; k < K; k++) {
        joint_w += lgamma(eta.bar(k)) - lgamma(beta_denom[k]);
      }

      double joint_d = 0.0;
      double lpd = 0.0;

      const double lg_alpha_bar = lgamma(alpha_bar);

      for (std::size_t d = 0; d < D; d++) {
        const int* Cd_d = &Cd[d * K];

        double denom = 0.0;
        for (std::size_t k = 0; k < K; k++) denom += Cd_d[k] + alpha[k];
        for (std::size_t k = 0; k < K; k++) {
          theta_prob[k] = (Cd_d[k] + alpha[k]) / denom;
          // Document half of the collapsed joint, same cancellation as above.
          if (Cd_d[k] != 0) joint_d += lgamma(Cd_d[k] + alpha[k]) - lg_alpha[k];
        }
        joint_d += lg_alpha_bar - lgamma(denom);

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
      // Row 3 is the COLLAPSED JOINT, log p(w, z | alpha, eta), with theta and
      // Phi analytically integrated out:
      //
      //   sum_d [ sum_k lgamma(Cd[d][k] + alpha_k) - lgamma(n_d + alpha_bar) ]
      //     + D [ lgamma(alpha_bar) - sum_k lgamma(alpha_k) ]
      //   + sum_k [ sum_v lgamma(Cv[v][k] + eta_kv) - lgamma(Ck[k] + eta_bar_k) ]
      //     + sum_k [ lgamma(eta_bar_k) - sum_v lgamma(eta_kv) ]
      //
      // This replaced a plug-in quantity -- row 2 plus Dirichlet log-densities
      // evaluated at theta_hat and beta_hat -- which was wrong twice over. It
      // had a sign error on both normalizers, and even corrected it was a
      // density, so it was unbounded above and routinely reported positive
      // numbers that meant nothing.
      //
      // The joint is a probability MASS: integrating out theta and Phi leaves a
      // discrete distribution over (w, z), so this is always <= 0. It is also
      // the unnormalized log target of a collapsed sampler, which makes it the
      // ordinary thing to watch for convergence, and marginalizing gives it an
      // Occam factor that row 2 lacks.
      //
      // It conditions on z, so it is a within-model diagnostic. Comparing it
      // ACROSS models is not valid -- that needs log p(w | alpha, eta) with z
      // marginalized too, which is intractable and wants the held-out
      // estimators of Wallach et al. (2009).
      log_likelihood.push_back(static_cast<double>(t));
      log_likelihood.push_back(lpd);
      log_likelihood.push_back(joint_d + joint_w);
    }

    if (verbose) progress.increment();
    RcppThread::checkUserInterrupt();  // main thread only, between iterations
  }

  // =========================================================================
  // Assemble the return value, matching fit_lda_c's contract exactly
  // =========================================================================

  // EXPORT ONLY WHAT new_tidylda() CAN READ. It takes Cd_mean/Cv_mean when
  // burnin > -1 and Cd/Cv otherwise, so exactly one of those pairs is reachable
  // on any given run. Copying out the other pair costs a full D*K and V*K R
  // allocation for nothing: at V = 81k, K = 1000 that is about 1.2 GB of peak
  // memory freed the moment the call returns.
  //
  // The averaging side was already gated this way; the raw side was not, which
  // is the whole of the asymmetry. Both are now.
  //
  // Cd_sum and Cv_sum used to be exported too. No R code has ever read them --
  // Cd_mean and Cv_mean are computed here from the same accumulators -- so they
  // are gone rather than gated.
  const bool raw_out = !averaging;                       // Cd/Cv reachable
  const bool cv_raw  = raw_out && !freeze_topics;
  const bool cv_mean = averaging && !freeze_topics;

  // C_word is current -- the word pass rebuilt and maintained every column and
  // nothing has moved since. C_doc is not: the word pass reassigned tokens after
  // the doc pass last touched it. The likelihood block above happens to leave a
  // fresh C_doc behind, but only when calc_likelihood is TRUE, so rebuild it
  // here. Only needed when Cd is actually going out.
  IntegerMatrix Cd_out(raw_out ? D : 0, raw_out ? K : 0);
  if (raw_out) {
    std::fill(Cd.begin(), Cd.end(), 0);
    for (std::size_t d = 0; d < D; d++) {
      for (token_t i = corpus.doc_begin(d); i < corpus.doc_end(d); i++) {
        Cd[d * K + corpus.old_z(i)]++;
      }
    }
    for (std::size_t d = 0; d < D; d++)
      for (std::size_t k = 0; k < K; k++) Cd_out(d, k) = Cd[d * K + k];
  }

  // D17: C^v goes out as V x K, the orientation the engine already holds it in.
  // The K x V transpose this used to do was a stopgap for the old R contract;
  // R's consumers were rewritten in Phase 6 to match instead. new_tidylda()
  // stores it as a dgCMatrix, and V x K is the orientation that makes each
  // word's topic counts one contiguous column of it.
  //
  // Frozen topics never build C^v, so it goes out empty rather than as a zero
  // matrix that could be mistaken for a real count.
  IntegerMatrix Cv_out(cv_raw ? V : 0, cv_raw ? K : 0);
  if (cv_raw)
    for (std::size_t v = 0; v < V; v++)
      for (std::size_t k = 0; k < K; k++) Cv_out(v, k) = Cv[v * K + k];

  NumericMatrix Cd_mean(averaging ? D : 0, averaging ? K : 0);
  NumericMatrix Cv_mean(cv_mean ? V : 0, cv_mean ? K : 0);

  if (averaging) {
    const double n_post = static_cast<double>(iterations - burnin);
    for (std::size_t d = 0; d < D; d++)
      for (std::size_t k = 0; k < K; k++)
        Cd_mean(d, k) = static_cast<double>(Cd_sum[d * K + k]) / n_post;

    if (cv_mean) {
      for (std::size_t v = 0; v < V; v++)
        for (std::size_t k = 0; k < K; k++)
          Cv_mean(v, k) = static_cast<double>(Cv_sum[v * K + k]) / n_post;
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
    _["log_likelihood"] = ll_out,
    _["alpha"]          = alpha,
    _["eta"]            = eta_in
  );
}
