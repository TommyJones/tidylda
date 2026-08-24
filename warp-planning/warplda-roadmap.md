---
output:
  pdf_document: default
  html_document: default
---

<!--
EDITING THIS FILE: USE ASCII ONLY, AND KEEP INLINE MATH LEGAL.

These documents render to PDF through pandoc -> pdflatex. Two things break it,
and both fail late and confusingly, so check before you commit:

  python3 warp-planning/check-markdown.py warp-planning/*.md

1. NON-ASCII CHARACTERS. pdflatex accepts only Unicode that inputenc's utf8
   option maps to a LaTeX command; anything else fails with "Unicode character
   ... not set up for use with LaTeX". LaTeX stops at the FIRST one, so a single
   reported error usually hides several more.

   Permitted: em dash, en dash, section sign. Nothing else.

   Write symbols as ASCII or as real math. Prefer ASCII for prose punctuation --
   "->" not an arrow glyph, "x" not a multiplication sign, "+/-" not a plus-minus
   sign. Reserve math delimiters for actual mathematics.

   The trap: U+2212 MINUS SIGN is visually near-identical to the ASCII hyphen and
   slips into negative numbers in results tables unnoticed.

2. PANDOC'S INLINE MATH RULES. Violate one and pandoc emits two literal dollar
   signs, putting a bare LaTeX command outside math mode -- "Missing $ inserted".
   Writing the delimiter as "D" below to avoid tripping the checker on this note:

     - an opening D must NOT be followed by whitespace
     - a closing D must NOT be preceded by whitespace
     - a closing D must NOT be followed immediately by a digit

   That last rule is why a plus-minus command wrapped in delimiters and butted
   straight against "5%" fails, while putting a space after the closing
   delimiter would not. It is also why decorative symbols next to numbers belong
   in ASCII: "1.6x" is safe, whereas the same thing written as math is one edit
   away from breaking.
-->
# warpLDA Project Roadmap

**Branch:** `warp` – **Companion document:**
`warp-planning/warplda-design-notes.md` (math, derivations, rationale)

---

## 1. How to use this document

**Read this before doing any work on the `warp` branch.** It is the project's
external memory, written because this work spans many sessions and conversation
context does not survive that long.

Three rules:

1. **Start here.** §2 tells you where the project stands and what to do next.
   §4 tells you what has already been decided.
2. **Update §2 at the end of every session.** A stale status ledger is worse
   than none — it will confidently send the next session in the wrong direction.
3. **§4 decisions are settled. Do not relitigate them.** Each was reached with
   effort and most have non-obvious reasoning. If one genuinely needs reopening,
   change it *and record why in the same edit* — a decision that silently
   reverses is how this project loses coherence.

Rationale for anything in §4 lives in `warp-planning/warplda-design-notes.md`,
cross-referenced by section. All paths in this document are relative to the
repository root.

---

## 2. Status ledger

> **Update this section at the end of every session.**

| | |
|---|---|
| **Last updated** | 2026-08-23 |
| **Branch** | `warp` |
| **Base commit** | `5abaa96` (Phase 0 fixes, merged from `main`) |
| **Current phase** | None. All scheduled phases are done. |
| **Last completed** | **Phase 7: the $O(N)$ word proposal.** Built, measured at 3-4x slower, reverted. §6.7. |
| **In flight** | Nothing. |

**Next action:** release prep for 0.1.0. All scheduled phases are done, and
Phase 7 closed as measured-and-declined without touching `src/`.

Release prep: `cran-comments.md`, a reverse-dependency check, and
submission. **Phase 7** (the $O(N)$ word proposal, §6.7) and **Phase 8** (memory
surgery, §6.8) are both unscheduled; Phase 7 is the cheaper of the two and should
land first.

**Where things stand.** `tidylda()`, `refit()` and `predict()` all call
`fit_lda_warp()`, which initializes, samples and threads. The engine is
statistically validated under both priors and reproducible under `set.seed()`
independent of thread count.

**What the whole project delivered, Phases 2 through 6:**

- **Statistically sound at every phase boundary**: scalar prior 8/8, matrix
  prior 8/8 unpaired and 8/8 paired.
- **About 3x faster than collapsed Gibbs single-threaded** at equal quality,
  and a further 7.6-7.8x on 12 physical cores -- roughly 23x end to end.
- **Reproducible under `set.seed()` at any thread count** (D12), achieved
  structurally rather than by argument.
- **Six pre-existing defects fixed**: `theta` under asymmetric `alpha`, an
  unstable `std::sort` shuffling exchangeable tokens, a per-token
  $O(K\log K)$ sort in initialization, two false statements in
  `tidylda()`'s documentation, and a redundant validation in
  `predict.tidylda()` that rejected valid `method` values.
- **D9, D14 and D20 revised on evidence; D12, D13, D16, D17 and D20 delivered;
  open questions 1 and 3 closed.**

**Known and knowingly accepted:** one tLDA cell is underpowered -- small/$K{=}50$
coherence, about 10 core-hours to close if a future run has budget.

**Background material already absorbed** — no need to re-read:
`ignore/parallel-rng-notes.md` (folded into D12 and D13).
`ignore/nuclear_option/` was examined only far enough to confirm it carries a
stale copy of the Phase 0 `denom` bug (design notes §9); it is a dead-end
experiment, is `.Rbuildignore`d, and needs no further reading.

---

## 3. The project in one page

**Goal.** Replace tidylda's collapsed Gibbs sampler with a warpLDA-based engine
written in Rcpp and parallelized with RcppThread, while preserving tLDA — the
matrix prior $\boldsymbol\eta$ that enables principled transfer learning, where
each topic gets its own Dirichlet prior over words, derived from the scaled
posterior of the model at $t-1$.

**Why.** The current sampler is $O(NK)$ per iteration: it builds a full
$K$-length probability vector for every token. warpLDA is $O(VK + N)$ — $O(K)$
per *word type*, then $O(1)$ per token. Since $V \ll N$ in any real corpus, this
is a large win before parallelism enters the picture at all.

**The non-negotiable.** Statistical parity with the current sampler must be
demonstrated before any of this merges. A previous parallel implementation was
abandoned precisely because it produced objectively worse topic models. That
attempt used a different, hackier approach lacking warpLDA's provable MH
guarantees, so this is not a repeat — but "it compiles and the likelihood curve
looks right" is not evidence of correctness, and the benchmarking harness exists
to make that concrete.

**Why the matrix prior is the interesting part.** warpLDA's efficiency rests on
Metropolis-Hastings proposals whose acceptance ratios cancel in a way that keeps
every memory access cache-resident. The central question for this project was
whether a topic-specific $\boldsymbol\eta$ breaks those cancellations. It does
not — see design notes §3. The change is close to mechanical.

---

## 4. Decision log

Settled. See §1 rule 3 before changing any of these.

| # | Decision | Rationale | Notes § |
|---|---|---|---|
| D1 | Work on branch `warp`, not a fork | Simpler; same repo | — |
| D2 | Port from text2vec `src/mcemlda/` (MIT licensed) as reference | Compact, faithful warpLDA; confirms the derivation line by line. Where it and the design notes disagree, **the derivation wins** — it is the reference for structure, not for correctness | §1 |
| D3 | Store $\boldsymbol\eta$ **dense**, no sparse representation | $\hat\beta^{(t-1)}$ is dense at $t=1$, so tLDA has no sparsity to exploit at any $t$ | §5.2 |
| D4 | Store $\boldsymbol\eta$ **column-major ($V \times K$)** | The word pass needs $\boldsymbol\eta_{\cdot v}$ contiguous beside $C^v_{\cdot v}$; this is the transpose of tidylda's current $K \times V$ | §3.4 |
| D5 | $\boldsymbol\eta$ as `float`, promoted to `double` for computation | Halves the largest allocation; keeps precision drift out of MH accept decisions | §5.4 |
| D6 | $C_k$ stays integer; $\bar\eta_k$ stays `double` | $C_k$ can exceed $2^{24}$ where `float` loses integer exactness | §5.4 |
| D7 | **Drop `optimize_alpha`** | A placeholder for unimplemented fixed-point estimation; carrying it forward adds real complexity for a hack. Consequence: $\boldsymbol\alpha$ is fixed, so its alias table is built once | §6.5 |
| D8 | Keep informed $\hat\beta\cdot\hat\theta$ initialization; discard warpLDA's uniform start | Required for seeded and tLDA models; uniform init also expected to slow convergence | §6.2 |
| D9 | `freeze_topics` (prediction): keep the algorithmic specializations — alias tables built once, $\boldsymbol\eta$/$C^v$/$C_k$ unused — but select them with a **runtime flag inside one kernel**, not a separate kernel | With topics frozen $q_w \propto \hat\beta_{kv}$ is fixed for the whole run, so per-word alias tables are built once instead of per word per iteration; that is where the speed comes from and it is what the flag guards. **Revised 2026-08-22.** The original said "separate specialization, not a branch in the hot loop". Measured at R's actual flags (`-O2`, where `-funswitch-loops` is *off*, so the compiler will not hoist it for you), 36M accept evaluations, min-of-9 interleaved: runtime flag 313.1/356.2 ms (training/frozen), `template<bool>` 313.5/354.3 ms, manual unswitch 314.0/373.3 ms. The branch is free — the loop is two FP divides plus scattered loads per evaluation, and a loop-invariant branch predicts perfectly against that. A second kernel would instead duplicate ~80 lines of MH bookkeeping to change four expressions, which is the code most likely to drift under later patches. Built in **Phase 3** | §6.3 |
| D10 | Burnin: accumulate `Cd_sum` during the doc pass, `Cv_sum` during the word pass | Each pass rebuilds its own matrix from `old_z`, so each is current in its own pass. No reconciliation sweep needed. Valid because these estimate *marginal* expectations, and totals stay exact. **Both matrices must be maintained through their pass's accept step, not merely rebuilt at the start of it** — a rebuilt-only matrix stops matching `old_z` at the first acceptance, and accumulating that biases the posterior mean | §6.4 |
| D11 | Likelihood evaluated every $n$-th iteration, parameterized. **The scheme is settled; the default value is not** — see §7 open question 2 | The likelihood is $O(\text{nnz}\cdot K + VK)$ and would otherwise dominate an $O(VK+N)$ sampler by ~$K$. **This is not thinning** — the chain advances every iteration and every post-burnin iteration still contributes to the count sums. Only a read-only diagnostic runs less often | §7 |
| D12 | RNG seeded **per work item**, not per thread: `seed = f(master, iteration, pass, index)`. **Built in Phase 2**, not Phase 5 | Gives reproducibility *independent of thread count*, not merely at a fixed count. Removes a confound from benchmarking and needs no caveat for CRAN. Master seed drawn from R's stream on the main thread so `set.seed()` governs. *Moved forward 2026-08-22: building it while single-threaded means Phase 5 changes only scheduling, so "Phase 5 at `threads = 1` reproduces Phase 2 bit for bit" becomes a real regression check separating a threading bug from an RNG-change bug* | §8.2 |
| D13 | Expand seeds through `splitmix64`. **Implemented in Phase 2** as `splitmix64` + `xoshiro256++`, inline in `src/warp_rng.h`; no `sitmo` dependency added | xorshift-family generators (including text2vec's `XOR128PLUS`) produce **correlated streams from nearby seeds** — a silent bias that looks like nothing until benchmarks come back subtly wrong. Verified: first draws of adjacent document seeds correlate at r = +0.002 over 200k pairs | §8.3 |
| D14 | Benchmark for **non-inferiority**, not model quality: $R^2$ and mean probabilistic coherence, across multiple seeds, no held-out data. **Unpaired** one-sided test, margin 5% | The question is whether MH is *no worse than* Gibbs on the same model and data, not whether either is good in the abstract. Both metrics ship with tidylda. Pass/fail criterion in §6.1 | §10 |
| D15 | Accept the $O(VK)$ word-proposal construction for now | The $O(N)$ alternative needs $V$ precomputed alias tables over $\boldsymbol\eta$ columns, costing ~2x $\boldsymbol\eta$ in permanent memory. **The code must carry a comment recording this alternative** — it is wanted downstream | §4.2 |
| D16 | The new engine subsumes **both** `create_lexicon()` and `fit_lda_c()`. **Done in Phase 4** | Building the CSR/CSC token structure *is* what `create_lexicon` does; fusing them eliminates the R/C++ round trip and the 16-bytes-per-token marshalling. `create_lexicon()` is now uncalled by the package and survives only as a reference for one test; it goes with `fit_lda_c()` in Phase 6 | §11 |
| D17 | Export $C^d$ and $C^v$ in the engine's own orientation — $C^d$ as $D \times K$, $C^v$ as $V \times K$ — with **no transpose on output**. Rewrite the R consumers to match. **Orientation done in Phase 6; the sparse half landed during release prep, for $C^v$ only** — measurement showed $C^d$ should stay dense. See the invariant in §5 | Supersedes an earlier plan to transpose $C^v$ to topic-major on every fit. Sparse storage shrinks the largest part of the returned object, and keeping the engine's orientation avoids transposing a $V \times K$ matrix on every run. Deferred because it touches the R surface rather than the sampler, and doing it early would churn code the engine work has not stabilised yet. Caveats and the consumer list are in §6.7 | §6.7 |
| D18 | MH steps configurable, default 1 | Default reproduces the reference exactly and costs nothing; the parameter is what allows experimentation with mixing under tLDA's sharper priors. Costs `mh_steps x 2` bytes per token above the default. Built in **Phase 2** | §11.1 |
| D19 | Alias table over $\boldsymbol\alpha$ in the doc-proposal draw — **binding**, Phase 2 | The reference's uniform-draw branch is only proportional to $\alpha_k$ when $\boldsymbol\alpha$ is symmetric; tidylda permits a vector. Omitting it yields code that runs fine and samples from the wrong prior. Costs one $O(K)$ setup, since D7 makes $\boldsymbol\alpha$ fixed | §3.5 |
| D20 | Scalar fast path for $\boldsymbol\eta$ — **done in Phase 6** | A memory win in the common non-transfer case, but an optimization rather than a correctness requirement. *Corrected 2026-08-23: this row previously said Phase 4, contradicting the §6 table. Phase 6 is right — a scalar path computes with a `double` $\eta$ where the matrix path uses the `float`-rounded value (D5), so it moves results and cannot ride along with a refactor whose whole value is being verifiable without a benchmark run.* In the event it was verified by diff after all: `Eta`'s scalar constructor rounds through `float` and *accumulates* $\bar\eta$ over $V$ additions rather than multiplying, which reproduces the matrix path bit for bit | §5.5 |

---

## 5. Invariants — must not break

| Invariant | Enforced by / at risk in |
|---|---|
| `counts$Cd` and `counts$Cv` remain usable as Dirichlet parameters and may hold fractional post-burnin means. **As of 0.1.0** (D17) both are topics-in-columns: $C^d$ is $D \times K$, $C^v$ is $V \times K$. $C^v$ is a sparse `dgCMatrix`; $C^d$ is a dense base matrix, **deliberately** — see §6.8 for the density measurement that split them. $C^v$ also carries dimnames, which is what lets `counts_cv()` identify a post-0.1.0 model when $k$ equals the vocabulary size. $C^v$ carries dimnames, which is what lets `counts_cv()` identify a post-0.1.0 model when $k$ equals the vocabulary size. Every read goes through `counts_cv()`, which transposes pre-0.1.0 saved models | `posterior.tidylda.R:100-133`, `refit.tidylda.R:223` |
| `set.seed()` reproducibility, including under parallelism — a CRAN requirement | D12, D13 |
| Public API unchanged: `tidylda()`, `refit()`, `predict()`, `posterior()`, `tidy`/`augment`/`glance` | all of `R/` |
| `refit`'s R-side vocabulary alignment and topic addition stay in R, untouched | `refit.tidylda.R:249-329` |
| Initialization remains informed, not uniform | D8 |
| R's RNG is touched only on the main thread; workers use `RcppThread::checkUserInterrupt()` | D12 |

---

## 6. Phase plan

| Phase | Deliverable | Exit criterion |
|---|---|---|
| **0** | Defect fixes on `main` | **Done — `5abaa96`** |
| **1** | Statistical benchmarking harness | Produces stable baseline distributions of $R^2$ and mean coherence for the current sampler, across multiple seeds and at least two corpora/$K$ settings |
| **2** | warpLDA engine, single-threaded, **scalar** prior. Includes D18 (`mh_steps`), D19 ($\alpha$ alias table) and — moved forward — D12/D13 (RNG) | **Done.** Matches CGS on the harness at a converged iteration count. Results in §6.3 |
| **3** | Generalize to matrix $\boldsymbol\eta$ (tLDA); D9's `freeze_topics`; `refit()` and `predict()` onto the new engine | **Done.** Whole public API on warpLDA; `fit_lda_c()` has no callers. Gate passes under both a scalar and a matrix prior — §6.3 |
| **4** | Fuse initialization into the engine; eliminate the R round trip | **Done.** One C++ entry point; `Docs`/`Zd` never materialized, so per-token marshalling drops from 16 bytes to zero. **Exit criterion revised** from "identical results to Phase 3" to "identical *initial state*" — see §6.3 |
| **4.5** | Replace `lsamp_one()`'s per-token $O(K\log K)$ sort with a constant-work draw | **Done.** Initialization 7.9x faster at $K{=}10$, 11.8x at $K{=}50$; gate re-run covers Phases 4 and 4.5 together. §6.3 |
| **5** | RcppThread parallelism | **Done.** Bit-identical at 1/2/4/8/16 threads; both gates pass; sampler efficiency 59-68% (about 78% clock-adjusted), 7.8x on 12 physical cores. §6.3 |
| **5.5** | Parallelize initialization | **Done.** Init 3.4x/6.1x/8.8x at $K=10/50/200$; end-to-end 7.6-7.8x on 12 physical cores and flat in $K$. Verified without the gate -- see 6.3 |
| **6** | Cleanup: D17 sparse column-major `counts` and its four consumers; D20 scalar $\boldsymbol\eta$ fast path; documentation, NEWS, CRAN preparation | **Done.** Old engine deleted and its tests rewritten against an independent R reference sampler; D17 and D20 landed; documentation, `NEWS.md` and `DESCRIPTION` rewritten at version 0.1.0; `devtools::check()` clean apart from the local compiler-flag NOTE. §6.4 |
| **6.5** | Triage the open GitHub issues against the 0.1.0 code; close what this project fixed or made moot — see §6.5 | **Done.** Ten issues read; 70 and 77 closed, six commented with status, 8/28/51 untouched. One open question surfaced: whether to surface the log posterior density (issue 30) |
| **6.6** | Replace the plug-in prior-inclusive likelihood with the collapsed joint; delete the sign-buggy row 3 — see §6.6 | **Done.** Agrees with an independent R implementation to 6e-10 (scalar) and 8e-12 (matrix); negative in every case; surfaced as a `log_joint` column. No gate run |
| **7** | The $O(N)$ word proposal — sparse count part plus a precomputed table over the prior — see §6.7 | **Built, measured, reverted.** 3-4x slower at every $K$; the term it targeted is flat in $K$ at `-O2`. No gate run needed — it never got that far |
| **8** | *Unscheduled.* Memory surgery for large corpora — see §6.8 | A separate project. Breaking, broadly |

**The `counts` documentation gap is closed.** `new_tidylda()`'s `@return` now
documents the slot in its topics-in-columns form, with a note that $C^v$ changed
orientation at 0.1.0. *Corrected during release prep: that entry described the
matrices as sparse `dgCMatrix`, which they never were.*

**Why this order.** The two historically risky things — statistical correctness
of the modified sampler, and parallel correctness — are isolated into separate
phases so neither arrives at the same time as anything else. Phase 2
deliberately uses a scalar prior so that a failure implicates the port and not
the matrix prior.

### 6.1 Phase 1 specification

Concrete enough to start from cold. This defines the comparison every later
phase is measured against, so the details are binding, not suggestions.

**Scope of Phase 1 specifically.** The warpLDA engine does not exist yet, so
Phase 1 ran the **Gibbs arm only — 240 fits** (per-$K$ seed counts, see the grid
note below). Phases 2, 3 and 5 run the same 240-fit grid for the warp engine and
compare it against this stored baseline.

**Files.** All under `warp-planning/benchmarks/` — tracked, and excluded from the
build by the `^warp-planning$` entry in `.Rbuildignore`. Not part of the package
and not run by `testthat`.

| File | Role |
|---|---|
| `bench-lib.R` | Shared helpers, sourced by the three scripts below. No side effects |
| `build-corpora.R` | Builds and persists both DTMs. Run once; refuses to overwrite without `--force` |
| `run-benchmark.R` | Runs the grid for one sampler, writes an `.rds`. Parallel across fits, resumable |
| `compare.R` | Loads two `.rds` files, runs the gate, reports pass/fail. `--self-test` validates the gate itself |
| `data/` | Persisted DTMs and a corpus manifest |
| `results/<engine>/` | One `.rds` per fit. Scratch — `.gitignore`d; the assembled `baseline-*.rds` is the artifact |

`bench-lib.R` is a deliberate addition to the three files originally specified
here, so that the schema of a result file is defined and consumed in one place
and `compare.R` does not duplicate the runner's knowledge of it.

**Corpora.** Two, at deliberately different scales:

| Name | Source | Size |
|---|---|---|
| `small` | `nih_sample_dtm` (ships with the package) | 100 docs |
| `medium` | fixed 1,000-row sample of `nih` | 997 docs (see below) |

`nih` is a 68,508 x 44 data frame in `data/nih.rda`, present in the repo but
`.Rbuildignore`d (`^data/nih.rda`), so it is available for benchmarking and
never shipped. Relevant columns are `APPLICATION_ID` and `ABSTRACT_TEXT`. Build
the DTM with `tidytext::unnest_tokens()` + `cast_sparse()`. The `unnest_tokens`
call at `tests/testthat/test-utils.R:7-11` is the pattern to follow; note that
block casts with `cast_dtm` on `nih_sample`, so only the tokenization transfers.

**The medium corpus is 997 documents, not 1,000** — three sampled abstracts
contribute no tokens and drop out at cast time.

> **The sample is drawn once and persisted**, never resampled. `set.seed(8675309)`,
> `dplyr::slice_sample(n = 1000)`, written to
> `warp-planning/benchmarks/data/nih-1000.rds`. If a later session resamples, its
> numbers are not comparable to this baseline and Phase 1 has to be redone.
> `build-corpora.R` refuses to overwrite an existing DTM without `--force` for
> exactly this reason.

**Vocabulary pruning** (fixes $V$, which the whole $O(VK)$ question depends on):
remove `tidytext::stop_words`, then drop terms appearing in fewer than 5
documents. Apply identically to both corpora and record the resulting $V$ in the
results table.

**Fitting hyperparameters.** Identical for every fit — $R^2$ and coherence are
not comparable across different iteration counts:

| Parameter | Value |
|---|---|
| `iterations` | 200 |
| `burnin` | 50 |
| `alpha` | 0.1 (package default) |
| `eta` | 0.05 (package default) |
| `calc_likelihood` | `TRUE` |

200/50 is the maintainer's standard setting for Gibbs on this package.
**Validated in Phase 1 and kept.** Averaged across seeds within each cell, the
fraction of total log-likelihood gain realized by iteration 50 is 0.89–0.94, by
100 is 0.97–0.98, and by 150 is 0.99; the final 25 iterations move the
likelihood by under 1% of total gain in every cell. `burnin = 50` already sits
at 89–94% of the gain, so the post-burnin averaging window draws from a chain
that has substantially converged.

**Grid.** $K \in \{10, 50\}$. Two $K$ values probe whether parity degrades as
topic count grows, which is where the $O(VK)$ term and MH mixing under sharp
priors would both first show up.

**Seeds per cell: 100 at $K=10$, 20 at $K=50$.** Calibrated from the Phase 1
baseline, not guessed — see §6.2. The original uniform 20 leaves coherence
underpowered at $K=10$ in both corpora (48 and 71 seeds required), while $R^2$
would be satisfied by 3 seeds anywhere. Per the margin note below, seeds were
added rather than the margin widened.

> **The statistical risk and the computational risk sit at opposite ends of the
> $K$ range, and this is worth holding onto.** The grid was designed expecting
> trouble at high $K$ — that is where $O(VK)$ and sharp-prior mixing bite. But
> between-seed spread in coherence roughly *doubles* going from $K=50$ to
> $K=10$ (CV 0.046/0.048 against 0.096/0.117), on both corpora, so it tracks
> topic count rather than corpus size. Misspecified $K$ is the high-variance
> regime, and since the margin is *relative* while the $K=10$ coherence mean is
> also lower, the margin shrinks at exactly the $K$ where the noise grows —
> a double penalty. **$K=10$ is therefore where a subtly-worse sampler is
> hardest to detect**, which is an argument for powering that cell properly, not
> for dropping it.

**Metrics.** `calc_lda_r2()`, and the mean of `calc_prob_coherence()` across all
topics — a single number per fit. Both from `R/utils.R`.

**Pass/fail — the merge gate for Phases 2, 3 and 5.** A **one-sided
non-inferiority test** (Welch, unpaired) per cell and per metric, at
$\alpha = 0.05$ with margin $\delta = 0.05\,\overline{\text{gibbs}}$:

$$\text{PASS} \iff \text{reject } H_0:\ \mu_{\text{warp}} - \mu_{\text{gibbs}} \le -\delta$$

Rejecting says warpLDA is not worse than Gibbs by more than the margin.

*Revised 2026-08-21 from the original paired TOST. Recorded per §1 rule 3.*
Three reasons, in order of weight:

1. **Pairing does not survive the phase plan.** Seed $i$ pairs the two engines
   only while they share an initialization path. That holds in Phases 2 and 3,
   but **Phase 4 fuses initialization into the engine** (D16), changing RNG
   consumption, after which seed $i$ no longer produces a common starting state.
   A paired test would silently degrade into two independent samples analyzed as
   paired — still valid, but less powerful, exactly when a power drop would be
   misread as a sampler regression. An unpaired test means the same thing in
   every phase.
2. **Only one direction is a failure.** warpLDA scoring *better* than Gibbs is
   not a reason to block a merge, so the upper bound has no business in the gate.
3. **The margin is what makes the claim testable.** A test of
   $H_0: \mu_{\text{warp}} - \mu_{\text{gibbs}} \le 0$ — "reject that MH is worse"
   — is a *superiority* test: a correct port has a true difference near zero and
   would fail it roughly 95% of the time. Non-inferiority against a margin is the
   formulation that means what it sounds like it means.

The upper side is still computed and **reported as a diagnostic**: a warp engine
that scores implausibly better than Gibbs is more likely one that is not really
sampling than a free win.

Two notes on the margin, both deliberate:

- *Scale.* A relative margin is safe here because average coherence for a
  working model runs around 0.1 or above, and negative topic coherences cluster
  near zero rather than going far negative. If the **Gibbs** baseline mean
  coherence comes out near zero, treat that as a broken baseline — not a margin
  problem — and stop. If **warpLDA** lands near zero, the gate has correctly
  caught a loud failure.
- *Power.* **Resolved by Phase 1.** 20 seeds does not resolve a 5% margin for
  coherence at $K=10$; seeds were added, not the margin widened. The per-cell
  numbers are in §6.2 and the resulting seed counts are above.

**One diagnostic the gate does not cover.** It compares means. A sampler that
mixes worse can match on the mean while being more seed-dependent, and that is
most likely at misspecified $K$, where between-seed spread is already largest.
`compare.R` therefore reports `sd_ratio` alongside each verdict and flags
anything above 1.5. It is not a pass/fail criterion — it is where to look first
if something downstream seems off.

**Persistence.** Write to
`warp-planning/benchmarks/baseline-5abaa96.rds`, named for the commit that
produced it, and fill in the table below so later phases can compare without
re-running and without trusting that a re-run reproduces.

The file stores more than the two metrics, because §6.1 asks for more than the
two metrics: **per-fit likelihood curves** (without which the 200/50 validation
below cannot be done), per-fit wall-clock, the corpus manifest, the
hyperparameters, `sessionInfo()`, and the actual git HEAD. `5abaa96` names the
last commit that touched the sampler; HEAD has since moved on with docs-only
commits, so the two are recorded separately.

### 6.2 Phase 1 baseline results

**Complete.** Gibbs sampler at `5abaa96`, produced on branch `warp` at
`dadf2cd`. Persisted to `warp-planning/benchmarks/baseline-5abaa96.rds`.
240 fits, all `status = ok`.

Corpora after pruning (stop words removed, then terms in fewer than 5 documents
dropped):

| Corpus | $D$ | $V$ | $N$ | nnz |
|---|---|---|---|---|
| small | 100 | 687 | 10,935 | 6,997 |
| medium | 997 | 4,443 | 181,758 | 119,120 |

Metrics, with the 5% margin and the minimum detectable difference at the seed
count actually used. **mdd < margin in every cell — the gate can resolve the
margin everywhere.**

| Corpus | $K$ | $n$ | mean $R^2$ (sd) | mean coherence (sd) | $R^2$ margin / mdd | coherence margin / mdd |
|---|---|---|---|---|---|---|
| small | 10 | 100 | 0.2189 (0.0056) | 0.1217 (0.0129) | 0.01095 / 0.00198 | 0.00609 / 0.00456 |
| small | 50 | 20 | 0.5096 (0.0078) | 0.1637 (0.0078) | 0.02548 / 0.00630 | 0.00819 / 0.00629 |
| medium | 10 | 100 | 0.1000 (0.0020) | 0.1119 (0.0098) | 0.00500 / 0.00071 | 0.00560 / 0.00345 |
| medium | 50 | 20 | 0.2056 (0.0026) | 0.1267 (0.0059) | 0.01028 / 0.00209 | 0.00633 / 0.00472 |

**Coherence is the binding metric.** Its mdd/margin ratio runs 0.62–0.77 against
$R^2$'s 0.14–0.25, so $R^2$ has power to spare everywhere (3 seeds would satisfy
it) while coherence sets the seed count. This is why the seed counts are
asymmetric — see the grid note in §6.1.

**The near-zero-coherence guard did not fire.** The lowest cell mean is 0.112,
and the smallest single-fit value across all 240 fits is well clear of zero.
There is no degenerate cell; in particular small/$K{=}50$, the cell most at risk
on a 100-document corpus, has the *highest* mean coherence of the four (0.1637).

**Timing**, single-threaded per fit, 20 fits in parallel: median 8.5 s
(small/$K{=}10$), 45.8 s (small/$K{=}50$), 143.9 s (medium/$K{=}10$), 781.8 s
(medium/$K{=}50$). The full 240-fit baseline is 8.7 core-hours, about 27 minutes
of wall clock at 20 workers. Phase 2 must budget the same again for the warp arm.

### 6.3 Phase 2 through 5 results — the warpLDA engine

**Both phases complete and passing.** The engine is in `src/warp_rng.h`,
`src/warp_alias.h`, `src/warp_corpus.h`, `src/warp_eta.h` and `src/warp_lda.cpp`.
As of Phase 3, `tidylda()`, `refit()` and `predict()` all dispatch to it.

Phase 2 (scalar prior) is below; Phase 3 (matrix prior, `freeze_topics`, the rest
of the API) follows it.

**Gate result: PASS on all eight cell x metric combinations.** Run
`compare.R baseline-5abaa96.rds run-warp.rds` to reproduce.

| Corpus | $K$ | metric | Gibbs | warp | diff | margin | sd ratio |
|---|---|---|---|---|---|---|---|
| medium | 10 | $R^2$ | 0.1000 | 0.1011 | +0.0011 | 0.0050 | 0.87 |
| medium | 10 | coherence | 0.1119 | 0.1131 | +0.0012 | 0.0056 | 0.82 |
| medium | 50 | $R^2$ | 0.2056 | 0.2064 | +0.0008 | 0.0103 | 0.99 |
| medium | 50 | coherence | 0.1267 | 0.1304 | +0.0038 | 0.0063 | 1.07 |
| small | 10 | $R^2$ | 0.2189 | 0.2230 | +0.0041 | 0.0109 | 0.89 |
| small | 10 | coherence | 0.1217 | 0.1200 | -0.0017 | 0.0061 | 0.91 |
| small | 50 | $R^2$ | 0.5096 | 0.5043 | -0.0052 | 0.0255 | 0.96 |
| small | 50 | coherence | 0.1637 | 0.1598 | -0.0039 | 0.0082 | 0.78 |

`sd_ratio` runs 0.78–1.07, so there is no variance inflation hiding behind a
matched mean — the failure mode that would show up first at misspecified $K$.
The one `check` flag (medium/$K{=}50$ coherence) means only that at $n = 20$ the
gate cannot rule out warp being better by more than 5%; the point estimate is
+3.0%.

#### Iteration counts differ by engine, deliberately

**Gibbs runs 200/50; warp runs 1200/300.** §6.1 originally pinned iterations
across the whole grid, on the sound reasoning that the metrics are not
comparable across iteration counts. That reasoning quietly assumes an iteration
means the same amount of work to both samplers. It does not: collapsed Gibbs
draws every token from its exact conditional, while warpLDA takes one
Metropolis-Hastings proposal per pass and rejects some fraction. Matching
iterations compares *mixing rates*; D14 asks about the *posterior*.

Measured convergence, warp against the Gibbs 200-iteration baseline (5 seeds):

| Corpus / $K$ | $R^2$ @200 -> @600 -> @1200 | coherence @200 -> @600 -> @1200 |
|---|---|---|
| small / 10 | -13.0% -> -2.6% -> **+0.8%** | -13.3% -> +3.4% -> **-1.5%** |
| small / 50 | -15.7% -> -5.1% -> **-0.9%** | -16.3% -> -10.2% -> **-4.3%** |
| medium / 10 | -9.6% -> -0.4% -> **+1.9%** | -23.1% -> -4.5% -> **+1.8%** |
| medium / 50 | -15.7% -> -3.2% -> **+0.5%** | -25.7% -> -6.6% -> **+1.3%** |

Cells inside +/-5% on both metrics: 0/4 at 200 iterations, 2/4 at 600, **4/4 at
1200**. Monotone approach from below in every cell, and warp's log-likelihood
plateaus at the same level as Gibbs (-61939 against -61945 on small/$K{=}10$).
That last point is the load-bearing one: a wrong $\bar\eta$ or a mis-derived
acceptance ratio still gives a valid MCMC chain, it just converges *somewhere
else*. Converging to the same likelihood is what separates those two cases.

`burnin` is held at 25% of iterations, as in Phase 1.

#### The speed win, which is the point of the project

Median seconds per fit, single-threaded:

| Corpus | $K$ | Gibbs @200 | warp @1200 | speedup |
|---|---|---|---|---|
| small | 10 | 8.5 | 5.1 | 1.7x |
| small | 50 | 45.8 | 11.2 | 4.1x |
| medium | 10 | 143.9 | 70.6 | 2.0x |
| medium | 50 | 781.8 | 115.5 | **6.8x** |

**Whole 240-fit grid: 8.73 core-hours for Gibbs, 2.80 for warp — 3.1x overall,
at equal or better quality, with six times the iterations.** Per *iteration* the
gap is far larger (9x to 43x, growing in $K$ exactly as $O(NK)$ against
$O(VK+N)$ predicts); most of it is spent buying convergence back. Parallelism in
Phase 5 multiplies this further.

#### Phase 3 — matrix $\boldsymbol\eta$, `freeze_topics`, and the rest of the API

**Complete.** The whole public API now runs on warpLDA; `fit_lda_c()` has no
callers and is deleted in Phase 6.

**Both gates pass.** Scalar prior, against the stored Gibbs baseline: PASS on all
eight cell x metric combinations, `sd_ratio` 0.83–1.15. Matrix prior
(`tlda-compare.R`): PASS on all eight, and PASS again on the paired test.

| | scalar $\eta$ (main grid) | matrix $\eta$ (tLDA) |
|---|---|---|
| unpaired gate | 8/8 PASS | 8/8 PASS |
| paired gate | n/a — see D14 | 8/8 PASS |
| worst `sd_ratio` | 1.15 | 1.03 |

**Speed**, median seconds per fit, warp @1200 against Gibbs @200:

| Corpus | $K$ | Gibbs | warp | speedup |
|---|---|---|---|---|
| small | 10 | 8.5 | 5.2 | 1.6x |
| small | 50 | 45.8 | 11.4 | 4.0x |
| medium | 10 | 143.9 | 71.6 | 2.0x |
| medium | 50 | 781.8 | 119.9 | **6.5x** |

Whole grid 8.73 -> 2.86 core-hours, **3.05x overall**. The tLDA grid shows lower
ratios (1.2x–5.1x, 2.20x overall) for a reason that has nothing to do with the
matrix prior: it fits **half** the corpus at the **full** vocabulary, so $N$
halves while $VK$ does not, and warp's $O(VK + N)$ cost falls by less than
Gibbs' $O(NK)$. Measured directly on the full corpus, the matrix prior costs
about **5%** over a scalar one (K=50: 5.85s vs 5.58s per 100 iterations).

**The tLDA comparison uses a paired test as well, and this is not a reversal of
D14.** `tlda-compare.R` runs both engines from one set of prepared inputs — same
tokens, same initial assignment, same prior — so the pairing is structural and
survives however initialization is later reorganized. That is a different
experiment from the main grid, where D14's reasoning against pairing stands. The
paired standard deviation is roughly a quarter of the unpaired one, which is
what makes 20 seeds sufficient at $K=50$ where the unpaired test needs far more.

**Seed counts in `tlda-compare.R` are 200 at $K=10$, 20 at $K=50$.** Coherence
means are lower here than in the main grid (half the documents, a strongly
informative transferred prior), so the *relative* margin is smaller and the
unpaired test needs more seeds. Three $K=50$ cells still have `mdd > margin` on
the unpaired test; they pass, and the paired test resolves them comfortably, but
raising $K=50$ to ~100 seeds would close it properly if this is ever re-run.

**Two implementation notes worth keeping.**

*The prior's layout is free.* D4 asks for $\boldsymbol\eta$ "column-major
($V \times K$)". R already stores a $K \times V$ matrix column-major, so all $K$
values for word $v$ are contiguous exactly as the word pass wants. `warp_eta.h`
copies straight through, downcasting to `float`; no data is transposed.

*Read $\boldsymbol\eta$ through a `double`.* It is stored `float` (D5), and the
word pass hoists a `const float*` for the current word. Every read must be cast
before it enters arithmetic: `Cv_w[k] + eta_w[k]` is `int + float` and evaluates
the **entire acceptance ratio in single precision**, which is precisely the drift
D5 exists to prevent. This was introduced and caught during Phase 3; it shifted
accept thresholds by ~2.6e-9 relative, passed every test and both gates, and was
visible only as results moving between runs when nothing that should affect them
had changed.

#### Phase 4 — initialization fused into the engine (D16)

**Complete.** `initialize_topic_counts()` now returns only `beta_initial` and
`Cd_start`; everything per-token happens in C++, in the same walk of the DTM the
sampler uses. `Docs` and `Zd` are never built, so the 16-bytes-per-token round
trip through R is gone entirely — at $N \sim 10^8$ that was ~1.6 GB allocated
twice for no purpose.

Also folded in, at no cost to results: document lengths now come from the
sparse column's nonzeros rather than a dense scan over the whole vocabulary with
a binary-search probe per word ($O(V\log \mathrm{nnz})$ per document, design
notes §9), and the `arma::vec` that `lsamp_one()` needs is allocated once
instead of being converted from a `std::vector<double>` on every token.

**The exit criterion changed, and the reason is a Phase 2 defect.**

The roadmap asked for *identical results to Phase 3*, so that this refactor
could be verified by diff rather than by a 5-core-hour benchmark run. That turned
out to be unreachable, and not because of anything Phase 4 did.

The Phase 2 `Corpus` constructor sorted each document's tokens by word:

```cpp
std::sort(order.begin(), order.end(),
          [&](std::size_t a, std::size_t b) { return docs[d][a] < docs[d][b]; });
```

That sort was **redundant** — `create_lexicon` already emitted tokens
word-ascending — and `std::sort` is **not stable**. Verified: on already-sorted
input with duplicate keys it is the identity at $n = 8$ but not at $n = 40$, 200
or 600, because libstdc++ switches from insertion sort to introsort above 16
elements and quicksort partitioning swaps equal elements. So it permuted tokens
within runs of the same word, each of which carries an independently sampled
topic.

Those tokens are exchangeable — same word, same document, identical
observations — so the shuffle changed nothing about the model. It changed only
which token slot holds which topic, and therefore which per-work-item RNG stream
(D12) reaches which token. Chains diverge from iteration 1.

**Revised exit criterion, which is the stronger claim:** the initialization
draws the same topics and every count matrix the sampler starts from is
identical. Verified against `create_lexicon` on all three:

```
Ck identical: TRUE    Cd identical: TRUE    Cv identical: TRUE
```

Locked by a test in `test-warp-engine.R`, using `iterations = 0` — which now
initializes and returns without sampling. That test calls `create_lexicon()`,
which Phase 6 deletes; keep it as a test-only reference or retire the test
deliberately, but do not let it vanish by accident.

End-state identity was only ever going to hold by accident of that sort, and
could not survive a compiler or libstdc++ change either. Initial-state identity
plus an untouched sampler isolates the difference to something provably
exchangeable.

**Statistical re-validation is still owed.** The gates were not re-run in this
phase, by agreement, because they contend with other work on this machine. The
argument above is why that is safe to defer rather than skip: run
`run-benchmark.R --engine=warp` and `tlda-compare.R` when the machine is free,
and expect both to pass.

#### Phase 4.5 — the initialization sampler

**Complete.** `src/warp_init_sample.h` replaces `lsamp_one()` on the
initialization path.

`lsamp_one()` drew one integer per token by sorting the whole $K$-vector
descending — twice, once for values and once for the index permutation — then
accumulating log-sum-exp across it with a `log_add_exp` per topic, each an `exp`
plus a `log1p`. Roughly three heap allocations and $O(K\log K)$ work per token.

The replacement subtracts the maximum, exponentiates once per topic into a
hoisted buffer while accumulating a running total, and walks that against
$u\cdot\text{total}$. One uniform, $O(K)$, no allocation.

| | $K=10$ | $K=50$ |
|---|---|---|
| initialization, before | 0.66 s | 3.90 s |
| initialization, after | 0.08 s | **0.33 s** |
| speedup | 7.9x | 11.8x |

Medium corpus. Initialization falls from 20% of a 200-iteration fit to about 2%,
and from 4% of a 1200-iteration fit to 0.3%.

**The sort was never load-bearing for correctness.** Accumulating in descending
order is a numerical-stability device; subtracting the maximum is the standard
and stronger one, because it bounds every exponent at zero regardless of order.
That matters most in exactly the tLDA case, where a transferred prior makes one
topic dominate its column by many orders of magnitude. Verified standalone: the
$\chi^2$ statistic is properly distributed over 200 replications (mean 48.84
against 49 expected, p95 65.1 against 66.3); log-weights of $\pm 700$ are handled
where a naive `exp` would overflow or flush to zero; and a dominant topic never
lets a negligible-weight category be drawn.

**A side benefit worth knowing.** The draw is now one `R::unif_rand()` per token
— a *fixed* consumption. `R::rexp()` uses Ahrens-Dieter, which consumes a
variable number of uniforms, so the stream position after initialization used to
depend on the data. It also makes initialization a candidate for the
per-work-item RNG scheme (D12) if it is ever parallelized.

**Test change, made deliberately.** Phase 4's exit test compared the fused
initialization against `create_lexicon()` by exact equality. Phase 4.5 breaks
that premise on purpose, so the test was retired and replaced with one that
checks the property D8 actually cares about and that no future sampler change
should break: initialization tracks the priors it was given rather than being
uniform. Reproducibility under `set.seed()` is tested separately.

**Both gates pass, covering Phases 4 and 4.5 together.**

| gate | result |
|---|---|
| main grid (scalar eta), 8/8 | PASS, `sd_ratio` 0.76-1.02 |
| tLDA (matrix eta), 8/8 unpaired | PASS, `sd_ratio` 0.84-1.07 |
| tLDA, 8/8 paired | PASS |

Differences are essentially zero throughout: r-squared diffs 0.0000 to 0.0007,
coherence -0.0002 to +0.0023.

**K=50 was raised from 20 to 100 seeds in the main grid**, and that resolved the
one cell that failed at n=20. small/K=50 coherence read -0.0061 with p = 0.168 at
20 seeds; at 100 seeds it reads -0.0029 with p = 0.0043 and mdd falls from 0.0063
to 0.0048. Adding seeds moved the estimate *toward* the baseline and sharpened
the test, which is the signature of a noisy estimate rather than a real shift --
a genuine regression would have held its size and become more significant. This
is the second time 6.1's "add seeds, do not widen the margin" has been the right
call.

**Pairing survived Phase 4 better than expected.** The two arms no longer share
an exact initial assignment -- warp builds its own token structure -- so the
paired test was expected to weaken. It did not much: sharing the scenario
(document split, base model, eta at t, Cd_start) carries most of the
correlation. The paired test is 1.9x to 36x sharper than the unpaired one, and
rescues three cells the unpaired test cannot resolve.

**One cell is underpowered and knowingly accepted.** small/K=50 coherence in the
tLDA grid has mdd 0.00542 against a 0.00496 margin even paired. It passes (diff
+0.0023, p = 0.0013) but the gate cannot fully certify it. The tLDA K=50 cells
stayed at 20 seeds because each row runs Gibbs at 200 iterations *plus* warp at
1200; topping up to 100 seeds would cost roughly 10 core-hours. Accepted as a
known limitation rather than paid for.

**Speed, tLDA arm** (both arms timed including initialization, single-threaded):

| Corpus | K | Gibbs | warp | speedup |
|---|---|---|---|---|
| small | 10 | 3.8 s | 3.1 s | 1.20x |
| small | 50 | 21.5 s | 8.3 s | 2.60x |
| medium | 10 | 62.3 s | 35.8 s | 1.74x |
| medium | 50 | 373.4 s | 70.2 s | **5.32x** |

**`src/sample.h` is dead code.** Nothing includes it. Its `log_sample_one()`
also reads `q[k - 1]` at `k = 0`, an out-of-bounds access. Not worth fixing —
delete it in Phase 6 alongside `fit_lda_c()` and `create_lexicon()`.

#### Phase 5 — RcppThread parallelism

**Complete, and the exit criterion holds.** `set.seed()` gives **bit-identical**
`Cd`, `Cv`, `Ck`, `Cd_mean` and log-likelihood curve at 1, 2, 4, 8 and 16
threads. `threads` had been validated in `tidylda_bridge()` and then silently
dropped since Phase 2; it now reaches the engine.

Both gates pass: main grid 8/8, tLDA 8/8 unpaired and 8/8 paired.

The design is in section 7 open question 1. In short: everything except $C_k$
partitions disjointly, $C_k$ is read-only within a pass with per-chunk deltas
merged afterwards, and that makes thread-count independence structural rather
than argued.

**Scaling of the sampler**, medium corpus, 12 physical cores:

| threads | speedup | efficiency |
|---|---|---|
| 4 | 2.79x | 70% |
| 8 | 5.31x | 66% |
| 12 | 7.83x | 65% |
| 24 | 6.44x | *slower than 12* |

Efficiency is flat rather than decaying, which rules out both Amdahl and false
sharing in the sampler -- either would worsen as threads are added. Most of the
shortfall is the processor: all-core turbo is 3556 MHz against 4300 MHz
single-core, so 17 points is clock the code cannot influence. Clock-adjusted
efficiency is about 78%. Twenty-four threads is slower than twelve because the
second SMT sibling on each core adds contention rather than throughput.

**A measurement trap, recorded because it is easy to fall into.** The first pass
at this reported efficiency *degrading* badly -- 56% down to 24% as the corpus
grew, and 36% down to 13% as $K$ grew -- which looked like a hard scaling limit.
It was an artifact. `fit_lda_warp()` does its own initialization, so timing the
whole call includes a **single-threaded** O(N*K) setup, and Amdahl then operates
on that rather than on the sampler. Timing initialization separately
(`iterations = 0`) and subtracting:

| $K$ | init, 1 thread | sampling 1t | sampling 12t | sampling efficiency |
|---|---|---|---|---|
| 10 | 0.98 s | 14.65 s | 1.89 s | 65% |
| 50 | 3.65 s | 16.27 s | 2.01 s | 68% |
| 200 | 16.97 s | 19.06 s | 2.68 s | 59% |

The sampler is stable at 59-68% and does **not** degrade with $N$ or $K$. The
lesson generalizes: do not time a function whose fixed setup cost scales with
the same parameter you are varying.

**What the trap was hiding is real, though, and becomes Phase 5.5.**
Initialization is O(N*K) and single-threaded. At $K = 200$ it is already 17 s
against 19 s of single-threaded sampling -- so perfect parallelism in the
sampler would still cap total speedup near 2x. It is serial because it draws
from R's RNG, which section 5's invariant confines to the main thread. D12's
per-work-item generator is exactly the machinery needed to fix it.

**Phase 5.5 can largely skip the full gate, and here is why.** Phase 4.5
replaced the initialization *algorithm*, which could have carried a bias, so it
was gated. Phase 5.5 changes only the *source of the uniform variate* --
`R::unif_rand()` for the xoshiro stream validated in Phase 2 -- with
`sample_log_weights()` untouched. Same distribution, different draw.

That is a change the existing baseline already covers. It carries 100 seeds per
cell at $K=10$, each with a different initial assignment, and its measured
spread *is* the seed-to-seed population. A new init stream is another draw from
that population, not a new condition.

So verify the two things that are actually new, cheaply:

1. **Identical across thread counts.** The only genuinely new risk, and the same
   property Phase 5 established for the sampler.
2. **The initialization distribution is unchanged.** Compare initial $C^d$ and
   $C^v$ between the old and new generators over ~100 seeds using
   `iterations = 0`, which costs about a second per run. This tests the claim
   directly rather than through 1200 iterations of chain.

Plus the existing informed-init and reproducibility tests and the full suite. If
either check is equivocal, fall back to the gate.

**Load imbalance was checked and is not a limiter.** The word pass splits over a
Zipf-distributed vocabulary, so equal-index chunks are unequal in work by up to
12x. That does not bind, because chunks are dynamically scheduled and there are
four per thread: the ceiling at 48 chunks is 97-100%, against 30-44% if chunks
equalled threads. Keep chunk count several times thread count. The one hard
floor is that a single word cannot be split, capping useful parallelism at the
reciprocal of the most frequent word's share -- measured at about 110 threads on
the medium corpus, and rising with vocabulary size.

#### Phase 5.5 — parallel initialization

**Complete.** Initialization now draws from D12's per-work-item generator rather
than R's, which is what allowed it to be threaded: R's RNG is main-thread-only
(section 5), so while initialization used it the O(N*K) setup stayed serial and
capped total speedup regardless of how well the sampler scaled.

| $K$ | init 1 thread | init 12 threads | speedup |
|---|---|---|---|
| 10 | 0.90 s | 0.27 s | 3.4x |
| 50 | 3.64 s | 0.59 s | 6.1x |
| 200 | 17.23 s | 1.95 s | **8.8x** |

Efficiency rises with $K$ because the remaining serial work -- one O(nnz) pass
for document lengths and `finalize()`'s O(N) counting sort -- is a smaller share
when there is more per-token work to divide. Those two are now the largest
serial pieces of initialization and are the next candidates if it ever shows up
in a profile again.

**End to end, 2.15M tokens, 12 physical cores:**

| $K$ | 1 thread | 12 threads | speedup | efficiency |
|---|---|---|---|---|
| 10 | 76.81 s | 10.09 s | 7.61x | 63% |
| 50 | 86.93 s | 11.10 s | 7.83x | 65% |
| 200 | 117.54 s | 15.30 s | 7.68x | 64% |

Flat in $K$, and matching the sampler's own efficiency -- initialization is no
longer a separate bottleneck. At $K = 200$ the old serial init would have made
the 12-thread total about 30.3 s, or 3.9x; it is now 7.7x.

**`Pass::init` was added to the RNG enum, and it is not decorative.**
Initialization runs conceptually at iteration 0, and the doc pass at iteration 0
already claims `work_item_rng(master, 0, Pass::doc, d)`. Sharing that stream
would have driven a token's starting topic and its very first proposal from the
same uniform -- a correlation with no reason to exist and no obvious symptom.

`Corpus` gained `begin_build()` / `set_token()` in place of the append-only
`add()` / `end_doc()`. Prefix-summing the document lengths first means every
document's slot range is known before any token is written, so documents fill
disjoint regions and need no coordination. The sparse matrix is read through raw
`col_ptrs` / `row_indices` / `values` after an explicit `sync()`, because arma's
sparse iterators can trigger lazy synchronization, which is not safe to invoke
from several threads at once.

**Verified without the full gate, deliberately.** Phase 4.5 replaced the
initialization *algorithm* and was gated because an algorithm can carry a bias.
This changed only the *source of the uniform*, with `sample_log_weights()`
untouched -- same distribution, different draw. That is variation the existing
baseline already characterizes, since it carries 100 seeds per cell at $K = 10$,
each with a different initial assignment. What was checked instead:

1. **Thread-count invariance** -- identical `Cd`, `Cv`, `Ck` at 1, 2, 4, 8 and
   16 threads. The only genuinely new risk.
2. **The draw matches its analytic expectation**, not the old code. The target
   is closed-form, so $E[C^d_{dk}] = \\sum_v n_{dv}\\,p(k \\mid d,v)$ is
   computable directly: correlation 0.999992 over 400 replications, per-cell
   z-scores with max 1.38 and none beyond 4. This is stronger than diffing two
   implementations, which could agree while both being wrong.

#### Notes for later phases

- `optimize_alpha` (D7) is accepted but ignored, with a once-per-session
  warning. Phase 6 turns that into a formal deprecation.
- The two Phase 2 skips are gone; the suite runs with no skips.
- The reference lives at `ignore/text2vec/` (pinned `0b31bdd8`, `.gitignore`d,
  not vendored per D2).

### 6.4 Phase 6 results -- cleanup, documentation, release preparation

No gate run. Nothing in this phase moves the sampler: A deletes dead code, B
changes only the container the counts are returned in, and C is bit-identical to
the matrix path by construction. Any difference here would be a bug, not a
result to re-validate.

**A. The old engine is gone.** `src/lda_gibbs2.cpp` (`fit_lda_c()`,
`create_lexicon()`), `src/sample_int.h`, `src/sample.h` and
`src/parallel_gibbs_utils.h` are deleted, along with `man/fit_lda_c.Rd` and
`man/create_lexicon.Rd`. `src/matrix_conversions.h` stays -- `warp_lda.cpp` uses
`mat_to_vec`. Eleven dangling `\link` targets in `R/utils.R` now point at
`fit_lda_warp`.

The two test files that depended on the old engine were rewritten rather than
retired. `test-cpp_funs.R` now exercises `fit_lda_warp()`: token conservation
across `Cd`, `Cv`, `Cd_mean` and `Cv_mean`, a finite log likelihood, and `alpha`
unchanged -- which now tests D7's removal of `optimize_alpha` rather than its
behaviour. `test-warp-engine.R`'s posterior-targeting test needed a reference
sampler, so it got **an independent collapsed Gibbs sampler written in R**,
about thirty lines over a deliberately tiny corpus. It shares no header, no RNG
and no data structure with the engine, so it can catch a warpLDA bug that a
C++-against-C++ comparison cannot, and it survives future engine changes. The
two agree on corpus log likelihood to 0.15%.

**B. D17 landed.** The engine emits $C^v$ as $V \times K$ and $C^d$ as
$D \times K$, both sparse, with no transpose on output. Four consumers updated,
per design notes §6.7. Every read goes through a new internal `counts_cv()`,
which detects orientation against `nrow(object$beta)` and transposes a
pre-0.1.0 saved model on the way in -- tidylda is on CRAN and those objects
exist in the wild. `test-counts-contract.R` covers the orientation, the
token-topic association through `posterior()`, and an old-format object round
trip.

One assertion in the first draft of that file was wrong and is worth recording:
it required `colSums(Cd_mean) == colSums(Cv_mean)`. Per D10 the two are
accumulated at different points in the iteration, so their topic marginals need
not agree; the test now asserts row marginals and the grand total.

**C. D20 landed.** `Eta` gained a scalar constructor that keeps one $K$-length
array instead of materializing $K \times V$, and `format_eta()` no longer
expands a scalar prior. Verified bit-identical to the matrix path on `beta`,
`theta`, `counts` and `log_likelihood`. At the test scale ($K{=}6$, $V{=}5210$)
the `eta` slot drops from 594.6 Kb to 56 bytes.

Two details make the bit-identity hold, both of which failed on the first
attempt. $\bar\eta$ is **accumulated** over $V$ additions rather than computed
as `v * eta`, because the matrix constructor reaches it by addition and the two
do not land on the same bits: 260.50000000002495 against 260.5. And the scalar
is **rounded through `float`**, matching D5's storage for the matrix path, so a
scalar fit and a matrix fit describing the identical model agree exactly.
Downstream, `new_tidylda()` tests `length(lda$eta) == 1` rather than
`is.matrix()`, since a scalar prior now returns from the engine as a 1x1 matrix.

**D. Documentation.** The sampler is named: **warpLDA**, with the arXiv 1510.08628
citation, once each in `DESCRIPTION`, the package-level doc, `README.md` and
`tidylda()`'s `@details`; "Metropolis-Hastings" thereafter. `predict()`'s
`method` argument is now `c("mh", "dot", "gibbs")` -- `"mh"` is the default,
`"gibbs"` still works and warns once per session. `optimize_alpha` is documented
as deprecated in both `tidylda()` and `refit.tidylda()`.

Three further defects surfaced during the read-through, none of them in the
engine:

- `predict.tidylda()` carried a redundant `stop("method must be one of 'gibbs'
  or 'dot'")` *after* `match.arg()` had already validated and normalized the
  argument -- so it rejected the new `"mh"`. Removed.
- `predict()`'s `threads` was documented as "currently ignored; only
  single-threaded prediction is implemented". It has been passed through to
  `fit_lda_warp()` since Phase 3.
- `vignettes/tidylda-intro.Rmd`'s YAML header ended
  `%\VignetteEncoding{UTF-8}---` on one line, gluing the closing fence onto the
  encoding directive.

`recover_counts_from_probs()` is deleted, along with its man page and the
commented-out call site that recorded it as returning wrong counts.

**E. Release checks.** `DESCRIPTION` is at `0.1.0` with a rewritten Description
field. `NEWS.md`'s top section is retitled and reorganized so the breaking items
lead, with the four Phase 0 entries preserved as its bug-fix section.
`devtools::check(vignettes = TRUE)` is clean: 0 errors, 0 warnings, and only the
local `-mno-omit-leaf-frame-pointer` NOTE. The spelling NOTE is gone --
`inst/WORDLIST` gained the new terms. Vignettes rebuild, which they had not been
doing here before: pandoc is present, just not on `PATH` (§8).

---

### 6.5 Phase 6.5 results -- retiring the issues this project closed

**Done.** All ten open issues read against the 0.1.0 code. Two closed, six
commented, two untouched. The rule was a recorded decision per issue, not a bulk
close, and applying it moved two issues out of the category the sketch had put
them in.

| # | Disposition |
|---|---|
| 70 | **Closed.** `recover_counts_from_probs()` is deleted, so the bug has no code path left |
| 77 | **Closed.** `predict()` parallelizes over documents as asked. Noted that parallel *training* also landed, contradicting the issue's premise that it was ill-advised -- true of the old batched Gibbs implementation, not of the warp engine |
| 30 | **Open**, against the sketch's expectation. See below |
| 34 | **Open.** Status note: `optimize_alpha` is deprecated, and a fixed $\boldsymbol\alpha$ is what lets D19's alias table be built once, so implementing fixed-point estimation is no longer the free addition it looked like when filed |
| 78 | **Open.** Re-scoped against the new baseline: the speedup minibatch training must justify is now roughly 23x larger. Noted that the broken batched implementation it would have built on is deleted |
| 63 | **Open**, but its implementation sketch was stale in two ways: it proposed using `recover_counts_from_probs()`, now deleted, and it says "remove rows from `Cv`" when D17 made $C^v$ tokens-by-topics |
| 68 | **Open**, and marginally stronger: `threads` is now functional and `mh_steps` is new, so there is more computational surface to move into a control object |
| 8, 28, 51 | **Open**, genuinely untouched by this project |

Both 34 and 78 carry a note that time to delivery will be long if it comes at
all, so nobody reads the status update as a commitment.

**Why 30 stayed open.** The sketch listed it as closeable on the grounds that
Phase 0 fixed the normalizer and D11 settled the evaluation interval. Reading the
issue against the code, that is not what it asks for. It has two parts, and the
answers differ:

- *"the likelihood calculation assumes that eta is a vector, not a matrix"* --
  **fixed.** The calculation reads `eta.column(v)` per word and normalizes with
  `eta.bar(k)`, so a tLDA matrix prior is handled correctly.
- *"I'd like to return both likelihoods... the latter has continually been
  problematic for me to calculate, given that its positive"* -- **computed but
  not surfaced.** The engine emits three rows: iteration, `lpd`, and
  `lpd + lp_alpha + lp_eta`. `new_tidylda()` reads the first two and drops the
  third.

**Row 3 also has a sign error, found while answering that question.**
`warp_lda.cpp:370-374` builds the Dirichlet normalizers as
$\sum_v \log\Gamma(\eta) - \log\Gamma(\bar\eta)$, which is $+\log B(\eta)$,
where a Dirichlet log density needs $-\log B(\eta)$. Same for
$\boldsymbol\alpha$. Verified by recomputing row 3 in R: with the code's sign it
matches to the cent, with the correct sign it differs by **81,393 nats** on a
40-document corpus. Nothing user-facing is affected, because `new_tidylda()`
drops the row before it reaches the model object.

An earlier version of this section, and the first comment posted to issue 30,
said the positivity was correct behaviour and stopped there. That was half right.
The quantity is *still* positive once corrected -- a Dirichlet density with
$\eta \ll 1$ is genuinely unbounded above -- but the reported number was
inflated by a defect, not only by the nature of a density. Corrected in both
places.

**The right fix is replacement, not repair.** Row 3 is a plug-in density: it
evaluates Dirichlet densities at point estimates and adds them. That is not what
the LDA literature reports when it reports a prior-inclusive likelihood. That is
the **collapsed joint** $\log p(w, z \mid \boldsymbol\alpha, \boldsymbol\eta)$,
with $\theta$ and $\Phi$ analytically integrated out -- a Dirichlet-multinomial.
Griffiths and Steyvers (2004) report it; Mallet prints it per iteration.

| | Quantity | Sign |
|---|---|---|
| Row 2, surfaced | $\log p(w \mid \hat\theta, \hat\beta)$, plug-in | $\le 0$ |
| Row 3, dropped | plug-in likelihood $+$ Dirichlet log-densities at the estimates | unbounded |
| Literature | $\log p(w, z \mid \boldsymbol\alpha, \boldsymbol\eta)$, marginalized | $\le 0$ |

Three properties make the third the right target:

- It is a probability **mass**, not a density. Integrating out $\theta$ and
  $\Phi$ leaves a discrete distribution over $(w, z)$, so it is bounded above by
  zero and cannot come out positive. That dissolves the original issue.
- It is the unnormalized log target of a collapsed sampler, so monitoring it is
  the ordinary MCMC practice of watching the log posterior up to a constant.
- Marginalizing gives it an automatic Occam factor, so it is comparable across
  $K$. A plug-in likelihood is not: it improves monotonically as topics are added
  and cannot be used to choose $K$ at all.

It is also cheap, because terms where a count is zero collapse into the
constant: $\log\Gamma(0 + \eta) - \log\Gamma(\eta) = 0$. So it evaluates
$\log\Gamma$ over the **nonzero counts alone**. On the model above: row 2 is
$-64{,}132$, row 3 as computed is $+111{,}516$, the collapsed joint is
$-75{,}958$.

Scheduled as **Phase 6.6** (§6.6), before the 0.1.0 release rather than after --
see that section for why.

---

### 6.6 Phase 6.6 results -- replacing the prior-inclusive likelihood

**Scheduled, and it belongs before the 0.1.0 release.** Three reasons, in order
of weight:

1. **0.1.0 already changes the reported log likelihood.** The Phase 0 normalizer
   fix moved every value, and `NEWS.md` already tells users the numbers differ
   from previous versions. Landing a second diagnostic change under the same
   heading costs one more bullet. Landing it in 0.2.0 makes it a fresh surprise
   in a release that otherwise would not have touched the metric.
2. **It ships a defect otherwise.** Row 3's sign error is inert only because
   `new_tidylda()` happens to drop the row. That is not a safety property, it is
   an accident of which index the R side reads.
3. **It is small and additive.** One new column on an existing tibble, no
   signature change, no change to what the sampler draws. Verifiable against an
   independent R implementation on a toy corpus, the same way the reference
   sampler in §6.4 was.

**What to do.**

- **Delete row 3.** Do not repair the sign. It is a plug-in Dirichlet density,
  which is not a quantity anyone asked for and not what the literature reports.
- **Add the collapsed joint** in its place:

$$\log p(w, z \mid \boldsymbol\alpha, \boldsymbol\eta) = \sum_d \left[ \sum_k \log\Gamma(C^d_{dk} + \alpha_k) - \log\Gamma(n_d + \bar\alpha) \right] + D\left[\log\Gamma(\bar\alpha) - \sum_k \log\Gamma(\alpha_k)\right]$$
$$+ \sum_k \left[ \sum_v \log\Gamma(C^v_{vk} + \eta_{kv}) - \log\Gamma(C_k + \bar\eta_k) \right] + \sum_k \left[\log\Gamma(\bar\eta_k) - \sum_v \log\Gamma(\eta_{kv})\right]$$

- **Keep row 2 as it is.** It is the current `log_likelihood` column and users
  may read it. Changing what an existing column means is worse than adding one.
- **Surface both** from `new_tidylda()`: `log_likelihood` unchanged, plus a new
  column for the joint.

**Implementation notes.**

- *The sparse trick is what makes it cheap.* Where $C^d_{dk} = 0$,
  $\log\Gamma(0 + \alpha_k) = \log\Gamma(\alpha_k)$, which is already part of
  the constant term. So accumulate
  $\log\Gamma(c + \alpha_k) - \log\Gamma(\alpha_k)$ over **nonzero counts
  only** and add the constant once. Same for $C^v$. That is $O(\mathrm{nnz})$
  rather than $O(DK + VK)$, and cheaper than row 2's $O(\mathrm{nnz}\cdot K)$.
- *The constant terms already exist.* `lgeta` and `lgalpha` at
  `warp_lda.cpp:365-374` are exactly the two bracketed constants above, with the
  sign flipped. Fix the sign there and they are reusable as-is -- which is also
  the tidiest way to ensure the bug cannot survive the edit.
- *Matrix $\boldsymbol\eta$ needs no special case.* `eta.column(v)` and
  `eta.bar(k)` supply $\eta_{kv}$ and $\bar\eta_k$ directly, and D20's scalar
  path already presents the same interface.
- *It needs a synchronized joint sample*, exactly like row 2, so it goes inside
  the existing likelihood block after the $C^d$ rebuild -- not in a new pass.

**Done, and it went as specified.** The joint agrees with an independent R
implementation -- written from the algebra, in the *direct* form with the full
constants, so agreeing also checks the cancellation rearrangement -- to
$6\times10^{-10}$ relative under a scalar prior and $8\times10^{-12}$ under a
matrix one. The residual is $\boldsymbol\eta$'s `float` storage (D5), not the
rearrangement. Negative in every case tested, as a probability mass must be.
Three cases cover scalar $\boldsymbol\eta$ with symmetric and asymmetric
$\boldsymbol\alpha$, and a matrix prior. No gate run: the sampler is untouched
and this is a read-only diagnostic.

`lgeta` and `lgalpha` are **gone rather than corrected**. The cancelling form
needs neither, since the two sums over $v$ that carried the sign error cancel
term by term against the count sums. Deleting them was the surest way to
guarantee the defect did not survive the edit, and it avoids the catastrophic
cancellation the old form courted: at $K{=}10$, $V{=}5210$, $\eta{=}0.05$ the
discarded constant was about $137{,}000$ nats against a result near
$-76{,}000$.

**Measured cost.** Both rows together are 6.8% of an iteration when evaluated
every iteration ($D{=}100$, $V{=}5210$, $N{=}21{,}585$, $K{=}50$), so 0.68% at
the default `likelihood_every = 10`. An earlier draft of this section claimed the
joint was cheaper than row 2. That was never measured and the asymptotics do not
support it cleanly -- $\mathrm{nnz}(C^v)$ itself grows with $K$, bounded by
$\min(VK, N)$ -- so the claim is withdrawn rather than defended. What is
measured is the combined block above, and that the nonzero form makes 4.2x fewer
$\log\Gamma$ calls than the dense form on the verification corpus.

**What this does not do.** The joint conditions on $z$, so it is a convergence
diagnostic and a within-model quantity, not a model-comparison statistic. Proper
comparison wants $\log p(w \mid \boldsymbol\alpha, \boldsymbol\eta)$ with $z$
marginalized too, which is intractable and needs the held-out estimators of
Wallach et al. (2009). Document the distinction rather than implying otherwise --
misusing this quantity across models is the standard error in the LDA
literature, and Griffiths and Steyvers reached for the harmonic-mean estimator to
avoid it, which has infinite variance and is worse.

---

### 6.7 Phase 7 results -- the $O(N)$ word proposal, built and declined

**Implemented, measured, and reverted.** The change made the sampler **3-4x
slower at every $K$ tested**, and the premise it rested on did not survive being
measured on an optimized build. `src/` is unchanged; the corpus and both
benchmark scripts are kept.

**What was built.** The word pass built a dense $K$-length weight vector and a
Vose alias table per word type per iteration. Phase 7 replaced that with the
mixture the doc pass has always used (`warp_lda.cpp:481-495`): with probability
$n_w/(n_w + \sum_k \eta_{kw})$ copy a uniformly chosen token of $w$, otherwise
draw from the prior -- uniform under a scalar prior, a per-word table built once
under a matrix prior. The mixture reproduces $q_w$ exactly, so this changed the
mechanism and not the distribution.

**The result**, both builds installed with `R CMD INSTALL` at `-O2` and timed
across $K$ (seconds per iteration, 12 threads):

| $K$ | large: before | after | | medium: before | after |
|---|---|---|---|---|---|
| 10 | 0.503 | 1.551 | | 0.0011 | 0.0048 |
| 50 | 0.415 | 1.571 | | 0.0014 | 0.0052 |
| 200 | 0.496 | 1.577 | | 0.0013 | 0.0048 |
| 500 | 0.361 | 1.564 | | | |
| 1000 | 0.524 | 1.487 | | | |

Runs interleaved before/after/before/after, so this is not thermal drift.

**Why it lost, and it is the same reason warpLDA is fast in the first place.**
The alias table is built immediately before it is used, so it sits in L1 for
every one of that word's draws. The mixture trades that cache-hot lookup for
`old_z(word_token(begin + random))` -- a scattered read into an $N$-element array,
once per token. At 14.7M tokens that is a cache miss per draw. The replacement
measured **flat in $K$ at about 1.56 s**, which is the signature of pure memory
latency: it does not care how much arithmetic it skipped.

The doc pass gets away with the same trick because a document's tokens are
**contiguous** in the doc-ordered array, so its random index is a local read. A
word's tokens are not contiguous -- they are reached through the CSC indirection.
The design notes already observed that this indirection is cheap *because it
ascends*; a random draw does not ascend, and that footnote turns out to be the
whole story.

**The premise was also wrong.** At `-O2`, per-iteration cost is **flat in $K$**:
0.503, 0.415, 0.496, 0.361, 0.524 from $K{=}10$ to $K{=}1000$. There is no
$O(VK)$ bottleneck to remove.

An earlier profile in this same phase found the opposite -- cost rising 1.4 to
2.5 s with $R^2 = 0.96$, $p = 0.004$ -- and it was wrong because it was built at
`-O0`. `pkgbuild::compile_dll()` defaults to `debug = TRUE` and
`devtools::load_all()` inherits it, so every timing taken that way describes
unoptimized code. That matters disproportionately here: `-O2` vectorizes the
dense $K$-loops the phase targeted, while the memory-bound MH accept work barely
benefits, so `-O0` inflates exactly the term under test. **Benchmark installed
builds, never `load_all()`.**

**D15 stands, now on evidence.** Its instruction was "revisit only if profiling
at high $K$ shows it dominating." Profiling at high $K$ shows it does not
dominate. The decision log entry needs no change; §4's amendment about D20
making the scalar case cheap is still true and still irrelevant, because the
scalar case is the one that got 3x slower.

**Kept from the phase:**

- `large` corpus (D=68,306, V=81,131, N=14.7M), gitignored but reproducible
  from `build-corpora.R`.
- `profile-vk.R` and `bench-phase7.R`, both carrying their methodology failures
  in the header so they are not repeated.
- A test robustness fix: `test-counts-contract.R` compared a 10-draw Monte Carlo
  posterior mean against the exact mean and demanded exact top-5 agreement,
  where ranks 5 and 6 differ by 11%. It passed only because the RNG happened to
  be in a favourable state at that point in the file; Phase 7's changed draw
  sequence exposed it. Now seeded, with 200 draws.

**Three measurement failures happened here, all mine**, and they are recorded
because each was invisible in its output:

1. Differencing `tidylda()` end to end at 3 and 8 iterations -- a 3 s difference
   between two 220 s runs. Produced a negative slope at $p = 0.39$.
2. Widening to a 20-iteration spread on the engine alone. Still put two of the
   five terms at *negative* cost.
3. Building at `-O0` throughout. Produced a strong, clean, entirely spurious
   $O(VK)$ signal.

The common thread: every one produced plausible-looking numbers. A wrong
measurement does not announce itself, and the only defenses that worked were an
independent second method and a physically impossible value showing up.

---

### 6.8 Phase 8 sketch — memory surgery for large corpora

**Unscheduled, and deliberately separate.** Recorded here so the reasoning is not
lost, not because it is next.

**Why it exists.** The point of this project is not only speed; it is to make
models on much larger corpora feasible. warpLDA addresses the *time* cost. It
does not address the fact that a fitted `tidylda` object is
$O(KV + DK)$ dense, and the $DK$ term grows with the corpus:

| Slot | Size | $K{=}500$, $V{=}10^5$, $D{=}10^4$ | same, $D{=}10^6$ |
|---|---|---|---|
| `beta`, `lambda`, `eta`, `counts$Cv` | $8KV$ each | 400 MB each | 400 MB each |
| `theta`, `counts$Cd` | $8DK$ each | 40 MB each | **4 GB each** |

D17's sparse half was **not built** (its orientation half was), so it is folded
into this phase rather than assumed done. It is also less of a win than it looks,
and the measurement matters because it changes what to build.

Density of the count matrices, measured at $K{=}20$ over 200 iterations on the
`nih_sample_dtm` corpus:

| | storage | $C^v$ | $C^d$ | sparse pays off below |
|---|---|---|---|---|
| `burnin = -1` | integer | **7.7%** | 38.4% | 33% |
| `burnin = 150` | double | **23.2%** | 81.3% | 67% |

`dgCMatrix` stores a double value plus an `int` row index per nonzero, so the
break-even against a dense matrix is 33% for integer counts and 67% for the
post-burn-in doubles.

**So sparsify $C^v$ and leave $C^d$ dense.** D17 prescribed both; only one of
them is a good idea, which is the kind of thing that only shows up when it is
measured rather than reasoned about.

**Done during release prep**, not deferred here. $C^v$ is a `dgCMatrix`: about
3x smaller, which is 19% off a whole fitted model at $K{=}200$ on a
4,443-token vocabulary, and a larger share as $V$ grows. Three consumers needed
adjusting, and they are worth recording because user code will hit exactly the
same three:

- `base::rowSums`/`colSums` do not dispatch on a `dgCMatrix` unless Matrix is
  attached, and tidylda *imports* it. `Matrix::colSums` in `refit()`.
- `t()` has the same problem and falls through to `t.default`. `Matrix::t()` in
  `counts_cv()`.
- `posterior()` builds `dir_par` as sparse-plus-dense, which Matrix promotes to a
  **dense `dgeMatrix`** rather than a sparse one; `generate_sample()`'s
  `as.data.frame()` cannot coerce that. Coerced with `as.matrix()` where
  `dir_par` is built.

The conversion happens at the point of storage rather than earlier, because
`beta` is computed from `t(Cv)` and would hit the `t.default` problem.

This remains a constant-factor improvement and does not change the asymptotics,
because `beta`, `lambda` and `theta` stay dense --- which is what this phase is
actually for.

**The design.** Store only what cannot be recomputed: the sparse counts, dense
$\boldsymbol\eta$, and $\boldsymbol\alpha$. Drop `beta`, `theta` and `lambda`
as stored matrices and derive them on demand:

$$\beta \propto C^v + \boldsymbol\eta, \qquad \theta \propto C^d + \boldsymbol\alpha, \qquad \lambda \text{ from Bayes' rule}$$

That removes the dense $DK$ term entirely — the term that grows with corpus
size — leaving $O(\text{nnz}(C^d) + \text{nnz}(C^v) + KV)$.

**Why it is a separate project.** It is breaking, and broadly so. `beta`,
`theta` and `lambda` are read by `summarize_topics()`, `print()`, `glance()`,
`tidy()`, `augment()`, `predict(method = "dot")` (via `lambda`), and
`refit()` (via `beta`, when constructing $\boldsymbol\eta^{(t)}$). Every one of
those becomes a computation rather than a lookup.

**One design question worth holding onto.** An S3 `$.tidylda` method or an active
binding could compute these on access and preserve `model$beta` syntax, making
the change non-breaking at the call site. The cost is that a user writing
`model$beta` inside a loop silently pays a full recomputation each time, which
argues for memoisation, which reintroduces the memory. Not resolved; just worth
not rediscovering from scratch.

**Prerequisites.** D17 (sparse counts) and D20 (scalar $\boldsymbol\eta$ fast
path) both point this direction and landed in Phase 6. Phase 7 (§6.7) should also
go first: it is independent, cheaper, and under a matrix prior the two pull in
opposite directions on memory.

---

## 7. Open questions

1. ~~**Work partitioning for parallelism.**~~ **Resolved in Phase 5.** Both
   passes partition on their own index, and everything except $C_k$ is disjoint:
   the doc pass owns slices of `Cd` and each document's tokens, the word pass the
   same for `Cv` and each word's tokens, and eta and the alias tables are
   read-only.

   $C_k$ is read-only within a pass. Each chunk accumulates its own delta and the
   deltas are summed in afterwards. That makes thread-count independence
   structural rather than argued: every work item sees the same $C_k$ however
   chunks land on threads; integer addition is associative and exact, so merge
   order cannot matter; and a work item's contribution depends only on its own
   state and the snapshot, so even the chunk count is free to be tuned for load
   balance. Verified bit identical at 1, 2, 4, 8 and 16 threads.

   Recorded because it cost effort to establish: **atomics do not solve this.**
   They remove the race but leave what each work item *reads* dependent on
   interleaving, so results would drift with thread count -- which is exactly
   what D12 forbids. The requirement is determinism, not merely safety.

   **What limits scaling, measured.** Efficiency is roughly flat at 65-70% from
   2 to 12 threads rather than decaying, which rules out both Amdahl and false
   sharing -- either would worsen as threads are added. Most of the shortfall is
   the processor: all-core turbo is 3556 MHz against 4300 MHz single-core, so 17
   points of the gap is clock the code cannot influence. Clock-adjusted
   efficiency is about 78%. Peak speedup is 7.8x on 12 physical cores; 24
   threads is *slower* than 12, because the second SMT sibling on each core adds
   contention rather than throughput.

   **Memory is not the limiter**, which was checked rather than assumed after an
   earlier draft of this entry claimed otherwise. Holding V and K fixed and
   growing the corpus until the token array crosses the 19.3 MB L3, per-token
   cost is flat -- 250.3, 223.8, 218.3 and 218.4 ns at 3.1, 6.1, 12.1 and 24.0
   MB. The word pass's one scattered access, the CSC indirection
   `z_[csc_token_[i]]` that the reference also has at `LDA.hpp:102`, is cheap
   because it *ascends* (the counting sort that builds `csc_token_` is stable, so
   within a word the indices increase) and because even assuming every access
   misses, N*64 bytes is about 380 MB and 19 ms at 6M tokens against roughly
   1.3 s of work per iteration. Streaming $C^v$ and eta is the same story at
   scale: about 400 MB per iteration, or 20 ms, against seconds of compute.

   **Report per-core and as-shipped numbers separately.** The Gibbs baseline
   was single-threaded, so any parallel comparison is parallel-warp against
   serial-Gibbs. That is the number users experience and is worth publishing,
   but it is not a per-core comparison of samplers and should not be presented
   as one.

2. **Likelihood evaluation interval.** Default of 10 proposed. **Phase 1
   evidence supports it**; confirm against the warpLDA curve in Phase 2 before
   fixing the default. Essentially all readable structure is in the first ~50
   iterations (89–94% of total gain). Past iteration 150 the curve is
   noise-dominated: iteration-to-iteration movement is comparable to the
   remaining drift, with drift/jitter ratios of 1.6–7.1. Sampling every 10th
   iteration loses no real structure and yields a *cleaner* curve than every
   iteration does.
3. ~~**Two expiring doc comments.**~~ **Both retired in Phase 5.**
   `R/tidylda-fit-methods.R`'s `@details` no longer claims the log likelihood
   returns positive numbers (it never did for the value `new_tidylda()` surfaces,
   which is `lpd` and negative), and no longer says parallelism is unimplemented.
   The block now documents the likelihood evaluation interval, `threads`, and
   that `optimize_alpha` is ignored. `man/tidylda.Rd` regenerated.

---

## 8. Environment and workflow

Things that cost time to discover once and should not cost it twice.

**Dependencies.** `RcppThread` (LinkingTo), `mvrsquared` and `tidytext`
(Imports) are required to build or load at all. `quanteda`, `tm`, `slam`,
`spelling` are Suggests exercised by the test suite — CI installs them, so
install them locally too or you verify less than CI does.

**Rendering these documents to PDF.** Both `.md` files here carry a
`pdf_document` header. Under `pdflatex` (what tinytex uses by default), only the
Unicode that `inputenc`'s `utf8` option maps to a LaTeX command will render;
anything else fails with *"Unicode character ... not set up for use with
LaTeX"*, and LaTeX stops at the first one, so a single reported error usually
hides several more.

Safe to use in prose: em dash, en dash and the section sign, which are in the
LaTeX core. **Not safe, and all previously present here** (named rather than
shown, since writing them would reintroduce the problem): U+2212 MINUS SIGN
(visually almost identical to the ASCII hyphen, and the one that surfaced
first), U+2192 RIGHTWARDS ARROW, any bare Greek letter such as U+03B7, and
U+00D7 MULTIPLICATION SIGN, U+00B7 MIDDLE DOT and U+00B1 PLUS-MINUS SIGN, which
need `textcomp` rather than the kernel.

All have been converted: minus signs to ASCII hyphens, the rest to math mode
(`\rightarrow`, `\eta`, `\times`, `\cdot`, `\pm` between dollar signs). Keep
it that way — write mathematical symbols as math, not as Unicode glyphs. To
check before rendering:

```sh
grep -nP '[^\x00-\x7F]' warp-planning/*.md | grep -vP '[\x{2013}\x{2014}\x{00A7}]'
```

**Local checks need `vignettes = FALSE`.** Pandoc is not installed on this
machine and `sudo` requires a password, so vignette building fails at the
*build* stage and aborts the whole check:

```r
devtools::check("/home/tommy/tidylda", document = FALSE, vignettes = FALSE)
```

CI has Pandoc and does build vignettes — confirmed `checking re-building of
vignette outputs ... OK`. Note `vignettes/tLDA.Rmd` has a malformed YAML header
(line 9 glues `---` onto the `%\VignetteEncoding{UTF-8}` line) which R tolerates.

**Known-constant check output.** One NOTE, a spelling diff flagging `ORCID` and
`tidylda's`. Locally there is a second NOTE for the non-portable compile flag
`-mno-omit-leaf-frame-pointer`; this comes from this Ubuntu R build's default
`CXXFLAGS` and does **not** appear on CI. Anything beyond these two is new and
is yours.

**CI.** `.github/workflows/R-CMD-check.yaml` runs `--as-cran` across five
OS/R combinations and installs Suggests. Jobs occasionally hang in `setup-r` and
hit GitHub's 6-hour timeout, surfacing as `cancelled` rather than `failure` —
that is infrastructure, not code.

**Git.** `origin` is SSH (`git@github.com:TommyJones/tidylda.git`). A global
`url."git@github.com:".insteadOf "https://github.com/"` rule rewrites HTTPS
GitHub URLs, so remotes created by GitHub Desktop work without intervention —
but `git remote -v` may still *display* an https URL. That is expected.

**Release conventions.** `NEWS.md` is prepended per version; `cran-comments.md`
is rewritten wholesale at release-prep time; `CRAN-SUBMISSION` is generated by
`devtools` and never hand-edited. `DESCRIPTION` reads `0.1.0` as of Phase 6 --
a **minor** bump, not the patch Phase 0 anticipated, because this release
replaces the sampler, changes the `counts` contract and deprecates two
arguments. `NEWS.md`'s Phase 0 entries live on under `0.1.0` as its bug-fix
section. **0.1.0 has not shipped**, so everything through Phase 7 lands in it and
`NEWS.md` edits cost nothing -- there is no released version whose notes are
being rewritten.

---

## 9. File map

**The engine (to be replaced).**

| Path | Role |
|---|---|
| `src/lda_gibbs2.cpp` | `create_lexicon()` (DTM to token chain, plus informed init) and `fit_lda_c()` (the sampler). Both subsumed by the new engine (D16) |
| `src/sample.h` | `sample_one()`, `log_sample_one()` — already RNG-agnostic, taking the variate as an argument. This is what makes D12 implementable without touching sampling logic |
| `src/matrix_conversions.h` | `mat_to_vec()` / `vec_to_mat()`; note `vec_to_mat` maps outer index to column |
| `src/parallel_gibbs_utils.h` | Batching helpers from the abandoned parallel attempt. Dead; goes away |

**The R surface (mostly unchanged).**

| Path | Role |
|---|---|
| `R/tidylda-fit-methods.R` | `tidylda()` and `tidylda_bridge()` — the main entry point |
| `R/refit.tidylda.R` | tLDA. Constructs $\boldsymbol\eta^{(t)}$ at `:225`; vocabulary and topic alignment at `:249-329` |
| `R/predict.tidylda.R` | `freeze_topics` path |
| `R/utils.R` | `format_eta()`, `format_alpha()`, `initialize_topic_counts()`, `new_tidylda()`, `calc_lda_r2()`, `calc_prob_coherence()` |
| `R/posterior.tidylda.R` | Reads `counts$Cd`/`counts$Cv` directly — see §5 |

**Reference and background.**

| Path | Role |
|---|---|
| `warp-planning/warplda-design-notes.md` | Derivations, cost analysis, constraints. The "why" behind §4 |
| `vignettes/tLDA.Rmd` | The tLDA model and the $a^{(t)} \rightarrow \omega_k^{(t)}$ derivation |
| text2vec `src/mcemlda/LDA.hpp` | Reference warpLDA implementation, MIT licensed. **Not vendored.** Clone `git@github.com:dselivanov/text2vec.git` and check out **`0b31bdd81f37baaf0bd2c8113cbaa578450c8730`** (2025-12-01) — every `LDA.hpp` line number cited in the design notes refers to that commit and will drift on any other |
| arXiv 1510.08628 | The warpLDA paper; also the source of tidylda's likelihood formula |
