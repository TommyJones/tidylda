---
output:
  pdf_document: default
  html_document: default
---
# warpLDA Project Roadmap

**Branch:** `warp` · **Companion document:**
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
| **Last updated** | 2026-08-22 |
| **Branch** | `warp` |
| **Base commit** | `5abaa96` (Phase 0 fixes, merged from `main`) |
| **Current phase** | Phase 3 — matrix $\boldsymbol\eta$ (tLDA) and `freeze_topics` |
| **Last completed** | **Phase 2: warpLDA engine built, single-threaded, scalar prior. Gate PASSES all eight cell × metric combinations.** Results, timings and three reference defects in §6.3. |
| **In flight** | Nothing. |

**Next action:** Phase 3 — generalize the engine to a matrix $\boldsymbol\eta$
(D3–D6 govern layout), add D9's `freeze_topics` specialization, and point
`refit()` and `predict()` at the new engine. Then re-run
`run-benchmark.R --engine=warp` and `compare.R`, and additionally exercise a
transferred prior, which the current grid does not.

**Where the engine is:** `src/warp_rng.h` (D12/D13), `src/warp_alias.h`,
`src/warp_corpus.h` (dual CSR/CSC view), `src/warp_lda.cpp` (`fit_lda_warp`).
`fit_lda_c()` is untouched and still serves `refit()` and `predict()`; it is
deleted in Phase 6 once nothing calls it.

**What Phase 2 settled:**

- **The port is statistically sound.** Gate PASSes everywhere, `sd_ratio`
  0.78–1.07, and warp's log-likelihood plateaus at the same level as Gibbs —
  the check that distinguishes "right posterior" from "some other stationary
  distribution."
- **Iteration counts now differ by engine: Gibbs 200/50, warp 1200/300.** An
  iteration is not the same unit of work to both samplers. §6.3 has the
  convergence table and the reasoning.
- **warp is 3.1× faster over the whole grid** (2.80 vs 8.73 core-hours) at equal
  quality, despite six times the iterations; 1.7×–6.8× per fit, growing in $K$.
- **The reference has a real bug in `beta_bar`** (§6.3(a)) that a faithful port
  would have inherited silently. Do not trust `LDA.hpp` over the design notes.
- **D9 now scheduled** (Phase 3); **D12/D13 moved forward** into Phase 2.

**Three things Phase 3 must not forget:**

1. **Remove the two `skip()`s** — `test-tidylda-fit-methods.R` (vector
   $\boldsymbol\eta$) and `test-refit.tidylda.R` (matrix $\boldsymbol\eta$).
   Both name Phase 3 in the skip message.
2. **Drop the scalar-$\boldsymbol\eta$ guard** in `tidylda_bridge()`, which
   currently errors on anything else.
3. **$\bar\eta_k$ becomes a $K$-vector**, computed once from
   $\sum_v \eta_{kv}$. It is a scalar today only because the prior is.

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
| D2 | Port from text2vec `src/mcemlda/` (MIT licensed) as reference | Compact, faithful warpLDA; confirms the derivation line by line. **But do not port it blindly** — three defects found in Phase 2, listed in §6.3 | §1 |
| D3 | Store $\boldsymbol\eta$ **dense**, no sparse representation | $\hat\beta^{(t-1)}$ is dense at $t=1$, so tLDA has no sparsity to exploit at any $t$ | §5.2 |
| D4 | Store $\boldsymbol\eta$ **column-major ($V \times K$)** | The word pass needs $\boldsymbol\eta_{\cdot v}$ contiguous beside $C^v_{\cdot v}$; this is the transpose of tidylda's current $K \times V$ | §3.4 |
| D5 | $\boldsymbol\eta$ as `float`, promoted to `double` for computation | Halves the largest allocation; keeps precision drift out of MH accept decisions | §5.4 |
| D6 | $C_k$ stays integer; $\bar\eta_k$ stays `double` | $C_k$ can exceed $2^{24}$ where `float` loses integer exactness | §5.4 |
| D7 | **Drop `optimize_alpha`** | A placeholder for unimplemented fixed-point estimation; carrying it forward adds real complexity for a hack. Consequence: $\boldsymbol\alpha$ is fixed, so its alias table is built once | §6.5 |
| D8 | Keep informed $\hat\beta\cdot\hat\theta$ initialization; discard warpLDA's uniform start | Required for seeded and tLDA models; uniform init also expected to slow convergence | §6.2 |
| D9 | `freeze_topics` (prediction) as a **separate specialization**, not a branch in the hot loop. **Built in Phase 3** | With topics frozen $q_w \propto \hat\beta_{kv}$ is fixed for the whole run, so alias tables build once; $\boldsymbol\eta$, $C^v$, $C_k$ go unused. *Phase assigned 2026-08-22: the §6 table never scheduled it, and it is what `predict()` needs. Pairing it with Phase 3's matrix $\boldsymbol\eta$ puts the whole public API on the new engine at the end of one phase* | §6.3 |
| D10 | Burnin: accumulate `Cd_sum` during the doc pass, `Cv_sum` during the word pass | Each pass rebuilds its own matrix from `old_z`, so each is current in its own pass. No reconciliation sweep needed. Valid because these estimate *marginal* expectations, and totals stay exact. **Implementation caveat (Phase 2):** the reference maintains `C_word` through the word pass's accept step (`LDA.hpp:116-117`) but **never updates `C_doc`** in the doc pass (`:168-178`), because it never reads it again. `C_doc[d]` therefore matches `old_z` only until the first acceptance. We add the two missing updates and accumulate at the end of each document; without them the natural placement silently averages a stale matrix | §6.4 |
| D11 | Likelihood evaluated every $n$-th iteration, parameterized. **The scheme is settled; the default value is not** — see §7 open question 2 | The likelihood is $O(\text{nnz}\cdot K + VK)$ and would otherwise dominate an $O(VK+N)$ sampler by ~$K$. **This is not thinning** — the chain advances every iteration and every post-burnin iteration still contributes to the count sums. Only a read-only diagnostic runs less often | §7 |
| D12 | RNG seeded **per work item**, not per thread: `seed = f(master, iteration, pass, index)`. **Built in Phase 2**, not Phase 5 | Gives reproducibility *independent of thread count*, not merely at a fixed count. Removes a confound from benchmarking and needs no caveat for CRAN. Master seed drawn from R's stream on the main thread so `set.seed()` governs. *Moved forward 2026-08-22: building it while single-threaded means Phase 5 changes only scheduling, so "Phase 5 at `threads = 1` reproduces Phase 2 bit for bit" becomes a real regression check separating a threading bug from an RNG-change bug* | §8.2 |
| D13 | Expand seeds through `splitmix64`. **Implemented in Phase 2** as `splitmix64` + `xoshiro256++`, inline in `src/warp_rng.h`; no `sitmo` dependency added | xorshift-family generators (including text2vec's `XOR128PLUS`) produce **correlated streams from nearby seeds** — a silent bias that looks like nothing until benchmarks come back subtly wrong. Verified: first draws of adjacent document seeds correlate at r = +0.002 over 200k pairs | §8.3 |
| D14 | Benchmark for **non-inferiority**, not model quality: $R^2$ and mean probabilistic coherence, across multiple seeds, no held-out data. **Unpaired** one-sided test, margin 5% | The question is whether MH is *no worse than* Gibbs on the same model and data, not whether either is good in the abstract. Both metrics ship with tidylda. Pass/fail criterion in §6.1 | §10 |
| D15 | Accept the $O(VK)$ word-proposal construction for now | The $O(N)$ alternative needs $V$ precomputed alias tables over $\boldsymbol\eta$ columns, costing ~2× $\boldsymbol\eta$ in permanent memory. **The code must carry a comment recording this alternative** — it is wanted downstream | §4.2 |
| D16 | The new engine subsumes **both** `create_lexicon()` and `fit_lda_c()` | Building the CSR/CSC token structure *is* what `create_lexicon` does; fusing them is what eliminates the R/C++ round trip and the 16-bytes-per-token marshalling | §11 |
| D17 | Export $C^d$ and $C^v$ as **sparse** matrices in the engine's own orientation — $C^d$ as $D \times K$, $C^v$ as $V \times K$ — with **no transpose on output**. Rewrite the R consumers to match. **Deferred to Phase 6**, after the engine works | Supersedes an earlier plan to transpose $C^v$ to topic-major on every fit. Sparse storage shrinks the largest part of the returned object, and keeping the engine's orientation avoids transposing a $V \times K$ matrix on every run. Deferred because it touches the R surface rather than the sampler, and doing it early would churn code the engine work has not stabilised yet. Caveats and the consumer list are in §6.7 | §6.7 |
| D18 | MH steps configurable, default 1 | Default reproduces the reference exactly and costs nothing; the parameter is what allows experimentation with mixing under tLDA's sharper priors. Costs `mh_steps × 2` bytes per token above the default. Built in **Phase 2** | §11.1 |
| D19 | Alias table over $\boldsymbol\alpha$ in the doc-proposal draw — **binding**, Phase 2 | The reference's uniform-draw branch is only proportional to $\alpha_k$ when $\boldsymbol\alpha$ is symmetric; tidylda permits a vector. Omitting it yields code that runs fine and samples from the wrong prior. Costs one $O(K)$ setup, since D7 makes $\boldsymbol\alpha$ fixed | §3.5 |
| D20 | Scalar fast path for $\boldsymbol\eta$ — **deferred to Phase 4**, not Phase 2 | A memory win in the common non-transfer case, but it is an optimization, not a correctness requirement, and Phase 2 has enough moving parts. `format_eta()` keeps materializing $K \times V$ until then | §5.5 |

---

## 5. Invariants — must not break

| Invariant | Enforced by / at risk in |
|---|---|
| `counts$Cd` and `counts$Cv` remain usable as Dirichlet parameters and may hold fractional post-burnin means. **Orientation changes in Phase 6** (D17): until then $C^v$ is $K \times V$ dense; after, $V \times K$ sparse. Whichever holds, `refit.tidylda.R:223`, `posterior.tidylda.R:100`/`:122` and `utils.R:642-646` must agree with it | `posterior.tidylda.R:100-133`, `refit.tidylda.R:223` |
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
| **2** | warpLDA engine, single-threaded, **scalar** prior. Includes D18 (`mh_steps`), D19 ($\alpha$ alias table) and — moved forward — D12/D13 (RNG) | **Done.** Matches CGS on the harness at a converged iteration count. Results and the three reference defects in §6.3 |
| **3** | Generalize to matrix $\boldsymbol\eta$ (tLDA); **plus D9's `freeze_topics` specialization**, which `predict()` needs. Remove the two Phase 2 skips in `test-refit.tidylda.R` and `test-tidylda-fit-methods.R`, and point `refit()`/`predict()` at the new engine | Parity preserved when a transferred prior is in play; `refit()` and `predict()` paths exercised; whole public API on the new engine |
| **4** | Fuse initialization into the engine; eliminate the R round trip | Identical results to Phase 3; one C++ entry point; per-token memory down from 16 bytes |
| **5** | RcppThread parallelism | Parity holds; `set.seed()` gives identical results across *different* thread counts (D12) |
| **6** | Cleanup: D17 sparse column-major `counts` and its four consumers; D20 scalar $\boldsymbol\eta$ fast path; documentation, NEWS, CRAN preparation | `devtools::check()` clean; `posterior()` and `refit()` verified against the new orientation; `counts` documented (see below); the expiring comments in §7 rewritten; `man/` regenerated |
| **7** | *Unscheduled.* Memory surgery for large corpora — see §6.4 | A separate project, after the engine lands |

**Documentation gap to close in Phase 6.** `counts` is a real slot on the
`tidylda` object but `new_tidylda()`'s `@return` never mentions it — it documents
`beta`, `theta`, `lambda`, `alpha`, `eta`, `summary`, `call`, `log_likelihood`
and `r2`, and stops. That is a documentation bug, not a deliberate signal that
the slot is private, and it should be fixed. Do it **after** D17 settles the
orientation, so the documentation describes the final sparse $V \times K$ form
rather than the interim one.

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

`nih` is a 68,508 × 44 data frame in `data/nih.rda`, present in the repo but
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

### 6.3 Phase 2 results — warpLDA engine, scalar prior

**Complete and passing.** The engine is in `src/warp_rng.h`, `src/warp_alias.h`,
`src/warp_corpus.h` and `src/warp_lda.cpp`. `tidylda()` dispatches to it;
`refit()` and `predict()` still call `fit_lda_c()` until Phase 3.

**Gate result: PASS on all eight cell × metric combinations.** Run
`compare.R baseline-5abaa96.rds run-warp.rds` to reproduce.

| Corpus | $K$ | metric | Gibbs | warp | diff | margin | sd ratio |
|---|---|---|---|---|---|---|---|
| medium | 10 | $R^2$ | 0.1000 | 0.1011 | +0.0011 | 0.0050 | 0.87 |
| medium | 10 | coherence | 0.1119 | 0.1131 | +0.0012 | 0.0056 | 0.82 |
| medium | 50 | $R^2$ | 0.2056 | 0.2064 | +0.0008 | 0.0103 | 0.99 |
| medium | 50 | coherence | 0.1267 | 0.1304 | +0.0038 | 0.0063 | 1.07 |
| small | 10 | $R^2$ | 0.2189 | 0.2230 | +0.0041 | 0.0109 | 0.89 |
| small | 10 | coherence | 0.1217 | 0.1200 | −0.0017 | 0.0061 | 0.91 |
| small | 50 | $R^2$ | 0.5096 | 0.5043 | −0.0052 | 0.0255 | 0.96 |
| small | 50 | coherence | 0.1637 | 0.1598 | −0.0039 | 0.0082 | 0.78 |

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

| Corpus / $K$ | $R^2$ @200 → @600 → @1200 | coherence @200 → @600 → @1200 |
|---|---|---|
| small / 10 | −13.0% → −2.6% → **+0.8%** | −13.3% → +3.4% → **−1.5%** |
| small / 50 | −15.7% → −5.1% → **−0.9%** | −16.3% → −10.2% → **−4.3%** |
| medium / 10 | −9.6% → −0.4% → **+1.9%** | −23.1% → −4.5% → **+1.8%** |
| medium / 50 | −15.7% → −3.2% → **+0.5%** | −25.7% → −6.6% → **+1.3%** |

Cells inside ±5% on both metrics: 0/4 at 200 iterations, 2/4 at 600, **4/4 at
1200**. Monotone approach from below in every cell, and warp's log-likelihood
plateaus at the same level as Gibbs (−61939 against −61945 on small/$K{=}10$).
That last point is the load-bearing one: a wrong $\bar\eta$ or a mis-derived
acceptance ratio still gives a valid MCMC chain, it just converges *somewhere
else*. Converging to the same likelihood is what separates those two cases.

`burnin` is held at 25% of iterations, as in Phase 1.

#### The speed win, which is the point of the project

Median seconds per fit, single-threaded:

| Corpus | $K$ | Gibbs @200 | warp @1200 | speedup |
|---|---|---|---|---|
| small | 10 | 8.5 | 5.1 | 1.7× |
| small | 50 | 45.8 | 11.2 | 4.1× |
| medium | 10 | 143.9 | 70.6 | 2.0× |
| medium | 50 | 781.8 | 115.5 | **6.8×** |

**Whole 240-fit grid: 8.73 core-hours for Gibbs, 2.80 for warp — 3.1× overall,
at equal or better quality, with six times the iterations.** Per *iteration* the
gap is far larger (9× to 43×, growing in $K$ exactly as $O(NK)$ against
$O(VK+N)$ predicts); most of it is spent buying convergence back. Parallelism in
Phase 5 multiplies this further.

#### Three defects in the reference implementation

D2 says to port from text2vec. It should not be ported blindly.

**(a) `beta_bar` is a sum over the wrong index — the serious one.**
`LDA.hpp:65` sets `beta_bar = n_topic * beta`. The collapsed conditional needs
$\bar\eta_k = \sum_{v=1}^{V}\eta_{kv}$, a sum over the **vocabulary** — so
$V\eta$, not $K\eta$. (`alpha_bar` really is a sum over topics, which is
probably how the two got conflated.) On the medium corpus at $K = 50$ that is
222.15 against 2.5. Since $\bar\eta$ only ever appears added to $C_k$, the
difference substantially changes how far each acceptance ratio is shrunk toward
1. Porting it faithfully yields a valid chain converging to the wrong posterior
— which "it compiles and the likelihood curve looks right" would never catch.
tidylda's own Gibbs sampler has this right (`lda_gibbs2.cpp:333`).

**(b) The doc pass never maintains `C_doc`.** See the D10 caveat.

**(c) `C_local` is written but never sized.** `init()` (`LDA.hpp:69-90`) resizes
`C_doc`, `C_word`, `C_all` and `C_local_diff` but not `C_local`, which is then
written at `:125-126` and read at `:203`. Out-of-bounds write. The design notes
suspected this; it is real. Neither `C_local` nor `C_local_diff` is needed here
— both serve text2vec's distributed bookkeeping.

A fourth item is not a defect in the reference but a trap when adapting it:
`init()` fills `new_z` with a uniform random topic (`:77`). An acceptance ratio
is a valid MH step only if the proposal came from the matching proposal
distribution, and at iteration 1 no pass has run. Under warpLDA's own uniform
initialization this is harmless; under tidylda's informed initialization (D8) it
would discard part of what that initialization bought. The engine sets
`new_z = old_z`, so the first resolve is a no-op via the reference's own skip.

#### Notes for later phases

- `optimize_alpha` (D7) is accepted but ignored, with a once-per-session
  warning. Phase 6 turns that into a formal deprecation.
- Two tests are `skip()`ped pending Phase 3, both naming it in the skip message:
  vector $\boldsymbol\eta$ in `test-tidylda-fit-methods.R` and matrix
  $\boldsymbol\eta$ in `test-refit.tidylda.R`.
- The reference lives at `ignore/text2vec/` (pinned `0b31bdd8`, `.gitignore`d,
  not vendored per D2).

### 6.4 Phase 7 sketch — memory surgery for large corpora

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

D17 is a constant-factor improvement against this — it sparsifies two of the six
slots. It does not change the asymptotics, because `beta`, `lambda` and `theta`
stay dense.

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
path) both point this direction and should land first.

---

## 7. Open questions

1. **Work partitioning for parallelism.** The RNG scheme is settled (D12/D13),
   but how documents and words divide across threads is not — including whether
   the two passes should partition differently, and how the shared $C_k$ vector
   is updated without contention. Phase 5.
2. **Likelihood evaluation interval.** Default of 10 proposed. **Phase 1
   evidence supports it**; confirm against the warpLDA curve in Phase 2 before
   fixing the default. Essentially all readable structure is in the first ~50
   iterations (89–94% of total gain). Past iteration 150 the curve is
   noise-dominated: iteration-to-iteration movement is comparable to the
   remaining drift, with drift/jitter ratios of 1.6–7.1. Sampling every 10th
   iteration loses no real structure and yields a *cleaner* curve than every
   iteration does.
3. **Two expiring doc comments**, both user-facing via `man/tidylda.Rd`, both
   Phase 6:
   - `R/tidylda-fit-methods.R:58-66` claims the log likelihood returns "positive
     numbers, rather than the expected negative numbers." Measured pre-fix values
     were $-36530 \rightarrow -33551$, so this describes neither current nor
     pre-fix behavior. Most likely it once referred to row 3 of the likelihood
     matrix (`lpd + lp_alpha + lp_eta`), which genuinely can be positive since a
     Dirichlet log-*density* is unbounded above, and which `new_tidylda()` no
     longer reads.
   - `R/tidylda-fit-methods.R:68-69` states parallelism is not implemented.
     Expires when Phase 5 lands.

---

## 8. Environment and workflow

Things that cost time to discover once and should not cost it twice.

**Dependencies.** `RcppThread` (LinkingTo), `mvrsquared` and `tidytext`
(Imports) are required to build or load at all. `quanteda`, `tm`, `slam`,
`spelling` are Suggests exercised by the test suite — CI installs them, so
install them locally too or you verify less than CI does.

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
`devtools` and never hand-edited. `DESCRIPTION` stays at `0.0.7.999` until
release prep. NEWS already carries a `0.0.8` section from Phase 0.

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
