---
title: "Replacing tidylda's Gibbs Sampler with warpLDA"
subtitle: "Design notes for the `warp` branch"
author: "Tommy Jones (with Claude)"
date: "2026-08-20"
output:
  pdf_document:
    toc: true
    number_sections: true
---

# Purpose and Scope

These are working notes for replacing tidylda's collapsed Gibbs sampler with a
warpLDA-based engine written in Rcpp, parallelized with RcppThread, while
preserving tLDA — the matrix prior $\boldsymbol\eta$ that enables principled
transfer learning.

Nothing here is committed to code yet. The purpose is to establish:

1. Whether warpLDA's efficiency guarantees survive the tLDA matrix prior (they do);
2. Where $\boldsymbol\eta$ lands in warpLDA's memory access pattern (one place, and it is the convenient one);
3. What in tidylda's existing feature set constrains the design;
4. What is already broken and should be fixed before any of this begins.

Reference implementation studied: `text2vec/src/mcemlda/` (Dmitriy Selivanov, MIT
licensed), specifically `LDA.hpp`. Line numbers below refer to that file, or to
`tidylda/src/lda_gibbs2.cpp` where noted.

**Status of the defects in section 9:** fixed and pushed to `main` in commit
`5abaa96`, and merged into `warp`. The log likelihood reported by the current
sampler is now verified correct, giving a trustworthy baseline for the
statistical benchmarking in section 10.

---

# Background: The Two Samplers

## What tidylda does now

tidylda uses a standard collapsed Gibbs sampler (CGS). For each token, it
discounts the current assignment, computes the full conditional over all $K$
topics, samples, and re-increments (`lda_gibbs2.cpp:439-484`). The sampling
distribution is

$$p(z_{d,n} = k \mid \cdot) \;\propto\; \frac{C^{v}_{kv} + \eta_{kv}}{C_k + \bar\eta_k}\cdot\frac{C^{d}_{dk} + \alpha_k}{n_d + \bar\alpha}$$

where $C^{v}$ is the topic-word count matrix, $C^{d}$ the document-topic count
matrix, $C_k = \sum_v C^v_{kv}$, and

$$\bar\eta_k \equiv \sum_{v=1}^{V}\eta_{kv}, \qquad \bar\alpha \equiv \sum_{k=1}^{K}\alpha_k .$$

$\bar\eta_k$ is already computed in the current code as `sum_eta[k]`
(`lda_gibbs2.cpp:330-334`). The $(n_d + \bar\alpha)$ denominator is constant in
$k$ and does no work in the sampler — it costs $K$ divisions per token for
nothing.

**This is $O(NK)$ per iteration**, where $N$ is the total token count. Correct,
simple, and slow. Section 4 compares this against warpLDA's cost directly; the
gap is the whole motivation for this project.

## What warpLDA does

warpLDA replaces the exact conditional with alternating
Metropolis-Hastings proposals, and — this is the actual innovation — reorders
the computation so that every memory access falls inside a small, cache-resident
working set.

Two proposals:

- **Doc-proposal:** $q_d(k) \propto C^d_{dk} + \alpha_k$
- **Word-proposal:** $q_w(k) \propto C^v_{kv} + \eta_{kv}$

Two passes per iteration, over the *same* token array viewed through two
different orderings (CSR by document, CSC by word):

- **Doc pass:** tokens grouped by $d$. Working set is $C^d_{d\cdot}$, a $K$-vector.
- **Word pass:** tokens grouped by $v$. Working set is $C^v_{\cdot v}$, a $K$-vector.

Each token carries two topic labels, `old_z` and `new_z` (`LDA.hpp:30-35`). A
pass *generates* proposals of its own type into `new_z`, and *resolves* the
acceptances left pending by the other pass. Each pass rebuilds its own count
matrix from `old_z` as it goes — per document in the doc pass
(`LDA.hpp:151-157`), per word in the word pass (`LDA.hpp:98-104`). That
rebuild-as-you-go property matters in section 6.4.

---

# The Metropolis-Hastings Derivation Under a Matrix Prior

## The target

Under tLDA the prior $\boldsymbol\eta$ is a $K \times V$ matrix, each row $k$ a
distinct Dirichlet prior over the vocabulary. The collapsed conditional is
unchanged in form:

$$p(k) \;\propto\; \underbrace{(C^d_{dk} + \alpha_k)}_{\text{document factor}} \cdot \underbrace{\frac{C^v_{kv} + \eta_{kv}}{C_k + \bar\eta_k}}_{\text{word factor}}$$

The only structural change from standard LDA is that $\eta_{kv}$ carries a topic
subscript, and consequently $\bar\eta_k$ does too. **Critically, $\boldsymbol\eta$
is fixed for the entire sampling run** — it is constructed once from the $t-1$
posterior and never updated. So $\bar\eta_k$ is a $K$-vector computed once at
initialization.

## Doc-proposal acceptance

Propose $k'$ from $q_d(k) \propto C^d_{dk} + \alpha_k$. The document factor of
$p$ is exactly $q_d$, so it cancels:

$$\pi_d = \min\left(1,\; \frac{(C^v_{k'v} + \eta_{k'v})\,/\,(C_{k'} + \bar\eta_{k'})}{(C^v_{kv} + \eta_{kv})\,/\,(C_k + \bar\eta_k)}\right)$$

This depends only on **word-side quantities**, all at the same $v$. Confirmed at
`LDA.hpp:111-113` (symmetric-prior special case).

## Word-proposal acceptance

Propose $k'$ from $q_w(k) \propto C^v_{kv} + \eta_{kv}$. Writing it out in full:

$$\pi_w = \frac{p(k')}{p(k)}\cdot\frac{q_w(k)}{q_w(k')}
= \underbrace{\frac{(C^d_{dk'}+\alpha_{k'})}{(C^d_{dk}+\alpha_k)}\cdot\frac{(C^v_{k'v}+\eta_{k'v})}{(C^v_{kv}+\eta_{kv})}\cdot\frac{(C_k+\bar\eta_k)}{(C_{k'}+\bar\eta_{k'})}}_{p(k')/p(k)}\cdot\underbrace{\frac{(C^v_{kv}+\eta_{kv})}{(C^v_{k'v}+\eta_{k'v})}}_{q_w(k)/q_w(k')}$$

The $(C^v + \eta)$ factors cancel **exactly**, leaving

$$\pi_w = \min\left(1,\; \frac{C^d_{dk'}+\alpha_{k'}}{C^d_{dk}+\alpha_k}\cdot\frac{C_k+\bar\eta_k}{C_{k'}+\bar\eta_{k'}}\right)$$

which is exactly `LDA.hpp:164-166`. **No $v$-indexed quantity survives.**
$\eta_{kv}$ vanishes entirely; only the global $K$-vectors $C_k$ and
$\bar\eta_k$ remain.

**Nothing needs to be arranged for this to happen.** The cancellation is a free
consequence of choosing $q_w$ to match the word factor of $p$. There is no
implementation work here and no risk, provided $q_w$ is left alone.

The forward-looking warning: if $q_w$ is ever replaced with a cheaper
approximation that does *not* exactly match the word factor — for instance one
that ignores $\eta$ to avoid touching the matrix — the residual $\eta$ terms
stop cancelling and reappear inside $\pi_w$, which is evaluated in the **doc
pass**, where $v$ varies token to token. That converts a cache-resident lookup
into random access across a $K \times V$ matrix, defeating the entire design.
The exact cancellation is therefore a **constraint to preserve**, not a problem
to solve.

## Which quantities need $v$-indexed data

This is a separate question from the cancellation, and conflating the two caused
confusion in an earlier draft of these notes.

| Quantity | Evaluated in | Needs $v$-indexed data? |
|---|---|---|
| $\pi_w$ (word-proposal acceptance) | doc pass | **no** |
| $\pi_d$ (doc-proposal acceptance) | word pass | yes — $C^v_{\cdot v}$ and $\boldsymbol\eta_{\cdot v}$ |
| $q_w$ construction | word pass | yes — $C^v_{\cdot v}$ and $\boldsymbol\eta_{\cdot v}$ |

Memory layout matters only for the bottom two rows, both of which live in the
word pass and both of which need only the single column belonging to the word
currently being processed.

**$C^v$ is already correct in the reference.** `C_word.at(w, t)` is a `DenseMat`
with $w$ as the row, so the $K$ counts for a word are contiguous. Nothing to
change.

**The layout requirement falls entirely on $\boldsymbol\eta$**, which tidylda
currently stores as $K \times V$. The engine needs it as $V \times K$ so that
$\boldsymbol\eta_{\cdot v}$ is contiguous and sits beside the $C^v_{\cdot v}$
vector already in cache. Getting this wrong costs a cache miss per token;
getting it right makes the matrix prior nearly free.

## The doc-proposal draw and asymmetric $\boldsymbol\alpha$

The doc-proposal is drawn in $O(1)$ without an alias table
(`LDA.hpp:183-191`): with probability $L_d/(L_d + \bar\alpha)$ copy the topic of
a uniformly random token in the document, otherwise draw a topic uniformly at
random.

**That second branch assumes $\boldsymbol\alpha$ is symmetric.** A uniform draw
over $K$ is only proportional to $\alpha_k$ when all $\alpha_k$ are equal.
Asymmetric $\boldsymbol\alpha$ requires an alias table over $\boldsymbol\alpha$
in that branch.

**This costs nothing and is not worth further thought.** $\boldsymbol\alpha$ is a
$K$-vector, fixed for the whole run now that `optimize_alpha` is being dropped
(section 6.5). One alias table is built at initialization: $O(K)$ time, a few
kilobytes. Every draw afterward is $O(1)$ — identical asymptotics to the uniform
draw it replaces. An alias draw (one uniform, one comparison, one indexed load)
is plausibly *faster* than `rng.sample() % n_topic`, since integer modulo is a
division. At worst it is a wash.

It is a silent correctness trap only in the sense that omitting it produces code
that runs fine and samples from the wrong prior.

---

# Computational Cost

## The headline comparison

| Sampler | Cost per iteration |
|---|---|
| tidylda CGS (current) | $O(NK)$ — a $K$-length `qz` built for **every token** |
| text2vec warpLDA | $O(VK + N)$ — $O(K)$ per **word type**, then $O(1)$ per token |
| warpLDA paper (theoretical) | $O(N)$ |

Since $V \ll N$ in any real corpus, **the reference implementation is already an
enormous improvement over the current sampler.** Chasing the theoretical $O(N)$
is not required to get the win this project is after.

## The $O(VK)$ term, and why we accept it for now

The reference builds a full $K$-length proposal vector and alias table for every
word type, every iteration (`LDA.hpp:132-141`):

```cpp
AliasUrn<decltype(rng)> urn(rng);
vector<double> prob;
for(topic_index_t t=0; t<n_topic; t++)
    prob.push_back(C_word.at(w,t) + beta);
urn.setup(prob);
```

This is a dense loop over all $K$ topics per word — accepted as the price of the
cache-friendly layout. A useful side effect for us: since the loop is dense
regardless, **sparsity in $\boldsymbol\eta$ would buy nothing, so its absence
costs nothing.** The tLDA change to this loop is simply

```cpp
prob[t] = C_word.at(w,t) + eta.at(w,t);   // eta stored word-major
```

plus `beta_bar` $\rightarrow$ `eta_bar[t]` in both acceptance ratios. That is
substantially the entire sampler-side change.

The paper reaches $O(N)$ by exploiting that $C^v_{\cdot v}$ has at most
$\min(n_v, K)$ nonzeros, splitting $q_w$ into a sparse count part plus a dense
prior part served by a shared table. Under tLDA that dense part is
$\boldsymbol\eta_{\cdot v}$ — word-specific — so recovering $O(1)$ would mean
precomputing $V$ alias tables over $\boldsymbol\eta$ columns. They would be built
once, since $\boldsymbol\eta$ is fixed, but would cost roughly $2\times$ the size
of $\boldsymbol\eta$ in permanent memory (~400MB at $V=10^5, K=500$) on top of
$\boldsymbol\eta$ itself. That is a poor trade against the memory goals in
section 5.

When the $VK$ term actually bites — roughly when $VK \gtrsim N$:

| $V$ | $K$ | $N$ | $VK$ vs $N$ | Verdict |
|---|---|---|---|---|
| $5\times10^3$ | 10 | $5\times10^4$ | $\approx N$ | fine |
| $5\times10^4$ | 100 | $5\times10^6$ | $\approx N$ | ~2× off ideal, fine |
| $2\times10^5$ | 1000 | $10^7$ | $20\times N$ | would hurt |

**Decision: ship the $O(VK)$ approach. Revisit only if profiling at high $K$
shows it dominating.**

> **Implementation requirement.** The word-proposal construction must carry a
> comment recording that this is the $O(VK)$ formulation, that the $O(N)$
> alternative exists via the sparse/dense split with precomputed per-column
> alias tables, and what that alternative costs in memory. This optimization is
> wanted downstream; the comment is how we avoid re-deriving it.

## Allocation in the hot path

The reference's real inefficiency is not the dense loop but **allocation
frequency**. `prob` is declared *inside* the per-word loop and grown with
`push_back` to length $K$, and a fresh `AliasUrn` is constructed per word. That
is $V \times \text{iterations}$ heap allocations, each with several reallocations
as the vector grows.

The fix is to hoist, not to change container:

```cpp
std::vector<double> prob(n_topic);          // allocated once, outside the loop
for (w ...) {
  for (t ...) prob[t] = C_word.at(w,t) + eta.at(w,t);   // overwrite, no alloc
  urn.setup(prob);
}
```

`std::vector<double>` remains the right choice: $K$ is a runtime value, so
`std::array` is out; VLAs are not standard C++; `alloca` is a bad idea. Once
allocated, a `vector<double>` is contiguous memory with identical access codegen
and cache behavior to a raw array — the heap costs you only at allocation time,
and hoisting removes that from the hot path entirely.

**Under threading this becomes per-thread state.** Each thread needs its own
`prob` buffer and its own `AliasUrn`, or there is a data race.

## Other notes on the reference

`C_local` is used at `LDA.hpp:125-126` and `:203` but never `resize`d in
`init()` (only `C_all` and `C_local_diff` are, `LDA.hpp:72-73`). Either it is
sized elsewhere in `warplda.cpp` or this is latent undefined behavior. Do not
copy the pattern.

---

# On $\boldsymbol\eta$ Becoming Dense

## The tLDA recursion

From the tLDA vignette, with $a^{(t)}$ the single tuning parameter:

$$\boldsymbol\eta_k^{(t)} = \omega_k^{(t)}\cdot\mathbb{E}\left[\boldsymbol\beta_k^{(t-1)}\right], \qquad \omega_k^{(t)} = a^{(t)}\sum_{v=1}^{V}\left(C^{v,(t-1)}_{kv} + \eta^{(t-1)}_{kv}\right)$$

At $a^{(t)} = 1$ this reduces to
$\boldsymbol\eta^{(t)}_k = \boldsymbol{C}^{v,(t-1)}_k + \boldsymbol\eta^{(t-1)}_k$,
which unrolls to a geometrically-weighted sum of every previous period's count
matrix plus the scaled base prior.

## Dense from $t=1$, not gradually

In the implementation (`refit.tidylda.R:225`):

```r
eta$eta <- prior_weight * w_star * object$beta
```

$\hat{\boldsymbol\beta}^{(t-1)}$ is **dense at $t=1$** — every word has nonzero
posterior probability because $\eta^{(0)} > 0$ everywhere. So $\boldsymbol\eta$
is fully dense immediately after the first refit. There is no sparsity to
exploit at any $t$, and any sparse-plus-shared-base decomposition is dead on
arrival for tLDA.

## Why this is fine

Per section 4.2, the word-proposal construction is dense over $K$ regardless, so
dense $\boldsymbol\eta$ costs no additional asymptotic compute. Further,
`format_eta()` (`utils.R:92-137`) **already materializes a dense $K \times V$
matrix even when the user passes a scalar**, so dense storage is the status quo,
not a regression.

**Conclusion: carry $\boldsymbol\eta$ dense, column-major ($V \times K$), and
stop worrying about it.**

## Precision

**Store $\boldsymbol\eta$ as `float`, compute in `double`.** The $V \times K$
matrix is where the memory is; single precision halves it (at $K=500$,
$V=10^5$: 400MB $\rightarrow$ 200MB) and is amply sufficient for a prior.
Promote to `double` on read for the acceptance ratios — free, since the FPU
works in doubles anyway, and it keeps any precision drift out of the MH accept
decision.

Keep $C_k$ integer and $\bar\eta_k$ `double`. Both are $K$-vectors, trivially
small, and $C_k$ can exceed $2^{24}$, above which `float` stops representing
integers exactly.

## Remaining memory levers

1. **Scalar fast path.** Propagate `eta_class` through to C++ and avoid
   materializing $K \times V$ at all when the user passed a scalar or vector.
   A memory *win* over current behavior in the common non-transfer case.
2. **Thresholding.** When $a^{(t)} < 1$ (recency bias), contributions from old
   periods decay geometrically, so pruning small entries is statistically
   principled and yields a bounded effective memory horizon. Arguably a modeling
   knob worth exposing regardless. Lowest priority.

---

# Constraints From tidylda's Existing Features

## Call flow

All three entry points converge on the same two C++ calls:

```
tidylda() / refit.tidylda() / predict.tidylda()
  |- format_alpha() / format_eta()   [R]   alpha: K-vector; eta: ALWAYS dense K x V
  |- initialize_topic_counts()       [R]   builds beta_initial, theta_initial, Cd_start
  |    \- create_lexicon()           [C++] DTM -> token chain + informed 1-pass init
  |- fit_lda_c()                     [C++] the sampler
  \- new_tidylda()                   [R]   counts -> beta / theta / lambda / summary
```

The R/C++ round-tripping to be eliminated is one trip out and one back in:
`create_lexicon` returns `Docs`, `Zd`, `Cd`, `Cv`, `Ck` as R objects, which are
immediately handed back to `fit_lda_c`. Those five objects are the entire seam.

Fusing them is also a significant memory win: `Docs` and `Zd` are
`vector<vector<size_t>>` — **16 bytes per token** marshalled through R — where
the reference packs the same information into a 4-byte
`Z{uint16 old_z, uint16 new_z}` plus CSR/CSC index arrays.

## Informed initialization (retain)

`create_lexicon` (`lda_gibbs2.cpp:120-131`) samples each token's initial topic
from $P(z_{d,n} = k) \propto \hat\beta_{kv}\cdot\hat\theta_{dk}$, in log space.
warpLDA's `init()` (`LDA.hpp:75-78`) assigns topics **uniformly at random**.
That difference is the whole ballgame for seeded and tLDA models, so warpLDA's
initialization must be replaced wholesale rather than adapted.

Rationale for retaining: it is the more principled choice when $\boldsymbol\eta$
is asymmetric, it is arguably required for the transfer case, and a
uniform-random start is expected to increase time to convergence.

Two useful properties: it samples each token independently from a distribution
held *fixed* within the document (`Cd[d]` is not updated during the loop),
making it embarrassingly parallel; and it is the natural place to build the
CSR/CSC structures in the same sweep.

## `freeze_topics` — the prediction path (retain, specialize)

With topics frozen, $p(k) \propto \hat\beta_{kv}(C^d_{dk} + \alpha_k)$. In
warpLDA terms this is *simpler* than training:

- the word-proposal $q_w(k) \propto \hat\beta_{kv}$ is **fixed for the entire
  run**, so per-word alias tables are built once rather than per iteration;
- $\pi_d = \hat\beta_{k'v}/\hat\beta_{kv}$;
- $\boldsymbol\eta$, $C^v$, and $C_k$ are not used at all.

Best implemented as a separate specialization rather than a branch in the hot
loop. Note the current code marshals a full $K \times V$ `eta` and a zero-filled
`Cv`/`Ck` into `fit_lda_c` on the prediction path
(`predict.tidylda.R:261-277`) where all three are dead weight.

## `burnin` — posterior averaging (retain; **resolved, and free**)

`Cd_sum`/`Cv_sum` accumulate for `t >= burnin` and are divided by
`iterations - burnin` at the end (`lda_gibbs2.cpp:523-541`, `650-668`). This is
common in Bayesian practice and uncommon in LDA implementations; it is a
differentiating feature and must be preserved.

The apparent difficulty was that **no point in warpLDA's cycle has both $C^d$
and $C^v$ current simultaneously.** Each pass rebuilds its own matrix from
`old_z` as it goes, so at the end of the doc pass $C^d$ is consistent and $C^v$
is stale, and at the end of the word pass the reverse.

**Resolution: accumulate each matrix during its own pass.** Add to `Cd_sum`
during the doc pass and to `Cv_sum` during the word pass, each at the moment
that matrix is current. No reconciliation sweep is ever needed.

This is valid because the averaged entries estimate **marginal** expectations,
and marginals need not come from a synchronized joint sample — averaging
destroys joint coherence in any case. The totals also stay honest:
$\sum_{d,k} C^d = \sum_{v,k} C^v = N$ per pass, since each rebuild is exact for
its own slice.

The governing principle, which also settles where to evaluate the likelihood:
**average the same quantity you would report as a point estimate.** A run with
`burnin = -1` returns the counts at the end of the cycle; a run with burnin
returns their post-burnin mean; the likelihood describes that same object.

## `optimize_alpha` — **being removed**

`lda_gibbs2.cpp:629-634` rescales $\boldsymbol\alpha$ each iteration
proportional to $C_k$. This was a placeholder pending implementation of
fixed-point estimation from the literature, and carrying it forward would add
complexity in service of a hack.

**Decision: drop it.** Consequence: $\boldsymbol\alpha$ is fixed for the whole
run, so the alias table of section 3.5 is built exactly once.

## `refit` vocabulary alignment and topic addition (untouched)

`refit.tidylda.R:249-329` is pure R matrix surgery on `eta`, `beta_initial`,
`alpha`, and `theta_initial` before anything reaches C++: new vocabulary gets a
flat prior at the 10th percentile of $\boldsymbol\eta$; new topics get
$\text{colMeans}(\boldsymbol\eta)$ scaled to `additional_eta_sum`. This is the
fiddly tLDA bookkeeping and it stays exactly where it is.

## The `counts` contract (constrains engine output)

`posterior.tidylda()` reads `x$counts$Cd` and `x$counts$Cv` directly and uses
them as Dirichlet parameters (`posterior.tidylda.R:100-133`):

```r
dir_par <- x$counts$Cv[which, ] + eta$eta[which, ]   # indexed by TOPIC
```

So the model object's `Cv` is **topic-major ($K \times V$)**, and when
`burnin > -1` it holds *fractional* post-burnin means, which `rdirichlet`
accepts.

**Implication:** the engine works internally in word-major ($V \times K$) but
must transpose $C^v$ on output. One transpose at the end of the run — cheap, but
easy to forget, and it would silently corrupt `posterior()` rather than error.

---

# The Log Likelihood Problem

This is the one genuinely unresolved cost issue, and it is separate from burnin.

Unlike the count averages, the likelihood needs $C^d$ and $C^v$ **at the same
point in the chain** — $\theta$ from one, $\beta$ from the other. So it does
require a synchronized snapshot, meaning an $O(N)$ sweep over `old_z` to rebuild
whichever matrix is stale.

**That sweep is not the problem.** It is $O(N)$ with sequential token access,
roughly a sixth of an iteration. Negligible.

**The likelihood computation itself is the problem.** As currently written it is
$O(N K + VK)$:

- `beta_prob` over all topics and words: $O(VK)$
- the per-token inner sum over $K$ (`lda_gibbs2.cpp:602-613`): $O(NK)$

Against the current $O(NK)$ sampler that roughly doubles the cost. Against
warpLDA's $O(VK + N)$ sampler **it dominates by a factor of about $K$** — the
diagnostic would cost far more than the model it is diagnosing.

Two mitigations, which compose:

1. **Iterate the DTM's nonzeros, not tokens.** The inner sum is currently taken
   per token occurrence; taking it per unique $(d,v)$ pair gives
   $O(\text{nnz}(X)\cdot K)$ instead of $O(NK)$. For a corpus averaging three
   occurrences per document-term pair that is a free 3× saving.
2. **Evaluate the likelihood every $n$-th iteration**, exposed as a parameter
   (`likelihood_every` or similar) defaulting to around 10. This preserves
   `calc_likelihood = TRUE` as the default — worth keeping — while cutting the
   cost by that factor. A convergence curve sampled every 10 iterations is just
   as readable.

> **This is not thinning.** Thinning conventionally means discarding draws from
> the chain, which changes the posterior average. Nothing of the kind happens
> here. The chain advances every iteration and **every post-burnin iteration
> still contributes to `Cd_sum` and `Cv_sum` exactly as it does today.** The
> only thing made less frequent is evaluation of a read-only diagnostic that
> feeds nothing the sampler uses. The posterior average is untouched.

> **Test impact.** `test-tidylda-fit-methods.R:41` asserts
> `nrow(log_likelihood) == tail(iteration, 1) + 1`, which assumes every
> iteration is recorded. A non-unit evaluation interval breaks that assertion
> and the test must be updated alongside.

---

# RNG and Reproducibility Under Parallelism

## The governing constraint

CRAN's position is that parallelism is not an excuse for irreproducibility:
`set.seed()` must still work. This is a hard requirement, not a nice-to-have,
and it shapes the parallel design more than performance does.

Two rules follow immediately:

1. **R's RNG may only be touched from the main thread.** `R::unif_rand()` and
   anything else reaching into R's API is not thread-safe. Worker threads get
   their own generators. (Likewise, use `RcppThread::checkUserInterrupt()`, never
   `Rcpp::checkUserInterrupt()`, from workers — the current code already does
   this at `lda_gibbs2.cpp:428`.)
2. **All randomness must ultimately derive from R's stream**, seeded on the main
   thread at initialization, so that `set.seed()` controls the entire run.

## Decision: seed per work item, not per thread

The obvious scheme — seed once from R, then draw `threads` integers per
iteration to seed each thread's generator — satisfies CRAN but only delivers
reproducibility **at a fixed thread count**. Running `threads = 4` and
`threads = 8` from the same seed gives different answers, because the work
partition and therefore the stream assignment differ. This is common practice
and usually just documented.

**We take the stronger option: derive each stream from the work item's index.**

$$\text{seed} = f(\text{master},\; \text{iteration},\; \text{pass},\; \text{index})$$

where *index* is the document id in the doc pass and the word id in the word
pass. Then document $d$ consumes the same random numbers no matter which thread
happens to process it, and results are reproducible **independent of thread
count**.

Rationale:

- No caveat needed in the documentation, and nothing surprising for a user who
  changes `threads` between runs.
- It removes a confound from the benchmarking harness (section 10): we want to
  compare *samplers*, not thread counts.
- It is defensible to CRAN without qualification.

This works because the number of draws consumed within a document is
deterministic given that document's state — the proposal draws are fixed per
token, and the accept/reject branch (`LDA.hpp:109-110`) is itself determined by
the state. Each item's stream is therefore self-consistent.

## Warning: stream correlation

**Do not seed many streams with sequential integers.** xorshift-family
generators — including the `XOR128PLUS` in text2vec's `qrand.hpp` — produce
correlated streams from nearby seeds, which would quietly bias the sampler in a
way that looks like nothing at all until the benchmarks come back subtly wrong.
(The previous parallel attempt failed for a different reason — it used a hackier
approach lacking warpLDA's MH guarantees — but this is the same *kind* of
silent, benchmark-only-visible failure, which is why the benchmarking harness
exists.)

Two acceptable fixes:

- Expand each seed through `splitmix64` before initializing the generator, the
  standard companion to the xoshiro/xorshift family; or
- Use a counter-based generator designed for parallel streams (Threefry or
  Philox). The `sitmo` package provides Threefry for Rcpp and is built for
  precisely this pattern.

## The groundwork already exists

`src/sample.h:33-35`:

> The function below does not use a Random number generator. Instead, you take a
> sample from a uniform(0,1) RV and pass it in the `u` argument. This lets you
> use any RNG you'd like.

Both `sample_one` and `log_sample_one` take the variate as an argument
(`double &u`, `double &e`) rather than drawing it internally. The samplers are
already RNG-agnostic, which means the "never call R's RNG from a worker thread"
rule falls out **structurally** rather than requiring discipline at every call
site. Whatever generator we choose plugs in without touching the sampling logic.

---

# Defects Found

Recorded for provenance. Five were identified; **four were fixed** in `5abaa96`
on `main` and merged to `warp` — (a), (c), (d) and (e) below. Item (b) and the
two items under "Not fixed, by decision" were deliberately left alone because
the rewrite removes the code they live in.

**(a) Log-likelihood `denom` never reset per topic.** `lda_gibbs2.cpp:551`
declared `denom` outside the topic loop at `:555`, so it accumulated across
topics and `beta_prob` rows shrank monotonically in $k$. Verified by recomputing
the corpus log probability independently in R from the model's own returned
counts: the corrected value now agrees to a relative difference of $4\times
10^{-16}$, where the previous value disagreed by roughly 4245 nats. Fitted
`beta` and `theta` are bit-identical before and after under a fixed seed.

*Note:* this was **not** the cause of the "returning positive numbers" comment at
`tidylda-fit-methods.R:58-66`. The surfaced metric was already negative. That
complaint could only ever have applied to row 3 of the likelihood matrix
(`lpd + lp_alpha + lp_eta`), which `new_tidylda()` never reads. **That comment
appears stale and is still in the code.**

**(b) Batch loop broken for `threads > 1`.** `lda_gibbs2.cpp:426` initialized
from a document index and bounded by a batch size. Dead code from the removed
parallel implementation; the batching infrastructure goes away with the rewrite.
Not separately fixed.

**(c) Hardcoded `k = 10`** at `tidylda-fit-methods.R:241`. Inert, since `k` is
never read in `initialize_topic_counts()`. Fixed.

**(d) Positional argument slip** at `predict.tidylda.R:281-287`, binding
`threads` to `new_tidylda()`'s `alpha`. Inert via the `is_prediction` early
return. Fixed.

**(e) Unguarded Suggests.** `test-utils.R` called `cast_dfm()` (quanteda) while
guarding only on `tm`. Fixed with `skip_if_not_installed("quanteda")`.

**Not fixed, by decision:** `create_lexicon`'s $O(N_d N_v)$ dense scan over a
sparse matrix (`lda_gibbs2.cpp:106-108`, `:116-133`) — replaced wholesale by
CSR/CSC construction. And `ignore/nuclear_option/nuclear_option.cpp` carries the
same `denom` bug at lines 147/339/387/490 but is `.Rbuildignore`d and not
compiled.

---

# Statistical Benchmarking

**The question is equivalence, not model quality.** We are not asking whether
either sampler produces good topic models in the abstract; we are asking whether
moving from CGS to MH gives *comparable* results on the same model and the same
data. That framing simplifies the harness considerably — no held-out split, no
perplexity, no synthetic ground-truth recovery.

**Metrics**, both of which ship with tidylda and both computed from
`beta`/`theta`, so they are directly comparable across samplers at fixed data
and $K$:

- `calc_lda_r2()` (wraps `mvrsquared`) — global reconstruction of the DTM
- `calc_prob_coherence()` — topic quality, averaged across topics; top-5 terms
  by default

**Design: paired comparison across multiple seeds.** Both samplers are
stochastic, so a single seed can easily hide a systematically worse sampler.
Run both engines over the same corpora, the same $K$, and a set of seeds, then
**compare distributions rather than point estimates.**

The bar for merging is that the warpLDA engine's distributions of $R^2$ and mean
coherence are not detectably worse than the CGS engine's on the same inputs.

---

# Scope of the Replacement

The warpLDA engine subsumes **both** `create_lexicon()` and `fit_lda_c()`:
constructing the CSR/CSC token structure *is* what `create_lexicon` does, and
doing it inside the new engine is what collapses the R/C++ round trip.

Retained unchanged: `format_alpha`, `format_eta` (modulo the scalar fast path),
all of `refit`'s alignment logic, `new_tidylda`, `posterior`,
`tidy`/`augment`/`glance`.

New engine responsibilities:

1. DTM $\rightarrow$ dual-view token array (CSR by doc, CSC by word), built once, in C++;
2. informed initialization from $\hat\beta \cdot \hat\theta$, replacing warpLDA's uniform start;
3. alternating doc-pass / word-pass MH sampling with matrix $\boldsymbol\eta$;
4. a `freeze_topics` specialization for prediction;
5. per-pass post-burnin count accumulation (section 6.4);
6. log-likelihood computation, evaluated at intervals — *not* thinned (section 7);
7. transpose $C^v$ to topic-major on output.

## MH steps

**Configurable, default 1.** Default 1 reproduces the reference exactly and
costs nothing; having the parameter in place from the start is what allows
experimentation with mixing under tLDA's more informative priors, where a
concentrated $\boldsymbol\eta$ makes $q_w$ sharper and may change mixing
behavior.

Design cost to be aware of: warpLDA *defers* acceptance to the next pass, so
`mh_steps > 1` means storing that many pending proposals per token rather than a
single `new_z` — `mh_steps × 2` bytes per token. The inline propose-accept loop
used by LightLDA does not fit the deferred structure, so the array is the right
shape. At the default it is exactly today's `Z{old_z, new_z}`.

---

# Open Questions

1. **Work partitioning for parallelism.** The RNG scheme is settled (section 8),
   but how documents and words are divided across threads is not — including
   whether the two passes should partition differently, and how the shared
   $C_k$ vector is updated without contention.
2. **Likelihood evaluation interval.** Default of 10 proposed in section 7; the
   right value depends on how noisy the curve looks in practice.
3. **Two expiring doc comments** in `tidylda()`'s `@details`, both user-facing
   via `man/tidylda.Rd`:
   - `tidylda-fit-methods.R:58-66` claims the log likelihood is "returning
     positive numbers, rather than the expected negative numbers." Measured
     pre-fix values were $-36530 \rightarrow -33551$, so this describes neither
     current nor pre-fix behavior. Most likely it once referred to row 3 of the
     likelihood matrix (`lpd + lp_alpha + lp_eta`), which genuinely can be
     positive since a Dirichlet log-*density* is unbounded above, and which
     `new_tidylda()` no longer reads. Needs rewording regardless, since 0.0.8
     changed the reported values.
   - `tidylda-fit-methods.R:68-69` states "Parallelism, is not currently
     implemented. The `threads` argument is a placeholder for planned
     enhancements." Expires when this project lands.

   Worth doing both in one documentation pass, with `devtools::document()` to
   regenerate `man/tidylda.Rd`.

---

# Proposed Sequencing

1. ~~Fix defects on `main`; push; merge into `warp`.~~ **Done — `5abaa96`.**
2. Stand up the statistical benchmarking harness (section 10) against the
   corrected `main` as baseline.
3. Port warpLDA single-threaded with **scalar** prior, validated against the
   reference implementation and the harness. Any deviation at this stage is
   unambiguously an implementation bug rather than a tLDA subtlety.
4. Generalize to matrix $\boldsymbol\eta$; re-benchmark for statistical parity,
   not just speed.
5. Fuse initialization into the engine; eliminate the R round trip.
6. Add RcppThread parallelism **last**, as a separate change re-validated against
   the harness.

The ordering isolates the two historically risky things — statistical
correctness of the modified sampler, and parallel correctness — so that neither
arrives simultaneously with anything else.
