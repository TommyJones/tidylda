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
| **Last updated** | 2026-08-21 |
| **Branch** | `warp` |
| **Base commit** | `5abaa96` (Phase 0 fixes, merged from `main`) |
| **Current phase** | Phase 1 — benchmarking harness |
| **Last completed** | Phase 0: four defect fixes on `main`, pushed and merged. Design notes and this roadmap written. |
| **In flight** | Nothing. |

**Next action:** Build the statistical benchmarking harness to the specification
in §6.1, and run it against the *current* Gibbs sampler to establish and persist
the baseline distributions.

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
| D2 | Port from text2vec `src/mcemlda/` (MIT licensed) as reference | Compact, faithful warpLDA; confirms the derivation line by line | §1 |
| D3 | Store $\boldsymbol\eta$ **dense**, no sparse representation | $\hat\beta^{(t-1)}$ is dense at $t=1$, so tLDA has no sparsity to exploit at any $t$ | §5.2 |
| D4 | Store $\boldsymbol\eta$ **column-major ($V \times K$)** | The word pass needs $\boldsymbol\eta_{\cdot v}$ contiguous beside $C^v_{\cdot v}$; this is the transpose of tidylda's current $K \times V$ | §3.4 |
| D5 | $\boldsymbol\eta$ as `float`, promoted to `double` for computation | Halves the largest allocation; keeps precision drift out of MH accept decisions | §5.4 |
| D6 | $C_k$ stays integer; $\bar\eta_k$ stays `double` | $C_k$ can exceed $2^{24}$ where `float` loses integer exactness | §5.4 |
| D7 | **Drop `optimize_alpha`** | A placeholder for unimplemented fixed-point estimation; carrying it forward adds real complexity for a hack. Consequence: $\boldsymbol\alpha$ is fixed, so its alias table is built once | §6.5 |
| D8 | Keep informed $\hat\beta\cdot\hat\theta$ initialization; discard warpLDA's uniform start | Required for seeded and tLDA models; uniform init also expected to slow convergence | §6.2 |
| D9 | `freeze_topics` (prediction) as a **separate specialization**, not a branch in the hot loop | With topics frozen $q_w \propto \hat\beta_{kv}$ is fixed for the whole run, so alias tables build once; $\boldsymbol\eta$, $C^v$, $C_k$ go unused | §6.3 |
| D10 | Burnin: accumulate `Cd_sum` during the doc pass, `Cv_sum` during the word pass | Each pass rebuilds its own matrix from `old_z`, so each is current in its own pass. No reconciliation sweep needed. Valid because these estimate *marginal* expectations, and totals stay exact | §6.4 |
| D11 | Likelihood evaluated every $n$-th iteration, parameterized. **The scheme is settled; the default value is not** — see §7 open question 2 | The likelihood is $O(\text{nnz}\cdot K + VK)$ and would otherwise dominate an $O(VK+N)$ sampler by ~$K$. **This is not thinning** — the chain advances every iteration and every post-burnin iteration still contributes to the count sums. Only a read-only diagnostic runs less often | §7 |
| D12 | RNG seeded **per work item**, not per thread: `seed = f(master, iteration, pass, index)` | Gives reproducibility *independent of thread count*, not merely at a fixed count. Removes a confound from benchmarking and needs no caveat for CRAN. Master seed drawn from R's stream on the main thread so `set.seed()` governs | §8.2 |
| D13 | Expand seeds through `splitmix64`, or use a counter-based generator (Threefry/Philox, e.g. `sitmo`) | xorshift-family generators (including text2vec's `XOR128PLUS`) produce **correlated streams from nearby seeds** — a silent bias that looks like nothing until benchmarks come back subtly wrong | §8.3 |
| D14 | Benchmark for **equivalence**, not model quality: $R^2$ and mean probabilistic coherence, paired across multiple seeds, no held-out data | The question is whether MH matches Gibbs on the same model and data, not whether either is good in the abstract. Both metrics ship with tidylda. Pass/fail criterion in §6.1 | §10 |
| D15 | Accept the $O(VK)$ word-proposal construction for now | The $O(N)$ alternative needs $V$ precomputed alias tables over $\boldsymbol\eta$ columns, costing ~2× $\boldsymbol\eta$ in permanent memory. **The code must carry a comment recording this alternative** — it is wanted downstream | §4.2 |
| D16 | The new engine subsumes **both** `create_lexicon()` and `fit_lda_c()` | Building the CSR/CSC token structure *is* what `create_lexicon` does; fusing them is what eliminates the R↔C++ round trip and the 16-bytes-per-token marshalling | §11 |
| D17 | Transpose $C^v$ to topic-major on output | The engine works word-major internally, but `posterior.tidylda()` indexes `counts$Cv` by topic. Cheap, but forgetting it corrupts `posterior()` silently rather than erroring | §6.7 |
| D18 | MH steps configurable, default 1 | Default reproduces the reference exactly and costs nothing; the parameter is what allows experimentation with mixing under tLDA's sharper priors. Costs `mh_steps × 2` bytes per token above the default. Built in **Phase 2** | §11.1 |
| D19 | Alias table over $\boldsymbol\alpha$ in the doc-proposal draw — **binding**, Phase 2 | The reference's uniform-draw branch is only proportional to $\alpha_k$ when $\boldsymbol\alpha$ is symmetric; tidylda permits a vector. Omitting it yields code that runs fine and samples from the wrong prior. Costs one $O(K)$ setup, since D7 makes $\boldsymbol\alpha$ fixed | §3.5 |
| D20 | Scalar fast path for $\boldsymbol\eta$ — **deferred to Phase 4**, not Phase 2 | A memory win in the common non-transfer case, but it is an optimization, not a correctness requirement, and Phase 2 has enough moving parts. `format_eta()` keeps materializing $K \times V$ until then | §5.5 |

---

## 5. Invariants — must not break

| Invariant | Enforced by / at risk in |
|---|---|
| `counts$Cd` and `counts$Cv` are usable as Dirichlet parameters, indexed by document and **topic** respectively, and may hold fractional post-burnin means | `posterior.tidylda.R:100-133` |
| `set.seed()` reproducibility, including under parallelism — a CRAN requirement | D12, D13 |
| Public API unchanged: `tidylda()`, `refit()`, `predict()`, `posterior()`, `tidy`/`augment`/`glance` | all of `R/` |
| `refit`'s R-side vocabulary alignment and topic addition stay in R, untouched | `refit.tidylda.R:249-329` |
| Initialization remains informed, not uniform | D8 |
| R's RNG is touched only on the main thread; workers use `RcppThread::checkUserInterrupt()` | D12 |

---

## 6. Phase plan

| Phase | Deliverable | Exit criterion |
|---|---|---|
| **0** | Defect fixes on `main` | ✅ **Done — `5abaa96`** |
| **1** | Statistical benchmarking harness | Produces stable baseline distributions of $R^2$ and mean coherence for the current sampler, across multiple seeds and at least two corpora/$K$ settings |
| **2** | warpLDA engine, single-threaded, **scalar** prior. Includes D18 (`mh_steps`) and D19 (α alias table). **Known breakage to fix here:** `tests/testthat/test-tidylda-fit-methods.R:41` asserts `nrow(log_likelihood) == tail(iteration,1)+1`, which a non-unit likelihood interval (D11) invalidates | Matches CGS on the harness. Any deviation here is unambiguously an implementation bug, not a tLDA subtlety |
| **3** | Generalize to matrix $\boldsymbol\eta$ (tLDA) | Parity preserved when a transferred prior is in play; `refit()` path exercised |
| **4** | Fuse initialization into the engine; eliminate the R round trip | Identical results to Phase 3; one C++ entry point; per-token memory down from 16 bytes |
| **5** | RcppThread parallelism | Parity holds; `set.seed()` gives identical results across *different* thread counts (D12) |
| **6** | Documentation, NEWS, CRAN preparation | `devtools::check()` clean; the expiring comments in §7 rewritten; `man/` regenerated |

**Why this order.** The two historically risky things — statistical correctness
of the modified sampler, and parallel correctness — are isolated into separate
phases so neither arrives at the same time as anything else. Phase 2
deliberately uses a scalar prior so that a failure implicates the port and not
the matrix prior.

### 6.1 Phase 1 specification

Concrete enough to start from cold. This defines the comparison every later
phase is measured against, so the details are binding, not suggestions.

**Scope of Phase 1 specifically.** The warpLDA engine does not exist yet, so
Phase 1 runs the **Gibbs half only — 80 fits**. The 160-fit grid below describes
the full paired comparison that Phases 2, 3 and 5 execute against this stored
baseline.

**Files.** All under `warp-planning/benchmarks/` — tracked, and excluded from the
build by the `^warp-planning$` entry in `.Rbuildignore`. Not part of the package
and not run by `testthat`.

| File | Role |
|---|---|
| `build-corpora.R` | Builds and persists both DTMs. Run once |
| `run-benchmark.R` | Runs the grid for one sampler, writes an `.rds` |
| `compare.R` | Loads two `.rds` files, runs TOST, reports pass/fail |
| `data/` | Persisted DTMs |

**Corpora.** Two, at deliberately different scales:

| Name | Source | Size |
|---|---|---|
| `small` | `nih_sample_dtm` (ships with the package) | 100 docs |
| `medium` | fixed 1,000-doc sample of `nih` | 1,000 docs |

`nih` is a 68,508 × 44 data frame in `data/nih.rda`, present in the repo but
`.Rbuildignore`d (`^data/nih.rda`), so it is available for benchmarking and
never shipped. Relevant columns are `APPLICATION_ID` and `ABSTRACT_TEXT`. Build
the DTM with `tidytext::unnest_tokens()` + `cast_sparse()`, following
`tests/testthat/test-utils.R:7-11`.

> **The 1,000-document sample is drawn once and persisted**, never resampled.
> Use `set.seed(8675309)`, `dplyr::slice_sample(n = 1000)`, and write the
> resulting DTM to `warp-planning/benchmarks/data/nih-1000.rds`. If a later
> session resamples, its numbers are not comparable to this baseline and Phase 1
> has to be redone.

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

200/50 is the maintainer's standard setting for Gibbs on this package. Validate
it in Phase 1 by inspecting the likelihood curves; if they have not plateaued by
200, raise both numbers and rerun the baseline — before Phase 2 depends on it.

**Grid.** $K \in \{10, 50\}$, 20 seeds per cell. Two $K$ values probe whether
parity degrades as topic count grows, which is where the $O(VK)$ term and MH
mixing under sharp priors would both first show up.

**Metrics.** `calc_lda_r2()`, and the mean of `calc_prob_coherence()` across all
topics — a single number per fit. Both from `R/utils.R`.

**Pass/fail — the merge gate for Phases 2, 3 and 5.** Paired per-seed
differences $d_i = \text{warp}_i - \text{gibbs}_i$, with a **TOST equivalence
test** at $\alpha = 0.05$ and a margin of **5% of the Gibbs baseline mean**,
applied to each metric separately:

$$\text{PASS} \iff \text{both one-sided tests reject } H_0: |\bar d| \ge 0.05\,\overline{\text{gibbs}}$$

TOST rather than a plain paired t-test because the claim is *equivalence*; a
t-test that fails to reject means only that the study was underpowered.

Two notes on the margin, both deliberate:

- *Scale.* A relative margin is safe here because average coherence for a
  working model runs around 0.1 or above, and negative topic coherences cluster
  near zero rather than going far negative. If the **Gibbs** baseline mean
  coherence comes out near zero, treat that as a broken baseline — not a margin
  problem — and stop. If **warpLDA** lands near zero, the gate has correctly
  caught a loud failure.
- *Power.* Whether 20 seeds resolves a 5% margin is empirical and Phase 1 will
  reveal it. If the baseline spread makes the margin unreachable, **add seeds
  rather than widen the margin**, and record the change here.

**Persistence.** Write to
`warp-planning/benchmarks/baseline-5abaa96.rds`, named for the commit that
produced it, and fill in the table below so later phases can compare without
re-running and without trusting that a re-run reproduces.

### 6.2 Phase 1 baseline results

*Populate when Phase 1 completes. Until then this table is the outstanding
deliverable.*

| Corpus | $V$ | $K$ | mean $R^2$ (sd) | mean coherence (sd) | $R^2$ margin | coherence margin |
|---|---|---|---|---|---|---|
| small | — | 10 | — | — | — | — |
| small | — | 50 | — | — | — | — |
| medium | — | 10 | — | — | — | — |
| medium | — | 50 | — | — | — | — |

---

## 7. Open questions

1. **Work partitioning for parallelism.** The RNG scheme is settled (D12/D13),
   but how documents and words divide across threads is not — including whether
   the two passes should partition differently, and how the shared $C_k$ vector
   is updated without contention. Phase 5.
2. **Likelihood evaluation interval.** Default of 10 proposed; the right value
   depends on how noisy the curve looks in practice. Phase 1 or 2.
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
| `src/lda_gibbs2.cpp` | `create_lexicon()` (DTM → token chain + informed init) and `fit_lda_c()` (the sampler). Both subsumed by the new engine (D16) |
| `src/sample.h` | `sample_one()`, `log_sample_one()` — already RNG-agnostic, taking the variate as an argument. This is what makes D12 implementable without touching sampling logic |
| `src/matrix_conversions.h` | `mat_to_vec()` / `vec_to_mat()`; note `vec_to_mat` maps outer index → column |
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
