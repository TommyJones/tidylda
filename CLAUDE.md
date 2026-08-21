# tidylda

An R package implementing Latent Dirichlet Allocation with tidyverse
conventions, plus tLDA — transfer learning via a matrix prior over words in
topics. Fitting is done with a collapsed Gibbs sampler in Rcpp.

## The `warp` branch

`warp` is a long-running project replacing the collapsed Gibbs sampler with a
**warpLDA** engine (Rcpp + RcppThread), while preserving tLDA's matrix prior.
It spans many sessions and its decisions are recorded externally, not in
conversation history.

**Before doing any work on this branch, read both of these:**

1. `warp-planning/warplda-roadmap.md` — current status, the next concrete
   action, the decision log, and the invariants that must not break. **Start
   here.**
2. `warp-planning/warplda-design-notes.md` — the mathematics, cost analysis,
   and rationale behind those decisions.

Two things the roadmap asks of you: update its status ledger (§2) at the end of
a session, and treat its decision log (§4) as settled — if a decision genuinely
needs reopening, change it and record why in the same edit.

## Building and testing

Pandoc is not installed here, so vignette building fails and aborts a normal
check. Use:

```r
devtools::check("/home/tommy/tidylda", document = FALSE, vignettes = FALSE)
```

Expect one NOTE (a spelling diff flagging `ORCID` and `tidylda's`) plus a
local-only NOTE about `-mno-omit-leaf-frame-pointer`. Anything else is new.
See roadmap §8 for dependencies and CI details.
