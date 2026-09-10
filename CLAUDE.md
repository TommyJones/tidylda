# tidylda

An R package implementing Latent Dirichlet Allocation with tidyverse
conventions, plus tLDA — transfer learning via a matrix prior over words in
topics. Fitting is done with a warpLDA Metropolis-Hastings sampler in Rcpp.

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

Pandoc ships with RStudio Server here but is not on the default `PATH`, so
vignette building fails and aborts a normal check unless you point at it:
`RSTUDIO_PANDOC=/usr/lib/rstudio-server/bin/quarto/bin/tools/x86_64`. For
routine checks, skip vignettes:

```r
devtools::check("/home/tommy/tidylda", document = FALSE, vignettes = FALSE)
```

Expect a local-only NOTE about `-mno-omit-leaf-frame-pointer` (it comes from
Ubuntu's R build flags, not from us) and a local-only WARNING that a complete
check needs `checkbashisms` — that script is not installed here, and CRAN's
Debian machines do run it against `configure` and `cleanup`. To check those
yourself, fetch `checkbashisms.pl` from devscripts and run
`perl checkbashisms -f configure cleanup`. Anything else is new.
See roadmap §8 for dependencies and CI details.
