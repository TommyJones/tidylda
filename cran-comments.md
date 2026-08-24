## Minor version 0.1.0

This release replaces the package's sampler. Model fitting now uses warpLDA
(Chen et al., 2016, <doi:10.48550/arXiv.1510.08628>), a Metropolis-Hastings
scheme, in place of the collapsed Gibbs sampler used through 0.0.7. Fit quality
was validated against the previous sampler across a grid of corpus sizes and
topic counts, under both scalar and matrix priors, before the change was made.

User-visible changes:

* One breaking change: the `counts$Cv` element of a fitted model is now tokens
    by topics rather than topics by tokens, and both count matrices are now
    sparse. This is documented in NEWS.md with the one-line fix. Models saved by
    earlier versions are detected and read correctly.
* Two deprecations, both backward compatible and warning once per session:
    `predict(method = "gibbs")` is now `method = "mh"`, and `optimize_alpha` is
    accepted but ignored.
* The parallelism promised by the `threads` argument is now implemented, with
    results independent of the number of threads.

## Test environments

* local: Ubuntu 24.04, R 4.6.0
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release,
    oldrel-1)
* win-builder: devel, release, oldrelease
* macOS builder: release

## R CMD check results

0 errors | 0 warnings | 0 notes

The local Ubuntu check reports one additional NOTE that does not appear on any
other platform:

> checking compilation flags used ... NOTE
>   Compilation used the following non-portable flag(s):
>     '-mno-omit-leaf-frame-pointer'

This flag is not set by the package. It comes from the Debian/Ubuntu build of R
itself, where it appears in `/usr/lib/R/etc/Makeconf` and is applied to every
package compiled on that system. `src/Makevars` sets only
`$(SHLIB_OPENMP_CXXFLAGS)` and `-DARMA_64BIT_WORD=1`.

## revdepcheck results

There are currently no downstream dependencies for this package, confirmed
against the CRAN package database with `tools::package_dependencies()`.
