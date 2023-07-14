## Patch release version 0.0.3
This version is a patch. In this version I have

* Added "class" and "distribution" options for `predict.tidylda` outputs
* Updated internal function `convert_dtm` to not use functions deprecated as of
  Matrix 1.4-2
* Updates for compatibility with R CMD check and tidy select variables
* Added vignettes to describe some of the novel features of tidylda.
* Fix a bug in `tidylda` where data not returned even if user specifies `return_data = TRUE`
* Patch a potential error caused in internal function `tidylda:::recover_counts_from_probs`
* Updated C++11 requirement consistent with current CRAN compilers

## Test environments
* macOS (on GitHub actions): release
* ubuntu 22.04.2 (on GitHub actions): release, devel, and oldrel
* win-builder: release, devel, and oldrel

## R CMD check results
There are no NOTES, WARNINGs or ERRORs.

## revdepcheck results
There are currently no downstream dependencies for this package.