# tidylda 0.0.7
* Additional checks added to refit.tidylda()

# tidylda 0.0.6
* Lifecycle is now stable. Removed references to experimental lifecycle.
* Minor updates to README and vignettes


# tidylda 0.0.5
* Fixed "Packages in Suggests should be used conditionally" issue flagged by CRAN
    Used `testthat::skip_if_not_installed('tm')` in offending test.
* Updated documentation to have valid NIH URL

# tidylda 0.0.4
* Fixed an issue flagged by CRAN related to RcppExports.cpp.
    See [here](https://github.com/RcppCore/Rcpp/issues/1287) for more info.

# tidylda 0.0.3

* Added "class" and "distribution" options for `predict.tidylda` outputs
* Updated internal function `convert_dtm` to not use functions deprecated as of
  Matrix 1.4-2
* Updates for compatibility with R CMD check and tidy select variables
* Added vignettes to describe some of the novel features of tidylda.
* Fix a bug in `tidylda` where data not returned even if user specifies `return_data = TRUE`
* Patch a potential error caused in internal function `tidylda:::recover_counts_from_probs`
* Updated C++11 requirement consistent with current CRAN compilers

# tidylda 0.0.2

* Fixed error encountered with call to `tidylda` with large data sets.
* Improved user experience for using `refit.tidylda` when fine tuning on only
  one document.
* Fixed miscalculation in `refit.tidylda` when beta from a previous model is used
  as the prior. Miscalculation only affected multiple sequential calls to `refit`
* Minor improvements to documentation.
* Model summary now displays top 5 terms per topic, instead of top 3.
* Removed all explicit dependencies on the `textmineR` package.

# tidylda 0.0.1
This is the first released version of tidylda!