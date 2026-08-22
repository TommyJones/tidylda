# tidylda 0.0.8
* Fixed the calculation of `theta` when `alpha` is asymmetric. The prior was
    added along the wrong axis of the document-topic count matrix, so instead of
    `alpha[k]` being added to topic `k`, the values were recycled diagonally
    across topics. Models fitted with a scalar `alpha` are unaffected, since
    every entry is then equal; models fitted with a vector `alpha` — including
    output from `refit.tidylda()` and any model fitted with
    `optimize_alpha = TRUE` — will now report different, correct values of
    `theta`.
* Fixed a miscalculation in the log likelihood reported by `tidylda()` and
    `refit.tidylda()`. A normalizing denominator was accumulated across topics
    instead of being reset for each topic, so the topic-word probabilities used
    in the calculation were incorrectly scaled. Note that this changes the log
    likelihood values reported for a given model relative to previous versions.
    Model fitting itself is unaffected.
* Removed a hardcoded number of topics in an internal call from `tidylda()` to
    `initialize_topic_counts()`. The value was unused, so model results are
    unaffected.
* Fixed an internal call in `predict.tidylda()` that passed an argument
    positionally into the wrong parameter of `new_tidylda()`. The argument was
    unused on the prediction path, so predictions are unaffected.
* Tests using the 'quanteda' package are now skipped when it is not installed,
    consistent with its status as a suggested package.

# tidylda 0.0.7
* Added additional checks to refit.tidylda()
* Fixed WARNINGS related to using a deprecated function from the C++ Armadillo
    library, as noted by CRAN checks on tidylda's page

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