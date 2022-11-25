# tidylda 0.0.3

* Updated internal function `convert_dtm` to not use functions deprecated as of
  Matrix 1.4-2

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