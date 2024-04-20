## Patch release version 0.0.5
Resubmission to fix one other instance of an invalid URL. All is good now.
Sorry for missing it the first time.

This version is a patch. In this version I have

* Fixed "Packages in Suggests should be used conditionally" issue flagged by CRAN
    Used `testthat::skip_if_not_installed('tm')` in offending test.
* Updated documentation to have valid NIH URL

## Test environments
* macOS (on GitHub actions): release
* ubuntu 22.04.3 (on GitHub actions): release, devel, and oldrel
* win-builder: release, devel, and oldrel

## R CMD check results
There are no NOTES, WARNINGs or ERRORs.

## revdepcheck results
There are currently no downstream dependencies for this package.