## Patch release version 0.0.5
This version is a patch. In this version I have

* Fixed "Packages in Suggests should be used conditionally" issue flagged by CRAN
    Used `testthat::skip_if_not_installed('tm')` in offending test.

## Test environments
* macOS (on GitHub actions): release
* ubuntu 22.04.3 (on GitHub actions): release, devel, and oldrel
* win-builder: release, devel, and oldrel

## R CMD check results
There are no NOTES, WARNINGs or ERRORs.

## revdepcheck results
There are currently no downstream dependencies for this package.