## Patch release version 0.0.3
This version is a patch. In this version I have

* Updated internal function `convert_dtm` to not use functions deprecated as of
  Matrix 1.4-2
* Updates for compatibility with R CMD check and tidy select variables

## Test environments
* macOS (on GitHub actions): release
* ubuntu 20.04 (on GitHub actions): release, devel, and oldrel
* win-builder: release, devel, and oldrel

## R CMD check results
There are no NOTES, WARNINGs or ERRORs.

## revdepcheck results
There are currently no downstream dependencies for this package.