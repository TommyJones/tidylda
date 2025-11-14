## Patch version 0.0.7
This version is a patch. In this version I have

* Added additional checks to refit.tidylda()
* Fixed WARNINGS related to using a deprecated function from the C++ Armadillo
    library, as noted by CRAN checks on tidylda's page

## Test environments
* local macOS install: release
* macOS (on GitHub actions): release
* ubuntu 24.04.3 (on GitHub actions): release
* win-builder: release, devel, and oldrel

## R CMD check results
One NOTE from win-builder running R version 4.4.3 (old release). 

> NOTE
  Author field differs from that derived from Authors@R
  Author:    'Tommy Jones [aut, cre] (ORCID: <https://orcid.org/0000-0001-6457-2452>), Brendan Knapp [ctb] (ORCID: <https://orcid.org/0000-0003-3284-4972>), Barum Park [ctb]'
  Authors@R: 'Tommy Jones [aut, cre] (<https://orcid.org/0000-0001-6457-2452>), Brendan Knapp [ctb] (<https://orcid.org/0000-0003-3284-4972>), Barum Park [ctb]'

This has not appeared on any other flavor of R CMD check.

## revdepcheck results
There are currently no downstream dependencies for this package.