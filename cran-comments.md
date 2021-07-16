## Resubmission
Thank you to Gregor Seyer for a quick review. I've included excerpts of his
comments and how I have addressed them below.

CRAN wrote: "If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file [...]""

* I have added links to three works. The "novel" aspects of this package are
  based on my own research which is ongoing and as-yet unpublished. As that
  research is published, I will add references. Note that this package does
  allow for a canonical use of LDA. The novel features are labeled as
  experimental in the README. The documentation of refit.tidylda includes
  implementation details of this experimental research.
  
CRAN wrote: "Please add \\value to .Rd files regarding exported methods and explain
the functions results in the documentation. Please write about the
structure of the output (class) and also what the output means. [...]"

* I have added the value tag to the following documentation, as requested by
  CRAN, describing the respective function's output
  - calc_lda_r2.Rd: \value
  - create_lexicon.Rd: \value
  - fit_lda_c.Rd: \value
  - print.tidylda.Rd: \value
  - recover_counts_from_probs.Rd: \value
  - tidy_triplet.Rd: \value
  - tidylda_bridge.Rd: \value
* I would like to note that all of these functions are internal. They are not
  exported in tidylda's NAMESPACE. So, I beg CRAN's leniency on the
  formality of the documentation.
  
CRAN wrote: "\\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \\dontrun{} adds the comment
("# Not run:") as a warning for the user. Does not seem necessary. [...]"

* I have modified the examples wrapped in \\dontrun{} so that they are able to
  be run and replaced \\dontrun{} with \\donttest{}

## Test environments
* local macOS install: release
* macOS (on GitHub actions): release
* ubuntu 20.04 (on GitHub actions): release
* win-builder: release, devel, and oldrel

## R CMD check results
There are no WARNINGs or ERRORs.

There is one NOTE about possibly misspelled words in the description. 

* These words {Blei, Wickham, al, et} are definitely spelled correctly. :)

## revdepcheck results
There are currently no downstream dependencies for this package.