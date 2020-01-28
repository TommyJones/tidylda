<!-- badges: start -->
  [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/TommyJones/tidylda?branch=master&svg=true)](https://ci.appveyor.com/project/TommyJones/tidylda)
[![Travis-CI Build Status](https://travis-ci.com/TommyJones/tidylda.svg?branch=master)](https://travis-ci.com/TommyJones/tidylda)
[![Coverage Status](https://img.shields.io/codecov/c/github/tommyjones/tidylda/master.svg)](https://codecov.io/github/tommyjones/tidylda?branch=master)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

Latent Dirichlet Allocation Using 'tidyverse' Conventions

Copyright 2020 by Thomas W. Jones

Implements an algorithim for Latent Dirichlet Allocation using style conventions from the [tidyverse](https://style.tidyverse.org/) and [tidymodels](https://tidymodels.github.io/model-implementation-principles/). 
    
In addition this implementation of LDA allows you to:

* use asymmetric prior parameters alpha and beta
* use a matrix prior parameter, beta to seed topics into a model
* use a previously-trained model as a prior for a new model
* apply LDA in a transfer-learning paradigm, updating a model's parameters with additional data (or additional iterations)

Note that the seeding of topics and transfer learning are **experimental** for now. They are almost-surely useful but their behaviors have not been optimized or well-studied. Caveat emptor!

