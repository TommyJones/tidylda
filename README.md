
Latent Dirichlet Allocation Using 'tidyverse' Conventions

Copyright 2019 by Thomas W. Jones

As of this writing this package is in development and will not run on its own.

Implements an algorithim for Latent Dirichlet Allocation using style conventions from the [tidyverse](https://style.tidyverse.org/) and [tidymodels](https://tidymodels.github.io/model-implementation-principles/index.html). 
    
This implementation of LDA differs from others you may have used. Specifically it allows you to:

* use asymmetric prior parameters $\boldsymbol\alpha$ and $\boldsymbol\beta$
* use a previously-trained model as a prior for a new model
* apply LDA in a transfer-learning paradigm, updating a model's parameters with additional data (or additional iterations)

[//]: # (the [WarpLDA](https://arxiv.org/abs/1510.08628) algorithm for )

This may deviate from tidy modeling principles in a few ways:

* _Parallelism by default_ - the tidy modeling style guide says to default to sequential execution and let users decide to use parallelism. However, text data is large enough and compute intensive enough that I believe defaulting to sequential execution would greatly dimish the package's usability.