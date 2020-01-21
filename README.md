
Note: As of this writing this package is in development and will not run on its own.
I'm going to fold some of this into the next update for [textmineR](https://github.com/tommyjones/textminer) before concentrating on making this its own package.

Latent Dirichlet Allocation Using 'tidyverse' Conventions

Copyright 2019 by Thomas W. Jones

Implements an algorithim for Latent Dirichlet Allocation using style conventions from the [tidyverse](https://style.tidyverse.org/) and [tidymodels](https://tidymodels.github.io/model-implementation-principles/). 
    
This implementation of LDA allows you to:

* use asymmetric prior parameters $\boldsymbol\alpha$ and $\boldsymbol\beta$
* use a matrix prior parameter, $\boldsymbol\beta$ to seed topics into a model
* use a previously-trained model as a prior for a new model
* apply LDA in a transfer-learning paradigm, updating a model's parameters with additional data (or additional iterations)


