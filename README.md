
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tidylda <img src='man/figures/logo.png' align="right" height="136.5" />

<!-- badges: start -->

[![DOI](https://joss.theoj.org/papers/10.21105/joss.06800/status.svg)](https://doi.org/10.21105/joss.06800)
[![R-CMD-check](https://GitHub.com/TommyJones/tidylda/actions/workflows/R-CMD-check.yaml/badge.svg)](https://GitHub.com/TommyJones/tidylda/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Codecov test
coverage](https://codecov.io/gh/TommyJones/tidylda/graph/badge.svg)](https://app.codecov.io/gh/TommyJones/tidylda)
<!-- badges: end -->

Latent Dirichlet Allocation Using ‘tidyverse’ Conventions

`tidylda` implements an algorithm for Latent Dirichlet Allocation using
style conventions from the [tidyverse](https://style.tidyverse.org/) and
[tidymodels](https://tidymodels.GitHub.io/model-implementation-principles/).

In addition this implementation of LDA allows you to:

- use asymmetric prior parameters alpha and eta
- use a matrix prior parameter, eta, to seed topics into a model
- use a previously-trained model as a prior for a new model
- apply LDA in a transfer-learning paradigm, updating a model’s
  parameters with additional data (or additional iterations)

Fitting uses [warpLDA](https://arxiv.org/abs/1510.08628) (Chen et al.,
2016), a Metropolis-Hastings sampler that alternates document-ordered
and word-ordered passes so each pass touches only a small,
cache-resident working set. It replaced the collapsed Gibbs sampler in
version 0.1.0 and is multithreaded, with results that do not depend on
the number of threads.

## Installation

You can install the latest CRAN release with:

``` r
install("tidylda")
```

You can install the development version from
[GitHub](https://GitHub.com/) with:

``` r
install.packages("remotes")

remotes::install_GitHub("tommyjones/tidylda")
```

For a list of dependencies see the DESCRIPTION file.

# Getting started

This package is still in its early stages of development. However, some
basic functionality is below. Here, we will use the `tidytext` package
to create a document term matrix, fit a topic model, predict topics of
unseen documents, and update the model with those new documents.

`tidylda` uses the following naming conventions for topic models:

- `theta` is a matrix whose rows are distributions of topics over
  documents, or P(topic\|document)
- `beta` is a matrix whose rows are distributions of tokens over topics,
  or P(token\|topic)
- `lambda` is a matrix whose rows are distributions of topics over
  tokens, or P(topic\|token) `lambda` is useful for making predictions
  with a computationally-simple and efficient dot product and it may be
  interesting to analyze in its own right.
- `alpha` is the prior that tunes `theta`
- `eta` is the prior that tunes `beta`

## Example

``` r
library(tidytext)
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2)
library(tidyr)
library(tidylda)
library(Matrix)
#> 
#> Attaching package: 'Matrix'
#> The following objects are masked from 'package:tidyr':
#> 
#>     expand, pack, unpack

### Initial set up ---
# load some documents
docs <- nih_sample 

# tokenize using tidytext's unnest_tokens
tidy_docs <- docs %>% 
  select(APPLICATION_ID, ABSTRACT_TEXT) %>% 
  unnest_tokens(output = word, 
                input = ABSTRACT_TEXT,
                stopwords = stop_words$word,
                token = "ngrams",
                n_min = 1, n = 2) %>% 
  count(APPLICATION_ID, word) %>% 
  filter(n>1) #Filtering for words/bigrams per document, rather than per corpus

tidy_docs <- tidy_docs %>% # filter words that are just numbers
  filter(! stringr::str_detect(tidy_docs$word, "^[0-9]+$"))

# append observation level data 
colnames(tidy_docs)[1:2] <- c("document", "term")


# turn a tidy tbl into a sparse dgCMatrix 
# note tidylda has support for several document term matrix formats
d <- tidy_docs %>% 
  cast_sparse(document, term, n)

# let's split the documents into two groups to demonstrate predictions and updates
d1 <- d[1:50, ]

d2 <- d[51:nrow(d), ]

# make sure we have different vocabulary for each data set to simulate the "real world"
# where you get new tokens coming in over time
d1 <- d1[, colSums(d1) > 0]

d2 <- d2[, colSums(d2) > 0]

### fit an intial model and inspect it ----
set.seed(123)

lda <- tidylda(
  data = d1,
  k = 10,
  iterations = 200, 
  burnin = 175,
  alpha = 0.1, # also accepts vector inputs
  eta = 0.05, # also accepts vector or matrix inputs
  optimize_alpha = FALSE, # experimental
  calc_likelihood = TRUE,
  calc_r2 = TRUE, # see https://arxiv.org/abs/1911.11061
  return_data = FALSE
)

# did the model converge?
# there are actual test stats for this, but should look like "yes"
qplot(x = iteration, y = log_likelihood, data = lda$log_likelihood, geom = "line") + 
    ggtitle("Checking model convergence")
#> Warning: `qplot()` was deprecated in ggplot2 3.4.0.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
```

<img src="man/figures/README-example-1.png" alt="" width="100%" />

``` r

# look at the model overall
glance(lda)
#> # A tibble: 1 × 5
#>   num_topics num_documents num_tokens iterations burnin
#>        <int>         <int>      <int>      <dbl>  <dbl>
#> 1         10            50       1524        200    175

print(lda)
#> A Latent Dirichlet Allocation Model of  10 topics,  50  documents, and  1524  tokens:
#> tidylda(data = d1, k = 10, iterations = 200, burnin = 175, alpha = 0.1, 
#>     eta = 0.05, optimize_alpha = FALSE, calc_likelihood = TRUE, 
#>     calc_r2 = TRUE, return_data = FALSE)
#> 
#> The model's R-squared is  0.2374 
#> The  5  most prevalent topics are:
#> # A tibble: 10 × 4
#>   topic prevalence coherence top_terms                                          
#>   <dbl>      <dbl>     <dbl> <chr>                                              
#> 1     7       12.7     0.198 cancer, dcis, breast, imaging, clinical, ...       
#> 2     8       12.0     0.276 mitochondrial, function, cmybp, sleep, studies, ...
#> 3     5       11.7     0.222 research, cancer, disparities, program, mhirt, ... 
#> 4     1       11.4     0.124 effects, imaging, v4, stiffening, wall, ...        
#> 5     6       10.3     0.14  risk, diabetes, sud, numeracy, redox, ...          
#> # ℹ 5 more rows
#> 
#> The  5  most coherent topics are:
#> # A tibble: 10 × 4
#>   topic prevalence coherence top_terms                                          
#>   <dbl>      <dbl>     <dbl> <chr>                                              
#> 1     4       8.17     0.466 cdk5, nmdar, nmdars, adhesion, calpain, ...        
#> 2     8      12.0      0.276 mitochondrial, function, cmybp, sleep, studies, ...
#> 3     5      11.7      0.222 research, cancer, disparities, program, mhirt, ... 
#> 4     7      12.7      0.198 cancer, dcis, breast, imaging, clinical, ...       
#> 5     3       8.43     0.142 cells, cell, specific, infection, lung, ...        
#> # ℹ 5 more rows

# it comes with its own summary matrix that's printed out with print(), above
lda$summary
#> # A tibble: 10 × 4
#>    topic prevalence coherence top_terms                                         
#>    <dbl>      <dbl>     <dbl> <chr>                                             
#>  1     1      11.4     0.124  effects, imaging, v4, stiffening, wall, ...       
#>  2     2       8.62    0.111  cns, research, extinction, lsds, drug, ...        
#>  3     3       8.43    0.142  cells, cell, specific, infection, lung, ...       
#>  4     4       8.17    0.466  cdk5, nmdar, nmdars, adhesion, calpain, ...       
#>  5     5      11.7     0.222  research, cancer, disparities, program, mhirt, ...
#>  6     6      10.3     0.14   risk, diabetes, sud, numeracy, redox, ...         
#>  7     7      12.7     0.198  cancer, dcis, breast, imaging, clinical, ...      
#>  8     8      12.0     0.276  mitochondrial, function, cmybp, sleep, studies, .…
#>  9     9       7.34    0.022  disease, data, ppg, pd, behavior, ...             
#> 10    10       9.36    0.0913 brain, plasticity, mechanisms, study, repair, ...


# inspect the individual matrices
tidy_theta <- tidy(lda, matrix = "theta")

tidy_theta
#> # A tibble: 500 × 3
#>    document topic   theta
#>    <chr>    <dbl>   <dbl>
#>  1 8574224      1 0.00238
#>  2 8574224      2 0.00238
#>  3 8574224      3 0.00238
#>  4 8574224      4 0.00238
#>  5 8574224      5 0.00238
#>  6 8574224      6 0.148  
#>  7 8574224      7 0.00238
#>  8 8574224      8 0.00238
#>  9 8574224      9 0.00238
#> 10 8574224     10 0.833  
#> # ℹ 490 more rows

tidy_beta <- tidy(lda, matrix = "beta")

tidy_beta
#> # A tibble: 15,240 × 3
#>    topic token             beta
#>    <dbl> <chr>            <dbl>
#>  1     1 adolescence  0.0000652
#>  2     1 age          0.0000652
#>  3     1 application  0.0000652
#>  4     1 depressive   0.0000652
#>  5     1 disorder     0.0000652
#>  6     1 emotionality 0.0000652
#>  7     1 information  0.00246  
#>  8     1 mdd          0.0000652
#>  9     1 onset        0.0000652
#> 10     1 onset mdd    0.0000652
#> # ℹ 15,230 more rows

tidy_lambda <- tidy(lda, matrix = "lambda")

tidy_lambda
#> # A tibble: 15,240 × 3
#>    topic token         lambda
#>    <dbl> <chr>          <dbl>
#>  1     1 adolescence  0.00777
#>  2     1 age          0.00929
#>  3     1 application  0.00784
#>  4     1 depressive   0.0204 
#>  5     1 disorder     0.0204 
#>  6     1 emotionality 0.0204 
#>  7     1 information  0.254  
#>  8     1 mdd          0.0114 
#>  9     1 onset        0.00789
#> 10     1 onset mdd    0.0204 
#> # ℹ 15,230 more rows

# append observation-level data
augmented_docs <- augment(lda, data = tidy_docs)
#> Joining with `by = join_by(document, term, n)`

augmented_docs
#> # A tibble: 4,566 × 4
#>    document term            n topic
#>    <chr>    <chr>       <int> <int>
#>  1 8574224  adolescence     1     6
#>  2 8646901  adolescence     1     6
#>  3 8689019  adolescence     1     6
#>  4 8705323  adolescence     1     6
#>  5 8574224  age             1    10
#>  6 8705323  age             1    10
#>  7 8757072  age             1    10
#>  8 8823186  age             1    10
#>  9 8574224  application     1     6
#> 10 8605875  application     1     6
#> # ℹ 4,556 more rows

### predictions on held out data ---
# two methods: mh (Metropolis-Hastings) is cleaner and more technically
# correct in the bayesian sense
p_mh <- predict(lda, new_data = d2[1, ], iterations = 100, burnin = 75)

# dot is faster, less prone to error (e.g. underflow), noisier, and frequentist
p_dot <- predict(lda, new_data = d2[1, ], method = "dot")

# pull both together into a plot to compare
tibble(topic = 1:ncol(p_mh), mh = p_mh[1,], dot = p_dot[1, ]) %>%
  pivot_longer(cols = mh:dot, names_to = "type") %>%
  ggplot() + 
  geom_bar(mapping = aes(x = topic, y = value, group = type, fill = type), 
           stat = "identity", position="dodge") +
  scale_x_continuous(breaks = 1:10, labels = 1:10) + 
  ggtitle("Metropolis-Hastings predictions vs. dot product predictions")
```

<img src="man/figures/README-example-2.png" alt="" width="100%" />

``` r

### Augment as an implicit prediction using the 'dot' method ----
# Aggregating over terms results in a distribution of topics over documents
# roughly equivalent to using the "dot" method of predictions.
augment_predict <- 
  augment(lda, tidy_docs, "prob") %>%
  group_by(document) %>% 
  select(-c(document, term)) %>% 
  summarise_all(function(x) sum(x, na.rm = T))
#> Joining with `by = join_by(document, term, n)`
#> Adding missing grouping variables: `document`

# reformat for easy plotting
augment_predict <- 
  as_tibble(t(augment_predict[, -c(1,2)]), .name_repair = "minimal")

colnames(augment_predict) <- unique(tidy_docs$document)

augment_predict$topic <- 1:nrow(augment_predict) %>% as.factor()

compare_mat <- 
  augment_predict %>%
  select(
    topic,
    augment = matches(rownames(d2)[1])
  ) %>%
  mutate(
    augment = augment / sum(augment), # normalize to sum to 1
    dot = p_dot[1, ]
  ) %>%
  pivot_longer(cols = c(augment, dot))

ggplot(compare_mat) + 
  geom_bar(aes(y = value, x = topic, group = name, fill = name), 
           stat = "identity", position = "dodge") +
  labs(title = "Prediction using 'augment' vs 'predict(..., method = \"dot\")'")
```

<img src="man/figures/README-example-3.png" alt="" width="100%" />

``` r

# Not shown: aggregating over documents results in recovering the "tidy" lambda.

### updating the model ----
# now that you have new documents, maybe you want to fold them into the model?
lda2 <- refit(
  object = lda, 
  new_data = d, # save me the trouble of manually-combining these by just using d
  iterations = 200, 
  burnin = 175,
  calc_likelihood = TRUE,
  calc_r2 = TRUE
)

# we can do similar analyses
# did the model converge?
qplot(x = iteration, y = log_likelihood, data = lda2$log_likelihood, geom = "line") +
  ggtitle("Checking model convergence")
```

<img src="man/figures/README-example-4.png" alt="" width="100%" />

``` r

# look at the model overall
glance(lda2)
#> # A tibble: 1 × 5
#>   num_topics num_documents num_tokens iterations burnin
#>        <int>         <int>      <int>      <dbl>  <dbl>
#> 1         10            99       2962        200    175

print(lda2)
#> A Latent Dirichlet Allocation Model of  10 topics,  99  documents, and  2962  tokens:
#> refit.tidylda(object = lda, new_data = d, iterations = 200, burnin = 175, 
#>     calc_likelihood = TRUE, calc_r2 = TRUE)
#> 
#> The model's R-squared is  0.1318 
#> The  5  most prevalent topics are:
#> # A tibble: 10 × 4
#>   topic prevalence coherence top_terms                                          
#>   <dbl>      <dbl>     <dbl> <chr>                                              
#> 1     5      14.1     0.172  research, cancer, program, disparities, dr, ...    
#> 2     3      13.4     0.129  cell, cells, lung, hiv, inflammation, ...          
#> 3    10      12.0     0.0980 microbiome, brain, plasticity, response, immune, .…
#> 4     7      11.0     0.118  cancer, tumor, clinical, cells, dcis, ...          
#> 5     9       9.67    0.0386 data, disease, study, develop, core, ...           
#> # ℹ 5 more rows
#> 
#> The  5  most coherent topics are:
#> # A tibble: 10 × 4
#>   topic prevalence coherence top_terms                                          
#>   <dbl>      <dbl>     <dbl> <chr>                                              
#> 1     4       4.06     0.483 cdk5, nmdar, nmdars, adhesion, calpain, ...        
#> 2     8       8.88     0.286 mitochondrial, function, cmybp, sleep, studies, ...
#> 3     5      14.1      0.172 research, cancer, program, disparities, dr, ...    
#> 4     6       8.72     0.144 health, risk, diabetes, sud, factors, ...          
#> 5     3      13.4      0.129 cell, cells, lung, hiv, inflammation, ...          
#> # ℹ 5 more rows


# how does that compare to the old model?
print(lda)
#> A Latent Dirichlet Allocation Model of  10 topics,  50  documents, and  1524  tokens:
#> tidylda(data = d1, k = 10, iterations = 200, burnin = 175, alpha = 0.1, 
#>     eta = 0.05, optimize_alpha = FALSE, calc_likelihood = TRUE, 
#>     calc_r2 = TRUE, return_data = FALSE)
#> 
#> The model's R-squared is  0.2374 
#> The  5  most prevalent topics are:
#> # A tibble: 10 × 4
#>   topic prevalence coherence top_terms                                          
#>   <dbl>      <dbl>     <dbl> <chr>                                              
#> 1     7       12.7     0.198 cancer, dcis, breast, imaging, clinical, ...       
#> 2     8       12.0     0.276 mitochondrial, function, cmybp, sleep, studies, ...
#> 3     5       11.7     0.222 research, cancer, disparities, program, mhirt, ... 
#> 4     1       11.4     0.124 effects, imaging, v4, stiffening, wall, ...        
#> 5     6       10.3     0.14  risk, diabetes, sud, numeracy, redox, ...          
#> # ℹ 5 more rows
#> 
#> The  5  most coherent topics are:
#> # A tibble: 10 × 4
#>   topic prevalence coherence top_terms                                          
#>   <dbl>      <dbl>     <dbl> <chr>                                              
#> 1     4       8.17     0.466 cdk5, nmdar, nmdars, adhesion, calpain, ...        
#> 2     8      12.0      0.276 mitochondrial, function, cmybp, sleep, studies, ...
#> 3     5      11.7      0.222 research, cancer, disparities, program, mhirt, ... 
#> 4     7      12.7      0.198 cancer, dcis, breast, imaging, clinical, ...       
#> 5     3       8.43     0.142 cells, cell, specific, infection, lung, ...        
#> # ℹ 5 more rows
```

There are several vignettes available in
[/vignettes](https://GitHub.com/TommyJones/tidylda/tree/main/vignettes).
They can be compiled using `knitr` or you can view [PDF versions on
CRAN](https://CRAN.R-project.org/package=tidylda).

See NEWS.md for a changelog, including changes from the CRAN release to
the development version on GitHub.

See the “Issues” tab on GitHub to see planned features as well as bug
fixes.

# Contributions

If you would like to contribute to this package, please start by opening
an issue on GitHub. Direct contributions via merge requests are
discouraged unless invited to do so.

If you have any suggestions for additional functionality, changes to
functionality, changes to arguments or other aspects of the API please
let me know by opening an issue on GitHub or sending me an email:
jones.thos.w at gmail.com.

## Code of Conduct

Please note that the tidylda project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
