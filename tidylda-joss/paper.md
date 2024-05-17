---
title: 'tidylda: An R Package for Latent Dirichlet Allocation Using ''tidyverse''
  Conventions'
tags:
- R
- topic models
- LDA
- natural language processing
- tidy data
date: "9 May 2024"
output: pdf_document
authors:
- name: Tommy Jones
  orcid: "0000-0001-6457-2452"
  equal-contrib: yes
  affiliation: 1
bibliography: paper.bib
affiliations:
- name: Foundation, USA
  index: 1
---

# Summary
_tidylda_ is a package for a topic model, Latent Dirichlet Allocation or LDA [@blei2002lda], that is natively compatible with the _tidyverse_ [@tidyverse]. _tidylda_'s Gibbs sampler is written in C++ for performance and offers several novel features. It also has methods for sampling from the posterior of a trained model, for more traditional Bayesian analyses.

_tidylda_ takes its syntactic cues from an ecosystem of R packages known as _the tidyverse_. The tidyverse's goal is to "facilitate a conversation between a human and computer about data" [@tidyverse]. Packages in---and adjacent to---the tidyverse share a common design philosophy and syntax based on "tidy data" principles [@wickham2014tidy]. Tidy data has each variable in a column, each observation in a row, and each observational unit in a table. Extensions include the _broom_ package [@broom] for "tidying" up outputs from statistical models and the in-development _tidymodels_ ecosystem [@tidymodels] which extends the tidyverse philosophy to statistical modeling and machine learning workflows. 

Silge et al. articulated a "tidy data" framework for text analyses---the _tidytext_ package [@tidytextjoss]. Their approach has "one row per document per token". The _tidytext_ package provides functionality to tokenize a corpus, transform it into this "tidy" format, and manipulate it in various ways, including preparing data for input into some of R's many topic modeling packages. The _tidytext_ package also provides tidying functions in the style of _broom_ to harmonize outputs from some of R's topic modeling packages. _tidylda_ manages inputs and outputs in the flavor of _tidytext_ but in one self contained package.

# Statement of Need
Packages that implement topic models in R and other languages are plentiful. Why do we need another? _tidylda_'s native compatibility with the _tidyverse_ makes it significantly more user friendly than other topic modeling packages. It also enables more traditional Bayesian analyses such as the ability to set more flexible priors, burn-in iterations, averaging over segments of the Gibbs sample chain, and sampling from the posterior that other packages lack. Finally, _tidylda_ implements a transfer learning algorithm developed in [@jones2023latent], unavailable in any other package and described in more detail in the following section.

# Latent Dirichlet Allocation and Notation
LDA is a Bayesian latent variable model for text [@blei2002lda]. It decomposes a data set of word counts, $\boldsymbol{X}$, whose row/column entries, $d,v$, represent the number of times word $v$ is found in document $d$, into two matrices: $\boldsymbol\Theta$ and $\boldsymbol{B}$. The former gives a distribution of (latent) topics over documents and the latter gives a distribution of words over topics. Formally, LDA is

\begin{align}
  z_{d_{n}}|\boldsymbol\theta_d &\sim 
  \text{Categorical}(\boldsymbol\theta_d)\\
  w_{d_{n}}|z_{k},\boldsymbol\beta_k^{(t)} &\sim
  \text{Categorical}(\boldsymbol\beta_k^{(t)}) \\
  \boldsymbol\theta_d &\sim
  \text{Dirichlet}(\boldsymbol\alpha)\\
  \boldsymbol\beta_k^{(t)} &\sim
  \text{Dirichlet}(\boldsymbol\eta)
\end{align}

where random variables $w_{d_{n}}$ and $z_{d_{n}}$ represent the word and topic of the $n$-th word of the $d$-th document. The user sets prior values for $\boldsymbol\alpha$ and $\boldsymbol\eta$ as well as specifying the number of topics, $K$. 

Posterior estimates of $\boldsymbol\Theta$ and $\boldsymbol{B}$ along with the data, $\boldsymbol{X}$, allow for the calculation of $\boldsymbol\Lambda$. Where $\beta_{k,v}$ is $P(\text{word}_v | \text{topic}_k)$, $\lambda_{k, v}$ is $P(\text{topic}_k | \text{word}_v)$. 

One of the common ways of estimating LDA is through collapsed Gibbs sampling. In the background, the sampler tracks the number of times topics are sampled with two matrices: $\boldsymbol{Cd}$ and $\boldsymbol{Cv}$. The former's row/column entries, $d,k$, are the number of times topic $k$ was sampled in document $d$. The latter's row/column entries, $k,v$ are the number of times topic $k$ was sampled for word $v$.

## Transfer LDA (tLDA)

Formally, tLDA modifies LDA in the following way:

\begin{align}
  \boldsymbol\beta_k^{(t)} &\sim
  \text{Dirichlet}(\omega_k^{(t)} \cdot \mathbb{E}\left[\boldsymbol\beta_k^{(t-1)}\right])
\end{align}

The above indicates that tLDA places a matrix prior for words over topics where $\eta_{k, v}^{(t)} = \omega_{k}^{(t)} \cdot \mathbb{E}\left[\beta_{k,v}^{(t-1)}\right] = \omega_{k}^{(t)} \cdot \frac{Cv_{k,v}^{(t-1)} + \eta_{k,v}^{(t-1)}}{\sum_{v=1}^V Cv_{k,v}^{(t-1)}}$. Because the posterior at time $t$ depends only on data at time $t$ and the state of the model at time $t-1$, tLDA models retain the Markov property. Rather than $K$ tuning weights, $\omega_k^{(t)}$, users tune a single parameter, $a$. 

When $a^{(t)} = 1$, fine tuning is equivalent to adding the data in $\boldsymbol{X}^{(t)}$ to $\boldsymbol{X}^{(t-1)}$. In other words, each word occurrence in $\boldsymbol{X}^{(t)}$ carries the same weight in the posterior as each word occurrence in $\boldsymbol{X}^{(t-1)}$. When $a^{(t)} < 1$, then the posterior has recency bias. When $a^{(t)} > 1$, then the posterior has precedent bias. Each word occurrence in $\boldsymbol{X}^{(t)}$ carries less weight than each word occurrence in $\boldsymbol{X}^{(t-1)}$. 

For more details on tLDA see [@jones2023latent].

# _tidylda_'s Novel Features

## Model Initialization and Gibbs Sampling

_tidylda_'s Gibbs sampler has several unique features, described below.

**Non-uniform initialization:** Most LDA Gibbs samplers initialize by assigning words to topics and topics to documents by sampling from a uniform distribution. This ensures initialization without incorporating any prior information. _tidylda_ incorporates the priors in its initialization. It begins by drawing $P(\text{topic}|\text{document})$ and $P(\text{word}|\text{topic})$ from Dirichlet distributions with parameters $\boldsymbol\alpha$ and $\boldsymbol\eta$, respectively. Then _tidylda_ uses the above probabilities to construct $P(\text{topic}|\text{word}, \text{document})$ and makes a single run of the Gibbs sampler to initialize two matrices tracking topics over documents and words over topics, denoted $\boldsymbol{Cd}$ and $\boldsymbol{Cv}$, respectively. 

This non-uniform initialization powers tLDA, described in Chapter 5, by starting a Gibbs run near where the previous run left off. For initial models, it uses the user's prior information to tune where sampling starts.

**Flexible priors:** _tidylda_ has multiple options for setting LDA priors. Users may set scalar values for $\boldsymbol\alpha$ and $\boldsymbol\eta$ to construct symmetric priors. Users may also choose to construct vector priors for both $\boldsymbol\alpha$ and $\boldsymbol\eta$ for a full specification of LDA. Additionally, _tidylda_ allows users to set a matrix prior for $\boldsymbol\eta$, enabled by its implementation of tLDA. This enables users to set priors over word-topic relationships informed by expert input. The best practices for encoding expert input in this manner are not yet well studied. Nevertheless, this capability makes _tidylda_ unique among LDA implementations. 

**Burn in iterations and posterior averaging:** Most LDA Gibbs samplers construct posterior estimates of $\boldsymbol\Theta$ and $\boldsymbol{B}$ from $\boldsymbol{Cd}$ and $\boldsymbol{Cv}$'s values of the final iteration of sampling, effectively using a single sample. This is inconsistent with best practices from Bayesian statistics, which is to average over many samples from a stable posterior. _tidylda_ enables averaging across multiple samples of the posterior with the `burnin` argument. When `burnin` is set to a positive integer, _tidylda_ averages the posterior across all iterations larger than `burnin`. For example, if `iterations` is 200 and `burnin` is 150, _tidylda_ will return a posterior estimate that is an average of the last 50 sampling iterations. This ensures that posterior estimates are more likely to be representative than any single sample.

**Transfer learning with tLDA:** Finally, and as discussed previously, _tidylda_'s Gibbs sampler enables transfer learning with tLDA. 

## Tidy Methods

_tidylda_'s construction follows _Conventions of R Modeling Packages_ [@tidymodelsbook]. In particular, it contains methods for `print`, `summary`, `glance`, `tidy`, and `augment`, consistent with other "tidy" packages. These methods are briefly described below.

* `print`, `summary`, and `glance` return various summaries of the contents of a _tidylda_ object, into which an LDA model trained with _tidylda_ is stored.
* `tidy` returns the contents of $\boldsymbol\Theta$, $\boldsymbol{B}$, or $\boldsymbol\Lambda$ (stored as `theta`, `beta`, and `lambda` respectively), as specified by the user, formatted as a tidy `tibble`, instead of a numeric matrix.
* `augment` appends model outputs to observational-level data. Taking the cue from _tidytext_ [@tidytextjss], "observational-level" data is one row per word per document. Therefore, the key statistic used by `augment` is $P(\text{topic}|\text{word}, \text{document})$. _tidylda_ calculates this as $\boldsymbol\Lambda \times P(\text{word}|\text{document})$, where $P(\text{word}|\text{document}_d) = \frac{\boldsymbol{x}_d}{\sum_{v=1}^V x_{d,v}}$. 

## Posterior Methods
_tidylda_ enables traditional Bayesian uncertainty quantification by sampling from the posterior. The posterior distribution for $\boldsymbol\theta_d$ is $\text{Dirichlet}(\boldsymbol{Cd}_d + \boldsymbol\alpha)$ and the posterior distribution for $\boldsymbol\beta_k$ is $\text{Dirichlet}(\boldsymbol{Cv}_k + \boldsymbol\eta)$ (or $\text{Dirichlet}(\boldsymbol{Cv}_k + \boldsymbol\eta_k)$ for tLDA). _tidylda_ enables a `posterior` method for _tidylda_ objects, allowing users to sample from the posterior to quantify uncertainty for estimates of estimated parameters.

_tidylda_ uses one of two calculations for predicting topic distributions (i.e., $\hat{\boldsymbol\theta}_d$) for new documents. The first, and default, is to run the Gibbs sampler, constructing a new $\boldsymbol{Cd}$ for the new documents but without updating topic-word distributions in $\boldsymbol{B}$. The second uses a dot product, $\boldsymbol{X}^{(new)} \cdot \boldsymbol\Lambda'$, where the rows of $\boldsymbol{X}^{(new)}$ are normalized to sum to $1$. _tidylda_ actually uses the dot product prediction combined with the _non-uniform initialization_---described above---to initialize $\boldsymbol{Cd}$ when predicting using the Gibbs sampler.

# Acknowledgements
Many people over the years have supported the development of _tidylda_. But most notably are

* Wil Doane, for making me a better programmer and giving ample good advice.
* Brendan Knapp, for helping with the C++ code.
* Barum Park, whose code formed the basis of the multinomial sampler in C++.
* My PhD committee, without whom _tidylda_ would be full of "good ideas", but not peer-reviewed research.

# References


