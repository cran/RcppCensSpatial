
<!-- README.md is generated from README.Rmd. Please edit that file -->

## The `RcppCensSpatial R` package

<!-- badges: start -->
<!-- badges: end -->

The `RcppCensSpatial R` package provides functions to estimate the
parameters in spatial models with censored and missing responses via EM
(Dempster, Laird, and Rubin 1977), SAEM (Delyon, Lavielle, and Moulines
1999), and MCEM (Wei and Tanner 1990) algorithm. These algorithms are
widely used to compute the maximum likelihood (ML) estimates in
incomplete data problems. The EM algorithm computes the ML estimates
when a closed expression for the conditional expectation of the
complete-data log-likelihood function is available. In the MCEM
algorithm, the conditional expectation is substituted by a Monte Carlo
(MC) approximation based on many independent simulations of the missing
data, while the SAEM algorithm splits the E-step into a simulation step
and an integration step. The SAEM algorithm was developed as an
alternative to the computationally intensive MCEM algorithm.

This package approximates the standard error of the estimates using the
method developed by (Louis 1982). It also has a function that performs
spatial prediction in a set of new locations. Besides the functions to
estimate parameters, this package allows calculating the covariance
matrix and the distance matrix.

### Functions

`RcppCensSpatial` package provides the following functions:

-   `CovMat`: Computes the spatial covariance matrix.
-   `dist2Dmatrix`: Computes the Euclidean distance matrix for a set of
    coordinates.
-   `EM.sclm`: Returns the ML estimates of the unknown parameters
    computed via the EM algorithm.
-   `MCEM.sclm`: Returns the ML estimates of the unknown parameters
    computed via the MCEM algorithm.
-   `SAEM.sclm`: Returns the ML estimates of the unknown parameters
    computed via the SAEM algorithm.
-   `predict.sclm`: Performs spatial prediction in a set of new
    locations.
-   `rCensSp`: Simulates censored spatial data for an established
    censoring rate.

Functions `print`, `summary`, `predict`, and `plot` also work for the
`sclm` class.

### Installation

You can install the released version of `RcppCensSpatial` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("RcppCensSpatial")
```

### References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-delyon1999convergence" class="csl-entry">

Delyon, Bernard, Marc Lavielle, and Eric Moulines. 1999. “Convergence of
a Stochastic Approximation Version of the EM Algorithm.” *Annals of
Statistics* 27 (1): 94–128. <https://www.jstor.org/stable/120120>.

</div>

<div id="ref-dempster1977maximum" class="csl-entry">

Dempster, Arthur P, Nan M Laird, and Donald B Rubin. 1977. “Maximum
Likelihood from Incomplete Data via the EM Algorithm.” *Journal of the
Royal Statistical Society: Series B (Methodological)* 39 (1): 1–22.
<https://www.jstor.org/stable/2984875>.

</div>

<div id="ref-louis1982finding" class="csl-entry">

Louis, Thomas. 1982. “Finding the Observed Information Matrix When Using
the EM Algorithm.” *Journal of the Royal Statistical Society: Series B
(Methodological)* 44 (2): 226–33.
<https://www.jstor.org/stable/2345828>.

</div>

<div id="ref-wei1990monte" class="csl-entry">

Wei, Greg CG, and Martin A Tanner. 1990. “A Monte Carlo Implementation
of the EM Algorithm and the Poor Man’s Data Augmentation Algorithms.”
*Journal of the American Statistical Association* 85 (411): 699–704.
<https://doi.org/10.1080/01621459.1990.10474930>.

</div>

</div>
