
<!-- README.md is generated from README.Rmd. Please edit that file -->

# natstrat

<!-- badges: start -->

[![R-CMD-check](https://github.com/kkbrum/natstrat/workflows/R-CMD-check/badge.svg)](https://github.com/kkbrum/natstrat/actions)
[![Codecov test
coverage](https://codecov.io/gh/kkbrum/natstrat/branch/master/graph/badge.svg)](https://codecov.io/gh/kkbrum/natstrat?branch=master)
<!-- badges: end -->

The goal of natstrat is to obtain unweighted natural strata that balance
many covariates. Natural strata fix a constant ratio of controls to
treated units within each stratum. This ratio need not be an integer.
The control units are chosen using randomized rounding of a linear
program that balances many covariates.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kkbrum/natstrat", build_vignettes = TRUE)
```

## Usage

To learn about how to use this package, please see the associated
vignette with:

``` r
browseVignettes(package = "natstrat")
```
