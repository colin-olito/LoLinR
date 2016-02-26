# Local Linear Regression for Estimating Monotonic Biological Rates in R (LoLinR)

[![Build Status](https://travis-ci.org/colin-olito/LoLinR.png?branch=master)](https://travis-ci.org/colin-olito/LoLinR)

## Documentation

We introduce the `LoLinR` package, which provides tools to implement local linear regression techniques for estimating monotonic rates from time-series or trace data in a statistically robust and reproducable fashion. The methods are a modification of traditional Loess regression techniques built around the wrapper function `rankLocReg()`. 

See the package documentation for `LoLinR` using `library(help=LoLinR)` and examples therein.


## Installation

The `LoLinR` package can be installed from github using the [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) package. 
**Using `devtools::install_github`**

If you do not yet have `devtools`, install with `install.packages("devtools")`

Then install `LoLinR` using the following:
`library(devtools)`
`install_github('colin-olito/LoLinR')`
`library(LoLinR)`

`LoLinR` does not currently have any dependencies on other packages





