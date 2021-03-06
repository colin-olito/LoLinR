# Local Linear Regression for Estimating Monotonic Biological Rates in R (LoLinR)

[![Build Status](https://travis-ci.org/colin-olito/LoLinR.png?branch=master)](https://travis-ci.org/colin-olito/LoLinR)

## Documentation

We introduce the `LoLinR` package, which provides tools to implement local linear regression techniques for estimating monotonic rates from time-series or trace data in a statistically robust and reproducible fashion. The methods are a modification of traditional Loess regression techniques built around the wrapper function `rankLocReg()`. 

See the package documentation for `LoLinR` through our [online](https://colin-olito.github.io/LoLinR/vignettes/LoLinR.html) `html_vignette`. You can also find the final version of the peer reviewed article describing the methods [here](http://jeb.biologists.org/content/early/2016/12/27/jeb.148775).

## Citing `LoLinR`

When using `LoLinR` for your research, please cite:

Olito, C., C.R. White, D.J. Marshall, and D.R. Barneche. 2017. Estimating monotonic rates from biological data using local linear regression. Journal of Experimental Biology. DOI: 10.1242/jeb.148775.

## Installation

The `LoLinR` package can be installed from github using the [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) package
using `devtools::install_github`.

If you do not yet have `devtools`, install with `install.packages("devtools")`.

Then install `LoLinR` using the following:  
`library(devtools)`  
`install_github('colin-olito/LoLinR')`  
`library(LoLinR)`

`LoLinR` does not currently have any functional dependencies on other packages. However, it loads `lmtest` for the purposes of unit testing.

## Contact & bug reporting

LoLinR is in the final stages of development, and will continue to get frequent updates until this process is finished. We currently need beta testing, and encourage users to test the package and report any bugs and/or problems by opening an issue on the `LoLinR` github webpage [here](https://github.com/colin-olito/LoLinR/issues). If you would like to report a bug/issue, and do not have a github account (and don't want to get one), please send a brief email detailing the problem you encountered to colin.olito<at>monash.edu.





