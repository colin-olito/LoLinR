context("Local Regressions w/ locReg()")

rm(list=ls())
source('R/functions.R')
	win      <-  c(sample(1:50,1), sample(51:100,1))
    x        <-  rnorm(100)
    y        <-  rnorm(100, mean=2, sd=0.3) + rnorm(100, mean=3, sd=0.5) * x
    x2       <-  rnorm(103)
    y2       <-  rnorm(103, mean=2, sd=0.3) + rnorm(103, mean=3, sd=0.5) * x2
    xNA      <-  x
    yNA      <-  y
    xNA[sample(c(1:100), 5)] <- NA
    yNA[sample(c(1:100), 5)] <- NA
    charVec  <-  letters[1:100]

locReg(wins=win, xall=x, yall=y, resids=FALSE)

test_that("Simple corner cases", {
    # returns correct behaviour regardless of input
    expect_is(breuschGodfrey(y, x), "list")

    # make sure output length matches expected behaviour
    expect_identical(length(breuschGodfrey(y, x)), 3L)

    # returns correct behaviour when lengths of x and y differ
    expect_error(breuschGodfrey(y2, x), "incompatible dimensions")

    # returns error if x is not numeric
    expect_error(breuschGodfrey(y, charVec), "NA/NaN/Inf in 'x'")
    expect_warning(breuschGodfrey(y, charVec), "NAs introduced by coercion")

    # returns correct behavoiur when x, y have NAs
    expect_error(breuschGodfrey(yNA, xNA), "NA/NaN/Inf in 'x'")
})
