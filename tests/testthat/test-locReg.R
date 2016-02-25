context("Local Regressions w/ locReg()")

win      <-  c(sample(1:50,1), sample(51:100,1))
x        <-  rnorm(100)
y        <-  rnorm(100, mean=2, sd=0.3) + rnorm(100, mean=3, sd=0.5) * x
x2       <-  rnorm(50)
y2       <-  rnorm(50, mean=2, sd=0.3) + rnorm(50, mean=3, sd=0.5) * x2
xNA      <-  x
yNA      <-  y
xNA[sample(c(1:100), 5)] <- NA
yNA[sample(c(1:100), 5)] <- NA
charVec  <-  letters[1:100]
locReg1  <-  locReg(wins=win, xall=x, yall=y, resids=FALSE)

test_that("Simple corner cases", {

    # returns correct behaviour regardless of input
    expect_is(locReg1, "data.frame")
    expect_identical(dim(locReg1), c(1L,9L))
    #expect_error(locReg(win, x, charVec), "NA/NaN/Inf in 'x'")
    #expect_warning(locReg(win, x, charVec), "NAs introduced by coercion")

    # returns correct behaviour when lengths of x and y differ
    expect_error(locReg(win, x, y2), "x and y must be of equal length")
})
