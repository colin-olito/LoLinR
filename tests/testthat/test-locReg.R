context("Local Regressions w/ locReg()")

win      <-  c(sample(1:50,1), sample(51:100,1))
x        <-  rnorm(100)
y        <-  rnorm(100, mean=runif(c(-2,2)), sd=0.3) + rnorm(100, mean=runif(c(-3,3)), sd=0.5) * x
x2       <-  rnorm(50)
y2       <-  rnorm(50, mean=2, sd=0.3) + rnorm(50, mean=3, sd=0.5) * x2
xNA      <-  x
yNA      <-  y
xNA[sample(c(1:100), 5)] <- NA
yNA[sample(c(1:100), 5)] <- NA
charVec  <-  letters[1:100]
locReg1  <-  locReg(wins=win, xall=x, yall=y, resids=FALSE)
lm1      <-  lm( y[win[1]:win[2]] ~ x[win[1]:win[2]])
lst      <-  list(win,x,y)

test_that("Simple corner cases", {

    # returns correct structure
    expect_is(locReg1, "data.frame")

    # returns correct behaviour regardless of input
    expect_identical(dim(locReg1), c(1L,9L))
    expect_error(locReg(win, x, charVec), "NA/NaN/Inf in 'y'")
    expect_error(locReg(win, xNA, yNA), "NA/NaN/Inf in 'x'")
    expect_error(locReg(win, x, yNA), "NA/NaN/Inf in 'y'")

    # returns correct behaviour when lengths of x and y differ
    expect_error(locReg(win, x, y2), "x and y must be of equal length")
  	
  	# Throws error when missing args
    expect_error(locReg(wins=win, xall=x), 'argument "yall" is missing, with no default')

    # throws error if user-provided windows are invalid
    expect_error(locReg(wins=rev(win), xall=x, yall=y), 'Left window boundary greater than or equal to Right boundary')
})

test_that("Basic Computations are Correct", {

	# returns same results as lm()
	expect_equivalent(locReg1$b0, lm1$coefficients[1])
	expect_equivalent(locReg1$b1, lm1$coefficients[2])
	expect_equivalent(locReg1$b1LoCI,confint(lm1)[2,1])
	expect_equivalent(locReg1$b1UpCI,confint(lm1)[2,2])
})
