context("Multiple Local Regressions w/ rankLocReg()")

x        <-  rnorm(100)
y        <-  rnorm(100, mean=runif(c(-2,2)), sd=0.3) + rnorm(100, mean=runif(c(-3,3)), sd=0.5) * x
x2       <-  rnorm(50)
y2       <-  rnorm(50, mean=2, sd=0.3) + rnorm(50, mean=3, sd=0.5) * x2
xNA      <-  x
yNA      <-  y
charVec  <-  letters[1:100]

xNA[sample(c(1:100), 5)] <- NA
yNA[sample(c(1:100), 5)] <- NA

allRegs1  <-  rankLocReg(xall=x, yall=y, alpha=0.8, method='eq')
allRegs2  <-  rankLocReg(xall=x[order(x)], yall=y[order(x)], alpha=0.8, method='eq')
allRegs3  <-  rankLocReg(xall=xNA[order(xNA)], yall=yNA[order(xNA)], alpha=0.8, method='eq')

test_that("Simple corner cases", {

    # returns correct structure
    expect_is(allRegs1, "rankLocReg")
    expect_is(allRegs2, "rankLocReg")
    expect_is(allRegs3, "rankLocReg")
    expect_identical(length(allRegs1), 6L)
    expect_identical(length(allRegs2), 6L)
    expect_identical(length(allRegs3), 6L)

    # make sure it removes NAs from inputs
    expect_identical(length(allRegs3$xall), length(xNA[!is.na(xNA) & !is.na(yNA)]))

    # returns correct behaviour regardless of input
    expect_error(rankLocReg(xall=charVec, yall=y, alpha=0.8, method='eq'), 'x must be numeric')
    expect_error(rankLocReg(xall=x, yall=charVec, alpha=0.8, method='eq'), 'x must be numeric')

    # returns correct behaviour when lengths of x and y differ
    expect_error(rankLocReg(xall=x, yall=y[-c(10)], alpha=0.8, method='eq'))
  	
  	# throws error when missing alpha
    expect_error(rankLocReg(xall=x[order(x)], yall=y[order(x)]), 'argument "alpha" is missing, with no default')

    # correct behaviour when method is missing, defaults to 'z'
    expect_is(rankLocReg(xall=x[order(x)], yall=y[order(x)], alpha=0.8), "rankLocReg")
    expect_identical(rankLocReg(xall=x[order(x)], yall=y[order(x)], alpha=0.8)$method, "z")

    # throws warning when x is unordered
    expect_warning(rankLocReg(xall=x, yall=y, alpha=0.8, method='eq'), "Dataset must be ordered by xall")
})