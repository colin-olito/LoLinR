context("Breusch-Godfrey Statistic")

x        <-  rnorm(100)
y        <-  rnorm(100, mean=runif(c(-2,2)), sd=0.3) + rnorm(100, mean=runif(c(-3,3)), sd=0.5) * x
x2       <-  rnorm(103)
y2       <-  rnorm(103, mean=2, sd=0.3) + rnorm(103, mean=3, sd=0.5) * x2
xNA      <-  x
yNA      <-  y
charVec  <-  letters[1:100]

xNA[sample(c(1:100), 5)] <- NA
yNA[sample(c(1:100), 5)] <- NA

test_that("Simple corner cases", {

    # returns correct structure
    expect_is(breuschGodfrey(y, x), "list")
    expect_identical(length(breuschGodfrey(y, x)), 3L)

    # returns correct behaviour when lengths of x and y differ
    expect_error(breuschGodfrey(y2, x), "incompatible dimensions")

    # returns error if x is not numeric
    expect_error(breuschGodfrey(y, charVec), "NA/NaN/Inf in 'x'")

    # returns correct behaviour when x, y have NAs
    expect_error(breuschGodfrey(yNA, xNA), "NA/NaN/Inf in 'x'")
})

test_that("Computations are correct", {
 
   # does calculated statistic match output from model function bgtest() from package lmtest?
    expect_equivalent(breuschGodfrey(y, x)$bg, bgtest(y ~ x, order=(length(x)-3))$statistic)

   # Check that bfN is being calculated as expected
    expect_equivalent(breuschGodfrey(y, x)$bg / length(x), breuschGodfrey(y, x)$bgN)
 
})
