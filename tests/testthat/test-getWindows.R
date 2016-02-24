context("Getting all possible windows")

test_that("Simple corner cases", {
    x        <-  rnorm(10)
    x2       <-  c(NA, x)
    charVec  <-  letters[1:10]

    # function breaks if alpha argument is missing
    expect_error(getWindows(x))
    
    # returns correct behaviour regardless of NA strings on vector
    expect_is(getWindows(x, alpha=.2), "matrix")
    expect_is(getWindows(x2, alpha=.2), "matrix")

    # function breaks if alpha is wrongly specified
    expect_error(getWindows(x, alpha=2), "alpha must take a value higher than 0 and lower or equal to 1")
    expect_error(getWindows(x, alpha=0), "alpha must take a value higher than 0 and lower or equal to 1")
    expect_error(getWindows(x, alpha=-1), "alpha must take a value higher than 0 and lower or equal to 1")
    expect_error(getWindows(x, alpha='a'), "alpha must take a value higher than 0 and lower or equal to 1")

    # returns error if x is not numeric
    expect_error(getWindows(charVec, alpha=.2), "x must be numeric")

    # returns error if length(x) == 1
    expect_error(getWindows(1, alpha=.2), "n < m")
})
