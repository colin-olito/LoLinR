context("Sample skewness")

test_that("Simple corner cases", {
    x        <-  rnorm(10)
    x2       <-  c(NA, x)
    charVec  <-  letters[1:10]

    # returns correct behaviour regardless of NA strings on vector
    expect_is(skew(x), "numeric")
    expect_is(skew(x2), "numeric")

    # make sure output length matches expected behaviour
    ## formula breaks if length(x) == 1 | length(x) == 2
    expect_true(is.na(skew(1)))
    expect_true(is.na(skew(1:2)))

    # returns error if x is not numeric
    expect_error(skew(charVec), "x must be numeric")

    # returns warning in case x has non-unique values
    expect_true(is.na(skew(x2, na.rm=FALSE)))
})
