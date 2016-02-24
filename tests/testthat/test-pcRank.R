context("Percentile ranking")

test_that("Simple corner cases", {
    x        <-  rnorm(10)
    x2       <-  c(NA, x)
    z        <-  seq_len(length(x))
    charVec  <-  letters[1:10]
    d1       <-  1:100
    p1       <-  sample(1:5, 100, replace=TRUE)
    e1       <-  sort(rexp(100, rate=.5))

    # returns correct behaviour regardless of NA strings on vector
    expect_is(pcRank(x), "numeric")
    expect_is(pcRank(x2), "numeric")

    # make sure output length matches expected behaviour
    expect_equal(length(x), length(pcRank(x)))
    expect_equal(length(x2[!is.na(x2)]), length(pcRank(x2)))

    # returns error if x is not numeric
    expect_error(pcRank(charVec), "x must be numeric")

    # make sure that the function ranks numeric vectors as expected
    expect_identical(pcRank(sort(x)), pcRank(z))
    expect_identical(pcRank(sort(d1)), pcRank(e1))

    # returns warning in case x has non-unique values
    expect_warning(pcRank(p1), 'input/output have ties')
})
