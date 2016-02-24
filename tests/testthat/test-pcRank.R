library(LoLinR)
context("Inputs are Correct")

test_that("inputs are correct", {
  expect_error(Skew(c("a","b","c")), "non-numeric argument")
  expect_equal(str_length("ab"), 2)
  expect_equal(str_length("abc"), 3)
})

test_that("str_length of factor is length of level", {
  expect_equal(str_length(factor("a")), 1)
  expect_equal(str_length(factor("ab")), 2)
  expect_equal(str_length(factor("abc")), 3)
})

test_that("str_length of missing is missing", {
  expect_equal(str_length(NA), NA_integer_)
  expect_equal(str_length(c(NA, 1)), c(NA, 1))
  expect_equal(str_length("NA"), 2)
})



FindLocLin  <-  function(yall, xall, alpha, ref.b1=FALSE, plots=TRUE, ...) {

    
