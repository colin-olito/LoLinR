context("Thinning data sets")

data(CormorantData)
thin2  <-  thinData(CormorantData, xy=c(1,2), by=2)
thin3  <-  thinData(CormorantData, xy=c(1,2), by=3)
N      <-  nrow(CormorantData)
len21  <-  length(seq(from=1, to=N, by=2))
len22  <-  length(seq(from=2, to=N, by=2))
len31  <-  length(seq(from=1, to=N, by=3))
len32  <-  length(seq(from=2, to=N, by=3))
len33  <-  length(seq(from=3, to=N, by=3))

x        <-  c(1:100)
y        <-  sort(runif(100, min=8, max=50))
y2       <-  y[1:90]
yNA      <-  y
yNA[sample(c(1:100), 5)] <- NA
yNAtest  <-  yNA[!is.na(yNA)]
charVEC  <-  sample(letters, 100, replace=TRUE)
testDF   <-  data.frame("x"=x, "y"=y, "yNA"=yNA, "charVEC"=as.character(charVEC))
thinNA   <- thinData(testDF, xy=c(1,3), by=2)

data(UrchinData)
thinTest  <-  thinData(UrchinData)

test_that("Simple corner cases", {
    # returns correct structure
    expect_is(thin2, "list")
    expect_identical(length(thin2), 2L)
    expect_identical(length(thin3), 3L)
    expect_identical(nrow(thin2$newData1), len21)
    expect_identical(nrow(thin2$newData2), len22)
    expect_identical(nrow(thin3$newData1), len31)
    expect_identical(nrow(thin3$newData2), len32)
    expect_identical(nrow(thin3$newData3), len33)

    # make sure it removes NAs from inputs
	expect_match(any(is.na(thin2$newData2)), "FALSE")
	expect_match(any(is.na(thin3$newData2)), "FALSE")
	expect_match(any(is.na(thin3$newData3)), "FALSE")
	expect_identical(length(thinNA$newData2$yNA),
					 length(yNAtest[seq(from=2, to=length(yNAtest), by=2)]))

    # returns correct behaviour regardless of input
    expect_error(thinData(testDF, xy=c(1,4)), 'data must contain only numeric variables')

	# correct behaviour when "thin" is missing, defaults to first two columns, by=2
    expect_identical(length(thinTest), 2L)
    expect_identical(names(UrchinData)[1:2], names(thinTest$newData1))
})