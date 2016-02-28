context("Thinning data sets")

data(thinned_cormorant_data)
head(thinned_cormorant_data)
thin2  <-  thinData(thinned_cormorant_data$Time,
					thinned_cormorant_data$VO2.ml.min, by=2)
thin3  <-  thinData(thinned_cormorant_data$Time,
					thinned_cormorant_data$VO2.ml.min, by=3)
N      <-  nrow(thinned_cormorant_data)
len21  <-  length(seq(from=1, to=N, by=2))
len22  <-  length(seq(from=2, to=N, by=2))
len31  <-  length(seq(from=1, to=N, by=3))
len32  <-  length(seq(from=2, to=N, by=3))
len33  <-  length(seq(from=3, to=N, by=3))

source('R/functions.R')

test_that("Simple corner cases", {

    # returns correct structure
    expect_is(thin2, "list")
    expect_identical(length(thin2), 2L)
    expect_identical(length(thin3), 3L)
    expect_identical(nrow(thin2[[1]]), len21)
    expect_identical(nrow(thin2[[2]]), len22)
    expect_identical(nrow(thin3[[1]]), len31)
    expect_identical(nrow(thin3[[2]]), len32)
    expect_identical(nrow(thin3[[3]]), len33)

    # make sure it removes NAs from inputs
    

    # returns correct behaviour regardless of input
    

    # returns correct behaviour when lengths of x and y differ
 	

	# correct behaviour when "thin" is missing, defaults to 2
    

	# Returns data frame with equal spacing of observations

#	expect_identical((nrow(thinned_cormorant_data) / nrow(newData)), 3)


    # returns correct structure
#    expect_is(allRegs1, "rankLocReg")
#    expect_is(allRegs2, "rankLocReg")
#    expect_is(allRegs3, "rankLocReg")
#    expect_identical(length(allRegs1), 6L)
#    expect_identical(length(allRegs2), 6L)
#    expect_identical(length(allRegs3), 6L)
#
#    # make sure it removes NAs from inputs
#    expect_identical(length(allRegs3$xall), length(xNA[!is.na(xNA) & !is.na(yNA)]))
#
#    # returns correct behaviour regardless of input
#    expect_error(rankLocReg(xall=charVec, yall=y, alpha=0.8, method='eq'), 'x must be numeric')
#    expect_error(rankLocReg(xall=x, yall=charVec, alpha=0.8, method='eq'), 'x must be numeric')
#
#    # returns correct behaviour when lengths of x and y differ
#    expect_error(rankLocReg(xall=x, yall=y[-c(10)], alpha=0.8, method='eq'))
# 	
#  	 # throws error when missing alpha
#    expect_error(rankLocReg(xall=x[order(x)], yall=y[order(x)]), 'argument "alpha" is missing, with no default')
#
#    # correct behaviour when method is missing, defaults to 'ns'
#    expect_is(rankLocReg(xall=x[order(x)], yall=y[order(x)], alpha=0.8), "rankLocReg")
#    expect_identical(rankLocReg(xall=x[order(x)], yall=y[order(x)], alpha=0.8)$method, "ns")
#
#    # throws warning when x is unordered
#    expect_warning(rankLocReg(xall=x, yall=y, alpha=0.8, method='eq'), "Dataset must be ordered by xall")
})