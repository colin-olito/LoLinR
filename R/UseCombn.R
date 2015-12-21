################################################################
#  Functions and Dependencies for the primary function FindLocLin(), used to find
#  the 'most linear' local regressions from a given dataset.
#
#  Run these functions



################################################################
#  Dependency -- Function TriCube():
#
# TriCube() is a function to calculate Tri-Cube weights for a given vector of values
# for an independent variable. This function is a dependency for FindLocLin(), and is
# used to calculate the weights for Weighted Least Squares Regression.
TriCube  <-  function(x, h) {
    j  <-  h + 1
    z  <-  abs(x - x[j]) / h
    ifelse(z < 1, (1 - z^3)^3, 0)
}


TriCube  <-  function(x, m, h) {
#browser()
    z  <-  abs(x - m) / h
    ifelse(z < 1, (1 - z^3)^3, 0)
}

TriCube(x=x,5.5,5)


################################################################
#  Dependency -- Function Skew():
#
# Skew() is a function to calculate the Fisher-Pearson Standardized Third Moment
# Coefficient for a given vector of numbers. This function is a dependency for
# FindLocLin(), and is used to calculate the skewness of standardized residuals.

# Note: Compare to moments package
Skew  <-  function(x) {
    n  <-  length(x)
    (n/((n-1)*(n-2))) * sum(((x-mean(x))/sd(x))^3)
}

################################################################
# Dependency -- nCm():
###############

# NOTE: This code is copied directly from the combn() function
#        from the R {utils} package. Not sure if this matters
#        as a dependency (i.e. whether we should strip this
#        function down to just what we need)... need to ask
#        Diego...

nCm <- function (x, m, FUN = NULL, simplify = TRUE, ...) 
{
    stopifnot(length(m) == 1L, is.numeric(m))
    if (m < 0) 
        stop("m < 0", domain = NA)
    if (is.numeric(x) && length(x) == 1L && x > 0 && trunc(x) == 
        x) 
        x <- seq_len(x)
    n <- length(x)
    if (n < m) 
        stop("n < m", domain = NA)
    x0 <- x
    if (simplify) {
        if (is.factor(x)) 
            x <- as.integer(x)
    }
    m <- as.integer(m)
    e <- 0
    h <- m
    a <- seq_len(m)
    nofun <- is.null(FUN)
    if (!nofun && !is.function(FUN)) 
        stop("'FUN' must be a function or NULL")
    len.r <- length(r <- if (nofun) x[a] else FUN(x[a], ...))
    count <- as.integer(round(choose(n, m)))
    if (simplify) {
        dim.use <- if (nofun) 
            c(m, count)
        else {
            d <- dim(r)
            if (length(d) > 1L) 
                c(d, count)
            else if (len.r > 1L) 
                c(len.r, count)
            else c(d, count)
        }
    }
    if (simplify) 
        out <- matrix(r, nrow = len.r, ncol = count)
    else {
        out <- vector("list", count)
        out[[1L]] <- r
    }
    if (m > 0) {
        i <- 2L
        nmmp1 <- n - m + 1L
        while (a[1L] != nmmp1) {
            if (e < n - h) {
                h <- 1L
                e <- a[m]
                j <- 1L
            }
            else {
                e <- a[m - h]
                h <- h + 1L
                j <- 1L:h
            }
            a[m - h + j] <- e + j
            r <- if (nofun) 
                x[a]
            else FUN(x[a], ...)
            if (simplify) 
                out[, i] <- r
            else out[[i]] <- r
            i <- i + 1L
        }
    }
    if (simplify) {
        if (is.factor(x0)) {
            levels(out) <- levels(x0)
            class(out) <- class(x0)
        }
        dim(out) <- dim.use
    }
    out
}

################################################################
# GET WINDOWS
#############
GetWindows <- function(y, alpha) {
    lenY <- length(y)
    minWin <- ceiling((alpha*lenY))
    allWindows <- nCm(lenY, 2)
    t(allWindows[, allWindows[2,] - allWindows[1,] >= minWin ])
}


################################################################
#  Dependency -- LocReg():
#############
# LocReg() is a wrapper for a function to fit a simple linear regression using either
# weighted or unweighted least squares. The function uses matrix notation to minimize
# memory usage relative to fitting any of the Base linear model functions in R. This
# function will consist of the 'guts' of the Fit Block... I imagine that much of the
# flexibility in the FindLocLin() function will come from minor modifications to this
# bit of code.

LocReg  <-  function(wins, xall, yall, ..., weights=TRUE, verbose=TRUE) {
#    if(verbose) {
#        cat(h, ' ', weights, '\n')
#    }
    #  Grab data window for local regression  #
browser()
    x <- xall[wins[1]:wins[2]]
    y <- yall[wins[1]:wins[2]]
    # Design Matrix #
    X  <-  matrix(cbind(1, x), ncol=2)
    if(weights) { # Use Weighted Least Squares Regression #
        w           <-  TriCube(x=x, m=mean(x), h=(x[wins[2]]-mean(x)))
        bHat        <-  (solve(t(X) %*% diag(w) %*% X)) %*% t(X) %*% diag(w) %*% y
        yHat        <-  X %*% bHat
        sigmaHatUb  <-  sum((y - (X %*% bHat))^2)/(length(y) - 2)
        stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
    } else { # Use Ordinary Least Squares Regression #
        bHat        <-  (solve(t(X) %*% X)) %*% t(X) %*% y
        yHat        <-  X %*% bHat
        sigmaHatUb  <-  sum((y - (X %*% bHat))^2) / (length(y) - 2)
        stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
    }
    data.frame(
        Lbound   = wins[1],
        Rbound   = wins[2],
        r2       =  1 - ((sum((y - yHat)^2)) / (sum((y - mean(y))^2))),
        alph     =  length(wins[1]:wins[2])/length(yall),
        b0       =  bHat[1],
        b1       =  bHat[2],
        skew     =  Skew(x = stdResid), stringsAsFactors=FALSE)
}


################################################################
FindLocLin  <-  function(yall, xall, alpha, plots=TRUE, ...) {
    #  Get windows # 
    windows <- GetWindows(y = yall, alpha) 
    #  Fit Local Regressions  #
    res  <-  apply(windows, 1, LocReg, xall=xall, yall=yall)
    res  <-  do.call(rbind.data.frame, res)
    #  Calculate combined metric (L) for linearity & fit  #
    res$L  <-  ((min(res$skew) + abs(res$skew))/sd(res$skew)) + ((max(res$r2) - res$r2)/sd(res$r2))    
    res    <-  res[with(res, order(L)), ][1:25,]
    #  Plots to accompany best results  #
    if(plots) {        
        toPdf(outputPlot(res, xall, yall), 'testPlots.pdf', height=15, width=15)
        toPdf(outputHist(res), 'testHist.pdf', height=7, width=7)
    }    
    res   
}


xtest <- c(1:20)
ytest <- 0.3*xtest + 1.5 + rnorm(20,sd=0.4)
    wins <- GetWindows(y=ytest, alpha=0.4)
    wins
    res  <-  apply(wins,1, LocReg, xall=xtest, yall=ytest, weights=TRUE)
    res  <-  do.call(rbind.data.frame, res)
    res

apply(wins,1, mean)

x <- x[wins[1,][1]:wins[1,][2]]
y <- y[wins[1,][1]:wins[1,][2]]
