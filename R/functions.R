################################################################
#  Functions and Dependencies for the primary function FindLocLin(), used to find
#  the 'most linear' local regressions from a given dataset.
#
#  Run these functions

library(lmtest)

################################################################
#  Dependency -- Function TriCube():
#
# TriCube() is a function to calculate Tri-Cube weights for a given vector of values
# for an independent variable. This function is a dependency for FindLocLin(), and is
# used to calculate the weights for Weighted Least Squares Regression.
TriCube <-  function(x, m, h) {
    z  <-  abs(x - m) / h
    ifelse(z < 1, (1 - z^3)^3, 0)
}

################################################################
#  Dependency -- perc.rank():
#
# perc.rank() is a function to calculate the percentile values of a vector x.
pc.rank <- function(x) trunc(rank(x,na.last = NA))/sum(!is.na(x))



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
# Dependency -- combn():
###############

# NOTE: This code is copied directly from the combn() function
#        from the R {utils} package. Not sure if this matters
#        as a dependency (i.e. whether we should strip this
#        function down to just what we need)... need to ask
#        Diego...

combn <- function (x, m, FUN = NULL, simplify = TRUE, ...) 
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
# Dependency -- simpleReg():
###############

# NOTE: simple function to run a linear regression... for use 
#       with BG()

simpleReg <- function(X, y) {
    if(!is.matrix(X))
        stop("X must be a model matrix")
    if(nrow(X) != length(y))
         stop("X must be a model matrix")
#    bHat        <-  (solve(t(X) %*% X)) %*% t(X) %*% y
    qr <- qr(X)$qr
    bHat <- chol2inv(qr) %*% t(X) %*% y
    yHat        <-  X %*% bHat
    sigmaHatUb  <-  sum((y - yHat)^2) / (length(y) - 2)
#    vcov        <-  sigmaHatUb * (solve(t(X) %*% X))
    vcov        <-  sigmaHatUb * chol2inv(qr)
    stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
     list(
        'bHat'       =  bHat,
        'yHat'       =  yHat,
        'sigmaHatUb' =  sigmaHatUb,
        'vcov'       =  vcov,
        'stdResid'   =  stdResid
        )
}

################################################################
# Dependency -- BG():
###############

# NOTE: This code is (very slightly) modified from bgtest.R
#        from the lmtest github repo. I have stripped it down
#       to minimal functionality for our purposes

BG <- function(y, x, order = FALSE, fill=0) {
    X  <-  matrix(cbind(1, x), ncol=2)
    n <- nrow(X) 
    k <- ncol(X)
    if(order)
        order <- 1:order
    else {
        order <- 1:(n-k-1)
    }
    m <- length(order)
    resids <- simpleReg(X, y)$stdResid
    Z <- sapply(order, function(x) c(rep(fill, length.out = x), resids[1:(n-x)]))
        if(any(na <- !complete.cases(Z))) {
            X <- X[!na, , drop = FALSE]
            Z <- Z[!na, , drop = FALSE]
            y <- y[!na]
            resids <- resids[!na]
            n <- nrow(X)
        }
    auxfit <- simpleReg(X=cbind(X, Z), resids)$yHat

    bg   <- n * sum(auxfit^2)/sum(resids^2)
    bg.n <- bg / n
    names(bg) <- "Breusch-Godfrey Statistic"
    names(bg.n) <- "Breusch-Godfrey Statistic / n"
    df <- m
    names(df) <- "df"
    list(
        bg    =  bg,
        bg.n  =  bg.n,
        df    =  df
        )
}


makedata <- function(n, b0=0.1, b1=2, sd.x = 1, sd.e = 0.1) {
    x <- rnorm(n, mean=10,sd=sd.x)
    e <- rnorm(n, sd=sd.e)
    y <- b0 + b1*x + e
    list(
        'y'          = y,
        'x'          = x
    )
}


## troubleshooting BG()
# d <- makedata(50)
# str(qr(cbind(1,d$x)))
# bgtest(d$y ~ d$x, order=(length(d$x)-3))
# BG(d$y, d$x)
# simpleReg(cbind(1,d$x), d$y)$bHat

################################################################
# GET WINDOWS
#############
GetWindows  <-  function(y, alpha) {
    lenY        <-  length(y)
    minWin      <-  floor((alpha*lenY))
    allWindows  <-  combn(lenY, 2)
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

LocReg  <-  function(wins, xall, yall, ..., weights=TRUE) {
    #  Grab data window for local regression  #
    x <- xall[wins[1]:wins[2]]
    y <- yall[wins[1]:wins[2]]
    # Design Matrix #
    X  <-  matrix(cbind(1, x), ncol=2)
    if(weights) { # Use Weighted Least Squares Regression #
        w           <-  TriCube(x=x, m=mean(x), h=(x[length(x)] - x[1])/2)
        bHat        <-  (solve(t(X) %*% diag(w) %*% X)) %*% t(X) %*% diag(w) %*% y
        yHat        <-  X %*% bHat
        sigmaHatUb  <-  sum((y - yHat)^2)/(length(y) - 2)
        vcov        <-  sigmaHatUb * (solve(t(X) %*% diag(w) %*% X))
        b1.CI       <-  bHat[2,] + qt(c(0.025,0.975),df=(length(x)-2))*sqrt(diag(vcov))[2]
        stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
        BGtest2      <-  BG(y, x, order = FALSE)
        bg.n2        <-  as.numeric(BGtest2$bg.n)
        BGtest      <-  bgtest(y ~ x, order = (nrow(X) - 3))
        bg.n        <-  as.numeric(BGtest$statistic) / nrow(X)
    } else { # Use Ordinary Least Squares Regression #
        bHat        <-  (solve(t(X) %*% X)) %*% t(X) %*% y
        yHat        <-  X %*% bHat
        sigmaHatUb  <-  sum((y - yHat)^2) / (length(y) - 2)
        vcov        <-  sigmaHatUb * (solve(t(X) %*% X))
        b1.CI       <-  bHat[2,] + qt(c(0.025,0.975),df=(length(x)-2))*sqrt(diag(vcov))[2]
        stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
        BGtest2      <-  BG(y, x, order = FALSE)
        bg.n2        <-  as.numeric(BGtest2$bg.n)
        BGtest      <-  bgtest(y ~ x, order = (nrow(X) - 3))
        bg.n        <-  as.numeric(BGtest$statistics) / nrow(X)
    }
    data.frame(
        weights  =  weights,
        Lbound   =  wins[1],
        Rbound   =  wins[2],
        alph     =  length(wins[1]:wins[2])/length(yall),
        b0       =  bHat[1],
        b1       =  bHat[2],
        b1.CI.lo =  b1.CI[1],
        b1.CI.hi =  b1.CI[2],
        skew     =  Skew(x = stdResid),
        bg.n     =  bg.n,
        bg.n2    =  bg.n2
     )
}

################################################################
#  THE MAIN WRAPPER FUNCTION -- FindLocLin():
#############
FindLocLin  <-  function(yall, xall, alpha, ref.b1=FALSE, method = c("ns","eq","pc"),
                         plots=TRUE, plot.name="testPlots.pdf.", all=FALSE, ...) {
    #  Get windows # 
    wins  <-  GetWindows(y = yall, alpha)
    #  Fit Local Regressions  #
    res   <-  apply(wins, 1, LocReg, xall=xall, yall=yall)
    res   <-  do.call(rbind.data.frame, res)
    #  Calculate combined metric (L) for linearity & fit  #
    res$CI.range <- res$b1.CI.hi - res$b1.CI.lo
    res  <-  res[, c('weights', 'Lbound', 'Rbound', 'alph', 'b0', 'b1', 'b1.CI.lo', 'b1.CI.hi', 'CI.range', 'skew', 'bg.n', 'bg.n2' )]
    res$L     <-  ((min(abs(res$skew)) + abs(res$skew))     / sd(res$skew)) +
                  ((res$bg.n           - min(res$bg.n) )    / sd(res$bg.n) ) +
                  ((res$CI.range       - min(res$CI.range)) / sd(res$CI.range))
    res$L.eq  <-  (((min(abs(res$skew)) + abs(res$skew))     / sd(res$skew))     / (max(((min(abs(res$skew)) + abs(res$skew))     / sd(res$skew))))) +
                  (((res$bg.n           - min(res$bg.n) )    / sd(res$bg.n) )    / (max(((res$bg.n           - min(res$bg.n) )    / sd(res$bg.n) )))) +
                  (((res$CI.range       - min(res$CI.range)) / sd(res$CI.range)) / (max(((res$CI.range       - min(res$CI.range)) / sd(res$CI.range)))))
    res$L.pc  <-  ((pc.rank(abs(res$skew))) +
                  (pc.rank((res$bg.n  - min(res$bg.n)))) +
                  (pc.rank((res$CI.range)))) / 3
    switch(match.arg(method),
           "ns" = {
               res   <- res[with(res, order(L)), ]
          },
           "eq" = {
               res   <- res[with(res, order(L.eq)), ]
           },
           "pc" = {
               res   <- res[with(res, order(L.pc)), ]
           }
          )
          numericB1  <-  is.numeric(ref.b1)
          if(numericB1) {
              res   <- res[with(res, order(abs(ref.b1 - res$b1))), ]
          } 
          if(!all) {
              res <- res[1:25,]
          }
    #  Plots to accompany best results  #
    if(plots) {        
        toPdf(outputPlot(res, xall, yall), filename=plot.name, height=15, width=15)
    }    
    nFits <- print(nrow(res))
    list(
#        'alpha'   =  alpha,
#        'weights'   =  weights,
        'nFits'   =  nFits,
        'res'     =  res
        )
}

####################
# PLOTTING FUNCTIONS
####################
toDev  <-  function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
  dev(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

toPdf  <-  function(expr, filename, ...) {
  toDev(expr, pdf, filename, ...)
}

outputPlot  <-  function(resultsTable, x, y) {
    col1 <- adjustcolor('#1B6889', alpha=0.5)
    par(mfrow=c(5,5))
    for(i in 1:nrow(resultsTable)) {
# Subset Data #
        ytemp <- y[c(resultsTable$Lbound[i]:resultsTable$Rbound[i])]
        xtemp <- x[c(resultsTable$Lbound[i]:resultsTable$Rbound[i])]
            # Plot #
        plot(y ~ x, pch=21, col='grey80', main=i)
        points(ytemp ~ xtemp, pch=21, bg=col1, ask=TRUE)
        abline(coef=c(resultsTable$b0[i],resultsTable$b1[i]), col=2)
    }
}

outputHist  <-  function(resultsTable) {
        hist(resultsTable$b1, breaks=25)
}






##############################
#  Dependency -- BestLocReg():
##############################
# for use with PlotBest() -- slight modification of LocReg, that
#  accepts the best window chosen by the user instead of subsetting
#  the dataset. Also uses a list() to output objects of different
#  length (e.g. bHat and stdResid).

BestLocReg  <-  function(bestwin, y, x, weights, ..., verbose=TRUE) {
    # Design Matrix #
    X  <-  matrix(cbind(1, x), ncol=2)
    if(weights) { # Use Weighted Least Squares Regression #
        w           <-  TriCube(x=x, m=mean(x), h=(x[length(x)] - x[1])/2)
        bHat        <-  (solve(t(X) %*% diag(w) %*% X)) %*% t(X) %*% diag(w) %*% y
        yHat        <-  X %*% bHat
        sigmaHatUb  <-  sum((y - yHat)^2)/(length(y) - 2)
        vcov        <-  sigmaHatUb * (solve(t(X) %*% diag(w) %*% X))
        b1.CI       <-  bHat[2,] + qt(c(0.025,0.975),df=(length(x)-2))*sqrt(diag(vcov))[2]
        stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
        BGtest      <-  BG(y, x, order = FALSE)
        bg.n        <-  as.numeric(BGtest$bg)
    } else { # Use Ordinary Least Squares Regression #
        bHat        <-  (solve(t(X) %*% X)) %*% t(X) %*% y
        yHat        <-  X %*% bHat
        sigmaHatUb  <-  sum((y - yHat)^2) / (length(y) - 2)
        vcov        <-  sigmaHatUb * (solve(t(X) %*% X))
        b1.CI       <-  bHat[2,] + qt(c(0.025,0.975),df=(length(x)-2))*sqrt(diag(vcov))[2]
        stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
        BGtest      <-  BG(y, x, order = FALSE)
        bg.n        <-  as.numeric(BGtest$bg)
    }
    list(
        'BestWindow' =  bestwin,
        'bHat'       =  bHat,
        'yHat'       =  yHat,
        'sigmaHatUb' =  sigmaHatUb,
        'stdResid'   =  stdResid,
        'skew'       =  Skew(x = stdResid),
        'b1.CI'      =  b1.CI,
        'bg'         =  bg.n
    )
}



###################################################
#  PlotBest():
#
#  Followup function for use with FindLocLin(). Generates residual
#    plots and stand alone scatterplot for a local linear regression
#    chosen by the user from the FindLocLin() output data frame.
#
#    Takes 5 arguments:
#      -- res:  FindLocLin() output data frame
#      -- yall:  Same yall input data as for FindLocLin()
#      -- xall:  Same xall input data as for FindLocLin()
#      -- best:  row number corresponding to the local linear
#                 regression from FindLocLIn output the user
#                 wishes to plot/inspect. Usually the 'best'
#                 linear regression.
#
###################################################

PlotBest <- function(res, yall, xall, best=1) {
    #  Recover data window for chosen local regression model  #
    weight=res$res$weights[1]
    bestwin <- c(res$res$Lbound[best],res$res$Rbound[best])
    y       <- yall[bestwin[1]:bestwin[2]]
    x       <- xall[bestwin[1]:bestwin[2]]
    #  Fit block  #
    LocFit <- BestLocReg(bestwin, y=y, x=x, weight)
    b1 <- as.numeric(LocFit$bHat[2,1])
    
    #  Residual Plots  #
    dev.new()
    # Standardized Residuals ~ x
    par(mfrow=c(2,2))
    plot(LocFit$stdResid ~ x,
         xlab="x", ylab="y", main="Std. Residuals ~ x")
    abline(h=0,col=1, lwd=2)
    abline(h=c(-2,2), lty=2)
    lf1 <- loess(LocFit$stdResid ~ x)
    points(x, lf1$fitted, type='l', col=2, lwd=2)
    # Standardized Residuals ~ Fitted Values
    plot(LocFit$stdResid ~ LocFit$yHat,
         xlab="Fitted Values",ylab="Standardized Residuals",main="Std. Residuals ~ Fitted Values")
    abline(h=0,col=1, lwd=2)
    abline(h=c(-2,2), lty=2)
    lf2 <- loess(LocFit$stdResid ~ LocFit$yHat)
    points(LocFit$yHat, lf2$fitted, type='l', col=2, lwd=2)
    # QQNorm Plot of Standardized Residuals #
    qqnorm(LocFit$stdResid, main="QQNorm plot of Std. Residuals")
    qqline(LocFit$stdResid,col=2)
    # Histogram of Standardized Residuals #
    hist(LocFit$stdResid, xlab="Standardized Residuals", ylab="Density",breaks=20, main="Density Plot of Std. Residuals")
    
    #  Overall Regression Plot  #
    dev.new()
    col1 <- adjustcolor('#1B6889', alpha=0.5)
    plot(yall ~ xall, pch=21, col='grey80', ask=TRUE,
         main=expression(paste("Best Local Regression: ", beta[1]," = ", b1)))
    points(y ~ x, pch=21, bg=col1, ask=TRUE)
    abline(coef=c(LocFit$bHat[1],LocFit$bHat[2]), col=1)    
}


PlotBeta1 <- function(results){
    c1 <- adjustcolor('#A67A01', alpha=0.75)
    c2 <- adjustcolor('#6D65FA', alpha=0.75)
    c3 <- adjustcolor('#B6084E', alpha=0.75)
#browser()
    dev.new()
    par(mfrow=c(1,1))
    #pdf(file="beta1.pdf", family="CM Roman", width=6, height=6)
    plot(density(results$res$b1), lwd=4, col=col1, xlab=expression(paste(beta[1])),
        main=expression(paste("Distribution of ", beta[1])), cex.main=2)
    abline(v=results$res$b1[results$res$L == min(results$res$L)], col=c1, lty=1, lwd=4)
    abline(v=results$res$b1[results$res$L.eq == min(results$res$L.eq)], col=c2, lty=2, lwd=4)
    abline(v=results$res$b1[results$res$L.pc == min(results$res$L.pc)], col=c3, lty=3, lwd=4)
    legend(x = min(res$b1) + (0.8 * (abs(range(results$res$b1)[2] - range(results$res$b1)[1]))),
        y = (0.95*max(density(res$b1)$y)),
        legend = c(expression(paste(italic(L))),
                   expression(paste(italic(L[eq]))),
                   expression(paste(italic(L["%"])))),
        lwd = 4,
        lty = c(1,2,3),
        col = c(c1, c2, c3),
        cex=1)
    #dev.off()
    #embed_fonts("beta1.pdf", outfile = "beta1.pdf")
}