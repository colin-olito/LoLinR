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
# Dependency -- nCm():
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
        BGtest      <-  bgtest(y ~ x, order=(length(x)-3))
        bg.df       <-  as.numeric(BGtest$statistic) / (BGtest$parameter)
    } else { # Use Ordinary Least Squares Regression #
        bHat        <-  (solve(t(X) %*% X)) %*% t(X) %*% y
        yHat        <-  X %*% bHat
        sigmaHatUb  <-  sum((y - yHat)^2) / (length(y) - 2)
        vcov        <-  sigmaHatUb * (solve(t(X) %*% X))
        b1.CI       <-  bHat[2,] + qt(c(0.025,0.975),df=(length(x)-2))*sqrt(diag(vcov))[2]
        stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
        BGtest      <-  bgtest(y ~ x, order=(length(x)-3))
        bg.df       <-  as.numeric(BGtest$statistic) / (BGtest$parameter)
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
        bg.df    =  bg.df
    )
}


################################################################
#  THE MAIN WRAPPER FUNCTION -- FindLocLin():
#############
FindLocLin  <-  function(yall, xall, alpha, ref.b1=FALSE, method = c("ns","sr","pc"),
                         plots=TRUE, plot.name="testPlots.pdf.", all=FALSE, ...) {
    #  Get windows # 
    wins  <-  GetWindows(y = yall, alpha)
    #  Fit Local Regressions  #
    res   <-  apply(wins, 1, LocReg, xall=xall, yall=yall)
    res   <-  do.call(rbind.data.frame, res)
    #  Calculate combined metric (L) for linearity & fit  #
    res$CI.range <- res$b1.CI.hi - res$b1.CI.lo
    res  <-  res[, c('weights', 'Lbound', 'Rbound', 'alph', 'b0', 'b1', 'b1.CI.lo', 'b1.CI.hi', 'CI.range', 'skew', 'bg.df')]
    res$L     <-  ((min(abs(res$skew)) + abs(res$skew))     / sd(res$skew)) +
                  ((res$bg.df          - min(res$bg.df))    / sd(res$bg.df)) +
                  ((res$CI.range       - min(res$CI.range)) / sd(res$CI.range))
    res$L.sr  <-  (((min(abs(res$skew)) + abs(res$skew))     / sd(res$skew))     / (max(((min(abs(res$skew)) + abs(res$skew))     / sd(res$skew))))) +
                  (((res$bg.df          - min(res$bg.df))    / sd(res$bg.df))    / (max(((res$bg.df          - min(res$bg.df))    / sd(res$bg.df))))) +
                  (((res$CI.range       - min(res$CI.range)) / sd(res$CI.range)) / (max(((res$CI.range       - min(res$CI.range)) / sd(res$CI.range)))))
    res$L.pc  <-  ((pc.rank(abs(res$skew))) +
                  (pc.rank((res$bg.df - min(res$bg.df)))) +
                  (pc.rank((res$CI.range)))) / 3
    switch(match.arg(method),
           "ns" = {
               res   <- res[with(res, order(L)), ]
          },
           "sr" = {
               res   <- res[with(res, order(L.sr)), ]
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
        acorr       <-  as.vector(acf(stdResid, plot=FALSE)$acf[-1])
        Xac         <-  matrix(cbind(1, seq(0,1,len=length(acorr))), ncol=2)
        b1.ac       <-  (solve(t(Xac) %*% Xac)) %*% t(Xac) %*% acorr 
        BGtest      <-  bgtest(y ~ x, order=(length(x)-3))
        bg          <-  as.numeric(BGtest$statistic) / as.numeric(BGtest$parameter)
    } else { # Use Ordinary Least Squares Regression #
        bHat        <-  (solve(t(X) %*% X)) %*% t(X) %*% y
        yHat        <-  X %*% bHat
        sigmaHatUb  <-  sum((y - yHat)^2) / (length(y) - 2)
        vcov        <-  sigmaHatUb * (solve(t(X) %*% X))
        b1.CI       <-  bHat[2,] + qt(c(0.025,0.975),df=(length(x)-2))*sqrt(diag(vcov))[2]
        stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
        acorr       <-  as.vector(acf(stdResid, plot=FALSE)$acf[-1])
        Xac         <-  matrix(cbind(1, seq(0,1,len=length(acorr))), ncol=2)
        b1.ac       <-  (solve(t(Xac) %*% Xac)) %*% t(Xac) %*% acorr 
        BGtest      <-  bgtest(y ~ x, order=(length(x)-3))
        bg          <-  as.numeric(BGtest$statistic) / as.numeric(BGtest$parameter)
    }
    list(
        'BestWindow' =  bestwin,
        'bHat'       =  bHat,
        'yHat'       =  yHat,
        'sigmaHatUb' =  sigmaHatUb,
        'stdResid'   =  stdResid,
        'b1.CI'      =  b1.CI,
        'b1.ac'      =  b1.ac,
        'bg'         =  bg,
        'skew'       =  Skew(x = stdResid)
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
