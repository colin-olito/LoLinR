## Main package functions

##' Wrapper - strips NA from input x and y
##'
##' @title Remove NA
##' @param x A numeric vector
##' @param y A numeric vector
##' @return a data.frame with complete.cases
stripNAs  <-  function(x, y) {
    dat  <-  data.frame(x, y)
    dat[complete.cases(dat), ]
}

##' Wrapper - check if vector is numeric
##'
##' @title Wrapper for \code{is.numeric}
##' @param x A numeric vector
##' @return Breaks function and returns error message if x 
##' is non-numeric
checkNumeric  <-  function(x) {
    xNumeric  <-  is.numeric(x)
    if(!xNumeric)
        stop('x must be numeric') 
}

##' Wrapper - Checks that lengths of inputs x and y are equal
##'
##' @title Check Equal
##' @param x A numeric vector
##' @param y A numeric vector
##' @return logical
checkEqualLength  <-  function(x, y) {
    equalLength  <-  length(x) == length(y)
    if(!equalLength)
        stop('x and y must be of equal length')
}

##' Calculate the percentile values of a vector x
##'
##' @title Calculate the percentile values of a vector x
##' @param x A numeric vector
##' @return A numeric vector of percentiles truncated between 0 and 1
##' @export
pcRank  <-  function(x) {
    checkNumeric(x)
    percentiles  <-  trunc(rank(x, na.last=NA)) / sum(!is.na(x))
    allUnique    <-  length(percentiles) == length(unique(percentiles))
    if(!allUnique)
        warning('input/output have ties')
    percentiles
}

##' Sample skewness
##'
##' @title Sample skewness (Fisher-Pearson Standardized Third Moment Coefficient)
##' @param x A numeric vector
##' @details This function is a dependency for \code{findLocLin}
##' where it is used to calculate the (sample) skewness of standardized residuals.
##' @return A numeric vector of length 1
##' @export
skew  <-  function(x, na.rm=TRUE) {
    checkNumeric(x)
    if(na.rm)
        x  <-  x[!is.na(x)]
    n  <-  length(x)
    (n/((n - 1) * (n - 2))) * sum(((x - mean(x)) / sd(x))^3)
}

################################################################
# Dependency -- breuschGodfrey():
###############



##' Breusch-Godfrey Statistic
##'
##' @title Modified Breusch-Godfrey Statistic ((n*R^2)/n)
##' @param x A numeric vector
##' @param y A numeric vector
##' @param order Order to which residuals are lagged. Defaults to \code{order <- (n - k - 1)},
##' where \code{n} is the number of observations, and \code{k} is the number of parameters
##' in the regression (2 by default). This represents the highest possible order given \code{n}.
##' @param fill Defaults to \code{fill = 0}, used to fill model matrix for lagged residuas in 
##' the auxillary regression.

##' @details NOTE: This function is a (very slightly) modified version 
##' of \code{bgtest()} from the \code{\link{lmtest}} package (available at 
##' \link{https://github.com/cran/lmtest/blob/master/R/bgtest.R}).
##' We have stripped it down to minimal functionality for our purposes.
##' All development credit goes to the authors of \code{\link{lmtest}}.
##'
##' This function is a dependency for \code{rankLocReg} where it is
##' used to calculate the Breusch-Godfrey statistic divided by the number 
##' of ovservations ((n*R^2)/n). For the purposes of \code{rankLocReg}, 
##' only the relative variance explained by the fitted values from the
##' auxillary regression and the residuals of the original regression is of
##' interest. Calculating ((n*R^2)/n) preserves this information, while 
##' avoiding the introduction of strong covariance between bgN and n; an 
##' undesirable behaviour for the linearity metric L. If desired, users can 
##' reproduce a standard Breusch-Godfrey Chi-squared test of significance 
##' by using running \code{qchisq()} with the output \code{bgN} and \code{df}.
##' However, we would recommend using the function \code{bgtest()} from the
##' package \code{lmtest}, as it is specifically designed for this purpose.
##' @return A list with: the standard BG statistic (bg), BG/n (bgN), and d.f. (df)
##' @export
breuschGodfrey  <-  function(y, x, order=FALSE, fill=0) {
    X  <-  matrix(cbind(1, x), ncol=2)
    n  <-  nrow(X) 
    k  <-  ncol(X)
    
    if(order) {
        order  <-  1:order
    } else {
        order  <-  1:(n - k - 1)
    }
    
    m       <-  length(order)
    resids  <-  lm.fit(X, y)$residuals
    Z       <-  sapply(order, function(x) c(rep(fill, length.out=x), resids[1:(n - x)]))
    na      <-  !complete.cases(Z)
    if(any(na)) {
        X       <-  X[!na, , drop=FALSE]
        Z       <-  Z[!na, , drop=FALSE]
        y       <-  y[!na]
        resids  <-  resids[!na]
        n       <-  nrow(X)
    }
    auxfit      <-  lm.fit(x=cbind(X, Z), y=resids)
    bg          <-  n * sum(auxfit$fitted^2) / sum(resids^2)
    bgN         <-  bg / n
    names(bg)   <-  'Breusch-Godfrey Statistic'
    names(bgN)  <-  'Breusch-Godfrey Statistic / n'
    df          <-  m
    names(df)   <-  'df'
    list(
        bg    =  bg,
        bgN   =  bgN,
        df    =  df
    )
}

##' Get all possible windows
##'
##' @title Get all possible windows between specified alpha 
##' and 1 
##' @param x A numeric vector
##' @details This function is a dependency for \code{findLocLin}
##' where it is used to extract all local windows
##' for local regressions. alpha must be higher than 0 and lower or equal to 1. 
##' @return A matrix of vector positions, with starting value on first column and ending value on second column.
##' @export
getWindows  <-  function(x, alpha) {
    checkNumeric(x)
    validAlpha  <-  alpha > 0 & alpha <= 1
    if(!validAlpha)
        stop('alpha must take a value higher than 0 and lower or equal to 1')
    lenX        <-  length(x)
    minWin      <-  floor((alpha * lenX))
    allWindows  <-  combn(lenX, 2)
    t(allWindows[, allWindows[2,] - allWindows[1,] >= minWin ])
}

################################################################
#  Dependency -- locReg():
#############
# locReg() is a wrapper for a function to fit a simple linear regression using either
# weighted or unweighted least squares. The function uses matrix notation to minimize
# memory usage relative to fitting any of the Base linear model functions in R. This
# function will consist of the 'guts' of the Fit Block... I imagine that much of the
# flexibility in the findLocLin() function will come from minor modifications to this
# bit of code.

locReg  <-  function(wins, xall, yall, resids=FALSE) {
    # check equal lengths of x and y
    checkEqualLength(xall, yall)

    #  grab data window for local regression
    x  <-  xall[wins[1]:wins[2]]
    y  <-  yall[wins[1]:wins[2]]
    
    # design matrix
    X  <-  matrix(cbind(1, x), ncol=2)
    
    # OLS regression
    lmFit    <-  lm.fit(x=X, y=y)
    bHat     <-  coefficients(lmFit)
    vc       <-  chol2inv(lmFit$qr$qr) * sum(lmFit$residuals^2) / lmFit$df.residual
    b1CI     <-  bHat[2] + qt(c(0.025, 0.975), df=(length(x) - 2)) * sqrt(diag(vc))[2]
    
    # breuschGodfrey test for serial correlation
	BGtest  <-  breuschGodfrey(y, x, order=FALSE)
    bgN     <-  as.numeric(BGtest$bgN)
    
    out     <-  data.frame(
                           Lbound   =  wins[1],
                           Rbound   =  wins[2],
                           alph     =  length(wins[1]:wins[2]) / length(yall),
                           b0       =  bHat[1],
                           b1       =  bHat[2],
                           b1LoCI   =  b1CI[1],
                           b1UpCI   =  b1CI[2],
                           skew     =  skew(x=(lmFit$residuals / lmFit$df.residual)),
                           bgN      =  bgN
                )
    if(resids) {
        sigmaHatUb     <-  sum((lmFit$residuals)^2) / (lmFit$df.residual)
        out$residuals  <-  lmFit$residuals / (sqrt(sigmaHatUb))
    }
    out
}

################################################################
#  THE MAIN WRAPPER FUNCTION -- findLocLin():
#############
rankLocLin  <-  function(xall, yall, alpha, method=c('ns', 'eq', 'pc'), plots=TRUE, verbose=TRUE) {
    if(is.unsorted(xall))
        warning("Dataset must be ordered by xall")
        
    # make sure that all NAs are dealt with
    dat   <-  stripNAs(xall, yall)
    xall  <-  dat$x
    yall  <-  dat$y

    #  get all possible windows
    wins  <-  getWindows(y=yall, alpha)
    
    #  fit local regressions
    allRegs   <-  apply(wins, 1, locReg, xall=xall, yall=yall)
    allRegs   <-  do.call(rbind.data.frame, allRegs)
    
    #  calculate combined metric (L) for linearity & fit
    allRegs$ciRange  <-  allRegs$b1UpCI - allRegs$b1LoCI
    allRegs          <-  allRegs[, c('Lbound', 'Rbound', 'alph', 'b0', 'b1', 'b1LoCI', 'b1UpCI', 'ciRange', 'skew', 'bgN')]
    allRegs$L        <-  ((min(abs(allRegs$skew)) + abs(allRegs$skew)) / sd(allRegs$skew)) + ((allRegs$bgN - min(allRegs$bgN)) / sd(allRegs$bgN)) + ((allRegs$ciRange - min(allRegs$ciRange)) / sd(allRegs$ciRange))
    allRegs$Leq      <-  (((min(abs(allRegs$skew)) + abs(allRegs$skew)) / sd(allRegs$skew)) / (max(((min(abs(allRegs$skew)) + abs(allRegs$skew)) / sd(allRegs$skew))))) + (((allRegs$bgN - min(allRegs$bgN)) / sd(allRegs$bgN)) / (max(((allRegs$bgN - min(allRegs$bgN)) / sd(allRegs$bgN))))) + (((allRegs$ciRange - min(allRegs$ciRange)) / sd(allRegs$ciRange)) / (max(((allRegs$ciRange - min(allRegs$ciRange)) / sd(allRegs$ciRange)))))
    allRegs$Lpc      <-  ((pcRank(abs(allRegs$skew))) + (pcRank((allRegs$bgN - min(allRegs$bgN)))) + (pcRank((allRegs$ciRange)))) / 3
    
    # choose weighting scheme for linearity metric L
    switch(match.arg(method),
        'ns' = {
            allRegs   <-  allRegs[with(allRegs, order(L)), ]
        },
        'eq' = {
            allRegs   <-  allRegs[with(allRegs, order(Leq)), ]
        },
        'pc' = {
            allRegs   <-  allRegs[with(allRegs, order(Lpc)), ]
        }
    )
    
    #  plots to accompany best local regression
    if(plots) {        
        dev.new(height=15, width=15)
        outputPlot(allRegs, xall, yall)
    }
    
    nFits  <-  nrow(allRegs)
    if(verbose)
        cat(sprintf('rankLocLin fitted %d local regressions', nFits), '\n')
    
    list(
        'nFits'    =  nFits,
        'allRegs'  =  allRegs
    )
}

####################
# PLOTTING FUNCTIONS
####################
outputPlot  <-  function(allRegs, x, y) {
    col1  <-  adjustcolor('#1B6889', alpha=0.5)
    par(mfrow=c(5,5))
    for(i in seq_len(nrow(allRegs))) {
        # Subset Data #
        ytemp  <-  y[c(allRegs$Lbound[i]:allRegs$Rbound[i])]
        xtemp  <-  x[c(allRegs$Lbound[i]:allRegs$Rbound[i])]
        # Plot #
        plot(y ~ x, pch=21, col='grey80', main=i)
        points(ytemp ~ xtemp, pch=21, bg=col1, ask=TRUE)
        abline(coef=c(allRegs$b0[i],allRegs$b1[i]), col=2)
    }
}

outputHist  <-  function(allRegs) {
    hist(allRegs$b1, breaks=25)
}

###################################################
#  plotBest():
#
#  Followup function for use with findLocLin(). Generates residual
#    plots and stand alone scatterplot for a local linear regression
#    chosen by the user from the findLocLin() output data frame.
#
#    Takes 5 arguments:
#      -- res:  findLocLin() output data frame
#      -- yall:  Same yall input data as for findLocLin()
#      -- xall:  Same xall input data as for findLocLin()
#      -- best:  row number corresponding to the local linear
#                 regression from findLocLIn output the user
#                 wishes to plot/inspect. Usually the 'best'
#                 linear regression.
#
###################################################

plotBest <- function(allRegs, xall, yall, best=1) {
    #  Recover data window for chosen local regression model
    bestwin  <-  c(allRegs$allRegs$Lbound[best], allRegs$allRegs$Rbound[best])
    y        <-  yall[bestwin[1]:bestwin[2]]
    x        <-  xall[bestwin[1]:bestwin[2]]
    
    #  fit block
    LocFit  <-  locReg(bestwin, xall, yall, resids=TRUE)
    b1      <-  LocFit$bHat[2]

    #  residual plots
    dev.new()
    
    # standardized residuals ~ x
    par(mfrow=c(2,2))
    plot(LocFit$stdResid ~ x, xlab='x', ylab='y', main='Std. Residuals ~ x')
    abline(h=0, col=1, lwd=2)
    abline(h=c(-2, 2), lty=2)
    lf1  <-  loess(LocFit$stdResid ~ x)
    points(x, lf1$fitted, type='l', col=2, lwd=2)
    
    # standardized residuals ~ fitted values
    plot(LocFit$stdResid ~ LocFit$yHat, xlab='Fitted Values',ylab='Standardized Residuals',main='Std. Residuals ~ Fitted Values')
    abline(h=0, col=1, lwd=2)
    abline(h=c(-2, 2), lty=2)
    lf2  <-  loess(LocFit$stdResid ~ LocFit$yHat)
    points(LocFit$yHat, lf2$fitted, type='l', col=2, lwd=2)
    
    # qqnorm plot of standardized residuals
    qqnorm(LocFit$stdResid, main='QQNorm plot of Std. Residuals')
    qqline(LocFit$stdResid, col=2)
    
    # histogram of standardized residuals
    hist(LocFit$stdResid, xlab='Standardized Residuals', ylab='Density', breaks=20, main='Density Plot of Std. Residuals')
    
    #  overall regression plot
    dev.new()
    col1  <-  adjustcolor('#1B6889', alpha=0.5)
    plot(yall ~ xall, pch=21, col='grey80', ask=TRUE, main=expression(paste('Best Local Regression: ', beta[1],' = ', b1)))
    points(y ~ x, pch=21, bg=col1, ask=TRUE)
    abline(coef=c(LocFit$bHat[1], LocFit$bHat[2]), col=1)
}

plotBeta1 <- function(results) {
    c1  <-  adjustcolor('#A67A01', alpha=0.75)
    c2  <-  adjustcolor('#6D65FA', alpha=0.75)
    c3  <-  adjustcolor('#B6084E', alpha=0.75)
    dev.new()
    par(mfrow=c(1, 1))
    plot(density(results$res$b1), lwd=4, col=col1, xlab=expression(paste(beta[1])), main=expression(paste('Distribution of ', beta[1])), cex.main=2)
    abline(v=results$res$b1[results$res$L == min(results$res$L)], col=c1, lty=1, lwd=4)
    abline(v=results$res$b1[results$res$Leq == min(results$res$Leq)], col=c2, lty=2, lwd=4)
    abline(v=results$res$b1[results$res$Lpc == min(results$res$Lpc)], col=c3, lty=3, lwd=4)
    legend(
          x       =  min(results$res$b1) + (0.8 * (abs(range(results$res$b1)[2] - range(results$res$b1)[1]))),
          y       =  0.95 * max(density(results$res$b1)$y),
          legend  =  c(expression(paste(italic(L))),
                      expression(paste(italic(L[eq]))),
                      expression(paste(italic(L['%'])))),
          lwd     =  4,
          lty     =  c(1, 2, 3),
          col     =  c(c1, c2, c3),
          cex     =  1
    )
}
