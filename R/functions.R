################################################################
#  Functions and Dependencies for the primary function findLocLin(), used to find
#  the 'most linear' local regressions from a given dataset.
#
#  Run these functions

library(lmtest)

################################################################
#  Dependency -- Function triCube():
#
# triCube() is a function to calculate Tri-Cube weights for a given vector of values
# for an independent variable. This function is a dependency for findLocLin(), and is
# used to calculate the weights for Weighted Least Squares Regression.
triCube <-  function(x, m, h) {
    z  <-  abs(x - m) / h
    ifelse(z < 1, (1 - z^3)^3, 0)
}

################################################################
#  Dependency -- perc.rank():
#
# pcRank() is a function to calculate the percentile values of a vector x.
pcRank <- function(x) {
    trunc(rank(x, na.last=NA)) / sum(!is.na(x))
}

################################################################
#  Dependency -- Function skew():
#
# skew() is a function to calculate the Fisher-Pearson Standardized Third Moment
# Coefficient for a given vector of numbers. This function is a dependency for
# findLocLin(), and is used to calculate the skewness of standardized residuals.

# Note: Compare to moments package
skew  <-  function(x) {
    n  <-  length(x)
    (n/((n - 1) * (n - 2))) * sum(((x - mean(x)) / sd(x))^3)
}

################################################################
# Dependency -- simpleReg():
###############

# NOTE: simple function to run a linear regression... for use 
#       with breuschGodfrey()
simpleReg <- function(X, y) {
    if(!is.matrix(X)) {
        stop('X must be a model matrix')
    }
    if(nrow(X) != length(y)) {
        stop('X must be a model matrix')
    }
    QR          <-  qr(X, tol=1e-10)$qr
    bHat        <-  chol2inv(QR) %*% t(X) %*% y
#    bHat        <-  solve(t(X) %*% X) %*% t(X) %*% y
    yHat        <-  X %*% bHat
    sigmaHatUb  <-  sum((y - yHat)^2) / (length(y) - 2)
    vcov        <-  sigmaHatUb * chol2inv(QR)
    stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
#    stdResid    <-  (y - yHat)
    list(
        'bHat'       =  bHat,
        'yHat'       =  yHat,
        'sigmaHatUb' =  sigmaHatUb,
        'vcov'       =  vcov,
        'stdResid'   =  stdResid
    )
}


d <- makedata(50)
X <- cbind(1,d$x)
solve(t(X) %*% X) %*% t(X) %*% d$y
chol2inv(qr(X, tol=1e-10)$qr) %*% t(X) %*% d$y
chol2inv(qr(X, tol=1e-10)$qr)
head(solve(t(X) %*% X))

head(qr(X, tol=1e-10)$qr)

lm1 <- lm.fit(X, d$y)
test <-  simpleReg(X, d$y)
test$sigmaHatUb

lm1$resid
as.vector(test$stdResid)
lm1$resid - as.vector(test$stdResid * sqrt(test$sigmaHatUb))

lm1$fitted - as.vector(test$yHat)




as.vector(simpleReg(X=cbind(X, Z), resids)$yHat) - lm.fit(x=cbind(X, Z), resids)$fitted

breuschGodfrey(d$y, d$x)
bgtest(d$y ~ d$x, order=(length(d$y)-3))

################################################################
# Dependency -- breuschGodfrey():
###############


# NOTE: This code is (very slightly) modified from bgtest.R
#        from the lmtest github repo. I have stripped it down
#       to minimal functionality for our purposes
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
    resids  <-  simpleReg(X, y)$stdResid
    Z       <-  sapply(order, function(x) c(rep(fill, length.out=x), resids[1:(n - x)]))
    na      <-  !complete.cases(Z)
    if(any(na)) {
        X       <-  X[!na, , drop=FALSE]
        Z       <-  Z[!na, , drop=FALSE]
        y       <-  y[!na]
        resids  <-  resids[!na]
        n       <-  nrow(X)
    }
    auxfit      <-  simpleReg(X=cbind(X, Z), resids)
    bg          <-  n * sum(auxfit$yHat^2) / sum(resids^2)
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

################################################################
# GET WINDOWS
#############
getWindows  <-  function(y, alpha) {
    lenY        <-  length(y)
    minWin      <-  floor((alpha * lenY))
    allWindows  <-  combn(lenY, 2)
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

locReg  <-  function(wins, xall, yall, weights=TRUE, ...) {
    #  Grab data window for local regression  #
    x  <-  xall[wins[1]:wins[2]]
    y  <-  yall[wins[1]:wins[2]]
    # Design Matrix #
    X  <-  matrix(cbind(1, x), ncol=2)
    if(weights) { # Use Weighted Least Squares Regression #
        w           <-  triCube(x=x, m=mean(x), h=(x[length(x)] - x[1]) / 2)
        bHat        <-  (solve(t(X) %*% diag(w) %*% X)) %*% t(X) %*% diag(w) %*% y
        yHat        <-  X %*% bHat
        sigmaHatUb  <-  sum((y - yHat)^2) / (length(y) - 2)
        vcov        <-  sigmaHatUb * (solve(t(X) %*% diag(w) %*% X))
        b1CI        <-  bHat[2, ] + qt(c(0.025, 0.975), df=(length(x) - 2)) * sqrt(diag(vcov))[2]
        stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
        BGtest2     <-  breuschGodfrey(y, x, order=FALSE)
        bgN2        <-  as.numeric(BGtest2$bgN)
        BGtest      <-  bgtest(y ~ x, order=(nrow(X) - 3))
        bgN         <-  as.numeric(BGtest$statistic) / nrow(X)
    } else { # Use Ordinary Least Squares Regression #
        bHat        <-  (solve(t(X) %*% X)) %*% t(X) %*% y
        yHat        <-  X %*% bHat
        sigmaHatUb  <-  sum((y - yHat)^2) / (length(y) - 2)
        vcov        <-  sigmaHatUb * (solve(t(X) %*% X))
        b1CI        <-  bHat[2, ] + qt(c(0.025, 0.975), df=(length(x) - 2)) * sqrt(diag(vcov))[2]
        stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
        BGtest2     <-  breuschGodfrey(y, x, order=FALSE)
        bgN2        <-  as.numeric(BGtest2$bgN)
        BGtest      <-  bgtest(y ~ x, order=(nrow(X) - 3))
        bgN         <-  as.numeric(BGtest$statistics) / nrow(X)
    }
    data.frame(
        weights  =  weights,
        Lbound   =  wins[1],
        Rbound   =  wins[2],
        alph     =  length(wins[1]:wins[2]) / length(yall),
        b0       =  bHat[1],
        b1       =  bHat[2],
        b1LoCI   =  b1CI[1],
        b1UpCI   =  b1CI[2],
        skew     =  skew(x=stdResid),
        bgN      =  bgN,
        bgN2     =  bgN2
    )
}

################################################################
#  THE MAIN WRAPPER FUNCTION -- findLocLin():
#############
findLocLin  <-  function(yall, xall, alpha, refB1=FALSE, method=c('ns', 'eq', 'pc'),
                         plots=TRUE, plotName='testPlots.pdf.', ...) {
    #  Get windows # 
    wins  <-  getWindows(y=yall, alpha)
    #  Fit Local Regressions  #
    res   <-  apply(wins, 1, locReg, xall=xall, yall=yall)
    res   <-  do.call(rbind.data.frame, res)
    #  Calculate combined metric (L) for linearity & fit  #
    res$ciRange  <-  res$b1UpCI - res$b1LoCI
    res          <-  res[, c('weights', 'Lbound', 'Rbound', 'alph', 'b0', 'b1', 'b1LoCI', 'b1UpCI', 'ciRange', 'skew', 'bgN', 'bgN2')]
    res$L        <-  ((min(abs(res$skew)) + abs(res$skew)) / sd(res$skew)) + ((res$bgN - min(res$bgN)) / sd(res$bgN)) + ((res$ciRange - min(res$ciRange)) / sd(res$ciRange))
    res$Leq      <-  (((min(abs(res$skew)) + abs(res$skew)) / sd(res$skew)) / (max(((min(abs(res$skew)) + abs(res$skew)) / sd(res$skew))))) + (((res$bgN - min(res$bgN)) / sd(res$bgN)) / (max(((res$bgN - min(res$bgN)) / sd(res$bgN))))) + (((res$ciRange - min(res$ciRange)) / sd(res$ciRange)) / (max(((res$ciRange - min(res$ciRange)) / sd(res$ciRange)))))
    res$Lpc      <-  ((pcRank(abs(res$skew))) + (pcRank((res$bgN - min(res$bgN)))) + (pcRank((res$ciRange)))) / 3
    
    switch(match.arg(method),
        'ns' = {
            res   <-  res[with(res, order(L)), ]
        },
        'eq' = {
            res   <-  res[with(res, order(Leq)), ]
        },
        'pc' = {
            res   <-  res[with(res, order(Lpc)), ]
        }
    )
    
    numericB1  <-  is.numeric(refB1)
    if(numericB1) {
        res  <-  res[with(res, order(abs(refB1 - res$b1))), ]
    } 

    #  Plots to accompany best results  #
    if(plots) {        
        dev.new(height=15, width=15)
        outputPlot(res, xall, yall)
    }
    nFits  <-  print(nrow(res))
    list(
        'nFits'  =  nFits,
        'res'    =  res
    )
}

####################
# PLOTTING FUNCTIONS
####################
outputPlot  <-  function(resultsTable, x, y) {
    col1  <-  adjustcolor('#1B6889', alpha=0.5)
    par(mfrow=c(5,5))
    for(i in seq_len(nrow(resultsTable))) {
        # Subset Data #
        ytemp  <-  y[c(resultsTable$Lbound[i]:resultsTable$Rbound[i])]
        xtemp  <-  x[c(resultsTable$Lbound[i]:resultsTable$Rbound[i])]
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
#  Dependency -- bestLocReg():
##############################
# for use with plotBest() -- slight modification of LocReg, that
#  accepts the best window chosen by the user instead of subsetting
#  the dataset. Also uses a list() to output objects of different
#  length (e.g. bHat and stdResid).

bestLocReg  <-  function(bestwin, y, x, weights, verbose=TRUE, ...) {
    # Design Matrix #
    X  <-  matrix(cbind(1, x), ncol=2)
    if(weights) { # Use Weighted Least Squares Regression #
        w           <-  triCube(x=x, m=mean(x), h=(x[length(x)] - x[1]) / 2)
        bHat        <-  (solve(t(X) %*% diag(w) %*% X)) %*% t(X) %*% diag(w) %*% y
        yHat        <-  X %*% bHat
        sigmaHatUb  <-  sum((y - yHat)^2)/(length(y) - 2)
        vcov        <-  sigmaHatUb * (solve(t(X) %*% diag(w) %*% X))
        b1CI        <-  bHat[2, ] + qt(c(0.025, 0.975), df=(length(x) - 2)) * sqrt(diag(vcov))[2]
        stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
        BGtest      <-  breuschGodfrey(y, x, order=FALSE)
        bgN         <-  as.numeric(BGtest$bg)
    } else { # Use Ordinary Least Squares Regression #
        bHat        <-  (solve(t(X) %*% X)) %*% t(X) %*% y
        yHat        <-  X %*% bHat
        sigmaHatUb  <-  sum((y - yHat)^2) / (length(y) - 2)
        vcov        <-  sigmaHatUb * (solve(t(X) %*% X))
        b1CI        <-  bHat[2, ] + qt(c(0.025, 0.975), df=(length(x) - 2)) * sqrt(diag(vcov))[2]
        stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
        BGtest      <-  breuschGodfrey(y, x, order=FALSE)
        bgN         <-  as.numeric(BGtest$bg)
    }
    list(
        'BestWindow' =  bestwin,
        'bHat'       =  bHat,
        'yHat'       =  yHat,
        'sigmaHatUb' =  sigmaHatUb,
        'stdResid'   =  stdResid,
        'skew'       =  skew(x = stdResid),
        'b1CI'       =  b1CI,
        'bg'         =  bgN
    )
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

plotBest <- function(res, yall, xall, best=1) {
    #  Recover data window for chosen local regression model  #
    weight   <-  res$res$weights[1]
    bestwin  <-  c(res$res$Lbound[best], res$res$Rbound[best])
    y        <-  yall[bestwin[1]:bestwin[2]]
    x        <-  xall[bestwin[1]:bestwin[2]]
    #  Fit block  #
    LocFit  <-  bestLocReg(bestwin, y=y, x=x, weight)
    b1      <-  as.numeric(LocFit$bHat[2,1])

    #  Residual Plots  #
    dev.new()
    # Standardized Residuals ~ x
    par(mfrow=c(2,2))
    plot(LocFit$stdResid ~ x, xlab='x', ylab='y', main='Std. Residuals ~ x')
    abline(h=0, col=1, lwd=2)
    abline(h=c(-2, 2), lty=2)
    lf1  <-  loess(LocFit$stdResid ~ x)
    points(x, lf1$fitted, type='l', col=2, lwd=2)
    # Standardized Residuals ~ Fitted Values
    plot(LocFit$stdResid ~ LocFit$yHat, xlab='Fitted Values',ylab='Standardized Residuals',main='Std. Residuals ~ Fitted Values')
    abline(h=0, col=1, lwd=2)
    abline(h=c(-2, 2), lty=2)
    lf2  <-  loess(LocFit$stdResid ~ LocFit$yHat)
    points(LocFit$yHat, lf2$fitted, type='l', col=2, lwd=2)
    # QQNorm Plot of Standardized Residuals #
    qqnorm(LocFit$stdResid, main='QQNorm plot of Std. Residuals')
    qqline(LocFit$stdResid, col=2)
    # Histogram of Standardized Residuals #
    hist(LocFit$stdResid, xlab='Standardized Residuals', ylab='Density', breaks=20, main='Density Plot of Std. Residuals')
    #  Overall Regression Plot  #
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
          x       =  min(res$b1) + (0.8 * (abs(range(results$res$b1)[2] - range(results$res$b1)[1]))),
          y       =  0.95 * max(density(res$b1)$y),
          legend  =  c(expression(paste(italic(L))),
                      expression(paste(italic(L[eq]))),
                      expression(paste(italic(L['%'])))),
          lwd     =  4,
          lty     =  c(1, 2, 3),
          col     =  c(c1, c2, c3),
          cex     =  1
    )
}