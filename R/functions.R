## Main package functions

transparentColor <- function(col, opacity=0.5) {
    if (length(opacity) > 1 && any(is.na(opacity))) {
        n <- max(length(col), length(opacity))
        opacity <- rep(opacity, length.out=n)
        col <- rep(col, length.out=n)
        ok <- !is.na(opacity)
        ret <- rep(NA, length(col))
        ret[ok] <- Recall(col[ok], opacity[ok])
        ret
    } else {
        tmp <- col2rgb(col)/255
        rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
    }
}

proportionalLabel <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
    usr  <-  par('usr')
    x.p  <-  usr[1] + px*(usr[2] - usr[1])
    y.p  <-  usr[3] + py*(usr[4] - usr[3])
    if(log=='x') {
        x.p<-10^(x.p)
    }
    if(log=='y') {
        y.p<-10^(y.p)
    }
    if(log=='xy') {
        x.p<-10^(x.p)
        y.p<-10^(y.p)
    }
    if(text){
        text(x.p, y.p, lab, adj=adj, ...)
    } else {
        points(x.p, y.p, ...)
    }
}

whiteGrid  <-  function(...) {
    proportionalLabel(rep(0.2, 2), c(0,1), text=FALSE, type='l', col='white', lwd=0.5, ...)
    proportionalLabel(rep(0.4, 2), c(0,1), text=FALSE, type='l', col='white', lwd=0.5, ...)
    proportionalLabel(rep(0.6, 2), c(0,1), text=FALSE, type='l', col='white', lwd=0.5, ...)
    proportionalLabel(rep(0.8, 2), c(0,1), text=FALSE, type='l', col='white', lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.2, 2), text=FALSE, type='l', col='white', lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.4, 2), text=FALSE, type='l', col='white', lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.6, 2), text=FALSE, type='l', col='white', lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.8, 2), text=FALSE, type='l', col='white', lwd=0.5, ...)
}

rounded  <-  function(value, precision=1, change=FALSE) {
  if(change) {
    value  <-  value * -1
  }
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}


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
    bgN     <-  BGtest$bgN
    
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
        sigmaHatUb  <-  sum((lmFit$residuals)^2) / (lmFit$df.residual)
        out         <-  list(
                            table      =  out, 
                            residuals  =  lmFit$residuals / (sqrt(sigmaHatUb)),
                            yHat       =  lmFit$fitted
                        )
    }
    out
}

################################################################
#  THE MAIN WRAPPER FUNCTION -- findLocLin():
#############
rankLocReg.default  <-  function(xall, yall, alpha, method=c('ns', 'eq', 'pc'), plots=TRUE, verbose=TRUE) {
    if(is.unsorted(xall))
        warning("Dataset must be ordered by xall")
        
    # make sure that all NAs are dealt with
    dat   <-  stripNAs(xall, yall)
    xall  <-  dat$x
    yall  <-  dat$y

    #  get all possible windows
    wins  <-  getWindows(x=yall, alpha)
    
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
        cat(sprintf('rankLocReg fitted %d local regressions', nFits), '\n')
    
    out  <-  list(
                 'nFits'    =  nFits,
                 'allRegs'  =  allRegs,
                 'xall'     =  xall,
                 'yall'     =  yall
             )

    class(out)  <-  'rankLocReg'
    out
}

rankLocReg <- function(x, ...) {
    UseMethod('rankLocReg')
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
plot.rankLocReg  <-  function(allRegs, rank=1) {
    #  recover data window for chosen local regression model
    bestwin  <-  c(allRegs$allRegs$Lbound[rank], allRegs$allRegs$Rbound[rank])
    y        <-  allRegs$yall[bestwin[1]:bestwin[2]]
    x        <-  allRegs$xall[bestwin[1]:bestwin[2]]
    
    #  fit block
    fit     <-  locReg(bestwin, allRegs$xall, allRegs$yall, resids=TRUE)
    locFit  <-  fit$table
    resids  <-  fit$residuals
    b1      <-  locFit$b1
    yHat    <-  fit$yHat

    #  residual plots
    dev.new(width=9, height=5)

    layout(matrix(c(
                    rep(c(rep(1, 4), rep(2, 2), rep(3, 2)), 2),
                    rep(c(rep(1, 4), rep(4, 2), rep(5, 2)), 2)
                   ), 
           nrow=4, ncol=8, byrow=TRUE)
    )
    
    #  overall regression plot
    outy  <-  allRegs$yall[c(1:(bestwin[1]-1), (bestwin[2]+1):length(allRegs$yall))]
    outx  <-  allRegs$xall[c(1:(bestwin[1]-1), (bestwin[2]+1):length(allRegs$yall))]

    par(mai=c(1.2, 0.8, 0.8, 0.4), cex=1)
    plot(allRegs$yall ~ allRegs$xall, axes=FALSE, type='n', xlab='Predictor', ylab='Response', cex.lab=1.2)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    axis(1, cex.axis=0.9)
    axis(2, las=1, cex.axis=0.9)
    points(outy ~ outx, pch=16, col=transparentColor('black', 0.2), cex=1.2)
    points(y ~ x, col='dodgerblue', cex=1.2)
    lines(x, locFit$b0 + locFit$b1*x, col='black', lwd=2, lty=2)
    proportionalLabel(c(0, 0.14), rep(1.1, 2), text=FALSE, xpd=NA, type='l', lwd=2, lty=2)
    proportionalLabel(0.15, 1.1, substitute('Rank '*pos*': '*italic(y) == a~sy~b%.%italic(x), list(pos=rank, a=rounded(locFit$b0, 2), sy=ifelse(b1 < 0, ' - ', ' + '), b=rounded(abs(b1), 2))), xpd=NA, adj=c(0, 0.5))

    # standardized residuals ~ x
    par(mai=c(0.6732, 0.5412, 0.5412, 0.2772), cex=0.8)
    yRange  <-  max(abs(c(floor(min(resids)), ceiling(max(resids)))))
    yRange  <-  c(-1*yRange, yRange)
    plot(resids ~ x, xlab='Predictor', ylab='Std. residuals', xpd=NA, ylim=yRange, type='n', axes=FALSE)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    axis(1, cex.axis=0.9)
    axis(2, las=1, cex.axis=0.9)
    points(resids ~ x, pch=16, col=transparentColor('dodgerblue', 0.5))
    abline(h=0, col=1, lwd=2)
    abline(h=c(-2, 2), lty=2)
    lf1  <-  loess(resids ~ x)
    lines(x, lf1$fitted, col='tomato', lwd=2)
    
    # standardized residuals ~ fitted values
    plot(resids ~ yHat, xlab='Fitted Values', ylab='Std. residuals', xpd=NA, ylim=yRange, type='n', axes=FALSE)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    axis(1, cex.axis=0.9)
    axis(2, las=1, cex.axis=0.9)
    points(resids ~ yHat, pch=16, col=transparentColor('dodgerblue', 0.5))
    abline(h=0, col=1, lwd=2)
    abline(h=c(-2, 2), lty=2)
    lf2  <-  loess(resids ~ yHat)
    lines(yHat, lf2$fitted, col='tomato', lwd=2)
    
    # qqnorm plot of standardized residuals
    par(mai=c(0.9732, 0.5412, 0.2412, 0.2772), cex=0.8)
    qqPlot  <-  qqnorm(resids, main='QQNorm plot of Std. Residuals', xpd=NA, plot=FALSE)
    plot(y ~ x, data=qqPlot, xlab='Theoretical quantiles', ylab='Sample quantiles', xpd=NA, ylim=yRange, xlim=yRange, type='n', axes=FALSE)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    axis(1, cex.axis=0.9)
    axis(2, las=1, cex.axis=0.9)
    points(qqPlot$y ~ qqPlot$x, pch=16, col=transparentColor('dodgerblue', 0.5))
    qqline(resids, col='tomato')
    
    # histogram of standardized residuals
    histPlot  <-  hist(resids, breaks=20, plot=FALSE)
    plot(NA, xlab='Std. Residuals', ylab='Density', xpd=NA, ylim=c(0, max(histPlot$density)), xlim=yRange, type='n', axes=FALSE)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    axis(1, cex.axis=0.9)
    axis(2, las=1, cex.axis=0.9)
    densities  <-  histPlot$density
    breakPts   <-  histPlot$breaks
    for(j in seq_along(densities)) {
        polygon(c(breakPts[j], breakPts[j+1], breakPts[j+1], breakPts[j], breakPts[j]), c(rep(usr[3], 2), rep(densities[j], 2), usr[3]), border='dodgerblue', col=transparentColor('dodgerblue', 0.5))
    }
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
