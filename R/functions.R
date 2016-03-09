######################### MAIN PACKAGE FUNCTIONS
########################

#' Re-rank an object of class \code{rankLocReg}.
#'
#' @title Re-rank \code{rankLocReg}
#' @param x An object of class \code{rankLocReg}.
#' @param newMethod Method with which object of class \code{rankLocReg} should be re-ranked.
#' @details This function also updates \code{x$call} and \code{x$method}.
#' @return An object of class \code{rankLocReg}.
#' @author Diego Barneche.
#' @export
#' @examples
#' # load sea urchin respirometry data
#' data(UrchinData)
#' x  <-  rankLocReg(xall=UrchinData$time, yall=UrchinData$D, alpha=0.8, method="eq", verbose=TRUE)
#' x2  <-  reRank(x, newMethod="pc")
#' # check outputs
#' x$call; x2$call; x$method; x2$method
reRank  <-  function(x, newMethod) {
    isRankLocReg  <-  class(x) == 'rankLocReg'
    if(!isRankLocReg)
        stop('x must be of class rankLocReg')
    isFinalStart  <-  x$method == newMethod
    if(isFinalStart)
        stop(paste0(newMethod, ' is already in place'))
    x$call$method  <-  newMethod
    x$method       <-  newMethod
    x$allRegs  <-  x$allRegs[order(x$allRegs[[paste0('L', newMethod)]]), ]
    x
}

#' Wrapper - strips NA from input x and y.
#'
#' @title Remove NA
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @return A data frame with complete cases.
#' @author Diego Barneche.
#' @export
stripNAs  <-  function(x, y) {
    dat  <-  data.frame(x, y)
    dat[complete.cases(dat), ]
}

#' Wrapper for \code{is.\link[base]{numeric}}.
#'
#' @title Check if vector is numeric
#' @param x A numeric vector.
#' @return Breaks function and returns error message if x
#' is non-numeric.
#' @author Diego Barneche.
#' @export
checkNumeric  <-  function(x) {
    xNumeric  <-  is.numeric(x)
    if(!xNumeric)
        stop('x must be numeric') 
}

#' Wrapper - Checks that lengths of inputs x and y are equal.
#'
#' @title Check for equal lengths
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @return Logical.
#' @author Colin Olito and Diego Barneche.
#' @export
checkEqualLength  <-  function(x, y) {
    equalLength  <-  length(x) == length(y)
    if(!equalLength)
        stop('x and y must be of equal length')
}

#' Thin large data sets to contain fewer observations
#'
#' @title Thin a large data set
#' @param data A data.frame.
#' @param xy A numeric vector designating the columns where x and y are located in data.
#' @param by An integer specifying the thinning rate. Analogous to the \code{by} argument 
#' from the \code{seq} command from the \pkg{base} package.
#' @details This function is used do reduce the size of large data sets to improve
#' overall computation time for \code{\link{rankLocReg}}. \code{\link{locReg}} takes
#' approximately 5 minutes to complete for data sets of ~500 observations, and 
#' computation time increases exponentially with larger data sets. Provided that the
#' data is of suitably high resolution, thinning larger data sets is a reasonable 
#' solution for reducing this computation time.
#' @return A list of data frames with thinned observations.
#' @author Colin Olito.
#' #' @examples
#' # load sea urchin respirometry data
#' data(UrchinData)
#' head(UrchinData)
#' # thin data for individual "D" by 1/2
#' newData  <-  thinData(UrchinData, xy=c(1,5), by=2)
#' head(newData[[1]])
#' @export
thinData  <-  function(data, xy=c(1,2), by=2) {
    # make sure input does not contain factors
    isfactor    <-  any(is.factor(data[,xy[1]]), is.factor(data[,xy[2]]))
    if(isfactor)
        stop('data must contain only numeric variables')

    # make sure that all NAs are dealt with
    dat         <-  stripNAs(data[,xy[1]], data[,xy[2]])
    names(dat)  <-  names(data)[xy]

    # thin data set
    keep  <-  seq(from=1, to=nrow(dat), by=by)
    
    # store thinned data set(s) as list
    newData    <-  list()
    listNames  <-  c()
    for(i in 1:by) {
        thin_dat             <-  dat[(keep + (i - 1)),]
        newData[[i]]         <-  thin_dat[complete.cases(thin_dat), ]
        listNames[i]         <-  paste("newData", i, sep="")
    }
    names(newData)  <-  listNames
    return(newData)
}

#' Calculate the percentile values of a vector x.
#'
#' @title Calculate the percentile values of x
#' @param x A numeric vector.
#' @return A numeric vector of percentiles truncated between 0 and 1.
#' @export
pcRank  <-  function(x) {
    checkNumeric(x)
    percentiles  <-  trunc(rank(x, na.last=NA)) / sum(!is.na(x))
    allUnique    <-  length(percentiles) == length(unique(percentiles))
    if(!allUnique)
        warning('input/output have ties')
    percentiles
}

#' Sample skewness (Fisher-Pearson Standardized Third Moment Coefficient).
#'
#' @title Sample skewness
#' @param x A numeric vector.
#' @param na.rm Logical. Should NAs be removed?
#' @details This function is a dependency for \code{\link{rankLocReg}}
#' where it is used to calculate the (sample) skewness of standardized residuals.
#' @return A numeric vector of length 1.
#' @author Colin Olito.
#' @export
skew  <-  function(x, na.rm=TRUE) {
    checkNumeric(x)
    if(na.rm)
        x  <-  x[!is.na(x)]
    n  <-  length(x)
    (n/((n - 1) * (n - 2))) * sum(((x - mean(x)) / sd(x))^3)
}

#' Get all possible windows between specified alpha.
#' and 1
#' 
#' @title Get all possible windows
#' 
#' @param x A numeric vector.
#' @param alpha Window size. Needs to be higher than 0 and lower or equal to 1.
#' @details This function is a dependency for \code{\link{rankLocReg}}
#' where it is used to extract all local windows
#' for local regressions. \code{alpha} must be higher than 0 and lower or equal to 1. 
#' @return A matrix of vector positions, with starting value on first column and ending value on second column.
#' @author Colin Olito and Diego Barneche.
#' @export
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

#' Modified Breusch-Godfrey Statistic \eqn{\frac{n * R^2}{n}}.
#'
#' @title Breusch-Godfrey Statistic
#' @param x A numeric vector.
#' @param y A numeric vector.
#' @param order Order to which residuals are lagged. Defaults to \code{order <- (n - k - 1)},
#' where \code{n} is the number of observations, and \code{k} is the number of parameters
#' in the regression (2 by default). This represents the highest possible order given \code{n}.
#' @param fill Defaults to \code{fill = 0}, used to fill model matrix for lagged residuals in 
#' the auxiliary regression.
#' @details NOTE: This function is a (very slightly) modified version 
#' of \code{\link[lmtest]{bgtest}} from the \pkg{lmtest} package (available at
#' \url{https://github.com/cran/lmtest/blob/master/R/bgtest.R}).
#' We have stripped it down to minimal functionality for our purposes.
#' All development credit goes to the authors of \pkg{lmtest}.
#'
#' This function is a dependency for \code{\link{rankLocReg}} where it is
#' used to calculate \code{bgN} (Breusch-Godfrey Statistic) divided by the number 
#' of observations \eqn{\frac{n * R^2}{n}}. For the purposes of \code{\link{rankLocReg}}, 
#' only the relative variance explained by the fitted values from the
#' auxiliary regression and the residuals of the original regression is of
#' interest. Calculating \eqn{\frac{n * R^2}{n}} preserves this information, while 
#' avoiding the introduction of strong covariance between bgN and n; an 
#' undesirable behaviour for the linearity metric L. If desired, users can 
#' reproduce a standard Breusch-Godfrey Chi-squared test of significance 
#' by running \code{\link[stats]{qchisq}} with the output \code{bgN} and \code{df}.
#' However, we would recommend using the function \code{\link[lmtest]{bgtest}} from the
#' package \pkg{lmtest}, as it is specifically designed for this purpose.
#' @author Colin Olito.
#' @import lmtest
#' @return A list with: the standard BG statistic (\code{bg}), BG/n (\code{bgN}), and d.f. (\code{df}).
#' @export
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

#' Local Regression Function.
#'
#' @title Perform local regression for specified window
#' 
#' @param wins A matrix of vector positions, where wins[,1] contains left boundary positions,
#' and wins[,2] contains right boundary positions for the desired local regressions. Default
#' behaviour is optimized for use with \code{\link{getWindows}} in \code{\link{rankLocReg}},
#' but can accept user-defined windows. 
#' @param xall A numeric vector.
#' @param yall A numeric vector.
#' @param resids Logical. Are residuals to be returned?
#' @return A data frame if \code{resids} is \code{FALSE}. A list if \code{resids} is \code{TRUE}.
#' @seealso \code{\link{getWindows}}, \code{\link{rankLocReg}}.
#' @author Colin Olito and Diego Barneche.
#' @export
#' @examples
#' # load sea urchin respirometry data
#' data(UrchinData)
#' # rank L metric by method 'eq'
#' wins  <-  getWindows(x=UrchinData$D, alpha=0.3)
#' reg   <-  locReg(wins[1, ], xall=UrchinData$time, yall=UrchinData$D, resids=TRUE)
locReg  <-  function(wins, xall, yall, resids=FALSE) {
    # check that window is ok
    badwin  <-  wins[1] >= wins[2]
    if(badwin)
        stop('Left window boundary greater than or equal to Right boundary')

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

#' Wrapper - calls function that creates class \code{rankLocReg}.
#'
#' @title Ranking local linear regressions according to linearity
#' @param xall A numeric vector.
#' @param yall A numeric vector.
#' @param alpha Minimum window size, expressed as a proportion of total data. Must be higher than 0 and less-than or equal to 1.
#' @param method Ranking method. See details.
#' @param verbose Logical. Should progress be printed?
#' @details \code{rankLocReg} is the main function around which \pkg{LoLinR} is built. Given independent and dependent variables
#' for a time series or trace data set, \code{rankLocReg} tries to find the 'most linear' ordered subset of the full data set. 
#' This is accomplished by fitting all possible local linear regressions with minimum window size \code{alpha}, and ranking them
#' according to the combined linearity metric $L$. $L$ quantifies linearity from 1) the skewness of the standardized residuals, 
#' 2) the range of the 95% confidence interval around the regression slope $\beta_1$, and 3) auto-correlation among the standardized
#' residuals (a modified Breusch-Godfrey $R^2$). These three components of $L$ can be weighted in 3 different ways: unweighted (\code{method="z"}), 
#' equal weights (\code{method="eq"}), and percentile ranks (\code{method="pc"}). If method is unspecified, default to \code{z}. 
#' For highly skewed, or otherwise ill-behaved data we strongly advise examining the relative behaviour of the different weighting
#' methods using \code{\link{plotBeta1}}.
#' 
#' For data sets with greater tha ~500 observations, we suggest thinning the number of observations using \code{\link{thinData}},
#' given that this does not compromise the resolution of the data for the question of interest.
#' @return A data frame with important output from all local regressions, ranked by metric \code{L} following raking method chosen by argument \code{method}.
#' @author Colin Olito and Diego Barneche.
#' @seealso \code{\link{locReg}}, \code{\link{thinData}}, \code{\link{plotBeta1}}.
#' @export
#' @examples
#' # load sea urchin respirometry data
#' data(UrchinData)
#' # rank L metric by method 'eq'
#' allRegs  <-  rankLocReg(xall=UrchinData$time, yall=UrchinData$D, alpha=0.3, method="eq", verbose=TRUE)
rankLocReg <- function(xall, yall, alpha, method=c('z', 'eq', 'pc'), verbose=TRUE) {
    UseMethod('rankLocReg')
}

#' Ranking local linear regressions.
#'
#' @title Ranking local linear regressions according to linearity
#' @param xall A numeric vector.
#' @param yall A numeric vector.
#' @param alpha Minimum window size, expressed as a proportion of total data. Must be higher than 0 and less-than or equal to 1.
#' @param method Ranking method. See details.
#' @param verbose Logical. Should progress be printed?
#' @details \code{rankLocReg} is the main function around which \pkg{LoLinR} is built. Given independent and dependent variables
#' for a time series or trace data set, \code{rankLocReg} tries to find the 'most linear' ordered subset of the full data set. 
#' This is accomplished by fitting all possible local linear regressions with minimum window size \code{alpha}, and ranking them
#' according to the combined linearity metric $L$. $L$ quantifies linearity from 1) the skewness of the standardized residuals, 
#' 2) the range of the 95% confidence interval around the regression slope $\beta_1$, and 3) auto-correlation among the standardized
#' residuals (a modified Breusch-Godfrey $R^2$). These three components of $L$ can be weighted in 3 different ways: unweighted (\code{method="z"}), 
#' equal weights (\code{method="eq"}), and percentile ranks (\code{method="pc"}). If method is unspecified, default to \code{z}. 
#' For highly skewed, or otherwise ill-behaved data we strongly advise examining the relative behaviour of the different weighting
#' methods using \code{\link{plotBeta1}}.
#' 
#' For data sets with greater tha ~500 observations, we suggest thinning the number of observations using \code{\link{thinData}},
#' given that this does not compromise the resolution of the data for the question of interest.
#' @return A data frame with important output from all local regressions, ranked by metric \code{L} following raking method chosen by argument \code{method}.
#' @author Colin Olito and Diego Barneche.
#' @seealso \code{\link{locReg}}, \code{\link{thinData}}, \code{\link{plotBeta1}}.
#' @export
#' @examples
#' # load sea urchin respirometry data
#' data(UrchinData)
#' # rank L metric by method 'eq'
#' allRegs  <-  rankLocReg(xall=UrchinData$time, yall=UrchinData$D, alpha=0.3, method="eq", verbose=TRUE)
rankLocReg.default  <-  function(xall, yall, alpha, method=c('z', 'eq', 'pc'), verbose=TRUE) {
    
    if(missing(method))
        method  <-  'z'
    
    # make sure input does not contain characters
    checkNumeric(c(xall, yall))

    # make sure that all NAs are dealt with
    dat   <-  stripNAs(xall, yall)
    xall  <-  dat$x
    yall  <-  dat$y

    # x needs to be sorted
    if(is.unsorted(xall))
        warning("Dataset must be ordered by xall")

    #  get all possible windows
    wins  <-  getWindows(x=yall, alpha)
    
    #  fit local regressions
    allRegs   <-  apply(wins, 1, locReg, xall=xall, yall=yall)
    allRegs   <-  do.call(rbind.data.frame, allRegs)
    
    #  calculate combined metric (L) for linearity & fit
    allRegs$ciRange  <-  allRegs$b1UpCI - allRegs$b1LoCI
    allRegs          <-  allRegs[, c('Lbound', 'Rbound', 'alph', 'b0', 'b1', 'b1LoCI', 'b1UpCI', 'ciRange', 'skew', 'bgN')]
    allRegs$Lz       <-  ((min(abs(allRegs$skew)) + abs(allRegs$skew)) / sd(allRegs$skew)) + ((allRegs$bgN - min(allRegs$bgN)) / sd(allRegs$bgN)) + ((allRegs$ciRange - min(allRegs$ciRange)) / sd(allRegs$ciRange))
    allRegs$Leq      <-  (((min(abs(allRegs$skew)) + abs(allRegs$skew)) / sd(allRegs$skew)) / (max(((min(abs(allRegs$skew)) + abs(allRegs$skew)) / sd(allRegs$skew))))) + (((allRegs$bgN - min(allRegs$bgN)) / sd(allRegs$bgN)) / (max(((allRegs$bgN - min(allRegs$bgN)) / sd(allRegs$bgN))))) + (((allRegs$ciRange - min(allRegs$ciRange)) / sd(allRegs$ciRange)) / (max(((allRegs$ciRange - min(allRegs$ciRange)) / sd(allRegs$ciRange)))))
    allRegs$Lpc      <-  ((pcRank(abs(allRegs$skew))) + (pcRank((allRegs$bgN - min(allRegs$bgN)))) + (pcRank((allRegs$ciRange)))) / 3
    
    # choose weighting scheme for linearity metric L
    switch(match.arg(method),
        'z' = {
            allRegs   <-  allRegs[with(allRegs, order(Lz)), ]
        },
        'eq' = {
            allRegs   <-  allRegs[with(allRegs, order(Leq)), ]
        },
        'pc' = {
            allRegs   <-  allRegs[with(allRegs, order(Lpc)), ]
        }
    )
    
    nFits  <-  nrow(allRegs)
    if(verbose)
        cat(sprintf('rankLocReg fitted %d local regressions', nFits), '\n')
    
    out  <-  list(
                 'nFits'    =  nFits,
                 'allRegs'  =  allRegs,
                 'xall'     =  xall,
                 'yall'     =  yall,
                 'call'     =  match.call(),
                 'method'   =  method
             )

    class(out)  <-  'rankLocReg'
    out
}

####################
# PLOTTING FUNCTIONS
####################

#' Plotting chosen local linear regression.
#'
#' @title Plotting chosen local linear regression
#' @param x An object of class \code{rankLocReg}.
#' @param ... Other parameters to be passed through to plotting functions.
#' @param rank Position, as in row number from input \code{x}, of local regression to be plotted.
#' @return Generates a scatterplot + residual-plot diagnostics for chosen local regression.
#' @seealso \code{\link{rankLocReg.default}}.
#' @export
#' @author Colin Olito and Diego Barneche.
#' @examples
#' # load sea urchin respirometry data
#' data(UrchinData)
#' # rank L metric by method 'eq'
#' allRegs  <-  rankLocReg(xall=UrchinData$time, yall=UrchinData$D, alpha=0.3, method="eq", verbose=TRUE)
#' # plot best local regression
#' plot(allRegs, rank=1)
plot.rankLocReg  <-  function(x, ..., rank=1) {
    #  recover data window for chosen local regression model
    bestwin  <-  c(x$allRegs$Lbound[rank], x$allRegs$Rbound[rank])
    y1       <-  x$yall[bestwin[1]:bestwin[2]]
    x1       <-  x$xall[bestwin[1]:bestwin[2]]
    
    #  fit block
    fit     <-  locReg(bestwin, x$xall, x$yall, resids=TRUE)
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
    outy  <-  x$yall[c(1:(bestwin[1]-1), (bestwin[2]+1):length(x$yall))]
    outx  <-  x$xall[c(1:(bestwin[1]-1), (bestwin[2]+1):length(x$yall))]

    par(mai=c(1.2, 0.8, 0.8, 0.4), cex=1)
    plot(x$yall ~ x$xall, axes=FALSE, type='n', xlab='Predictor', ylab='Response', cex.lab=1.2, ylim=c(min(x$yall), (max(x$yall) + 0.1*(max(x$yall) - min(x$yall)))))
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    axis(1, cex.axis=0.9)
    axis(2, las=1, cex.axis=0.9)
    points(outy ~ outx, pch=16, col=transparentColor('black', 0.2), cex=1.2)
    points(y1 ~ x1, col='dodgerblue', cex=1.2)
    lines(x1, locFit$b0 + locFit$b1*x1, col='black', lwd=2, lty=2)
    proportionalLabel(c(0, 0.14), rep(1.1, 2), text=FALSE, xpd=NA, type='l', lwd=2, lty=2)
    proportionalLabel(0.15, 1.1, substitute('Rank '*pos*': '*italic(y) == a~sy~b%.%italic(x), list(pos=rank, a=rounded(locFit$b0, 2), sy=ifelse(b1 < 0, ' - ', ' + '), b=rounded(abs(b1), 2))), xpd=NA, adj=c(0, 0.5))
    proportionalLabel(c(0, 0.14), rep(1.1, 2), text=FALSE, xpd=NA, type='l', lwd=2, lty=2)
    proportionalLabel(0.95, 0.95, paste0('n = ', length(y1)), xpd=NA, adj=c(1, 0.5), font=3, col='dodgerblue')

    # standardized residuals ~ x
    par(mai=c(0.6732, 0.5412, 0.5412, 0.2772), cex=0.8)
    yRange  <-  max(abs(c(floor(min(resids)), ceiling(max(resids)))))
    yRange  <-  c(-1*yRange, yRange)
    plot(resids ~ x1, xlab='Predictor', ylab='Std. residuals', xpd=NA, ylim=yRange, type='n', axes=FALSE)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    axis(1, cex.axis=0.9)
    axis(2, las=1, cex.axis=0.9)
    points(resids ~ x1, pch=16, col=transparentColor('dodgerblue', 0.5))
    abline(h=0, col=1, lwd=2)
    abline(h=c(-2, 2), lty=2)
    lf1  <-  loess(resids ~ x1)
    lines(x1, lf1$fitted, col='tomato', lwd=2)
    
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
    plot(y1 ~ x1, data=qqPlot, xlab='Theoretical quantiles', ylab='Sample quantiles', xpd=NA, ylim=yRange, xlim=yRange, type='n', axes=FALSE)
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

#' Plotting 25 best local linear regressions.
#'
#' @title Plotting 25 best local linear regressions
#' @param allRegs An object of class \code{rankLocReg}.
#' @return Generates scatterplots for 25 best local regressions.
#' @author Colin Olito and Diego Barneche.
#' @seealso \code{rankLocReg.default}
#' @export
#' @examples
#' # load sea urchin respirometry data
#' data(UrchinData)
#' # rank L metric by method 'eq'
#' allRegs  <-  rankLocReg(xall=UrchinData$time, yall=UrchinData$D, alpha=0.3, method="eq", verbose=TRUE)
#' outputRankLocRegPlot(allRegs)
outputRankLocRegPlot  <-  function(allRegs) {
    dev.new(width=7, height=7)
    par(mfrow=c(5,5), omi=rep(1, 4), mai=rep(0,4), cex=1)
    locFit  <-  allRegs$allRegs

    for(i in 1:25) {
        # subset data
        outy  <-  allRegs$yall[c(1:(locFit$Lbound[i]-1), (locFit$Rbound[i]+1):length(allRegs$yall))]
        outx  <-  allRegs$xall[c(1:(locFit$Lbound[i]-1), (locFit$Rbound[i]+1):length(allRegs$yall))]
        y     <-  allRegs$yall[locFit$Lbound[i]:locFit$Rbound[i]]
        x     <-  allRegs$xall[locFit$Lbound[i]:locFit$Rbound[i]]

        # plot
        plot(allRegs$yall ~ allRegs$xall, axes=FALSE, type='n', xlab='Predictor', ylab='Response', cex.lab=1.2, ylim=c(min(allRegs$yall), (max(allRegs$yall) + 0.25*(max(allRegs$yall) - min(allRegs$yall)))))
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
        whiteGrid()
        box()

        # check whether axes and labels are to be plotted
        if(i %in% seq(1, 21, 5))
            axis(2, las=1, cex.axis=0.5)
        if(i %in% 21:25)
            axis(1, cex.axis=0.5, mgp=c(3, 0.5, 0))
        
        points(outy ~ outx, pch=16, col=transparentColor('black', 0.2), cex=1.2)
        points(y ~ x, col='dodgerblue', cex=0.8)
        lines(x, locFit$b0[i] + locFit$b1[i]*x, col='black', lwd=2, lty=2)
        proportionalLabel(0.03, 0.9, substitute(italic(z)*italic(y) == a~sy~b%.%italic(x), list(z=paste0(i, ';   '), a=rounded(locFit$b0[i], 2), sy=ifelse(locFit$b1[i] < 0, ' - ', ' + '), b=rounded(abs(locFit$b1[i]), 2))), adj=c(0, 0.5), cex=0.5)
    }
    mtext('Response', side=2, line=2.5, outer=TRUE)
    mtext('Predictor', side=1, line=2.5, outer=TRUE)
}

#' Distribution of all local slopes.
#'
#' @title Distribution of all local slopes
#' @param allRegs An object of class \code{rankLocReg}.
#' @details Generates a distribution of all local regression slopes.
#' @seealso \code{\link{rankLocReg.default}}
#' @author Colin Olito and Diego Barneche.
#' @export
#' @examples
#' # load sea urchin respirometry data
#' data(UrchinData)
#' # rank L metric by method 'eq'
#' allRegs  <-  rankLocReg(xall=UrchinData$time, yall=UrchinData$D, alpha=0.3, method="eq", verbose=TRUE)
#' plotBeta1(allRegs)
plotBeta1 <- function(allRegs) {

    c1  <-  'tomato'
    c2  <-  'darkolivegreen'
    c3  <-  'dodgerblue4'
    
    dev.new(width=7, height=7)
    par(omi=rep(0.5, 4), cex=1)
    locFit     <-  allRegs$allRegs
    b1Density  <-  density(locFit$b1)

    plot(NA, xlab=expression(paste(beta[1])), type='n', axes=FALSE, ylab='Density', cex.lab=1.2, xlim=range(b1Density$x), ylim=c(0, (max(b1Density$y)+0.05*max(b1Density$y))), yaxs='i')
    proportionalLabel(0.5, 1.1, expression(paste('Distribution of ', beta[1])), xpd=NA, adj=c(0.5, 0.5), font=3, cex=2)
    usr  <-  par('usr')
    rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
    whiteGrid()
    box()
    polygon(c(b1Density$x), c(b1Density$y), col=transparentColor('dodgerblue2', 0.5), border='dodgerblue2')
    axis(1)
    axis(2, las=1)

    abline(v=locFit$b1[locFit$Lz  == min(locFit$Lz)], col=c1, lty=1, lwd=3)
    abline(v=locFit$b1[locFit$Leq == min(locFit$Leq)], col=c2, lty=2, lwd=3)
    abline(v=locFit$b1[locFit$Lpc == min(locFit$Lpc)], col=c3, lty=3, lwd=3)
    legend(
          x       =  usr[2],
          y       =  usr[4],
          legend  =  c(expression(paste(italic(L))),
                      expression(paste(italic(L[eq]))),
                      expression(paste(italic(L['%'])))),
          lwd     =  4,
          lty     =  c(1, 2, 3),
          col     =  c(c1, c2, c3),
          cex     =  1,
          xjust   =  1,
          yjust   =  1,
          bty     =  'n',
          border  =  NA
    )
}

#' Creates transparent colours
#'
#' @title Creates transparent colours
#' @param col Colour.
#' @param opacity Relative y-axis position (in proportion) where character is to be plotted.
#' @export
transparentColor <- function(col, opacity=0.5) {
    if (length(opacity) > 1 && any(is.na(opacity))) {
        n        <-  max(length(col), length(opacity))
        opacity  <-  rep(opacity, length.out=n)
        col      <-  rep(col, length.out=n)
        ok       <-  !is.na(opacity)
        ret      <-  rep(NA, length(col))
        ret[ok]  <-  Recall(col[ok], opacity[ok])
        ret
    } else {
        tmp  <-  col2rgb(col)/255
        rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
    }
}

#' Plot text or points according to relative axis position.
#'
#' @title Plot text or points according to relative axis position
#' @param px Relative x-axis position (in proportion) where character is to be plotted.
#' @param py Relative y-axis position (in proportion) where character is to be plotted.
#' @param lab Plotted text. Works if argument \code{\link[graphics]{text}} is \code{TRUE}.
#' @param adj See argument of same name in R base function \code{\link[graphics]{par}}.
#' @param text Logical. Should text or points be plotted?
#' @param log Used if the original plot uses the argument log, e.g. \code{log='x'}, \code{log='y'} or \code{log='xy'}.
#' @param ... Additional arguments to R base function \code{\link[graphics]{text}}.
#' @export
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

#' Draw equally-spaced white lines on plot window.
#'
#' @title Equally-spaced white lines on plot window
#' @param ... Additional arguments to internal function \code{\link{proportionalLabel}}.
#' @author Diego Barneche
#' @export
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

#' Internal. Create nice rounded numbers for plotting.
#'
#' @title Rounded numbers for plotting
#' @param value A numeric vector.
#' @param precision Number of rounding digits.
#' @return A character vector.
#' @author Diego Barneche.
#' @export
rounded  <-  function(value, precision=1) {
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}

######################
# AUXILLIARY FUNCTIONS
######################

#' Summary for object of class \code{rankLocReg}.
#'
#' @title Summary for object of class \code{rankLocReg}
#' @param object An object of class \code{rankLocReg}.
#' @param ... Additional arguments affecting the summary produced.
#' @return A summary list with main features calculated by function \code{\link{rankLocReg}}.
#' @seealso \code{\link{rankLocReg}}.
#' @author Diego Barneche.
#' @export
#' @examples
#' # load sea urchin respirometry data
#' data(UrchinData)
#' # rank L metric by method 'eq'
#' allRegs  <-  rankLocReg(xall=UrchinData$time, yall=UrchinData$D, alpha=0.3, method="eq", verbose=TRUE)
#' summary(allRegs)
summary.rankLocReg <- function(object, ...) {
    out <- list(
                call          =  object$call,
                data          =  summary(data.frame(xall=object$xall, yall=object$yall)),
                summaryTable  =  head(object$allRegs),
                nFits         =  object$nFits,
                method        =  object$method
           )

    class(out) <- 'summary.rankLocReg'
    out
}

#' Wrapper summary for object of class \code{rankLocReg}.
#'
#' @title Wrapper summary for object of class \code{rankLocReg}
#' @param x An object of class \code{rankLocReg}.
#' @param ... Further arguments passed to or from other methods.
#' @return A summary list with main features calculated by function \code{\link{rankLocReg}}.
#' @seealso \code{\link{rankLocReg}}, \code{\link{summary.rankLocReg}}
#' @author Diego Barneche.
#' @export
#' @examples
#' # load sea urchin respirometry data
#' data(UrchinData)
#' # rank L metric by method 'eq'
#' allRegs  <-  rankLocReg(xall=UrchinData$time, yall=UrchinData$D, alpha=0.3, method="eq", verbose=TRUE)
#' summary(allRegs)
print.summary.rankLocReg <- function(x, ...) {
    cat('Call:\n')
    print(x$call)
    cat('\n')

    cat('Number of fitted local regressions:\n')
    print(x$nFits)
    cat('\n')

    cat('Used dataset:\n')
    print(x$data)
    cat('\n')

    cat('Best ', nrow(x$summaryTable), ' local regressions (L-metric ranking) fitted with method ', x$method, '\n')
    print(x$summaryTable)
    cat('\n')
}
