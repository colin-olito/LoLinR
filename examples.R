rm(list=ls())
source('R/functions.R')
col1 <- adjustcolor('#1B6889', alpha=0.5)
col2 <- adjustcolor('#1B6889', alpha=0.2)




###################################################
##  SEA URCHIN CLOSED-CHAMBER RESPIROMETRY DATA  ##
###################################################

# Import test VO2 data #
data     <-  read.csv("data/TestO2data.csv", header=TRUE, stringsAsFactors=FALSE)

##  Test new ref.b1 option  ##
results  <-  FindLocLin(yall=data$D, xall=data$time, alpha=0.3, ref.b1 = FALSE,
                        plots=FALSE, weights=TRUE, all=FALSE, verbose=TRUE)
results



##  Using all=TRUE to examine distributions of L  ##
##   and relations b/w different parts of metric  ##
results  <-  FindLocLin(yall=data$D, xall=data$time, alpha=0.3, ref.b1 = FALSE,
                        plots=FALSE, weights=TRUE, all=TRUE, verbose=TRUE)

xrange <- (results$res$Rbound - results$res$Lbound)
CIrange <- results$res$b1.CI.hi - results$res$b1.CI.lo

##  Density plots of Metric components  ##
par(mfrow=c(2,3), cex.lab=1.5)
plot(density(results$res$L), lwd=4, col=col1, xlab=expression(paste(italic(L))),
     main=expression(paste("Distribution of ", italic(L))), cex.main=2)
plot(density(results$res$b1), lwd=4, col=col1, xlab=expression(paste(beta[1])),
     main=expression(paste("Distribution of ", beta[1])), cex.main=2)
plot(density(results$res$skew), lwd=4, col=col1, xlab=expression(paste(italic(Skew))),
     main=expression(paste("Distribution of Skew")), cex.main=2)
plot(density(results$res$b1.ac), lwd=4, col=col1, xlab=expression(paste(beta[1~acorr])),
     main=expression(paste("Distribution of ", beta[1~acorr])), cex.main=2)
plot(density(CIrange), lwd=4, col=col1, xlab=expression(paste(italic(L))),
     main=expression(paste("Distribution of C.I. range")), cex.main=2)

##  L ~ component parts  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(results$res$L ~ CIrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab="C.I. range")
plot(results$res$L ~  abs(results$res$skew), pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab="Skew")
plot(results$res$L ~  results$res$b1.ac, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab=expression(paste(beta[1~acorr])))

##  L,components ~ sample size  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(results$res$L ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab=expression(paste(italic(n))))
plot(CIrange ~ xrange, pch=21, bg=col2, col=NA,
     ylab="C.I. range", xlab=expression(paste(italic(n))))
plot(abs(results$res$skew) ~ xrange, pch=21, bg=col2, col=NA,
     ylab="Skew", xlab=expression(paste(italic(n))))
plot(results$res$b1.ac ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1~acorr])), xlab=expression(paste(italic(n))))

##  Correlations among L components  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(abs(results$res$skew) ~ results$res$b1.ac, pch=21, bg=col2, col=NA,
     xlab=expression(paste(beta[1~acorr])), ylab="|Skew|")
plot(abs(results$res$skew) ~ CIrange, pch=21, bg=col2, col=NA,
     xlab="C.I. range", ylab="|Skew|")
plot(CIrange ~ results$res$b1.ac, pch=21, bg=col2, col=NA,
     xlab=expression(paste(beta[1~acorr])), ylab="C.I. range")



#########################################################
##  COMPUTATION TIME BENCHMARKS USING SEA URCHIN DATA  ##
#########################################################

system.time({
results  <-  FindLocLin(yall=data$D, xall=data$time, alpha=0.2, plots=FALSE, weights=TRUE, verbose=FALSE)
})[1]

results <- FindLocLin(yall=data$D, xall=data$time, alpha=0.2, plots=FALSE, weights=TRUE, verbose=FALSE)


## Benchmark -- increasing dataset size ##
size <- seq(50,500, by=25)
time <- c()
nRegs <- c()

for (i in 1:length(size)) {
# Generate test data for benchmarks #
xtest <- seq(size[i])
set.seed(123)
e <- rnorm(length(xtest), sd=3)
ytest <- sort(0.1*xtest + 2 + e)
# Benchmark #
time[i] <- system.time({ results <- FindLocLin(yall=ytest, xall=xtest, alpha=0.2, plots=FALSE,
                                       weights=TRUE, verbose=FALSE)})[1]
nRegs[i] <- results$nFits
}
BenchOut <- data.frame(size=size,time=time, nFits=nRegs)
#write.table(BenchOut, file="BenchOut.txt")

bench <- read.table("data/BenchOut.txt", header=TRUE, stringsAsFactors=FALSE)
head(bench)

par(mfrow=c(2,2))
plot(time/60 ~ size,  data=bench, ylab='minutes', xlab='nObs')
plot(time/60 ~ nFits, data=bench, ylab='minutes')
plot(nFits   ~ size,  data=bench, xlab='nObs')








#####################################################
##  Test CORMORANT FLOW-THROUGH RESPIROMETRY DATA  ##
#####################################################

data     <-  read.csv("data/thinned_cormorant_data.csv", header=TRUE, stringsAsFactors=FALSE)
head(data)
plot(data$Vo2..ml.min. ~ data$Time..h., pch=21, col=1, bg=col1)

results  <-  FindLocLin(yall=data$Vo2..ml.min., xall=data$Time..h., alpha=0.3, ref.b1 = FALSE,
                        plots=TRUE, weights=FALSE, verbose=FALSE)
results

PlotBest(res=results, yall=data$Vo2..ml.min., xall=data$Time..h., best=1)



##  Using all=TRUE to examine distributions of L  ##
##   and relations b/w different parts of metric  ##
results  <-  FindLocLin(yall=data$Vo2..ml.min., xall=data$Time..h., alpha=0.3, ref.b1 = FALSE,
                        plots=FALSE, weights=TRUE, all=TRUE, verbose=TRUE)

xrange <- (results$res$Rbound - results$res$Lbound)
CIrange <- results$res$b1.CI.hi - results$res$b1.CI.lo

##  Density plots of Metric components  ##
par(mfrow=c(2,3), cex.lab=1.5)
plot(density(results$res$L), lwd=4, col=col1, xlab=expression(paste(italic(L))),
     main=expression(paste("Distribution of ", italic(L))), cex.main=2)
plot(density(results$res$b1), lwd=4, col=col1, xlab=expression(paste(beta[1])),
     main=expression(paste("Distribution of ", beta[1])), cex.main=2)
plot(density(results$res$skew), lwd=4, col=col1, xlab=expression(paste(italic(Skew))),
     main=expression(paste("Distribution of Skew")), cex.main=2)
plot(density(results$res$b1.ac), lwd=4, col=col1, xlab=expression(paste(beta[1~acorr])),
     main=expression(paste("Distribution of ", beta[1~acorr])), cex.main=2)
plot(density(CIrange), lwd=4, col=col1, xlab=expression(paste(italic(L))),
     main=expression(paste("Distribution of C.I. range")), cex.main=2)

##  L ~ component parts  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(results$res$L ~ CIrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab="C.I. range")
plot(results$res$L ~  abs(results$res$skew), pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab="Skew")
plot(results$res$L ~  results$res$b1.ac, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab=expression(paste(beta[1~acorr])))

##  L,components ~ sample size  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(results$res$L ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab=expression(paste(italic(n))))
plot(CIrange ~ xrange, pch=21, bg=col2, col=NA,
     ylab="C.I. range", xlab=expression(paste(italic(n))))
plot(abs(results$res$skew) ~ xrange, pch=21, bg=col2, col=NA,
     ylab="Skew", xlab=expression(paste(italic(n))))
plot(results$res$b1.ac ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1~acorr])), xlab=expression(paste(italic(n))))

##  Correlations among L components  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(abs(results$res$skew) ~ results$res$b1.ac, pch=21, bg=col2, col=NA,
     xlab=expression(paste(beta[1~acorr])), ylab="|Skew|")
plot(abs(results$res$skew) ~ CIrange, pch=21, bg=col2, col=NA,
     xlab="C.I. range", ylab="|Skew|")
plot(CIrange ~ results$res$b1.ac, pch=21, bg=col2, col=NA,
     xlab=expression(paste(beta[1~acorr])), ylab="C.I. range")




###################################
##  Test TOAD RESPIROMETRY DATA  ##
###################################

data     <-  read.csv("data/thinned_toad_data.csv", header=TRUE, stringsAsFactors=FALSE)
head(data)
plot(data$Fo2 ~ data$Time..s., pch=21, col=1, bg=col1)

res  <-  FindLocLin(yall=data$Fo2, xall=data$Time..s., alpha=0.3, ref.b1 = FALSE,
                        plots=TRUE, weights=FALSE, verbose=FALSE)
res
res2  <-  FindLocLin(yall=data$Fo2, xall=data$Time..s., alpha=0.3, ref.b1 = 0.0,
                        plots=TRUE, weights=FALSE, verbose=FALSE)
res2

PlotBest(res=res, yall=data$Fo2, xall= data$Time..s., best=1)
PlotBest(res=res2, yall=data$Fo2, xall= data$Time..s., best=1)




##  Using all=TRUE to examine distributions of L  ##
##   and relations b/w different parts of metric  ##
results  <-  FindLocLin(yall=data$Fo2, xall= data$Time..s., alpha=0.3, ref.b1 = FALSE,
                        plots=FALSE, weights=TRUE, all=TRUE, verbose=TRUE)

xrange <- (results$res$Rbound - results$res$Lbound)
CIrange <- results$res$b1.CI.hi - results$res$b1.CI.lo

##  Density plots of Metric components  ##
par(mfrow=c(2,3), cex.lab=1.5)
plot(density(results$res$L), lwd=4, col=col1, xlab=expression(paste(italic(L))),
     main=expression(paste("Distribution of ", italic(L))), cex.main=2)
plot(density(results$res$b1), lwd=4, col=col1, xlab=expression(paste(beta[1])),
     main=expression(paste("Distribution of ", beta[1])), cex.main=2)
plot(density(results$res$skew), lwd=4, col=col1, xlab=expression(paste(italic(Skew))),
     main=expression(paste("Distribution of Skew")), cex.main=2)
plot(density(results$res$b1.ac), lwd=4, col=col1, xlab=expression(paste(beta[1~acorr])),
     main=expression(paste("Distribution of ", beta[1~acorr])), cex.main=2)
plot(density(CIrange), lwd=4, col=col1, xlab=expression(paste(italic(L))),
     main=expression(paste("Distribution of C.I. range")), cex.main=2)

##  L ~ component parts  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(results$res$L ~ CIrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab="C.I. range")
plot(results$res$L ~  abs(results$res$skew), pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab="Skew")
plot(results$res$L ~  results$res$b1.ac, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab=expression(paste(beta[1~acorr])))

##  L,components ~ sample size  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(results$res$L ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab=expression(paste(italic(n))))
plot(CIrange ~ xrange, pch=21, bg=col2, col=NA,
     ylab="C.I. range", xlab=expression(paste(italic(n))))
plot(abs(results$res$skew) ~ xrange, pch=21, bg=col2, col=NA,
     ylab="Skew", xlab=expression(paste(italic(n))))
plot(results$res$b1.ac ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1~acorr])), xlab=expression(paste(italic(n))))

##  Correlations among L components  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(abs(results$res$skew) ~ results$res$b1.ac, pch=21, bg=col2, col=NA,
     xlab=expression(paste(beta[1~acorr])), ylab="|Skew|")
plot(abs(results$res$skew) ~ CIrange, pch=21, bg=col2, col=NA,
     xlab="C.I. range", ylab="|Skew|")
plot(CIrange ~ results$res$b1.ac, pch=21, bg=col2, col=NA,
     xlab=expression(paste(beta[1~acorr])), ylab="C.I. range")



########################################
##  Test COCKROACH RESPIROMETRY DATA  ##
########################################

data     <-  read.csv("data/thinned_cockroach_data.csv", header=TRUE, stringsAsFactors=FALSE)
head(data)
plot(data$X.CO2..ppm ~ data$Time..s., pch=21, col=1, bg=col1)

res  <-  FindLocLin(yall=data$X.CO2..ppm, xall=data$Time..s., alpha=0.3, ref.b1 = FALSE,
                        plots=TRUE, weights=FALSE, verbose=FALSE)
res
res2  <-  FindLocLin(yall=data$X.CO2..ppm, xall=data$Time..s., alpha=0.3, ref.b1 = 0.0,
                        plots=TRUE, weights=FALSE, verbose=FALSE)
res2

PlotBest(res=res, yall=data$X.CO2..ppm , xall= data$Time..s., best=1)
PlotBest(res=res2, yall=data$X.CO2..ppm , xall= data$Time..s., best=1)






##  Using all=TRUE to examine distributions of L  ##
##   and relations b/w different parts of metric  ##
results  <-  FindLocLin(yall=data$X.CO2..ppm, xall=data$Time..s., alpha=0.3, ref.b1 = 0.0,
                        plots=FALSE, weights=TRUE, all=TRUE, verbose=TRUE)
PlotBest(res=results, yall=data$X.CO2..ppm , xall= data$Time..s., best=1)

xrange <- (results$res$Rbound - results$res$Lbound)
CIrange <- results$res$b1.CI.hi - results$res$b1.CI.lo

##  Density plots of Metric components  ##
par(mfrow=c(2,3), cex.lab=1.5)
plot(density(results$res$L), lwd=4, col=col1, xlab=expression(paste(italic(L))),
     main=expression(paste("Distribution of ", italic(L))), cex.main=2)
plot(density(results$res$b1), lwd=4, col=col1, xlab=expression(paste(beta[1])),
     main=expression(paste("Distribution of ", beta[1])), cex.main=2)
plot(density(results$res$skew), lwd=4, col=col1, xlab=expression(paste(italic(Skew))),
     main=expression(paste("Distribution of Skew")), cex.main=2)
plot(density(results$res$b1.ac), lwd=4, col=col1, xlab=expression(paste(beta[1~acorr])),
     main=expression(paste("Distribution of ", beta[1~acorr])), cex.main=2)
plot(density(CIrange), lwd=4, col=col1, xlab=expression(paste(italic(L))),
     main=expression(paste("Distribution of C.I. range")), cex.main=2)

##  L ~ component parts  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(results$res$L ~ CIrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab="C.I. range")
plot(results$res$L ~  abs(results$res$skew), pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab="Skew")
plot(results$res$L ~  results$res$b1.ac, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab=expression(paste(beta[1~acorr])))

##  L,components ~ sample size  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(results$res$L ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab=expression(paste(italic(n))))
plot(CIrange ~ xrange, pch=21, bg=col2, col=NA,
     ylab="C.I. range", xlab=expression(paste(italic(n))))
plot(abs(results$res$skew) ~ xrange, pch=21, bg=col2, col=NA,
     ylab="Skew", xlab=expression(paste(italic(n))))
plot(results$res$b1.ac ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1~acorr])), xlab=expression(paste(italic(n))))

##  Correlations among L components  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(abs(results$res$skew) ~ results$res$b1.ac, pch=21, bg=col2, col=NA,
     xlab=expression(paste(beta[1~acorr])), ylab="|Skew|")
plot(abs(results$res$skew) ~ CIrange, pch=21, bg=col2, col=NA,
     xlab="C.I. range", ylab="|Skew|")
plot(CIrange ~ results$res$b1.ac, pch=21, bg=col2, col=NA,
     xlab=expression(paste(beta[1~acorr])), ylab="C.I. range")






########################################################################################
##  Testing logistic growth data, and usefulness of ref.b1 (as suggested by Martino)  ##
########################################################################################

logistic <- function(k, x0, n, xmin, xmax, sd.e) {
    x <- seq(0,10,length.out=200)
    e <- rnorm(n,sd=sd.e)
    l <- 1/(1+ exp(-k*(x - x0) + e))
    list(
        'x' = x,
        'l' = l
    )
}

lg <- logistic(k=2, x0=5, n=100, xmin=2, xmax=10, sd.e=0.7)
plot(lg$l~lg$x)

res1  <-  FindLocLin(yall=lg$l, xall=lg$x, alpha=0.2, ref.b1 = 0.5,
                        plots=FALSE, weights=TRUE, verbose=FALSE)
res1
res2  <-  FindLocLin(yall=lg$l, xall=lg$x, alpha=0.2, ref.b1 = FALSE,
                        plots=FALSE, weights=TRUE, verbose=FALSE)
res2


PlotBest(res=res1, yall=lg$l, xall=lg$x, best=1)
PlotBest(res=res2, yall=lg$l, xall=lg$x, best=1)





##  Using all=TRUE to examine distributions of L  ##
##   and relations b/w different parts of metric  ##
results  <-  FindLocLin(yall=lg$l, xall=lg$x, alpha=0.3, ref.b1 = 0.5,
                        plots=FALSE, weights=TRUE, all=TRUE, verbose=TRUE)

xrange <- (results$res$Rbound - results$res$Lbound)
CIrange <- results$res$b1.CI.hi - results$res$b1.CI.lo

##  Density plots of Metric components  ##
par(mfrow=c(2,3), cex.lab=1.5)
plot(density(results$res$L), lwd=4, col=col1, xlab=expression(paste(italic(L))),
     main=expression(paste("Distribution of ", italic(L))), cex.main=2)
plot(density(results$res$b1), lwd=4, col=col1, xlab=expression(paste(beta[1])),
     main=expression(paste("Distribution of ", beta[1])), cex.main=2)
plot(density(results$res$skew), lwd=4, col=col1, xlab=expression(paste(italic(Skew))),
     main=expression(paste("Distribution of Skew")), cex.main=2)
plot(density(results$res$b1.ac), lwd=4, col=col1, xlab=expression(paste(beta[1~acorr])),
     main=expression(paste("Distribution of ", beta[1~acorr])), cex.main=2)
plot(density(CIrange), lwd=4, col=col1, xlab=expression(paste(italic(L))),
     main=expression(paste("Distribution of C.I. range")), cex.main=2)

##  L ~ component parts  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(results$res$L ~ CIrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab="C.I. range")
plot(results$res$L ~  abs(results$res$skew), pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab="Skew")
plot(results$res$L ~  results$res$b1.ac, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab=expression(paste(beta[1~acorr])))

##  L,components ~ sample size  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(results$res$L ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab=expression(paste(italic(n))))
plot(CIrange ~ xrange, pch=21, bg=col2, col=NA,
     ylab="C.I. range", xlab=expression(paste(italic(n))))
plot(abs(results$res$skew) ~ xrange, pch=21, bg=col2, col=NA,
     ylab="Skew", xlab=expression(paste(italic(n))))
plot(results$res$b1.ac ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1~acorr])), xlab=expression(paste(italic(n))))

##  Correlations among L components  ##
par(mfrow=c(2,2), cex.lab=1.5)
plot(abs(results$res$skew) ~ results$res$b1.ac, pch=21, bg=col2, col=NA,
     xlab=expression(paste(beta[1~acorr])), ylab="|Skew|")
plot(abs(results$res$skew) ~ CIrange, pch=21, bg=col2, col=NA,
     xlab="C.I. range", ylab="|Skew|")
plot(CIrange ~ results$res$b1.ac, pch=21, bg=col2, col=NA,
     xlab=expression(paste(beta[1~acorr])), ylab="C.I. range")





