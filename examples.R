rm(list=ls())
source('R/functions.R')
col1 <- adjustcolor('#1B6889', alpha=0.5)
col2 <- adjustcolor('#1B6889', alpha=0.2)
library(extrafont)
#col1 <- adjustcolor('dodgerblue', alpha=0.5)
#col2 <- adjustcolor('dodgerblue', alpha=0.2)



###################################################
##  SEA URCHIN CLOSED-CHAMBER RESPIROMETRY DATA  ##
###################################################

# Import test VO2 data #
data     <-  read.csv("data/TestO2data.csv", header=TRUE, stringsAsFactors=FALSE)

##  Using all=TRUE to examine distributions of L  ##
##   and relation b/w different parts of metric  ##
results  <-  FindLocLin(yall=data$D, xall=data$time, alpha=0.3, method="sr", ref.b1 = FALSE,
                        plots=FALSE, weights=TRUE, all=TRUE, verbose=TRUE)
#toPdf(PlotBest(res=results, yall=data$D, xall=data$time, best=1), filename='BG.Lpc.pdf')
PlotBest(res=results, yall=data$D, xall=data$time, best=1)
head(results$res)

xrange <- (results$res$Rbound - results$res$Lbound)

## Component Metrics  ##
L.skew  <-  ((min(abs(results$res$skew)) + abs(results$res$skew)) / sd(results$res$skew))
L.CI    <- ((results$res$CI.range - min(results$res$CI.range)) / sd(results$res$CI.range))
L.BG    <- ((results$res$bg.df - min(results$res$bg.df)) / sd(results$res$bg.df))

L.skew.sr  <-  (((min(abs(results$res$skew)) + abs(results$res$skew)) / sd(results$res$skew)) /
                (max(((min(abs(results$res$skew)) + abs(results$res$skew)) / sd(results$res$skew)))))
L.CI.sr    <- (((results$res$CI.range - min(results$res$CI.range)) / sd(results$res$CI.range))/
               (max(((results$res$CI.range - min(results$res$CI.range)) / sd(results$res$CI.range)))))
L.BG.sr    <- (((results$res$bg.df - min(results$res$bg.df)) / sd(results$res$bg.df)) /
               (max(((results$res$bg.df - min(results$res$bg.df))/ sd(results$res$bg.df)))))

L.skew.pc <- pc.rank(abs(results$res$skew))
L.CI.pc <- (pc.rank((results$res$CI.range)))
L.BG.pc <- (pc.rank((results$res$bg.df - min(results$res$bg.df))))



##  Histograms of Metric Components
par(mfrow=c(2,3))
hist(L.skew, breaks=50)
hist(L.CI, breaks=50)
hist(L.BG, breaks=50)
hist(L.skew.sr, breaks=50)
hist(L.CI.sr, breaks=50)
hist(L.BG.sr, breaks=50)


##  Effects of different standardizations on each component metric  ##

par(mfrow=c(3,3))
plot(L.skew ~ L.skew.sr, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(Skew))), xlab=expression(paste(italic(Skew[sr]))))
plot(L.skew.pc ~ L.skew, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(Skew["%"]))),xlab=expression(paste(italic(Skew))))
plot(L.skew.pc ~ L.skew.sr, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(BG["%"]))), xlab=expression(paste(italic(BG[sr]))))
plot(L.BG ~ L.BG.sr, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(BG))), xlab=expression(paste(italic(BG[sr]))))
plot(L.BG.pc ~ L.BG, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(BG["%"]))),xlab=expression(paste(italic(BG))))
plot(L.BG.pc ~ L.BG.sr, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(BG["%"]))), xlab=expression(paste(italic(BG[sr]))))
plot(L.CI ~ L.CI.sr, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(CI))), xlab=expression(paste(italic(CI[sr]))))
plot(L.CI.pc ~ L.CI, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(BG["%"]))),xlab=expression(paste(italic(CI))))
plot(L.CI.pc ~ L.CI.sr, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(CI["%"]))), xlab=expression(paste(italic(CI[sr]))))


##  Density plots of Metric components  ##
par(mfrow=c(3,4), cex.lab=1.5)
plot(density(results$res$L), lwd=4, col=col1, xlab=expression(paste(italic(L))),
     main=expression(paste("Distribution of ", italic(L))), cex.main=2)
plot(density(results$res$L.sr), lwd=4, col=col1, xlab=expression(paste(italic(L[sr]))),
     main=expression(paste("Distribution of ", italic(L[sr]))), cex.main=2)
plot(density(results$res$L.pc), lwd=4, col=col1, xlab=expression(paste(italic(L["%"]))),
     main=expression(paste("Distribution of ", italic(L["%"]))), cex.main=2)
plot(density(results$res$b1), lwd=4, col=col1, xlab=expression(paste(beta[1])),
     main=expression(paste("Distribution of ", beta[1])), cex.main=2)
plot(density(L.skew), lwd=4, col=col1, xlab=expression(paste(italic(Skew))),
     main=expression(paste("Distribution of Skew")), cex.main=2)
plot(density(L.CI), lwd=4, col=col1, xlab=expression(paste(italic(C.I.~range))),
     main=expression(paste("Distribution of C.I. range")), cex.main=2)
plot(density(L.BG), lwd=4, col=col1, xlab=expression(paste(BG)),
     main=expression(paste("Distribution of ", BG)), cex.main=2)
plot(NA)
plot(density(L.skew.sr), lwd=4, col=col1, xlab=expression(paste(italic(Skew[sr]))),
     main=expression(paste("Distribution of ", Skew[sr])), cex.main=2)
plot(density(L.CI.sr), lwd=4, col=col1, xlab=expression(paste(italic(C.I.~range[sr]))),
     main=expression(paste("Distribution of C.I. ",range[sr])), cex.main=2)
plot(density(L.BG.sr), lwd=4, col=col1, xlab=expression(paste(BG[sr])),
     main=expression(paste("Distribution of ", BG[sr])), cex.main=2)



##  L ~ component parts  ##
par(mfrow=c(3,3), cex.lab=1.5)
plot(results$res$L ~ L.CI, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab=expression(paste(italic(C.I.~range))))
plot(results$res$L ~  abs(results$res$skew), pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab=expression(paste(italic(Skew))))
plot(results$res$L ~  L.BG, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab=expression(paste(BG)))
plot(results$res$L.sr ~ L.CI.sr, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L[sr]))), xlab=expression(paste(italic(C.I.~range[sr]))))
plot(results$res$L.sr ~  L.skew.sr, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L[sr]))), xlab=expression(paste(italic(C.I.[sr]))))
plot(results$res$L.sr ~  L.BG.sr, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L[sr]))), xlab=expression(paste(BG[sr])))
plot(results$res$L.pc ~ L.CI.pc, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L["%"]))), xlab=expression(paste(italic(C.I.~range["%"]))))
plot(results$res$L.pc ~  L.skew.pc, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L["%"]))), xlab=expression(paste(italic(Skew["%"]))))
plot(results$res$L.pc ~  L.BG.pc, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L["%"]))), xlab=expression(paste(BG["%"])))



##  L,components ~ sample size  ##
par(mfrow=c(2,3), cex.lab=1.5)
plot(results$res$L ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L))), xlab=expression(paste(italic(n))))
plot(results$res$L.sr ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L[sr]))), xlab=expression(paste(italic(n))))
plot(results$res$L.pc ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(L["%"]))), xlab=expression(paste(italic(n))))

par(mfrow=c(3,3), cex.lab=1.5)
plot(L.skew ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(Skew))), xlab=expression(paste(italic(n))))
plot(L.CI ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(C.I.~range))), xlab=expression(paste(italic(n))))
plot(L.BG ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(BG))), xlab=expression(paste(italic(n))))
plot(L.skew.sr ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(Skew[sr]))), xlab=expression(paste(italic(n))))
plot(L.CI.sr ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(C.I.~range[sr]))), xlab=expression(paste(italic(n))))
plot(L.BG.sr ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(BG[sr]))), xlab=expression(paste(italic(n))))
plot(L.skew.pc ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(Skew["%"]))), xlab=expression(paste(italic(n))))
plot(L.CI.pc ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(C.I.~range["%"]))), xlab=expression(paste(italic(n))))
plot(L.BG.pc ~ xrange, pch=21, bg=col2, col=NA,
     ylab=expression(paste(italic(BG["%"]))), xlab=expression(paste(italic(n))))



##  Correlations among L components  ##
par(mfrow=c(3,3), cex.lab=1.5)
plot(L.skew ~ L.BG, pch=21, bg=col2, col=NA,
     xlab=expression(paste(BG)), ylab=expression(paste(Skew)))
plot(L.skew ~ L.CI, pch=21, bg=col2, col=NA,
     xlab=expression(paste(C.I.~range)), ylab=expression(paste(Skew)))
plot(L.CICorrelations among L components  ~ L.BG, pch=21, bg=col2, col=NA,
     xlab=expression(paste(BG)), ylab=expression(paste(C.I.~range)))
plot(L.skew.sr ~ L.BG.sr, pch=21, bg=col2, col=NA,
     xlab=expression(paste(BG[sr])), ylab=expression(paste(Skew[sr])))
plot(L.skew.sr ~ L.CI.sr, pch=21, bg=col2, col=NA,
     xlab=expression(paste(C.I.~range[sr])), ylab=expression(paste(Skew[sr])))
plot(L.CI.sr ~ L.BG.sr, pch=21, bg=col2, col=NA,
     xlab=expression(paste(BG[sr])), ylab=expression(paste(C.I.~range[sr])))
plot(L.skew.pc ~ L.BG.pc, pch=21, bg=col2, col=NA,
     xlab=expression(paste(BG["%"])), ylab=expression(paste(Skew["%"])))
plot(L.skew.pc ~ L.CI.pc, pch=21, bg=col2, col=NA,
     xlab=expression(paste(C.I.~range["%"])), ylab=expression(paste(Skew["%"])))
plot(L.CI.pc ~ L.BG.pc, pch=21, bg=col2, col=NA,
     xlab=expression(paste(BG["%"])), ylab=expression(paste(C.I.~range["%"])))




##  Beta[1] as a function of L, components ##
par(mfrow=c(1,3), cex.lab=1.5)
plot(results$res$b1 ~ results$res$L, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1])), xlab=expression(paste(L)))
plot(results$res$b1 ~ results$res$L.sr, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1])), xlab=expression(paste(L[sr])))
plot(results$res$b1 ~ results$res$L.pc, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1])), xlab=expression(paste(L["%"])))

par(mfrow=c(3,3), cex.lab=1.5)
plot(results$res$b1 ~ L.skew, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1])), xlab=expression(paste(Skew)))
plot(results$res$b1 ~ L.CI, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1])), xlab=expression(paste(C.I.~range)))
plot(results$res$b1 ~ L.BG, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1])), xlab=expression(paste(BG)))
plot(results$res$b1 ~ L.skew.sr, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1])), xlab=expression(paste(Skew[sr])))
plot(results$res$b1 ~ L.CI.sr, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1])), xlab=expression(paste(C.I.~range[sr])))
plot(results$res$b1 ~ L.BG.sr, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1])), xlab=expression(paste(BG[sr])))
plot(results$res$b1 ~ L.skew.pc, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1])), xlab=expression(paste(Skew["%"])))
plot(results$res$b1 ~ L.CI.pc, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1])), xlab=expression(paste(C.I.~range["%"])))
plot(results$res$b1 ~ L.BG.pc, pch=21, bg=col2, col=NA,
     ylab=expression(paste(beta[1])), xlab=expression(paste(BG["%"])))

c1 <- adjustcolor('#A67A01', alpha=0.75)
c2 <- adjustcolor('#6D65FA', alpha=0.75)
c3 <- adjustcolor('#B6084E', alpha=0.75)

par(mfrow=c(1,1))
#pdf(file="beta1.pdf", family="CM Roman", width=6, height=6)
plot(density(results$res$b1), lwd=4, col=col1, xlab=expression(paste(beta[1])),
     main=expression(paste("Distribution of ", beta[1])), cex.main=2)
abline(v=results$res$b1[results$res$L == min(results$res$L)], col=c1, lty=1, lwd=4)
abline(v=results$res$b1[results$res$L.sr == min(results$res$L.sr)], col=c2, lty=2, lwd=4)
abline(v=results$res$b1[results$res$L.pc == min(results$res$L.pc)], col=c3, lty=3, lwd=4)
legend(x = min(results$res$b1) + (0.8 * (abs(range(results$res$b1)[2] - range(results$res$b1)[1]))),
       y = (0.95*max(density(results$res$b1)$y)),
       legend = c(expression(paste(italic(L))),
                  expression(paste(italic(L[sr]))),
                  expression(paste(italic(L["%"])))),
       lwd = 4,
       lty = c(1,2,3),
       col = c(c1, c2, c3),
       cex=1)
#dev.off()
#embed_fonts("beta1.pdf", outfile = "beta1.pdf")


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





