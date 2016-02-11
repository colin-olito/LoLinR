rm(list=ls())
source('R/functions.R')
col1 <- adjustcolor('#1B6889', alpha=0.5)


# Import test VO2 data #
data     <-  read.csv("data/TestO2data.csv", header=TRUE, stringsAsFactors=FALSE)


##  Test new ref.b1 option  ##
results  <-  FindLocLin(yall=data$D, xall=data$time, alpha=0.3, ref.b1 = FALSE,
                        plots=TRUE, weights=TRUE, verbose=FALSE)
results

PlotBest(res=results, yall=data$D, xall=data$time, best=1)





## Benchmarks -- VO2 data ##
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









##  Test CORMORANT FLOW-THROUGH RESPIROMETRY DATA  ##
data     <-  read.csv("data/thinned_cormorant_data.csv", header=TRUE, stringsAsFactors=FALSE)
head(data)
plot(data$Vo2..ml.min. ~ data$Time..h., pch=21, col=1, bg=col1)

results  <-  FindLocLin(yall=data$Vo2..ml.min., xall=data$Time..h., alpha=0.3, ref.b1 = FALSE,
                        plots=TRUE, weights=FALSE, verbose=FALSE)
results

PlotBest(res=results, yall=data$Vo2..ml.min., xall=data$Time..h., best=1)



##  Not really sure what to do with this... I am pretty sure that our function
##  is doing what it is supposed to do... but I'm having a hard time understanding
##  why it is picking the local regressions that it is.  Basically, it won't give a
##  reasonable answer unless alpha = 0.4... but I'm confused as to why it's returning
##  bad answers with lower alpha.  Also, working with this data set makes me doubt
##  the utility of having ref.b1... it just seems to give strange results.

#  alpha = 0.3, ref.b1 = FALSE:  292    426
# r2        skew       L
0.3519084  0.23391666   -0.5772622

#  alpha = 0.4, ref.b1 = FALSE:  199    390
# r2        skew       L
0.3088907  0.2051832    -0.4559385

summary.lm
summary(lm(data$Vo2..ml.min.[292:426] ~ data$Time..h.[292:426]))
summary(lm(data$Vo2..ml.min.[199:390] ~ data$Time..h.[199:390]))
summary(lm(data$Vo2..ml.min.[292:426] ~ data$Time..h.[292:426], weights=wts))
summary(lm(data$Vo2..ml.min.[199:390] ~ data$Time..h.[199:390], weights=wts2))

wts <- TriCube(x=data$Time..h.[292:426],
        m=mean(data$Time..h.[292:426]),
        h=(data$Time..h.[292:426][length(data$Time..h.[292:426])] - data$Time..h.[292:426][1])/2)
plot(wts ~ data$Time..h.[292:426])



wts2 <- TriCube(x=data$Time..h.[199:390],
        m=mean(data$Time..h.[199:390]),
        h=(data$Time..h.[199:390][length(data$Time..h.[199:390])] - data$Time..h.[199:390][1])/2)









##  Test TOAD RESPIROMETRY DATA  ##
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







##  Test COCKROACH RESPIROMETRY DATA  ##
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










##  Here I am testing with made up data to see if our ref.b1 is useful in
##  the context of finding maximum slopes. (as suggested by Martino)

logistic <- function(k, x0, n, xmin, xmax, sd.e) {
    x <- seq(0,10,length.out=200)
    e <- rnorm(n,sd=sd.e)
    l <- 1/(1+ exp(-k*(x - x0) + e))
    list(
        'x' = x,
        'l' = l
    )
}

lg <- logistic(k=2, x0=5, n=100, xmin=2, xmax=10, sd.e=0.5)
plot(lg$l~lg$x)

res1  <-  FindLocLin(yall=lg$l, xall=lg$x, alpha=0.2, ref.b1 = 0.5,
                        plots=FALSE, weights=TRUE, verbose=FALSE)
res1
res2  <-  FindLocLin(yall=lg$l, xall=lg$x, alpha=0.2, ref.b1 = FALSE,
                        plots=FALSE, weights=TRUE, verbose=FALSE)
res2


PlotBest(res=res, yall=lg$l, xall=lg$x, best=1)
PlotBest(res=res2, yall=lg$l, xall=lg$x, best=1)
