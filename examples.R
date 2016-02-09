rm(list=ls())
source('R/functions.R')

# Import test VO2 data #
data     <-  read.csv("data/TestO2data.csv", header=TRUE, stringsAsFactors=FALSE)


##  Test new ref.b1 option  ##
results  <-  FindLocLin(yall=data$D, xall=data$time, alpha=0.2, ref.b1 = FALSE, plots=TRUE, weights=TRUE, verbose=FALSE)
results

PlotBest(res=results, yall=data$D, xall=data$time, best=1)

## Benchmark -- VO2 data ##
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







col1 <- adjustcolor('#1B6889', alpha=0.5)



##  Test CORMORANT FLOW-THROUGH RESPIROMETRY DATA  ##
data     <-  read.csv("data/thinned_cormorant_data.csv", header=TRUE, stringsAsFactors=FALSE)
head(data)
plot(data$Vo2..ml.min. ~ data$Time..h., pch=21, col=1, bg=col1)

results  <-  FindLocLin(yall=data$Vo2..ml.min., xall=data$Time..h., alpha=0.3, ref.b1 = 0.0,
                        plots=TRUE, weights=TRUE, verbose=FALSE)
results

PlotBest(res=results, yall=data$Vo2..ml.min., xall=data$Time..h., best=1)
















