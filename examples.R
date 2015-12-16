rm(list=ls())
source('R/functions2.R')

data     <-  read.csv("data/TestO2data.csv", header=TRUE, stringsAsFactors=FALSE)
system.time({
results  <-  FindLocLin(yall=data$D, xall=data$time, alpha=0.2, plots=FALSE, weights=TRUE, verbose=FALSE)
})

PlotBest(res=results, best=1, yall=data$D, xall=data$time)
