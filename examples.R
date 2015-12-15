source('R/functions.R')

###################################################
#  Examples:
###################################################

data     <-  read.csv("data/TestO2data.csv", header=TRUE, stringsAsFactors=FALSE)
results  <-  FindLocLin(yall=data$D, xall=data$time, alpha=0.2, weights=TRUE, plots=FALSE)
PlotBest(res=results, best=1, yall=data$D, xall=data$time)
