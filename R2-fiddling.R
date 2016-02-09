makedata <- function(n, b0=0.1, b1=2, sd.x = 1, sd.e = 0.1) {
    x <- rnorm(n, mean=10,sd=sd.x)
    e <- rnorm(n, sd=sd.e)
    y <- b0 + b1*x + e

    list(
        'y'          = y,
        'x'          = x
    )
}

doReg <- function(y, x) {
    X           <-  matrix(cbind(1, x), ncol=2)
    bHat        <-  (solve(t(X) %*% X)) %*% t(X) %*% y
    yHat        <-  X %*% bHat
    sigmaHatUb  <-  sum((y - (X %*% bHat))^2) / (length(y) - 2)
    stdResid    <-  (y - yHat) / sqrt(sigmaHatUb)
    r2          <- 1 - ((sum((y - yHat)^2)) / (sum((y - mean(y))^2)))
    r2adj       <- r2 - (1 - r2) * (2/(length(x) - 2 - 1))
    list(
        'y'          = y,
        'x'          = x,
        'bHat'       = bHat,
        'yHat'       = yHat,
        'sigmaHatUb' = sigmaHatUb,
        'stdResid'   = stdResid,
        'r2'         = r2,
        'r2adj'      = r2adj
    )
}






col1 <- adjustcolor('blue', alpha=0.2)
col2 <- adjustcolor('red', alpha=0.1)


r <- c()
radj <- c()

for (j in 0:100) {
    testdat <- makedata(n=200, sd.x = 2, sd.e = 1)
    ns <- c(5:length(testdat$x))
    for (i in 1:length(ns)) {
        reg <- doReg(testdat$y[1:ns[i]], testdat$x[1:ns[i]])
        r[(j*length(ns))+i]    <- reg$r2
        radj[(j*length(ns))+i] <- reg$r2adj
    }
}

par(mfrow=c(2,2))
plot(r~rep(ns,times=101), ylim=c(0,1), pch=21, col=NA, bg=col1, xlab="Sample size")
points(radj~rep(ns,times=101), pch=21, col=NA, bg=col2)
legend(x=(length(ns)*0.8), y=0.8,c("r2", "r2.adj"), cex=1.25, pch=21, col=NA, fill=c(col1,col2))
plot(density(r), type='l', lwd=4, col="blue")
lines(density(radj), lwd=4, col="red")
legend(x=0.6, y=(max(density(radj)$y *0.9)),c("r2", "r2.adj"), cex=1.25, lty=1, lwd=4, col=c("blue","red"))
plot(testdat$y~testdat$x, pch=21, bg="grey80")





par(mfrow=c(2,2))
r <- c()
radj <- c()
testdat <- makedata(n=300, sd.e = 1)
ns <- c(5:length(testdat$x))
for (i in 1:length(ns)) {
    reg <- doReg(testdat$y[1:ns[i]], testdat$x[1:ns[i]])
    r[i]    <- reg$r2
    radj[i] <- reg$r2adj
}


plot(r, pch=21, bg='grey80')
points(radj, pch=21, bg='red')
hist(r, breaks=40)
plot(testdat$y ~ testdat$x, pch=21, bg="grey80")
