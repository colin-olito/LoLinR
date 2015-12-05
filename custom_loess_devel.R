#=============================================
# Colin Olito. Created 27/10/2015.
#
#  Development code for custom LOESS function for
#  finding most linear subset of O2 consumption data
#
#  NOTES: 





###################################
#  Step 1: OLS & WLS Regression machinery #
#
#  Notes: I'm stripping code from an old
#  stats class here so I can keep the
#  memory usage for the regression fits
#  to a minimum.
###################################

# Import test data #

data=read.csv("taxes.csv",header=TRUE)
data=data.frame(cbind(data$price, data$tax))
data=na.omit(data)
length(data[,1])  #should equal 107
is(data)
str(data)
######
# OLS #

# Design Matix #

y=data$tax
X=matrix(,nrow=length(y),ncol=2)
X[,1]=d
X[,2]=data$price
head(X)


# Fit Model #

b.hat = (solve(t(X) %*% X)) %*% t(X) %*% y
b.hat
y.hat = X %*% b.hat

hat.mat = X %*% (solve(t(X) %*% X)) %*% t(X) %*% y
sigma.hat = sum((y-(X %*% b.hat))^2) / length(y)
sigma.hat.unbiased = sum((y - (X %*% b.hat))^2) / (length(y) - 2)


# Diagnostic Plots #

par(mfrow=c(2,2))
std.resid = (y - y.hat)/sqrt(y.hat)
plot(std.resid ~ data$price,
       xlab="Price", ylab="Standardized Residuals", main="Std. Residuals ~ Price")
abline(h=0,col=2)
plot(std.resid ~ y.hat,
       xlab="Fitted Values",ylab="Standardized Residuals",main="Std. Residuals ~ Fitted Values")
abline(h=0,col=2)
qqnorm(std.resid, main="QQNorm plot of Std. Residuals")
qqline(std.resid,col=2)
hist(std.resid, xlab="Standardized Residuals", ylab="Density",breaks=10, main="Density Plot of Std. Residuals")



#  TriCube is Tukey tricubed weight function
#  accepts independent variable
#  centerpoint, and proportion
#  of data used to calculate weights

TriCube<-function(x, i, r) {
  h <- r/2
  x2 <- (x - x[i]) / sum(x - x[i])
  z <- abs(x2-x2[i])/h
  ifelse(z < 1, (1 - z^3)^3, 0)
}
testweights <- TriCube(X[,2], 50 , 0.2)


Skew <- function(x){
    n <- length(x)
    G <- (n/((n-1)*(n-2))) * sum(((x-mean(x))/sd(x))^3)
    return(G)
}


# Fit Model #
wb.hat = (solve(t(X) %*% diag(testweights) %*% X)) %*% t(X) %*% diag(testweights) %*% y
wb.hat
wy.hat = X %*% wb.hat
what.mat = X %*% (solve(t(X) %*% diag(testweights) %*% X)) %*% t(X) %*% diag(testweights) %*% y
wsigma.hat = sum((y - (X %*% wb.hat))^2)/length(y)
wsigma.hat.unbiased = sum((y-(X %*% wb.hat))^2)/(length(y) - 2)


# Diagnostic Plots #
wstd.resid = (y - y.hat)/sqrt(y.hat)
plot(wstd.resid ~ y,
       xlab="Price", ylab="Weighted Residuals", main="Weighted Residuals ~ Price")
abline(h=0,col=2)
plot(wstd.resid ~ wy.hat,
       xlab="Fitted Values", ylab="Weighted Residuals",main="Weighted Residuals ~ Fitted Values")
abline(h=0,col=2)
qqnorm(wstd.resid, main="QQNorm plot of Std. Residuals")
qqline(wstd.resid,col=2)
hist(wstd.resid, xlab="Standardized Residuals", ylab="Density",breaks=10, main="Density Plot of Std. Residuals")


# Calculate Skew of residuals
Skew(wstd.resid)

# Confirm that WLS is working correctly

wb.hat
coef(lm(data$tax ~ data$price, weights=testweights))

diff(range(y))












###################################
#  Step 2: Iteration machinery #
#
#  Notes: Basic idea is to iterate through
#           2 variables:
#              1) data range used for local regression (r).
#                  Start w/ minimum, iterate through larger
#                  proportions of data.
#
#              2) centerpoint for local regression (i), moving
#                 left to right.
#
#           at each step, need to record:
#             -- centerpoint (i)
#             -- range (r)
#             -- regression coefficient (beta)
#             -- r^2 value
#             -- skewness of residuals
###################################

head(y)
length(y)

res <- array(,dim=c(1,4))
alpha <- 0.2
h <- ceiling((alpha*length(y))/2)

while((2*h)+1 <= length(y)){

    for (i in 1:length(y)) {
        if (i-h < 1 | i+h > length(y))
            next
        else
            range <- c(y[i-h], y[i+h])

        res <- rbind(res,c(i,h,range))

        # FIT BLOCK GOES HERE  #
        
    }
    h <- h+1
}
res










###################################
#  Step 3: Complete Function #
#
#  Notes: Basic idea is to iterate through
#           2 variables:
#              1) data range used for local regression (r).
#                  Start w/ minimum, iterate through larger
#                  proportions of data.
#
#              2) centerpoint for local regression (i), moving
#                 left to right.
#
#           at each step, need to record:
#             -- centerpoint (i)
#             -- range (r)
#             -- regression coefficient (beta)
#             -- r^2 value
#             -- skewness of residuals
###################################
str(data)
xall <- as.vector(sort(as.numeric(as.character(data$X1))))
yall <- as.vector(sort(as.numeric(as.character(data$X2))))
length(x1)
length(y1)
is(x1)
is(y1)

alpha=0.2
yall <- data$A
xall <- data$time
head(data)

FindLocLin <- function(yall, xall, alpha, weights = TRUE, plots = TRUE) {

#  Initialize results  #
res <- array(,dim=c(1,7))

#  Initialize minimum data window  #
h <- ceiling((alpha*length(yall))/2)

#  Start Loops  #
while((2*h)+1 <= length(yall)) {
    for (i in 1:length(yall)) {

        if (i-h < 1 | i+h > length(yall))
            next
        else {
            y <- yall[c((i-h) : (i+h))]
            x <- xall[c((i-h) : (i+h))]
        }

        
##  FIT BLOCK  ##        

# Design Matrix #
X <- matrix(,nrow=length(y),ncol=2)
X[,1] <- 1
X[,2] <- x
        
        if(weights == FALSE){ # Use Ordinary Least Squares Regression #

            # Fit Model #
            b.hat = (solve(t(X) %*% X)) %*% t(X) %*% y
            y.hat = X %*% b.hat
            hat.mat = X %*% (solve(t(X) %*% X)) %*% t(X) %*% y
            sigma.hat = sum((y-(X %*% b.hat))^2) / length(y)
            sigma.hat.ub = sum((y - (X %*% b.hat))^2) / (length(y) - 2)
            std.resid = (y - y.hat)/sqrt(y.hat)
        }

        if(weights == TRUE){ # Use Weighted Least Squares Regression #

            # Calcualate weights #
            w <- TriCube(x=x, h=h)

            # Fit Model #
            b.hat = (solve(t(X) %*% diag(w) %*% X)) %*% t(X) %*% diag(w) %*% y
            y.hat = X %*% b.hat
            hat.mat = X %*% (solve(t(X) %*% diag(w) %*% X)) %*% t(X) %*% diag(w) %*% y
            sigma.hat = sum((y - (X %*% b.hat))^2)/length(y)
            sigma.hat.ub = sum((y-(X %*% b.hat))^2)/(length(y) - 2)
            std.resid = (y - y.hat)/sqrt(y.hat)
        }

        
##  Calculate Statistics of Interest  ##

    r2 <- 1 - ((sum((y - y.hat)^2)) / (sum((y - mean(y))^2)))
    alph <- (2*h) / length(yall)
    res <- rbind(res,c(i, h, alph, b.hat[1,], b.hat[2,], Skew(x = std.resid), r2))

    }
h <- h+1
}

##  Compile and Sort Results  ##
res <- data.frame(res[-1,])

#  Calculate combined metric (L) for linearity & fit  #
L <- ((min(abs(res[,6])) + abs(res[,6]))/sd(res[,6])) + ((max(res[,7]) - res[,7])/sd(res[,7]))
res <- cbind(res,L)

res <- data.frame(res[order(res[,8]),][1:25,])
names(res) <- c("i","h", "alpha", "b0","b1","skew","r2", "L")


    
################################
##  Plots to accompany best results  ##

if(plots==TRUE) {

pdf(file="testplots.pdf", height=15, width=15)
par(mfrow=c(5,5))

for(i in 1:nrow(res)) {
ytemp <- yall[c((res$i[i]-res$h[i]) : (res$i[i]+res$h[i]))]
xtemp <- xall[c((res$i[i]-res$h[i]) : (res$i[i]+res$h[i]))]

plot(yall ~ xall, pch=21, col='grey80', main=i)
points(ytemp ~ xtemp, pch=21, bg=1, col=2,ask=TRUE)
abline(coef=c(res$b0[i],res$b1[i]), col=2)
}
graphics.off()

dev.new()
hist(res$b1, breaks=25)

}


##  Return Results  ##    
return(res)

}  #*** END OF FUNCTION


test <- FindLocLin(yall = data$X1, xall = data$X2, alpha=0.2, weights=TRUE, plots=TRUE)
test


PlotBest(1,res=test,yall=yall, xall=xall,weights=TRUE)









###################################################
#  Step 4: Followup Function to inspect best linear regression  #
#
#  Notes: Accepts
#             -- DataFrame produced by FindLocLin, and
#             -- Choice for best regression (number between 1 and 25)
#             -- range (r)
#             -- regression coefficient (beta)
#             -- r^2 value
#             -- skewness of residuals
###################################################

PlotBest <- function(best, res, yall, xall, weights = TRUE) {

    yall <- as.vector(as.numeric(as.character(yall)))
    xall <- as.vector(as.numeric(as.character(xall)))

    y <- yall[c((res$i[best] - res$h[best]) : (res$i[best] + res$h[best]))]
    x <- xall[c((res$i[best] - res$h[best]) : (res$i[best] + res$h[best]))]

#  FIT BLOCK  #        
            # Design Matix #
            X <- matrix(,nrow=length(y),ncol=2)
            X[,1] <- 1
            X[,2] <- x
#            I=matrix(c(1,0,0,1),nrow=2, ncol=2)
        
        if(weights == FALSE){ # Use Ordinary Least Squares Regression #

            # Fit Model #
            b.hat = (solve(t(X) %*% X)) %*% t(X) %*% y
            y.hat = X %*% b.hat
            hat.mat = X %*% (solve(t(X) %*% X)) %*% t(X) %*% y
            sigma.hat = sum((y-(X %*% b.hat))^2) / length(y)
            sigma.hat.ub = sum((y - (X %*% b.hat))^2) / (length(y) - 2)
            std.resid = (y - y.hat)/sqrt(y.hat)
        }

        if(weights == TRUE){ # Use Weighted Least Squares Regression #

            # Calcualate weights #
            w <- TriCube(X[,2], (h+1) , h)

            # Fit Model #
            b.hat = (solve(t(X) %*% diag(w) %*% X)) %*% t(X) %*% diag(w) %*% y
            y.hat = X %*% b.hat
            hat.mat = X %*% (solve(t(X) %*% diag(w) %*% X)) %*% t(X) %*% diag(w) %*% y
            sigma.hat = sum((y - (X %*% b.hat))^2)/length(y)
            sigma.hat.ub = sum((y-(X %*% b.hat))^2)/(length(y) - 2)
            std.resid = (y - y.hat)/sqrt(y.hat)
        }

    
##  Residual Plots  ##

#pdf(file="residplots.pdf", height=10, width=10)
dev.new()
par(mfrow=c(2,2))
plot(std.resid ~ x,
       xlab="x", ylab="y", main="Std. Residuals ~ x")
abline(h=0,col=2)
plot(std.resid ~ y.hat,
       xlab="Fitted Values",ylab="Standardized Residuals",main="Std. Residuals ~ Fitted Values")
abline(h=0,col=2)
qqnorm(std.resid, main="QQNorm plot of Std. Residuals")
qqline(std.resid,col=2)
hist(std.resid, xlab="Standardized Residuals", ylab="Density",breaks=20, main="Density Plot of Std. Residuals")
#graphics.off()
    
##  Overall Regression Plot  ##
dev.new()
plot(yall ~ xall, pch=21, col='grey80', ask=TRUE,
                  main=expression(paste("Best Local Regression: ",b[o], " = ", b.hat[1,], " = ", b.hat[2,])))
points(y ~ x, pch=21, bg=1, col=2,ask=TRUE)
abline(coef=c(b.hat[1,],b.hat[2,]), col=1)
    
}  #END OF FUNCTION



PlotBest(18,res=test,yall=yall, xall=xall,weights=TRUE)







#  TriCube is Tukey tricubed weight function
#  accepts independent variable
#  centerpoint, and proportion
#  of data used to calculate weights


TriCube<-function(x, i, h) {
  x2 <- (x - x[i]) / sum(x - x[i])
  z <- abs(x2-x2[i])/h
  ifelse(z < 1, (1 - z^3)^3, 0)
}


Skew <- function(x){
    n <- length(x)
    G <- (n/((n-1)*(n-2))) * sum(((x-mean(x))/sd(x))^3)
    return(G)
}











###################################################
#  Step 5: Test with real 02 data  #
#
###################################################



#  Import Data  #
?read.csv
data <- read.csv("TestO2data.csv", header=TRUE, stringsAsFactors=FALSE)
head(data)
nrow(data)

results2 <- FindLocLin(yall=data$D, xall=data$time, alpha=0.2, weights=TRUE, plots=TRUE)
results2
PlotBest(res=results2, best=1, yall=data$A, xall=data$time, weights=FALSE)





















#############################################3
##  Legacy function, before using combined metric L


FindLocLin <- function(yall, xall, alpha, weights = TRUE, plots = TRUE) {

#  Initialize whole dataset  #
yall <- as.vector(as.numeric(as.character(yall)))
xall <- as.vector(as.numeric(as.character(xall)))

#  Initialize results  #
res <- array(,dim=c(1,7))

#  Initialize data window  #
h <- ceiling((alpha*length(yall))/2)

#  Start Loop  #
while((2*h)+1 <= length(yall)) {
    for (i in 1:length(yall)) {

        if (i-h < 1 | i+h > length(yall))
            next
        else {
            y <- yall[c((i-h) : (i+h))]
            x <- xall[c((i-h) : (i+h))]
        }

        
#  FIT BLOCK  #        
            # Design Matix #
            X <- matrix(,nrow=length(y),ncol=2)
            X[,1] <- 1
            X[,2] <- x
#            I=matrix(c(1,0,0,1),nrow=2, ncol=2)
        
        if(weights == FALSE){ # Use Ordinary Least Squares Regression #

            # Fit Model #
            b.hat = (solve(t(X) %*% X)) %*% t(X) %*% y
            y.hat = X %*% b.hat
            hat.mat = X %*% (solve(t(X) %*% X)) %*% t(X) %*% y
            sigma.hat = sum((y-(X %*% b.hat))^2) / length(y)
            sigma.hat.ub = sum((y - (X %*% b.hat))^2) / (length(y) - 2)
            std.resid = (y - y.hat)/sqrt(y.hat)
        }

        if(weights == TRUE){ # Use Weighted Least Squares Regression #

            # Calcualate weights #
            w <- TriCube(x, h)

            # Fit Model #
            b.hat = (solve(t(X) %*% diag(w) %*% X)) %*% t(X) %*% diag(w) %*% y
            y.hat = X %*% b.hat
            hat.mat = X %*% (solve(t(X) %*% diag(w) %*% X)) %*% t(X) %*% diag(w) %*% y
            sigma.hat = sum((y - (X %*% b.hat))^2)/length(y)
            sigma.hat.ub = sum((y-(X %*% b.hat))^2)/(length(y) - 2)
            std.resid = (y - y.hat)/sqrt(y.hat)
        }

# Calculate Statistics of Interest   #

r2 <- 1 - ((sum((y - y.hat)^2)) / (sum((y - mean(y))^2)))
alph <- (2*h) / length(yall)
res <- rbind(res,c(i, h, alph, b.hat[1,], b.hat[2,], Skew(x = std.resid), r2))
    }
    h <- h+1
}

# Compile and Sort Results   #
res <- res[-1,]
res <- data.frame(res[order(abs(res[,6]), -res[,7]),][1:25,])

#res <- data.frame(res[1:25,])
names(res) <- c("i","h", "alpha", "b0","b1","skew","r2")



    
################################
##  Plots to accompany best results  ##

if(plots==TRUE) {

pdf(file="testplots.pdf", height=15, width=15)
par(mfrow=c(5,5))

for(i in 1:nrow(res)) {
ytemp <- yall[c((res$i[i]-res$h[i]) : (res$i[i]+res$h[i]))]
xtemp <- xall[c((res$i[i]-res$h[i]) : (res$i[i]+res$h[i]))]

plot(yall ~ xall, pch=21, col='grey80', main=i)
points(ytemp ~ xtemp, pch=21, bg=1, col=2,ask=TRUE)
abline(coef=c(res$b0[i],res$b1[i]), col=2)
}
graphics.off()

dev.new()
hist(res$b1, breaks=25)
}


##  Return Results  ##    
return(res)

}  #*** END OF FUNCTION
