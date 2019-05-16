# title: "Prediction of the thickness of the sedimentary value"
# author: "Jingmei Lu"
# date: "May, 2014"

#load package
library(tseries)

# load data
sedimentary <- read.table("D:\\SJSU\\Math265\\project\\data\\deposits.txt")
head(sedimentary)
dim(sedimentary)
sedimentary <- sedimentary$V1
# check missing value
summary(sedimentary)
sum(is.na(sedimentary))
# no missing value

# visualize data
plot(sedimentary, type = 'b')

# data transformation
# log transformation
sum(sedimentary <= 0) # all value > 0
plot(log(sedimentary), type = 'b')
adf.test(log(sedimentary))
# p value is 0.05595, non-stationary

# difference data
data.trans <- diff(log(sedimentary), lag = 1, differences = 1)
plot(data.trans, type = 'b')
adf.test(data.trans)
# p value is 0.01, stationary

# generate acf and pacf plot
par(mfrow=c(2,1))
acf(data.trans, type = c("correlation"), ylab="Sample ACF")
acf(data.trans, type = c("partial"), ylab="Sample PACF")
par(mfrow=c(1,1))
# The sample autocorrelation cut off after lag1. And the sample partial 
# autocorrelation might cuts off after lag5. These two plots might imply 
# an ARIMA(5, 1, 1) model. 

# create periodogram plot of the differenced log data
spec.pgram(data.trans, taper = .1)
# The plot might have one peak at f=0.5, which indicates an AR(1) process. 
# There might be two dips at about frequency 0 and 0.2, which indicate MA(1) or MA(3).
# So the candidate model might be ARIMA(5, 1, 1) , ARIMA(1, 1, 1), ARIMA(1, 1, 3).

# separate data into train and test set
n <- length(sedimentary)
train <- log(sedimentary[1:564])
test <- log(sedimentary[565:624])

#  Maybe fit "all" ARIMA models and look at AIC
my.aic.values <- matrix(0, nrow=10, ncol=10)

for(p in 1:8)
  for(q in 1:6)  
  {
    my.fit <- arima(train, order = c(p-1, 1, q-1))
    my.aic.values[p,q] <- my.fit$aic
  }
my.aic.values
# from AIC value, choose candidate models as follow
# arima(0,1,2), arima(0,1,3), arima(1,1,1), arima(1,1,2), arima(1,1,3)
# arima(2,1,1), arima(3,1,3), arima(5,1,1), arima(6,1,3)

# check mse, estimated coeffient and prediction of all candidate models
fit.1 <- arima(train, order = c(0,1,2))
fit.2 <- arima(train, order = c(0,1,3))
fit.3 <- arima(train, order = c(1,1,1))
fit.4 <- arima(train, order = c(1,1,2))
fit.5 <- arima(train, order = c(1,1,3))
fit.6 <- arima(train, order = c(2,1,1))
fit.7 <- arima(train, order = c(3,1,3))
fit.8 <- arima(train, order = c(5,1,1))
fit.9 <- arima(train, order = c(6,1,3))
fit.auto <- auto.arima(log(sedimentary))
# by reviewing significance of estimated coefficients of all models, 
# only estimated coefficients of arima(0, 1, 2), arima(1,1,1) are all significant.

# prediction of test data
pred.1 <- predict(fit.1, n.ahead = 60, se.fit = TRUE)
pred.3 <- predict(fit.3, n.ahead = 60, se.fit = TRUE)

library(DMwR)
acc1 <- regr.eval(test, pred.1$pred)
acc3 <- regr.eval(test, pred.3$pred)
# By comparing the MSE of train data, aic and MSE of test data, the final model
# is arima(1,1,1)

# put original data, prediction data of fit.1 and fit.3 in one plot
pred.1.orig<- exp(pred.1$pred)
pred.3.orig<- exp(pred.3$pred)
plot(549:564, sedimentary[549:564], xlim=c(549,624), ylim=c(1, 180),  type="b")
lines(565:624, pred.1.orig, type="b", col='green')
lines(565:624, pred.3.orig, type="b", col='red')
points(565:624, sedimentary[565:624])
# fit.1 and fit.3 almost generate same predictions.

# create spectrum plot
##  A function for the TRUE SPECTRUM (on the decibel scale) 
##  of an ARMA process.

my.spectrum <- function(phi.of.b, theta.of.b, variance=1)
{
  p <- length(phi.of.b)
  q <- length(theta.of.b)
  
  omega <- seq(from=0, to=pi, by=.001)
  
  phi.of.e.minus.i.omega <- 1
  phi.of.e.i.omega <- 1
  
  if(p>1)
  {   for(i in 2:p)
  {
    phi.of.e.minus.i.omega <-  phi.of.e.minus.i.omega + phi.of.b[i]*exp(complex(imaginary = -(i-1))*omega)
    phi.of.e.i.omega <-  phi.of.e.i.omega + phi.of.b[i]*exp(complex(imaginary = (i-1))*omega)
  }
  }
  
  theta.of.e.minus.i.omega <- 1
  theta.of.e.i.omega <- 1
  
  if(q>1)
  {
    for(i in 2:q)
    {
      theta.of.e.minus.i.omega <-  theta.of.e.minus.i.omega + theta.of.b[i]*exp(complex(imaginary = -(i-1))*omega)
      theta.of.e.i.omega <-  theta.of.e.i.omega + theta.of.b[i]*exp(complex(imaginary = (i-1))*omega)
    }
  }
  
  my.spectrum <- (variance/(2*pi))*Re(theta.of.e.minus.i.omega*theta.of.e.i.omega/(phi.of.e.minus.i.omega*phi.of.e.i.omega))
  
  plot(omega, 10*log10(my.spectrum), ylab="spectrum (in decibels)", type="l")   
}

# put sample periodogram and theoretical spetrum plot together
par(mfrow=c(2,1))
spec.pgram(data.trans, taper = .1)  ##  tapr = .1 by default
my.spectrum(phi.of.b=c(1, -fit.3$coef[1]), theta.of.b=c(1, fit.3$coef[2]), 
            variance=fit.3$sigma2)
par(mfrow=c(1,1))

# model adequcy check
tsdiag(fit.3)
qqnorm(fit.3$residuals)

# generate prediction value and 95% prediction interval
lower.limits  <- pred.3$pred - 1.96*pred.3$se
upper.limits <- pred.3$pred + 1.96*pred.3$se
# Put back on the "Y" scale of original data
pred.orig<- exp(pred.3$pred)
lower.limit.orig <- exp(lower.limits)
upper.limit.orig <- exp(upper.limits)

cbind(lower.limit.orig, pred.orig, upper.limit.orig)

plot(549:564, sedimentary[549:564], xlim=c(549,624), ylim=c(1, 180),  type="b")
lines(565:624, pred.orig, type="b", col='red')
lines(565:624, lower.limit.orig, type="l", col='red')
lines(565:624, upper.limit.orig, type="l", col='red')
points(565:624, sedimentary[565:624], col='green')

# forecast next 12 new values
fit.final <- arima(log(sedimentary), order = c(1,1,1)) # use entire data to fit model
forecast.final <- predict(fit.final, n.ahead = 12, se.fit = TRUE)
forecast.orig <- exp(forecast.final$pred)
upper <- forecast.final$pred + 1.96*forecast.final$se
lower <- forecast.final$pred - 1.96*forecast.final$se
upper.orig <- exp(upper)
lower.orig <- exp(lower)
cbind(lower.orig, forecast.orig, upper.orig)
plot(600:624, sedimentary[600:624], xlim=c(600,640), ylim=c(1, 180),  type="b")
lines(625:636, forecast.orig, type="b", col='red')
lines(625:636, lower.orig, type="l", col='red')
lines(625:636, upper.orig, type="l", col='red')
