library(tseries)
library(ggplot2)
library(roxygen2)
load_pkgload("C:/Users/user/Documents/Graduation Project/code2")
path <- "C:/Users/user/Documents/Graduation Project/code/Multikink R code and read data/data"
setwd(path)
data <- read.table("AirTemp.txt")
y <- data[,3]
x <- seq(1850,2020,0.0833)
tempmonth <- as.data.frame(cbind(x,y))
n <- length(y)
z <- NULL
tau1 <- 0.1
tau2 <- 0.3
tau3 <- 0.5
tau4 <- 0.7
tau5 <- 0.9
obj <- kinkTest(y, x, z, tau, NB = 500,sparsity="nid",bandwidth_type="Hall-Sheather")
fit1 <- mkqr.bea(y, x, z, tau1,Cn=log(n))
fit2 <- mkqr.bea(y, x, z, tau2,Cn=log(n))
fit3 <- mkqr.bea(y, x, z, tau3,Cn=log(n))
fit4 <- mkqr.bea(y, x, z, tau4,Cn=log(n))
fit5 <- mkqr.bea(y, x, z, tau5,Cn=log(n))
k <- length(fit1$psi.est)
windows()
p <- ggplot(data = tempmonth,mapping = aes(x = x,y = y))
X <- matrix(rep(x,k),nrow=n)
PSI1 <- matrix(rep(fit1$psi.est,rep(n,k)),ncol=k)
xzu1 <- cbind(1,x,pmax(X-PSI1,0))
yt1 <- xzu1%*%fit1$bet.est
d1 <- data.frame(x=x,y=yt1)
PSI2 <- matrix(rep(fit2$psi.est,rep(n,k)),ncol=k)
xzu2 <- cbind(1,x,pmax(X-PSI2,0))
yt2 <- xzu2%*%fit2$bet.est
d2 <- data.frame(x=x,y=yt2)
PSI3 <- matrix(rep(fit3$psi.est,rep(n,k)),ncol=k)
xzu3 <- cbind(1,x,pmax(X-PSI3,0))
yt3 <- xzu3%*%fit3$bet.est
d3 <- data.frame(x=x,y=yt3)
PSI4 <- matrix(rep(fit4$psi.est,rep(n,k)),ncol=k)
xzu4 <- cbind(1,x,pmax(X-PSI4,0))
yt4 <- xzu4%*%fit4$bet.est
d4 <- data.frame(x=x,y=yt4)
PSI5 <- matrix(rep(fit5$psi.est,rep(n,k)),ncol=k)
xzu5 <- cbind(1,x,pmax(X-PSI5,0))
yt5 <- xzu5%*%fit5$bet.est
d5 <- data.frame(x=x,y=yt5)
pic <- p + geom_point() + xlab("age(year)") + ylab("TempMonthly") +
  geom_line(data = d1,aes(x,yt1),color = 'red',linetype = 'longdash') +
  geom_line(data = d2,aes(x,yt2),color = 'blue',linetype = 'longdash') +
  geom_line(data = d3,aes(x,yt3),color = 'green',linetype = 'longdash') +
  geom_line(data = d4,aes(x,yt4),color = 'purple',linetype = 'longdash') +
  geom_line(data = d5,aes(x,yt5),color = 'orange',linetype = 'longdash')
pic

sink("TempBRISQ.txt")
cat("BRISQ,tau=0.1,Est1:","\n")
print(fit1)
cat("BRISQ,tau=0.3,Est2:","\n")
print(fit2)
cat("BRISQ,tau=0.5,Est3:","\n")
print(fit3)
cat("BRISQ,tau=0.7,Est4:","\n")
print(fit4)
cat("BRISQ,tau=0.9,Est5:","\n")
print(fit5)
sink()

