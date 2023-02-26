library(tseries)
library(ggplot2)
library(roxygen2)
load_pkgload("C:/Users/user/Documents/Graduation Project/code2")
path <- "C:/Users/user/Documents/Graduation Project/code/Multikink R code and read data/data"
setwd(path)
data <- read.table("AirTemp.txt")
path1 <- "C:/Users/user/Documents/Graduation Project/code/Multikink R code and read data/BentCable"
setwd(path1)
source("BentCableIter2.R")
source("BentCableVarEst.R")
source("LinearPSIIter.R")
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
k <- 3

psi <- quantile(x,seq(0,1,l=(k+2)))[-c(1,k+2)]
# psi <- c(3)
X <- matrix(rep(x,k),nrow=n)
PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
pz <- ifelse(is.null(z), 0, ncol(as.matrix(z)))

if(missing(z) || is.null(z)){
  p <- 2; XREG <- matrix(x,nrow=n,ncol=1);xreg.names <- c("x")
}else{
  XREG <- cbind(z,x); p <- ncol(XREG)+1
  if(p==3) xreg.names <- c("x","z") else xreg.names <- c(paste("x",1:ncol(x),sep=""),"z")#将自变量与协变量组成矩阵形式并对每一列进行命名
}

U <- pmax((X - PSI), 0)#取每个X-PSI与0之间的较大值(n个)，若X位于纽结点位置右侧返回X与纽结点位置的差值，若位于左侧直接返回0
V <- ifelse((X > PSI), -1, 0)#对应U的位置为正值取-1，否则取0.U和V可以用来理解指代自变量x的值位于哪一段的折线上(哪两个纽结点之间)
#U指文献第三页2.4式当中的U值，而V指2.4式中的V值

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
XZ <- cbind(x,z,U,V)
rownames(XZ) <- NULL
obj <- rq(y~XZ, tau = tau3,  method ="br")

obj1 <- rq(y~XZ, tau = tau1,  method ="br")
obj2 <- rq(y~XZ, tau = tau2,  method ="br")
obj3 <- rq(y~XZ, tau = tau3,  method ="br")
obj4 <- rq(y~XZ, tau = tau4,  method ="br")
obj5 <- rq(y~XZ, tau = tau5,  method ="br")

bet.ini1 <- obj1$coefficients[1:(2+pz+k)]
bet.ini2 <- obj2$coefficients[1:(2+pz+k)]
bet.ini3 <- obj3$coefficients[1:(2+pz+k)]
bet.ini4 <- obj4$coefficients[1:(2+pz+k)]
bet.ini5 <- obj5$coefficients[1:(2+pz+k)]

bp.ini <- psi
bet.ini <- obj$coefficients[1:(2+pz+k)]
h <- sd(obj$residuals)*n^(-0.26)
yta <- 0.1
a <- 0.5
b <- 2

Est1 <- BentCableIter2(x,y,z,tau1,h,bet.ini1, bp.ini,tol=1e-8, max.iter=5000,yta,a,b)
VarEst1 <- BentCableVar(y,x,z,tau1,h,Est1$coefficients,Est1$psi)
Est2 <- BentCableIter2(x,y,z,tau2,h,bet.ini2, bp.ini,tol=1e-8, max.iter=5000,yta,a,b)
VarEst2 <- BentCableVar(y,x,z,tau2,h,Est2$coefficients,Est2$psi)
Est3 <- BentCableIter2(x,y,z,tau3,h,bet.ini3, bp.ini,tol=1e-8, max.iter=5000,yta,a,b)
VarEst3 <- BentCableVar(y,x,z,tau3,h,Est3$coefficients,Est3$psi)
Est4 <- BentCableIter2(x,y,z,tau4,h,bet.ini4, bp.ini,tol=1e-8, max.iter=5000,yta,a,b)
VarEst4 <- BentCableVar(y,x,z,tau4,h,Est4$coefficients,Est4$psi)
Est5 <- BentCableIter2(x,y,z,tau5,h,bet.ini5, bp.ini,tol=1e-8, max.iter=5000,yta,a,b)
VarEst5 <- BentCableVar(y,x,z,tau5,h,Est5$coefficients,Est5$psi)

fit1.coef <- Est1$coefficients
fit2.coef <- Est2$coefficients
fit3.coef <- Est3$coefficients
fit4.coef <- Est4$coefficients
fit5.coef <- Est5$coefficients

X <- matrix(rep(x,k),nrow=n)
PSI1 <- matrix(rep(Est1$psi,rep(n,k)),ncol=k)
xzu1 <- cbind(1,x,(X-PSI1)*pnorm((X-PSI1)/h))
yt1 <- xzu1%*%fit1.coef
d1 <- data.frame(x=x,y=yt1)
PSI2 <- matrix(rep(Est2$psi,rep(n,k)),ncol=k)
xzu2 <- cbind(1,x,(X-PSI2)*pnorm((X-PSI2)/h))
yt2 <- xzu2%*%fit2.coef
d2 <- data.frame(x=x,y=yt2)
PSI3 <- matrix(rep(Est3$psi,rep(n,k)),ncol=k)
xzu3 <- cbind(1,x,(X-PSI3)*pnorm((X-PSI3)/h))
yt3 <- xzu3%*%fit3.coef
d3 <- data.frame(x=x,y=yt3)
PSI4 <- matrix(rep(Est4$psi,rep(n,k)),ncol=k)
xzu4 <- cbind(1,x,(X-PSI4)*pnorm((X-PSI4)/h))
yt4 <- xzu4%*%fit4.coef
d4 <- data.frame(x=x,y=yt4)
PSI5 <- matrix(rep(Est5$psi,rep(n,k)),ncol=k)
xzu5 <- cbind(1,x,(X-PSI5)*pnorm((X-PSI5)/h))
yt5 <- xzu5%*%fit5.coef
d5 <- data.frame(x=x,y=yt5)

windows()
p <- ggplot(data = tempmonth,mapping = aes(x = x,y = y))
pic <- p + geom_point() + xlab("age(year)") + ylab("TempMonthly") +
  geom_line(data = d1,aes(x,yt1),color = 'red',linetype = 'longdash') +
  geom_line(data = d2,aes(x,yt2),color = 'blue',linetype = 'longdash') +
  geom_line(data = d3,aes(x,yt3),color = 'green',linetype = 'longdash') +
  geom_line(data = d4,aes(x,yt4),color = 'yellow',linetype = 'longdash') +
  geom_line(data = d5,aes(x,yt5),color = 'orange',linetype = 'longdash')

pic

sink("TempBentCable.txt")
cat("BentCable,tau=0.1,Est1:","\n")
print(Est1)
print(VarEst1)
cat("BentCable,tau=0.3,Est2:","\n")
print(Est2)
print(VarEst2)
cat("BentCable,tau=0.5,Est3:","\n")
print(Est3)
print(VarEst3)
cat("BentCable,tau=0.7,Est4:","\n")
print(Est4)
print(VarEst4)
cat("BentCable,tau=0.9,Est5:","\n")
print(Est5)
print(VarEst5)
sink()