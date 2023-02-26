#' @title BentCable method of parameters stand error estimation for multi-kink quantile regression 
#' 
#' @description This function gives BentCable method of parameters stand error estimation for multi-kink
#' quantile regression.
#'
#' @rdname BentCableVarEst
#' @param y A vector of response variable
#' @param x A numeric variable with change point
#' @param z A vector of covariates (without change point)
#' @param tau quantile level
#' @param h bandwidth for kernel function
#' @param bet  Coefficients estimation results of multi-kink quantile regression
#' @param psi  Change points estimation of multi-kink quantile regression
#'
#' @return A list with the elements
#' \item{bet.se}{The estimated standard error of the regression coefficients.}
#' \item{psi.se}{The estimated standard error of the change points.}
#'
#' @author Yu Yue
#' @keywords BentCable stand error
#' 
#' @importFrom quantreg rq
#' @export
#'
#' @example {
#' 
#' ##simulated data
#'  n <- 1000
#'  bet <- c(1,1,1,-3,4)
#'  psi <- c(-1,2)
#'  tau <- 0.25
#'  data <- MKsimdat(n,bet,psi,tau1,etype = "normal",scedast = "homo")
#'  y <- data$y
#'  x <- data$x
#'  z <- data$z
#'  h <- data$h
#'  bet.ini <- data$bet.ini
#'  bp.ini <- data$bp.ini
#'  tol <- 1e-8
#'  max.iter <- 5000
#'  yta <- 0.1
#'  a <- 0.5
#'  b <- 2
#'  Result <- BentCableIter(x,y,z,tau,h,bet.ini, bp.ini,tol, max.iter,yta,a,b)
#'  VarEst <- BentCableVar(x,y,z,tau,h,Result$coefficients,Result$psi)
#' }
#'


CableProduce <- function(x,psi,h){
  k <- length(psi)
  n <- length(x)
  Sn <- matrix(0,n,k)
  for(j in 1:k){
    for(i in 1:n){
      if(x[i]<psi[j]-h){
        Sn[i,j] <- 0
      }else if(x[i]>psi[j]+h){
        Sn[i,j] <- x[i]-psi[j]
      }else{
        Sn[i,j] <- (x[i]-psi[j]+h)^2/(4*h)
      }
    }
  }
  
  return(Sn)
}

CableDerive <- function(x,psi,h){
  k <- length(psi)
  n <- length(x)
  qn <- matrix(0,n,k)
  for(j in 1:k){
    for(i in 1:n){
      if(x[i]<psi[j]-h){
        qn[i,j] <- 0
      }else if(x[i]>psi[j]+h){
        qn[i,j] <- -1
      }else{
        qn[i,j] <- -(x[i]-psi[j]+h)/(2*h)
      }
    }
  }
  
  return(qn)
}


BentCableVar <- function(y,x,z,tau,h,bet,psi){
  n <- length(y)
  k <- length(psi)
  bet.psi <- bet[3:(2+k)]
  p <- ifelse(is.null(z), 0, ncol(as.matrix(z)))#协??量??????
  XX <- matrix(rep(x,k),nrow=n)
  PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
  BET.psi <- matrix(rep(bet.psi,rep(n,k)),ncol=k)
  
  x.psi <- pmax((XX-PSI),0)
  M <- cbind(1,x,x.psi,z)
  Sn <- CableProduce(x,psi,h)
  qn <- CableDerive(x,psi,h)
  
  if(p > 0){
    H <- cbind(1,x,Sn,z,BET.psi*qn)
  }else{
    H <- cbind(1,x,Sn,BET.psi*qn)
  }
  
  
  density_fun <- function(y,M,tau,bandwidth_type){#?????芏裙兰坪???
    eps <- .Machine$double.eps^(2/3)#??示R???源?????小?杀?示????值??(2/3)?畏???为??始eps
    
    if (bandwidth_type == "Bofinger") {
      bandwidth <- n^(-1/5) * ((9/2 * dnorm(qnorm(tau))^4)/(2 * qnorm(tau)^2 + 1)^2)^(1/5)
    } else if (bandwidth_type == "Chamberlain") {
      alpha <- 0.05
      bandwidth <- qnorm(1 - tau/2) * sqrt(alpha * (1 - alpha)/n)
    } else if (bandwidth_type == "Hall-Sheather") {
      alpha <- 0.05
      bandwidth <- n^(-1/3) * qnorm(1 - alpha/2)^(2/3) *
        ((3/2 * dnorm(qnorm(tau))^2)/(2 * qnorm(tau)^2 + 1))^(1/3)
    } else {#???诓?同?拇?????????????同?拇???
      stop("Not a valid bandwith method!")
    }
    # Compute the density
    # Hendricks and Koenker (1992)
    bu <- suppressWarnings(fitted.values(rq(y~M[,-1],tau+bandwidth,method="br")))#fitted.values????为摘录??模?偷?????值?暮???
    bl <- suppressWarnings(fitted.values(rq(y~M[,-1],tau-bandwidth,method="br")))#?舜?为?玫??苑?位???蛹?????为??????位?????蟹?位???毓??玫???????值????????为????值??围
    density <- pmax(0, 2 * bandwidth/((bu - bl) - eps))#?玫???应??每???员?量位?玫??芏?值
    
    return(diag(density))#??????mkqr.fit???械??芏裙兰频???式
  }
  
  D.achieve <- function(H,DD,tau){
    -tau*(1-tau)*t(H)%*%DD%*%H/n
  }
  
  C.achieve <- function(H,tau){
    tau*(1-tau)*t(H)%*%H/n
  }
  
  DD <- density_fun(y,M,tau,bandwidth_type = "Hall-Sheather")
  D <- D.achieve(H,DD,tau)
  C <- C.achieve(H,tau)
  
  psi.Sigma <- n^(-1)*tau*(1-tau)*solve(D)%*%C%*%t(solve(D))
  bet.Sigma <- n^(-1)*tau*(1-tau)*solve(D)%*%C%*%t(solve(D))
  
  psi.se.vector <- as.vector(sqrt(diag(psi.Sigma)))
  bet.se.vector <- as.vector(sqrt(diag(bet.Sigma)))
  
  return(list(bet.se = bet.se.vector[1:(length(bet.se.vector)-k)],psi.se = psi.se.vector[(length(psi.se.vector)-k+1):length(psi.se.vector)]))
  
}