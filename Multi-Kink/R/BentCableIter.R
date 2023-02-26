#' @title BentCable method of parameters estimation for multi-kink quantile regression 
#' 
#' @description This function gives BentCable method of parameters estimation for multi-kink
#' quantile regression.
#'
#' @rdname BentCableIter
#' @param y A vector of response variable
#' @param x A numeric variable with change point
#' @param z A vector of covariates (without change point)
#' @param tau quantile level
#' @param h bandwidth for kernel function
#' @param bet.ini  A initial vector of regression coefficients
#' @param bp.ini  A initial value vector of change points
#' @param tol iterative tolerance value, 1e-8 for default
#' @param max.iter the maximum iteration steps, 5000 for default
#' @param yta learning rate for iteration,0.1 for default
#' @param a learning rate adjustment parameters (reduce learning rate),0.5 for default
#' @param b learning rate adjustment parameters (increase learning rate),2 for default
#'
#' @return A list with the elements
#' \item{psi}{The estimated change points.}
#' \item{coefficients}{The estimated regression coefficients with intercept.}
#' \item{h}{The bandwidth results after iterations.}
#' \item{iter}{The iteration steps of algorithm.}
#' \item{Loss.result}{End value of loss function.}
#' \item{Sn}{BentCable kink function of model}
#'
#' @author Yu Yue
#' @keywords BentCable
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
#' 
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

hDerive <- function(x,psi,h){
  k <- length(psi)
  n <- length(x)
  h.der <- matrix(0,n,k)
  for(j in 1:k){
    for (i in 1:n) {
      if(abs(x[i]-psi[j]) <= h){
        h.der[i,j] <- (h^2-(x[i]-psi[j])^2)/(4*h^2)
      }else{
        h.der[i,j] <- 0
      }
    }
  }
  
  return(h.der)
}

BentCableIter <- function(x,y,z,tau,h,bet.ini, bp.ini,tol=1e-8, max.iter=5000,yta,a,b){
  
  if(is.null(z)){
    z <- NULL
  }else{
    z <- as.matrix(z)
  }
  n <- length(y)
  k <- length(bp.ini)
  
  checkfun <- function(u){
    u*(tau-ifelse(u<0, 1, 0))#u<0ȡ1????ȡ0
  }
  
  wfun <- function(u, tau) tau-(u<0)#??Ӧ?????е?psi_\tau????
  
  CableModel <- function(x,z,psi,h,bet){
    X <- matrix(rep(x,k),nrow=n)
    PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
    Sn.ca <- CableProduce(x,psi,h)
    if(is.null(z)){
      xzu <- cbind(1,x,Sn.ca)
    }else{
      xzu <- cbind(1,x,z,Sn.ca)
    }
    return(xzu%*%bet)
  }
  
  psi.result <- psi.est(x,y,z,bp.ini,tau,tol,max.iter,yta)
  
  p <- ifelse(is.null(z), 0, ncol(as.matrix(z)))
  
  Sn <- CableProduce(x,psi.result,h)
  qn <- CableDerive(x,psi.result,h)
  hn <- hDerive(x,psi.result,h)
  
  Coef.new <- bet0 <- bet.ini
  h.new <- h.old <- h
  Intercept.old <- Intercept <- bet0[1]
  Coef.x.old <- Coef.x <- bet0[2]
  if(p>0){
    Coef.z.old <- Coef.z <- bet0[3:(2+p)]
    Coef.z.new <- rep(0,length(Coef.z.old))
  }else{
    Coef.z.old <- NULL
    Coef.z.new <- NULL
  }
  Coef.Ker.old <- Coef.Ker <- bet0[(3+p):(2+p+k)]
  Coef.Ker.new <- rep(0,length(Coef.Ker.old))
  
  iter <- 1
  
  while (iter <= max.iter) {
    PSI <- matrix(rep(psi.result,rep(n,k)),ncol=k)
    ModelEqua <- CableModel(x,z,psi.result,h,Coef.new)
    DiffValue <- as.vector(y-ModelEqua)
    DerivePart <- wfun(DiffValue,tau)
    Sn.new <- CableProduce(x,psi.result,h.new)
    qn.new <- CableDerive(x,psi.result,h.new)
    hn.new <- hDerive(x,psi.result,h.new)
    
    Intercept.new <- Intercept.old + yta*mean(DerivePart)
    Coef.x.new <- Coef.x.old + yta*mean(DerivePart*x)
    for(i in 1:k){
      Coef.Ker.new[i] <- Coef.Ker.old[i] + yta*mean(DerivePart*Sn.new[,i])
    }
    
    if(p > 0){
      for(i in 1:p){
        Coef.z.new[i] <- Coef.z.old[i] + yta*mean(DerivePart*z[,i])
      }
    }else{
      Coef.z.new <- NULL
    }
    
    for(i in 1:k){
      h.new <- h.old + yta*mean(DerivePart*hn.new[i])
    }
    
    if(p > 0){
      Coef.old <- c(Intercept.old,Coef.x.old,Coef.z.old,Coef.Ker.old)
      Coef.new <- c(Intercept.new,Coef.x.new,Coef.z.new,Coef.Ker.new)
    }else{
      Coef.old <- c(Intercept.old,Coef.x.old,Coef.Ker.old)
      Coef.new <- c(Intercept.new,Coef.x.new,Coef.Ker.new)
    }
    
    Loss.old <- mean(checkfun(y - CableModel(x,z,psi.result,h.old,Coef.old)))
    Loss.new <- mean(checkfun(y - CableModel(x,z,psi.result,h.new,Coef.new)))
    Sn.ca <- CableProduce(x,psi.result,h)
    
    if(Loss.old < Loss.new){
      yta <- a*yta
      iter <- iter + 1
    }else if(Loss.old - Loss.new > tol){
      Intercept.old <- Intercept.new
      Coef.x.old <- Coef.x.new
      Coef.Ker.old <- Coef.Ker.new
      Coef.z.old <- Coef.z.new
      h.old <- h.new
      iter <- iter + 1
      yta <- b*yta
    }else{
      break
    }
  }
  
  return(list(psi = psi.result,coefficients = Coef.new,h = h.new,iter = iter,Loss.result = Loss.new,Sn = Sn.ca))
}

