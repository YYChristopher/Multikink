#' @title BSQreg method of parameters estimation for multi-kink quantile regression 
#' 
#' @description This function gives BSQreg method of parameters estimation for multi-kink
#' quantile regression.
#'
#' @rdname BSQregIter
#' @param y A vector of response variable
#' @param x A numeric variable with change point
#' @param z A vector of covariates (without change point)
#' @param tau quantile level
#' @param h bandwidth for kernel function
#' @param bet.ini  A initial vector of regression coefficients
#' @param bp.ini  A initial value vector of change points
#' @param alpha parameter of new loss function
#' @param tol iterative tolerance value, 1e-8 for default
#' @param max.iter the maximum iteration steps, 5000 for default
#' @param yta learning rate for iteration,0.1 for default
#' @param a learning rate adjustment parameters (reduce learning rate),0.5 for default
#' @param b learning rate adjustment parameters (increase learning rate),2 for default
#'
#' @return A list with the elements
#' \item{psi}{The estimated change points.}
#' \item{coefficients}{The estimated regression coefficients with intercept.}
#' \item{iter}{The iteration steps of algorithm.}
#' \item{Loss.result}{End value of loss function.}
#'
#' @author Yu Yue
#' @keywords BSQreg Boosting
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
#'  alpha <- 0.1
#'  tol <- 1e-8
#'  max.iter <- 5000
#'  yta <- 0.1
#'  a <- 0.5
#'  b <- 2
#'  Result <- BS.QReg(x,y,z,tau,h,bet.ini, bp.ini,alpha,tol, max.iter,yta,a,b)
#' 
#' }
#'




S.func <- function(x,alpha,tau){
  fun <- tau*x + alpha*log(1 + exp(-x/alpha))
  return(fun)
}

S.derive <- function(x,alpha,tau){
  deri <- tau - 1/(1+exp(x/alpha))
  return(deri)
}

reg <- function(X,y) {
  X <- qr(X)#QR factorization
  as.matrix(qr.coef(X,y))
}

BS.QReg <- function(x,y,z,tau,h,bet.ini, bp.ini,alpha,tol, max.iter,yta,a,b){
  
  if(is.null(z)){
    z <- NULL
  }else{
    z <- as.matrix(z)
  }
  n <- length(y)
  k <- length(bp.ini)
  
  p <- ifelse(is.null(z), 0, ncol(as.matrix(z)))
  
  KernelModel <- function(x,z,psi,h,bet){
    X <- matrix(rep(x,k),nrow=n)
    PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
    if(is.null(z)){
      xzu <- cbind(1,x,pmax(X-PSI,0))
    }else{
      xzu <- cbind(1,x,z,pmax(X-PSI,0))
    }
    return(xzu%*%bet)
  }
  
  psi.result <- psi.est(x,y,z,bp.ini,tau,tol,max.iter,yta)
  
  ols.derive <- function(x,z,psi,h,U){
    X <- matrix(rep(x,k),nrow=n)
    PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
    if(is.null(z)){
      xzu <- cbind(1,x,pmax(X-PSI,0))
    }else{
      xzu <- cbind(1,x,z,pmax(X-PSI,0))
    }
    bt <- reg(xzu,U)
    et = U - xzu%*%bt
    return(list(coef = bt,residual = et))
  }
  
  Coef.new <- bet0 <- bet.ini
  
  psi.new <- psi.old <- psi0 <- bp.ini
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
  Coef.new <- Coef.old <- c(Intercept.old,Coef.x.old,Coef.z.old,Coef.Ker.old)
  
  iter <- 1
  
  while(iter <= max.iter){
    
    PSI <- matrix(rep(psi.result,rep(n,k)),ncol=k)
    ModelEqua <- KernelModel(x,z,psi.result,h,Coef.old)
    DiffValue <- as.vector(y - ModelEqua)
    DerivePart <- S.derive(DiffValue,alpha,tau)
    
    Loss.old <- mean(S.func(DiffValue,alpha,tau))
    Coef.add <- ols.derive(x,z,psi.result,h,DerivePart)
    
    Coef.new <- Coef.old + yta*(Coef.add$coef)
    ModelEqua.new <- KernelModel(x,z,psi.result,h,Coef.new)
    DiffValue.new <- as.vector(y - ModelEqua.new)
    Loss.new <- mean(S.func(DiffValue.new,alpha,tau))
    
    if(Loss.old < Loss.new){
      yta <- a*yta
      iter <- iter + 1
    }else if(Loss.old - Loss.new > tol){
      Coef.old <- Coef.new
      iter <- iter + 1
      yta <- b*yta
    }else{
      break
    }
  }
  
  return(list(psi = psi.result,coefficients = Coef.new,iter = iter,Loss.result = Loss.new))
}