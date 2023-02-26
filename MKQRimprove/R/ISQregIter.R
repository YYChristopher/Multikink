#' @title Induced Smooth method of parameters estimation for multi-kink quantile regression
#'
#' @description This function gives Induced Smooth method of parameters estimation for multi-kink
#' quantile regression.
#'
#' @rdname SmoothIter
#' @param x A numeric variable with change point
#' @param y A vector of response variable
#' @param z A vector of covariates (without change point)
#' @param tau quantile level
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
#' \item{Varest}{The estimated standard error of the regression coefficients.}
#' \item{psi.se}{The estimated standard error of the change points.}
#' \item{iter}{The iteration steps of algorithm.}
#' \item{Loss.result}{End value of loss function.}
#'
#' @author Yu Yue
#' @keywords Induced Smooth
#'
#' @importFrom quantreg rq
#' @export
#'
#' @example
#'
#' \dontrun{
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
#'  bet.ini <- data$bet.ini
#'  bp.ini <- data$bp.ini
#'  tol <- 1e-8
#'  max.iter <- 5000
#'  yta <- 0.1
#'  a <- 0.5
#'  b <- 2
#'  Result <- ISQReg(x,y,z,tau,bet.ini, bp.ini,tol, max.iter,yta,a,b)
#'
#' }
#'




ISQReg <- function(x,y,z,tau,bet.ini, bp.ini,tol, max.iter,yta,a,b){
  #bet.ini顺序:1,x,z,pmax((X-PSI),0)

  n <- length(y)
  k <- length(bp.ini)
  X <- matrix(rep(x,k),nrow=n)
  PSI <- matrix(rep(bp.ini,rep(n,k)),ncol=k)

  if(is.null(z)){
    z <- NULL
    M <- cbind(1,x,pmax((X-PSI),0))
  }else{
    z <- as.matrix(z)
    M <- cbind(1,x,z,pmax((X-PSI),0))
  }

  pz <- ifelse(is.null(z), 0, ncol(as.matrix(z)))

  checkfun <- function(u,tau){
    u*(tau-ifelse(u<0, 1, 0))
  }

  LinearModel <- function(x,z,psi,bet){
    X <- matrix(rep(x,k),nrow=n)
    PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
    if(is.null(z)){
      xzu <- cbind(1,x,pmax((X-PSI),0))
    }else{
      xzu <- cbind(1,x,z,pmax((X-PSI),0))
    }
    return(xzu%*%bet)
  }

  psi.result <- psi.est(x,y,z,bp.ini,tau,tol,max.iter,yta)

  m <- k + 2 + pz
  H <- diag(x = 1/n,nrow = m,ncol = m)

  U.produce <- function(x,y,z,bet,psi,tau,H){
    X <- matrix(rep(x,k),nrow=n)
    PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
    if(is.null(z)){
      pwt <- cbind(1,x,pmax((X-PSI),0))
    }else{
      pwt <- cbind(1,x,z,pmax((X-PSI),0))
    }
    U.trans <-matrix(rep(tau - pnorm((pwt%*%bet - y)/sqrt(diag(pwt%*%H%*%t(pwt)))),ncol(pwt)),nrow = n)
    U <- as.vector(colMeans(pwt*U.trans))
    return(U)
  }

  A.produce <- function(x,y,z,bet,psi,tau,H){
    X <- matrix(rep(x,k),nrow=n)
    PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
    if(is.null(z)){
      pwt <- cbind(1,x,pmax((X-PSI),0))
    }else{
      pwt <- cbind(1,x,z,pmax((X-PSI),0))
    }
    ratio <- as.vector(-(dnorm((pwt%*%bet - y)/sqrt(diag(pwt%*%H%*%t(pwt))))/sqrt(diag(pwt%*%H%*%t(pwt)))))
    Ratio <- diag(x = ratio,nrow = n,ncol = n)
    A <- as.matrix(t(pwt)%*%Ratio%*%pwt)
    return(A)
  }

  density_fun <- function(y,M,tau,bandwidth_type){
    eps <- .Machine$double.eps^(2/3)

    if (bandwidth_type == "Bofinger") {
      bandwidth <- n^(-1/5) * ((9/2 * dnorm(qnorm(tau))^4)/(2 * qnorm(tau)^2 + 1)^2)^(1/5)
    } else if (bandwidth_type == "Chamberlain") {
      alpha <- 0.05
      bandwidth <- qnorm(1 - tau/2) * sqrt(alpha * (1 - alpha)/n)
    } else if (bandwidth_type == "Hall-Sheather") {
      alpha <- 0.05
      bandwidth <- n^(-1/3) * qnorm(1 - alpha/2)^(2/3) *
        ((3/2 * dnorm(qnorm(tau))^2)/(2 * qnorm(tau)^2 + 1))^(1/3)
    } else {
      stop("Not a valid bandwith method!")
    }
    # Compute the density
    # Hendricks and Koenker (1992)
    bu <- suppressWarnings(fitted.values(rq(y~M[,-1],tau+bandwidth,method="br")))
    bl <- suppressWarnings(fitted.values(rq(y~M[,-1],tau-bandwidth,method="br")))
    density <- pmax(0, 2 * bandwidth/((bu - bl) - eps))

    return(diag(density))
  }

  Coef.old <- Coef.new <- bet0 <- bet.ini
  H.new <- H.old <- H

  iter <- 1

  while(iter <= max.iter){
    PSI <- matrix(rep(psi.result,rep(n,k)),ncol=k)
    ModelEqua <- LinearModel(x,z,psi.result,Coef.old)
    DiffValue <- as.vector(y - ModelEqua)
    Loss.old <- mean(checkfun(DiffValue,tau))

    U.iter <- U.produce(x,y,z,Coef.old,psi.result,tau,H.old)
    A.iter <- A.produce(x,y,z,Coef.old,psi.result,tau,H.old)
    Coef.new <- Coef.old + solve(-A.iter)%*%U.iter

    bb <- Coef.new[3:(2+k)]
    BB <- matrix(rep(bb),nrow=n,ncol=k,byrow=T)
    DD <- density_fun(y,M,tau,bandwidth_type = "Hall-Sheather")
    ht <- as.matrix(cbind(1,x,pmax((X - PSI),0),z,BB*ifelse(X>PSI,-1,0)))
    Cn <- tau*(1-tau)/n*t(ht)%*%ht
    Dn <- t(ht)%*%DD%*%ht/n
    Dn.inv <- solve(Dn+1e-8)
    Gamma <- n^(-1) * Dn.inv %*% Cn %*% Dn.inv
    H.new <- Gamma[1:(2+pz+k),1:(2+pz+k)]
    psi.var <- diag(Gamma)[(3+pz+k):(2+pz+2*k)]

    ModelEqua.new <- LinearModel(x,z,psi.result,Coef.new)
    DiffValue.new <- as.vector(y - ModelEqua.new)
    Loss.new <- mean(checkfun(DiffValue.new,tau))

    if(Loss.old < Loss.new){
      yta <- a*yta
      iter <- iter + 1
    }else if(Loss.old - Loss.new > tol){
      Coef.old <- Coef.new
      H.old <- H.new
      iter <- iter + 1
      yta <- b*yta
    }else{
      break
    }
  }

  return(list(psi = psi.result,coefficients = Coef.new,Varest = sqrt(diag(H.new)),psi.se = sqrt(psi.var),iter = iter,Loss.result = Loss.new))

}
