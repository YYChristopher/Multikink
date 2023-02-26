#' @title SmoothKernel method of parameters stand error estimation for multi-kink quantile regression
#'
#' @description This function gives SmoothKernel method of parameters stand error estimation for multi-kink
#' quantile regression.
#'
#' @rdname SmoothVarEst
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
#' @keywords Smooth Kernel stand error
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
#'  h <- data$h
#'  bet.ini <- data$bet.ini
#'  bp.ini <- data$bp.ini
#'  tol <- 1e-8
#'  max.iter <- 5000
#'  yta <- 0.1
#'  a <- 0.5
#'  b <- 2
#'  Result <- SmoothMKQR(x,y,z,tau,h,bet.ini, bp.ini,tol, max.iter,yta,a,b)
#'  VarEst <- SmoothVar(x,y,z,tau,h,Result$coefficients,Result$psi)
#' }
#'




SmoothVar <- function(x,y,z,tau,h,bet,psi){

  n <- length(y)
  k <- length(psi)
  bet.psi <- bet[3:(2+k)]
  p <- ifelse(is.null(z), 0, ncol(as.matrix(z)))
  X <- matrix(rep(x,k),nrow=n)
  PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
  BET.psi <- matrix(rep(bet.psi,rep(n,k)),ncol=k)

  x.psi <- pmax((X-PSI),0)
  M <- cbind(1,x,x.psi,z)

  if(p > 0){
    Q <- cbind(1,x,(X-PSI)*pnorm((X-PSI)/h),z,-BET.psi*(pnorm((X-PSI)/h)+(X-PSI)*dnorm((X-PSI)/h)/h))
  }else{
    Q <- cbind(1,x,(X-PSI)*pnorm((X-PSI)/h),-BET.psi*(pnorm((X-PSI)/h)+(X-PSI)*dnorm((X-PSI)/h)/h))
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

  D.achieve <- function(Q,DD,tau){
    -tau*(1-tau)*t(Q)%*%DD%*%Q/n
  }

  C.achieve <- function(Q,tau){
    tau*(1-tau)*t(Q)%*%Q/n
  }

  DD <- density_fun(y,M,tau,bandwidth_type = "Hall-Sheather")
  D <- D.achieve(Q,DD,tau)
  C <- C.achieve(Q,tau)

  psi.Sigma <- n^(-1)*tau*(1-tau)*solve(D)%*%C%*%t(solve(D))
  bet.Sigma <- n^(-1)*tau*(1-tau)*solve(D)%*%C%*%t(solve(D))

  psi.se.vector <- as.vector(sqrt(diag(psi.Sigma)))
  bet.se.vector <- as.vector(sqrt(diag(bet.Sigma)))

  return(list(bet.se = bet.se.vector[1:(length(bet.se.vector)-k)],psi.se = psi.se.vector[(length(psi.se.vector)-k+1):length(psi.se.vector)]))

}
