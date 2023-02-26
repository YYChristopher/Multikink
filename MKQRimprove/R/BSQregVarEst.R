#' @title BSQreg method of parameters stand error estimation for multi-kink quantile regression
#'
#' @description This function gives BSQreg method of parameters stand error estimation for multi-kink
#' quantile regression.
#'
#' @rdname BSQregVarEst
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
#' @keywords BSQreg stand error
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
#'  alpha <- 0.1
#'  tol <- 1e-8
#'  max.iter <- 5000
#'  yta <- 0.1
#'  a <- 0.5
#'  b <- 2
#'  Result <- BS.QReg(x,y,z,tau,h,bet.ini, bp.ini,alpha,tol, max.iter,yta,a,b)
#'  VarEst <- BSQVar(x,y,z,tau,h,Result$coefficients,Result$psi)
#' }
#'




BSQVar <- function(y,x,z,tau,h,bet,psi){
  n <- length(y)
  k <- length(psi)
  bet.psi <- bet[3:(2+k)]
  p <- ifelse(is.null(z), 0, ncol(as.matrix(z)))#协变量的列数
  XX <- matrix(rep(x,k),nrow=n)
  PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
  BET.psi <- matrix(rep(bet.psi,rep(n,k)),ncol=k)

  x.psi <- pmax((XX-PSI),0)
  M <- cbind(1,x,x.psi,z)

  if(p > 0){
    Q <- cbind(1,x,(XX-PSI)*pnorm((XX-PSI)/h),z,-BET.psi*(pnorm((XX-PSI)/h)+(XX-PSI)*dnorm((XX-PSI)/h)/h))
  }else{
    Q <- cbind(1,x,(XX-PSI)*pnorm((XX-PSI)/h),-BET.psi*(pnorm((XX-PSI)/h)+(XX-PSI)*dnorm((XX-PSI)/h)/h))
  }

  density_fun <- function(y,M,tau,bandwidth_type){#定义密度估计函数
    eps <- .Machine$double.eps^(2/3)#表示R中自带的最小可表示的正值的(2/3)次方作为初始eps

    if (bandwidth_type == "Bofinger") {
      bandwidth <- n^(-1/5) * ((9/2 * dnorm(qnorm(tau))^4)/(2 * qnorm(tau)^2 + 1)^2)^(1/5)
    } else if (bandwidth_type == "Chamberlain") {
      alpha <- 0.05
      bandwidth <- qnorm(1 - tau/2) * sqrt(alpha * (1 - alpha)/n)
    } else if (bandwidth_type == "Hall-Sheather") {
      alpha <- 0.05
      bandwidth <- n^(-1/3) * qnorm(1 - alpha/2)^(2/3) *
        ((3/2 * dnorm(qnorm(tau))^2)/(2 * qnorm(tau)^2 + 1))^(1/3)
    } else {#基于不同的带宽种类给定不同的带宽
      stop("Not a valid bandwith method!")
    }
    # Compute the density
    # Hendricks and Koenker (1992)
    bu <- suppressWarnings(fitted.values(rq(y~M[,-1],tau+bandwidth,method="br")))#fitted.values函数为摘录出模型的拟合值的函数
    bl <- suppressWarnings(fitted.values(rq(y~M[,-1],tau-bandwidth,method="br")))#此处为得到以分位数加减带宽为真正分位数进行分位数回归得到的拟合值结果，即为拟合值范围
    density <- pmax(0, 2 * bandwidth/((bu - bl) - eps))#得到对应的每个自变量位置的密度值

    return(diag(density))#类似于mkqr.fit当中的密度估计的形式
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
