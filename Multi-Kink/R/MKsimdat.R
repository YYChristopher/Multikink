#' Simulated data from the multi-kink quantile regression
#'
#' The function generates data from the multi-kink quantile regression model.
#'
#' The bent line quantile regression model:
#' Q_tau (Y|x,z) = [1, x, z,(x-\delta_1)_+,\cdots,(x-\delta_K)_+] * bet0. 
#' 
#' @rdname MKsimdat
#' @param n sample size.
#' @param bet0 the vector of real multi-kink regression coefficients.
#' @param psi0 the real location vector of the change points.
#' @param tau the quantile level.
#' @param etype type of error, where 
#'  etype = "normal" for N(0,1)
#'  etype = "student" for t(3)
#'  etype = "laplace" for LA(0,1)
#' @param scedast type of scedasticity, where
#'  scedast = "homo" for homoscedasticity
#'  scedast = "hetero" for heteroscedasticity
#'
#' @return A list with the elements
#' \item{y}{The response variable.}
#' \item{x}{The scalar covariate with change points.}
#' \item{z}{A vector of covariates.}
#'
#' @author Yu Yue 
#' @keywords MKsimdat
#' @importFrom stats runif rbinom rnorm qnorm  rt qt  rmutil 
#' @export
#'
#' @examples
#' /dontrun{
#' 
#'  ## simulated data
#'  n <- 1000
#'  bet <- c(1,1,1,-3,4)
#'  psi <- c(-1,2)
#'  tau <- 0.25
#'  data <- MKsimdat(n,bet,psi,tau1,etype = "normal",scedast = "homo")
#' }
#'


MKsimdat <- function(n,bet0,psi0,tau,etype=c("normal","student","laplace"),scedast=c("homo","hetero")){
  library(quantreg)
  library(rmutil)
  etype <- match.arg(etype,c("normal","student","laplace"))
  scedast <- match.arg(scedast,c("homo","hetero"))
  
  k <- length(psi0)
  x <- runif(n, -5, 5)
  z <- rnorm(n, 1, 1)
  X <- matrix(rep(x,k),nrow=n)
  PSI <- matrix(rep(psi0,rep(n,k)),ncol=k)
  XZ <- cbind(1,z,x,pmax((X-PSI),0))
  if(etype == "normal"){
    e <- rnorm(n,0,1) - qnorm(tau,0,1)
  }else if(etype == "student"){
    e <- rt(n,3) - qt(tau,3)
  }else if(etype == "laplace"){
    e <- rlaplace(n,0,1) - qlaplace(tau,0,1)
  }else{
    stop("Not a valid error type!")
  }
  
  if(scedast == "homo"){
    err <- e
  }else if(scedast == "hetero"){
    err <- e*(1 + 0.2*x)
  }else{
    stop("Not a valid scedasticity type!")
  }
  y <- as.vector(XZ %*% bet0) + err
  
  psi <- quantile(x,seq(0,1,l=(k+2)))[-c(1,k+2)]#Initiation of change points
  X <- matrix(rep(x,k),nrow=n)
  PSI <- matrix(rep(psi,rep(n,k)),ncol=k)
  pz <- ifelse(is.null(z), 0, ncol(as.matrix(z)))
  U <- pmax((X - PSI), 0)
  V <- ifelse((X > PSI), -1, 0)
  XZ <- cbind(x,z,U,V)
  obj <- rq(y~XZ, tau = tau,method ="br")
  
  bp.ini <- psi
  bet.ini <- obj$coefficients[1:(2+pz+k)]#initiation of coefficients
  h <- sd(e)*n^(-0.26)#initiation of bandwidth
  
  return(list(x=x,y=y,z=z,h=h,bp.ini=bp.ini,bet.ini=bet.ini))
}

