#' @title Change points estimation for multi-kink quantile regression 
#' 
#' @description This function gives change points estimation for multi-kink quantile regression.
#' 
#'
#' @rdname LinearPSIIter
#' @param x A numeric variable with change point
#' @param y A vector of response variable
#' @param z A vector of covariates (without change point)
#' @param psi.ini  A initial value vector of change points
#' @param tau quantile level
#' @param tol iterative tolerance value, 1e-8 for default
#' @param max.iter the maximum iteration steps, 5000 for default
#' @param yta learning rate for iteration,0.1 for default
#'
#' @return A list with the elements
#' \item{psi.new}{The estimated change points.}
#'
#' @author Yu Yue
#' @keywords Change points estimation
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
#'  psi.ini <- data$bp.ini
#'  tol <- 1e-8
#'  max.iter <- 5000
#'  yta <- 0.1
#'  psi.Result <- psi.est(x,y,z,psi.ini,tau,tol,max.iter,yta)
#' 
#' }
#'


psi.est <- function(x,y,z,psi.ini,tau,tol,max.iter,yta){
  
  if(is.null(z)){
    z <- NULL
  }else{
    z <- as.matrix(z)
  }
  
  k <- length(psi.ini)
  n <- length(y)
  psi.new <- psi.old <- psi.ini
  
  X <- matrix(rep(x,k),nrow=n)
  
  iter <- 1
  while(iter <= max.iter){
    PSI <- matrix(rep(psi.old,rep(n,k)),ncol=k)
    pz <- ifelse(is.null(z), 0, ncol(as.matrix(z)))
    
    if(missing(z) || is.null(z)){
      p <- 2; XREG <- matrix(x,nrow=n,ncol=1);xreg.names <- c("x")
    }else{
      XREG <- cbind(z,x); p <- ncol(XREG)+1
      if(p==3) xreg.names <- c("x","z") else xreg.names <- c(paste("x",1:ncol(x),sep=""),"z")#将自变量与协变量组成矩阵形式并对每一列进行命名
    }
    
    U <- pmax((X - PSI), 0)#取每个X-PSI与0之间的较大值(n个)，若X位于纽结点位置右侧返回X与纽结点位置的差值，若位于左侧直接返回0
    V <- ifelse((X > PSI), -1, 0)
    
    XZ <- cbind(x,z,U,V)
    rownames(XZ) <- NULL
    obj <- rq(y~XZ, tau = tau,  method ="br")
    
    beta.c <- coef(obj)[(2+pz+1):(2+pz+k)]#取beta为obj中XZU项的对应系数
    gamma.c <-  coef(obj)[(3+pz+k):(2+pz+2*k)]
    
    psi.new <- psi.old + yta*gamma.c/beta.c
    
    if(sum(gamma.c^2) < tol && sum((psi.new - psi.old)^2)<tol){
      break
    }else{
      psi.old <- psi.new
      iter <- iter + 1
    }
  }
  
  return(psi.new)
}