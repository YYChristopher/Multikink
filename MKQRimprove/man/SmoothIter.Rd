% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ISQregIter.R, R/SmoothIter.R
\name{ISQReg}
\alias{ISQReg}
\alias{SmoothMKQR}
\title{Induced Smooth method of parameters estimation for multi-kink quantile regression}
\usage{
ISQReg(x, y, z, tau, bet.ini, bp.ini, tol, max.iter, yta, a, b)

SmoothMKQR(
  x,
  y,
  z,
  tau,
  h,
  bet.ini,
  bp.ini,
  tol = 1e-08,
  max.iter = 5000,
  yta = 0.1,
  a,
  b
)
}
\arguments{
\item{x}{A numeric variable with change point}

\item{y}{A vector of response variable}

\item{z}{A vector of covariates (without change point)}

\item{tau}{quantile level}

\item{bet.ini}{A initial vector of regression coefficients}

\item{bp.ini}{A initial value vector of change points}

\item{tol}{iterative tolerance value, 1e-8 for default}

\item{max.iter}{the maximum iteration steps, 5000 for default}

\item{yta}{learning rate for iteration,0.1 for default}

\item{a}{learning rate adjustment parameters (reduce learning rate),0.5 for default}

\item{b}{learning rate adjustment parameters (increase learning rate),2 for default}

\item{h}{bandwidth for kernel function}
}
\value{
A list with the elements
\item{psi}{The estimated change points.}
\item{coefficients}{The estimated regression coefficients with intercept.}
\item{Varest}{The estimated standard error of the regression coefficients.}
\item{psi.se}{The estimated standard error of the change points.}
\item{iter}{The iteration steps of algorithm.}
\item{Loss.result}{End value of loss function.}

A list with the elements
\item{psi}{The estimated change points.}
\item{coefficients}{The estimated regression coefficients with intercept.}
\item{iter}{The iteration steps of algorithm.}
\item{Loss.result}{End value of loss function.}
}
\description{
This function gives Induced Smooth method of parameters estimation for multi-kink
quantile regression.

This function gives SmoothKernel method of parameters estimation for multi-kink
quantile regression.
}
\author{
Yu Yue

Yu Yue
}
\keyword{Induced}
\keyword{Kernel}
\keyword{Smooth}
