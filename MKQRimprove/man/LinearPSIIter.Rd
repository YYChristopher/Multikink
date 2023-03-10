% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LinearPSIIter.R
\name{psi.est}
\alias{psi.est}
\title{Change points estimation for multi-kink quantile regression}
\usage{
psi.est(x, y, z, psi.ini, tau, tol, max.iter, yta)
}
\arguments{
\item{x}{A numeric variable with change point}

\item{y}{A vector of response variable}

\item{z}{A vector of covariates (without change point)}

\item{psi.ini}{A initial value vector of change points}

\item{tau}{quantile level}

\item{tol}{iterative tolerance value, 1e-8 for default}

\item{max.iter}{the maximum iteration steps, 5000 for default}

\item{yta}{learning rate for iteration,0.1 for default}
}
\value{
A list with the elements
\item{psi.new}{The estimated change points.}
}
\description{
This function gives change points estimation for multi-kink quantile regression.
}
\author{
Yu Yue
}
\keyword{Change}
\keyword{estimation}
\keyword{points}
