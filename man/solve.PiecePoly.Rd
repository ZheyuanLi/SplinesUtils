\name{solve.PiecePoly}

\alias{solve.PiecePoly}

\title{
``solve'' method for ``PiecePoly'' class.
}

\description{
The function solves equation \eqn{f(x) = b} and \eqn{f'(x) = b}, where \eqn{f(x)} is a spline represented as piecewise polynomials. The function can be used for finding roots and extrema of \eqn{f(x)} for example.
}

\usage{
\method{solve}{PiecePoly}(a, b = 0, deriv = 0L, \dots)
}

\arguments{
  \item{a}{A ``PiecePoly'' object representing a spline \eqn{f(x)} as piecewise polynomials.}
  \item{b}{The RHS of equations \eqn{f(x) = b} and \eqn{f'(x) = b}.}
  \item{deriv}{Derivatives applied to \eqn{f(x)}. Only \code{deriv = 0} and \code{deriv = 1} are implemented.}
  \item{\dots}{Not used by the method.}
}

\value{
A vector of \eqn{x}-coordinate values that solve the equation.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

