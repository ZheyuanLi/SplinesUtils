\name{predict.PiecePoly}

\alias{predict.PiecePoly}

\title{
``predict'' method for ``PiecePoly'' class.
}

\description{
The function evaluates piecewise polynomials or its derivatives at a vector of \eqn{x}-coordinate values.
}

\usage{
\method{predict}{PiecePoly}(object, newx, deriv = 0L, \dots)
}

\arguments{
  \item{object}{A ``PiecePoly'' object.}
  \item{newx}{A vector of x-coordinate values where the piecewise polynomials or its derivatives is evaluated.}
  \item{deriv}{Derivatives of the piecewise polynomials can be evaluated. Much interest lies in the 1st or 2nd derivative, so while a higher derivative can be computed, it is not performed by the function.}
  \item{\dots}{Not used by the method.}
}

\value{
A vector of evaluated \eqn{y}-coordinate values.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

