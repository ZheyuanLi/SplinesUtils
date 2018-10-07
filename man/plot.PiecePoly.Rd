\name{plot.PiecePoly}

\alias{plot.PiecePoly}

\title{
``plot'' method for ``PiecePoly'' class.
}

\description{
The function plots piecewise polynomials or its derivatives.
}

\usage{
\method{plot}{PiecePoly}(x, spread = 3L, deriv = 0L, show.knots = FALSE, \dots)
}

\arguments{
  \item{x}{A ``PiecePoly'' object.}
  \item{spread}{A graphical parameter determining how dense the evaluation grid is set up for plotting. Note that this is a multiplier factor. For a spline of degree \eqn{d}, there will be \code{(spread * d)} evenly spaced points between any two nearby knots.}
  \item{deriv}{Derivatives of the spline (here represented as piecewise polynomials) can be plotted. Much interest lies in the 1st and 2nd derivative, so while higher derivative can be computed, it is not performed by the functions.}
  \item{show.knots}{If \code{TRUE}, vertical dotted lines will be drawn through knots of the spline. Note that most spline constructor functions would place knots by quantiles of the \eqn{x}-coordinate values.}
  \item{\dots}{Not used by the method.}
}

\value{
The function invisibly returns \eqn{(x, y)} data used for plotting in a list.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

