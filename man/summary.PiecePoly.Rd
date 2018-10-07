\name{summary.PiecePoly}

\alias{summary.PiecePoly}

\title{
``summary'' method for ``PiecePoly'' class.
}

\description{
The function works through columns of the coefficient matrix in a \link{PiecePoly} object, exporting equation of each piecewise polynomial as a formatted string.
}

\usage{
\method{summary}{PiecePoly}(object, no.eqn = NULL, \dots)
}

\arguments{
  \item{object}{A ``PiecePoly'' object.}
  \item{no.eqn}{Number of piecewise polynomials to export. Normally not used by users (hence left as \code{NULL} and later set as the number of piecewise polynomials). A small value of 6 is internally used by \code{\link{print.PiecePoly}} for a sample display of a ``PiecePoly'' object.}
  \item{\dots}{Not used by the method.}
}

\value{
The function \emph{invisibly} returns a list of strings, where \eqn{x} is always used as the variable name.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

