\name{print.PiecePoly}

\alias{print.PiecePoly}

\title{
``print'' method for ``PiecePoly'' class.
}

\description{
The function hides contents of a ``PiecePoly'' object when printing it. It calls \code{\link{summary.PiecePoly}} to export up to 6 piecewise polynomial equations as formatted strings and prints them onto standard output for a sample display.
}

\usage{
\method{print}{PiecePoly}(x, \dots)
}

\arguments{
  \item{x}{A ``PiecePoly'' object.}
  \item{\dots}{Not used by the method.}
}

\value{
The function returns nothing.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

