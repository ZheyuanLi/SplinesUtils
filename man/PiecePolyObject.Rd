\name{PiecePolyObject}

\alias{PiecePolyObject}
\alias{PiecePoly}
\alias{SplinesUtils}

\title{
Ojbect of ``PiecePoly'' class.
}

\description{
\code{\link{CubicInterpSplineAsPiecePoly}}, \code{\link{SmoothSplineAsPiecePoly}} and \code{\link{RegBsplineAsPiecePoly}} return objects of ``PiecePoly'' class, for which S3 methods \code{\link[=summary.PiecePoly]{summary}}, \code{\link[=print.PiecePoly]{print}}, \code{\link[=plot.PiecePoly]{plot}}, \code{\link[=predict.PiecePoly]{predict}} and \code{\link[=solve.PiecePoly]{solve}} are defined.
}

\details{
A spline \eqn{f(x)} of degree \eqn{d} and \eqn{k} knots \eqn{x_1,\ x_2, \ldots,\ x_k}, is number of piecewise polynomials of degree \eqn{d} that are \eqn{(d - 1)} times differentiable at the knots. For the \eqn{i}-th piece on interval \eqn{[x_i,\ x_{i + 1}]}, the polynomial \eqn{f_i(x)} is often parametrized in a shifted form \deqn{f_i(x) = b_{0i} + b_{1i}(x - x_i) + b_{2i}(x - x_i) ^ 2 + \ldots + b_{di}(x - x_i) ^ d.} This form is appealing because \eqn{b_{0i}} is just the value of the spline at knot \eqn{x_i}. It also facilitates integrating the spline. However, an unshifted version can also be exported \deqn{f_i(x) = a_{0i} + a_{1i}x + a_{2i}x ^ 2 + \ldots + a_{di}x ^ d.}

Piecewise polynomial representation is conveniently used in interpolation problems by \code{stats::spline} and \code{stats::splinefun}. However, B-splines representation is often used in regression and smoothing problems, for example by \code{splines::bs}, \code{splines::ns} and \code{stats::smooth.spline}. Coefficients for B-spline basis are difficult to interpret, so it might be helpful to reparametrize a fitted regression spline or smoothing spline into its conventional piecewise polynomial form.

Package \link{SplinesUtils} conducts such reparametrization by exact least squares fitting. \eqn{d + 1} evenly spaced points are taken on \eqn{[x_i,\ x_{i + 1}]}, where spline values are predicted. Given those \eqn{(x,\ y)} data, a polynomial of degree \eqn{d} in either shifted or non-shifted form can be exactly fitted. The resulting piecewise polynomial coefficients \eqn{b_{0i},\ b_{1i}, \ldots,\ b_{di}} or \eqn{a_{0i},\ a_{1i}, \ldots,\ a_{di}} are stored in the \eqn{i}-th column of a \eqn{(d + 1) \times (k - 1)} coefficient matrix.

The piecewise polynomial representation makes evaluation of splines and its derivatives very easy. It is also straightforward to solve equations like \eqn{f(x) = y} and \eqn{f'(x) = y} for \eqn{x} values given a \eqn{y} value. In particular, for \eqn{y = 0} the solutions to two equations respectively give the roots and extrema of the spline.
}

\value{
An object of ``PiecePoly'' class is a list with at least the following two components:

  \item{PiecePoly}{A list with two components: \code{coef} and \code{shift}. The former is a \eqn{(d + 1) \times (k - 1)} matrix, and the latter is a logical value.}
  \item{knots}{A vector of knots of the spline.}

\code{\link{CubicInterpSplineAsPiecePoly}} and \code{\link{SmoothSplineAsPiecePoly}} add an extra ``CubicInterpSpline'' class, while \code{\link{RegBsplineAsPiecePoly}} adds and extra ``RegBspline'' class.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

