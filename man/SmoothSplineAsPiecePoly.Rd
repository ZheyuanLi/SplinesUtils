\name{SmoothSplineAsPiecePoly}

\alias{SmoothSplineAsPiecePoly}

\title{
Export a fitted smoothing spline as piecewise cubic polynomials.
}

\description{
The function extracts knots from a fitted smoothing spline object, evaluate the smoothing spline at knots, then calls \code{\link{CubicInterpSplineAsPiecePoly}} to construct a natural cubic interpolation spline.
}

\usage{
SmoothSplineAsPiecePoly(SmoothSpline)
}

\arguments{
  \item{SmoothSpline}{A fitted smoothing spline as returned by \code{stats::smooth.spline}.}
}

\value{
An object of ``PiecePoly'' and ``CubicInterpSpline'' class. Note that the piecewise polynomial is always parametrized in shifted form.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{

\dontrun{
## a toy dataset
set.seed(0)
x <- 1:100 + runif(100, -0.1, 0.1)
y <- poly(x, 9) %*% rnorm(9)
y <- y + rnorm(length(y), 0, 0.2 * sd(y))

## fit a smoothing spline
sm <- smooth.spline(x, y)

## coerce "smooth.spline" object to "PiecePoly" object
oo <- SmoothSplineAsPiecePoly(sm)

## print the "PiecePoly"
oo

## plot the "PiecePoly"
plot(oo)

## find all the roots
xr <- solve(oo)
points(xr, rep.int(0, length(xr)))
abline(h = 0, lty = 2)

## find all stationary / saddle points
xs <- solve(oo, deriv = 1)

## predict the "PiecePoly" at stationary / saddle points
ys <- predict(oo, xs)
points(xs, ys, pch = 19)
}
}
