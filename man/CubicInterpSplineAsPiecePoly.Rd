\name{CubicInterpSplineAsPiecePoly}

\alias{CubicInterpSplineAsPiecePoly}

\title{
Construct a cubic interpolation spline as piecewise cubic polynomials.
}

\description{
The function is a thin wrapper of \code{stats::splinefun}. The latter readily constructs the interpolation spline as piecewise cubic polynomials but stores construction data like knots and piecewise polynomial coefficients in an environment. The Function reorganizes such information into a "PiecePoly" object.
}

\usage{
CubicInterpSplineAsPiecePoly(x, y, method)
}

\arguments{
  \item{x}{A vector of \eqn{x}-coordinate values of the data points to interpolate.}
  \item{y}{A vector of \eqn{y}-coordinate values of the data points to interpolate.}
  \item{method}{Methods for cubic spline interpolation, including ``fmm'', ``natural'', ``periodic'' and ``hyman'' in \code{stats::spline} and \code{stats::splinefun}.}
}

\value{
An object of class ``PiecePoly'' and ``CubicInterpSpline'' classes. See \link{PiecePolyObject}. Note that the piecewise polynomial is always parametrized in shifted form.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
\dontrun{
## a toy dataset
set.seed(0)
## data for "fmm" and "natural" method
x <- 1:10 + runif(10, -0.1, 0.1)
y1 <- rnorm(10, 3, 1)
## periodic data for "periodic" method
x2 <- c(x, x[10] + runif(1))
y2 <- c(y1, y1[1])
## monotone data for "hyman" method
y3 <- sort(y1)

## perform interpolation
fmm <- CubicInterpSplineAsPiecePoly(x, y1, "fmm")
natural <- CubicInterpSplineAsPiecePoly(x, y1, "natural")
periodic <- CubicInterpSplineAsPiecePoly(x2, y2, "periodic")
hyman <- CubicInterpSplineAsPiecePoly(x, y3, "hyman")

## test print method
fmm
natural
periodic
hyman

## export all polynomial equations as formatted strings
fmm_eqn <- summary(fmm)
natural_eqn <- summary(natural)
periodic_eqn <- summary(periodic)
hyman_eqn <- summary(hyman)

## plot all splines
par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
plot(fmm, show.knots = TRUE); points(x, y1)
plot(natural, show.knots = TRUE); points(x, y1)
plot(periodic, show.knots = TRUE); points(x2, y2)
plot(hyman, show.knots = TRUE); points(x, y3)

## plot 1st derivatives of all splines
par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
plot(fmm, deriv = 1)
plot(natural, deriv = 1)
plot(periodic, deriv = 1)
plot(hyman, deriv = 1)

## backsolve the spline given y = 2.85
xr1 <- solve(fmm, b = 2.85)
xr2 <- solve(natural, b = 2.85)
xr3 <- solve(periodic, b = 2.85)
xr4 <- solve(hyman, b = 2.85)
par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
plot(fmm); abline(h = 2.85, lty = 2); points(xr1, rep.int(2.85, length(xr1)))
plot(natural); abline(h = 2.85, lty = 2); points(xr2, rep.int(2.85, length(xr2)))
plot(periodic); abline(h = 2.85, lty = 2); points(xr3, rep.int(2.85, length(xr3)))
plot(hyman); abline(h = 2.85, lty = 2); points(xr4, rep.int(2.85, length(xr4)))

## find all extrema
## "hyman" excluded as it is monotonic spline
xs1 <- solve(fmm, deriv = 1)
xs2 <- solve(natural, deriv = 1)
xs3 <- solve(periodic, deriv = 1)
par(mfrow = c(2, 2), mar = c(4, 4, 1, 1))
ys1 <- predict(fmm, xs1)
ys2 <- predict(natural, xs2)
ys3 <- predict(periodic, xs3)
plot(fmm); points(xs1, ys1, pch = 19)
plot(natural); points(xs2, ys2, pch = 19)
plot(periodic); points(xs3, ys3, pch = 19)
}
}
