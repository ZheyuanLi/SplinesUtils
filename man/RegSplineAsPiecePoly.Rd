\name{RegSplineAsPiecePoly}

\alias{RegSplineAsPiecePoly}

\title{
Reparametrize a fitted bs() or ns() term in a linear model ("lm"), generalized linear model ("glm"), or linear mixed model ("lme") as piecewise polynomials.
}

\description{
Reparametrize a single fitted bs() or ns() term in a linear model ("lm"), generalized linear model ("glm"), or linear mixed model ("lme") as piecewise polynomials. Do not interpret the spline as fitted values or linear predictors of the model.
}

\usage{
RegSplineAsPiecePoly(RegModel, SplineTerm, shift = TRUE)
}

\arguments{
  \item{RegModel}{A fitted "lm", "glm" or "lme" that has at least one spline term constructed by \code{splines::bs} or \code{splines:ns}.}
  \item{SplineTerm}{A character string giving the name of the spline term to be exported.}
  \item{shift}{A logical value determining whether the piecewise polynomial should be parametrized in shifted form.}
}

\value{
An object of ``PiecePoly'' and ``RegSpline'' classes. See \link{PiecePolyObject}.
}

\author{
Zheyuan Li \email{zheyuan.li@bath.edu}
}

\examples{
\dontrun{
## a toy data set
set.seed(0)
## a curve of 'x' (note that my 'x' are not uniformly sampled)
x <- sort(rnorm(400, 0, pi))
fx <- sin(x) + 0.2 * x + cos(abs(x))
x_range <- range(x)
## a 2-level grouping variable (1st group with mean 0 and 2nd group with mean 1)
g <- sample(gl(2, 200, labels = letters[1:2]))
fg <- c(0, 1)[g]
## a linear effect
u <- runif(400, x_range[1], x_range[2])
fu <- -0.05 * u
## a spurious covariate
z <- runif(400, x_range[1], x_range[2])
## an additive model
y <- fx + fg + fu + rnorm(400, 0, 0.5)

## splines are constructed so that boundary knots are clamped
## df = length(knots) + degree + has.intercept
## for 'ns()', degree = 3L is fixed

## fit a linear model
fit <- lm(y ~ g + bs(x, df = 20) + u + ns(z, df = 10))

## export spline term "bs(x, df = 20)"
spl_bs <- RegSplineAsPiecePoly(fit, "bs(x, df = 20)")

## predict the spline at observed points
y_bs <- predict(spl_bs, x)
## not that this is not the fitted values
all.equal(y_bs, fitted(fit), check.attributes = FALSE)
#[1] "Mean relative difference: 0.8498324"
## it is identical to the following, up to a constant vertical shift
## (`predict.lm` applies term-wise centering, see https://stackoverflow.com/q/37963904/4891738)
y_bs2 <- predict(fit, type = "terms", terms = "bs(x, df = 20)")
## the standard error of the difference is 0, i.e., the difference is constant
c(sd(y_bs - y_bs2))
#[1] 9.727983e-15

## plot the spline
plot(spl_bs)

## find all its local extrema
xs <- solve(spl_bs, b = 0, deriv = 1)
ys <- predict(spl_bs, xs)
points(xs, ys, pch = 19)
}
}

