mystop <- function (msg) stop("\n  ", msg, call. = FALSE)

########################################################
# cubic interpolation splines in R core as "PiecePoly" #
########################################################

CubicInterpSplineAsPiecePoly <- function (x, y, method = c("fmm", "natural", "periodic", "hyman")) {
  ## method validation
  if (!(method %in% c("fmm", "natural", "periodic", "hyman")))
    mystop("'method' must be one of the following: 'fmm', 'natural', 'periodic', 'hyman'!")
  ## use `splinefun` for cubic spline interpolation
  CubicInterpSpline <- stats::splinefun(x, y, method)
  ## extract construction info
  construction_info <- environment(CubicInterpSpline)$z
  ## export as an "PiecePoly" object
  pieces <- seq_len(length(construction_info$x) - 1L)
  PiecePolyCoef <- with(construction_info, rbind(y[pieces], b[pieces], c[pieces], d[pieces], deparse.level = 0L))
  structure(list(PiecePoly = list(coef = PiecePolyCoef, shift = TRUE),
                 knots = construction_info$x), method = method,
                 class = c("PiecePoly", "CubicInterpSpline"))
  }

###################################################
# cubic smoothing spline in R core as "PiecePoly" #
###################################################

## represent a fitted smoothing spline as an interpolation spline
SmoothSplineAsPiecePoly <- function (SmoothSpline) {
  ## input validation
  if (!inherits(SmoothSpline, "smooth.spline"))
    mystop("This function only works with models that inherit 'smooth.spline' class!")
  ## knots of the smoothing spline
  kx <- with(SmoothSpline$fit, knot * range + min)
  kx <- kx[4:(length(kx) - 3)]
  ky <- predict(SmoothSpline, kx, 0L)[[2]]  ## deriv = 0L
  ## natural cubic spline interpolation over the knots
  CubicInterpSplineAsPiecePoly(kx, ky, method = "natural")
  }

##############################################################
# `bs` or `ns` term in a "lm", "glm" or "lme" as "PiecePoly" #
##############################################################

ExtractSplineTerm <- function (termsObject, SplineTerm) {
  ## term labels
  tl <- attr(termsObject, "term.labels")  ## tl <- labels(termsObject)
  ## is SplineTerm found in termsObject?
  if (!(SplineTerm %in% tl)) {
    ## header
    h <- "\n  Required SplineTerm not found! Available terms are:\n"
    ## ordered list
    ol <- sprintf("    * %s\n", tl)
    ## stop
    stop(h, ol, call. = FALSE)
  }
  ## is SplineTerm a bs() or ns()?
  is_bs <- grepl("bs(", SplineTerm, fixed = TRUE)
  is_ns <- grepl("ns(", SplineTerm, fixed = TRUE)
  if (!(is_bs || is_ns))
    mystop("Provided SplineTerm is neither 'bs()' nor 'ns()' from package 'splines'!")
  ## position of SplineTerm in termsObject
  pos <- match(SplineTerm, tl)
  ## extract predict call
  predvars <- attr(termsObject, "predvars")  ## try: as.list(predvars)
  SplineCall <- predvars[[2L + pos]]  ## try: as.list(SplineCall)
  ## change variable name in SplineCall to x
  SplineCall[[2]] <- quote(x)
  ## extract degree and knots from SplineCall
  if (is_bs) {
    degree <- SplineCall[[3]]
    interior_knots <- unname(SplineCall[[4]])
    boundary_knots <- SplineCall[[5]]
    SplineCall[[4]] <- interior_knots
  } else {
    degree <- 3L
    interior_knots <- unname(SplineCall[[3]])
    boundary_knots <- SplineCall[[4]]
    SplineCall[[3]] <- interior_knots
  }
  x <- c(boundary_knots[1], interior_knots, boundary_knots[2])
  ## return
  list(pos = pos, degree = degree, knots = x, call = SplineCall)
}

PiecePolyRepara <- function (SplineTerm, SplineCoef, shift = TRUE) {
  x <- SplineTerm$knots
  degree <- SplineTerm$degree
  SplineCall <- SplineTerm$call
  n_pieces <- length(x) - 1L
  ## initialize piecewise polynomial coefficients (matrix)
  PiecePolyCoef <- matrix(0, degree + 1L, n_pieces)
  ## loop through pieces
  i <- 1L
  while (i <= n_pieces) {
    ## an evenly-spaced grid between x[i] and x[i + 1]
    xg <- seq.int(x[i], x[i + 1L], length.out = degree + 1L)
    ## spline values on the grid
    yg <- c(eval(SplineCall, data.frame(x = xg)) %*% SplineCoef)
    ## basis matrix on the grid
    Xg <- outer(xg - shift * x[i], 0:degree, "^")
    ## estimate linear regression yg ~ Xg + 0 by solving normal equation
    A <- base::crossprod(Xg)
    b <- base::crossprod(Xg, yg)
    U <- chol.default(A)
    PiecePolyCoef[, i] <- base::backsolve(U, base::forwardsolve(t.default(U), b))
    ## next piece
    i <- i + 1L
  }
  ## return
  PiecePolyCoef
}

## export a `bs` or `ns` term in a "lm", "glm" or "lme" as piecewise polynomials
RegSplineAsPiecePoly <- function (RegModel, SplineTerm, shift = TRUE) {
  
  ## input validation
  if (!inherits(RegModel, c("lm", "glm", "lme")))
    mystop("This function only works with models that inherit 'lm', 'glm' or 'lme' class!")

  ## extract "terms" on the RHS of the model formula
  RegSpline <- ExtractSplineTerm(terms(RegModel), SplineTerm)
  pos <- RegSpline$pos
  
  ## regression coefficients associated with SplineTerm
  if (inherits(RegModel, c("lm", "glm"))) {
    RegSplineCoef <- with(RegModel, coefficients[assign == pos])
  } else {
    ind <- attr(RegModel$fixDF, "assign")[[SplineTerm]]
    RegSplineCoef <- RegModel$coefficients$fixed[ind]  ## fixed.effects(RegModel)
  }
  RegSplineCoef <- unname(RegSplineCoef)

  ## is there `NA` in coefficients?
  na <- is.na(RegSplineCoef)
  if (any(na)) {
    warning("NA coefficients found for SplineTerm; Replacing NA by 0")
    RegSplineCoef[na] <- 0
    }

  ## it is not possible to use CubicInterpSplineAsPiecePoly for reparametrization
  ## it only supports ns() term but not bs() term
  PiecePolyCoef <- PiecePolyRepara(RegSpline, RegSplineCoef, shift)

  ## return coefficient matrix, shift, knots, SplineCall and RegSplineCoef
  structure(list(PiecePoly = list(coef = PiecePolyCoef, shift = shift),
                 knots = RegSpline$knots), class = c("PiecePoly", "RegSpline"))
}

RegBsplineAsPiecePoly <- function (lm_glm, SplineTerm, shift = TRUE) {
  warning("RegBsplineAsPiecePoly is now obsolete and likely to be removed in the next update; use RegSplineAsPiecePoly")
  RegSplineAsPiecePoly(lm_glm, SplineTerm, shift)
}

##############################
# S3 methods for "PiecePoly" #
##############################

## print a polynomial equation given polynomial coefficients
## the workhorse function for `summary.PiecePoly`
PolyEqn <- function (pc, shift, knot) {

  ## split coefficients into three parts: pc0, pc1, pc2
  pc0 <- pc[1L]       ## coefficient for 1
  pc1 <- pc[2L]       ## coefficient for x
  pc2 <- pc[-c(1:2)]  ## coefficients for higher power term

  ###############################################
  # print polynomial terms as formatted strings #
  ###############################################

  ## reverse sign of knot
  rev_knot_sgn <- if (knot > 0) " - " else " + "
  knot <- abs(knot)

  ## 0-th order term
  xterm0 <- sprintf("%.3g", pc0)
  ## 1-st order term: sign, absolute coefficient, base
  coef_sgn <- if (pc1 > 0) " + " else " - "
  if (shift) {
    xterm1 <- sprintf("%s%.3g * (x%s%.3g)", coef_sgn, abs(pc1), rev_knot_sgn, knot)
    } else {
    xterm1 <- sprintf("%s%.3g * x", coef_sgn, abs(pc1))
    }
  ## higher order term: sign, absolute coefficient, base, power
  if (length(pc2)) {
    sgn <- rep.int(" - ", length(pc2))
    sgn[pc2 > 0] <- " + "
    if (shift) {
      xterm2 <- sprintf("%s%.3g * (x%s%.3g) ^ %d", sgn, abs(pc2), rev_knot_sgn, knot, seq_along(pc2) + 1L)
      } else {
      xterm2 <- sprintf("%s%.3g * x ^ %d", sgn, abs(pc2), seq_along(pc2) + 1L)
      }
    xterm2 <- paste0(xterm2, collapse = "")
    } else xterm2 <- ""

  ## concatenate all terms to a single polynomial equation
  paste0(xterm0, xterm1, xterm2, collapse = "")
  }

## `summary` method for "PiecePoly"
##   summary
##   function (object, ...)
## export piecewise polynomial equations as formatted strings
summary.PiecePoly <- function (object, no.eqn = NULL, ...) {
  ## change symbol
  PiecePolyObject <- object
  ## extract piecewise polynomial info
  PiecePolyCoef <- PiecePolyObject[[c(1, 1)]]
  shift <- PiecePolyObject[[c(1, 2)]]
  knots <- PiecePolyObject$knots
  n_pieces <- dim(PiecePolyCoef)[[2]]
  cat(sprintf("%d piecewise polynomials of degree %d are constructed!\n", n_pieces, dim(PiecePolyCoef)[1L] - 1L))
  ## how many pieces to summary?
  if (is.null(no.eqn)) no.eqn <- n_pieces
  else no.eqn <- min(no.eqn, n_pieces)
  ## call `PolyEqn`
  EqnLst <- vector("list", no.eqn)
  i <- 1L
  while (i <= no.eqn) {
    EqnLst[[i]] <- PolyEqn(PiecePolyCoef[, i], shift, knots[i])
    i <- i + 1L
    }
  ## invisibly return equation list
  output <- EqnLst
  }

## `print` method for "PiecePoly"
##    print
##    function (x, ...) 
## it summarizes and prints at most 6 piecewise polynomial equations
print.PiecePoly <- function (x, ...) {
  ## change symbol
  PiecePolyObject <- x
  ## summarize at most 6 equations
  Head <- summary.PiecePoly(PiecePolyObject, 6L)
  ## helpful message
  cat(sprintf("Use 'summary' to export all of them.\nThe first %d are printed below.\n", length(Head)))
  for (i in seq_along(Head)) {
    cat(Head[[i]]); cat("\n")
    }
  }

## `plot` method for "PiecePoly"
##    plot
##    function (x, y, ...)
plot.PiecePoly <- function (x, spread = 3L, deriv = 0L, show.knots = FALSE, ...) {
  ## change symbol
  PiecePolyObject <- x
  ## extract piecewise polynomial coefficients
  PiecePolyCoef <- PiecePolyObject$PiecePoly$coef
  shift <- PiecePolyObject$PiecePoly$shift
  n_pieces <- dim(PiecePolyCoef)[2L]
  ## get degree and power
  degree <- dim(PiecePolyCoef)[1L] - 1L
  ## extract knots
  x <- PiecePolyObject$knots
  ## deriv validation
  if (deriv > degree) {
    plot.default(x, rep.int(0, length(x)), "l", ylab = "y")
    return(NULL)
    }
  ## get power
  power <- 0:(degree - deriv)
  ## evaluation grid
  xg <- vector("list", n_pieces)  ## x-values of spline
  yg <- vector("list", n_pieces)  ## y-values of spline
  ## loop through pieces
  i <- 1L
  while (i <= n_pieces) {
    ## an evenly-spaced grid between x[i] and x[i + 1]
    xg[[i]] <- seq.int(x[i], x[i + 1], length.out = spread * (degree + 1L))
    ## evaluate the polynomial with computed coefficients
    pc <- PiecePolyCoef[, i]
    if (deriv > 0) pc <- pc[-seq_len(deriv)] * choose(deriv:degree, deriv) * factorial(deriv)
    yg[[i]] <- c(outer(xg[[i]] - shift * x[i], power, "^") %*% pc)
    ## next piece
    i <- i + 1L
    }
  xg <- unlist(xg)
  yg <- unlist(yg)
  ## plot the spline and mark knots locations
  plot.default(xg, yg, "l", xlab = "x", ylab = "y")
  if (show.knots) abline(v = x, lty = 3, col = 8)
  ## invisibly return data for plotting
  xgyg <- list(x = xg, y = yg)
  ## interior knots are included twice in 'xg', but this is actually good
  ## in case 'deriv' = 'degree', we don't need 'type = "s"' inside 'plot.default'
  }

## `predict` method for "PiecePoly"
##    predict
##    function (object, ...)
predict.PiecePoly <- function (object, newx, deriv = 0L, ...) {
  ## change symbol
  PiecePolyObject <- object
  ## extract piecewise polynomial coefficients
  PiecePolyCoef <- PiecePolyObject$PiecePoly$coef
  shift <- PiecePolyObject$PiecePoly$shift
  ## get degree
  degree <- dim(PiecePolyCoef)[1L] - 1L
  ## deriv validation
  if (deriv > degree) return(numeric(length(newx)))
  ## get power
  power <- 0:(degree - deriv)
  ## extract knots
  x <- PiecePolyObject$knots
  ## which piece?
  piece_id <- findInterval(newx, x, TRUE)
  ind <- split.default(seq_len(length(newx)), piece_id)
  unique_piece_id <- as.integer(names(ind))
  n_pieces <- length(unique_piece_id)
  ## loop through pieces
  y <- numeric(length(newx))
  i <- 1L
  while (i <= n_pieces) {
    ii <- unique_piece_id[i]
    xi <- newx[ind[[i]]] - shift * x[ii]
    pc <- PiecePolyCoef[, ii]
    if (deriv > 0) pc <- pc[-seq_len(deriv)] * choose(deriv:degree, deriv) * factorial(deriv)
    y[ind[[i]]] <- c(outer(xi, power, "^") %*% pc)
    i <- i + 1L
    }
  y
  }

## `solve` method for "PiecePoly"
##    solve
##    function (a, b, ...)
##  1. backsolve 'x' value given a 'y' value on the spline
##  2. find extrema of the spline
solve.PiecePoly <- function (a, b = 0, deriv = 0L, ...) {
  ## change symbol
  PiecePolyObject <- a
  y <- b
  ## helpful message (in case someone used `y = y0` than `b = y0` to give RHS which returns misleading results)
  cat(sprintf("solving equation for RHS value %.7g\n", y))
  ## extract piecewise polynomial coefficients
  PiecePolyCoef <- PiecePolyObject$PiecePoly$coef
  shift <- PiecePolyObject$PiecePoly$shift
  n_pieces <- dim(PiecePolyCoef)[2L]
  ## get degree
  degree <- dim(PiecePolyCoef)[1L] - 1L
  ## extract knots
  x <- PiecePolyObject$knots
  ## deriv validation
  if (deriv >= degree) mystop("'deriv' can not exceed 'degree'!")
  ## list of roots on each piece
  xr <- vector("list", n_pieces)
  ## loop through pieces
  i <- 1L
  while (i <= n_pieces) {
    ## polynomial coefficient
    pc <- PiecePolyCoef[, i]
    ## take derivative
    if (deriv > 0) pc <- pc[-seq_len(deriv)] * choose(deriv:degree, deriv) * factorial(deriv)
    pc[1] <- pc[1] - y
    ## complex roots
    croots <- base::polyroot(pc)
    ## real roots (be careful when testing 0 for floating point numbers)
    rroots <- Re(croots)[round(Im(croots), 10) == 0]
    ## is shifting needed?
    if (shift) rroots <- rroots + x[i]
    ## real roots in (x[i], x[i + 1])
    xr[[i]] <- rroots[(rroots >= x[i]) & (rroots <= x[i + 1])]
    ## next piece
    i <- i + 1L
    }
  ## collapse list to atomic vector and return
  unlist(xr)
  }
