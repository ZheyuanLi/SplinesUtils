# SplinesUtils

## Overview

This is version 0.2 of R package SplinesUtils. The previous version is 0.1-1, committed on 2018-10-08. The new version fixes a few code and documentation bugs, also with enhanced functionality (see story below). This package is now making its way to CRAN.

## Alert

Function `RegBsplineAsPiecePoly` is now obsolete, and superseded by `RegSplineAsPiecePoly`. The older function can still be used for now, but will pop up a warning.

## Description

Splines that are used in regression or smoothing models often adopt basis representation. While such representation is recommended for numerical computations, it is neither intuitive to interpret, nor easy to derive mathematical property of the spline. Package "SplinesUtils" provides functions that reparametrize univariate cubic interpolation splines, smoothing splines and regression splines (used in "lm", "glm", "lme") to piecewise polynomials, by exporting piecewise polynomial coefficients as a matrix and print piecewise polynomial equations as formatted strings. The package also provides generic functions like 1) `plot` and `predict` for plotting and predicting a spline or its derivatives; 2) `solve` for solving a spline or its derivatives for a given value, specially useful for identifying local extrama of the spline.

## Installation

First install R package `devtools`, then install `SplinesUtils` with R command:

```
devtools::install_github("ZheyuanLi/SplinesUtils")
```

## Story

The package (version 0.1-1) was specially motivated by [this Stack Overflow thread](https://stackoverflow.com/q/44739192/4891738). The package was then used by some researchers and I was occasionally approached for [help on Stack Overflow](https://stackoverflow.com/q/57609337/4891738) or [small bug report](https://github.com/ZheyuanLi/SplinesUtils/issues/1).

Recently another researcher contacted me, asking whether regression splines in a linear mixed model can be exported as well. This is definitely possible but extra investigation is needed. The researcher then opened [a thread on Stack Overflow](https://stackoverflow.com/q/65303707/4891738) offering a good, reproducible example, which leads to version 0.2 of the package.

