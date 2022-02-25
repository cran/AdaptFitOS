\name{aspFormula}
\alias{aspFormula}
\title{
An asp formula
}
\description{
A formula to be used in \code{\link{asp2}}. 

Dummies for categorical covariates are constructed automatically if a variable is given as factor (with contrasts as set by \code{options("contasts")} or specified by a list in argument \code{contrasts}). Note that only parametric interactions are supported and that
interacting covariates have to be multiplied beforehand and given as a new variable in the formula.
Smooth terms are given by

f(x,basis="os",degree=3,knots,var.knots,var.basis,var.degree,adap=TRUE)

with the following arguments:
}

\arguments{
  \item{x}{the covariate}
  \item{basis}{the spline basis function to be used. "trunc.poly" for truncated polynomials, "tps" for thin plate splines and "os" for B-splines (default).}
  \item{degree}{the degree of the basis. In the case of B-splines also a vector of the form c(p,q)
  with p the B-spline degree and q the penalty order (the integrated q-th squared derivative is penalized, see references).
  If only a scalar is given q is chosen such that p=2*q-1. Defaults are degree=3 (basis="tps"), degree=1 (basis="trunc.poly") and degree=c(3,2) (basis="os"), respectively.}
  \item{knots}{the knots to be used. Using e.g. kn=default.knots(x,40) beforehand leads to 40 quantile based knots in the case of "tps" and "trunc.poly" bases.
  In the case of B-splines ("os"), knots are always equidistant and are automatically generated with the number equal to the length of the vector of knots given plus boundary knots.
  If no knots are given the number of knots is automatically chosen to be equal to floor(n/max(4, floor(n/35)) - 1).}
  \item{adap}{TRUE for spatially adaptive smoothing parameter}
  \item{var.knots}{the knots for the spline basis for adaptively estimating the smoothing parameter. If missing the number of knots is automatically chosen to be equal to floor(knots/max(4, floor(knots/35)) - 1).}
 \item{var.basis}{spline basis function for adaptive smoothing parameter estimation. If missing, the same basis as for estimation of f is used.}
 \item{var.degree}{spline degree for adaptive smoothing parameter estimation. If missing, the same degree as for estimation of f is used.}
  \item{spar}{the smoothing parameter if desired. Usually this is left unspecified, such that the smoothing parameter is estimated by restricted maximum likelihood (see references). Currently doesn't work for basis="os".}
  \item{df}{the number of degrees of freedom corresponding to the REML choice of smoothing parameter if desired. Usually this is left unspecified, such that the smoothing parameter is estimated by restricted maximum likelihood (see references). Currently doesn't work for basis="os".}
}
