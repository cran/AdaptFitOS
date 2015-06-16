\name{default.knots}
\alias{default.knots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute default knots for a given x vector
}
\description{
Computes default knots for a given x vector.}
\usage{
default.knots(x, num.knots, knotchoice="quantiles")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  The covariate. Note that for B-splines, only the range of \code{x} is considered.
}
  \item{num.knots}{
  The number of knots. Defaults to \code{floor(n/max(4, floor(n/35)) - 1)}. 
}
  \item{knotchoice}{
  Either "equidistant" or "quantiles" for equidistant and quantile based knots, respectively. Note that in case of B-splines, knots are always equidistant.
}
}
