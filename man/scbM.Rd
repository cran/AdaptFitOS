\name{scbM}
\alias{scbM}
\title{Calculate simultaneous confidence bands for penalized splines}
\description{
Calculates simultaneous (uniform) confidence bands for the mixed model representation of penalized
splines based on volume-of-tube formula.
}
\usage{
scbM(object,select=NULL,drv=0,level=0.95,div=1000,calc.stdev=TRUE,bayes=FALSE)
}
\arguments{
  \item{object}{an asp object.}
  \item{select}{vector specifying which curves in an additive model should be considered. If NULL, all curves are considered.}
  \item{drv}{the derivative order. Defaults to 0, i.e. the estimated function itself is plotted. First and second derivatives are supported.}
  \item{level}{level of confidence.}
\item{div}{precision for the integral used for calculation of the length of the curve, default is 1000.}
%\item{mixed}{FALSE for the confidence bands with the (approximate) frequentist coverage probability, TRUE  for the (approximate) Bayesian coverage probability}
\item{calc.stdev}{TRUE to compute standard deviation and confidence bands for each value of the covariates. Computationally intensive for large data sets. Use plot.scbm() or plot.asp() to compute standard deviation and bounds only for a grid.
If FALSE only critical values are computed.}
\item{bayes}{FALSE for confidence bands with (approximate) frequentist coverage probability, TRUE for (approximate) Bayesian coverage probability. See Krivobokova et al. (2010) for details.}
}
\details{Returns a \code{scbm} object and prints critical values.
The resulting confidence bands have (approximate) frequentist coverage probabilities with automatic bias correction (see references).
Makes use of the \code{libtube} library by Catherine Loader (see package \code{\link[locfit]{locfit}}).
}

\value{
A list object of class \code{scbm} containing
\item{aspobject}{an asp object.}
\item{drv}{the derivative order.}
\item{crit}{a list of critical values.}
\item{sigma2}{the variance of the residuals.}
\item{cov.coef}{a list of covariance matrices of spline coefficients in the mixed model framework.}
\item{Stdev}{the standard deviations of estimates. Only given if calc.stdev=TRUE.}
\item{fitted}{a list of fitted values. Only given if calc.stdev=TRUE.}
\item{lcb}{a list of lower bounds of confidence bands. Only given if calc.stdev=TRUE.}
\item{ucb}{a list of upper bounds of confidence bands. Only given if calc.stdev=TRUE.}
\item{\dots}{ further}
}
\references{
Krivobokova, T., Kneib, T., and Claeskens, G. (2010)\cr
Simultaneous confidence bands for penalized spline estimators.
\emph{Journal of the American Statistical Association}, 105(490):852-863.

Wiesenfarth, M., Krivobokova, T., Klasen, S., Sperlich, S. (2012).\cr
Direct Simultaneous Inference in Additive Models and its Application to Model Undernutrition.
\emph{Journal of the American Statistical Association}, accepted.
}
\author{
Manuel Wiesenfarth \email{mwiesen at gwdg de},
Tatyana Krivobokova \email{tkrivob at gwdg de}
}

\examples{
\dontrun{
beta=function(l,m,x) return(gamma(l+m)*(gamma(l)*gamma(m))^(-1)*x^(l-1)*(1-x)^(m-1))
f1 = function(x) return((0.6*beta(30,17,x)+0.4*beta(3,11,x))*1/0.958)
f2 = function(z) return((sin(2*pi*(z-0.5))^2)*1/.3535)
f3 = function(z) return((exp(-400*(z-0.6)^2)+5/3*
                           exp(-500*(z-0.75)^2)+2*exp(-500*(z-0.9)^2))*1/0.549)
center=function(x) return(x-mean(x))

set.seed(1)
N <- 500
x1 = runif(N,0,1)
x2 = runif(N,0,1)
x3 = runif(N,0,1)


kn1 <- default.knots(x1,40)
kn2 <- default.knots(x2,40)
kn3 <- default.knots(x3,40)
kn.var3 <- default.knots(kn3,5)

y <- f1(x1)+f2(x2)+f3(x3)+0.3*rnorm(N)

# fit model with last smoothing parameter adaptive
  fit2b=  asp2(y~f(x1, basis="os", degree=3, knots=kn1, adap=FALSE)
                +f(x2, basis="os", degree=3, knots=kn2, adap=FALSE)
                +f(x3, basis="os", degree=3, knots=kn3, adap=TRUE,
                       var.knots=kn.var3, var.basis="os", var.degree=3),
                niter = 20, niter.var = 200)

  # compute 95% simultaneous confidence bands.
  # You could skip this and use "fit2b" indstead of "scb2b" later on, however,
  # if N is large, computing the SCBs various times can take some time
  # if you don't need fitted values and bounds for all covariate points
  # (can be computationally intensive due to large matrix dimensions),
  # set calc.stdev=F such that these are not computed.
    scb2b<- scbM(fit2b,calc.stdev=FALSE)
    plot(scb2b,pages=1)

  # plot first derivative of f(x1)
    scb2bdrv<- scbM(fit2b,drv=1,calc.stdev=FALSE)
    plot(scb2bdrv,select=1)
    #the following would give the same result
    #plot(fit2b,select=1,drv=1)
    #different style
    plot(scb2bdrv,select=1,scb.lty="blank", shade=TRUE,shade.col="steelblue")
}
}
\keyword{smooth}
\keyword{dplot}

