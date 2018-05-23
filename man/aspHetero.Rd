\name{aspHetero}
\alias{aspHetero}
\title{ Estimate varying residual variance  }
\description{Estimates a varying residual variance on basis of an asp object. Resulting object can be plotted with simultaneous confidence bands corrected for heteroscedasticity.}
\usage{
aspHetero(object, xx, nknots=5, knots=NULL, basis="os", 
		 degree=c(3,2), tol=1e-8, niter=100, niter.var=250)
}
\arguments{
\item{object}{an asp object.}
\item{xx}{the covariate.}
\item{nknots}{the number of knots. Does not apply when \code{knots} are given.}
\item{knots}{the knots. Does not apply if basis=="os". Otherwise, if \code{NULL} \code{nknots} equidistant knots are used.}
\item{basis}{the spline basis: "os" (default), "trunc.poly" or "tps".}
\item{degree}{the spline degree (and penalty order in case of B-splines). Defaults to c(3,2).}
\item{tol}{tolerance for the convergence criterion. Default is 1e-8.}
\item{niter}{a maximum number of iterations for residual variance function estimation, default is 100.}
\item{niter.var}{a maximum number of iterations for the variance of
    random effects estimation within the residual variance function estimation routine, default is 250.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
An object of class asp with varying variances, with additional element \code{sigmax} including information on the spline of the varying variance.}

\references{
Wiesenfarth, M., Krivobokova, T., Klasen, S., Sperlich, S. (2012).\cr
Direct Simultaneous Inference in Additive Models and its Application to Model Undernutrition.
\emph{Journal of the American Statistical Association}, 107(500): 1286-1296.
}
\author{
Manuel Wiesenfarth \email{m.wiesenfarth at dkfz de}
}

\examples{
attach(mcycle)

y=accel
kn1 <- default.knots(times,20)
# fit model with constant residual variance
  fit= asp2(accel~f(times,basis="os",degree=3,knots=kn1,adap=FALSE),
  								niter = 20, niter.var = 200)


# fit model with varying residual variance
fith=aspHetero(fit,times,tol=1e-8)
op <- par(mfrow = c(1,3))
plot(fit);plot(fith)
#sigma() returns the fitted varying residual variance
plot(sort(times),sigma(fith)[order(times)],type="l")
par(op)
}
\keyword{smooth}
