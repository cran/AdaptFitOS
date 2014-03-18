\name{asp2}
\alias{asp2}
\title{Fit a semiparametric regression model with spatially adaptive penalized splines}
\description{
Fits semiparametric
regression models using the mixed model
representation of penalized splines with spatially adaptive
penalties based on the \code{\link[AdaptFit]{asp}} function of the AdaptFit-package (which is based on \code{\link{spm}} function of \code{SemiPar} package). Particular differences to AdaptFit
include the availability of simultaneous confidence bands and B-spline basis functions and different functionality of the \code{plot} function. However, random effects, autocorrelations and interaction surfaces are not supported.
Further, only Gaussian responses are supported. Note that in contrast to AdaptFit and SemiPar, estimated curves are centered to have zero mean.
}
\usage{
asp2(form, spar.method = "REML", contrasts=NULL, omit.missing = NULL, returnFit=FALSE,
     niter = 20, niter.var = 50, tol=1e-6, tol.theta=1e-6, control=NULL)
}
\arguments{
  \item{form}{a formula describing the model to be fitted. See \code{\link{aspFormula}} for further information.
          Note, that an intercept is always included, whether 
          given in the formula or not.}
% \item{adap}{\code{TRUE} if an adaptive fit should be performed, otherwise the fit is identical to that of function \code{spm} (only that curves are centered to have zero mean).}
 \item{spar.method}{method for automatic smoothing parameter
                     selection. May be "REML" (restricted maximum
                     likelihood) or "ML" (maximum likelihood).}
 \item{contrasts}{an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.}
  \item{omit.missing}{a logical value indicating whether
                fields with missing values are to be
                      omitted.}
  \item{niter}{a maximum number of iterations for the mean estimation, default is 20.}
  \item{niter.var}{a maximum number of iterations for the variance of
    random effects estimation, default is 50.}
  \item{tol}{tolerance for the convergence criterion. Default is 1e-6.     }
  \item{tol.theta}{tolerance for the convergence criterion (smoothing parameter function routine). Default is 1e-6.     }
  \item{returnFit}{a logical value indicating whether the fitted object
          should be returned when the maximum number of iterations is
          reached without convergence of the algorithm. Default is
          FALSE.}
\item{control}{see lmeControl in the documentation to
  \code{nlme}.}

}
\value{
  A list object of class \code{asp} containing the fitted model.
  The components are:
\item{fitted}{fitted values.}
\item{coef.mean}{estimated mean coefficients.}
\item{design.matrices}{design matrices both for knots und subknots.}
\item{x}{x values.}
\item{knots}{knots.}
\item{y.cov}{estimated covariance matrix of the response.}
\item{random.var}{estimated covariance matrix of the random effects.}
\item{subknots}{subknots.}
\item{coef.random}{estimated spline coefficients of the covariance
  matrix of the random effects.}
\item{var.random.var}{estimated variance of the spline coefficients of
  the covariance matrix of the random effects.}
  \item{fit}{mimics fit object of lme().}
  \item{info}{information about the inputs.}
  \item{aux}{auxiliary information such as variability
   estimates.}
}
\references{
Krivobokova, T., Crainiceanu, C.M. and Kauermann, G. (2008)\cr
Fast Adaptive Penalized Splines.\emph{ Journal of Computational and
Graphical Statistics. 17(1) 1-20.}

Ruppert, D., Wand, M.P. and Carroll, R.J. (2003)\cr
\emph{Semiparametric Regression} Cambridge University Press.\cr
\url{http://stat.tamu.edu/~carroll/semiregbook/}

Wiesenfarth, M., Krivobokova, T., Klasen, S., Sperlich, S. (2012).\cr
Direct Simultaneous Inference in Additive Models and its Application to Model Undernutrition.
\emph{Journal of the American Statistical Association}, accepted.
}
\author{
Manuel Wiesenfarth \email{mwiesen at gwdg de},
Tatyana Krivobokova \email{tkrivob at gwdg.de}
}
\seealso{
\code{\link{gam}} (in package `mgcv'),
\code{\link{asp}} (in package `AdaptFit'),
\code{\link{lme}} (in package `nlme')
}
\examples{
\dontrun{
############
# Examples as in package AdaptFit
  ## scatterplot smoothing
    x <- 1:1000/1000
    mu <- exp(-400*(x-0.6)^2)+5*exp(-500*(x-0.75)^2)/3+2*exp(-500*(x-0.9)^2)
    y <- mu+0.5*rnorm(1000)

  #fit with default knots
    y.fit <- asp2(y~f(x,adap=TRUE))
    plot(y.fit,residuals=TRUE,lwd=2,scb.lwd=2,scb.lty="dashed")
    # with shaded confidence region. Use scb.lty="blank" to plot the shades only.
      plot(y.fit,residuals=TRUE,shade=TRUE,scb.lty="blank")

  ## one more scatterplot smoothing with specified knots and subknots
    x <- 1:400/400
    mu <- sqrt(x*(1-x))*sin((2*pi*(1+2^((9-4*6)/5)))/(x+2^((9-4*6)/5)))
    y <- mu+0.2*rnorm(400)

    kn <- default.knots(x,80)
    kn.var <- default.knots(kn,20)

    y.fit <- asp2(y~f(x,knots=kn))
    y.fit2 <- asp2(y~f(x,knots=kn,var.knots=kn.var,adap=TRUE))
    op <- par(mfrow = c(1, 2))
    plot(y.fit)
    plot(y.fit2)
    par(op)


  ## additive models
    x1 <- 1:300/300
    x2 <- runif(300)
    mu1 <- exp(-400*(x1-0.6)^2)+5*exp(-500*(x1-0.75)^2)/3+2*exp(-500*(x1-0.9)^2)
    mu2 <- sin(2*pi*x2)
    y2 <- mu1+mu2+0.3*rnorm(300)

    y2.fit <- asp2(y2~f(x1,adap=TRUE)+f(x2,adap=TRUE))
    # switch off adaptive fitting for the first function
      y21.fit <- asp2(y2~f(x1,adap=FALSE)+f(x2,adap=TRUE))
    op <- par(mfrow = c(2, 2))
    plot(y2.fit)
    plot(y21.fit)
    par(op)

# Use predict to avoid centering of smooths in case of scatterplot smoothing
  library(SemiPar)
  data(fossil)
  attach(fossil)
  fit <- asp2(strontium.ratio~f(age,basis="tps",adap=FALSE))
  newdata.age <-data.frame(age=seq(90,130,length.out=50) )
  preds <- predict(fit,newdata=newdata.age,se=FALSE)
  plot(age,strontium.ratio)
  lines(newdata.age$age,preds,col="red")


##################
#more examples
  beta=function(l,m,x) return(gamma(l+m)*(gamma(l)*gamma(m))^(-1)*x^(l-1)*(1-x)^(m-1))
  f1 = function(x) return((0.6*beta(30,17,x)+0.4*beta(3,11,x))*1/0.958)
  f2 = function(z) return((sin(2*pi*(z-0.5))^2)*1/.3535)
  f3 = function(z) return((exp(-400*(z-0.6)^2)+5/3*
                       exp(-500*(z-0.75)^2)+2*exp(-500*(z-0.9)^2))*1/0.549)

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

  # semiparametric model
    fit1=  asp2(y~x1+f(x2,basis="os",degree=3,knots=kn2,adap=FALSE)
                    +f(x3,basis="os",degree=3,knots=kn3,var.knots=kn.var3,adap=FALSE),
                    niter = 20, niter.var = 200)
    summary(fit1)
    plot(fit1,pages=1)


  # all effects flexible
  # fit model with all smoothing parameters constant
    fit2a= asp2(y~f(x1,basis="os",degree=3,knots=kn1,adap=FALSE)
                  +f(x2,basis="os",degree=3,knots=kn2,adap=FALSE)
                  +f(x3,basis="os",degree=3,knots=kn3,adap=FALSE),
                  niter = 20, niter.var = 200)
    plot(fit2a,pages=1)

  # fit model with last smoothing parameter adaptive
    fit2b= asp2(y~f(x1,basis="os",degree=3,knots=kn1,adap=FALSE)
                  +f(x2,basis="os",degree=3,knots=kn2,adap=FALSE)
                  +f(x3,basis="os",degree=3,knots=kn3,adap=TRUE,
                     var.knots=kn.var3,var.basis="os",var.degree=3),
                 niter = 20, niter.var = 200)

  # plot smoothing parameter function for covariate x3.
  # Note that in the case of B-splines additional knots are added, see references.
    plot(seq(0,1,length.out=42), fit2b$y.cov/fit2b$random.var[85:126],
                ylab=expression(lambda(x3)),xlab="x3",type="l",lwd=3)

  # compute 95% simultaneous confidence bands.
  # You could skip this and use "fit2b" indstead of "scb2b" later on, however,
  # if N is large, computing the SCBs various times can take some time
  # if you don't need fitted values and bounds for all covariate points
  # (can be computationally intensive due to large matrix dimensions),
  # set calc.stdev=F such that these are not computed.
    scb2b<- scbM(fit2b,calc.stdev=FALSE)
    plot(scb2b,pages=1)

  # plot only f(x2).
    plot(scb2b,select=2,mfrow=c(1,1),lwd=3,ylab="f(x2)",xlab="x2")
  # plot.scbm (and plot.asp) returns fitted values and confidence limits,
  # if you only need the returned object set plot=FALSE
    pscb=plot(scb2b,plot=FALSE)
  # add pointwise confidence intervals to the plot
    polygon(c(pscb$grid.x[[2]], rev(pscb$grid.x[[2]])),
            c(pscb$fitted[[2]]+1.96*pscb$Stdev[[2]],
              rev(pscb$fitted[[2]]-1.96*pscb$Stdev[[2]])),
            col = grey(0.85), border = NA)
    lines(pscb$grid.x[[2]],pscb$lcb[[2]],lty="dotted",lwd=3)
    lines(pscb$grid.x[[2]],pscb$fitted[[2]],lwd=3)
    lines(pscb$grid.x[[2]],pscb$ucb[[2]],lty="dotted",lwd=3)

  # plot first derivative of f(x1)
    scb2bdrv<- scbM(fit2b,drv=1,calc.stdev=FALSE)
    plot(scb2bdrv,select=1)
    #the following would give the same result
    #x11();plot(fit2b,select=1,drv=1)
  # different style
    plot(scb2bdrv,select=1,scb.lty="blank", shade=TRUE,shade.col="steelblue")
}
}
\keyword{nonlinear} \keyword{models} \keyword{smooth} \keyword{regression}\keyword{adaptive}

