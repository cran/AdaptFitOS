\name{summary.asp}
\alias{summary.asp}
\title{Summaries and hypothesis tests}
\description{Takes a fitted \code{asp} object produced by
\code{asp2()} and summarises the fit, including tests for significance of nonparametric effects as well as their deviation from a parametric fit. 
}
\usage{
\method{summary}{asp}(object,test1=FALSE,test2=FALSE,signif=0.05,...)
}
\arguments{ 
  
\item{object}{a fitted \code{asp} object as produced by \code{asp2()}.}
\item{test1}{\code{TRUE} in order to include a test for significance of a nonparametrically estimated effect. The test correpsonds to checking whether the zero line is entirely inside the simultaneous confidence band.}
\item{test2}{\code{TRUE} in order to inlcude the nonparametric specification test proposed in Wiesenfarth et al. (2012). Only works with B-splines. The function under the null hypothesis is a polynomial of degree q-1 where q is the penalty order.}
\item{signif}{the significance level.}

\item{...}{other arguments.}

}
\details{Produces tables for the linear (parametric)
and non-linear (nonparametric) components. The linear
table provides coefficient estimates, standard errors
and p-values. The non-linear table provides degrees of
freedom values and other information including tests for significance of nonparametric effects as well as their deviation from a parametric fit.
See Wiesenfarth et al (2011, 2012) and Wiesenfarth (2012) for details on the hypothesis tests. 
}

\value{The function generates summary tables.}

\references{
Ruppert, D., Wand, M.P. and Carroll, R.J. (2003)\cr
\emph{Semiparametric Regression} Cambridge University Press.\cr
\url{https://web.stat.tamu.edu/~carroll/semiregbook/}

Wiesenfarth, M., Krivobokova, T., & Sperlich, S. (2011)\cr A Volume-of-tube based Test for Penalized Splines Estimators. \emph{Int. Statistical Inst.: Proc. 58th World Statistical Congress, 2011, Dublin.} \emph{http://www.2011.isiproceedings.org/papers/950754.pdf}


Wiesenfarth, M., Krivobokova, T., Klasen, S., Sperlich, S. (2012).\cr
Direct Simultaneous Inference in Additive Models and its Application to Model Undernutrition.
\emph{Journal of the American Statistical Association}, 107(500): 1286-1296.

Wiesenfarth, M. (2012). Estimation and Inference in Special Nonparametric Models. \emph{Doctoral dissertation, Goettingen, Georg-August Universitaet, Diss., 2012}. \emph{http://d-nb.info/104297182X/34}

}

\seealso{
\code{\link{plot.asp}}, \code{\link{predict.asp}}

}

\examples{
data(onions,package="SemiPar")
attach(onions)
log.yield <- log(yield)
fit <- asp2(log.yield~location+f(dens, degree=c(3,2)))
summary(fit,test1=TRUE,test2=TRUE)
}
\keyword{models} \keyword{smooth} \keyword{regression} 


