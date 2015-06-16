\name{summary.asp}
\alias{summary.asp}
\title{Semiparametric regression summary}
\description{Takes a fitted \code{asp} object produced by
\code{asp2()} and summarises the fit. This function is extensively based on the summary function in package \code{SemiPar}.}
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
freedom values and other information.}

\value{The function generates summary tables.}

\references{
Ruppert, D., Wand, M.P. and Carroll, R.J. (2003)\cr
\emph{Semiparametric Regression} Cambridge University Press.\cr
\url{http://stat.tamu.edu/~carroll/semiregbook/}

Wiesenfarth, M., Krivobokova, T., Klasen, S., Sperlich, S. (2012).\cr
Direct Simultaneous Inference in Additive Models and its Application to Model Undernutrition.
\emph{Journal of the American Statistical Association}, 107(500): 1286-1296.
}

\seealso{
\code{\link{plot.asp}}, \code{\link{predict.asp}}

\code{\link{asp}} (package AdaptFit)

}

\examples{
#\dontrun{
library(SemiPar)
data(onions)
attach(onions)
log.yield <- log(yield)
fit <- asp2(log.yield~location+f(dens, degree=c(3,2)))
summary(fit,test1=TRUE,test2=TRUE)
#}
}
\keyword{models} \keyword{smooth} \keyword{regression} 


