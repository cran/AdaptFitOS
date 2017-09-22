\name{residuals.asp}
\alias{residuals.asp}
\title{Residuals for semiparametric regression.}
\description{Extracts residuals from a semiparametric regression fit object.}
\usage{
\method{residuals}{asp}(object,...)
}
\arguments{ 
  
\item{object}{a fitted \code{asp} object as produced by \code{asp2()}.}

\item{...}{other possible arguments.}

}
\details{Extracts residuals from a semiparametric regression fit object.
The residuals are defined to be the difference between the response
variable and the fitted values.}

\value{The vector of residuals.}


\seealso{
\code{\link{plot.asp}}, \code{\link{predict.asp}}, \code{\link{summary.asp}}, \code{\link{fitted.asp}}

\code{\link{asp}} (package AdaptFit)
}

\examples{
data(fossil,package="SemiPar")
attach(fossil)
fit <- asp2(strontium.ratio~f(age))
plot(age,residuals(fit))
abline(0,0)
}
\keyword{models} \keyword{smooth} \keyword{regression} 


