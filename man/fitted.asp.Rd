\name{fitted.asp}
\alias{fitted.asp}
\title{Fitted values for semiparametric regression.}
\description{Extracts fitted values from a semiparametric 
regression fit object. }
\usage{
\method{fitted}{asp}(object,...)
}
\arguments{ 
  
\item{object}{a fitted \code{asp} object as produced by \code{asp2()}.}

\item{...}{other possible arguments.}

}
\details{Extracts fitted from a semiparametric regression fit object.
The fitted are defined to be the set of values obtained
when the predictor variable data are substituted into the
fitted regression model.}

\value{The vector of fitted.}


\seealso{
\code{\link{plot.asp}},
\code{\link{predict.asp}},
\code{\link{summary.asp}},
\code{\link{residuals.asp}},
\code{\link{asp}} (package AdaptFit)
}

\examples{
#\dontrun{
library(SemiPar)
data(fossil)
attach(fossil)
fit <- asp2(strontium.ratio~f(age))
plot(fit, bands=FALSE)
points(age,fitted(fit)-fit$coef[1],col="red")
#}
}
\keyword{models} \keyword{smooth} \keyword{regression} 


