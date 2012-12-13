\name{fitted.asp}
\alias{fitted.asp}
\title{Fitted values for semiparametric regression.}
\description{Extracts fitted values from a semiparametric 
regression fit object. This function is identical to the respective function in package \code{SemiPar}.}
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

\references{
     Ruppert, D., Wand, M.P. and Carroll, R.J. (2003)\cr
     \emph{Semiparametric Regression} Cambridge University Press.\cr
     \url{http://stat.tamu.edu/~carroll/semiregbook/}


     Ganguli, B. and Wand, M.P. (2005)\cr
     \emph{SemiPar 1.0 Users' Manual}.\cr
     \url{http://www.uow.edu.au/~mwand/SPmanu.pdf}
}
\author{M.P. Wand 
\email{mwand@uow.edu.au}
(other contributors listed in SemiPar Users' Manual).}

\seealso{
\code{\link{plot.asp}},
\code{\link{predict.asp}},
\code{\link{summary.asp}},
\code{\link{residuals.asp}},
\code{\link{asp}} (package AdaptFit)
}

\examples{
\dontrun{
library(SemiPar)
data(fossil)
attach(fossil)
fit <- asp2(strontium.ratio~f(age))
plot(fit, bands=FALSE)
points(age,fitted(fit)-fit$coef[1],col="red")
}
}
\keyword{models} \keyword{smooth} \keyword{regression} 


