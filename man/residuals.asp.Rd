\name{residuals.asp}
\alias{residuals.asp}
\title{Residuals for semiparametric regression.}
\description{Extracts residuals from a semiparametric regression fit object.
This function is nearly identical to the respective function in package \code{SemiPar}. }
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
\code{\link{plot.asp}}, \code{\link{predict.asp}}, \code{\link{summary.asp}}, \code{\link{fitted.asp}}

\code{\link{asp}} (package AdaptFit)
}

\examples{
library(SemiPar)
data(fossil)
attach(fossil)
fit <- asp2(strontium.ratio~f(age))
plot(age,residuals(fit))
abline(0,0)
}
\keyword{models} \keyword{smooth} \keyword{regression} 


