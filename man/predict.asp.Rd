\name{predict.asp}
\alias{predict.asp}
\title{Semiparametric regression prediction.}
\description{Takes a fitted \code{asp} object produced by
\code{asp2} and obtains predictions at new data values.}
\usage{
\method{predict}{asp}(object,newdata,se,...)
}
\arguments{ 
  
\item{object}{a fitted \code{asp} object as produced by \code{asp2()}. Does not work with \code{basis="os"}.}

\item{newdata}{a data frame containing the values of the predictors at
               which predictions are required. The columns should
               have the same name as the predictors. Further, minima and maxima should currently coincide with those of the predictors.}

\item{se}{when this is TRUE standard error estimates are
          returned for each prediction. The default is FALSE.}

\item{...}{other arguments.}

}
\details{Takes a fitted \code{asp} object produced by
\code{asp2()} and obtains predictions at new data values
as specified by the `newdata' argument. If `se=TRUE' then
standard error estimates are also obtained.}

\value{If se=FALSE then a vector of predictions at `newdata' is returned.
       If se=TRUE then a list with components named `fit' and `se' is
       returned. The `fit' component contains the predictions. 
       The `se' component contains standard error estimates.}

\author{Manuel Wiesenfarth, based on implementation of M.P. Wand (package \code{SemiPar}).
}

\seealso{
\code{\link{plot.asp}}, \code{\link{summary.asp}}
}


\examples{
data(fossil,package="SemiPar")
attach(fossil)
fit <- asp2(strontium.ratio~f(age, basis="tps"))
newdata.age <- data.frame(age=c(90,100,110,120,130))
preds <-  predict(fit,newdata=newdata.age,se=TRUE)
print(preds)

# Use predict to avoid centering of smooths in case of scatterplot 
# smoothing
  fit <- asp2(strontium.ratio~f(age,basis="tps"))
  newdata.age <-data.frame(age=seq(90,130,length.out=50) )
  preds <- predict(fit,newdata=newdata.age,se=TRUE)
  plot(age,strontium.ratio)
  lines(newdata.age$age,preds$fit,col="red")
  lines(unlist(newdata.age),preds$fit+2*preds$se,col="blue")
  lines(unlist(newdata.age),preds$fit-2*preds$se,col="green")
}
\keyword{models} \keyword{smooth} \keyword{regression} 


