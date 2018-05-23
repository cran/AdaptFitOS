\name{plot.asp}
\alias{plot.asp}
\alias{plot.scbm}
\title{Plots fitted curves or their derivatives}
\description{
Plots fitted curves or their derivatives together with simultaneous confidence bands.}
\usage{
\method{plot}{asp}(x, select=NULL, drv=0, bands=TRUE, level=0.95, grid=50, pages=0, 
				plot=TRUE, ylim=NULL, xlab=NULL, ylab=NULL,
     			scb.lwd=1, scb.lty="dotted", shade=FALSE, shade.col=grey(0.85),
     			residuals=FALSE, residuals.col="steelblue", 
     			bayes=FALSE, rug=TRUE,...)
}
\arguments{
  \item{x}{an \code{asp} object created by \code{\link{asp}} or \code{\link{aspHetero}}}
  \item{select}{vector specifying which curves in an additive model should be plotted. If NULL, all curves are plotted.}
   \item{drv}{the derivative order. Defaults to 0, i.e. the estimated curves themselves are plotted. First and second derivatives are supported. Does not apply to objects created by \code{\link{scbM}}.}
  \item{bands}{TRUE in order to include simultaneous confidence bands.}
  \item{grid}{number of points used for the plot, default value 50.}
  \item{plot}{if FALSE no plot is given}
  \item{ylim}{vector or list of vecotrs of limits on y axes. If NULL limits are automatically chosen.
   If multiple curves are plotted and a two-dimensional vector is given, y axes for all curves will be equal.
    A list with length equal to the number of smooth curves in the model can be given to specify different y-axes for each smooth.
  }
\item{pages}{The number of pages over which to spread the output as in package \code{\link{mgcv}}.
For example, if pages=1 then all terms will be plotted on one page in an automatic fashion.
If pages=0 (default) all graphics settings are left as they are. }
%  \item{mfrow}{the layout in the case of multiple functions as in par(mfrow=c(a,b)).
%  If NULL and not specified using \code{par()}, all curves are plotted in one window in an automatic fashion. Iff you want to plot the curves in an additive model on single pages e.g. in a pdf, use mfrow=c(1,1).}
  \item{level}{the level of confidence (does not apply to objects created by \code{\link{scbM}}).}
  \item{xlab}{label for the x axis. A list with length equal to the number of smooth curves in the model can be given to specify different labels for each smooth.}
  \item{ylab}{label for the y axis. A list with length equal to the number of smooth curves in the model can be given to specify different labels for each smooth.}
  \item{scb.lwd}{line width for simultaneous confidence bands}
  \item{scb.lty}{line type for simultaneous confidence bands. Use scb.lty="blank", if you only want to plot the shades.}
  \item{shade}{set to TRUE to produce shaded regions as simultaneous confidence bands for smooths}
  \item{shade.col}{define the color used for shading confidence bands}
  \item{residuals}{if TRUE, partial residuals are added to the plot}
  \item{residuals.col}{color of partial residuals}
  \item{rug}{adds a rug representation (1-d plot) of the data to the plot.}
\item{bayes}{FALSE for simultaneous confidence bands with (approximate) frequentist coverage probability, TRUE for (approximate) Bayesian coverage probability. See Krivobokova et al. (2010) for details.}
  \item{\dots}{further arguments to be passed to plot()}
}

\details{
plot.asp() first calls \code{\link{scbM}} and then plot.scbm() to plot an asp object.
If plotting takes long (because of a large data set) and you want to plot multiple times with different settings, use \code{\link{scbM}} and then plot the resulting \code{scbm} object with \code{plot.scbm()}.
Different to packages SemiPar and AdaptFit, estimated fits are centred to have zero mean.
The simultaneous confidence bands have (approximate) frequentist coverage probabilities with automatic bias correction (see references).
}

\value{
\item{grid.x }{list of the grid values used}
\item{fitted }{list of the fitted values on the grid}
\item{lcb }{list of the lower bounds of the confidence bands}
\item{ucb }{list of the upper bounds of the confidence bands}
\item{drv }{the derivative order}
\item{Stdev.fit }{the standard deviations on the grid}
\item{ylim }{list of ylim used for plotting}
\item{residuals }{the partial residuals.}
}

\references{
Krivobokova, T., Kneib, T., and Claeskens, G. (2010)\cr
Simultaneous confidence bands for penalized spline estimators.
\emph{Journal of the American Statistical Association}, 105(490):852-863.

Wiesenfarth, M., Krivobokova, T., Klasen, S., Sperlich, S. (2012).\cr
Direct Simultaneous Inference in Additive Models and its Application to Model Undernutrition.
\emph{Journal of the American Statistical Association}, 107(500): 1286-1296.

}

\author{
Manuel Wiesenfarth \email{m.wiesenfarth at dkfz de}
}

\seealso{
\code{\link[SemiPar]{plot.spm}} in package SemiPar
}
\examples{
# see asp2()
}
\keyword{aplot}
\keyword{dplot}
