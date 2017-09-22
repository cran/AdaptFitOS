\name{AdaptFitOS-package}
\alias{AdaptFitOS-package}
\alias{AdaptFitOS}
\docType{package}
\title{
Adaptive Semiparametric Regression with Simultaneous Confidence Bands
}
\description{
Based on package \code{AdaptFit}, fits semiparametric regression models with spatially adaptive penalized splines and computes simultaneous confidence bands.
%Based on the function \code{\link[AdaptFit]{asp}} of the AdaptFit package (which is based on \code{\link{spm}} function of \code{SemiPar} package)
Particular differences to AdaptFit
include the availability of simultaneous confidence bands and B-spline basis functions and different functionality of the \code{plot} function. However, random effects, autocorrelations and interaction surfaces are not supported.
Further, only Gaussian responses are supported. Note that in contrast to AdaptFit, estimated curves are centered to have zero mean and dummies for categorical covariates are constructed automatically if
a variable is given as factor.
For computation of the critical value for simultaneous confidence bands based on Hotelling's volume-of-tube formula, some functions of the \code{libtube} library by Catherine Loader (see package \code{\link[locfit]{locfit}}) are used.
See the references for details on the construction of the confidence bands.
}
\details{
\tabular{ll}{
Package: \tab AdaptFitOS\cr
Version: \tab 0.42\cr
Date: \tab 2012-06-03\cr
Depends: \tab MASS, nlme, cluster\cr
}

Index:
\preformatted{
asp2             Fit a semiparametric regression model with
                 spatially adaptive penalized
aspFormula       An asp formula
aspHetero        Estimate varying residual variance
fitted.asp       Fitted values for semiparametric regression.
plot.asp         Plots fitted curves or their derivatives
                 including simultaneous confidence bands
plot.scbm        Plots fitted curves in a scbm object
                 including simultaneous confidence bands
predict.asp      Semiparametric regression prediction.
residuals.asp    Residuals for semiparametric regression.
scbM             Calculate simultaneous confidence bands for
                 penalized splines
summary.asp      Semiparametric regression summary
}
The function \code{asp2()} is used to fit the model. Using the resulting \code{asp} object, fitted curves or their derivatives can be plotted with \code{\link{plot.asp}} and information on the parametric effects can be printed using \code{\link{summary}}.

See Wiesenfarth et al (2012) for technical details and Wiesenfarth (2012, Chapter 5.1) for some more details on the use of the package (including a demonstration on how plots in Wiesenfarth et al are obtained).

}
\author{
Manuel Wiesenfarth and Tatyana Krivobokova

Maintainer: Manuel Wiesenfarth <m.wiesenfarth at dkfz de>
}
\references{

Krivobokova, T., Crainiceanu, C.M. and Kauermann, G. (2008)\cr
Fast Adaptive Penalized Splines.
\emph{ Journal of Computational and Graphical Statistics}, 17(1):1-20.

Krivobokova, T., Kneib, T., and Claeskens, G. (2010)\cr
Simultaneous confidence bands for penalized spline estimators.
\emph{Journal of the American Statistical Association}, 105(490):852-863.

Wiesenfarth, M., Krivobokova, T., Klasen, S., Sperlich, S. (2012).\cr
Direct Simultaneous Inference in Additive Models and its Application to Model Undernutrition.
\emph{Journal of the American Statistical Association}, 107(500): 1286-1296.

Wiesenfarth, M. (2012). Estimation and Inference in Special Nonparametric Models. \emph{Doctoral dissertation, Goettingen, Georg-August Universitaet, Diss., 2012}. \emph{http://d-nb.info/104297182X/34}

}


\keyword{ package }
\seealso{
\code{\link[SemiPar]{spm}} (package SemiPar),
\code{\link[AdaptFit]{asp}} (package AdaptFit)
}
