% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AucCurve.R
\name{auc.curve}
\alias{auc.curve}
\title{Area under the ROC curve (AUC) by numerical integration}
\usage{
\method{auc}{curve}(curve.obj, lo.lim = 0, hi.lim = 1, rescale = TRUE, stop.err = FALSE, ...)
}
\arguments{
\item{curve.obj}{Curve object. For example, from \code{smooth} function.}

\item{lo.lim}{lower limit of integration on x-axis (1-Specificity). Default is \code{lo.lim=0}.}

\item{hi.lim}{Upper limit of integration on x-axis (1-Specificity). Default is \code{hi.lim=1}.}

\item{rescale}{Flag for rescaling the partial ACU by the maximum area in rectangle of interest. Default is \code{TRUE}.}

\item{stop.err}{Flag for integrate function to "stop.on.error". Default is \code{FALSE}.}
}
\value{
Returns the area under the smoothed ROC curve (AUC). Partial AUC if limits of integration are set:

\describe{
\item{\code{auc}}{Area under the ROC curve (AUC).}
}
}
\description{
This function computes the area under an estimated ROC curve by numerical integration. If the curve can be drawn, this function will numerically integrate it.
}
\details{
Can compute the partial area under the ROC curve (pAUC) by setting \code{lo.lim} and/or \code{hi.lim}. This is done by numerical integration using R's \code{integrate} function and the \code{approx} function for interpolation between points on the curve. Denser curve objects are less dependent on the interpolation.
}
\examples{

## Fake ROC data for examples 
y.0 <- round(rnorm(100, mean=0, sd=1)   , 1) 
y.1 <- round(rnorm(100, mean=1, sd=1.1) , 1)
myroc <- set.roc(score=y.0, grp=y.1, as.groups=TRUE)

## Area under smoothed ROC curve
auc(smooth(myroc))

## Area under empirical ROC curve
auc(smooth(myroc, adj=0.01))

## Partial AUC
auc(smooth(myroc), lo.lim=0, hi.lim=0.2)
auc(smooth(myroc), lo.lim=0, hi.lim=0.2, rescale=FALSE)
auc(smooth(myroc, adj=0.01), lo.lim=0, hi.lim=0.2)

}
\references{
Add references here
}
\seealso{
\url{www.statisticalevidence.com} \code{\link{auc.roc}}
}
\author{
Jeffrey D Blume, \email{j.blume@vanderbilt.edu}
}
\keyword{AUC}
\keyword{ROC}
\keyword{curve,}
