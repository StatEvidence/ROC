################################################
## Function:	auc.other
## Version:		1.0 
## Purpose:		Numerically Integrate ROC curve 
## 
## Author: 		JD Blume
## Date: 		May 2020
################################################
#'
#' Area under the ROC curve (AUC) by numerical integration 
#'
#' @description This function computes the area under an estimated ROC curve by numerical integration. If the curve can be drawn, this function will numerically integrate it.
#'
#' @param curve.obj Curve object. For example, from \code{smooth} function.
#' @param lo.lim lower limit of integration on x-axis (1-Specificity). Default is \code{lo.lim=0}.
#' @param hi.lim Upper limit of integration on x-axis (1-Specificity). Default is \code{hi.lim=1}.
#' @param rescale Flag for rescaling the partial ACU by the maximum area in rectangle of interest. Default is \code{TRUE}.
#' @param stop.err Flag for integrate function to "stop.on.error". Default is \code{FALSE}.
#'
#' @details Can compute the partial area under the ROC curve (pAUC) by setting \code{lo.lim} and/or \code{hi.lim}. This is done by numerical integration using R's \code{integrate} function and the \code{approx} function for interpolation between points on the curve. Denser curve objects are less dependent on the interpolation.
#'
#' @return Returns the area under the smoothed ROC curve (AUC). Partial AUC if limits of integration are set:
#' 
#' \describe{
#' \item{\code{auc}}{Area under the ROC curve (AUC).}
#' }
#' 
#' @author Jeffrey D Blume, \email{j.blume@@vanderbilt.edu}
#' @references Add references here
#' @seealso	\url{www.statisticalevidence.com} \code{\link{auc.roc}}
#' @keywords ROC curve, AUC
#' 
#' @export
#'
#' @examples
#'
#' ## Fake ROC data for examples 
#' y.0 <- round(rnorm(100, mean=0, sd=1)   , 1) 
#' y.1 <- round(rnorm(100, mean=1, sd=1.1) , 1)
#' myroc <- set.roc(score=y.0, grp=y.1, as.groups=TRUE)
#' 
#' ## Area under smoothed ROC curve
#' auc(smooth(myroc))
#'
#' ## Area under empirical ROC curve
#' auc(smooth(myroc, adj=0.01))
#'
#' ## Partial AUC
#' auc(smooth(myroc), lo.lim=0, hi.lim=0.2)
#' auc(smooth(myroc), lo.lim=0, hi.lim=0.2, rescale=FALSE)
#' auc(smooth(myroc, adj=0.01), lo.lim=0, hi.lim=0.2)
#' 

auc.other <- function(fit.obj, lo.lim=0, hi.lim=1, rescale=TRUE, stop.err=FALSE, level=0.95, ...) {

#### Errors / Warnings

	if (!missing(fit.obj)) {
		if (is(fit.obj, "fitobj")!=TRUE) {stop("Object not FITOBJ class.")}
		}
		
	if (lo.lim>=hi.lim) {stop("Lower limit of integration must be less than higher limit")}

#### Create curve function

	roc.curve <- function(new.pt, actual.pts) { 
		ycoordinate=approx(actual.pts$fpf, actual.pts$tpf, xout=new.pt, ties=min)
		return(ycoordinate$y)
	}

#### Get curve object

	crv.obj <- curve(fit.obj, level=level)

#### Organize curve points

	smooth.pts.fit <- unique(data.frame( cbind(fpf=c(0,crv.obj$x,1), tpf=c(0,crv.obj$y,1)) ) )
	smooth.pts.lo  <- unique(data.frame( cbind(fpf=c(0,crv.obj$x,1), tpf=c(0,crv.obj$ci.lo,1)) ) )
	smooth.pts.hi  <- unique(data.frame( cbind(fpf=c(0,crv.obj$x,1), tpf=c(0,crv.obj$ci.hi,1)) ) )

#### Integrate curve and obtain AUC

	int.auc.fit <- integrate(roc.curve, lower=lo.lim, upper=hi.lim, actual.pts=smooth.pts.fit, stop.on.error = stop.err)$value
	int.auc.lo  <- integrate(roc.curve, lower=lo.lim, upper=hi.lim, actual.pts=smooth.pts.lo, stop.on.error = stop.err)$value
	int.auc.hi  <- integrate(roc.curve, lower=lo.lim, upper=hi.lim, actual.pts=smooth.pts.hi, stop.on.error = stop.err)$value

#### Rescale partial AUC
	
	if (rescale==TRUE) { total.area <- (hi.lim - lo.lim) } else {total.area <- 1}

	int.auc.fit <- int.auc.fit/total.area
	int.auc.lo  <- int.auc.lo/total.area
	int.auc.hi  <- int.auc.hi/total.area

#### Results table

	auc <- c("Other", format(round(c(int.auc.fit, int.auc.lo, int.auc.hi), 4), nsmall = 2))

	tab <- data.frame(rbind(auc))
	colnames(tab) <- c("AUC       ", "Estimate ", "CI Lo  ", "CI Hi")

#### Return AUC

	outlist <- list("auc"=int.auc.fit, "auc.lo"=int.auc.lo, "auc.hi"= int.auc.hi, tab=tab, level=crv.obj$level)
	class(outlist) <- c("auc", class(outlist))

return(outlist) 
}

###
##
#