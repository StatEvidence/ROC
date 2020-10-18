################################################
## Function:	hist.roc
## Version:		2.0		
## Purpose:		Plot histograms displaying data
##				and density curves (if supplied)
## 
## Author: 		JD Blume
## Date: 		April 2020
################################################
#'
#' Plot histograms of ROC data and estimated densities
#'
#' @description This function plots histograms of the ROC and associated smooth or fitted densities. 
#'
#' @param roc.obj ROC object from \code{set.roc} function.
#' @param curve.obj Curve object. For example, from \code{smooth} function.
#' @param breaks Histogram breaks. Default is \code{"Scott"} but \code{"FD"} is also a good option.
#' @param hcols Distogram colors for group 0 and group 1. Default is \code{c("steelblue1", "darkorange")}.
#' @param dcols Density line colors for group 0 and group 1. Default is \code{c("dodgerblue3", "chocolate3")}.
#' @param lwds  Density line widths for group 0 and group 1. Default is \code{c(3,3)}.
#' @param show.legend Flag to display default legend. Default is \code{TRUE}.
#' @param title Title
#' @param xlab Lable for x-axis.
#' @param blty Line type for smoothed densities when plotted with fitted densities (solid line). 
#' @param blwd Line width for smoothed densities when plotted with fitted densities (solid line).
#'
#' @details Details?
#'
#' @return 
#' 
#' @author Jeffrey D Blume, \email{j.blume@@vanderbilt.edu}
#' @references Add references here
#' @seealso	\url{www.statisticalevidence.com} \code{\link{plot.roc}}
#' @keywords ROC curve, Histogram
#' 
#' @export
#'
#' @examples
#'
#' ## Fake ROC data for examples 
#' y.0 <- round(rnorm(100, mean=0, sd=1)   , 1) 
#' y.1 <- round(rnorm(100, mean=1.75, sd=1.1) , 1)
#' myroc <- set.roc(score=y.0, grp=y.1, as.groups=TRUE)
#' 
#' myroc.smooth <- smooth(myroc)
#' 
#' ## Histogram of data and smoothed densities
#' hist(myroc)
#' hist(myroc, myroc.smooth)
#' 
#' ## More smooth, alternate colors, forced breaks
#' hist(myroc, smooth(myroc,adj=1), hcols=c("dodgerblue3","firebrick"), dcols=c("darkblue","darkred"), seq(-4, 6, length.out = 21))
#' 
#'

hist.roc <- function(roc.obj, curve.obj=NULL, fit.obj=NULL,
						breaks="Scott",
						hcols = c("steelblue1", "darkorange"),
						dcols = c("dodgerblue3", "chocolate3"),
						lwds = c(3, 3) ,
						show.legend=TRUE,
						title = "Density Estimate", xlab="Score", blty=3, blwd=2,...){

#### Errors / Warnings

	if ( !is(roc.obj, "roc") ) { stop("ROC object not ROC class.") }
	
	if ( !missing(fit.obj) ) { 
			
		if ( !is(fit.obj, "fitobj") ) { stop("Object not FITOBJ class.") }

	}
	
	if ( !missing(curve.obj) ) { 
	
		if ( !is(curve.obj, "curve") ) { 
		
			if ( !is(curve.obj, "fitobj") | !missing(fit.obj) ) {stop("Curve Object not CURVE class.") } 
				
			fit.obj <- curve.obj 
			curve.obj <- substitute()
	
		}
		
	}	

#### Identify data for each group

	y.0 <- roc.obj$data.0
	y.1 <- roc.obj$data.1
	e	<- 0.3 * diff(range(roc.obj$data.all$score))
	top <- max(hist(y.0, breaks=breaks)$density, hist(y.1, breaks=breaks)$density)
	
	x.lo <- min(c(y.0,y.1))-e
	x.hi <- max(c(y.0,y.1))+e
	
	xrange <- seq(x.lo, x.hi, 0.01)
	
#### Compute fitted densities
	
	if (!missing(fit.obj)) {
			
		a  <- fit.obj$a$coef
		b  <- fit.obj$b$coef
		t0 <- fit.obj$t0$coef
		t1 <- fit.obj$t1$coef

		mu.0 <- -t0/t1
		sd.0 <- sqrt( (1/t1)^2 )

		mu.1 <- -t0/t1 + a/(b*t1)
		sd.1 <- sqrt( 1/(b^2*t1^2) )

		fn.0 <- dnorm(xrange, mean= mu.0, sd=sd.0)
		fn.1 <- dnorm(xrange, mean= mu.1, sd=sd.1)
			
		}
	
#### Plot Histograms

	tcolor <- c(col2rgb(hcols[2],alpha=TRUE)/255)

	hist(y.1, breaks=breaks, freq=FALSE,
			xlim=c( x.lo, x.hi ),
			ylim=c(0, top),
			col=hcols[2],
			main=title, xlab="") 
			mtext(xlab, 1, line=2.25)

	hist(y.0, breaks=breaks, freq=FALSE, add=TRUE,
			col=hcols[1]) 

	hist(y.1, breaks=breaks, freq=FALSE,
			add=TRUE,
			col=rgb(tcolor[1], tcolor[2], tcolor[3], 0.5))   # 50% transparent 

#### Add smooth density estimates	

	if (!missing(curve.obj) & missing(fit.obj)) {
		lines(curve.obj$density.0, lwd=3, col=dcols[1]) 
		lines(curve.obj$density.1, lwd=3, col=dcols[2])
		}

	if (missing(curve.obj) & !missing(fit.obj)) {
		lines(xrange, fn.0, lwd=3, col=dcols[1])
		lines(xrange, fn.1, lwd=3, col=dcols[2])
		}
	
	if (!missing(curve.obj) & !missing(fit.obj)) {
		lines(curve.obj$density.0, lwd=blwd, col=dcols[1], lty=blty) 
		lines(curve.obj$density.1, lwd=blwd, col=dcols[2], lty=blty)
	
		lines(xrange, fn.0, lwd=3, col=dcols[1])
		lines(xrange, fn.1, lwd=3, col=dcols[2])
		}
	
	if (show.legend==TRUE) {

		legend("bottom", c("Group 0", "Group 1"), xpd = TRUE, 
			horiz = TRUE, inset = c(min(c(y.0, y.1))-e, -0.19), 
			bty = "n", lty=1, lwd=3, col = c(dcols[1], dcols[2]), 
			cex = 1)
		}
}

###
##
#