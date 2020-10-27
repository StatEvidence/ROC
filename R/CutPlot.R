################################################
## Function:	cutplot
## Version:		2.0		
## Purpose:		Plot CDFs and pdfs of model
## 
## Author: 		JD Blume
## Date: 		October 2020
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
#'

cutplot <- function(roc.obj, fit.obj=NULL, at,
						adj=0.5, ref=TRUE,
						cols = c("dodgerblue3", "chocolate3"),
						lwds = c(2.5, 2.5), 
						show.legend=TRUE, show.auc=TRUE,
						xlab="Score", ...){

##########
## Binormal model curves change to logistic based on model
##########


#### Errors / Warnings

	if ( missing(roc.obj) ) { stop("ROC object required.") }

	if ( !is(roc.obj, "roc") ) { stop("ROC object not ROC class.") }
	
	if ( !missing(fit.obj) ) { 
			
		if ( !is(fit.obj, "fitobj") ) { stop("Object not FITOBJ class.") }

	 	}

	
#### Identify data for each group

	y.0 <- roc.obj$data.0
	y.1 <- roc.obj$data.1
	
	e	<- 0.3 * diff(range(roc.obj$data.all$score))
		
	cut.lo <- min(c(y.0,y.1))-e
	cut.hi <- max(c(y.0,y.1))+e
	
	cutrange <- seq(cut.lo, cut.hi, 0.01)
	
	auc.emp <- auc(roc.obj)$auc
	auc.fit <- NA
	
#### Compute empirical densities	

	sroc <- smooth.roc(roc.obj, adj=adj)

	eFn.0 <- ecdf(y.0)
	eFn.1 <- ecdf(y.1)
	
	upper.lim <- max(sroc$density.0$y, sroc$density.1$y)

	
#### Extract ROC (a, b) and density (t0, t1) parameters from fit
	
	if (!missing(fit.obj)) {
			
		a  <- fit.obj$a$coef['(Intercept)']
		b  <- fit.obj$b$coef['(Intercept)']
		t0 <- fit.obj$t0$coef['(Intercept)']
		t1 <- fit.obj$t1$coef['(Intercept)']
		
		auc.fit <- ifelse(fit.obj$model$model=="binormal", 
								auc(fit.obj)$auc, 
									auc.other(fit.obj)$auc)

		}


#### Compute predicted ROC (a, b) and density (t0, t1) parameters from fit

	if (!missing(at)) {
	
		keep <- predict(fit.obj, at=at)
	
		a  <- keep$parms$a
		b  <- keep$parms$b
		t0 <- keep$parms$t0
		t1 <- keep$parms$t1
		
		auc.fit <- keep$auc$est
		}
		
		
#### Compute fitted density parameters

	if (!missing(fit.obj)) {
		
		mu.0 <- -t0/t1
		sd.0 <- sqrt( (1/t1)^2 )

		mu.1 <- -t0/t1 + a/(b*t1)
		sd.1 <- sqrt( 1/(b^2*t1^2) )
	
		}

#### Map back to densities 

	if (!missing(fit.obj)) {

	## Class 0
	if (fit.obj$model$dist.0=="normal") 	{ 
			fn.0 <- dnorm(cutrange, mean=mu.0, sd=sd.0)
			Fn.0 <- pnorm(cutrange, mean=mu.0, sd=sd.0)  
			}

	if (fit.obj$model$dist.0=="logistic") 	{ 
			fn.0 <- dlogis(cutrange, location=mu.0, scale=sd.0)
			Fn.0 <- plogis(cutrange, location=mu.0, scale=sd.0) 
			} 
		
	if (fit.obj$model$dist.0=="gamma") 		{
			fn.0 <- dgamma(cutrange, shape=mu.0, scale=sd.0)
			Fn.0 <- pgamma(cutrange, shape=mu.0, scale=sd.0)
			}
	
	## Class 1	
	if (fit.obj$model$dist.1=="normal") 	{ 
			fn.1 <- dnorm(cutrange, mean=mu.1, sd=sd.1)
			Fn.1 <- pnorm(cutrange, mean=mu.1, sd=sd.1)  
			}

	if (fit.obj$model$dist.1=="logistic") 	{ 
			fn.1 <- dlogis(cutrange, location=mu.1, scale=sd.1)
			Fn.1 <- plogis(cutrange, location=mu.1, scale=sd.1) 
			} 

	if (fit.obj$model$dist.1=="gamma") 		{
			fn.1 <- dgamma(cutrange, shape=mu.1, scale=sd.1)
			Fn.1 <- pgamma(cutrange, shape=mu.1, scale=sd.1)
			}
	
	upper.lim <- max(fn.0, fn.1, upper.lim)
	
		}
	
#### Plot cutpoint distribution as CDFs and pdfs 

	par(mfrow=c(2,1),
		mar=c(0.1, 0.75, 0.75, 0.5), 
		oma = c(4, 4, 0.2, 0.2)+0.1 
		)

	layout(matrix(1:2,ncol=1),heights=c(1,2))
	
	## Top plot (PDFs)
	plot(0, 0, 	ylim=c(0, upper.lim), type="n", 
				xlim=c(min(cutrange), max(cutrange)), 
				xlab=" ", ylab=" ", xaxt="n", las=1 
		)
	
	mtext(side=2, line=3, "Density")	
				
	grid (NULL, NULL, lty = 3, col = "grey81") 	
		
	if ( ref==TRUE ) { 			
		lines(sroc$density.0, lty=2, lwd=1.2, col="black") 
		lines(sroc$density.1, lty=2, lwd=1.2, col="black")
		}

	if ( !missing(fit.obj) ) { 	
		lines(cutrange, fn.0, lwd=lwds[1], col=cols[1])
		lines(cutrange, fn.1, lwd=lwds[2], col=cols[2])
		}
		
		
	if (show.auc==TRUE) {

		legend("topright", c(	paste("AUC emp:", round(auc.emp,4)),
								paste("AUC fit    :", round(auc.fit,4)) ), 
								bty="n")
		}	
		
	## Bottom plot (CDFs)
	plot(0, 0, 	ylim=c(0, 1), type="n", 
				xlim=c(min(cutrange), max(cutrange)), 
				xlab=" ", ylab=" ", las=1 
		)

	grid (NULL, NULL, lty = 3, col = "grey81") 
	
	if ( ref==TRUE ) { 			
		lines(cutrange, eFn.0(cutrange),   lty=2, lwd=1.2, col="black")
		lines(cutrange, 1-eFn.1(cutrange), lty=2, lwd=1.2, col="black")
		}

	if ( !missing(fit.obj) ) {
		lines(cutrange, Fn.0,   lwd=lwds[1], col=cols[1])
		lines(cutrange, 1-Fn.1, lwd=lwds[2], col=cols[2])
		}
		
	mtext(side=2, line=3, "Sensitivity / Specificity")	
	mtext(side=1, line=2.15, xlab)

	if (show.legend==TRUE & !missing(fit.obj)) {
 	
		legend("right", c("Group 0", "(Specificity)", "Group 1", "(Sensitivity)"), bty="n",
		 	lty=c(1, NA, 1, NA), lwd=c(3, NA, 3, NA) , col = c(cols[1], NA, cols[2], NA), 
			cex = 1)

		}

	par(mfrow=c(1,1))

	outlist <- list(
		graph= data.frame(score=cutrange, 
					sens.emp=eFn.0(cutrange), spec.emp=1-eFn.1(cutrange), 
					sens.fit=Fn.0, spec.fit=1-Fn.1), auc.emp=auc.emp, auc.fit=auc.fit)
	invisible(outlist)

}


###
##
#