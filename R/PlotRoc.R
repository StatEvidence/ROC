################################################
## Function:	plot.roc
## Version:		2.1		(old version = ROC.plot)
## Purpose:		Set ROC plot; plot empirical pts
##				plot CI region for empirical pts
##
## Author: 		JD Blume
## Date: 		April 2020
################################################
#'
#' Plot Empirical ROC points and confidence region
#'
#' @description This function plots empirical ROC points and their pointwise confidence region. 
#'
#' @param roc.obj ROC object from \code{set.roc} function.
#' @param show.plot Flag for new plot. Default is \code{TRUE}.
#' @param show.points Flag to plot empirical ROC points. Default is \code{TRUE}.
#' @param roc.line Options for connecting empirical ROC points are \code{'step.fun'}, \code{'linear'}, \code{'none'}. . Default is \code{'none'}. See details.
#' @param ci.region flag for plotting pointwise CI region based on empirical pts. Default is \code{FALSE}.
#' @param flip.xaxis Flag to flip x-axis. Will show specificity instead of false positive rate (1-specificity). Default is \code{FALSE}.
#' @param points.col Color for points. Default is \code{'black'}.
#' @param points.pch Plotting symbol for points. Default is \code{pch=20}.
#' @param line.col color for ROC line. Default is \code{'black'}.
#' @param line.type line style for line for ROC line. Default is \code{lty=1}.
#' @param line.lwd line width for line for ROC line. Default is \code{lwd=2}.
#' @param ci.border flag for drawing border lines of CI region. Default is \code{NA}.
#' @param ci.color set color of CI region. Default is \code{'aliceblue'}.
#'
#' @details The ROC line can be set to either a step function, linear interpolation between the points (connect the dots), or none. The linear interpolation corresponds with empirical estimates of the area under the ROC curve, but can be misleading when there are few categories or many ties in the data.
#'
#' The confidence interval region is essentially the union of the confidence regions over all empirical points. Recall that each point has uncertainty along both the y- and x-axes. The CI limits are obtained when the ROC object is created with the \code{set.roc} function. If these were not computed, the function will not accept the \code{ci.region=TRUE} setting.
#' 
#' @return 
#' 
#' @author Jeffrey D Blume, \email{j.blume@@vanderbilt.edu}
#' @references Add references here
#' @seealso	\url{www.statisticalevidence.com} \code{\link{smooth.roc}}
#' @keywords ROC curve, AUC
#' 
#' @export
#'
#' @examples
#'
#' ## Fake ROC data for examples (with ties)
#' y.0 <- round( rnorm(100, mean=0, sd=1)   , 1)
#' y.1 <- round( rnorm(100, mean=1, sd=1.1) , 1)
#' myroc <- set.roc(score=y.0, grp=y.1, as.groups=TRUE)
#'
#' ## Various plots
#' plot(myroc)
#' plot(myroc, flip.xaxis=TRUE, ci.region=TRUE)
#' plot(myroc, roc.line="linear")
#' plot(myroc, roc.line="step.fun", line.col="green", show.plot=FALSE)
#'
#' ## CI region shows the span of the CIs around the empirical ROC points
#' plot(myroc, roc.line="linear", line.col="blue", points.col="blue", show.points=TRUE, ci.region=TRUE, ci.border="black")
#' points(myroc$fp.ci[,1], myroc$tp, col='hotpink',pch=20) 		## lower CI limit on FPR
#' points(myroc$fp, myroc$tp.ci[,2], col='purple',pch=20) 		## Upper CI limit on TPR
#' points(myroc$fp.ci[,2], myroc$tp, col='hotpink',pch=20)		## Upper CI limit on FPR
#' points(myroc$fp, myroc$tp.ci[,1], col='purple',pch=20)		## lower CI limit on TPR
#' legend('bottomright', c("Empirical ROC", "Edge of CI region", "FP CI limit", "TP CI limit"), col=c("blue", "black", "hotpink", "purple"), pch=c(20,NA,20,20), lty=c(1,1,NA,NA), bty='n')
#'

plot.roc <- function(roc.obj, show.plot=TRUE, show.points=TRUE, roc.line='none',
				ci.region=FALSE, flip.xaxis = FALSE,
				points.col='black', points.pch=20, 
				line.col='black', line.type=1, line.lwd=1.5,
				ci.border=NA, ci.color='aliceblue',...){

#### Errors / Warnings

	if (is(roc.obj,"roc")!=TRUE) {stop("Object not ROC class.")}

#### Get CI limits for ROC curve
	if (ci.region==TRUE) {

		## 	 Compute maximum CI region from CI limits above the cruve

		top <- data.frame(rbind(cbind(x=roc.obj$fp.ci[,1],y=roc.obj$tp),
		cbind(x=roc.obj$fp,y=roc.obj$tp.ci[,2])))

		top 	 <-	top[order(top[,'x'],top[,'y']),]
		top$ymax <- cummax(top[,'y'])

		##   Compute minimum CI region from CI limits below the cruve

		bot 	 <- data.frame(rbind(cbind(x=roc.obj$fp.ci[,2],y=roc.obj$tp),
		cbind(x=roc.obj$fp,y=roc.obj$tp.ci[,1])))
		bot 	 <- bot[order(-bot[,'x'],-bot[,'y']),]
		bot$ymin <- cummin(bot[,'y'])

		##   Compute ROC CI region using polygon function

		region <- data.frame(x=c(top$x, bot$x), y=c(top$ymax, bot$ymin))
		}

#### Set plot

	if (show.plot==TRUE) {
		plot(roc.obj$fp, roc.obj$tp, type="n", 
				pty="s", xaxt="n", yaxt="n",
				ylab="Sensitivity", xlab="", 
				xlim=c(0,1), ylim=c(0,1))
		
		axis(side=2, at=seq(0,1,0.2), 
				labels=format(seq(0,1,0.2),digits=2),
				las=1)
		
		if (flip.xaxis==TRUE) {
			axis(side=1, at=seq(0,1,0.2), labels=format(seq(1,0,-0.2), digits=2))
			mtext(side=1, line=3, "Specificity")
			}		
		else{ axis(side=1); mtext(side=1, line=3, "1-Specificity")}

	}

#### Draw CI region 
		
	if (ci.region==TRUE) {
		polygon(region$x, region$y, col=ci.color, border=ci.border)
	}

#### Draw reference line 

	if (show.plot==TRUE) {abline(0, 1, lty=2, lwd=0.5)}

#### Add empirical ROC points
	
	if (show.points==TRUE) {
		points(roc.obj$fp, roc.obj$tp, pch=points.pch, col=points.col)
	}

#### Add ROC line

	if (roc.line=='linear') {
		lines(roc.obj$fp, roc.obj$tp, lty=line.type, col=line.col, lwd=line.lwd)
	}
	
	if (roc.line=="step.fun") {
		
		step.fp <- sort(c(roc.obj$fp, roc.obj$fp), decreasing=FALSE)
		step.tp <- sort(c(roc.obj$tp, roc.obj$tp[-which.max(roc.obj$tp)],0), decreasing=FALSE)

		lines(step.fp, step.tp, col=line.col, lwd=line.lwd, lty=line.type)
	}
	
}

###
##
#