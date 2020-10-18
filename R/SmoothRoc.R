################################################
## Function:	smooth.roc
## Version:		1.0 
## Purpose:		Fit kernal density smoothed ROC
## 
## Author: 		JD Blume
## Date: 		April 2020
################################################
#'
#' Smooth an empirical ROC with kernal regression smoothers
#'
#' @description This function smooth the empirical ROC curve by smoothing the underlying density estimates for the two group. Kernal regression smoothers are used for the smoothing. 
#'
#' @param roc.obj ROC object from \code{set.roc} function.
#' @param adj Bandwith for kernal smoothers (used for both groups). Default is \code{0.8}.
#' @param acc Resolution along the x-axis (1-Specificity). Default is \code{0.01}.
#'
#' @details The bandwith for kernal smoothers is equal in the two groups.
#'
#' @return Returns a smoothed ROC curve (list) with the following elements:
#' 
#' \describe{
#' \item{\code{x}}{Survivor function for group 0.}
#' \item{\code{y}}{Survivor function for group 1.}
#' \item{\code{density.0}}{Kernal density estimate for group 0.}
#' \item{\code{density.1}}{Kernal density estimate for group 1.}
#' }
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
#' ## Fake ROC data for examples 
#' y.0 <- round( rnorm(100, mean=0, sd=1)   , 1) 
#' y.1 <- round( rnorm(100, mean=1, sd=1.1) , 1)
#' myroc <- set.roc(score=y.0, grp=y.1, as.groups=TRUE)
#'
#' ## Let's smooth
#' myroc.smooth <- smooth(myroc)
#'
#' ## Plotting
#' plot(myroc)
#' lines(myroc.smooth, col="purple", lwd=2)
#' lines(smooth(myroc,adj=0.5), col="hotpink", lwd=2)
#' 
#' ## Small bandwidth recovers empirical ROC curve
#' plot(myroc, roc.line="linear", line.col='red')
#' lines(smooth(myroc, adj=0.01), col="black")  
#'

smooth.roc <- function(roc.obj, adj=0.8, acc=0.1,...) {

#### Errors / Warnings

	if (is(roc.obj,"roc")!=TRUE) {stop("Object not ROC class.")}

#### Identify data for each group
	
	y.0 <- roc.obj$data.0
	y.1 <- roc.obj$data.1

#### Extend range of k-density beyond data
	
	e.0 <- 0.3 * diff(range(y.0))
	e.1 <- 0.3 * diff(range(y.1)) 
	e   <- 0.3 * diff(range(roc.obj$data.all$score))

#### Smoothed k-densities
	kdens.0 <- density(y.0, adjust=adj, from=min(y.0)-e.0, to=max(y.0) + e.0)
	kdens.1 <- density(y.1, adjust=adj, from=min(y.1)-e.1, to=max(y.1) + e.1)

#### Compute Survivor functions
	z = seq(min(y.0,y.1)-e, max(y.0,y.1) + e, by=acc)

	mat.0 <- outer(X=z, Y=y.0, function(X,Y) {pnorm(X, mean=Y, sd=kdens.0$bw)} ) 
	mat.1 <- outer(X=z, Y=y.1, function(X,Y) {pnorm(X, mean=Y, sd=kdens.1$bw)} ) 
	
	s.0 <- 1-rowMeans(mat.0)
	s.1 <- 1-rowMeans(mat.1)

	outlist <- list("x"=s.0, "y"=s.1, "density.0"=kdens.0, "density.1"=kdens.1)

class(outlist) <- "curve"
return(outlist) 
}

###
##
#