################################################
## Function:	smooth
## Version:		1.0 
## Purpose:		Smooth empirical ROC by kernal density 
## 
## Author: 		JD Blume
## Date: 		May 2020
################################################
#'
#' Generic function for smoothing ROC curves
#'
#' @description 
#'
#' @param x An ROC object. 
#'
#' @details 
#'
#' @return 
#' 
#' @author Jeffrey D Blume, \email{j.blume@@vanderbilt.edu}
#' @references Add references here
#' @seealso	\url{smooth.default} \code{\link{smooth.roc}}
#' @keywords ROC curve, smooth, density
#' 
#' @export
#'
#' @examples
#' 

smooth <- function(roc.obj,...) { UseMethod("smooth") }

smooth.default <- function (...) stats::ksmooth(...)

###
##
#