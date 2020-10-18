################################################
## Function:	curve
## Version:		1.0 
## Purpose:		Obtain ROC curve from fitted object 
## 
## Author: 		JD Blume
## Date: 		May 2020
################################################
#'
#' Generic function for obtaining ROC curves from fitted ROC objects
#'
#' @description 
#'
#' @param x An ROC fitted object. 
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

curve <- function(fit.obj,...) { UseMethod("curve") }

curve.default <- function (...) graphics::curve(...)

###
##
#