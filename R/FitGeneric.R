################################################
## Function:	fit
## Version:		1.0 
## Purpose:		Parent function for Fitting computations 
## 
## Author: 		JD Blume
## Date: 		May 2020
################################################
#'
#' Generic function for ROC model fitting
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
#' @seealso	\url{fit.roc} 
#' @keywords ROC curve, AUC
#' 
#' @export
#'
#' @examples
#' 

fit <- function(roc.obj,...) { UseMethod("fit") }

###
##
#