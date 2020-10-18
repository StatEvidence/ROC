################################################
## Function:	auc
## Version:		1.0 
## Purpose:		Parent function for AUC computations 
## 
## Author: 		JD Blume
## Date: 		May 2020
################################################
#'
#' Generic function for AUC computations
#'
#' @description 
#'
#' @param x An ROC or Curve object. 
#'
#' @details 
#'
#' @return 
#' 
#' @author Jeffrey D Blume, \email{j.blume@@vanderbilt.edu}
#' @references Add references here
#' @seealso	\url{auc.curve} \code{\link{auc.roc}}
#' @keywords ROC curve, AUC
#' 
#' @export
#'
#' @examples
#' 

auc <- function(x,...) { UseMethod("auc") }

###
##
#